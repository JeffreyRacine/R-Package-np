# An active entry owns metadata only, never result payloads or worker closures.
# Cancellation is cooperative at existing rank-local task/chunk boundaries.
# This protocol does not recover a lost rank or a failed native MPI primitive.
.npRmpi_fanout_state <- new.env(parent = emptyenv())
.npRmpi_bootstrap_task_progress <- new.env(parent = emptyenv())

# Preserve the original condition class, including message-less interrupts.
.npRmpi_fanout_rethrow <- function(condition) {
  if (inherits(condition, "interrupt")) {
    signalCondition(condition)
    invokeRestart("abort")
  }
  stop(condition)
}

# Optional, call-scoped completion within one task. Consumers report a
# cumulative count only after completing their numerical batch. The owner
# converts it to deltas and credits any unreported tail exactly once.
.npRmpi_bootstrap_task_checkpoint <- function(done) {
  checkpoint <- .npRmpi_bootstrap_task_progress$checkpoint
  if (is.function(checkpoint)) checkpoint(done)
  invisible(NULL)
}

.npRmpi_bootstrap_task_evaluate <- function(expr, weight, report = NULL) {
  if (is.null(report)) return(force(expr))
  state <- new.env(parent = emptyenv())
  state$done <- 0L
  previous <- .npRmpi_bootstrap_task_progress$checkpoint
  on.exit(.npRmpi_bootstrap_task_progress$checkpoint <- previous, add = TRUE)
  checkpoint <- function(done) {
    if (length(done) != 1L || is.na(done) || !is.finite(done) ||
        done != as.integer(done) || done < state$done || done > weight)
      stop("invalid internal bootstrap completion count", call. = FALSE)
    delta <- as.integer(done - state$done)
    if (delta > 0L) report(delta)
    state$done <- as.integer(done)
    invisible(NULL)
  }
  .npRmpi_bootstrap_task_progress$checkpoint <- checkpoint
  value <- force(expr)
  if (!inherits(value, "try-error") && !inherits(value, "npRmpi_local_failure"))
    checkpoint(weight)
  value
}

.npRmpi_fanout_native <- function(action, comm, arg1 = NULL, arg2 = NULL) {
  comm <- as.integer(comm)
  # Literal entry points keep registration/arity checking visible to R CMD
  # check. This dispatcher accepts only the five private protocol actions.
  switch(action,
    begin = .Call("np_mpi_fanout_begin", comm, arg1, PACKAGE = "npRmpi"),
    send = .Call("np_mpi_fanout_send", comm, arg1, arg2, PACKAGE = "npRmpi"),
    poll = .Call("np_mpi_fanout_poll", comm, PACKAGE = "npRmpi"),
    finish = .Call("np_mpi_fanout_finish", comm, PACKAGE = "npRmpi"),
    owner = .Call("np_mpi_fanout_owner", comm, PACKAGE = "npRmpi"),
    stop("unknown private MPI fan-out control action", call. = FALSE))
}

.npRmpi_fanout_retained <- function(comm) {
  key <- as.character(comm)
  if (exists(key, envir = .npRmpi_fanout_state, inherits = FALSE))
    return(get(key, envir = .npRmpi_fanout_state, inherits = FALSE))
  if (!isTRUE(getOption("npRmpi.mpi.initialized", FALSE))) return(NULL)
  .npRmpi_fanout_native("owner", comm)
}

.npRmpi_fanout_assert_idle <- function(comm, activation = FALSE) {
  tx <- .npRmpi_fanout_retained(comm)
  if (is.null(tx) || tx$phase %in% c("prepared", "quiescent")) return(invisible(TRUE))
  if (activation && isTRUE(tx$owner.active) && identical(tx$phase, "activating"))
    return(invisible(TRUE))
  stop("MPI pool has active or quarantined fan-out work; no new work was dispatched. For an abandoned call, retry npRmpi.quit(force = TRUE) to complete cooperative cleanup.",
       call. = FALSE)
}

.npRmpi_fanout_notice <- function(text) {
  # Safety disposition is not ordinary estimator progress, and warn=2 or a
  # user message handler must not replace the original condition.
  tryCatch(.np_message(text), error = function(e) NULL, interrupt = function(e) NULL)
  invisible(NULL)
}

# Only the canonical bootstrap caller supplies this data-only, call-scoped
# owner reference. It is never part of a transaction or a worker payload.
.npRmpi_fanout_recovery_context <- new.env(parent = emptyenv())
.npRmpi_fanout_recovery_context$current <- NULL

.npRmpi_with_fanout_recovery <- function(owner, code) {
  previous <- .npRmpi_fanout_recovery_context$current
  on.exit({ .npRmpi_fanout_recovery_context$current <- previous }, add = TRUE)
  .npRmpi_fanout_recovery_context$current <- owner
  force(code)
}

.npRmpi_fanout_recovery_owner <- function(comm) {
  owner <- .npRmpi_fanout_recovery_context$current
  if (is.list(owner) && identical(owner[["comm"]], comm)) owner else NULL
}

.npRmpi_fanout_recovery_step <- function(owner, first = FALSE) {
  if (is.null(owner)) return(invisible(NULL))
  disable <- function(e) {
    state <- owner[["ref"]][[owner[["slot"]]]]
    if (is.list(state)) {
      state$enabled <- FALSE
      owner[["ref"]][[owner[["slot"]]]] <- state
    }
    invisible(NULL)
  }
  tryCatch({
    state <- owner[["ref"]][[owner[["slot"]]]]
    if (!isTRUE(getOption("np.messages", TRUE)) || !is.list(state) ||
        !isTRUE(state[["enabled"]]) || !isTRUE(state[["visible"]]) ||
        isTRUE(state[["message_muffled"]]) ||
        !identical(state[["id"]], .np_progress_registry$active_id))
      return(invisible(NULL))
    if (first) {
      state$label <- "Waiting for MPI cleanup"
      state$known_total <- FALSE
      state$total <- NULL
      state$last_done <- NULL
      state$last_emitted_done <- NULL
      state$unknown_total_fields <- NULL
      state$started <- .np_progress_now()
      state$start_note_pending <- FALSE
      state$fanout_recovery <- TRUE
      owner[["ref"]][[owner[["slot"]]]] <- state
    }
    if (!isTRUE(state[["fanout_recovery"]])) return(invisible(NULL))
    owner[["ref"]][[owner[["slot"]]]] <- .np_progress_step_at(
      state, now = .np_progress_now(), force = first)
    invisible(NULL)
  }, error = disable, interrupt = disable)
}

.npRmpi_fanout_new <- function(comm, workers, n, scheduler, weights = NULL) {
  key <- as.character(comm)
  if (!is.null(.npRmpi_fanout_retained(comm)))
    stop("an earlier MPI fan-out has not reached quiescence", call. = FALSE)
  .npRmpi_protocol_rank_tag("manual_bcast_base", as.integer(workers),
                           where = "MPI fan-out")
  .npMacMseriesAccelerateOptionValue()
  .npRmpi_fanout_metadata(comm, workers, n, scheduler, weights,
    .npRmpi_lease_ensure_session(), .npRmpi_lease_next_id("fanout"),
    .npRmpi_protocol_tag("fanout_control"))
}

.npRmpi_fanout_metadata <- function(comm, workers, n, scheduler, weights,
                                    session, operation, control) {
  tx <- new.env(parent = emptyenv())
  tx$comm <- comm
  tx$key <- as.character(comm)
  tx$session <- session
  tx$id <- operation
  tx$scheduler <- scheduler
  tx$control <- control
  tx$n <- as.integer(n)
  tx$workers <- as.integer(workers)
  tx$phase <- "prepared"
  tx$owner.active <- FALSE
  tx$native <- FALSE
  tx$cancelled <- FALSE
  tx$cleanup.started <- NA_real_
  tx$cleanup.budget <- 0
  tx$rank <- rep.int("initial", workers)
  tx$assigned <- rep(list(integer()), workers)
  tx$tags <- integer(workers)
  tx$scheduled <- 0L
  tx$next.task <- integer(workers)
  tx$seen <- rep.int(FALSE, n)
  tx$weights <- if (is.null(weights)) rep.int(1L, n) else weights
  tx$progress <- integer(workers)
  tx$raw.messages <- 0L
  tx$raw.bytes <- 0
  tx
}

.npRmpi_fanout_header <- function(tx) {
  list(session = tx$session, operation = tx$id, scheduler = tx$scheduler,
       control = tx$control, comm = tx$comm)
}

.npRmpi_fanout_envelope <- function(header, kind, tasks = integer(), value = NULL) {
  list(session = header$session, operation = header$operation,
       kind = kind, tasks = as.integer(tasks), value = value)
}

.npRmpi_fanout_identity <- function(tx, message, source) {
  is.list(message) && !anyDuplicated(names(message)) &&
    all(c("session", "operation", "kind", "tasks", "value") %in% names(message)) &&
    identical(message[["session"]], tx$session) &&
    identical(message[["operation"]], tx$id) && length(source) == 1L &&
    !is.na(source) && source >= 1L && source <= tx$workers
}

.npRmpi_fanout_protocol_error <- function() {
  stop("MPI fan-out received an unexpected operation/source/task reply", call. = FALSE)
}

.npRmpi_fanout_accept <- function(tx, message, source, tag, discard = FALSE) {
  if (!.npRmpi_fanout_identity(tx, message, source)) {
    if (discard) return(NULL)
    .npRmpi_fanout_protocol_error()
  }
  kind <- message[["kind"]]
  tasks <- message[["tasks"]]
  assigned <- tx$assigned[[source]]
  if (identical(kind, "closed")) {
    if (!identical(as.integer(tag), tx$control) || length(tasks) ||
        tx$rank[[source]] != "terminal") .npRmpi_fanout_protocol_error()
    tx$rank[[source]] <- "closed"
    return(list(kind = kind, source = source))
  }
  if (identical(kind, "terminal")) {
    if (!identical(as.integer(tag), tx$control) || length(tasks) ||
        tx$rank[[source]] %in% c("initial", "terminal", "closed"))
      .npRmpi_fanout_protocol_error()
    # A valid terminal receipt discharges transport even if the result is
    # incomplete and must now raise a collector error.
    was.stopping <- tx$rank[[source]] == "stopping"
    tx$rank[[source]] <- "terminal"
    .npRmpi_fanout_native("send", tx$comm, as.integer(source), 2L)
    if (!discard && (any(!tx$seen[assigned]) ||
        (identical(tx$scheduler, "dynamic") && !was.stopping)))
      .npRmpi_fanout_protocol_error()
    return(list(kind = kind, source = source))
  }
  if (identical(kind, "ready")) {
    if (!identical(tx$scheduler, "dynamic") ||
        !identical(as.integer(tag), tx$control) ||
        !identical(tasks, assigned) || tx$rank[[source]] != "busy")
      .npRmpi_fanout_protocol_error()
    tx$rank[[source]] <- "ready"
    if (!discard && any(!tx$seen[tasks])) .npRmpi_fanout_protocol_error()
    return(list(kind = kind, source = source))
  }
  if (tx$rank[[source]] != "busy" || !identical(as.integer(tag), tx$tags[[source]]))
    .npRmpi_fanout_protocol_error()
  if (identical(kind, "progress")) {
    boot <- message[["value"]]
    if (!identical(tx$scheduler, "bundle") || length(tasks) ||
        !is.integer(boot) || length(boot) != 1L || is.na(boot) || boot < 1L ||
        tx$progress[[source]] + boot > sum(tx$weights[assigned]))
      .npRmpi_fanout_protocol_error()
    tx$progress[[source]] <- tx$progress[[source]] + boot
    return(list(kind = kind, source = source, boot = boot))
  }
  if (!identical(kind, "result") || !is.integer(tasks) || !length(tasks) ||
      anyNA(tasks) || anyDuplicated(tasks) || any(!(tasks %in% assigned)) ||
      any(tx$seen[tasks]) || !is.list(message[["value"]]) ||
      length(message[["value"]]) != length(tasks)) .npRmpi_fanout_protocol_error()
  tx$seen[tasks] <- TRUE
  if (identical(tx$scheduler, "dynamic")) {
    # Preserve the incumbent result-arrival assignment order even when READY
    # control receipts from different ranks arrive in a different order.
    tx$scheduled <- tx$scheduled + 1L
    tx$next.task[[source]] <- tx$scheduled
  }
  list(kind = kind, source = source, tasks = tasks,
       value = if (discard) NULL else message[["value"]],
       boot = sum(tx$weights[tasks]))
}

.npRmpi_fanout_receive <- function(tx, discard = FALSE, poll = TRUE, sleep = 0.0005,
                                  timeout = 0, started = 0, what = "fan-out",
                                  recovery = NULL) {
  if (poll) {
    while (!isTRUE(mpi.iprobe(mpi.any.source(), mpi.any.tag(), tx$comm))) {
      if (discard) .npRmpi_fanout_cleanup_boundary(tx)
      if (!discard && timeout > 0 &&
          unname(proc.time()[["elapsed"]]) - started > timeout)
        stop(sprintf("MPI %s dispatch timeout waiting on worker results (timeout=%.3fs)",
                     what, timeout), call. = FALSE)
      if (discard) .npRmpi_fanout_recovery_step(recovery)
      Sys.sleep(sleep)
    }
  } else {
    mpi.probe(mpi.any.source(), mpi.any.tag(), tx$comm)
  }
  suspendInterrupts({
    status <- mpi.get.sourcetag()
    bytes <- .npRmpi_recv_raw_probed(status, comm = tx$comm)
    tx$raw.messages <- tx$raw.messages + 1L
    tx$raw.bytes <- tx$raw.bytes + length(bytes)
    # Raw consumption, decoding and control-state commit are interrupt-atomic;
    # callbacks/publication run only after this section. A data decode error
    # still leaves its separate READY/terminal control receipt to be drained.
    decode <- function() {
      message <- unserialize(bytes)
      .npRmpi_fanout_accept(tx, message, as.integer(status[[1L]]),
                            as.integer(status[[2L]]), discard)
    }
    if (discard) tryCatch(decode(), error = function(e) NULL) else decode()
  })
}

.npRmpi_fanout_send <- function(tx, rank, object, tag, tasks = integer()) {
  bytes <- serialize(object, NULL)
  suspendInterrupts({
    mpi.send(bytes, type = 4L, dest = rank, tag = tag, comm = tx$comm)
    tx$assigned[[rank]] <- as.integer(tasks)
    if (length(tasks)) tx$scheduled <- max(tx$scheduled, tasks)
    tx$tags[[rank]] <- as.integer(tag)
    tx$rank[[rank]] <- if (length(tasks)) "busy" else "stopping"
  })
  invisible(NULL)
}

.npRmpi_fanout_stop_ready <- function(tx) {
  if (!identical(tx$scheduler, "dynamic") &&
      !identical(tx$scheduler, "bundle")) return(invisible(NULL))
  for (rank in which(tx$rank %in% c("initial", "ready")))
    .npRmpi_fanout_send(tx, rank, 0L, tx$stop.tag)
  invisible(NULL)
}

.npRmpi_fanout_drain <- function(tx, recovery = NULL) {
  if (identical(tx$phase, "prepared") || identical(tx$phase, "quiescent"))
    return(invisible(TRUE))
  # A failed native activation/collective cannot be repaired by an R handler.
  if (tx$phase %in% c("activating", "uncertain"))
    stop("MPI fan-out activation did not complete; session is not reusable", call. = FALSE)
  tx$phase <- "draining"
  tx$cleanup.started <- unname(proc.time()[["elapsed"]])
  tx$cleanup.budget <- .npRmpi_session_recv_timeout()
  .npRmpi_fanout_notice("Cancelling MPI work at task boundaries; draining replies. Interrupt again to return with this pool quarantined.")
  .npRmpi_fanout_recovery_step(recovery, first = TRUE)
  suspendInterrupts({
    tx$cancelled <- TRUE
    for (rank in which(!(tx$rank %in% c("terminal", "closed"))))
      .npRmpi_fanout_native("send", tx$comm, as.integer(rank), 1L)
  })
  .npRmpi_fanout_stop_ready(tx)
  while (any(tx$rank != "closed")) {
    .npRmpi_fanout_cleanup_boundary(tx)
    .npRmpi_fanout_receive(tx, discard = TRUE, recovery = recovery)
    .npRmpi_fanout_stop_ready(tx)
    .npRmpi_fanout_recovery_step(recovery)
  }
  while (!identical(.npRmpi_fanout_native("poll", tx$comm), 1L)) {
    .npRmpi_fanout_cleanup_boundary(tx)
    Sys.sleep(0.0005)
  }
  suspendInterrupts({
    .npRmpi_fanout_native("finish", tx$comm)
    tx$native <- FALSE
    tx$phase <- "quiescent"
  })
  invisible(TRUE)
}

.npRmpi_fanout_cleanup_boundary <- function(tx) {
  # Serviced under both idle polling and continuous incoming traffic.
  if (tx$cleanup.budget > 0 &&
      unname(proc.time()[["elapsed"]]) - tx$cleanup.started >= tx$cleanup.budget)
    stop(structure(list(message = "MPI cooperative cleanup budget exhausted",
                        call = NULL),
                   class = c("npRmpi_cleanup_pending", "error", "condition")))
  invisible(NULL)
}

.npRmpi_fanout_cleanup_attempt <- function(tx, recovery = NULL) {
  failure <- tryCatch({ .npRmpi_fanout_drain(tx, recovery); NULL },
                      error = identity, interrupt = identity)
  if (!is.null(failure)) {
    tx$phase <- if (inherits(failure, c("interrupt", "npRmpi_cleanup_pending")))
      "quarantined" else "uncertain"
    if (identical(tx$phase, "uncertain")) .npRmpi_lease_poison()
    .npRmpi_fanout_notice("MPI cleanup is incomplete; this pool is quarantined, not closed. No partial result was returned. An explicit npRmpi.quit(force = TRUE) may resume cooperative cleanup; native transport failures cannot be retried safely.")
  }
  failure
}

.npRmpi_fanout_forget <- function(tx) {
  if (tx$phase %in% c("prepared", "quiescent") &&
      exists(tx$key, envir = .npRmpi_fanout_state, inherits = FALSE))
    rm(list = tx$key, envir = .npRmpi_fanout_state)
  invisible(NULL)
}

.npRmpi_fanout_run <- function(tx, code, recovery = .npRmpi_fanout_recovery_owner(tx$comm)) {
  force(recovery)
  assign(tx$key, tx, envir = .npRmpi_fanout_state)
  tx$owner.active <- TRUE
  on.exit({
    tx$owner.active <- FALSE
    if (!(tx$phase %in% c("prepared", "quiescent", "uncertain", "quarantined")))
      .npRmpi_fanout_cleanup_attempt(tx, recovery)
    .npRmpi_fanout_forget(tx)
  }, add = TRUE)
  failure <- new.env(parent = emptyenv())
  failure$condition <- NULL
  value <- tryCatch(force(code), error = function(e) {
    failure$condition <- e
    NULL
  }, interrupt = function(e) {
    failure$condition <- e
    NULL
  })
  if (!is.null(failure$condition)) {
    .npRmpi_fanout_cleanup_attempt(tx, recovery)
    suspendInterrupts({
      .npRmpi_fanout_forget(tx)
      .npRmpi_fanout_rethrow(failure$condition)
    })
  }
  if (any(tx$rank != "closed"))
    stop("MPI fan-out returned before worker quiescence", call. = FALSE)
  while (!identical(.npRmpi_fanout_native("poll", tx$comm), 1L)) Sys.sleep(0.0005)
  suspendInterrupts({
    .npRmpi_fanout_native("finish", tx$comm)
    tx$native <- FALSE
    tx$phase <- "quiescent"
  })
  value
}

.npRmpi_fanout_quiesce <- function(comm = 1L, resume = FALSE) {
  tx <- .npRmpi_fanout_retained(comm)
  if (is.null(tx)) return(invisible(TRUE))
  # Ordinary reentrant lifecycle calls must unwind the active owner rather
  # than close its communicator and resume its R frame. Finalizers likewise
  # must not start a blocking cleanup behind the caller's back.
  if (isTRUE(tx$owner.active))
    stop("MPI lifecycle cannot begin inside an active fan-out", call. = FALSE)
  if (!(tx$phase %in% c("prepared", "quiescent"))) {
    if (!resume || isTRUE(.npRmpi_exit_finalizer_state$running))
      stop("MPI pool is quarantined; implicit cleanup was not attempted. Use npRmpi.quit(force = TRUE) explicitly.", call. = FALSE)
    failure <- .npRmpi_fanout_cleanup_attempt(tx)
    if (!is.null(failure)) {
      if (inherits(failure, "interrupt")) {
        signalCondition(failure)
        invokeRestart("abort")
      }
      stop(failure)
    }
  }
  .npRmpi_fanout_forget(tx)
  invisible(TRUE)
}

.npRmpi_fanout_worker_terminal <- function(header, terminal, closed) {
  suspendInterrupts({
    mpi.send(terminal, type = 4L, dest = 0L, tag = header$control, comm = header$comm)
    while (.npRmpi_fanout_native("poll", header$comm) < 2L)
      tryCatch(Sys.sleep(0.0005), interrupt = function(e) NULL)
    .npRmpi_fanout_native("finish", header$comm)
    mpi.send(closed, type = 4L, dest = 0L, tag = header$control, comm = header$comm)
  })
  invisible(NULL)
}

.npRmpi_fanout_worker_owner <- function(header, code) {
  # Allocate retirement payloads before entering task code. The outer owner
  # covers ordinary failures in serialization/progress as well as FUN.
  terminal <- serialize(.npRmpi_fanout_envelope(header, "terminal"), NULL)
  closed <- serialize(.npRmpi_fanout_envelope(header, "closed"), NULL)
  suspendInterrupts(.npRmpi_fanout_native("begin", header$comm, NULL))
  on.exit(.npRmpi_fanout_worker_terminal(header, terminal, closed), add = TRUE)
  tryCatch(force(code), error = function(e) invisible(NULL),
           interrupt = function(e) invisible(NULL))
  invisible(NULL)
}

.npRmpi_fanout_worker_cancelled <- function(header) {
  .npRmpi_fanout_native("poll", header$comm) %% 2L == 1L
}

.npRmpi_fanout_worker_result <- function(header, tasks, parts, tag) {
  mpi.send.Robj(.npRmpi_fanout_envelope(header, "result", tasks, parts),
                 0L, tag, header$comm)
  if (identical(header$scheduler, "dynamic"))
    mpi.send.Robj(.npRmpi_fanout_envelope(header, "ready", tasks),
                   0L, header$control, header$comm)
  invisible(NULL)
}

.npRmpi_fanout_worker_call <- function(FUN, args) {
  tryCatch(do.call(FUN, args), error = function(e)
    structure(conditionMessage(e), class = "try-error", condition = e))
}

.npRmpi_fanout_worker_header <- function(shared, scheduler, comm) {
  header <- if (is.list(shared)) shared[["transaction"]] else NULL
  scalar.text <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  if (!is.list(shared) || anyDuplicated(names(shared)) ||
      !is.list(header) || anyDuplicated(names(header)) ||
      !all(c("session", "operation", "scheduler", "control", "comm") %in% names(header)) ||
      !scalar.text(header[["session"]]) || !scalar.text(header[["operation"]]) ||
      !identical(header[["scheduler"]], scheduler) ||
      !identical(header[["control"]], .npRmpi_protocol_tag("fanout_control")) ||
      !is.numeric(header[["comm"]]) || length(header[["comm"]]) != 1L ||
      is.na(header[["comm"]]) || header[["comm"]] != comm ||
      !is.function(shared[["FUN"]]) || !is.list(shared[["dot.arg"]]))
    stop("MPI fan-out worker received a missing or malformed transaction header", call. = FALSE)
  invisible(header)
}

.npRmpi_fanout_worker_apply <- function(tmpfunarg, n, tag, comm) {
  header <- tmpfunarg[["transaction"]]
  .npRmpi_fanout_worker_owner(header, {
  x <- mpi.scatter.Robj(root = 0L, comm = comm)
  rank <- mpi.comm.rank(comm)
  if (rank <= n && !.npRmpi_fanout_worker_cancelled(header)) {
    value <- .npRmpi_fanout_worker_call(tmpfunarg$FUN, c(list(x), tmpfunarg$dot.arg))
    .npRmpi_fanout_worker_result(header, as.integer(rank), list(value), tag)
  }
  })
}

.npRmpi_fanout_worker_dynamic <- function(tmpfunarg, n, comm) {
  header <- tmpfunarg[["transaction"]]
  .npRmpi_fanout_worker_owner(header, {
  repeat {
    request <- mpi.recv.Robj(0L, mpi.any.tag(), comm)
    tag <- mpi.get.sourcetag()[[2L]]
    # Consume the incumbent request/stop handshake before cancellation. A
    # task already sent by the master must not remain queued for a later call.
    if (tag > n || .npRmpi_fanout_worker_cancelled(header)) break
    args <- if (is.list(request)) request$data.arg else NULL
    value <- if (is.null(args)) {
      structure("mpi.applyLB worker received malformed task payload", class = "try-error")
    } else .npRmpi_fanout_worker_call(tmpfunarg$FUN, c(args, tmpfunarg$dot.arg))
    .npRmpi_fanout_worker_result(header, as.integer(tag), list(value), tag)
  }
  })
}

.npRmpi_fanout_worker_bundle <- function(tmpfunarg, n, comm) {
  header <- tmpfunarg[["transaction"]]
  .npRmpi_fanout_worker_owner(header, {
  request <- mpi.recv.Robj(0L, mpi.any.tag(), comm)
  tag <- mpi.get.sourcetag()[[2L]]
  if (tag <= n) {
  tasks <- request$task_indices
  parts <- vector("list", length(tasks))
  boot <- 0L
  completed <- 0L
  for (i in seq_along(tasks)) {
    if (.npRmpi_fanout_worker_cancelled(header)) break
    value <- .npRmpi_bootstrap_task_evaluate(
      .npRmpi_fanout_worker_call(tmpfunarg$FUN,
        c(list(request$tasks[[i]]), tmpfunarg$dot.arg)),
      weight = request$tasks[[i]]$bsz,
      report = if (isTRUE(tmpfunarg$progress.enabled) &&
                   isTRUE(request$tasks[[i]]$report.internal)) function(delta) {
        mpi.send.Robj(.npRmpi_fanout_envelope(header, "progress", value = delta),
                       0L, tag, comm)
      })
    if (isTRUE(tmpfunarg$stream.results)) {
      .npRmpi_fanout_worker_result(header, tasks[[i]], list(value), tag)
    } else {
      parts[i] <- list(value)
      if (isTRUE(tmpfunarg$progress.enabled) &&
          !isTRUE(request$tasks[[i]]$report.internal)) {
        boot <- boot + request$tasks[[i]]$bsz
        if (i %% tmpfunarg$progress.stride == 0L || i == length(tasks)) {
          mpi.send.Robj(.npRmpi_fanout_envelope(header, "progress",
                         value = as.integer(boot)), 0L, tag, comm)
          boot <- 0L
        }
      }
    }
    completed <- i
  }
  if (!isTRUE(tmpfunarg$stream.results) && completed == length(tasks))
    .npRmpi_fanout_worker_result(header, tasks, parts, tag)
  }
  })
}

.npRmpi_fanout_apply <- function(X, FUN, dot.arg, comm, dynamic = FALSE,
                                poll = FALSE, sleep = 0.01) {
  n <- length(X)
  workers <- mpi.comm.size(comm) - 1L
  if (!dynamic && n > workers) stop("data length must be at most total slave size")
  if (!is.function(FUN)) stop("FUN is not a function")
  force(dot.arg)
  # Retain mpi.apply/iapply's incumbent RNG draw; transaction identity itself
  # is deterministic and does not consume statistical RNG.
  tag <- if (dynamic) NA_integer_ else floor(runif(1, 1, 1000))
  tx <- .npRmpi_fanout_new(comm, workers, n, if (dynamic) "dynamic" else "scatter")
  tx$tags <- rep.int(as.integer(tag), workers)
  tx$stop.tag <- as.integer(n + 1L)
  header <- .npRmpi_fanout_header(tx)
  shared <- serialize(list(FUN = FUN, dot.arg = dot.arg, transaction = header), NULL)
  scatter <- if (!dynamic) {
    # Preserve the incumbent c() dispatch before as.list() for classed X.
    if (n < workers) X <- c(X, as.list(integer(workers - n)))
    .npRmpi_scatter_prepare(c(list("master"), as.list(X)), comm)
  } else NULL
  out <- as.list(integer(n))
  .npRmpi_fanout_run(tx, {
    suspendInterrupts({
      tx$phase <- "activating"
      if (dynamic) mpi.bcast.cmd(.mpi.worker.applyLB, n = n, comm = comm)
      else mpi.bcast.cmd(.mpi.worker.apply, n = n, tag = tag, comm = comm)
      .npRmpi_bcast_prepared(shared, rank = 0L, comm = comm)
      .npRmpi_fanout_native("begin", comm, tx)
      tx$native <- TRUE
      if (!dynamic) {
        .npRmpi_scatter_prepared(scatter, root = 0L, comm = comm)
        for (rank in seq_len(workers)) {
          tx$assigned[[rank]] <- if (rank <= n) as.integer(rank) else integer()
          tx$rank[[rank]] <- "busy"
        }
      }
      tx$phase <- "active"
    })
    sent <- 0L
    if (dynamic) {
      for (rank in seq_len(workers)) {
        sent <- sent + 1L
        .npRmpi_fanout_send(tx, rank, list(data.arg = list(X[[sent]])), sent, sent)
      }
    }
    while (any(tx$rank != "closed")) {
      message <- .npRmpi_fanout_receive(tx, poll = poll, sleep = sleep)
      if (identical(message$kind, "result"))
        out[[message$tasks]] <- message$value[[1L]]
      if (identical(message$kind, "ready")) {
        next.task <- tx$next.task[[message$source]]
        if (next.task <= n)
          .npRmpi_fanout_send(tx, message$source, list(data.arg = list(X[[next.task]])), next.task, next.task)
        else .npRmpi_fanout_send(tx, message$source, 0L, tx$stop.tag)
      }
    }
    out
  })
}

.npRmpi_fanout_bootstrap <- function(tasks, worker, dot.arg, workers, comm,
                                    master.local, ncol.out, progress.enabled,
                                    progress.step, what) {
  n <- length(tasks)
  weights <- vapply(tasks, function(task) as.integer(task$bsz), integer(1L))
  tx <- .npRmpi_fanout_new(comm, workers, n,
    if (master.local) "bundle" else "dynamic", weights)
  tx$tags <- integer(workers)
  tx$stop.tag <- as.integer(if (master.local) workers + 1L else n + 1L)
  local.idx <- if (master.local) which((seq_len(n) - 1L) %% (workers + 1L) == 0L) else integer()
  worker.idx <- if (master.local) lapply(seq_len(workers), function(rank)
    which((seq_len(n) - 1L) %% (workers + 1L) == rank)) else NULL
  stream <- master.local && .npRmpi_bootstrap_stream_bundle_results(worker.idx, tasks, ncol.out)
  # Internal-reporting tasks can retain completion beacons even when their
  # result payloads are streamed. Other consumers keep the incumbent policy.
  beacons <- master.local && progress.enabled && (!stream ||
    all(vapply(tasks, function(task) isTRUE(task$report.internal), logical(1L))))
  stride <- if (beacons) .npRmpi_bootstrap_bundle_progress_stride(worker.idx) else 0L
  shared <- serialize(list(FUN = worker, dot.arg = dot.arg,
    transaction = .npRmpi_fanout_header(tx), stream.results = stream,
    progress.enabled = beacons, progress.stride = stride), NULL)
  state <- new.env(parent = emptyenv())
  state$out <- vector("list", n)
  state$done <- 0L
  timeout <- .npRmpi_bootstrap_dispatch_timeout_sec()
  started <- unname(proc.time()[["elapsed"]])
  receive <- function() {
    message <- .npRmpi_fanout_receive(tx, timeout = timeout, started = started, what = what)
    if (is.null(message)) return(invisible(NULL))
    if (identical(message$kind, "progress")) state$done <- state$done + message$boot
    if (identical(message$kind, "result")) {
      state$out[message$tasks] <- message$value
      if (!beacons) state$done <- state$done + message$boot
      .npRmpi_bootstrap_transport_trace(what, "fanout.recv",
        list(src = message$source, tag = tx$tags[[message$source]]))
    }
    if (message$kind %in% c("progress", "result")) progress.step(state$done)
    if (identical(message$kind, "ready")) {
      next.task <- tx$next.task[[message$source]]
      if (next.task <= n) {
        .npRmpi_fanout_send(tx, message$source,
          list(data.arg = list(tasks[[next.task]])), next.task, next.task)
        .npRmpi_bootstrap_transport_trace(what, "fanout.send.next",
          list(dest = message$source, tag = next.task, task_idx = next.task))
      } else {
        .npRmpi_fanout_send(tx, message$source, 0L, tx$stop.tag)
        .npRmpi_bootstrap_transport_trace(what, "fanout.send.stop",
          list(dest = message$source, tag = tx$stop.tag))
      }
    }
    invisible(NULL)
  }
  drain.available <- function() {
    if (!identical(tx$phase, "active") || isTRUE(state$polling))
      return(invisible(NULL))
    state$polling <- TRUE
    on.exit(state$polling <- FALSE, add = TRUE)
    # Rank-local native computation temporarily maps comm[1] to MPI_COMM_SELF.
    # Receipt polling must use the saved pool, then restore the numerical
    # owner's local mode even if collection or rendering signals a condition.
    old.mode <- .Call("C_np_set_local_regression_mode", FALSE, PACKAGE = "npRmpi")
    on.exit(.Call("C_np_set_local_regression_mode", old.mode, PACKAGE = "npRmpi"),
            add = TRUE)
    while (any(tx$rank != "closed") &&
           isTRUE(mpi.iprobe(mpi.any.source(), mpi.any.tag(), comm))) receive()
    invisible(NULL)
  }
  if (master.local && progress.enabled &&
      any(vapply(tasks, function(task) isTRUE(task$report.internal), logical(1L)))) {
    previous.forward <- .np_progress_runtime$fit_forward
    on.exit(.np_progress_runtime$fit_forward <- previous.forward, add = TRUE)
    .np_progress_runtime$fit_forward <- function() {
      # Existing native heartbeats service completed worker receipts while
      # the master's numerical tile is still active. Otherwise blocking sends
      # can stall workers until the master's next whole-tile boundary.
      # The native progress bridge deliberately tolerates cosmetic R errors.
      # Retain transport/interrupt conditions here and rethrow at the next
      # ordinary R-owned task boundary, where the transaction can unwind.
      if (is.null(state$pump.failure)) {
        state$pump.failure <- tryCatch({ drain.available(); NULL },
          error = identity, interrupt = identity)
      }
      if (is.function(previous.forward)) previous.forward()
      else progress.step(state$done)
      invisible(NULL)
    }
  }
  .npRmpi_fanout_run(tx, {
    suspendInterrupts({
      tx$phase <- "activating"
      if (master.local) mpi.bcast.cmd(.npRmpi_bootstrap_worker_bundle, n = workers, comm = comm)
      else mpi.bcast.cmd(.mpi.worker.applyLB, n = n, comm = comm)
      .npRmpi_bcast_prepared(shared, rank = 0L, comm = comm)
      .npRmpi_fanout_native("begin", comm, tx)
      tx$native <- TRUE
      tx$phase <- "active"
    })
    .npRmpi_bootstrap_transport_trace(what, "fanout.master_assist.start",
      list(n_remote = if (master.local) sum(lengths(worker.idx) > 0L) else n,
           slave_num = workers, local_n = length(local.idx),
           scheduler = if (master.local) "static_bundle" else "dynamic_worker_only",
           stream_results = stream, progress_beacons = beacons, progress_stride = stride))
    for (rank in seq_len(workers)) {
      ids <- if (master.local) as.integer(worker.idx[[rank]]) else if (rank <= n) as.integer(rank) else integer()
      if (length(ids)) {
        payload <- if (master.local) list(task_indices = ids, tasks = tasks[ids]) else list(data.arg = list(tasks[[rank]]))
        .npRmpi_fanout_send(tx, rank, payload, rank, ids)
        .npRmpi_bootstrap_transport_trace(what, "fanout.send.initial",
          list(dest = rank, tag = rank, task_count = length(ids),
               task_first = ids[[1L]], task_last = ids[[length(ids)]]))
      } else {
        .npRmpi_fanout_send(tx, rank, 0L, tx$stop.tag)
        .npRmpi_bootstrap_transport_trace(what, "fanout.send.stop.initial",
          list(dest = rank, tag = tx$stop.tag))
      }
    }
    for (task.local.idx in local.idx) {
      .npRmpi_bootstrap_transport_trace(what, "fanout.master_local_chunk.start",
        list(task_idx = task.local.idx))
      internal.progress <- progress.enabled &&
        isTRUE(tasks[[task.local.idx]]$report.internal)
      state$out[task.local.idx] <- list(.npRmpi_bootstrap_task_evaluate(
        .npRmpi_fanout_worker_call(worker,
          c(list(tasks[[task.local.idx]]), dot.arg)),
        weight = weights[[task.local.idx]],
        report = if (internal.progress) function(delta) {
          if (!is.null(state$pump.failure))
            .npRmpi_fanout_rethrow(state$pump.failure)
          state$done <- state$done + delta
          drain.available()
          progress.step(state$done)
        }))
      if (!is.null(state$pump.failure))
        .npRmpi_fanout_rethrow(state$pump.failure)
      if (!internal.progress)
        state$done <- state$done + weights[[task.local.idx]]
      tx$seen[[task.local.idx]] <- TRUE
      progress.step(state$done)
      .npRmpi_bootstrap_transport_trace(what, "fanout.master_local_chunk.done",
        list(task_idx = task.local.idx, bsz = weights[[task.local.idx]]))
      while (any(tx$rank != "closed") && isTRUE(mpi.iprobe(mpi.any.source(), mpi.any.tag(), comm))) receive()
    }
    while (any(tx$rank != "closed")) receive()
    .npRmpi_bootstrap_transport_trace(what, "fanout.master_assist.done",
      list(done = sum(tx$rank == "closed"), n_remote = workers, local_done = length(local.idx)))
    state$out
  })
}
