.npRmpi_session_apply_options <- function(autodispatch = NULL,
                                          np.messages = NULL,
                                          autodispatch.verify.options = NULL,
                                          autodispatch.option.sync = NULL) {
  if (!is.null(autodispatch))
    options(npRmpi.autodispatch = isTRUE(autodispatch))
  if (!is.null(np.messages))
    options(np.messages = isTRUE(np.messages))
  if (!is.null(autodispatch.verify.options))
    options(npRmpi.autodispatch.verify.options = isTRUE(autodispatch.verify.options))
  if (!is.null(autodispatch.option.sync))
    options(npRmpi.autodispatch.option.sync = as.character(autodispatch.option.sync)[1L])
  invisible(TRUE)
}

.npRmpi_session_reset_spmd_state <- function() {
  # SPMD sequencing must restart for each fresh session world.
  options(npRmpi.spmd.seq_id = 0L)
  options(npRmpi.autodispatch.remote.counter = 0L)
  invisible(TRUE)
}

.npRmpi_session_reset_spmd_state_allranks <- function(comm = 1L,
                                                      strict = TRUE,
                                                      where = "npRmpi session state reset") {
  cmd <- quote({
    fn <- get(".npRmpi_session_reset_spmd_state", envir = asNamespace("npRmpi"), inherits = FALSE)
    fn()
    invisible(TRUE)
  })
  out <- try(
    .npRmpi_bcast_cmd_expr(expr = cmd, comm = comm, caller.execute = TRUE),
    silent = TRUE
  )
  if (inherits(out, "try-error")) {
    msg <- paste(as.character(out), collapse = " ")
    if (isTRUE(strict))
      stop(sprintf("%s failed across ranks: %s", where, msg), call. = FALSE)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

.npRmpi_safe_int <- function(expr) {
  tryCatch(as.integer(expr), error = function(e) NA_integer_)
}

.npRmpi_safe <- function(expr, fallback = NULL) {
  tryCatch(expr, error = function(e) fallback)
}

.npRmpi_session_ensure_comm <- function(comm = 1L) {
  size <- .npRmpi_safe_int(mpi.comm.size(comm))
  if (!is.na(size) && as.integer(size) >= 1L)
    return(invisible(as.integer(size)))
  invisible(mpi.comm.dup(0, comm))
}

.npRmpi_session_has_active_pool <- function(comm = 1L) {
  size <- .npRmpi_safe_int(mpi.comm.size(comm))
  rank <- .npRmpi_safe_int(mpi.comm.rank(comm))
  if (is.na(size) || is.na(rank))
    return(FALSE)
  as.integer(size) >= 2L
}

.npRmpi_autodispatch_prime_options <- function() {
  keys <- tryCatch(.npRmpi_autodispatch_option_keys(), error = function(e) character(0))
  if (!length(keys))
    return(invisible(FALSE))

  snap <- tryCatch(.npRmpi_autodispatch_option_snapshot(keys), error = function(e) NULL)
  if (!is.list(snap))
    return(invisible(FALSE))

  options(npRmpi.autodispatch.option.snapshot = snap)

  invisible(TRUE)
}

.npRmpi_session_recv_timeout <- function() {
  opt <- getOption("npRmpi.session.recv.timeout", NULL)
  if (!is.null(opt)) {
    if (is.numeric(opt) && length(opt) == 1L && is.finite(opt) && opt > 0)
      return(as.numeric(opt))
    return(0)
  }

  envv <- suppressWarnings(as.numeric(Sys.getenv("NP_RMPI_SESSION_RECV_TIMEOUT_SEC", "0")))
  if (!is.finite(envv) || envv <= 0)
    return(0)
  as.numeric(envv)
}

.npRmpi_attach_state_reset <- function() {
  options(npRmpi.attach.close.state = "closed")
  options(npRmpi.attach.session.id = NA_integer_)
  invisible(TRUE)
}

.npRmpi_attach_next_session_id <- function() {
  sid <- suppressWarnings(as.integer(getOption("npRmpi.attach.session.counter", 0L)))
  if (!is.finite(sid) || sid < 0L)
    sid <- 0L
  sid <- as.integer(sid + 1L)
  if (!is.finite(sid) || sid < 1L)
    sid <- 1L
  options(npRmpi.attach.session.counter = sid)
  sid
}

.npRmpi_attach_close_ack_tag <- function(rank) {
  rank <- suppressWarnings(as.integer(rank)[1L])
  if (!is.finite(rank) || rank < 1L)
    rank <- 1L
  as.integer(62000L + rank)
}

.npRmpi_attach_close_release_tag <- function(rank) {
  rank <- suppressWarnings(as.integer(rank)[1L])
  if (!is.finite(rank) || rank < 1L)
    rank <- 1L
  as.integer(62100L + rank)
}

.npRmpi_attach_close_ack_timeout <- function() {
  opt <- getOption("npRmpi.attach.close.ack.timeout", NULL)
  if (!is.null(opt)) {
    if (is.numeric(opt) && length(opt) == 1L && is.finite(opt) && opt > 0)
      return(as.numeric(opt))
    return(5)
  }
  envv <- suppressWarnings(as.numeric(Sys.getenv("NP_RMPI_ATTACH_CLOSE_ACK_TIMEOUT_SEC", "5")))
  if (!is.finite(envv) || envv <= 0)
    return(5)
  as.numeric(envv)
}

.npRmpi_attach_send_close_ack <- function(comm = 1L,
                                          session_id = getOption("npRmpi.attach.session.id", NA_integer_)) {
  rank <- .npRmpi_safe_int(mpi.comm.rank(comm))
  if (is.na(rank) || rank < 1L)
    return(invisible(FALSE))
  payload <- list(
    kind = "attach_close_ack",
    rank = as.integer(rank),
    session_id = suppressWarnings(as.integer(session_id)[1L]),
    timestamp = as.numeric(Sys.time())
  )
  sent <- try(
    mpi.send.Robj(payload, dest = 0L, tag = .npRmpi_attach_close_ack_tag(rank), comm = comm),
    silent = TRUE
  )
  invisible(!inherits(sent, "try-error"))
}

.npRmpi_attach_send_close_release <- function(comm = 1L,
                                              ranks,
                                              session_id = getOption("npRmpi.attach.session.id", NA_integer_),
                                              barrier = TRUE) {
  ranks <- as.integer(ranks)
  ranks <- sort(unique(ranks[is.finite(ranks) & ranks >= 1L]))
  if (!length(ranks))
    return(list(ok = TRUE, sent = integer(0), failed = integer(0)))

  sid <- suppressWarnings(as.integer(session_id)[1L])
  payload <- list(
    kind = "attach_close_release",
    session_id = sid,
    barrier = isTRUE(barrier),
    timestamp = as.numeric(Sys.time())
  )

  sent <- integer(0)
  failed <- integer(0)
  for (rk in ranks) {
    ok <- try(
      mpi.send.Robj(payload, dest = rk, tag = .npRmpi_attach_close_release_tag(rk), comm = comm),
      silent = TRUE
    )
    if (inherits(ok, "try-error"))
      failed <- c(failed, rk)
    else
      sent <- c(sent, rk)
  }

  list(
    ok = !length(failed),
    sent = sort(unique(as.integer(sent))),
    failed = sort(unique(as.integer(failed)))
  )
}

.npRmpi_attach_wait_close_release <- function(comm = 1L,
                                              session_id = getOption("npRmpi.attach.session.id", NA_integer_),
                                              timeout_sec = .npRmpi_attach_close_ack_timeout(),
                                              poll_sleep = 0.01) {
  rank <- .npRmpi_safe_int(mpi.comm.rank(comm))
  if (is.na(rank) || rank < 1L)
    return(list(ok = FALSE, barrier = FALSE, payload = NULL))

  timeout_sec <- if (is.numeric(timeout_sec) && length(timeout_sec) == 1L && is.finite(timeout_sec) && timeout_sec > 0)
    as.numeric(timeout_sec) else 5
  poll_sleep <- if (is.numeric(poll_sleep) && length(poll_sleep) == 1L && is.finite(poll_sleep) && poll_sleep > 0)
    as.numeric(poll_sleep) else 0.01
  deadline <- as.numeric(Sys.time()) + timeout_sec
  tag <- .npRmpi_attach_close_release_tag(rank)
  sid <- suppressWarnings(as.integer(session_id)[1L])

  repeat {
    if (isTRUE(.npRmpi_safe(mpi.iprobe(source = 0L, tag = tag, comm = comm), fallback = FALSE))) {
      msg <- .npRmpi_safe(mpi.recv.Robj(source = 0L, tag = tag, comm = comm), fallback = NULL)
      if (is.list(msg) &&
          identical(as.character(msg$kind)[1L], "attach_close_release")) {
        msg.sid <- suppressWarnings(as.integer(msg$session_id)[1L])
        if (is.na(sid) || is.na(msg.sid) || identical(msg.sid, sid)) {
          return(list(
            ok = TRUE,
            barrier = isTRUE(msg$barrier),
            payload = msg
          ))
        }
      }
      return(list(ok = FALSE, barrier = FALSE, payload = msg))
    }
    if (as.numeric(Sys.time()) >= deadline)
      break
    Sys.sleep(poll_sleep)
  }

  list(ok = FALSE, barrier = FALSE, payload = NULL)
}

.npRmpi_attach_collect_close_acks <- function(comm = 1L,
                                              session_id = getOption("npRmpi.attach.session.id", NA_integer_),
                                              expected = NULL,
                                              timeout_sec = .npRmpi_attach_close_ack_timeout(),
                                              poll_sleep = 0.01) {
  if (is.null(expected)) {
    size <- .npRmpi_safe_int(mpi.comm.size(comm))
    size <- if (is.na(size)) 0L else as.integer(size)
    expected <- if (size >= 2L) seq_len(size - 1L) else integer(0)
  }
  expected <- as.integer(expected)
  expected <- expected[is.finite(expected) & expected >= 1L]
  expected <- sort(unique(expected))

  if (!length(expected)) {
    return(list(ok = TRUE, acked = integer(0), missing = integer(0), payloads = list()))
  }

  sid <- suppressWarnings(as.integer(session_id)[1L])
  payloads <- vector("list", length(expected))
  names(payloads) <- as.character(expected)
  acked <- integer(0)
  timeout_sec <- if (is.numeric(timeout_sec) && length(timeout_sec) == 1L && is.finite(timeout_sec) && timeout_sec > 0)
    as.numeric(timeout_sec) else 5
  poll_sleep <- if (is.numeric(poll_sleep) && length(poll_sleep) == 1L && is.finite(poll_sleep) && poll_sleep > 0)
    as.numeric(poll_sleep) else 0.01
  deadline <- as.numeric(Sys.time()) + timeout_sec

  repeat {
    pending <- expected[!(expected %in% acked)]
    if (!length(pending))
      break
    if (as.numeric(Sys.time()) >= deadline)
      break
    for (rk in pending) {
      tag <- .npRmpi_attach_close_ack_tag(rk)
      if (!isTRUE(.npRmpi_safe(mpi.iprobe(source = rk, tag = tag, comm = comm), fallback = FALSE)))
        next
      msg <- .npRmpi_safe(mpi.recv.Robj(source = rk, tag = tag, comm = comm), fallback = NULL)
      payloads[[as.character(rk)]] <- msg
      if (!is.list(msg))
        next
      if (!identical(as.character(msg$kind)[1L], "attach_close_ack"))
        next
      if (!identical(as.integer(msg$rank)[1L], as.integer(rk)))
        next
      msg.sid <- suppressWarnings(as.integer(msg$session_id)[1L])
      if (!is.na(sid) && !is.na(msg.sid) && !identical(msg.sid, sid))
        next
      acked <- c(acked, rk)
    }
    if (length(expected[!(expected %in% acked)]) > 0L)
      Sys.sleep(poll_sleep)
  }

  acked <- sort(unique(as.integer(acked)))
  missing <- expected[!(expected %in% acked)]
  list(ok = !length(missing), acked = acked, missing = missing, payloads = payloads)
}

.npRmpi_worker_normalize_message <- function(msg) {
  if (!is.call(msg) || length(msg) < 4L)
    return(msg)
  if (!identical(msg[[1L]], as.name("do.call")))
    return(msg)
  if (!is.environment(msg[[4L]]))
    return(msg)
  as.call(list(
    as.name("do.call"),
    msg[[2L]],
    msg[[3L]],
    quote = FALSE,
    envir = msg[[4L]]
  ))
}

.npRmpi_worker_is_break_message <- function(msg) {
  if (is.character(msg) && length(msg) >= 1L && identical(msg[[1L]], "kaerb"))
    return(TRUE)
  if (!is.call(msg) || length(msg) < 2L)
    return(FALSE)
  hd <- msg[[1L]]
  payload <- msg[[2L]]
  if (!is.character(payload) || length(payload) < 1L || !identical(payload[[1L]], "kaerb"))
    return(FALSE)
  if (is.symbol(hd) && identical(as.character(hd), ".npRmpi_with_manual_bcast_context"))
    return(TRUE)
  if (is.function(hd))
    return(TRUE)
  FALSE
}

.npRmpi_worker_loop <- function(comm = 1L,
                                nonblock = TRUE,
                                sleep = 0.1,
                                recv.timeout = 0,
                                loop.label = "worker",
                                timeout.remediation = "") {
  repeat {
    if (recv.timeout > 0)
      base::setTimeLimit(elapsed = recv.timeout, transient = TRUE)
    msg <- try(mpi.bcast.cmd(rank = 0, comm = comm, nonblock = nonblock, sleep = sleep),
               silent = TRUE)
    if (recv.timeout > 0)
      base::setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    if (inherits(msg, "try-error")) {
      msg.text <- as.character(msg)
      if (any(grepl("reached elapsed time limit", msg.text, fixed = TRUE))) {
        rem <- if (is.character(timeout.remediation) && nzchar(timeout.remediation))
          paste("Remediation:", timeout.remediation) else ""
        stop(
          paste(
            sprintf("npRmpi %s receive timeout waiting for master command.", loop.label),
            "Possible blocked MPI route or transport deadlock.",
            rem
          ),
          call. = FALSE
        )
      }
      stop(msg.text, call. = FALSE)
    }
    if (.npRmpi_worker_is_break_message(msg))
      break

    msg <- .npRmpi_worker_normalize_message(msg)
    res <- try(.npRmpi_eval_scmd(msg, envir = .GlobalEnv), silent = TRUE)
    if (inherits(res, "try-error")) {
      cmd <- paste(utils::capture.output(print(msg)), collapse = " ")
      err <- as.character(res)
      rank <- .npRmpi_safe_int(mpi.comm.rank(comm))
      base::cat(
        sprintf("\n[%s rank %d] CMD: %s\n[%s rank %d] ERROR: %s\n",
                loop.label, rank, cmd, loop.label, rank, paste(err, collapse = " ")),
        file = stderr()
      )
      flush(stderr())
      stop(err, call. = FALSE)
    }
  }
  invisible(TRUE)
}

.npRmpi_session_attach_worker_loop <- function(comm = 1L,
                                               nonblock = TRUE,
                                               sleep = 0.1,
                                               session_id = getOption("npRmpi.attach.session.id", NA_integer_)) {
  options(echo = FALSE)
  recv.timeout <- .npRmpi_session_recv_timeout()
  .npRmpi_worker_loop(
    comm = comm,
    nonblock = nonblock,
    sleep = sleep,
    recv.timeout = recv.timeout,
    loop.label = "attach slave",
    timeout.remediation = paste(
      "verify attach launch wiring (mpiexec -n <master+slaves>,",
      "clear R_PROFILE_USER/R_PROFILE for attach), confirm FI_TCP_IFACE,",
      "and relaunch with a fresh MPI world."
    )
  )
  .npRmpi_attach_send_close_ack(comm = comm, session_id = session_id)
  release <- .npRmpi_attach_wait_close_release(
    comm = comm,
    session_id = session_id,
    timeout_sec = .npRmpi_attach_close_ack_timeout()
  )
  if (isTRUE(release$ok) && isTRUE(release$barrier))
    .npRmpi_safe(mpi.barrier(0), fallback = NULL)
  .npRmpi_safe(if (comm != 0L) mpi.comm.free(comm), fallback = NULL)
  mpi.quit()
  invisible(FALSE)
}

npRmpi.init <- function(...,
                         nslaves = 1,
                         comm = 1,
                         mode = c("auto", "spawn", "attach"),
                         autodispatch = TRUE,
                         autodispatch.verify.options = FALSE,
                         autodispatch.option.sync = c("onchange", "always", "never"),
                         np.messages = NULL,
                         nonblock = TRUE,
                         sleep = 0.1,
                         quiet = FALSE) {
  .npRmpi_abort_if_rmpi_attached(where = "npRmpi.init()")

  nslaves <- npValidateNonNegativeInteger(nslaves, "nslaves")
  if (nslaves < 1L) {
    stop(
      "'nslaves' must be >= 1 for npRmpi; use package 'np' for serial workflows",
      call. = FALSE
    )
  }
  comm <- npValidatePositiveInteger(comm, "comm")
  autodispatch <- npValidateScalarLogical(autodispatch, "autodispatch")
  autodispatch.verify.options <- npValidateScalarLogical(autodispatch.verify.options, "autodispatch.verify.options")
  if (!is.null(np.messages))
    np.messages <- npValidateScalarLogical(np.messages, "np.messages")
  nonblock <- npValidateScalarLogical(nonblock, "nonblock")
  quiet <- npValidateScalarLogical(quiet, "quiet")
  sleep <- npValidatePositiveFiniteNumeric(sleep, "sleep")

  if (!isTRUE(getOption("npRmpi.mpi.initialized", FALSE))) {
    stop(
      paste(
        "MPI runtime is not initialized in this R session.",
        "Restart R and load npRmpi without NP_RMPI_SKIP_INIT,",
        "or call this only in an MPI-initialized subprocess."
      ),
      call. = FALSE
    )
  }

  mode <- match.arg(mode)
  autodispatch.option.sync <- match.arg(autodispatch.option.sync)
  world.size <- .npRmpi_safe_int(mpi.comm.size(0))
  world.size <- if (is.na(world.size)) 1L else as.integer(world.size)

  if (identical(mode, "auto")) {
    mode <- if (world.size > 1L) "attach" else "spawn"
  }

  effective.autodispatch <- autodispatch

  .npRmpi_session_apply_options(
    autodispatch = effective.autodispatch,
    np.messages = np.messages,
    autodispatch.verify.options = autodispatch.verify.options,
    autodispatch.option.sync = autodispatch.option.sync
  )
  if (isTRUE(effective.autodispatch))
    .npRmpi_autodispatch_prime_options()
  else
    options(npRmpi.autodispatch.option.snapshot = NULL)
  .npRmpi_session_reset_spmd_state()
  .npRmpi_attach_state_reset()

  {
    clear.fun <- tryCatch(
      get(".npRmpi_profile_clear", envir = asNamespace("npRmpi"), inherits = FALSE),
      error = function(e) NULL
    )
    if (is.function(clear.fun))
      tryCatch(clear.fun(), error = function(e) NULL)
  }

  if (identical(mode, "spawn")) {
    if (.npRmpi_session_has_active_pool(comm = comm)) {
      .npRmpi_session_reset_spmd_state_allranks(
        comm = comm,
        strict = TRUE,
        where = "npRmpi.init() reused-slave SPMD reset"
      )
      options(npRmpi.master.only = FALSE)
      if (!quiet)
        mpi.hostinfo(comm)
      return(invisible(TRUE))
    }
    mpi.spawn.Rslaves(..., nslaves = nslaves, comm = comm, quiet = quiet, nonblock = nonblock, sleep = sleep)
    mpi.bcast.cmd(np.mpi.initialize(), caller.execute = TRUE, comm = comm)
    options(npRmpi.master.only = FALSE)
    return(invisible(TRUE))
  }

  if (world.size < 2L)
    stop("attach mode requires a pre-launched MPI world (e.g. mpiexec -n <master+slaves>)")

  # Attach mode standardizes on communicator 1 used by the C core.
  # Ensure it exists under mpiexec-launched sessions without requiring
  # user-side mpi.comm.dup(...) boilerplate in scripts.
  comm <- 1L
  .npRmpi_safe(mpi.comm.dup(0, comm), fallback = NULL)
  .npRmpi_session_ensure_comm(comm = comm)
  np.mpi.initialize()
  mpi.barrier(0)
  options(npRmpi.master.only = FALSE)

  rank <- .npRmpi_safe_int(mpi.comm.rank(comm))
  rank <- if (is.na(rank)) 0L else as.integer(rank)
  sid.local <- if (rank == 0L) .npRmpi_attach_next_session_id() else 0L
  sid.msg <- .npRmpi_safe(mpi.bcast(as.integer(sid.local), type = 1, rank = 0, comm = comm),
                          fallback = as.integer(sid.local))
  sid <- suppressWarnings(as.integer(sid.msg)[1L])
  if (!is.finite(sid) || sid < 1L)
    sid <- as.integer(sid.local)
  options(npRmpi.attach.session.id = sid)
  options(npRmpi.attach.close.state = "open")

  if (rank == 0L) {
    if (!quiet) slave.hostinfo(comm)
    return(invisible(TRUE))
  }

  .npRmpi_session_attach_worker_loop(comm = comm, nonblock = nonblock, sleep = sleep, session_id = sid)
}

npRmpi.quit <- function(force = FALSE,
                        dellog = TRUE,
                        comm = 1,
                        mode = c("auto", "spawn", "attach")) {
  force <- npValidateScalarLogical(force, "force")
  dellog <- npValidateScalarLogical(dellog, "dellog")
  comm <- npValidatePositiveInteger(comm, "comm")
  mode <- match.arg(mode)
  size.comm <- .npRmpi_safe_int(mpi.comm.size(comm))
  size.comm <- if (is.na(size.comm)) 0L else as.integer(size.comm)
  size.world <- .npRmpi_safe_int(mpi.comm.size(0))
  size.world <- if (is.na(size.world)) 1L else as.integer(size.world)

  if (identical(mode, "auto")) {
    mode <- if (size.world > 1L && size.comm > 1L) "attach" else "spawn"
  }

  {
    clear.fun <- tryCatch(
      get(".npRmpi_profile_clear", envir = asNamespace("npRmpi"), inherits = FALSE),
      error = function(e) NULL
    )
    if (is.function(clear.fun))
      tryCatch(clear.fun(), error = function(e) NULL)
  }
  options(npRmpi.autodispatch.option.snapshot = NULL)

  if (size.comm < 2L)
  {
    .npRmpi_session_reset_spmd_state()
    .npRmpi_attach_state_reset()
    options(npRmpi.master.only = FALSE)
    return(invisible(FALSE))
  }

  if (identical(mode, "attach")) {
    comm <- 1L
    attach.state <- as.character(getOption("npRmpi.attach.close.state", "open"))[1L]
    if (identical(attach.state, "closed")) {
      .npRmpi_session_reset_spmd_state()
      .npRmpi_attach_state_reset()
      options(npRmpi.master.only = FALSE)
      return(invisible(TRUE))
    }
    if (identical(attach.state, "closing")) {
      .np_warning("attach-mode close already in progress; ignoring duplicate npRmpi.quit(mode='attach') call.",
              call. = FALSE)
      return(invisible(FALSE))
    }
    options(npRmpi.attach.close.state = "closing")

    size.attach <- .npRmpi_safe_int(mpi.comm.size(comm))
    size.attach <- if (is.na(size.attach)) 0L else as.integer(size.attach)
    expected <- if (size.attach >= 2L) seq_len(size.attach - 1L) else integer(0)
    sid <- suppressWarnings(as.integer(getOption("npRmpi.attach.session.id", NA_integer_))[1L])

    recv.timeout <- .npRmpi_session_recv_timeout()
    if (recv.timeout > 0)
      base::setTimeLimit(elapsed = recv.timeout, transient = TRUE)
    shut <- try(mpi.bcast.cmd(cmd = "kaerb", rank = 0, comm = comm, caller.execute = FALSE),
                silent = TRUE)
    if (recv.timeout > 0)
      base::setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    if (inherits(shut, "try-error")) {
      msg <- as.character(shut)
      .np_warning(
        paste(
          "attach-mode worker shutdown broadcast failed or timed out;",
          "continuing communicator teardown.",
          paste(msg, collapse = " ")
        ),
        call. = FALSE
      )
    }

    ack.info <- .npRmpi_attach_collect_close_acks(
      comm = comm,
      session_id = sid,
      expected = expected,
      timeout_sec = .npRmpi_attach_close_ack_timeout()
    )
    release.info <- .npRmpi_attach_send_close_release(
      comm = comm,
      ranks = ack.info$acked,
      session_id = sid,
      barrier = isTRUE(ack.info$ok)
    )
    if (!isTRUE(ack.info$ok)) {
      .np_warning(
        sprintf(
          "attach-mode worker shutdown ACK incomplete for session %s; missing ranks: %s",
          if (is.na(sid)) "NA" else as.character(sid),
          paste(ack.info$missing, collapse = ",")
        ),
        call. = FALSE
      )
    }
    if (!isTRUE(release.info$ok)) {
      .np_warning(
        sprintf(
          "attach-mode worker shutdown release incomplete for session %s; failed ranks: %s",
          if (is.na(sid)) "NA" else as.character(sid),
          paste(release.info$failed, collapse = ",")
        ),
        call. = FALSE
      )
    }
    if (isTRUE(ack.info$ok) && isTRUE(release.info$ok))
      .npRmpi_safe(mpi.barrier(0), fallback = NULL)

    if (comm != 0L) {
      .npRmpi_safe(mpi.comm.free(comm), fallback = NULL)
    }
    .npRmpi_session_reset_spmd_state()
    .npRmpi_attach_state_reset()
    options(npRmpi.master.only = FALSE)
    return(invisible(TRUE))
  }

  mpi.close.Rslaves(dellog = dellog, comm = comm, force = force)
  .npRmpi_session_reset_spmd_state()
  .npRmpi_attach_state_reset()
  options(npRmpi.master.only = FALSE)
  invisible(TRUE)
}

npRmpi.session.info <- function(comm=1){
  np_ver <- .npRmpi_safe(utils::packageVersion("npRmpi"), fallback = NA)
  rmpi_ver <- getOption("npRmpi.embedded.backend.version", NA_character_)
  reuse <- isTRUE(getOption("npRmpi.reuse.slaves", FALSE))
  env_no_reuse <- Sys.getenv("NP_RMPI_NO_REUSE_SLAVES", unset="")

  mpi_ver <- .npRmpi_safe(mpi.get.version(), fallback = NA)
  mpi_initialized <- isTRUE(getOption("npRmpi.mpi.initialized", FALSE))
  if (mpi_initialized) {
    comm_size <- .npRmpi_safe_int(mpi.comm.size(comm))
    comm_rank <- .npRmpi_safe_int(mpi.comm.rank(comm))
    proc <- .npRmpi_safe(mpi.get.processor.name(), fallback = NA_character_)
  } else {
    comm_size <- NA_integer_
    comm_rank <- NA_integer_
    proc <- NA_character_
  }

  info <- list(
    npRmpi = if (is.na(np_ver)[1L]) NA_character_ else as.character(np_ver),
    Rmpi = if (is.na(rmpi_ver)[1L]) NA_character_ else as.character(rmpi_ver)[1L],
    platform = R.version$platform,
    sysname = Sys.info()[["sysname"]],
    release = Sys.info()[["release"]],
    reuse_slaves = reuse,
    NP_RMPI_NO_REUSE_SLAVES = env_no_reuse,
    mpi_version = if (is.na(mpi_ver)[1L]) NA else mpi_ver,
    comm = comm,
    comm_rank = comm_rank,
    comm_size = comm_size,
    nslaves = if (is.na(comm_size)) NA_integer_ else max(as.integer(comm_size) - 1L, 0L),
    master_only = isTRUE(getOption("npRmpi.master.only", FALSE)),
    processor = proc
  )

  cat("npRmpi session info\n")
  cat("  npRmpi:", info$npRmpi, "\n")
  cat("  Rmpi  :", info$Rmpi, "\n")
  cat("  OS    :", info$sysname, info$release, "|", info$platform, "\n")
  cat("  reuse :", info$reuse_slaves, "(NP_RMPI_NO_REUSE_SLAVES=", info$NP_RMPI_NO_REUSE_SLAVES, ")\n", sep="")
  cat("  comm  :", info$comm, "rank", info$comm_rank, "of", info$comm_size, "(nslaves=", info$nslaves, ")\n", sep=" ")
  cat("  mode  :", if (isTRUE(info$master_only)) "master-only" else "slave-pool", "\n")
  if (!is.na(info$processor))
    cat("  host  :", info$processor, "\n")
  if (!all(is.na(info$mpi_version)))
    cat("  mpi   :", paste(info$mpi_version, collapse=" "), "\n")

  invisible(info)
}
