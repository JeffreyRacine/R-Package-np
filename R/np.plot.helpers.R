## the idea is that you have done bandwidth selection
## you just need to supply training data and the bandwidth
## this tool will help you visualize the result

.np_seed_enter <- function(random.seed = 42L) {
  save.seed <- NULL
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    save.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    exists.seed <- TRUE
  } else {
    exists.seed <- FALSE
  }

  set.seed(random.seed)
  list(exists.seed = exists.seed, save.seed = save.seed)
}

.np_seed_exit <- function(state, remove_if_absent = FALSE) {
  if (isTRUE(state$exists.seed)) {
    assign(".Random.seed", state$save.seed, envir = .GlobalEnv)
  } else if (isTRUE(remove_if_absent) &&
             exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(NULL)
}

.np_with_seed <- function(random.seed = 42L, code) {
  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state), add = TRUE)
  force(code)
}

.np_plot_progress_enabled <- function() {
  if (!interactive() && !isTRUE(getOption("np.plot.progress.noninteractive", FALSE)))
    return(FALSE)
  isTRUE(getOption("np.messages")) &&
    isTRUE(getOption("np.plot.progress", TRUE))
}

.np_plot_progress_interval_sec <- function() {
  val <- suppressWarnings(as.numeric(getOption("np.plot.progress.interval.sec", 0.5))[1L])
  if (!is.finite(val) || is.na(val) || val < 0)
    val <- 0.5
  val
}

.np_plot_progress_begin <- function(total, label) {
  total <- as.integer(total)
  if (is.na(total) || total < 1L || !.np_plot_progress_enabled())
    return(NULL)

  list(
    total = total,
    label = as.character(label)[1L],
    started = unname(as.double(proc.time()[["elapsed"]])),
    last = -Inf,
    interval = .np_plot_progress_interval_sec(),
    console = newLineConsole()
  )
}

.np_plot_progress_tick <- function(state, done, force = FALSE) {
  if (is.null(state))
    return(state)

  done <- as.integer(done)
  if (is.na(done))
    done <- 0L
  done <- max(0L, min(state$total, done))

  now <- unname(as.double(proc.time()[["elapsed"]]))
  if (!isTRUE(force) &&
      done < state$total &&
      (now - state$last) < state$interval) {
    return(state)
  }

  elapsed <- max(0, now - state$started)
  rate <- if (elapsed > 0) done / elapsed else 0
  remain <- max(0L, state$total - done)
  eta <- if (rate > 0) remain / rate else Inf
  eta.txt <- if (is.finite(eta)) sprintf("%.1fs", eta) else "NA"
  pct <- 100 * done / state$total

  msg <- sprintf(
    "%s %d/%d (%.1f%%, elapsed %.1fs, eta %s)",
    state$label, done, state$total, pct, elapsed, eta.txt
  )

  state$console <- printClear(state$console)
  state$console <- printPush(msg = msg, console = state$console)
  state$last <- now
  state
}

.np_plot_progress_end <- function(state) {
  if (is.null(state))
    return(invisible(NULL))
  state <- .np_plot_progress_tick(state = state, done = state$total, force = TRUE)
  state$console <- printClear(state$console)
  invisible(NULL)
}

.np_plot_restore_par <- function(oldpar) {
  if (is.null(oldpar))
    return(invisible(NULL))
  if (is.list(oldpar) && !is.null(oldpar[["new"]]))
    oldpar[["new"]] <- FALSE
  suppressWarnings(try(par(oldpar), silent = TRUE))
  invisible(NULL)
}

.np_mammen_draws <- function(n, B) {
  a <- (1 - sqrt(5)) / 2
  p.a <- (sqrt(5) + 1) / (2 * sqrt(5))
  u <- matrix(stats::runif(n * B), nrow = n, ncol = B)
  out <- matrix(1 - a, nrow = n, ncol = B)
  out[u <= p.a] <- a
  out
}

.np_rademacher_draws <- function(n, B) {
  u <- matrix(stats::runif(n * B), nrow = n, ncol = B)
  out <- matrix(1.0, nrow = n, ncol = B)
  out[u <= 0.5] <- -1.0
  out
}

.np_wild_draws <- function(n, B, wild = c("mammen", "rademacher")) {
  if (length(wild) > 1L)
    wild <- wild[1L]
  wild <- match.arg(wild, c("mammen", "rademacher"))
  if (identical(wild, "mammen")) {
    return(.np_mammen_draws(n = n, B = B))
  }
  .np_rademacher_draws(n = n, B = B)
}

.np_wild_chunk_size <- function(n, B) {
  chunk.opt <- getOption("np.plot.wild.chunk.size")
  if (!is.null(chunk.opt)) {
    chunk.opt <- as.integer(chunk.opt)
    if (length(chunk.opt) != 1L || is.na(chunk.opt) || chunk.opt < 1L)
      stop("option 'np.plot.wild.chunk.size' must be a positive integer")
    return(min(B, chunk.opt))
  }

  if (n < 1L || B < 1L)
    return(1L)

  # Keep the temporary n x chunk response matrix in a moderate memory range.
  target.bytes <- 64 * 1024 * 1024
  chunk <- as.integer(floor(target.bytes / (8 * n)))
  if (!is.finite(chunk) || is.na(chunk) || chunk < 1L)
    chunk <- 1L

  # In active MPI sessions, very large wild chunks have triggered allocator
  # instability on some stacks; keep a conservative default cap.
  if (isTRUE(getOption("npRmpi.mpi.initialized", FALSE))) {
    mpi.cap <- suppressWarnings(as.integer(getOption("np.plot.wild.chunk.max.mpi", 64L))[1L])
    if (is.na(mpi.cap) || mpi.cap < 1L)
      mpi.cap <- 64L
    chunk <- min(chunk, mpi.cap)
  }

  min(B, chunk)
}

.np_wild_boot_t <- function(H, fit.mean, residuals, B, wild = c("mammen", "rademacher")) {
  B <- as.integer(B)
  n <- length(residuals)
  if (length(fit.mean) != n)
    stop("length mismatch between fitted means and residuals for wild bootstrap")
  if (B < 1L)
    stop("argument 'plot.errors.boot.num' must be a positive integer")
  .npRmpi_bootstrap_transport_trace(
    what = "wild",
    event = "wild.entry",
    fields = list(
      n = n,
      h_rows = nrow(H),
      h_cols = ncol(H),
      B = B
    )
  )

  t.mpi <- .npRmpi_wild_boot_t_parallel(
    H = H,
    fit.mean = fit.mean,
    residuals = residuals,
    B = B,
    wild = wild,
    comm = 1L
  )
  if (is.matrix(t.mpi)) {
    .npRmpi_bootstrap_transport_trace(
      what = "wild",
      event = "wild.return",
      fields = list(
        t_rows = nrow(t.mpi),
        t_cols = ncol(t.mpi)
      )
    )
    return(t.mpi)
  }
  .npRmpi_bootstrap_fail_or_fallback(
    msg = "wild bootstrap fan-out did not return matrix output",
    what = "wild"
  )
}

.np_plot_is_wild_method <- function(method) {
  isTRUE(length(method) == 1L && !is.na(method) && method == "wild")
}

.np_plot_reject_wild_unsupervised <- function(method, where) {
  if (.np_plot_is_wild_method(method)) {
    stop(sprintf("plot.errors.boot.method='wild' is not supported for %s; use one of 'inid', 'fixed', or 'geom'", where))
  }
  invisible(NULL)
}

.np_plot_inid_fastpath_enabled <- function() {
  TRUE
}

.np_plot_block_fastpath_enabled <- function() {
  TRUE
}

.np_plot_require_bws <- function(bws, where) {
  if (is.null(bws))
    stop(sprintf("required argument 'bws' is missing or NULL in %s", where))
  invisible(TRUE)
}

# Policy hook: fastpath remains default-on; keep as helper so policy can be
# centralized without touching call sites.
.npRmpi_plot_inid_ksum_fastpath_enabled <- function() {
  TRUE
}

.npRmpi_bootstrap_worker_count <- function(comm = 1L) {
  if (!isTRUE(getOption("npRmpi.mpi.initialized", FALSE)))
    return(0L)
  size <- tryCatch(as.integer(mpi.comm.size(comm = comm)), error = function(e) NA_integer_)
  if (is.na(size) || size <= 1L)
    return(0L)
  size - 1L
}

.npRmpi_bootstrap_tune_chunk_size <- function(B,
                                              chunk.size,
                                              comm = 1L,
                                              include.master = TRUE) {
  B <- as.integer(B)
  chunk.size <- as.integer(chunk.size)
  if (is.na(B) || B < 1L)
    return(1L)
  if (is.na(chunk.size) || chunk.size < 1L)
    chunk.size <- 1L

  workers <- .npRmpi_bootstrap_worker_count(comm = comm)
  slots <- workers + if (isTRUE(include.master)) 1L else 0L
  if (slots > 1L && B >= slots) {
    max.chunk <- max(1L, as.integer(floor(B / slots)))
    chunk.size <- min(chunk.size, max.chunk)
  }

  max(1L, min(B, chunk.size))
}

.npRmpi_bootstrap_fail_or_fallback <- function(msg, what = "bootstrap") {
  # Fail-fast by design: MPI-selected routes must not silently fallback.
  stop(sprintf("MPI %s %s", what, msg), call. = FALSE)
}

.npRmpi_bootstrap_dispatch_timeout_sec <- function() {
  val.env <- Sys.getenv("NP_RMPI_BOOTSTRAP_DISPATCH_TIMEOUT_SEC", unset = "")
  val.opt <- getOption("npRmpi.bootstrap.dispatch.timeout.sec", NA_real_)
  val <- if (nzchar(val.env)) val.env else val.opt
  val <- suppressWarnings(as.numeric(val)[1L])
  if (!is.finite(val) || is.na(val) || val <= 0)
    return(0.0)
  val
}

.npRmpi_bootstrap_phase_trace_path <- function() {
  path.opt <- getOption("npRmpi.bootstrap.phase.file", "")
  path.env <- Sys.getenv("NP_RMPI_BOOTSTRAP_PHASE_FILE", unset = "")
  path <- if (nzchar(path.env)) path.env else path.opt
  path <- as.character(path)[1L]
  if (is.na(path) || !nzchar(path))
    return("")
  path
}

.npRmpi_bootstrap_phase_trace_append <- function(what, phase, where) {
  path <- .npRmpi_bootstrap_phase_trace_path()
  if (!nzchar(path))
    return(invisible(FALSE))

  dirpath <- dirname(path)
  if (!identical(dirpath, ".") && !dir.exists(dirpath)) {
    ok <- tryCatch({
      dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
    }, error = function(e) FALSE)
    if (!isTRUE(ok) && !dir.exists(dirpath))
      return(invisible(FALSE))
  }

  line <- sprintf(
    "%s\tpid=%d\twhat=%s\tphase=%s\twhere=%s\n",
    format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6%z"),
    Sys.getpid(),
    as.character(what)[1L],
    as.character(phase)[1L],
    if (is.na(where) || !nzchar(where)) "" else as.character(where)[1L]
  )

  tryCatch({
    cat(line, file = path, append = TRUE)
    TRUE
  }, error = function(e) FALSE)
}

.npRmpi_bootstrap_transport_trace_path <- function() {
  path.opt <- getOption("npRmpi.bootstrap.transport.trace.file", "")
  path.env <- Sys.getenv("NP_RMPI_BOOTSTRAP_TRANSPORT_TRACE_FILE", unset = "")
  path.fallback <- Sys.getenv("NP_RMPI_TRANSPORT_TRACE_FILE", unset = "")
  path <- if (nzchar(path.env)) {
    path.env
  } else if (nzchar(path.fallback)) {
    path.fallback
  } else {
    path.opt
  }
  path <- as.character(path)[1L]
  if (is.na(path) || !nzchar(path))
    return("")
  path
}

.npRmpi_bootstrap_transport_trace <- function(what, event, fields = list()) {
  path <- .npRmpi_bootstrap_transport_trace_path()
  if (!nzchar(path))
    return(invisible(FALSE))

  dirpath <- dirname(path)
  if (!identical(dirpath, ".") && !dir.exists(dirpath)) {
    ok <- tryCatch({
      dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
    }, error = function(e) FALSE)
    if (!isTRUE(ok) && !dir.exists(dirpath))
      return(invisible(FALSE))
  }

  if (is.null(names(fields)))
    names(fields) <- paste0("f", seq_along(fields))
  if (length(fields) > 0)
    fields <- fields[!is.na(names(fields)) & nzchar(names(fields))]
  kv <- if (length(fields) > 0) {
    paste(
      paste0(
        names(fields),
        "=",
        vapply(fields, function(v) {
          vv <- as.character(v)[1L]
          if (is.na(vv)) "NA" else vv
        }, character(1))
      ),
      collapse = "\t"
    )
  } else ""

  line <- paste(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6%z"),
    paste0("pid=", Sys.getpid()),
    paste0("what=", as.character(what)[1L]),
    paste0("event=", as.character(event)[1L]),
    kv,
    sep = "\t"
  )

  tryCatch({
    cat(paste0(line, "\n"), file = path, append = TRUE)
    TRUE
  }, error = function(e) FALSE)
}

.npRmpi_bootstrap_phase_mark <- function(what = "bootstrap",
                                         phase = NA_character_,
                                         where = NA_character_) {
  phase <- as.character(phase)[1L]
  if (is.na(phase) || !nzchar(phase))
    return(invisible(FALSE))

  label <- paste0("phase:", what, ":", phase)
  where <- as.character(where)[1L]
  if (!is.na(where) && nzchar(where))
    label <- paste0(label, ":", where)

  .npRmpi_profile_add_comm_elapsed(elapsed_sec = 0.0, where = label)
  .npRmpi_bootstrap_phase_trace_append(what = what, phase = phase, where = where)
  invisible(TRUE)
}

.npRmpi_bootstrap_assert_bindings <- function(required.bindings = NULL,
                                              what = "bootstrap") {
  if (is.null(required.bindings))
    return(invisible(TRUE))

  if (!is.list(required.bindings))
    stop("required.bindings must be a named list")

  nms <- names(required.bindings)
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms)))
    stop("required.bindings must be a named list with non-empty names")

  missing <- nms[vapply(required.bindings, is.null, logical(1))]
  if (length(missing)) {
    .npRmpi_bootstrap_fail_or_fallback(
      msg = sprintf("missing required worker bindings (%s)",
                    paste(missing, collapse = ", ")),
      what = what
    )
  }

  if (isTRUE(getOption("npRmpi.bootstrap.binding.serialize.check", FALSE))) {
    bad <- nms[!vapply(
      required.bindings,
      function(obj) isTRUE(tryCatch({
        serialize(obj, NULL)
        TRUE
      }, error = function(e) FALSE)),
      logical(1)
    )]
    if (length(bad)) {
      .npRmpi_bootstrap_fail_or_fallback(
        msg = sprintf("non-serializable worker bindings (%s)",
                      paste(bad, collapse = ", ")),
        what = what
      )
    }
  }

  invisible(TRUE)
}

.npRmpi_bootstrap_prepare_worker <- function(worker,
                                             required.bindings = NULL) {
  if (!is.function(worker))
    stop("worker must be a function")

  bind.env <- new.env(parent = asNamespace("npRmpi"))
  if (!is.null(required.bindings)) {
    nms <- names(required.bindings)
    if (is.null(nms))
      stop("required.bindings must be named when provided")
    for (nm in nms)
      assign(nm, required.bindings[[nm]], envir = bind.env)
  }

  worker.prepared <- worker
  environment(worker.prepared) <- bind.env
  worker.prepared
}

.npRmpi_bootstrap_fanout_enabled <- function(comm = 1L,
                                             n = NA_integer_,
                                             B = NA_integer_,
                                             chunk.size = NA_integer_,
                                             what = "bootstrap") {
  if (isTRUE(.npRmpi_autodispatch_called_from_bcast()))
    .npRmpi_bootstrap_fail_or_fallback(
      msg = "cannot run inside mpi.bcast.cmd context; invoke plot/bootstrap from master context so work can be fanned out across workers",
      what = what
    )
  if (!isTRUE(.npRmpi_has_active_slave_pool(comm = comm)))
    if (!isTRUE(.npRmpi_master_only_mode(comm = comm)))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "requires an active MPI slave pool; call npRmpi.init(...) first",
        what = what
      )
  workers <- .npRmpi_bootstrap_worker_count(comm = comm)
  if (!isTRUE(.npRmpi_master_only_mode(comm = comm)) && workers < 1L)
    .npRmpi_bootstrap_fail_or_fallback(
      msg = "requires at least one active MPI worker",
      what = what
    )
  TRUE
}

.npRmpi_bootstrap_chunk_tasks <- function(B, chunk.size) {
  B <- as.integer(B)
  chunk.size <- as.integer(chunk.size)
  if (B < 1L || chunk.size < 1L)
    stop("invalid chunk configuration")

  starts <- seq.int(1L, B, by = chunk.size)
  lens <- pmin(chunk.size, B - starts + 1L)
  seeds <- sample.int(.Machine$integer.max, length(starts))

  lapply(seq_along(starts), function(i) {
    list(
      start = as.integer(starts[i]),
      bsz = as.integer(lens[i]),
      seed = as.integer(seeds[i])
    )
  })
}

.npRmpi_bootstrap_collect_chunks <- function(parts, tasks, ncol.out, what = "bootstrap") {
  if (!is.list(parts) || length(parts) != length(tasks)) {
    .npRmpi_bootstrap_fail_or_fallback(
      msg = "fan-out returned malformed chunk results",
      what = what
    )
    return(NULL)
  }

  has.try.error <- vapply(parts, function(x) inherits(x, "try-error"), logical(1))
  if (any(has.try.error)) {
    first.err <- parts[[which(has.try.error)[1L]]]
    first.msg <- tryCatch(as.character(first.err)[1L], error = function(e) NA_character_)
    if (is.na(first.msg) || !nzchar(first.msg))
      first.msg <- "worker error"
    .npRmpi_bootstrap_fail_or_fallback(
      msg = sprintf("fan-out worker error detected (%s)", first.msg),
      what = what
    )
    return(NULL)
  }

  total.rows <- sum(vapply(tasks, function(tt) as.integer(tt$bsz), integer(1)))
  out <- matrix(NA_real_, nrow = total.rows, ncol = as.integer(ncol.out))
  rowi <- 1L

  for (i in seq_along(parts)) {
    bsz <- as.integer(tasks[[i]]$bsz)
    chunk <- parts[[i]]
    if (!is.matrix(chunk))
      chunk <- as.matrix(chunk)

    if (!identical(dim(chunk), c(bsz, as.integer(ncol.out)))) {
      if (length(chunk) != (bsz * as.integer(ncol.out))) {
        .npRmpi_bootstrap_fail_or_fallback(
          msg = "fan-out chunk dimension mismatch",
          what = what
        )
        return(NULL)
      }
      chunk <- matrix(as.numeric(chunk), nrow = bsz, ncol = as.integer(ncol.out))
    }

    out[rowi:(rowi + bsz - 1L), ] <- chunk
    rowi <- rowi + bsz
  }

  out
}

.npRmpi_bootstrap_run_fanout <- function(tasks,
                                         worker,
                                         ncol.out,
                                         what = "bootstrap",
                                         profile.where = NA_character_,
                                         comm = 1L,
                                         required.bindings = NULL,
                                         ...) {
  total.boot <- sum(vapply(tasks, function(tt) as.integer(tt$bsz), integer(1)))
  progress <- .np_plot_progress_begin(
    total = total.boot,
    label = sprintf("Plot bootstrap %s", what)
  )
  on.exit({
    .np_plot_progress_end(progress)
  }, add = TRUE)

  .npRmpi_bootstrap_phase_mark(
    what = what,
    phase = "preflight",
    where = profile.where
  )
  .npRmpi_bootstrap_assert_bindings(
    required.bindings = required.bindings,
    what = what
  )
  worker.exec <- .npRmpi_bootstrap_prepare_worker(
    worker = worker,
    required.bindings = required.bindings
  )

  .npRmpi_bootstrap_phase_mark(
    what = what,
    phase = "dispatch",
    where = profile.where
  )
  workers <- .npRmpi_bootstrap_worker_count(comm = comm)
  use.master.local <- !isTRUE(.npRmpi_has_active_slave_pool(comm = comm)) &&
    isTRUE(.npRmpi_master_only_mode(comm = comm))
  .npRmpi_bootstrap_transport_trace(
    what = what,
    event = "fanout.start",
    fields = list(
      workers = workers,
      tasks = length(tasks),
      master_local = isTRUE(use.master.local),
      comm = comm
    )
  )
  t.comm <- proc.time()
  parts <- if (isTRUE(use.master.local)) {
    tryCatch({
      parts.local <- vector("list", length(tasks))
      done.boot <- 0L
      .npRmpi_bootstrap_transport_trace(
        what = what,
        event = "fanout.master_local.start",
        fields = list(tasks = length(tasks))
      )
      for (ii in seq_along(tasks)) {
        task <- tasks[[ii]]
        parts.local[[ii]] <- do.call(worker.exec, c(list(task), list(...)))
        done.boot <- done.boot + as.integer(task$bsz)
        progress <- .np_plot_progress_tick(state = progress, done = done.boot)
      }
      .npRmpi_bootstrap_transport_trace(
        what = what,
        event = "fanout.master_local.done",
        fields = list(tasks = length(tasks))
      )
      parts.local
    }, error = function(e) e)
  } else if (workers >= 1L && length(tasks) >= 2L) {
    tryCatch({
      n.tasks <- length(tasks)
      remote.idx <- seq_len(n.tasks)
      parts.out <- vector("list", n.tasks)

      if (length(remote.idx) == 0L) {
        parts.out
      } else {
        n.remote <- length(remote.idx)
        slave.num <- workers
        mpi.anysource <- mpi.any.source()
        mpi.anytag <- mpi.any.tag()
        dispatch.timeout <- .npRmpi_bootstrap_dispatch_timeout_sec()
        dispatch.started <- unname(as.double(proc.time()[["elapsed"]]))

        mpi.bcast.cmd(.mpi.worker.applyLB, n = n.remote, comm = comm)
        mpi.bcast.Robj(list(FUN = worker.exec, dot.arg = list(...)), rank = 0, comm = comm)
        .npRmpi_bootstrap_transport_trace(
          what = what,
          event = "fanout.master_assist.start",
          fields = list(n_remote = n.remote, slave_num = slave.num, local_n = 0L)
        )

        init <- min(slave.num, n.remote)
        if (init > 0L) {
          for (i in seq_len(init)) {
            task.idx <- remote.idx[i]
            mpi.send.Robj(list(data.arg = list(tasks[[task.idx]])), dest = i, tag = i, comm = comm)
            .npRmpi_bootstrap_transport_trace(
              what = what,
              event = "fanout.send.initial",
              fields = list(dest = i, tag = i, task_idx = task.idx)
            )
          }
        }

        if (init < slave.num) {
          stop.tag <- as.integer(n.remote + 1L)
          for (i in seq.int(init + 1L, slave.num)) {
            mpi.send.Robj(as.integer(0), dest = i, tag = stop.tag, comm = comm)
            .npRmpi_bootstrap_transport_trace(
              what = what,
              event = "fanout.send.stop.initial",
              fields = list(dest = i, tag = stop.tag)
            )
          }
        }

        sent <- init
        done <- 0L
        done.boot <- 0L

        while (done < n.remote) {
          progressed <- FALSE

          if (done < n.remote && isTRUE(mpi.iprobe(mpi.anysource, mpi.anytag, comm))) {
            srctag <- mpi.get.sourcetag()
            src <- srctag[1L]
            tag <- srctag[2L]
            res <- mpi.recv.Robj(source = src, tag = tag, comm = comm)
            .npRmpi_bootstrap_transport_trace(
              what = what,
              event = "fanout.recv",
              fields = list(src = src, tag = tag, done_next = done + 1L, n_remote = n.remote)
            )
            done <- done + 1L
            task.idx <- remote.idx[tag]
            parts.out[[task.idx]] <- res
            done.boot <- done.boot + as.integer(tasks[[task.idx]]$bsz)
            progress <- .np_plot_progress_tick(state = progress, done = done.boot)

            sent <- sent + 1L
            if (sent <= n.remote) {
              next.idx <- remote.idx[sent]
              mpi.send.Robj(list(data.arg = list(tasks[[next.idx]])), dest = src, tag = sent, comm = comm)
              .npRmpi_bootstrap_transport_trace(
                what = what,
                event = "fanout.send.next",
                fields = list(dest = src, tag = sent, task_idx = next.idx)
              )
            } else {
              mpi.send.Robj(as.integer(0), dest = src, tag = as.integer(n.remote + 1L), comm = comm)
              .npRmpi_bootstrap_transport_trace(
                what = what,
                event = "fanout.send.stop",
                fields = list(dest = src, tag = as.integer(n.remote + 1L))
              )
            }
            progressed <- TRUE
          }

          if (!progressed && done < n.remote) {
            if (dispatch.timeout > 0) {
              elapsed.wait <- unname(as.double(proc.time()[["elapsed"]])) - dispatch.started
              if (is.finite(elapsed.wait) && elapsed.wait > dispatch.timeout) {
                .npRmpi_bootstrap_fail_or_fallback(
                  msg = sprintf(
                    "dispatch timeout waiting on worker results (done=%d/%d, timeout=%.3fs)",
                    done,
                    n.remote,
                    dispatch.timeout
                  ),
                  what = what
                )
              }
            }
            Sys.sleep(0.0005)
          }
        }
        .npRmpi_bootstrap_transport_trace(
          what = what,
          event = "fanout.master_assist.done",
          fields = list(done = done, n_remote = n.remote, local_done = 0L)
        )
      }

      parts.out
    }, error = function(e) e)
  } else {
    .npRmpi_bootstrap_transport_trace(
      what = what,
      event = "fanout.applylb.fallback",
      fields = list(tasks = length(tasks), workers = workers)
    )
    tryCatch(
      mpi.applyLB(tasks, worker.exec, ..., comm = comm),
      error = function(e) e
    )
  }
  .npRmpi_profile_add_comm_elapsed(
    elapsed_sec = unname(as.double((proc.time() - t.comm)[["elapsed"]])),
    where = if (!is.na(profile.where) && nzchar(profile.where))
      if (isTRUE(use.master.local))
        paste0(profile.where, ":master-local")
      else if (workers >= 1L && length(tasks) >= 2L)
        paste0(profile.where, ":master-assist")
      else
        profile.where
    else if (isTRUE(use.master.local))
      paste0("master.local:", what)
    else if (workers >= 1L && length(tasks) >= 2L)
      paste0("master.assist:", what)
    else
      paste0("mpi.applyLB:", what)
  )
  if (inherits(parts, "error")) {
    .npRmpi_bootstrap_transport_trace(
      what = what,
      event = "fanout.error",
      fields = list(message = conditionMessage(parts))
    )
    .npRmpi_bootstrap_phase_mark(
      what = what,
      phase = "dispatch_error",
      where = profile.where
    )
    .npRmpi_bootstrap_fail_or_fallback(
      msg = sprintf("fan-out failed (%s)", conditionMessage(parts)),
      what = what
    )
  }

  progress <- .np_plot_progress_tick(state = progress, done = total.boot, force = TRUE)

  .npRmpi_bootstrap_phase_mark(
    what = what,
    phase = "collect",
    where = profile.where
  )
  out <- .npRmpi_bootstrap_collect_chunks(
    parts = parts,
    tasks = tasks,
    ncol.out = ncol.out,
    what = what
  )
  .npRmpi_bootstrap_phase_mark(
    what = what,
    phase = "done",
    where = profile.where
  )
  .npRmpi_bootstrap_transport_trace(
    what = what,
    event = "fanout.done",
    fields = list(parts = length(parts), total_boot = total.boot)
  )
  out
}

.npRmpi_wild_boot_t_parallel <- function(H, fit.mean, residuals, B, wild, comm = 1L) {
  n <- length(residuals)
  p <- nrow(H)
  chunk.size <- .np_wild_chunk_size(n = n, B = B)
  chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
    B = B,
    chunk.size = chunk.size,
    comm = comm,
    include.master = TRUE
  )
  .npRmpi_bootstrap_fanout_enabled(
    comm = comm,
    n = n,
    B = B,
    chunk.size = chunk.size,
    what = "wild"
  )
  tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
  if (!length(tasks))
    .npRmpi_bootstrap_fail_or_fallback(
      msg = "wild fan-out produced no tasks",
      what = "wild"
    )

  H <- as.matrix(H)
  H.vec <- as.double(H)
  n.h <- ncol(H)
  p.h <- nrow(H)
  fit.mean <- as.double(fit.mean)
  residuals <- as.double(residuals)
  wild <- match.arg(if (length(wild) > 1L) wild[1L] else wild,
                    c("mammen", "rademacher"))

  worker <- function(task, n, p, H.vec, fit.mean, residuals, wild.method) {
    H <- matrix(as.double(H.vec), nrow = p, ncol = n)
    set.seed(as.integer(task$seed))
    bsz <- as.integer(task$bsz)
    u <- matrix(stats::runif(n * bsz), nrow = n, ncol = bsz)
    if (identical(wild.method, "mammen")) {
      a <- (1 - sqrt(5)) / 2
      p.a <- (sqrt(5) + 1) / (2 * sqrt(5))
      draws <- matrix(1 - a, nrow = n, ncol = bsz)
      draws[u <= p.a] <- a
    } else {
      draws <- matrix(1, nrow = n, ncol = bsz)
      draws[u <= 0.5] <- -1
    }
    ystar <- matrix(fit.mean, nrow = n, ncol = bsz) +
      matrix(residuals, nrow = n, ncol = bsz) * draws
    t(H %*% ystar)
  }

  .npRmpi_bootstrap_run_fanout(
    tasks = tasks,
    worker = worker,
    ncol.out = p,
    what = "wild",
    profile.where = "mpi.applyLB:wild",
    n = n.h,
    p = p.h,
    H.vec = H.vec,
    fit.mean = fit.mean,
    residuals = residuals,
    wild.method = wild,
    comm = comm
  )
}

.np_inid_chunk_size <- function(n, B) {
  chunk.opt <- getOption("np.plot.inid.chunk.size")
  if (!is.null(chunk.opt)) {
    chunk.opt <- as.integer(chunk.opt)
    if (length(chunk.opt) != 1L || is.na(chunk.opt) || chunk.opt < 1L)
      stop("option 'np.plot.inid.chunk.size' must be a positive integer")
    return(min(B, chunk.opt))
  }

  if (n < 1L || B < 1L)
    return(1L)

  target.bytes <- 64 * 1024 * 1024
  chunk <- as.integer(floor(target.bytes / (8 * n)))
  if (!is.finite(chunk) || is.na(chunk) || chunk < 1L)
    chunk <- 1L
  min(B, chunk)
}

.np_inid_counts_matrix <- function(n, B, counts = NULL) {
  n <- as.integer(n)
  B <- as.integer(B)
  if (n < 1L || B < 1L)
    stop("invalid inid bootstrap dimensions")

  if (!is.null(counts)) {
    counts <- as.matrix(counts)
    if (!is.numeric(counts) ||
        nrow(counts) != n ||
        ncol(counts) != B)
      stop("counts must be an n x B numeric matrix")
    return(counts)
  }

  stats::rmultinom(n = B, size = n, prob = rep.int(1 / n, n))
}

.np_block_counts_drawer <- function(n,
                                    B,
                                    blocklen,
                                    sim = c("fixed", "geom"),
                                    n.sim = n,
                                    endcorr = TRUE) {
  sim <- match.arg(sim)
  n <- as.integer(n)
  B <- as.integer(B)
  n.sim <- as.integer(n.sim)
  blocklen <- as.integer(blocklen)

  if (n < 1L || B < 1L || n.sim < 1L)
    stop("invalid block bootstrap dimensions")
  if (length(blocklen) != 1L || is.na(blocklen) || blocklen < 1L || blocklen > n)
    stop("invalid block length for block bootstrap")

  ts.array <- utils::getFromNamespace("ts.array", "boot")
  make.ends <- utils::getFromNamespace("make.ends", "boot")
  ts.draws <- ts.array(
    n = n,
    n.sim = n.sim,
    R = B,
    l = blocklen,
    sim = sim,
    endcorr = isTRUE(endcorr)
  )

  starts <- as.matrix(ts.draws$starts)
  lengths <- ts.draws$lengths

  function(start, stopi) {
    start <- as.integer(start)
    stopi <- as.integer(stopi)
    if (start < 1L || stopi < start || stopi > B)
      stop("invalid block bootstrap chunk bounds")

    idx <- seq.int(start, stopi)
    out <- matrix(0.0, nrow = n, ncol = length(idx))

    for (jj in seq_along(idx)) {
      rr <- idx[jj]
      ends <- if (identical(sim, "geom")) {
        cbind(starts[rr, ], lengths[rr, ])
      } else {
        cbind(starts[rr, ], lengths)
      }
      inds <- apply(ends, 1L, make.ends, n)
      inds <- if (is.list(inds)) {
        as.integer(unlist(inds)[seq_len(n.sim)])
      } else {
        as.integer(inds)[seq_len(n.sim)]
      }
      out[, jj] <- tabulate(inds, nbins = n)
    }

    out
  }
}

.np_inid_lc_boot_from_hat <- function(H, ydat, B, counts = NULL, counts.drawer = NULL) {
  H <- as.matrix(H)
  ydat <- as.double(ydat)
  B <- as.integer(B)
  n <- length(ydat)
  if (B < 1L)
    stop("argument 'plot.errors.boot.num' must be a positive integer")
  if (ncol(H) != n)
    stop("hat matrix columns must match length of ydat")

  W <- t(H)
  Wy <- W * ydat
  t0 <- as.vector(H %*% ydat)

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    den <- crossprod(counts.mat, W)
    num <- crossprod(counts.mat, Wy)
    return(list(
      t = num / pmax(den, .Machine$double.eps),
      t0 = t0
    ))
  }

  chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
    B = B,
    chunk.size = .np_inid_chunk_size(n = n, B = B),
    comm = 1L,
    include.master = TRUE
  )
  if (!is.null(counts.drawer)) {
    tmat <- matrix(NA_real_, nrow = B, ncol = nrow(H))
    W.local <- W
    Wy.local <- Wy
    counts.drawer.local <- counts.drawer
    n.local <- n

    if (.npRmpi_bootstrap_fanout_enabled(
      comm = 1L,
      n = n,
      B = B,
      chunk.size = chunk.size,
      what = "inid-lc-block"
    )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        start <- as.integer(task$start)
        stopi <- start + as.integer(task$bsz) - 1L
        counts.chunk <- .np_inid_counts_matrix(
          n = n.local,
          B = as.integer(task$bsz),
          counts = counts.drawer.local(start, stopi)
        )
        den <- crossprod(counts.chunk, W.local)
        num <- crossprod(counts.chunk, Wy.local)
        num / pmax(den, .Machine$double.eps)
      }

      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = nrow(H),
        what = "inid-lc-block",
        profile.where = "mpi.applyLB:inid-lc-block",
        comm = 1L,
        required.bindings = list(
          W.local = W.local,
          Wy.local = Wy.local,
          counts.drawer.local = counts.drawer.local,
          n.local = n.local
        )
      )
    }

    if (anyNA(tmat)) {
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-lc-block fan-out returned incomplete results",
        what = "inid-lc-block"
      )
    }

    return(list(t = tmat, t0 = t0))
  }

  .npRmpi_bootstrap_fanout_enabled(
    comm = 1L,
    n = n,
    B = B,
    chunk.size = chunk.size,
    what = "inid-lc"
  )
  tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
  prob <- rep.int(1 / n, n)

  worker <- function(task) {
    set.seed(as.integer(task$seed))
    bsz <- as.integer(task$bsz)
    counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
    den <- crossprod(counts.chunk, W)
    num <- crossprod(counts.chunk, Wy)
    num / pmax(den, .Machine$double.eps)
  }

  t.mpi <- .npRmpi_bootstrap_run_fanout(
    tasks = tasks,
    worker = worker,
    ncol.out = nrow(H),
    what = "inid",
    profile.where = "mpi.applyLB:inid",
    comm = 1L,
    required.bindings = list(
      n = n,
      prob = prob,
      W = W,
      Wy = Wy
    )
  )

  list(t = t.mpi, t0 = t0)
}

.np_inid_lp_unpack_sym_row <- function(mrow, p) {
  A <- matrix(0.0, nrow = p, ncol = p)
  idx <- 1L
  for (a in seq_len(p)) {
    for (b in a:p) {
      A[a, b] <- mrow[idx]
      A[b, a] <- mrow[idx]
      idx <- idx + 1L
    }
  }
  A
}

.np_inid_lp_predict_row <- function(A, z, rhs, ridge.grid) {
  ridge.grid <- as.double(ridge.grid)
  if (!length(ridge.grid) || anyNA(ridge.grid) || any(!is.finite(ridge.grid)) || any(ridge.grid < 0))
    stop("argument 'ridge.grid' must be a non-empty non-negative finite numeric vector")
  z <- as.double(z)
  if (!length(z))
    return(NA_real_)
  z1 <- z[1L]

  for (ridge in ridge.grid) {
    Ar <- A
    zr <- z
    if (ridge > 0)
      diag(Ar) <- diag(Ar) + ridge
    if (ridge > 0)
      zr[1L] <- z1 + ridge * z1 / NZD(Ar[1L, 1L])
    beta <- tryCatch(
      drop(solve(Ar, matrix(zr, ncol = 1L))),
      error = function(e) NULL
    )
    if (!is.null(beta) && all(is.finite(beta)))
      return(sum(rhs * beta))
  }

  NA_real_
}

.np_inid_lp_predict_chunk_general <- function(Mvals, Zvals, rhs, ridge.grid) {
  Mvals <- as.matrix(Mvals)
  Zvals <- as.matrix(Zvals)
  rhs <- as.double(rhs)

  bsz <- nrow(Mvals)
  p <- ncol(Zvals)
  out <- numeric(bsz)

  for (ii in seq_len(bsz)) {
    A <- .np_inid_lp_unpack_sym_row(mrow = Mvals[ii, ], p = p)
    out[ii] <- .np_inid_lp_predict_row(
      A = A,
      z = as.double(Zvals[ii, ]),
      rhs = rhs,
      ridge.grid = ridge.grid
    )
  }

  out
}

.np_inid_lp_predict_chunk <- function(Mvals, Zvals, rhs, ridge.grid) {
  Mvals <- as.matrix(Mvals)
  Zvals <- as.matrix(Zvals)
  rhs <- as.double(rhs)

  bsz <- nrow(Mvals)
  p <- ncol(Zvals)
  out <- rep(NA_real_, bsz)

  if (p == 1L) {
    den <- as.double(Mvals[, 1L])
    out <- rhs[1L] * as.double(Zvals[, 1L]) / pmax(den, .Machine$double.eps)
    return(out)
  }

  if (p == 2L) {
    a <- as.double(Mvals[, 1L])
    b <- as.double(Mvals[, 2L])
    c <- as.double(Mvals[, 3L])
    u <- as.double(Zvals[, 1L])
    v <- as.double(Zvals[, 2L])

    det <- a * c - b * b
    good <- is.finite(det) & (abs(det) > .Machine$double.eps)
    if (any(good)) {
      invdet <- 1 / det[good]
      beta1 <- (c[good] * u[good] - b[good] * v[good]) * invdet
      beta2 <- (a[good] * v[good] - b[good] * u[good]) * invdet
      out[good] <- rhs[1L] * beta1 + rhs[2L] * beta2
    }
  } else if (p == 3L) {
    a <- as.double(Mvals[, 1L])
    b <- as.double(Mvals[, 2L])
    c <- as.double(Mvals[, 3L])
    d <- as.double(Mvals[, 4L])
    e <- as.double(Mvals[, 5L])
    f <- as.double(Mvals[, 6L])
    u <- as.double(Zvals[, 1L])
    v <- as.double(Zvals[, 2L])
    w <- as.double(Zvals[, 3L])

    det <- a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d)
    good <- is.finite(det) & (abs(det) > .Machine$double.eps)
    if (any(good)) {
      c11 <- d[good] * f[good] - e[good] * e[good]
      c12 <- c[good] * e[good] - b[good] * f[good]
      c13 <- b[good] * e[good] - c[good] * d[good]
      c22 <- a[good] * f[good] - c[good] * c[good]
      c23 <- b[good] * c[good] - a[good] * e[good]
      c33 <- a[good] * d[good] - b[good] * b[good]
      invdet <- 1 / det[good]

      beta1 <- (c11 * u[good] + c12 * v[good] + c13 * w[good]) * invdet
      beta2 <- (c12 * u[good] + c22 * v[good] + c23 * w[good]) * invdet
      beta3 <- (c13 * u[good] + c23 * v[good] + c33 * w[good]) * invdet
      out[good] <- rhs[1L] * beta1 + rhs[2L] * beta2 + rhs[3L] * beta3
    }
  }

  bad <- which(!is.finite(out))
  if (length(bad)) {
    p <- ncol(Zvals)
    for (ii in bad) {
      A <- .np_inid_lp_unpack_sym_row(mrow = Mvals[ii, ], p = p)
      out[ii] <- .np_inid_lp_predict_row(
        A = A,
        z = as.double(Zvals[ii, ]),
        rhs = rhs,
        ridge.grid = ridge.grid
      )
    }
  }

  out
}

.np_counts_to_indices <- function(counts.col) {
  rep.int(seq_along(counts.col), as.integer(counts.col))
}

.np_inid_boot_from_regression_exact <- function(xdat,
                                                exdat,
                                                bws,
                                                ydat,
                                                B,
                                                counts = NULL,
                                                counts.drawer = NULL,
                                                gradients = FALSE,
                                                gradient.order = 1L,
                                                slice.index = 1L) {
  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  ydat <- as.double(ydat)
  B <- as.integer(B)

  n <- nrow(xdat)
  neval <- nrow(exdat)
  if (length(ydat) != n)
    stop("length of ydat must match training rows in exact regression bootstrap helper")
  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid exact regression bootstrap dimensions")

  use.local.direct <- isTRUE(identical(bws$type, "generalized_nn"))

  fit_one <- function(x.train, y.train) {
    fit <- .np_regression_direct(
      bws = bws,
      txdat = x.train,
      tydat = y.train,
      exdat = exdat,
      gradients = isTRUE(gradients),
      gradient.order = gradient.order,
      local.mode = use.local.direct
    )
    if (isTRUE(gradients))
      as.vector(fit$grad[, slice.index])
    else
      as.vector(fit$mean)
  }

  t0 <- fit_one(x.train = xdat, y.train = ydat)
  nout <- length(t0)
  chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
    B = B,
    chunk.size = .np_inid_chunk_size(n = n, B = B),
    comm = 1L,
    include.master = TRUE
  )
  tmat <- matrix(NA_real_, nrow = B, ncol = nout)

  compute_chunk <- function(counts.chunk) {
    counts.chunk <- as.matrix(counts.chunk)
    bsz <- ncol(counts.chunk)
    out <- matrix(NA_real_, nrow = bsz, ncol = nout)
    for (jj in seq_len(bsz)) {
      idx <- .np_counts_to_indices(counts.chunk[, jj])
      out[jj, ] <- fit_one(
        x.train = xdat[idx, , drop = FALSE],
        y.train = ydat[idx]
      )
    }
    out
  }

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    if (.npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-regression-exact-counts"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        start <- as.integer(task$start)
        stopi <- start + as.integer(task$bsz) - 1L
        compute_chunk(counts.chunk = counts.mat[, start:stopi, drop = FALSE])
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = nout,
        what = "inid-regression-exact-counts",
        profile.where = "mpi.applyLB:inid-regression-exact-counts",
        comm = 1L,
        required.bindings = list(
          counts.mat = counts.mat,
          compute_chunk = compute_chunk,
          .np_ksum_unconditional_eval_exact = .np_ksum_unconditional_eval_exact
        )
      )
    }
    if (anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-regression-exact-counts fan-out returned incomplete results",
        what = "inid-regression-exact-counts"
      )
  } else {
    if (!is.null(counts.drawer) &&
        .npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-regression-exact-block"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        start <- as.integer(task$start)
        stopi <- start + as.integer(task$bsz) - 1L
        counts.chunk <- .np_inid_counts_matrix(
          n = n,
          B = as.integer(task$bsz),
          counts = counts.drawer(start, stopi)
        )
        compute_chunk(counts.chunk = counts.chunk)
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = nout,
        what = "inid-regression-exact-block",
        profile.where = "mpi.applyLB:inid-regression-exact-block",
        comm = 1L,
        required.bindings = list(
          n = n,
          counts.drawer = counts.drawer,
          compute_chunk = compute_chunk,
          .np_ksum_unconditional_eval_exact = .np_ksum_unconditional_eval_exact
        )
      )
    }

    if (!is.null(counts.drawer) && anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-regression-exact-block fan-out returned incomplete results",
        what = "inid-regression-exact-block"
      )

    prob <- rep.int(1 / n, n)
    if (anyNA(tmat) &&
        .npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-regression-exact"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        set.seed(as.integer(task$seed))
        bsz <- as.integer(task$bsz)
        counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
        compute_chunk(counts.chunk = counts.chunk)
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = nout,
        what = "inid-regression-exact",
        profile.where = "mpi.applyLB:inid-regression-exact",
        comm = 1L,
        required.bindings = list(
          n = n,
          prob = prob,
          compute_chunk = compute_chunk,
          .np_ksum_unconditional_eval_exact = .np_ksum_unconditional_eval_exact
        )
      )
    }

    if (anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-regression-exact fan-out returned incomplete results",
        what = "inid-regression-exact"
      )
  }

  if (any(!is.finite(t0)) || any(!is.finite(tmat)))
    stop("inid regression exact helper path produced non-finite values")

  list(t = tmat, t0 = t0)
}

.np_inid_boot_from_regression <- function(xdat,
                                          exdat,
                                          bws,
                                          ydat,
                                          B,
                                          counts = NULL,
                                          counts.drawer = NULL,
                                          ridge = 1.0e-12,
                                          gradients = FALSE,
                                          gradient.order = 1L,
                                          slice.index = 1L) {
  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  ydat <- as.double(ydat)
  B <- as.integer(B)

  n <- nrow(xdat)
  neval <- nrow(exdat)
  if (length(ydat) != n)
    stop("length of ydat must match training rows")
  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid inid regression bootstrap dimensions")
  if (isTRUE(gradients) || !identical(bws$type, "fixed")) {
    return(.np_inid_boot_from_regression_exact(
      xdat = xdat,
      exdat = exdat,
      bws = bws,
      ydat = ydat,
      B = B,
      counts = counts,
      counts.drawer = counts.drawer,
      gradients = gradients,
      gradient.order = gradient.order,
      slice.index = slice.index
    ))
  }
  ridge.grid <- npRidgeSequenceFromBase(n.train = n, ridge.base = ridge, cap = 1.0)

  regtype <- if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  ncon <- bws$ncon

  if (isTRUE(gradients)) {
    fit0 <- .np_regression_direct(
      bws = bws,
      txdat = xdat,
      tydat = ydat,
      exdat = exdat,
      gradients = TRUE,
      gradient.order = gradient.order
    )
    t0 <- fit0$grad[, slice.index]

    nout <- length(t0)
    chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
      B = B,
      chunk.size = .np_inid_chunk_size(n = n, B = B),
      comm = 1L,
      include.master = TRUE
    )
    tmat <- matrix(NA_real_, nrow = B, ncol = nout)

    compute_grad_chunk <- function(counts.chunk) {
      counts.chunk <- as.matrix(counts.chunk)
      bsz <- ncol(counts.chunk)
      out <- matrix(NA_real_, nrow = bsz, ncol = nout)
      for (jj in seq_len(bsz)) {
        idx <- rep.int(seq_len(n), as.integer(counts.chunk[, jj]))
        fit.b <- .np_regression_direct(
          bws = bws,
          txdat = xdat[idx, , drop = FALSE],
          tydat = ydat[idx],
          exdat = exdat,
          gradients = TRUE,
          gradient.order = gradient.order,
          local.mode = use.local.direct
        )
        out[jj, ] <- fit.b$grad[, slice.index]
      }
      out
    }

    if (!is.null(counts)) {
      counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)

      if (.npRmpi_bootstrap_fanout_enabled(
            comm = 1L,
            n = n,
            B = B,
            chunk.size = chunk.size,
            what = "inid-regression-grad-counts"
          )) {
        tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
        worker <- function(task) {
          start <- as.integer(task$start)
          stopi <- start + as.integer(task$bsz) - 1L
          compute_grad_chunk(counts.chunk = counts.mat[, start:stopi, drop = FALSE])
        }
        tmat <- .npRmpi_bootstrap_run_fanout(
          tasks = tasks,
          worker = worker,
          ncol.out = nout,
          what = "inid-regression-grad-counts",
          profile.where = "mpi.applyLB:inid-regression-grad-counts",
          comm = 1L,
          required.bindings = list(
            counts.mat = counts.mat,
            compute_grad_chunk = compute_grad_chunk
          )
        )
      }

      if (anyNA(tmat)) {
        .npRmpi_bootstrap_fail_or_fallback(
          msg = "inid-regression-grad-counts fan-out returned incomplete results",
          what = "inid-regression-grad-counts"
        )
      }
    } else {
      if (!is.null(counts.drawer) &&
          .npRmpi_bootstrap_fanout_enabled(
            comm = 1L,
            n = n,
            B = B,
            chunk.size = chunk.size,
            what = "inid-regression-grad-block"
          )) {
        tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
        worker <- function(task) {
          start <- as.integer(task$start)
          stopi <- start + as.integer(task$bsz) - 1L
          counts.chunk <- .np_inid_counts_matrix(
            n = n,
            B = as.integer(task$bsz),
            counts = counts.drawer(start, stopi)
          )
          compute_grad_chunk(counts.chunk = counts.chunk)
        }
        tmat <- .npRmpi_bootstrap_run_fanout(
          tasks = tasks,
          worker = worker,
          ncol.out = nout,
          what = "inid-regression-grad-block",
          profile.where = "mpi.applyLB:inid-regression-grad-block",
          comm = 1L,
          required.bindings = list(
            n = n,
            counts.drawer = counts.drawer,
            compute_grad_chunk = compute_grad_chunk
          )
        )
      }

      if (!is.null(counts.drawer) && anyNA(tmat)) {
        .npRmpi_bootstrap_fail_or_fallback(
          msg = "inid-regression-grad-block fan-out returned incomplete results",
          what = "inid-regression-grad-block"
        )
      }

      prob <- rep.int(1 / n, n)
      if (anyNA(tmat) &&
          .npRmpi_bootstrap_fanout_enabled(
            comm = 1L,
            n = n,
            B = B,
            chunk.size = chunk.size,
            what = "inid-regression-grad"
          )) {
        tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
        worker <- function(task) {
          set.seed(as.integer(task$seed))
          bsz <- as.integer(task$bsz)
          counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
          compute_grad_chunk(counts.chunk = counts.chunk)
        }
        tmat <- .npRmpi_bootstrap_run_fanout(
          tasks = tasks,
          worker = worker,
          ncol.out = nout,
          what = "inid-regression-grad",
          profile.where = "mpi.applyLB:inid-regression-grad",
          comm = 1L,
          required.bindings = list(
            n = n,
            prob = prob,
            compute_grad_chunk = compute_grad_chunk
          )
        )
      }

      if (anyNA(tmat)) {
        .npRmpi_bootstrap_fail_or_fallback(
          msg = "inid-regression-grad fan-out returned incomplete results",
          what = "inid-regression-grad"
        )
      }
    }

    if (any(!is.finite(t0)) || any(!is.finite(tmat)))
      stop("inid regression helper gradient path produced non-finite values")

    return(list(t = tmat, t0 = as.vector(t0)))
  }

  if (identical(regtype, "lc")) {
    H <- suppressWarnings(
      tryCatch(
        npreghat.rbandwidth(
          bws = bws,
          txdat = xdat,
          exdat = exdat,
          s = 0L,
          output = "matrix"
        ),
        error = function(e) NULL
      )
    )
    if (!is.null(H)) {
      if (!is.matrix(H))
        H <- matrix(as.double(H), nrow = neval, ncol = n)
      if (nrow(H) == neval && ncol(H) == n) {
        return(.np_inid_lc_boot_from_hat(
          H = H,
          ydat = ydat,
          B = B,
          counts = counts,
          counts.drawer = counts.drawer
        ))
      }
    }
  }

  degree <- if (identical(regtype, "lc")) {
    rep.int(0L, ncon)
  } else if (identical(regtype, "ll")) {
    rep.int(1L, ncon)
  } else {
    npValidateGlpDegree(
      regtype = "lp",
      degree = bws$degree,
      ncon = ncon
    )
  }

  basis <- npValidateLpBasis(
    regtype = "lp",
    basis = if (is.null(bws$basis)) "glp" else bws$basis
  )
  bernstein.basis <- npValidateGlpBernstein(
    regtype = "lp",
    bernstein.basis = isTRUE(bws$bernstein.basis)
  )

  kw <- .np_plot_kernel_weights_direct(
    bws = bws,
    txdat = xdat,
    exdat = exdat,
    operator = "normal",
    where = "direct regression kernel weights"
  )
  if (nrow(kw) != n || ncol(kw) != neval)
    stop("kernel-weight matrix shape mismatch")

  W <- W.lp(
    xdat = xdat,
    degree = degree,
    basis = basis,
    bernstein.basis = bernstein.basis
  )
  W.eval <- W.lp(
    xdat = xdat,
    exdat = exdat,
    degree = degree,
    basis = basis,
    bernstein.basis = bernstein.basis
  )
  W <- as.matrix(W)
  W.eval <- as.matrix(W.eval)

  if (nrow(W) != n || nrow(W.eval) != neval || ncol(W.eval) != ncol(W))
    stop("regression moment design matrix shape mismatch")

  p <- ncol(W)
  mcols <- p * (p + 1L) / 2L
  rhs <- W.eval
  ones <- matrix(1.0, nrow = n, ncol = 1L)

  Mfeat <- vector("list", neval)
  Zfeat <- vector("list", neval)
  t0 <- numeric(neval)

  for (i in seq_len(neval)) {
    k <- as.double(kw[, i])
    WK <- W * k
    Zfeat[[i]] <- WK * ydat

    mf <- matrix(0.0, nrow = n, ncol = mcols)
    idx <- 1L
    for (a in seq_len(p)) {
      for (b in a:p) {
        mf[, idx] <- WK[, a] * W[, b]
        idx <- idx + 1L
      }
    }
    Mfeat[[i]] <- mf

    M0 <- crossprod(ones, mf)
    Z0 <- crossprod(ones, Zfeat[[i]])
    t0[i] <- if (p > 3L) {
      .np_inid_lp_predict_chunk_general(
        Mvals = M0,
        Zvals = Z0,
        rhs = rhs[i, ],
        ridge.grid = ridge.grid
      )[1L]
    } else {
      .np_inid_lp_predict_chunk(
        Mvals = M0,
        Zvals = Z0,
        rhs = rhs[i, ],
        ridge.grid = ridge.grid
      )[1L]
    }
  }

  compute_chunk <- function(counts.chunk) {
    counts.chunk <- as.matrix(counts.chunk)
    bsz <- ncol(counts.chunk)
    out <- matrix(NA_real_, nrow = bsz, ncol = neval)
    for (i in seq_len(neval)) {
      Mvals <- crossprod(counts.chunk, Mfeat[[i]])
      Zvals <- crossprod(counts.chunk, Zfeat[[i]])
      out[, i] <- if (p > 3L) {
        .np_inid_lp_predict_chunk_general(
          Mvals = Mvals,
          Zvals = Zvals,
          rhs = rhs[i, ],
          ridge.grid = ridge.grid
        )
      } else {
        .np_inid_lp_predict_chunk(
          Mvals = Mvals,
          Zvals = Zvals,
          rhs = rhs[i, ],
          ridge.grid = ridge.grid
        )
      }
    }
    out
  }

  chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
    B = B,
    chunk.size = .np_inid_chunk_size(n = n, B = B),
    comm = 1L,
    include.master = TRUE
  )
  tmat <- matrix(NA_real_, nrow = B, ncol = neval)

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)

    if (.npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-regression-counts"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        start <- as.integer(task$start)
        stopi <- start + as.integer(task$bsz) - 1L
        compute_chunk(counts.chunk = counts.mat[, start:stopi, drop = FALSE])
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = neval,
        what = "inid-regression-counts",
        profile.where = "mpi.applyLB:inid-regression-counts",
        comm = 1L,
        required.bindings = list(
          counts.mat = counts.mat,
          compute_chunk = compute_chunk,
          .np_ksum_conditional_eval_exact = .np_ksum_conditional_eval_exact
        )
      )
    }

    if (anyNA(tmat)) {
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-regression-counts fan-out returned incomplete results",
        what = "inid-regression-counts"
      )
    }
  } else {
    if (!is.null(counts.drawer) &&
        .npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-regression-block"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        start <- as.integer(task$start)
        stopi <- start + as.integer(task$bsz) - 1L
        counts.chunk <- .np_inid_counts_matrix(
          n = n,
          B = as.integer(task$bsz),
          counts = counts.drawer(start, stopi)
        )
        compute_chunk(counts.chunk = counts.chunk)
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = neval,
        what = "inid-regression-block",
        profile.where = "mpi.applyLB:inid-regression-block",
        comm = 1L,
        required.bindings = list(
          n = n,
          counts.drawer = counts.drawer,
          compute_chunk = compute_chunk,
          .np_ksum_conditional_eval_exact = .np_ksum_conditional_eval_exact
        )
      )
    }

    if (!is.null(counts.drawer) && anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-regression-block fan-out returned incomplete results",
        what = "inid-regression-block"
      )

    prob <- rep.int(1 / n, n)

    if (anyNA(tmat) &&
        .npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-regression"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        set.seed(as.integer(task$seed))
        bsz <- as.integer(task$bsz)
        counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
        compute_chunk(counts.chunk = counts.chunk)
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = neval,
        what = "inid-regression",
        profile.where = "mpi.applyLB:inid-regression",
        comm = 1L,
        required.bindings = list(
          n = n,
          prob = prob,
          compute_chunk = compute_chunk,
          .np_ksum_conditional_eval_exact = .np_ksum_conditional_eval_exact
        )
      )
    }

    if (anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-regression fan-out returned incomplete results",
        what = "inid-regression"
      )
  }

  if (any(!is.finite(t0)) || any(!is.finite(tmat)))
    stop("inid regression helper path produced non-finite values")

  list(t = tmat, t0 = t0)
}

.np_inid_scoef_numeric_y <- function(ydat, bws) {
  if (is.factor(ydat)) {
    if (is.null(bws$ydati))
      stop("factor response requires bws$ydati for smooth coefficient inid helper")
    yadj <- adjustLevels(data.frame(ydat), bws$ydati)
    return((bws$ydati$all.dlev[[1L]])[as.integer(yadj[, 1L])])
  }
  as.double(ydat)
}

.np_inid_scoef_predict_row <- function(mrow, zrow, rhs, ridge.grid) {
  A <- .np_inid_lp_unpack_sym_row(mrow = mrow, p = length(rhs))
  tyw <- as.double(zrow)
  nc <- ncol(A)
  ridge.grid <- as.double(ridge.grid)
  if (!length(ridge.grid) || anyNA(ridge.grid) || any(!is.finite(ridge.grid)) || any(ridge.grid < 0))
    stop("argument 'ridge.grid' must be a non-empty non-negative finite numeric vector")

  maxPenalty <- sqrt(.Machine$double.xmax)
  coef.ii <- rep(maxPenalty, nc)
  for (ridge in ridge.grid) {
    ridge.val <- ridge * tyw[1L] / NZD(A[1L, 1L])
    coef.try <- tryCatch(
      solve(
        A + diag(rep(ridge, nc)),
        tyw + c(ridge.val, rep(0, nc - 1L))
      ),
      error = function(e) e
    )
    if (!(inherits(coef.try, "error") || any(!is.finite(coef.try))))
      return(sum(as.double(rhs) * as.double(coef.try)))
  }

  sum(as.double(rhs) * coef.ii)
}

.np_inid_scoef_predict_chunk <- function(Mvals, Zvals, rhs) {
  Mvals <- as.matrix(Mvals)
  Zvals <- as.matrix(Zvals)
  rhs <- as.double(rhs)

  bsz <- nrow(Mvals)
  p <- ncol(Zvals)
  out <- rep(NA_real_, bsz)

  if (p == 1L) {
    den <- as.double(Mvals[, 1L])
    good <- is.finite(den) & (abs(den) > .Machine$double.eps)
    out[good] <- rhs[1L] * as.double(Zvals[good, 1L]) / den[good]
    return(out)
  }

  if (p == 2L) {
    a <- as.double(Mvals[, 1L])
    b <- as.double(Mvals[, 2L])
    c <- as.double(Mvals[, 3L])
    u <- as.double(Zvals[, 1L])
    v <- as.double(Zvals[, 2L])
    det <- a * c - b * b
    good <- is.finite(det) & (abs(det) > .Machine$double.eps)
    if (any(good)) {
      invdet <- 1 / det[good]
      beta1 <- (c[good] * u[good] - b[good] * v[good]) * invdet
      beta2 <- (a[good] * v[good] - b[good] * u[good]) * invdet
      out[good] <- rhs[1L] * beta1 + rhs[2L] * beta2
    }
    return(out)
  }

  if (p == 3L) {
    a <- as.double(Mvals[, 1L])
    b <- as.double(Mvals[, 2L])
    c <- as.double(Mvals[, 3L])
    d <- as.double(Mvals[, 4L])
    e <- as.double(Mvals[, 5L])
    f <- as.double(Mvals[, 6L])
    u <- as.double(Zvals[, 1L])
    v <- as.double(Zvals[, 2L])
    w <- as.double(Zvals[, 3L])

    det <- a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d)
    good <- is.finite(det) & (abs(det) > .Machine$double.eps)
    if (any(good)) {
      c11 <- d[good] * f[good] - e[good] * e[good]
      c12 <- c[good] * e[good] - b[good] * f[good]
      c13 <- b[good] * e[good] - c[good] * d[good]
      c22 <- a[good] * f[good] - c[good] * c[good]
      c23 <- b[good] * c[good] - a[good] * e[good]
      c33 <- a[good] * d[good] - b[good] * b[good]
      invdet <- 1 / det[good]

      beta1 <- (c11 * u[good] + c12 * v[good] + c13 * w[good]) * invdet
      beta2 <- (c12 * u[good] + c22 * v[good] + c23 * w[good]) * invdet
      beta3 <- (c13 * u[good] + c23 * v[good] + c33 * w[good]) * invdet
      out[good] <- rhs[1L] * beta1 + rhs[2L] * beta2 + rhs[3L] * beta3
    }
    return(out)
  }

  out
}

.np_inid_boot_from_scoef <- function(txdat,
                                     ydat,
                                     tzdat,
                                     exdat,
                                     ezdat,
                                     bws,
                                     B,
                                     counts = NULL,
                                     counts.drawer = NULL,
                                     leave.one.out = FALSE) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)
  B <- as.integer(B)

  miss.z <- missing(tzdat) || is.null(tzdat)
  if (miss.z) {
    tzdat <- txdat
    ezdat <- exdat
  } else {
    tzdat <- toFrame(tzdat)
    ezdat <- toFrame(ezdat)
  }

  if (nrow(txdat) != nrow(tzdat))
    stop("smooth coefficient inid helper requires aligned txdat/tzdat rows")
  if (nrow(exdat) != nrow(ezdat))
    stop("smooth coefficient inid helper requires aligned exdat/ezdat rows")
  if (ncol(txdat) != ncol(exdat))
    stop("smooth coefficient inid helper requires matching txdat/exdat columns")
  if (nrow(txdat) < 1L || nrow(exdat) < 1L || B < 1L)
    stop("invalid smooth coefficient inid helper dimensions")

  if (length(ydat) != nrow(txdat))
    stop("length of ydat must match training rows in smooth coefficient inid helper")

  txdat <- adjustLevels(txdat, bws$xdati)
  exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
  if (!miss.z) {
    tzdat <- adjustLevels(tzdat, bws$zdati)
    ezdat <- adjustLevels(ezdat, bws$zdati, allowNewCells = TRUE)
  }

  y.num <- .np_inid_scoef_numeric_y(ydat = ydat, bws = bws)
  X.train <- toMatrix(txdat)
  X.eval <- toMatrix(exdat)
  W.train <- as.matrix(data.frame(1, X.train))
  W.eval <- as.matrix(data.frame(1, X.eval))

  kw <- .np_plot_kernel_weights_direct(
    bws = bws,
    txdat = tzdat,
    exdat = ezdat,
    operator = "normal",
    where = "direct smooth-coefficient kernel weights"
  )

  n <- nrow(W.train)
  neval <- nrow(W.eval)
  if (nrow(kw) != n || ncol(kw) != neval)
    stop("smooth coefficient inid helper kernel-weight matrix shape mismatch")

  p <- ncol(W.train)
  mcols <- p * (p + 1L) / 2L
  ones <- matrix(1.0, nrow = n, ncol = 1L)
  ridge.grid <- npRidgeSequenceAdditive(n.train = n, cap = 1.0)

  Mfeat <- vector("list", neval)
  Zfeat <- vector("list", neval)
  t0 <- numeric(neval)

  for (i in seq_len(neval)) {
    k <- as.double(kw[, i])
    WK <- W.train * k
    zf <- WK * y.num

    mf <- matrix(0.0, nrow = n, ncol = mcols)
    idx <- 1L
    for (a in seq_len(p)) {
      for (b in a:p) {
        mf[, idx] <- WK[, a] * W.train[, b]
        idx <- idx + 1L
      }
    }

    Mfeat[[i]] <- mf
    Zfeat[[i]] <- zf

    M0 <- crossprod(ones, mf)
    Z0 <- crossprod(ones, zf)
    t0i <- .np_inid_scoef_predict_chunk(Mvals = M0, Zvals = Z0, rhs = W.eval[i, ])[1L]
    if (!is.finite(t0i)) {
      t0i <- .np_inid_scoef_predict_row(
        mrow = M0[1L, ],
        zrow = Z0[1L, ],
        rhs = W.eval[i, ],
        ridge.grid = ridge.grid
      )
    }
    t0[i] <- t0i
  }

  compute_chunk <- function(counts.chunk) {
    counts.chunk <- as.matrix(counts.chunk)
    bsz <- ncol(counts.chunk)
    out.chunk <- matrix(NA_real_, nrow = bsz, ncol = neval)
    for (i in seq_len(neval)) {
      Mvals <- crossprod(counts.chunk, Mfeat[[i]])
      Zvals <- crossprod(counts.chunk, Zfeat[[i]])
      if (bsz == 1L) {
        Mvals <- matrix(Mvals, nrow = 1L)
        Zvals <- matrix(Zvals, nrow = 1L)
      }
      out <- .np_inid_scoef_predict_chunk(Mvals = Mvals, Zvals = Zvals, rhs = W.eval[i, ])
      bad <- which(!is.finite(out))
      if (length(bad)) {
        for (bb in bad) {
          out[bb] <- .np_inid_scoef_predict_row(
            mrow = Mvals[bb, ],
            zrow = Zvals[bb, ],
            rhs = W.eval[i, ],
            ridge.grid = ridge.grid
          )
        }
      }
      out.chunk[, i] <- out
    }
    out.chunk
  }

  chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
    B = B,
    chunk.size = .np_inid_chunk_size(n = n, B = B),
    comm = 1L,
    include.master = TRUE
  )
  tmat <- matrix(NA_real_, nrow = B, ncol = neval)

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)

    if (.npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-scoef-counts"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        start <- as.integer(task$start)
        stopi <- start + as.integer(task$bsz) - 1L
        compute_chunk(counts.chunk = counts.mat[, start:stopi, drop = FALSE])
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = neval,
        what = "inid-scoef-counts",
        profile.where = "mpi.applyLB:inid-scoef-counts",
        comm = 1L,
        required.bindings = list(
          counts.mat = counts.mat,
          compute_chunk = compute_chunk
        )
      )
    }

    if (anyNA(tmat)) {
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-scoef-counts fan-out returned incomplete results",
        what = "inid-scoef-counts"
      )
    }
  } else {
    if (!is.null(counts.drawer) &&
        .npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-scoef-block"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        start <- as.integer(task$start)
        stopi <- start + as.integer(task$bsz) - 1L
        counts.chunk <- .np_inid_counts_matrix(
          n = n,
          B = as.integer(task$bsz),
          counts = counts.drawer(start, stopi)
        )
        compute_chunk(counts.chunk = counts.chunk)
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = neval,
        what = "inid-scoef-block",
        profile.where = "mpi.applyLB:inid-scoef-block",
        comm = 1L,
        required.bindings = list(
          n = n,
          counts.drawer = counts.drawer,
          compute_chunk = compute_chunk
        )
      )
    }

    if (!is.null(counts.drawer) && anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-scoef-block fan-out returned incomplete results",
        what = "inid-scoef-block"
      )

    prob <- rep.int(1 / n, n)

    if (anyNA(tmat) &&
        .npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-scoef"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        set.seed(as.integer(task$seed))
        bsz <- as.integer(task$bsz)
        counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
        compute_chunk(counts.chunk = counts.chunk)
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = neval,
        what = "inid-scoef",
        profile.where = "mpi.applyLB:inid-scoef",
        comm = 1L,
        required.bindings = list(
          n = n,
          prob = prob,
          compute_chunk = compute_chunk
        )
      )
    }

    if (anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-scoef fan-out returned incomplete results",
        what = "inid-scoef"
      )
  }

  if (any(!is.finite(t0)) || any(!is.finite(tmat)))
    stop("inid smooth coefficient helper produced non-finite values")

  list(t = tmat, t0 = t0)
}

.np_plreg_numeric_x_matrix <- function(txdat, exdat, bws) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)

  p <- ncol(txdat)
  if (p < 1L)
    stop("plreg inid helper path requires at least one linear regressor")
  if (ncol(exdat) != p)
    stop("training/evaluation linear regressor dimensions do not match")

  x.train.num <- matrix(0.0, nrow = nrow(txdat), ncol = p)
  x.eval.num <- matrix(0.0, nrow = nrow(exdat), ncol = p)

  for (j in seq_len(p)) {
    if (is.factor(txdat[[j]])) {
      trj <- adjustLevels(txdat[, j, drop = FALSE], bws$bw[[j + 1L]]$ydati)
      evj <- adjustLevels(exdat[, j, drop = FALSE], bws$bw[[j + 1L]]$ydati, allowNewCells = TRUE)
      lev <- bws$bw[[j + 1L]]$ydati$all.dlev[[1L]]
      x.train.num[, j] <- lev[as.integer(trj[, 1L])]
      x.eval.num[, j] <- lev[as.integer(evj[, 1L])]
    } else {
      x.train.num[, j] <- as.double(txdat[[j]])
      x.eval.num[, j] <- as.double(exdat[[j]])
    }
  }

  list(train = x.train.num, eval = x.eval.num)
}

.np_plreg_weighted_coef <- function(X, y, w, ridge = 1.0e-12) {
  X <- as.matrix(X)
  y <- as.double(y)
  w <- as.double(w)
  if (nrow(X) != length(y) || length(w) != length(y))
    stop("weighted plreg solve dimension mismatch")

  w <- pmax(w, 0.0)

  XtWX <- crossprod(X, X * w)
  XtWy <- drop(crossprod(X, y * w))
  XtWy0 <- XtWy[1L]
  ridge.grid <- npRidgeSequenceFromBase(
    n.train = nrow(X),
    ridge.base = max(0.0, as.double(ridge)),
    cap = 1.0
  )

  for (ridge.try in ridge.grid) {
    A <- XtWX
    ridge <- ridge.try
    z <- XtWy
    if (ridge > 0)
      diag(A) <- diag(A) + ridge
    if (ridge > 0)
      z[1L] <- XtWy0 + ridge * XtWy0 / NZD(A[1L, 1L])
    beta <- tryCatch(
      drop(solve(A, matrix(z, ncol = 1L))),
      error = function(e) NULL
    )
    if (!is.null(beta) && all(is.finite(beta)))
      return(as.double(beta))
  }

  stop("plreg weighted solve failed")
}

.np_inid_boot_from_plreg <- function(txdat,
                                     ydat,
                                     tzdat,
                                     exdat,
                                     ezdat,
                                     bws,
                                     B,
                                     counts = NULL,
                                     counts.drawer = NULL,
                                     ridge = 1.0e-12) {
  txdat <- toFrame(txdat)
  tzdat <- toFrame(tzdat)
  exdat <- toFrame(exdat)
  ezdat <- toFrame(ezdat)
  B <- as.integer(B)

  n <- nrow(txdat)
  neval <- nrow(exdat)
  p <- ncol(txdat)
  if (nrow(tzdat) != n)
    stop("plreg inid helper path requires aligned txdat/tzdat rows")
  if (nrow(ezdat) != neval)
    stop("plreg inid helper path requires aligned exdat/ezdat rows")
  if (n < 1L || neval < 1L || p < 1L || B < 1L)
    stop("invalid plreg inid helper path dimensions")

  y.num <- if (is.factor(ydat)) {
    ty <- adjustLevels(data.frame(ydat), bws$bw$yzbw$ydati)
    bws$bw$yzbw$ydati$all.dlev[[1L]][as.integer(ty[, 1L])]
  } else {
    as.double(ydat)
  }
  if (length(y.num) != n)
    stop("length of ydat must match training rows")

  x.num <- .np_plreg_numeric_x_matrix(txdat = txdat, exdat = exdat, bws = bws)
  x.train.num <- x.num$train
  x.eval.num <- x.num$eval

  counts.mat <- if (!is.null(counts)) {
    .np_inid_counts_matrix(n = n, B = B, counts = counts)
  } else if (!is.null(counts.drawer)) {
    .np_inid_counts_matrix(n = n, B = B, counts = counts.drawer(1L, B))
  } else {
    .np_inid_counts_matrix(n = n, B = B)
  }

  y.train <- .np_inid_boot_from_regression(
    xdat = tzdat,
    exdat = tzdat,
    bws = bws$bw$yzbw,
    ydat = y.num,
    B = B,
    counts = counts.mat,
    ridge = ridge
  )
  y.eval <- .np_inid_boot_from_regression(
    xdat = tzdat,
    exdat = ezdat,
    bws = bws$bw$yzbw,
    ydat = y.num,
    B = B,
    counts = counts.mat,
    ridge = ridge
  )

  x.train <- vector("list", p)
  x.eval <- vector("list", p)
  for (j in seq_len(p)) {
    x.train[[j]] <- .np_inid_boot_from_regression(
      xdat = tzdat,
      exdat = tzdat,
      bws = bws$bw[[j + 1L]],
      ydat = x.train.num[, j],
      B = B,
      counts = counts.mat,
      ridge = ridge
    )
    x.eval[[j]] <- .np_inid_boot_from_regression(
      xdat = tzdat,
      exdat = ezdat,
      bws = bws$bw[[j + 1L]],
      ydat = x.train.num[, j],
      B = B,
      counts = counts.mat,
      ridge = ridge
    )
  }

  xres.train0 <- matrix(0.0, nrow = n, ncol = p)
  xres.eval0 <- matrix(0.0, nrow = neval, ncol = p)
  for (j in seq_len(p)) {
    xres.train0[, j] <- x.train.num[, j] - as.double(x.train[[j]]$t0)
    xres.eval0[, j] <- x.eval.num[, j] - as.double(x.eval[[j]]$t0)
  }
  yres0 <- y.num - as.double(y.train$t0)
  beta0 <- .np_plreg_weighted_coef(
    X = xres.train0,
    y = yres0,
    w = rep.int(1.0, n),
    ridge = ridge
  )
  t0 <- as.double(y.eval$t0) + as.vector(xres.eval0 %*% beta0)

  tmat <- matrix(NA_real_, nrow = B, ncol = neval)
  xres.train.b <- matrix(0.0, nrow = n, ncol = p)
  xres.eval.b <- matrix(0.0, nrow = neval, ncol = p)

  for (b in seq_len(B)) {
    for (j in seq_len(p)) {
      xres.train.b[, j] <- x.train.num[, j] - x.train[[j]]$t[b, ]
      xres.eval.b[, j] <- x.eval.num[, j] - x.eval[[j]]$t[b, ]
    }
    yres.b <- y.num - y.train$t[b, ]
    beta.b <- .np_plreg_weighted_coef(
      X = xres.train.b,
      y = yres.b,
      w = counts.mat[, b],
      ridge = ridge
    )
    tmat[b, ] <- y.eval$t[b, ] + as.vector(xres.eval.b %*% beta.b)
  }

  if (any(!is.finite(t0)) || any(!is.finite(tmat)))
    stop("plreg inid helper path produced non-finite values")

  list(t = tmat, t0 = t0)
}

.np_boot_matrix_from_ksum <- function(ksum, B, nout, where = "ksum helper path") {
  if (is.null(dim(ksum))) {
    if (B == 1L && length(ksum) == nout)
      return(matrix(as.double(ksum), nrow = 1L))
    stop(sprintf("%s returned unexpected vector shape", where))
  }

  km <- as.matrix(ksum)
  if (nrow(km) == B && ncol(km) == nout)
    return(km)
  if (nrow(km) == nout && ncol(km) == B)
    return(t(km))
  if (B == 1L && length(km) == nout)
    return(matrix(as.double(km), nrow = 1L))

  stop(sprintf("%s returned unexpected matrix shape", where))
}

.np_ksum_unconditional_eval_exact <- function(xdat, exdat, bws, operator) {
  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  kw <- .np_plot_kernel_weights_direct(
    bws = bws,
    txdat = xdat,
    exdat = exdat,
    operator = operator,
    where = "direct unconditional exact kernel weights"
  )
  colSums(kw) / nrow(xdat)
}

.np_ksum_kernel_weights_matrix <- function(kernel.weights, ntrain, neval, where = "ksum helper path") {
  if (is.null(kernel.weights))
    stop(sprintf("%s did not return kernel weights", where))

  kw <- as.matrix(kernel.weights)
  if (nrow(kw) == ntrain && ncol(kw) == neval)
    return(kw)
  if (nrow(kw) == neval && ncol(kw) == ntrain)
    return(t(kw))
  if (length(kw) == (ntrain * neval))
    return(matrix(as.double(kw), nrow = ntrain, ncol = neval))

  stop(sprintf("%s returned kernel weights with unexpected shape", where))
}

.np_ksum_extract_kernel_weights <- function(npksum.out, where = "ksum helper path") {
  if (is.null(npksum.out) || !is.list(npksum.out))
    stop(sprintf("%s did not return a valid npksum object", where))

  kw <- npksum.out$kw
  if (is.null(kw))
    kw <- npksum.out$kernel.weights
  if (is.null(kw))
    stop(sprintf("%s did not return kernel weights", where))
  kw
}

.np_plot_kernel_weights_direct <- function(bws,
                                           txdat,
                                           exdat,
                                           operator,
                                           where = "plot kernel weights direct") {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)

  if (!(txdat %~% exdat))
    stop(sprintf("%s requires similar training/evaluation frames", where))

  if (!isa(bws, "kbandwidth"))
    bws <- kbandwidth(bws)

  if (length(bws$bw) != ncol(txdat))
    stop(sprintf("%s bandwidth length/data width mismatch", where))

  txdat <- adjustLevels(txdat, bws$xdati, allowNewCells = TRUE)
  exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
  npKernelBoundsCheckEval(exdat, bws$icon, bws$ckerlb, bws$ckerub, argprefix = "cker")

  txm <- toMatrix(txdat)
  exm <- toMatrix(exdat)
  tuno <- txm[, bws$iuno, drop = FALSE]
  tcon <- txm[, bws$icon, drop = FALSE]
  tord <- txm[, bws$iord, drop = FALSE]
  euno <- exm[, bws$iuno, drop = FALSE]
  econ <- exm[, bws$icon, drop = FALSE]
  eord <- exm[, bws$iord, drop = FALSE]

  operator <- as.character(operator)
  if (length(operator) == 1L)
    operator <- rep.int(operator, ncol(txdat))
  if (length(operator) != ncol(txdat))
    stop(sprintf("%s operator length must equal number of columns", where))
  op.int <- as.integer(ALL_OPERATORS[operator])
  if (anyNA(op.int))
    stop(sprintf("%s encountered unsupported operator(s): %s",
                 where, paste(unique(operator[is.na(op.int)]), collapse = ", ")))

  tnrow <- nrow(txdat)
  enrow <- nrow(exdat)
  nkw <- tnrow * enrow

  myopti <- list(
    num_obs_train = tnrow,
    num_obs_eval = enrow,
    num_uno = bws$nuno,
    num_ord = bws$nord,
    num_con = bws$ncon,
    int_LARGE_SF = SF_ARB,
    BANDWIDTH_reg_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN
    ),
    int_MINIMIZE_IO = if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
    kerneval = switch(bws$ckertype,
      gaussian = CKER_GAUSS + bws$ckerorder / 2 - 1,
      epanechnikov = CKER_EPAN + bws$ckerorder / 2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS
    ),
    ukerneval = switch(bws$ukertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR
    ),
    okerneval = switch(bws$okertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_LR,
      nliracine = OKER_NLR,
      racineliyan = OKER_RLY
    ),
    miss.ex = FALSE,
    leave.one.out = FALSE,
    bandwidth.divide = TRUE,
    mcv.numRow = attr(bws$xmcv, "num.row"),
    wncol = 0L,
    yncol = 0L,
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
    return.kernel.weights = TRUE,
    permutation.operator = PERMUTATION_OPERATORS[["none"]],
    compute.score = FALSE,
    compute.ocg = FALSE,
    suppress.parallel = TRUE
  )

  cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])
  asDouble <- function(data) if (is.null(data)) as.double(0.0) else as.double(data)

  out <- .Call(
    "C_np_kernelsum",
    asDouble(tuno), asDouble(tord), asDouble(tcon),
    as.double(0.0), as.double(0.0),
    asDouble(euno), asDouble(eord), asDouble(econ),
    as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
    as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
    as.integer(op.int),
    as.integer(myopti), as.double(1.0),
    as.integer(enrow),
    as.integer(0L),
    as.integer(nkw),
    as.double(cker.bounds.c$lb),
    as.double(cker.bounds.c$ub),
    PACKAGE = "npRmpi"
  )

  kw <- out[["kernel.weights"]]
  if (is.null(kw))
    kw <- out[["kw"]]
  .np_ksum_kernel_weights_matrix(
    kernel.weights = kw,
    ntrain = tnrow,
    neval = enrow,
    where = where
  )
}

.np_inid_boot_from_ksum_unconditional_exact <- function(xdat,
                                                        exdat,
                                                        bws,
                                                        B,
                                                        operator,
                                                        counts = NULL,
                                                        counts.drawer = NULL) {
  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  B <- as.integer(B)
  n <- nrow(xdat)
  neval <- nrow(exdat)
  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid unconditional exact bootstrap dimensions")

  fit_one <- function(x.train) {
    .np_ksum_unconditional_eval_exact(
      xdat = x.train,
      exdat = exdat,
      bws = bws,
      operator = operator
    )
  }

  t0 <- fit_one(x.train = xdat)
  nout <- length(t0)
  chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
    B = B,
    chunk.size = .np_inid_chunk_size(n = n, B = B),
    comm = 1L,
    include.master = TRUE
  )
  tmat <- matrix(NA_real_, nrow = B, ncol = nout)

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    if (.npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-ksum-unconditional-exact-counts"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        start <- as.integer(task$start)
        stopi <- start + as.integer(task$bsz) - 1L
        counts.chunk <- as.matrix(counts.mat[, start:stopi, drop = FALSE])
        bsz <- ncol(counts.chunk)
        out <- matrix(NA_real_, nrow = bsz, ncol = nout)
        for (jj in seq_len(bsz)) {
          idx <- rep.int(seq_len(n), counts.chunk[, jj])
          out[jj, ] <- .np_ksum_unconditional_eval_exact(
            xdat = xdat[idx, , drop = FALSE],
            exdat = exdat,
            bws = bws,
            operator = operator
          )
        }
        out
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = nout,
        what = "inid-ksum-unconditional-exact-counts",
        profile.where = "mpi.applyLB:inid-ksum-unconditional-exact-counts",
        comm = 1L,
        required.bindings = list(
          counts.mat = counts.mat,
          xdat = xdat,
          exdat = exdat,
          bws = bws,
          operator = operator,
          n = n,
          nout = nout,
          .np_ksum_unconditional_eval_exact = .np_ksum_unconditional_eval_exact
        )
      )
    }
    if (anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-ksum-unconditional-exact-counts fan-out returned incomplete results",
        what = "inid-ksum-unconditional-exact-counts"
      )
  } else {
    if (!is.null(counts.drawer) &&
        .npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-ksum-unconditional-exact-block"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        start <- as.integer(task$start)
        stopi <- start + as.integer(task$bsz) - 1L
        counts.chunk <- .np_inid_counts_matrix(
          n = n,
          B = as.integer(task$bsz),
          counts = counts.drawer(start, stopi)
        )
        bsz <- ncol(counts.chunk)
        out <- matrix(NA_real_, nrow = bsz, ncol = nout)
        for (jj in seq_len(bsz)) {
          idx <- rep.int(seq_len(n), counts.chunk[, jj])
          out[jj, ] <- .np_ksum_unconditional_eval_exact(
            xdat = xdat[idx, , drop = FALSE],
            exdat = exdat,
            bws = bws,
            operator = operator
          )
        }
        out
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = nout,
        what = "inid-ksum-unconditional-exact-block",
        profile.where = "mpi.applyLB:inid-ksum-unconditional-exact-block",
        comm = 1L,
        required.bindings = list(
          n = n,
          counts.drawer = counts.drawer,
          xdat = xdat,
          exdat = exdat,
          bws = bws,
          operator = operator,
          nout = nout,
          .np_ksum_unconditional_eval_exact = .np_ksum_unconditional_eval_exact
        )
      )
    }

    if (!is.null(counts.drawer) && anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-ksum-unconditional-exact-block fan-out returned incomplete results",
        what = "inid-ksum-unconditional-exact-block"
      )

    prob <- rep.int(1 / n, n)
    if (anyNA(tmat) &&
        .npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-ksum-unconditional-exact"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        set.seed(as.integer(task$seed))
        bsz <- as.integer(task$bsz)
        counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
        out <- matrix(NA_real_, nrow = bsz, ncol = nout)
        for (jj in seq_len(bsz)) {
          idx <- rep.int(seq_len(n), counts.chunk[, jj])
          out[jj, ] <- .np_ksum_unconditional_eval_exact(
            xdat = xdat[idx, , drop = FALSE],
            exdat = exdat,
            bws = bws,
            operator = operator
          )
        }
        out
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = nout,
        what = "inid-ksum-unconditional-exact",
        profile.where = "mpi.applyLB:inid-ksum-unconditional-exact",
        comm = 1L,
        required.bindings = list(
          n = n,
          prob = prob,
          xdat = xdat,
          exdat = exdat,
          bws = bws,
          operator = operator,
          nout = nout,
          .np_ksum_unconditional_eval_exact = .np_ksum_unconditional_eval_exact
        )
      )
    }

    if (anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-ksum-unconditional-exact fan-out returned incomplete results",
        what = "inid-ksum-unconditional-exact"
      )
  }

  list(t = tmat, t0 = t0)
}

.np_inid_boot_from_ksum_unconditional <- function(xdat,
                                                  exdat,
                                                  bws,
                                                  B,
                                                  operator,
                                                  counts = NULL,
                                                  counts.drawer = NULL) {
  if (!identical(bws$type, "fixed")) {
    return(.np_inid_boot_from_ksum_unconditional_exact(
      xdat = xdat,
      exdat = exdat,
      bws = bws,
      B = B,
      operator = operator,
      counts = counts,
      counts.drawer = counts.drawer
    ))
  }

  xdat <- toFrame(xdat)
  exdat <- toFrame(exdat)
  B <- as.integer(B)
  n <- nrow(xdat)
  neval <- nrow(exdat)

  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid unconditional inid bootstrap dimensions")

  # Use direct local kernel weights to avoid nested MPI-aware npksum() calls in helper fanout.
  kw <- .np_plot_kernel_weights_direct(
    bws = bws,
    txdat = xdat,
    exdat = exdat,
    operator = operator,
    where = "direct unconditional kernel weights"
  )
  t0 <- colSums(kw) / n

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    return(list(t = crossprod(counts.mat, kw) / n, t0 = t0))
  }

  chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
    B = B,
    chunk.size = .np_inid_chunk_size(n = n, B = B),
    comm = 1L,
    include.master = TRUE
  )
  prob <- rep.int(1 / n, n)
  tmat <- matrix(NA_real_, nrow = B, ncol = neval)
  # Prebind closure inputs for worker serialization across MPI slaves.
  kw.local <- kw
  counts.drawer.local <- counts.drawer
  n.local <- n
  prob.local <- prob

  if (!is.null(counts.drawer) &&
      .npRmpi_bootstrap_fanout_enabled(
        comm = 1L,
        n = n,
        B = B,
        chunk.size = chunk.size,
        what = "inid-ksum-unconditional-block"
      )) {
    tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
    worker <- function(task) {
      start <- as.integer(task$start)
      stopi <- start + as.integer(task$bsz) - 1L
      counts.chunk <- .np_inid_counts_matrix(
        n = n.local,
        B = as.integer(task$bsz),
        counts = counts.drawer.local(start, stopi)
      )
      crossprod(counts.chunk, kw.local) / n.local
    }
    tmat <- .npRmpi_bootstrap_run_fanout(
      tasks = tasks,
      worker = worker,
      ncol.out = neval,
      what = "inid-ksum-unconditional-block",
      profile.where = "mpi.applyLB:inid-ksum-unconditional-block",
      comm = 1L,
      required.bindings = list(
        kw.local = kw.local,
        counts.drawer.local = counts.drawer.local,
        n.local = n.local
      )
    )
  }

  if (!is.null(counts.drawer) && anyNA(tmat)) {
    .npRmpi_bootstrap_fail_or_fallback(
      msg = "inid-ksum-unconditional-block fan-out returned incomplete results",
      what = "inid-ksum-unconditional-block"
    )
  }

  if (anyNA(tmat) &&
      .npRmpi_bootstrap_fanout_enabled(
        comm = 1L,
        n = n,
        B = B,
        chunk.size = chunk.size,
        what = "inid-ksum-unconditional"
      )) {
    tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
    worker <- function(task) {
      set.seed(as.integer(task$seed))
      bsz <- as.integer(task$bsz)
      counts.chunk <- stats::rmultinom(n = bsz, size = n.local, prob = prob.local)
      crossprod(counts.chunk, kw.local) / n.local
    }
    tmat <- .npRmpi_bootstrap_run_fanout(
      tasks = tasks,
      worker = worker,
      ncol.out = neval,
      what = "inid-ksum-unconditional",
      profile.where = "mpi.applyLB:inid-ksum-unconditional",
      comm = 1L,
      required.bindings = list(
        kw.local = kw.local,
        prob.local = prob.local,
        n.local = n.local
      )
    )
  }

  if (anyNA(tmat))
    .npRmpi_bootstrap_fail_or_fallback(
      msg = "inid-ksum-unconditional fan-out returned incomplete results",
      what = "inid-ksum-unconditional"
    )

  list(t = tmat, t0 = t0)
}

.np_con_inid_ksum_eligible <- function(bws) {
  isTRUE(identical(bws$cxkertype, bws$cykertype)) &&
    isTRUE(identical(bws$cxkerorder, bws$cykerorder)) &&
    isTRUE(identical(bws$cxkerbound, bws$cykerbound)) &&
    isTRUE(identical(bws$uxkertype, bws$uykertype)) &&
    isTRUE(identical(bws$oxkertype, bws$oykertype))
}

.np_con_make_kbandwidth_x <- function(bws, xdat) {
  xdat <- toFrame(xdat)
  kbandwidth.numeric(
    bw = bws$xbw,
    bwscaling = FALSE,
    # npksum helper constructors require raw bandwidths; bwscaling flags are
    # non-fit-defining here and are intentionally normalized to FALSE.
    bwtype = bws$type,
    ckertype = bws$cxkertype,
    ckerorder = bws$cxkerorder,
    ckerbound = bws$cxkerbound,
    ckerlb = if (!is.null(bws$cxkerlb)) bws$cxkerlb else NULL,
    ckerub = if (!is.null(bws$cxkerub)) bws$cxkerub else NULL,
    ukertype = bws$uxkertype,
    okertype = bws$oxkertype,
    nobs = nrow(xdat),
    xdati = untangle(xdat),
    xnames = names(xdat)
  )
}

.np_con_make_kbandwidth_xy <- function(bws, xdat, ydat) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  xydat <- data.frame(xdat, ydat)
  ckerlb <- c(if (is.null(bws$cxkerlb)) numeric(0) else bws$cxkerlb,
              if (is.null(bws$cykerlb)) numeric(0) else bws$cykerlb)
  ckerub <- c(if (is.null(bws$cxkerub)) numeric(0) else bws$cxkerub,
              if (is.null(bws$cykerub)) numeric(0) else bws$cykerub)

  kbandwidth.numeric(
    bw = c(bws$xbw, bws$ybw),
    bwscaling = FALSE,
    # npksum helper constructors require raw bandwidths; bwscaling flags are
    # non-fit-defining here and are intentionally normalized to FALSE.
    bwtype = bws$type,
    ckertype = bws$cxkertype,
    ckerorder = bws$cxkerorder,
    ckerbound = bws$cxkerbound,
    ckerlb = if (length(ckerlb)) ckerlb else NULL,
    ckerub = if (length(ckerub)) ckerub else NULL,
    ukertype = bws$uxkertype,
    okertype = bws$oxkertype,
    nobs = nrow(xydat),
    xdati = untangle(xydat),
    xnames = names(xydat)
  )
}

.np_ksum_conditional_eval_exact <- function(xdat,
                                            ydat,
                                            exdat,
                                            eydat,
                                            kbx,
                                            kbxy,
                                            cdf) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)
  xop <- rep.int("normal", ncol(xdat))
  yop <- rep.int(if (cdf) "integral" else "normal", ncol(ydat))
  xyop <- c(xop, yop)
  xydat <- data.frame(xdat, ydat)
  exydat <- data.frame(exdat, eydat)

  den <- colSums(.np_plot_kernel_weights_direct(
    bws = kbx,
    txdat = xdat,
    exdat = exdat,
    operator = xop,
    where = "direct conditional exact denominator kernel weights"
  )) / nrow(xdat)
  num <- colSums(.np_plot_kernel_weights_direct(
    bws = kbxy,
    txdat = xydat,
    exdat = exydat,
    operator = xyop,
    where = "direct conditional exact numerator kernel weights"
  )) / nrow(xydat)

  num / pmax(den, .Machine$double.eps)
}

.np_inid_boot_from_ksum_conditional_exact <- function(xdat,
                                                      ydat,
                                                      exdat,
                                                      eydat,
                                                      bws,
                                                      B,
                                                      cdf,
                                                      counts = NULL,
                                                      counts.drawer = NULL) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)
  B <- as.integer(B)
  n <- nrow(xdat)
  neval <- nrow(exdat)

  if (nrow(ydat) != n || nrow(eydat) != neval)
    stop("conditional exact bootstrap helper requires aligned x/y training and evaluation rows")
  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid conditional exact bootstrap dimensions")
  if (!.np_con_inid_ksum_eligible(bws))
    return(NULL)

  kbx <- tryCatch(.np_con_make_kbandwidth_x(bws = bws, xdat = xdat),
                  error = function(e) NULL)
  kbxy <- tryCatch(.np_con_make_kbandwidth_xy(bws = bws, xdat = xdat, ydat = ydat),
                   error = function(e) NULL)
  if (is.null(kbx) || is.null(kbxy))
    return(NULL)

  fit_one <- function(x.train, y.train) {
    .np_ksum_conditional_eval_exact(
      xdat = x.train,
      ydat = y.train,
      exdat = exdat,
      eydat = eydat,
      kbx = kbx,
      kbxy = kbxy,
      cdf = cdf
    )
  }

  t0 <- fit_one(x.train = xdat, y.train = ydat)
  nout <- length(t0)
  chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
    B = B,
    chunk.size = .np_inid_chunk_size(n = n, B = B),
    comm = 1L,
    include.master = TRUE
  )
  tmat <- matrix(NA_real_, nrow = B, ncol = nout)

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    if (.npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-ksum-conditional-exact-counts"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        start <- as.integer(task$start)
        stopi <- start + as.integer(task$bsz) - 1L
        counts.chunk <- as.matrix(counts.mat[, start:stopi, drop = FALSE])
        bsz <- ncol(counts.chunk)
        out <- matrix(NA_real_, nrow = bsz, ncol = nout)
        for (jj in seq_len(bsz)) {
          idx <- rep.int(seq_len(n), counts.chunk[, jj])
          out[jj, ] <- .np_ksum_conditional_eval_exact(
            xdat = xdat[idx, , drop = FALSE],
            ydat = ydat[idx, , drop = FALSE],
            exdat = exdat,
            eydat = eydat,
            kbx = kbx,
            kbxy = kbxy,
            cdf = cdf
          )
        }
        out
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = nout,
        what = "inid-ksum-conditional-exact-counts",
        profile.where = "mpi.applyLB:inid-ksum-conditional-exact-counts",
        comm = 1L,
        required.bindings = list(
          counts.mat = counts.mat,
          xdat = xdat,
          ydat = ydat,
          exdat = exdat,
          eydat = eydat,
          kbx = kbx,
          kbxy = kbxy,
          cdf = cdf,
          n = n,
          nout = nout,
          .np_ksum_conditional_eval_exact = .np_ksum_conditional_eval_exact
        )
      )
    }
    if (anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-ksum-conditional-exact-counts fan-out returned incomplete results",
        what = "inid-ksum-conditional-exact-counts"
      )
  } else {
    if (!is.null(counts.drawer) &&
        .npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-ksum-conditional-exact-block"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        start <- as.integer(task$start)
        stopi <- start + as.integer(task$bsz) - 1L
        counts.chunk <- .np_inid_counts_matrix(
          n = n,
          B = as.integer(task$bsz),
          counts = counts.drawer(start, stopi)
        )
        bsz <- ncol(counts.chunk)
        out <- matrix(NA_real_, nrow = bsz, ncol = nout)
        for (jj in seq_len(bsz)) {
          idx <- rep.int(seq_len(n), counts.chunk[, jj])
          out[jj, ] <- .np_ksum_conditional_eval_exact(
            xdat = xdat[idx, , drop = FALSE],
            ydat = ydat[idx, , drop = FALSE],
            exdat = exdat,
            eydat = eydat,
            kbx = kbx,
            kbxy = kbxy,
            cdf = cdf
          )
        }
        out
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = nout,
        what = "inid-ksum-conditional-exact-block",
        profile.where = "mpi.applyLB:inid-ksum-conditional-exact-block",
        comm = 1L,
        required.bindings = list(
          n = n,
          counts.drawer = counts.drawer,
          xdat = xdat,
          ydat = ydat,
          exdat = exdat,
          eydat = eydat,
          kbx = kbx,
          kbxy = kbxy,
          cdf = cdf,
          nout = nout,
          .np_ksum_conditional_eval_exact = .np_ksum_conditional_eval_exact
        )
      )
    }

    if (!is.null(counts.drawer) && anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-ksum-conditional-exact-block fan-out returned incomplete results",
        what = "inid-ksum-conditional-exact-block"
      )

    prob <- rep.int(1 / n, n)
    if (anyNA(tmat) &&
        .npRmpi_bootstrap_fanout_enabled(
          comm = 1L,
          n = n,
          B = B,
          chunk.size = chunk.size,
          what = "inid-ksum-conditional-exact"
        )) {
      tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
      worker <- function(task) {
        set.seed(as.integer(task$seed))
        bsz <- as.integer(task$bsz)
        counts.chunk <- stats::rmultinom(n = bsz, size = n, prob = prob)
        out <- matrix(NA_real_, nrow = bsz, ncol = nout)
        for (jj in seq_len(bsz)) {
          idx <- rep.int(seq_len(n), counts.chunk[, jj])
          out[jj, ] <- .np_ksum_conditional_eval_exact(
            xdat = xdat[idx, , drop = FALSE],
            ydat = ydat[idx, , drop = FALSE],
            exdat = exdat,
            eydat = eydat,
            kbx = kbx,
            kbxy = kbxy,
            cdf = cdf
          )
        }
        out
      }
      tmat <- .npRmpi_bootstrap_run_fanout(
        tasks = tasks,
        worker = worker,
        ncol.out = nout,
        what = "inid-ksum-conditional-exact",
        profile.where = "mpi.applyLB:inid-ksum-conditional-exact",
        comm = 1L,
        required.bindings = list(
          n = n,
          prob = prob,
          xdat = xdat,
          ydat = ydat,
          exdat = exdat,
          eydat = eydat,
          kbx = kbx,
          kbxy = kbxy,
          cdf = cdf,
          nout = nout,
          .np_ksum_conditional_eval_exact = .np_ksum_conditional_eval_exact
        )
      )
    }

    if (anyNA(tmat))
      .npRmpi_bootstrap_fail_or_fallback(
        msg = "inid-ksum-conditional-exact fan-out returned incomplete results",
        what = "inid-ksum-conditional-exact"
      )
  }

  list(t = tmat, t0 = t0)
}

.np_inid_boot_from_ksum_conditional <- function(xdat,
                                                ydat,
                                                exdat,
                                                eydat,
                                                bws,
                                                B,
                                                cdf,
                                                counts = NULL,
                                                counts.drawer = NULL) {
  if (!identical(bws$type, "fixed")) {
    return(.np_inid_boot_from_ksum_conditional_exact(
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      eydat = eydat,
      bws = bws,
      B = B,
      cdf = cdf,
      counts = counts,
      counts.drawer = counts.drawer
    ))
  }

  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)
  B <- as.integer(B)
  n <- nrow(xdat)
  neval <- nrow(exdat)

  if (nrow(ydat) != n || nrow(eydat) != neval)
    stop("conditional inid helper path requires aligned x/y training and evaluation rows")
  if (n < 1L || neval < 1L || B < 1L)
    stop("invalid conditional inid bootstrap dimensions")
  if (!.np_con_inid_ksum_eligible(bws))
    return(NULL)
  # Use direct local kernel weights to avoid nested MPI-aware npksum() calls in helper fanout.
  kbx <- tryCatch(.np_con_make_kbandwidth_x(bws = bws, xdat = xdat),
                  error = function(e) NULL)
  kbxy <- tryCatch(.np_con_make_kbandwidth_xy(bws = bws, xdat = xdat, ydat = ydat),
                   error = function(e) NULL)
  if (is.null(kbx) || is.null(kbxy))
    return(NULL)

  xop <- rep.int("normal", ncol(xdat))
  yop <- rep.int(if (cdf) "integral" else "normal", ncol(ydat))
  xyop <- c(xop, yop)

  xydat <- data.frame(xdat, ydat)
  exydat <- data.frame(exdat, eydat)
  den.kw <- .np_plot_kernel_weights_direct(
    bws = kbx,
    txdat = xdat,
    exdat = exdat,
    operator = xop,
    where = "direct conditional denominator kernel weights"
  )
  num.kw <- .np_plot_kernel_weights_direct(
    bws = kbxy,
    txdat = xydat,
    exdat = exydat,
    operator = xyop,
    where = "direct conditional numerator kernel weights"
  )
  den0 <- colSums(den.kw) / n
  num0 <- colSums(num.kw) / n
  t0 <- num0 / pmax(den0, .Machine$double.eps)

  if (!is.null(counts)) {
    counts.mat <- .np_inid_counts_matrix(n = n, B = B, counts = counts)
    den <- crossprod(counts.mat, den.kw) / n
    num <- crossprod(counts.mat, num.kw) / n
    return(list(t = num / pmax(den, .Machine$double.eps), t0 = t0))
  }

  chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
    B = B,
    chunk.size = .np_inid_chunk_size(n = n, B = B),
    comm = 1L,
    include.master = TRUE
  )
  prob <- rep.int(1 / n, n)
  tmat <- matrix(NA_real_, nrow = B, ncol = neval)
  # Prebind closure inputs for worker serialization across MPI slaves.
  den.kw.local <- den.kw
  num.kw.local <- num.kw
  counts.drawer.local <- counts.drawer
  n.local <- n
  prob.local <- prob

  if (!is.null(counts.drawer) &&
      .npRmpi_bootstrap_fanout_enabled(
        comm = 1L,
        n = n,
        B = B,
        chunk.size = chunk.size,
        what = "inid-ksum-conditional-block"
      )) {
    tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
    worker <- function(task) {
      start <- as.integer(task$start)
      stopi <- start + as.integer(task$bsz) - 1L
      counts.chunk <- .np_inid_counts_matrix(
        n = n.local,
        B = as.integer(task$bsz),
        counts = counts.drawer.local(start, stopi)
      )
      den <- crossprod(counts.chunk, den.kw.local) / n.local
      num <- crossprod(counts.chunk, num.kw.local) / n.local
      num / pmax(den, .Machine$double.eps)
    }
    tmat <- .npRmpi_bootstrap_run_fanout(
      tasks = tasks,
      worker = worker,
      ncol.out = neval,
      what = "inid-ksum-conditional-block",
      profile.where = "mpi.applyLB:inid-ksum-conditional-block",
      comm = 1L,
      required.bindings = list(
        den.kw.local = den.kw.local,
        num.kw.local = num.kw.local,
        counts.drawer.local = counts.drawer.local,
        n.local = n.local
      )
    )
  }

  if (!is.null(counts.drawer) && anyNA(tmat)) {
    .npRmpi_bootstrap_fail_or_fallback(
      msg = "inid-ksum-conditional-block fan-out returned incomplete results",
      what = "inid-ksum-conditional-block"
    )
  }

  if (anyNA(tmat) &&
      .npRmpi_bootstrap_fanout_enabled(
        comm = 1L,
        n = n,
        B = B,
        chunk.size = chunk.size,
        what = "inid-ksum-conditional"
      )) {
    tasks <- .npRmpi_bootstrap_chunk_tasks(B = B, chunk.size = chunk.size)
    worker <- function(task) {
      set.seed(as.integer(task$seed))
      bsz <- as.integer(task$bsz)
      counts.chunk <- stats::rmultinom(n = bsz, size = n.local, prob = prob.local)
      den <- crossprod(counts.chunk, den.kw.local) / n.local
      num <- crossprod(counts.chunk, num.kw.local) / n.local
      num / pmax(den, .Machine$double.eps)
    }
    tmat <- .npRmpi_bootstrap_run_fanout(
      tasks = tasks,
      worker = worker,
      ncol.out = neval,
      what = "inid-ksum-conditional",
      profile.where = "mpi.applyLB:inid-ksum-conditional",
      comm = 1L,
      required.bindings = list(
        den.kw.local = den.kw.local,
        num.kw.local = num.kw.local,
        prob.local = prob.local,
        n.local = n.local
      )
    )
  }

  if (anyNA(tmat))
    .npRmpi_bootstrap_fail_or_fallback(
      msg = "inid-ksum-conditional fan-out returned incomplete results",
      what = "inid-ksum-conditional"
    )

  list(t = tmat, t0 = t0)
}

gen.label = function(label, altlabel){
  paste(if (is.null(label)) altlabel else label)
}

gen.tflabel = function(condition, tlabel, flabel){
  paste(if (isTRUE(condition)) tlabel else flabel)
}

draw.error.bands = function(ex, ely, ehy, lty = 2, col = par("col")){
  lines(ex,ely,lty=lty,col=col)
  lines(ex,ehy,lty=lty,col=col)
}

draw.error.bars = function(ex, ely, ehy, hbar = TRUE, hbarscale = 0.3, lty = 2, col = par("col")){
  yy = double(3*length(ex))
  jj = seq_along(ex)*3

  yy[jj-2] = ely
  yy[jj-1] = ehy
  yy[jj] = NA
  
  xx = double(3*length(ex))
  xx[jj-2] = ex
  xx[jj-1] = ex
  xx[jj] = NA

  lines(xx,yy,lty=lty,col=col)

  if (hbar){
    ## hbars look silly if they are too wide in relation to their height
    ## this only matters in the limit of few points, since that is when
    ## hbardist may get relatively large

    golden = (1+sqrt(5))/2
    hbardist = abs(max(ex) - min(ex))/length(ex)*hbarscale

    yg = abs(yy[jj-2]-yy[jj-1])/golden
    htest = (hbardist >= yg)
    
    hdelta = pmin(yg, hbardist)/2
    xx[jj-2] = ex - hdelta
    xx[jj-1] = ex + hdelta
    
    ty = yy[jj-1]
    yy[jj-1] = yy[jj-2]

    lines(xx,yy,col=col)

    yy[jj-2] = ty
    yy[jj-1] = ty

    lines(xx,yy,col=col)
  }
}

draw.errors =
  function(ex, ely, ehy,
           plot.errors.style,
           plot.errors.bar,
           plot.errors.bar.num,
           lty,
           col = par("col")){
    if (plot.errors.style == "bar"){
      ei = seq(1,length(ex),length.out = min(length(ex),plot.errors.bar.num))
      draw.error.bars(ex = ex[ei],
                      ely = ely[ei],
                      ehy = ehy[ei],
                      hbar = (plot.errors.bar == "I"),
                      lty = lty,
                      col = col)
    } else if (plot.errors.style == "band") {
      draw.error.bands(ex = ex,
                       ely = ely,
                       ehy = ehy,
                       lty = lty,
                       col = col)
    }
  }

draw.all.error.types <- function(ex, center, all.err,
                                 plot.errors.style = "band",
                                 plot.errors.bar = "|",
                                 plot.errors.bar.num = min(length(ex), 25),
                                 lty = 2, add.legend = TRUE, legend.loc = "topleft",
                                 xi.factor = FALSE){
  if (is.null(all.err)) return(invisible(NULL))

  if (xi.factor) {
    plot.errors.style <- "bar"
    plot.errors.bar <- "I"
  }

  draw_one <- function(err, col) {
    if (is.null(err)) return(invisible(NULL))
    lower <- center - err[,1]
    upper <- center + err[,2]
    good <- complete.cases(ex, lower, upper)
    if (!any(good)) return(invisible(NULL))
    draw.errors(ex = ex[good], ely = lower[good], ehy = upper[good],
                plot.errors.style = plot.errors.style,
                plot.errors.bar = plot.errors.bar,
                plot.errors.bar.num = plot.errors.bar.num,
                lty = lty, col = col)
  }

  draw_one(all.err$pointwise, "red")
  draw_one(all.err$simultaneous, "green3")
  draw_one(all.err$bonferroni, "blue")

  if (add.legend) {
    legend(legend.loc,
           legend = c("Pointwise","Simultaneous","Bonferroni"),
           lty = 2, col = c("red","green3","blue"), lwd = 2, bty = "n")
  }
}

plotFactor <- function(f, y, ...){
  dot.args <- list(...)
  dot.names <- names(dot.args)
  has.user.lty <- !is.null(dot.names) && any(dot.names == "lty")

  if (has.user.lty) {
    do.call(plot, c(list(x = f, y = y), dot.args))
  } else {
    plot(x = f, y = y, lty = "blank", ...)
  }

  l.f = rep(f, each=3)
  l.f[3*seq_along(f)] = NA

  l.y = unlist(lapply(y, function (p) { c(0,p,NA) }))

  lines(x = l.f, y = l.y, lty = 2)
  points(x = f, y = y)
}

.np_plot_conditional_eval <- function(bws,
                                      xdat,
                                      ydat,
                                      exdat,
                                      eydat,
                                      cdf = FALSE,
                                      gradients = FALSE,
                                      proper = FALSE,
                                      proper.method = NULL,
                                      proper.control = list()) {
  fit.start <- proc.time()[3]
  proper <- npValidateScalarLogical(proper, "proper")

  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  eydat <- toFrame(eydat)

  if (nrow(xdat) != nrow(ydat))
    stop("conditional plot helper requires aligned training rows")
  if (nrow(exdat) != nrow(eydat))
    stop("conditional plot helper requires aligned evaluation rows")
  if (!xdat %~% exdat)
    stop("'xdat' and 'exdat' are not similar data frames!")
  if (!ydat %~% eydat)
    stop("'ydat' and 'eydat' are not similar data frames!")

  xdat <- adjustLevels(xdat, bws$xdati)
  ydat <- adjustLevels(ydat, bws$ydati)
  exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
  eydat <- adjustLevels(eydat, bws$ydati, allowNewCells = TRUE)
  npKernelBoundsCheckEval(exdat, bws$ixcon, bws$cxkerlb, bws$cxkerub, argprefix = "cxker")
  npKernelBoundsCheckEval(eydat, bws$iycon, bws$cykerlb, bws$cykerub, argprefix = "cyker")

  txeval <- exdat
  tyeval <- eydat
  tnrow <- nrow(xdat)
  enrow <- nrow(exdat)

  ydat <- toMatrix(ydat)
  tyuno <- ydat[, bws$iyuno, drop = FALSE]
  tycon <- ydat[, bws$iycon, drop = FALSE]
  tyord <- ydat[, bws$iyord, drop = FALSE]

  xdat <- toMatrix(xdat)
  txuno <- xdat[, bws$ixuno, drop = FALSE]
  txcon <- xdat[, bws$ixcon, drop = FALSE]
  txord <- xdat[, bws$ixord, drop = FALSE]

  eydat <- toMatrix(eydat)
  eyuno <- eydat[, bws$iyuno, drop = FALSE]
  eycon <- eydat[, bws$iycon, drop = FALSE]
  eyord <- eydat[, bws$iyord, drop = FALSE]

  exdat <- toMatrix(exdat)
  exuno <- exdat[, bws$ixuno, drop = FALSE]
  excon <- exdat[, bws$ixcon, drop = FALSE]
  exord <- exdat[, bws$ixord, drop = FALSE]

  reg.engine <- if (is.null(bws$regtype.engine)) {
    if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  } else {
    as.character(bws$regtype.engine)
  }
  basis.engine <- if (is.null(bws$basis.engine)) {
    if (is.null(bws$basis)) "glp" else bws$basis
  } else {
    bws$basis.engine
  }
  degree.engine <- if (is.null(bws$degree.engine)) {
    if (bws$xncon > 0L) {
      if (identical(reg.engine, "lc")) rep.int(0L, bws$xncon) else npValidateGlpDegree(
        regtype = "lp",
        degree = bws$degree,
        ncon = bws$xncon
      )
    } else {
      integer(0)
    }
  } else {
    as.integer(bws$degree.engine)
  }
  bernstein.engine <- if (is.null(bws$bernstein.basis.engine)) {
    isTRUE(bws$bernstein.basis)
  } else {
    isTRUE(bws$bernstein.basis.engine)
  }

  reg.c <- npRegtypeToC(
    regtype = if (identical(reg.engine, "lp")) "lp" else "lc",
    degree = degree.engine,
    ncon = bws$xncon,
    context = if (isTRUE(cdf)) "npcdist" else "npcdens"
  )
  degree.c <- if (bws$xncon > 0L) {
    as.integer(if (is.null(reg.c$degree)) rep.int(0L, bws$xncon) else reg.c$degree)
  } else {
    integer(0)
  }
  basis.code <- as.integer(npLpBasisCode(basis.engine))

  myopti <- list(
    num_obs_train = tnrow,
    num_obs_eval = enrow,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
    BANDWIDTH_den_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN
    ),
    int_MINIMIZE_IO = if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
    xkerneval = switch(bws$cxkertype,
      gaussian = CKER_GAUSS + bws$cxkerorder / 2 - 1,
      epanechnikov = CKER_EPAN + bws$cxkerorder / 2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS
    ),
    ykerneval = switch(bws$cykertype,
      gaussian = CKER_GAUSS + bws$cykerorder / 2 - 1,
      epanechnikov = CKER_EPAN + bws$cykerorder / 2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS
    ),
    uxkerneval = switch(bws$uxkertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR
    ),
    uykerneval = switch(bws$uykertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR
    ),
    oxkerneval = switch(bws$oxkertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_NLR,
      racineliyan = OKER_RLY
    ),
    oykerneval = switch(bws$oykertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_NLR,
      racineliyan = OKER_RLY
    ),
    num_yuno = bws$ynuno,
    num_yord = bws$ynord,
    num_ycon = bws$yncon,
    num_xuno = bws$xnuno,
    num_xord = bws$xnord,
    num_xcon = bws$xncon,
    no.exy = FALSE,
    gradients = gradients,
    ymcv.numRow = attr(bws$ymcv, "num.row"),
    xmcv.numRow = attr(bws$xmcv, "num.row"),
    densOrDist = if (isTRUE(cdf)) NP_DO_DIST else NP_DO_DENS,
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO
  )

  cxker.bounds.c <- npKernelBoundsMarshal(bws$cxkerlb[bws$ixcon], bws$cxkerub[bws$ixcon])
  cyker.bounds.c <- npKernelBoundsMarshal(bws$cykerlb[bws$iycon], bws$cykerub[bws$iycon])

  myout <- .np_plot_with_local_compiled_eval(.Call(
    "C_np_density_conditional",
    as.double(tyuno), as.double(tyord), as.double(tycon),
    as.double(txuno), as.double(txord), as.double(txcon),
    as.double(eyuno), as.double(eyord), as.double(eycon),
    as.double(exuno), as.double(exord), as.double(excon),
    as.double(c(
      bws$xbw[bws$ixcon], bws$ybw[bws$iycon],
      bws$ybw[bws$iyuno], bws$ybw[bws$iyord],
      bws$xbw[bws$ixuno], bws$xbw[bws$ixord]
    )),
    as.double(bws$ymcv), as.double(attr(bws$ymcv, "pad.num")),
    as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
    as.double(bws$nconfac), as.double(bws$ncatfac), as.double(bws$sdev),
    as.integer(myopti),
    as.integer(enrow),
    as.integer(bws$xndim),
    as.double(cxker.bounds.c$lb),
    as.double(cxker.bounds.c$ub),
    as.double(cyker.bounds.c$lb),
    as.double(cyker.bounds.c$ub),
    as.integer(reg.c$code),
    as.integer(degree.c),
    as.integer(bernstein.engine),
    basis.code,
    PACKAGE = "npRmpi"
  ))

  if (isTRUE(cdf))
    names(myout)[1L] <- "condist"

  if (isTRUE(gradients)) {
    xidx <- seq_len(bws$xndim)
    rorder <- numeric(bws$xndim)
    rorder[c(xidx[bws$ixcon], xidx[bws$ixuno], xidx[bws$ixord])] <- xidx
    myout$congrad <- matrix(data = myout$congrad, nrow = enrow, ncol = bws$xndim, byrow = FALSE)
    myout$congrad <- myout$congrad[, rorder, drop = FALSE]
    myout$congerr <- matrix(data = myout$congerr, nrow = enrow, ncol = bws$xndim, byrow = FALSE)
    myout$congerr <- myout$congerr[, rorder, drop = FALSE]
  } else {
    myout$congrad <- NA
    myout$congerr <- NA
  }

  fit.elapsed <- proc.time()[3] - fit.start
  optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
  total.time <- fit.elapsed + if (is.na(optim.time)) 0.0 else optim.time

  if (isTRUE(cdf)) {
    out <- condistribution(
      bws = bws,
      xeval = txeval,
      yeval = tyeval,
      condist = myout$condist,
      conderr = myout$conderr,
      congrad = myout$congrad,
      congerr = myout$congerr,
      ntrain = tnrow,
      trainiseval = FALSE,
      gradients = gradients,
      rows.omit = integer(0),
      timing = bws$timing,
      total.time = total.time,
      optim.time = optim.time,
      fit.time = fit.elapsed
    )

    if (isTRUE(proper)) {
      .np_condist_finalize_proper_object(
        object = out,
        proper = TRUE,
        proper.method = proper.method,
        proper.control = proper.control,
        where = "plot()"
      )
    } else {
      out
    }
  } else {
    out <- condensity(
      bws = bws,
      xeval = txeval,
      yeval = tyeval,
      condens = myout$condens,
      conderr = myout$conderr,
      congrad = myout$congrad,
      congerr = myout$congerr,
      ll = myout$log_likelihood,
      ntrain = tnrow,
      trainiseval = FALSE,
      gradients = gradients,
      rows.omit = integer(0),
      timing = bws$timing,
      total.time = total.time,
      optim.time = optim.time,
      fit.time = fit.elapsed
    )

    if (isTRUE(proper)) {
      .np_condens_finalize_proper_object(
        object = out,
        proper = TRUE,
        proper.method = proper.method,
        proper.control = proper.control,
        where = "plot()"
      )
    } else {
      out
    }
  }
}

.np_plot_panel_fun <- function(plot.bootstrap, plot.bxp) {
  if (plot.bootstrap && plot.bxp) bxp else plotFactor
}

.np_plot_resolve_xydat <- function(bws, xdat, ydat, miss.xy) {
  if (any(miss.xy) && !all(miss.xy))
    stop("one of, but not both, xdat and ydat was specified")

  if (all(miss.xy) && !is.null(bws$formula)) {
    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1, m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    ydat <- model.response(tmf)
    xdat <- tmf[, attr(attr(tmf, "terms"), "term.labels"), drop = FALSE]
    return(list(xdat = xdat, ydat = ydat))
  }

  if (all(miss.xy) && !is.null(bws$call)) {
    xdat <- data.frame(.np_eval_bws_call_arg(bws, "xdat"))
    ydat <- .np_eval_bws_call_arg(bws, "ydat")
  }

  xdat <- toFrame(xdat)
  goodrows <- seq_len(nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (!is.null(rows.omit))
    goodrows[rows.omit] <- 0

  if (all(goodrows == 0))
    stop("Data has no rows without NAs")

  goodrows <- goodrows[goodrows != 0]
  list(xdat = xdat[goodrows, , drop = FALSE],
       ydat = ydat[goodrows])
}

.np_plot_regression_eval <- function(bws,
                                     xdat,
                                     ydat,
                                     exdat,
                                     gradients = FALSE,
                                     gradient.order = 1L,
                                     need.asymptotic = FALSE) {
  if (isTRUE(need.asymptotic)) {
    return(npreg(
      txdat = xdat,
      tydat = ydat,
      exdat = exdat,
      bws = bws,
      gradients = gradients,
      gradient.order = gradient.order,
      warn.glp.gradient = FALSE
    ))
  }

  fit <- .np_regression_direct(
    bws = bws,
    txdat = xdat,
    tydat = ydat,
    exdat = exdat,
    gradients = gradients,
    gradient.order = gradient.order,
    local.mode = identical(bws$type, "generalized_nn")
  )

  neval <- length(fit$mean)
  fit$merr <- rep(NA_real_, neval)
  if (isTRUE(gradients))
    fit$gerr <- matrix(NA_real_, nrow = neval, ncol = NCOL(fit$grad))

  fit
}

.np_plot_unconditional_eval <- function(xdat,
                                        exdat,
                                        bws,
                                        cdf = FALSE,
                                        need.asymptotic = FALSE) {
  if (isTRUE(need.asymptotic)) {
    return(if (isTRUE(cdf)) {
      npudist(tdat = xdat, edat = exdat, bws = bws)
    } else {
      npudens(tdat = xdat, edat = exdat, bws = bws)
    })
  }

  est <- .np_ksum_unconditional_eval_exact(
    xdat = xdat,
    exdat = exdat,
    bws = bws,
    operator = if (isTRUE(cdf)) "integral" else "normal"
  )

  if (isTRUE(cdf)) {
    list(dist = est, derr = rep(NA_real_, length(est)))
  } else {
    list(dens = est, derr = rep(NA_real_, length(est)))
  }
}

.npRmpi_with_local_bootstrap <- function(expr) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  options(npRmpi.autodispatch.context = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
  force(expr)
}

.np_plot_with_local_compiled_eval <- function(expr) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.local <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old.local), add = TRUE)
  old.mode <- .Call("C_np_set_local_regression_mode", TRUE, PACKAGE = "npRmpi")
  on.exit(.Call("C_np_set_local_regression_mode", old.mode, PACKAGE = "npRmpi"), add = TRUE)
  force(expr)
}

.npRmpi_plot_behavior_for_rank <- function(plot.behavior) {
  if (!.npRmpi_autodispatch_called_from_bcast())
    return(plot.behavior)

  rank <- tryCatch(mpi.comm.rank(), error = function(e) NA_integer_)
  if (is.na(rank) || rank == 0L)
    return(plot.behavior)

  "data"
}

.npRmpi_profile_env <- local({
  env <- new.env(parent = emptyenv())
  env$last <- NULL
  env$history <- list()
  env$active_id <- NULL
  env$active <- list()
  env
})

.npRmpi_profile_enabled <- function() {
  !isFALSE(getOption("npRmpi.profile.enable", TRUE))
}

.npRmpi_profile_level <- function() {
  lvl <- as.character(getOption("npRmpi.profile.level", "basic"))[1L]
  if (is.na(lvl) || !(lvl %in% c("basic", "detailed")))
    lvl <- "basic"
  lvl
}

.npRmpi_profile_history_limit <- function() {
  lim <- suppressWarnings(as.integer(getOption("npRmpi.profile.history.max", 200L))[1L])
  if (is.na(lim) || lim < 1L)
    lim <- 200L
  lim
}

.np_nrows_safe <- function(x) {
  if (is.null(x))
    return(NA_integer_)
  if (is.atomic(x) && is.null(dim(x)))
    return(as.integer(length(x)))
  xf <- tryCatch(toFrame(x), error = function(e) NULL)
  if (is.null(xf))
    return(NA_integer_)
  as.integer(nrow(xf))
}

.npRmpi_profile_bootstrap_begin <- function(where, method = NA_character_,
                                            B = NA_integer_, ntrain = NA_integer_,
                                            neval = NA_integer_) {
  if (!isTRUE(.npRmpi_profile_enabled()))
    return(NULL)

  rank <- tryCatch(as.integer(mpi.comm.rank()), error = function(e) NA_integer_)
  size <- tryCatch(as.integer(mpi.comm.size()), error = function(e) NA_integer_)
  active <- list(
    where = where,
    level = .npRmpi_profile_level(),
    method = if (length(method)) as.character(method)[1L] else NA_character_,
    B = suppressWarnings(as.integer(B)[1L]),
    ntrain = suppressWarnings(as.integer(ntrain)[1L]),
    neval = suppressWarnings(as.integer(neval)[1L]),
    rank = rank,
    size = size,
    via_bcast = isTRUE(.npRmpi_autodispatch_called_from_bcast()),
    comm_elapsed_sec = 0.0,
    comm_calls = 0L,
    comm_notes = character(0),
    start_proc = proc.time(),
    start_wall = Sys.time()
  )
  active$id <- paste0(format(active$start_wall, "%Y%m%d%H%M%OS6"), "-", sample.int(1e8, 1L))
  .npRmpi_profile_env$active_id <- active$id
  .npRmpi_profile_env$active[[active$id]] <- active
  active
}

.npRmpi_profile_add_comm_elapsed <- function(elapsed_sec, where = NA_character_) {
  id <- .npRmpi_profile_env$active_id
  if (is.null(id))
    return(invisible(FALSE))
  cur <- .npRmpi_profile_env$active[[id]]
  if (is.null(cur))
    return(invisible(FALSE))

  elapsed <- suppressWarnings(as.double(elapsed_sec)[1L])
  if (is.na(elapsed) || elapsed < 0)
    elapsed <- 0.0

  cur$comm_elapsed_sec <- as.double(cur$comm_elapsed_sec) + elapsed
  cur$comm_calls <- as.integer(cur$comm_calls) + 1L
  if (!is.na(where) && nzchar(as.character(where)[1L]) &&
      identical(cur$level, "detailed")) {
    cur$comm_notes <- c(cur$comm_notes, as.character(where)[1L])
  }

  .npRmpi_profile_env$active[[id]] <- cur
  invisible(TRUE)
}

.npRmpi_profile_bootstrap_end <- function(value, ctx) {
  if (is.null(ctx))
    return(value)

  id <- ctx$id
  active <- if (!is.null(id)) .npRmpi_profile_env$active[[id]] else NULL
  if (is.list(active)) {
    ctx$comm_elapsed_sec <- as.double(active$comm_elapsed_sec)
    ctx$comm_calls <- as.integer(active$comm_calls)
    if (identical(active$level, "detailed"))
      ctx$comm_notes <- active$comm_notes
  }

  dt <- proc.time() - ctx$start_proc
  wall <- unname(as.double(dt[["elapsed"]]))
  comm <- suppressWarnings(as.double(ctx$comm_elapsed_sec)[1L])
  if (is.na(comm) || comm < 0)
    comm <- 0.0
  if (!is.finite(wall) || wall <= 0)
    wall <- 0.0
  compute <- max(0.0, wall - comm)
  denom <- comm + compute
  ratio <- if (denom > 0) min(1.0, max(0.0, comm / denom)) else NA_real_

  record <- list(
    where = ctx$where,
    level = ctx$level,
    method = ctx$method,
    B = ctx$B,
    ntrain = ctx$ntrain,
    neval = ctx$neval,
    rank = ctx$rank,
    size = ctx$size,
    via_bcast = ctx$via_bcast,
    wall_elapsed_sec = wall,
    comm_elapsed_sec = comm,
    compute_elapsed_sec = compute,
    comm_ratio = ratio,
    comm_calls = suppressWarnings(as.integer(ctx$comm_calls)[1L]),
    user_sec = unname(as.double(dt[["user.self"]])),
    system_sec = unname(as.double(dt[["sys.self"]])),
    timestamp_start = ctx$start_wall,
    timestamp_end = Sys.time()
  )
  if (identical(ctx$level, "detailed")) {
    record$comm_notes <- if (length(ctx$comm_notes)) ctx$comm_notes else character(0)
  }

  .npRmpi_profile_env$last <- record
  .npRmpi_profile_env$history <- c(.npRmpi_profile_env$history, list(record))
  if (!is.null(id)) {
    .npRmpi_profile_env$active[[id]] <- NULL
    if (identical(.npRmpi_profile_env$active_id, id))
      .npRmpi_profile_env$active_id <- NULL
  }
  keep <- .npRmpi_profile_history_limit()
  if (length(.npRmpi_profile_env$history) > keep) {
    .npRmpi_profile_env$history <- tail(.npRmpi_profile_env$history, keep)
  }

  if (is.list(value))
    value$timing.profile <- record
  value
}

.npRmpi_profile_finalize_bootstrap <- function(boot.err, bxp, boot.all.err, ctx) {
  out <- list(boot.err = boot.err, bxp = bxp, boot.all.err = boot.all.err)
  .npRmpi_profile_bootstrap_end(out, ctx)
}

.npRmpi_profile_last <- function() {
  .npRmpi_profile_env$last
}

.npRmpi_profile_history <- function(n = NULL) {
  h <- .npRmpi_profile_env$history
  if (is.null(n))
    return(h)
  n <- suppressWarnings(as.integer(n)[1L])
  if (is.na(n) || n < 1L)
    return(list())
  tail(h, n)
}

.npRmpi_profile_clear <- function() {
  .npRmpi_profile_env$last <- NULL
  .npRmpi_profile_env$history <- list()
  .npRmpi_profile_env$active_id <- NULL
  .npRmpi_profile_env$active <- list()
  invisible(NULL)
}

## Rank-based simultaneous confidence set helper, vendored from
## MCPAN::SCSrank (MCPAN 1.1-21, GPL-2; Schaarschmidt, Gerhard, Sill).
np.plot.SCSrank <- function(x, conf.level = 0.95, alternative = "two.sided", ...) {
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))

  DataMatrix <- x
  N <- nrow(DataMatrix)
  k <- round(conf.level * N, 0)
  RankDat <- apply(DataMatrix, 2, rank)

  switch(alternative,
    "two.sided" = {
      W1 <- apply(RankDat, 1, max)
      W2 <- N + 1 - apply(RankDat, 1, min)

      Wmat <- cbind(W1, W2)
      w <- apply(Wmat, 1, max)
      tstar <- round(sort(w)[k], 0)

      SCI <- function(x) {
        sortx <- sort(x)
        cbind(sortx[N + 1 - tstar], sortx[tstar])
      }

      SCS <- t(apply(DataMatrix, 2, SCI))
    },
    "less" = {
      W1 <- apply(RankDat, 1, max)
      tstar <- round(sort(W1)[k], 0)

      SCI <- function(x) {
        sortx <- sort(x)
        cbind(-Inf, sortx[tstar])
      }

      SCS <- t(apply(DataMatrix, 2, SCI))
    },
    "greater" = {
      W2 <- N + 1 - apply(RankDat, 1, min)
      tstar <- round(sort(W2)[k], 0)

      SCI <- function(x) {
        sortx <- sort(x)
        cbind(sortx[N + 1 - tstar], Inf)
      }

      SCS <- t(apply(DataMatrix, 2, SCI))
    }
  )

  colnames(SCS) <- c("lower", "upper")

  attr(SCS, which = "k") <- k
  attr(SCS, which = "N") <- N
  OUT <- list(conf.int = SCS, conf.level = conf.level, alternative = alternative)
  return(OUT)
}

compute.bootstrap.quantile.bounds <- function(boot.t, alpha, band.type, warn.coverage = TRUE) {
  B <- nrow(boot.t)
  neval <- ncol(boot.t)

  .np_plot_bootstrap_tail_warning <- function(B, alpha, band.type, neval) {
    m <- if (identical(band.type, "pointwise")) 1L else max(1L, as.integer(neval))
    min.B <- ceiling((2.0 * m) / alpha - 1.0)
    if (B >= min.B)
      return(invisible(NULL))

    m.desc <- if (m == 1L) {
      "m=1 (pointwise tails)"
    } else {
      sprintf("m=n.eval=%d (Bonferroni-conservative tails)", neval)
    }
    warning(sprintf(
      paste0("plot.errors.boot.num=%d is too small for plot.errors.type='%s' ",
             "(alpha=%g). Minimum recommended is %d using ",
             "B >= ceiling(2*m/alpha - 1), with %s. ",
             "For 2D perspective plots on a full neval x neval grid, m=neval^2."),
      B, band.type, alpha, min.B, m.desc
    ), call. = FALSE)
  }

  if (isTRUE(warn.coverage) && band.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
    warn.type <- if (identical(band.type, "all")) "bonferroni/simultaneous/all" else band.type
    .np_plot_bootstrap_tail_warning(B = B, alpha = alpha, band.type = warn.type, neval = neval)
  }

  if (band.type == "pointwise") {
    probs <- c(alpha / 2.0, 1.0 - alpha / 2.0)
    return(t(apply(boot.t, 2, quantile, probs = probs)))
  }

  if (band.type == "bonferroni") {
    probs <- c(alpha / (2.0 * neval), 1.0 - alpha / (2.0 * neval))
    return(t(apply(boot.t, 2, quantile, probs = probs)))
  }

  if (band.type == "simultaneous") {
    return(np.plot.SCSrank(boot.t, conf.level = 1.0 - alpha)$conf.int)
  }

  if (band.type == "all") {
    return(list(
      pointwise = compute.bootstrap.quantile.bounds(boot.t, alpha, "pointwise", warn.coverage = FALSE),
      bonferroni = compute.bootstrap.quantile.bounds(boot.t, alpha, "bonferroni", warn.coverage = FALSE),
      simultaneous = compute.bootstrap.quantile.bounds(boot.t, alpha, "simultaneous", warn.coverage = FALSE)
    ))
  }

  stop("'band.type' must be one of pointwise, bonferroni, simultaneous, all")
}

.np_plot_asymptotic_error_from_se <- function(se, alpha, band.type, m = length(se)) {
  se <- as.numeric(se)
  n <- length(se)
  m <- max(1L, as.integer(m[1L]))

  if (!length(se))
    return(list(err = matrix(numeric(0), nrow = 0L, ncol = 2L), all.err = NULL))

  make_err <- function(mult) {
    cbind(mult * se, mult * se)
  }

  if (band.type == "all") {
    err.pointwise <- make_err(qnorm(alpha / 2.0, lower.tail = FALSE))
    err.bonf <- make_err(qnorm(alpha / (2.0 * m), lower.tail = FALSE))
    err.sim <- matrix(NA_real_, nrow = n, ncol = 2L)
    return(list(
      err = err.pointwise,
      all.err = list(
        pointwise = err.pointwise,
        simultaneous = err.sim,
        bonferroni = err.bonf
      )
    ))
  }

  if (band.type == "simultaneous")
    return(list(err = matrix(NA_real_, nrow = n, ncol = 2L), all.err = NULL))

  mult <- if (band.type == "bonferroni") {
    qnorm(alpha / (2.0 * m), lower.tail = FALSE)
  } else if (band.type %in% c("pmzsd", "pointwise")) {
    qnorm(alpha / 2.0, lower.tail = FALSE)
  } else {
    stop("unsupported asymptotic interval type")
  }

  list(err = make_err(mult), all.err = NULL)
}

compute.all.error.range <- function(center, all.err) {
  if (is.null(all.err)) {
    return(c(NA_real_, NA_real_))
  }
  lower <- c(center - all.err$pointwise[,1],
             center - all.err$simultaneous[,1],
             center - all.err$bonferroni[,1])
  upper <- c(center + all.err$pointwise[,2],
             center + all.err$simultaneous[,2],
             center + all.err$bonferroni[,2])
  rng <- c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE))
  if (all(is.finite(rng)))
    return(rng)

  center <- center[is.finite(center)]
  if (!length(center))
    return(c(NA_real_, NA_real_))
  range(center, finite = TRUE)
}

compute.default.error.range <- function(center, err) {
  lower <- c(center - err[,1], err[,3] - err[,1])
  upper <- c(center + err[,2], err[,3] + err[,2])
  rng <- c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE))
  if (all(is.finite(rng)))
    return(rng)

  center <- center[is.finite(center)]
  if (!length(center))
    return(c(NA_real_, NA_real_))
  range(center, finite = TRUE)
}

.np_plot_normalize_common_options <- function(plot.behavior,
                                             plot.errors.method,
                                             plot.errors.boot.method,
                                             plot.errors.boot.wild = c("rademacher", "mammen"),
                                             plot.errors.boot.blocklen,
                                             plot.errors.center,
                                             plot.errors.type,
                                             plot.errors.alpha,
                                             plot.errors.style,
                                             plot.errors.bar,
                                             xdat,
                                             common.scale,
                                             ylim,
                                             allow_asymptotic_quantile = TRUE) {
  scalar_choice <- function(value, default) {
    if (is.null(value) || length(value) < 1L || is.na(value[1L])) default else value[1L]
  }

  plot.behavior <- match.arg(
    scalar_choice(plot.behavior, "plot"),
    c("plot", "plot-data", "data")
  )
  plot.errors.method <- match.arg(
    scalar_choice(plot.errors.method, "none"),
    c("none", "bootstrap", "asymptotic")
  )
  plot.errors.boot.method <- match.arg(
    scalar_choice(plot.errors.boot.method, "wild"),
    c("wild", "inid", "fixed", "geom")
  )
  plot.errors.boot.wild <- match.arg(
    scalar_choice(plot.errors.boot.wild, "rademacher"),
    c("rademacher", "mammen")
  )
  plot.errors.center <- match.arg(
    scalar_choice(plot.errors.center, "estimate"),
    c("estimate", "bias-corrected")
  )
  plot.errors.type <- match.arg(
    scalar_choice(plot.errors.type, "simultaneous"),
    c("simultaneous", "pointwise", "bonferroni", "pmzsd", "all")
  )

  if (identical(plot.errors.method, "bootstrap") &&
      isTRUE(.npRmpi_autodispatch_called_from_bcast()))
    stop(
      "cannot run bootstrap plot paths inside mpi.bcast.cmd context; invoke plot(...) from master context with npRmpi.autodispatch=TRUE",
      call. = FALSE
    )

  if (!is.numeric(plot.errors.alpha) || length(plot.errors.alpha) != 1 ||
      is.na(plot.errors.alpha) || plot.errors.alpha <= 0 || plot.errors.alpha >= 0.5)
    stop("the tail probability plot.errors.alpha must lie in (0,0.5)")

  plot.errors.style <- match.arg(
    scalar_choice(plot.errors.style, "band"),
    c("band", "bar")
  )
  plot.errors.bar <- match.arg(
    scalar_choice(plot.errors.bar, "|"),
    c("|", "I")
  )

  common.scale <- common.scale | (!is.null(ylim))

  if (plot.errors.method == "none" && plot.errors.type == "all") {
    warning("plot.errors.type='all' requires bootstrap errors; setting plot.errors.method='bootstrap'")
    plot.errors.method <- "bootstrap"
  }

  if (allow_asymptotic_quantile && plot.errors.method == "asymptotic") {
    if (plot.errors.type == "simultaneous")
      warning("asymptotic simultaneous confidence bands are unavailable here; returning NA interval limits")
    if (plot.errors.type == "all")
      warning("asymptotic simultaneous confidence bands are unavailable here; 'all' returns pointwise/bonferroni and NA simultaneous")

    if (plot.errors.center == "bias-corrected") {
      warning("no bias corrections can be calculated with asymptotics, centering on estimate")
      plot.errors.center <- "estimate"
    }
  }

  if (is.element(plot.errors.boot.method, c("fixed", "geom")) &&
      is.null(plot.errors.boot.blocklen))
    plot.errors.boot.blocklen <- b.star(xdat, round = TRUE)[1,1]

  list(
    plot.behavior = plot.behavior,
    plot.errors.method = plot.errors.method,
    plot.errors.boot.method = plot.errors.boot.method,
    plot.errors.boot.wild = plot.errors.boot.wild,
    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    plot.errors.style = plot.errors.style,
    plot.errors.bar = plot.errors.bar,
    common.scale = common.scale,
    plot.errors = (plot.errors.method != "none")
  )
}


compute.bootstrap.errors = function(...,bws){
  UseMethod("compute.bootstrap.errors",bws)
}

compute.bootstrap.errors.rbandwidth =
  function(xdat, ydat,
           exdat,
           gradients,
           gradient.order,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.rbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.rbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    boot.out <- NULL

    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL
    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method == "inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))

    if (is.wild.hat && gradients) {
      cont.idx <- which(bws$xdati$icon)
      if (is.na(match(slice.index, cont.idx))) {
        stop("plot.errors.boot.method='wild' supports gradients only for continuous slices in npRmpi; no serial fallback is permitted", call. = FALSE)
      }
    }

    inid.helper.ok <- isTRUE(.np_plot_inid_fastpath_enabled())
    block.helper.ok <- isTRUE(.np_plot_block_fastpath_enabled())

    if (is.inid && !isTRUE(inid.helper.ok))
      stop("inid regression helper unavailable for this configuration in npRmpi; no serial fallback is permitted", call. = FALSE)
    if (is.block && !isTRUE(block.helper.ok))
      stop("fixed/geom regression helper unavailable for this configuration in npRmpi; no serial fallback is permitted", call. = FALSE)

    if ((is.inid && isTRUE(inid.helper.ok)) || (is.block && isTRUE(block.helper.ok))) {
      counts.drawer <- if (is.block) {
        .np_block_counts_drawer(
          n = nrow(xdat),
          B = plot.errors.boot.num,
          blocklen = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method
        )
      } else {
        NULL
      }
      boot.out <- tryCatch(
        .np_inid_boot_from_regression(
          xdat = xdat,
          exdat = exdat,
          bws = bws,
          ydat = ydat,
          B = plot.errors.boot.num,
          counts.drawer = counts.drawer,
          gradients = gradients,
          gradient.order = gradient.order,
          slice.index = slice.index
        ),
        error = function(e) {
          stop(sprintf("%s regression helper failed in compute.bootstrap.errors.rbandwidth (%s)",
                       if (is.block) plot.errors.boot.method else "inid",
                       conditionMessage(e)),
               call. = FALSE)
        }
      )
    }

    if (is.null(boot.out) && is.wild.hat) {
      if (length(plot.errors.boot.wild) > 1L)
        plot.errors.boot.wild <- plot.errors.boot.wild[1L]
      plot.errors.boot.wild <- match.arg(plot.errors.boot.wild, c("mammen", "rademacher"))
      .npRmpi_bootstrap_transport_trace(
        what = "rbandwidth.wild",
        event = "wild.fit.start",
        fields = list(
          slice = slice.index,
          factor_slice = isTRUE(slice.index > 0L && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])),
          n_eval = nrow(exdat)
        )
      )

      fit.train <- suppressWarnings(npreg.rbandwidth(
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        gradients = FALSE,
        warn.glp.gradient = FALSE
      ))
      .npRmpi_bootstrap_transport_trace(
        what = "rbandwidth.wild",
        event = "wild.fit.done",
        fields = list(
          slice = slice.index,
          n_train = length(fit.train$mean)
        )
      )

      s.vec <- NULL
      if (gradients) {
        cont.idx <- which(bws$xdati$icon)
        cpos <- match(slice.index, cont.idx)
        gorder <- if (length(gradient.order) == 1L) {
          rep.int(as.integer(gradient.order), length(cont.idx))
        } else {
          as.integer(gradient.order)
        }
        if (length(gorder) != length(cont.idx))
          gorder <- rep.int(1L, length(cont.idx))
        s.vec <- integer(length(cont.idx))
        s.vec[cpos] <- gorder[cpos]
      }
      .npRmpi_bootstrap_transport_trace(
        what = "rbandwidth.wild",
        event = "wild.hat.start",
        fields = list(
          slice = slice.index,
          gradients = gradients
        )
      )

      H <- suppressWarnings(npreghat.rbandwidth(
        bws = bws,
        txdat = xdat,
        exdat = exdat,
        s = s.vec,
        output = "matrix"
      ))
      .npRmpi_bootstrap_transport_trace(
        what = "rbandwidth.wild",
        event = "wild.hat.done",
        fields = list(
          slice = slice.index,
          h_rows = nrow(H),
          h_cols = ncol(H)
        )
      )

      t0 <- as.vector(H %*% as.double(ydat))
      eps <- as.double(ydat - fit.train$mean)
      n <- length(eps)
      B <- plot.errors.boot.num
      .npRmpi_bootstrap_transport_trace(
        what = "rbandwidth.wild",
        event = "wild.boot.start",
        fields = list(
          slice = slice.index,
          n = n,
          B = B
        )
      )

      boot.out <- list(
        t = .np_wild_boot_t(
          H = H,
          fit.mean = fit.train$mean,
          residuals = eps,
          B = B,
          wild = plot.errors.boot.wild
        ),
        t0 = t0
      )
      .npRmpi_bootstrap_transport_trace(
        what = "rbandwidth.wild",
        event = "wild.boot.done",
        fields = list(
          slice = slice.index,
          t_rows = nrow(boot.out$t),
          t_cols = ncol(boot.out$t)
        )
      )
    }

    if (is.null(boot.out))
      stop("no MPI helper path available for this regression bootstrap configuration in npRmpi; no serial fallback is permitted", call. = FALSE)

    all.bp <- list()

    if (slice.index > 0 && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- bws$xdati$all.lev[[slice.index]]
      rm(boot.frame)
    }
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise",
          warn.coverage = FALSE)
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.scbandwidth =
  function(xdat, ydat, zdat,
           exdat, ezdat,
           gradients,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.scbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.scbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    miss.z <- missing(zdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method == "inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))
    boot.out <- NULL

    if (is.inid) {
      if (!isTRUE(.np_plot_inid_fastpath_enabled()))
        stop("inid bootstrap requires fastpath-enabled helper for smooth coefficient plots", call. = FALSE)
      if (isTRUE(gradients))
        stop("inid bootstrap for smooth coefficient gradients is not supported in helper mode", call. = FALSE)
      boot.out <- .npRmpi_with_local_bootstrap({
        tryCatch(
          .np_inid_boot_from_scoef(
            txdat = xdat,
            ydat = ydat,
            tzdat = if (miss.z) NULL else zdat,
            exdat = exdat,
            ezdat = if (miss.z) NULL else ezdat,
            bws = bws,
            B = plot.errors.boot.num,
            leave.one.out = FALSE
          ),
          error = function(e) {
            stop(sprintf("inid smooth coefficient helper failed in compute.bootstrap.errors.scbandwidth (%s)",
                         conditionMessage(e)),
                 call. = FALSE)
          }
        )
      })
    }
    if (is.null(boot.out) && is.block) {
      if (!isTRUE(.np_plot_block_fastpath_enabled()))
        stop("fixed/geom bootstrap requires fastpath-enabled helper for smooth coefficient plots", call. = FALSE)
      if (isTRUE(gradients))
        stop("fixed/geom bootstrap for smooth coefficient gradients is not supported in helper mode", call. = FALSE)
      counts.drawer <- .np_block_counts_drawer(
        n = nrow(xdat),
        B = plot.errors.boot.num,
        blocklen = plot.errors.boot.blocklen,
        sim = plot.errors.boot.method
      )
      boot.out <- .npRmpi_with_local_bootstrap({
        tryCatch(
          .np_inid_boot_from_scoef(
            txdat = xdat,
            ydat = ydat,
            tzdat = if (miss.z) NULL else zdat,
            exdat = exdat,
            ezdat = if (miss.z) NULL else ezdat,
            bws = bws,
            B = plot.errors.boot.num,
            counts.drawer = counts.drawer,
            leave.one.out = FALSE
          ),
          error = function(e) {
            stop(sprintf("%s smooth coefficient helper failed in compute.bootstrap.errors.scbandwidth (%s)",
                         plot.errors.boot.method,
                         conditionMessage(e)),
                 call. = FALSE)
          }
        )
      })
    }

    if (is.null(boot.out) && is.wild.hat) {
      stop("plot.errors.boot.method='wild' is unsupported for smooth coefficient bootstrap in npRmpi canonical SPMD mode; use 'inid', 'fixed', or 'geom'", call. = FALSE)
    }

    if (is.null(boot.out))
      stop("no MPI helper path available for this smooth coefficient bootstrap configuration in npRmpi; no serial fallback is permitted", call. = FALSE)

    all.bp <- list()

    if ((slice.index > 0) && (((slice.index <= ncol(xdat)) && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])) ||
                              ((slice.index > ncol(xdat)) && (bws$zdati$iord[slice.index-ncol(xdat)] || bws$zdati$iuno[slice.index-ncol(xdat)])))) {
      boot.frame <- as.data.frame(boot.out$t)

      if(slice.index <= ncol(xdat))
          u.lev <- bws$xdati$all.ulev[[slice.index]]
      else
          u.lev <- bws$zdati$all.ulev[[slice.index-ncol(xdat)]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))

      if(slice.index <= ncol(xdat))
          all.bp$names <- bws$xdati$all.lev[[slice.index]]
      else
          all.bp$names <- bws$zdati$all.lev[[slice.index-ncol(xdat)]]
      rm(boot.frame)
    }
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise",
          warn.coverage = FALSE)
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.plbandwidth =
  function(xdat, ydat, zdat,
           exdat, ezdat,
           gradients,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.plbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.plbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method == "inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))

    if (is.wild.hat) {
      if (length(plot.errors.boot.wild) > 1L)
        plot.errors.boot.wild <- plot.errors.boot.wild[1L]
      plot.errors.boot.wild <- match.arg(plot.errors.boot.wild, c("mammen", "rademacher"))

      fit.mean <- as.vector(npplreghat(
        bws = bws,
        txdat = xdat,
        tzdat = zdat,
        exdat = xdat,
        ezdat = zdat,
        y = ydat,
        output = "apply"
      ))
      H <- npplreghat(
        bws = bws,
        txdat = xdat,
        tzdat = zdat,
        exdat = exdat,
        ezdat = ezdat,
        output = "matrix"
      )

      t0 <- as.vector(H %*% as.double(ydat))
      eps <- as.double(ydat - fit.mean)
      n <- length(eps)
      B <- plot.errors.boot.num

      boot.out <- .npRmpi_with_local_bootstrap({
        list(
          t = .np_wild_boot_t(
            H = H,
            fit.mean = fit.mean,
            residuals = eps,
            B = B,
            wild = plot.errors.boot.wild
          ),
          t0 = t0
        )
      })
    } else {
      boot.out <- NULL
      if (is.inid) {
        if (!isTRUE(.np_plot_inid_fastpath_enabled()))
          stop("inid bootstrap requires fastpath-enabled helper for partially linear plots", call. = FALSE)
        boot.out <- tryCatch(
          .np_inid_boot_from_plreg(
            txdat = xdat,
            ydat = ydat,
            tzdat = zdat,
            exdat = exdat,
            ezdat = ezdat,
            bws = bws,
            B = plot.errors.boot.num
          ),
          error = function(e) {
            stop(sprintf("inid plreg helper failed in compute.bootstrap.errors.plbandwidth (%s)",
                         conditionMessage(e)),
                 call. = FALSE)
          }
        )
      } else if (is.block) {
        if (!isTRUE(.np_plot_block_fastpath_enabled()))
          stop("fixed/geom bootstrap requires fastpath-enabled helper for partially linear plots", call. = FALSE)
        counts.drawer <- .np_block_counts_drawer(
          n = nrow(xdat),
          B = plot.errors.boot.num,
          blocklen = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method
        )
        boot.out <- tryCatch(
          .np_inid_boot_from_plreg(
            txdat = xdat,
            ydat = ydat,
            tzdat = zdat,
            exdat = exdat,
            ezdat = ezdat,
            bws = bws,
            B = plot.errors.boot.num,
            counts.drawer = counts.drawer
          ),
          error = function(e) {
            stop(sprintf("%s plreg helper failed in compute.bootstrap.errors.plbandwidth (%s)",
                         plot.errors.boot.method,
                         conditionMessage(e)),
                 call. = FALSE)
          }
        )
      }

      if (is.null(boot.out))
        stop("no MPI helper path available for this partially linear bootstrap configuration in npRmpi; no serial fallback is permitted", call. = FALSE)
    }

    all.bp <- list()

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$zdati
      ti <- slice.index - bws$xndim
    }
    
    if (slice.index > 0 && (tdati$iord[ti] || tdati$iuno[ti])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- tdati$all.lev[[ti]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise",
          warn.coverage = FALSE)
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.bandwidth =
  function(xdat, 
           exdat,
           cdf,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.bandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.bandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    .np_plot_reject_wild_unsupervised(plot.errors.boot.method, "unconditional density/distribution estimators")
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))
    fast.inid <- isTRUE(.np_plot_inid_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.inid)
    fast.block <- isTRUE(.np_plot_block_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.block)

    if (is.inid && !isTRUE(fast.inid))
      stop("inid unconditional helper unavailable for this configuration in npRmpi; no serial fallback is permitted", call. = FALSE)
    if (is.block && !isTRUE(fast.block))
      stop("fixed/geom unconditional helper unavailable for this configuration in npRmpi; no serial fallback is permitted", call. = FALSE)

    boot.out <- NULL
    if (fast.inid || fast.block) {
      op <- if (cdf) "integral" else "normal"
      counts.drawer <- if (fast.block) {
        .np_block_counts_drawer(
          n = nrow(xdat),
          B = plot.errors.boot.num,
          blocklen = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method
        )
      } else {
        NULL
      }
      boot.out <- .npRmpi_with_local_bootstrap({
        tryCatch(
          .np_inid_boot_from_ksum_unconditional(
            xdat = xdat,
            exdat = exdat,
            bws = bws,
            B = plot.errors.boot.num,
            operator = op,
            counts.drawer = counts.drawer
          ),
          error = function(e) {
            stop(sprintf("%s ksum helper path failed in compute.bootstrap.errors.bandwidth (%s)",
                         if (fast.block) plot.errors.boot.method else "inid",
                         conditionMessage(e)),
                 call. = FALSE)
          }
        )
      })
    }

    if (is.null(boot.out))
      stop("no MPI helper path available for this unconditional bootstrap configuration in npRmpi; no serial fallback is permitted", call. = FALSE)

    all.bp <- list()

    if (slice.index > 0 && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- bws$xdati$all.lev[[slice.index]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise",
          warn.coverage = FALSE)
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.dbandwidth =
  function(xdat, 
           exdat,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.dbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.dbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    .np_plot_reject_wild_unsupervised(plot.errors.boot.method, "unconditional density/distribution estimators")
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    is.inid = plot.errors.boot.method=="inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))
    fast.inid <- isTRUE(.np_plot_inid_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.inid) &&
      isTRUE(identical(bws$type, "fixed"))
    fast.block <- isTRUE(.np_plot_block_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.block) &&
      isTRUE(identical(bws$type, "fixed"))

    boot.out <- NULL
    if (fast.inid || fast.block) {
      counts.drawer <- if (fast.block) {
        .np_block_counts_drawer(
          n = nrow(xdat),
          B = plot.errors.boot.num,
          blocklen = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method
        )
      } else {
        NULL
      }
      boot.out <- .npRmpi_with_local_bootstrap({
        tryCatch(
          .np_inid_boot_from_ksum_unconditional(
            xdat = xdat,
            exdat = exdat,
            bws = bws,
            B = plot.errors.boot.num,
            operator = "integral",
            counts.drawer = counts.drawer
          ),
          error = function(e) {
            stop(sprintf("%s ksum helper path failed in compute.bootstrap.errors.dbandwidth (%s)",
                         if (fast.block) plot.errors.boot.method else "inid",
                         conditionMessage(e)),
                 call. = FALSE)
          }
        )
      })
    }

    if (is.null(boot.out))
      stop("no MPI helper path available for this unconditional distribution bootstrap configuration in npRmpi; no serial fallback is permitted", call. = FALSE)

    all.bp <- list()

    if (slice.index > 0 && (bws$xdati$iord[slice.index] || bws$xdati$iuno[slice.index])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- bws$xdati$all.ulev[[slice.index]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- bws$xdati$all.lev[[slice.index]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise",
          warn.coverage = FALSE)
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.conbandwidth =
  function(xdat, ydat,
           exdat, eydat,
           cdf,
           quantreg,
           tau,
           gradients,
           gradient.index,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.conbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.conbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    exdat = toFrame(exdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    tboo =
      if(quantreg) "quant"
      else if (cdf) "dist"
      else "dens"

    if (!identical(tboo, "quant")) {
      .np_plot_reject_wild_unsupervised(plot.errors.boot.method, "conditional density/distribution estimators")
    }

    is.inid = plot.errors.boot.method=="inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))
    fast.inid <- isTRUE(.np_plot_inid_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.inid) &&
      isTRUE(!quantreg) &&
      isTRUE(!gradients) &&
      isTRUE(.np_con_inid_ksum_eligible(bws))
    fast.block <- isTRUE(.np_plot_block_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.block) &&
      isTRUE(!quantreg) &&
      isTRUE(!gradients) &&
      isTRUE(.np_con_inid_ksum_eligible(bws))

    if (is.inid && !isTRUE(fast.inid))
      stop("inid conditional helper unavailable for this configuration in npRmpi; no serial fallback is permitted", call. = FALSE)
    if (is.block && !isTRUE(fast.block))
      stop("fixed/geom conditional helper unavailable for this configuration in npRmpi; no serial fallback is permitted", call. = FALSE)

    boot.out <- NULL
    if (fast.inid || fast.block) {
      counts.drawer <- if (fast.block) {
        .np_block_counts_drawer(
          n = nrow(xdat),
          B = plot.errors.boot.num,
          blocklen = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method
        )
      } else {
        NULL
      }
      boot.out <- .npRmpi_with_local_bootstrap({
        tryCatch(
          .np_inid_boot_from_ksum_conditional(
            xdat = xdat,
            ydat = ydat,
            exdat = exdat,
            eydat = eydat,
            bws = bws,
            B = plot.errors.boot.num,
            cdf = cdf,
            counts.drawer = counts.drawer
          ),
          error = function(e) {
            stop(sprintf("%s ksum helper failed in compute.bootstrap.errors.conbandwidth (%s)",
                         if (fast.block) plot.errors.boot.method else "inid",
                         conditionMessage(e)),
                 call. = FALSE)
          }
        )
      })
    }

    if (is.null(boot.out))
      stop("no MPI helper path available for this conditional bootstrap configuration in npRmpi; no serial fallback is permitted", call. = FALSE)

    all.bp <- list()

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$ydati
      ti <- slice.index - bws$xndim
    }
    
    if (slice.index > 0 && (tdati$iord[ti] || tdati$iuno[ti])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- tdati$all.lev[[ti]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise",
          warn.coverage = FALSE)
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.condbandwidth =
  function(xdat, ydat,
           exdat, eydat,
           cdf,
           quantreg,
           tau,
           gradients,
           gradient.index,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.condbandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.condbandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(exdat)
    )
    exdat = toFrame(exdat)
    boot.err = matrix(data = NA, nrow = dim(exdat)[1], ncol = 3)
    boot.all.err <- NULL

    tboo =
      if(quantreg) "quant"
      else if (cdf) "dist"
      else "dens"

    if (!identical(tboo, "quant")) {
      .np_plot_reject_wild_unsupervised(plot.errors.boot.method, "conditional density/distribution estimators")
    }

    is.inid = plot.errors.boot.method=="inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))
    fast.inid <- isTRUE(.np_plot_inid_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.inid) &&
      isTRUE(!quantreg) &&
      isTRUE(!gradients) &&
      isTRUE(.np_con_inid_ksum_eligible(bws))
    fast.block <- isTRUE(.np_plot_block_fastpath_enabled()) &&
      isTRUE(.npRmpi_plot_inid_ksum_fastpath_enabled()) &&
      isTRUE(is.block) &&
      isTRUE(!quantreg) &&
      isTRUE(!gradients) &&
      isTRUE(.np_con_inid_ksum_eligible(bws))
    if (is.inid && !isTRUE(fast.inid))
      stop("inid conditional helper unavailable for this configuration in npRmpi; no serial fallback is permitted", call. = FALSE)
    if (is.block && !isTRUE(fast.block))
      stop("fixed/geom conditional helper unavailable for this configuration in npRmpi; no serial fallback is permitted", call. = FALSE)

    boot.out <- NULL
    if (fast.inid || fast.block) {
      counts.drawer <- if (fast.block) {
        .np_block_counts_drawer(
          n = nrow(xdat),
          B = plot.errors.boot.num,
          blocklen = plot.errors.boot.blocklen,
          sim = plot.errors.boot.method
        )
      } else {
        NULL
      }
      boot.out <- .npRmpi_with_local_bootstrap({
        tryCatch(
          .np_inid_boot_from_ksum_conditional(
            xdat = xdat,
            ydat = ydat,
            exdat = exdat,
            eydat = eydat,
            bws = bws,
            B = plot.errors.boot.num,
            cdf = cdf,
            counts.drawer = counts.drawer
          ),
          error = function(e) {
            stop(sprintf("%s ksum helper failed in compute.bootstrap.errors.condbandwidth (%s)",
                         if (fast.block) plot.errors.boot.method else "inid",
                         conditionMessage(e)),
                 call. = FALSE)
          }
        )
      })
    }

    if (is.null(boot.out))
      stop("no MPI helper path available for this conditional bootstrap configuration in npRmpi; no serial fallback is permitted", call. = FALSE)

    all.bp <- list()

    if (slice.index <= bws$xndim){
      tdati <- bws$xdati
      ti <- slice.index
    } else {
      tdati <- bws$ydati
      ti <- slice.index - bws$xndim
    }
    
    if (slice.index > 0 && (tdati$iord[ti] || tdati$iuno[ti])){
      boot.frame <- as.data.frame(boot.out$t)
      u.lev <- tdati$all.ulev[[ti]]

      ## if we are bootstrapping a factor, there should be one
      ## set of replications for each level
      stopifnot(length(u.lev)==ncol(boot.frame))
      
      all.bp$stats <- matrix(data = NA, nrow = 5, ncol = length(u.lev))
      all.bp$conf <- matrix(data = NA, nrow = 2, ncol = length(u.lev))

      for (i in seq_along(u.lev)){
        t.bp <- boxplot.stats(boot.frame[,i])
        all.bp$stats[,i] <- t.bp$stats
        all.bp$conf[,i] <- t.bp$conf
        all.bp$out <- c(all.bp$out,t.bp$out)
        all.bp$group <- c(all.bp$group, rep.int(i,length(t.bp$out)))
      }
      all.bp$n <- rep.int(plot.errors.boot.num, length(u.lev))
      all.bp$names <- tdati$all.lev[[ti]]
      rm(boot.frame)
    }

    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise",
          warn.coverage = FALSE)
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = all.bp,
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }

compute.bootstrap.errors.sibandwidth =
  function(xdat, ydat,
           gradients,
           plot.errors.boot.method,
           plot.errors.boot.wild = c("rademacher", "mammen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           ...,
           bws){
    .np_plot_require_bws(bws = bws, where = "compute.bootstrap.errors.sibandwidth")
    prof.ctx <- .npRmpi_profile_bootstrap_begin(
      where = "compute.bootstrap.errors.sibandwidth",
      method = plot.errors.boot.method,
      B = plot.errors.boot.num,
      ntrain = .np_nrows_safe(xdat),
      neval = .np_nrows_safe(xdat)
    )

    boot.err = matrix(data = NA, nrow = nrow(xdat), ncol = 3)
    boot.all.err <- NULL

    is.wild.hat <- .np_plot_is_wild_method(plot.errors.boot.method)
    is.inid <- plot.errors.boot.method=="inid"
    is.block <- is.element(plot.errors.boot.method, c("fixed", "geom"))

    if (is.wild.hat) {
      if (length(plot.errors.boot.wild) > 1L)
        plot.errors.boot.wild <- plot.errors.boot.wild[1L]
      plot.errors.boot.wild <- match.arg(plot.errors.boot.wild, c("mammen", "rademacher"))

      fit.train <- npindex.sibandwidth(
        txdat = xdat,
        tydat = ydat,
        bws = bws,
        gradients = FALSE
      )
      H <- npindexhat(
        bws = bws,
        txdat = xdat,
        exdat = xdat,
        output = "matrix",
        s = if (gradients) 1L else 0L
      )

      t0 <- as.vector(H %*% as.double(ydat))
      eps <- as.double(ydat - as.vector(fit.train$mean))
      n <- length(eps)
      B <- plot.errors.boot.num

      boot.out <- .npRmpi_with_local_bootstrap({
        list(
          t = .np_wild_boot_t(
            H = H,
            fit.mean = as.vector(fit.train$mean),
            residuals = eps,
            B = B,
            wild = plot.errors.boot.wild
          ),
          t0 = t0
        )
      })
    } else if (is.inid) {
      inid.helper.ok <- isTRUE(.np_plot_inid_fastpath_enabled()) &&
        !isTRUE(gradients)
      if (!isTRUE(inid.helper.ok)) {
        stop("inid single-index helper unavailable for this configuration in npRmpi; no serial fallback is permitted", call. = FALSE)
      } else {
        boot.out <- .npRmpi_with_local_bootstrap({
          tryCatch({
            tx.index <- data.frame(index = as.vector(toMatrix(xdat) %*% bws$beta))
            rbw <- .np_indexhat_rbw(bws = bws, idx.train = tx.index)
            .np_inid_boot_from_regression(
              xdat = tx.index,
              exdat = tx.index,
              bws = rbw,
              ydat = ydat,
              B = plot.errors.boot.num
            )
          }, error = function(e) {
            stop(sprintf("inid single-index helper failed in compute.bootstrap.errors.sibandwidth (%s)",
                         conditionMessage(e)),
                 call. = FALSE)
          })
        })
      }
    } else if (is.block) {
      block.helper.ok <- isTRUE(.np_plot_block_fastpath_enabled()) &&
        !isTRUE(gradients)
      if (!isTRUE(block.helper.ok)) {
        stop("fixed/geom single-index helper unavailable for this configuration in npRmpi; no serial fallback is permitted", call. = FALSE)
      } else {
        boot.out <- .npRmpi_with_local_bootstrap({
          tryCatch({
            tx.index <- data.frame(index = as.vector(toMatrix(xdat) %*% bws$beta))
            rbw <- .np_indexhat_rbw(bws = bws, idx.train = tx.index)
            .np_inid_boot_from_regression(
              xdat = tx.index,
              exdat = tx.index,
              bws = rbw,
              ydat = ydat,
              B = plot.errors.boot.num,
              counts.drawer = .np_block_counts_drawer(
                n = nrow(tx.index),
                B = plot.errors.boot.num,
                blocklen = plot.errors.boot.blocklen,
                sim = plot.errors.boot.method
              )
            )
          }, error = function(e) {
            stop(sprintf("%s single-index helper failed in compute.bootstrap.errors.sibandwidth (%s)",
                         plot.errors.boot.method,
                         conditionMessage(e)),
                 call. = FALSE)
          })
        })
      }
    }

    if (is.null(boot.out))
      stop("no MPI helper path available for this single-index bootstrap configuration in npRmpi; no serial fallback is permitted", call. = FALSE)
    
    if (plot.errors.type == "pmzsd") {
      boot.err[,1:2] = qnorm(plot.errors.alpha/2, lower.tail = FALSE)*sqrt(diag(cov(boot.out$t)))
    }
    else if (plot.errors.type %in% c("pointwise", "bonferroni", "simultaneous", "all")) {
      if (plot.errors.type == "all") {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "pointwise",
          warn.coverage = FALSE)
        boot.all.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = "all")
        boot.all.err <- lapply(boot.all.bounds, function(bb) {
          cbind(boot.out$t0 - bb[,1], bb[,2] - boot.out$t0)
        })
      } else {
        boot.bounds <- compute.bootstrap.quantile.bounds(
          boot.t = boot.out$t,
          alpha = plot.errors.alpha,
          band.type = plot.errors.type)
      }
      boot.err[,1] <- boot.out$t0 - boot.bounds[,1]
      boot.err[,2] <- boot.bounds[,2] - boot.out$t0
    }
    if (plot.errors.center == "bias-corrected")
      boot.err[,3] <- 2*boot.out$t0-colMeans(boot.out$t)
    .npRmpi_profile_finalize_bootstrap(
      boot.err = boot.err,
      bxp = list(),
      boot.all.err = boot.all.err,
      ctx = prof.ctx
    )
  }


uocquantile <- function(x, prob) {
  if(anyNA(prob)) stop("'prob' contains missing values")
  if(any(prob < 0 | prob > 1, na.rm = TRUE)) stop("'prob' outside [0,1]")
  if(anyNA(x)) stop("missing values and NaN's not allowed")
  if (is.ordered(x)){
    x <- droplevels(x)
    tq = unclass(table(x))
    tq = tq / sum(tq)
    tq[length(tq)] <- 1.0
    bscape <- levels(x)
    tq <- cumsum(tq)
    j <- sapply(prob, function(p){ which(tq >= p)[1] })
    bscape[j]
  } else if (is.factor(x)) {
    ## just returns mode
    x <- droplevels(x)
    tq = unclass(table(x))
    j = which(tq == max(tq))[1]
    levels(x)[j]
  } else {
    quantile(x, probs = prob)
  }
}


trim.quantiles <- function(dat, trim){
  if (sign(trim) == sign(-1)){
    trim <- abs(trim)
    tq <- quantile(dat, probs = c(0.0, 0.0+trim, 1.0-trim,1.0))
    tq <- c(2.0*tq[1]-tq[2], 2.0*tq[4]-tq[3])
  }
  else {
    tq <- quantile(dat, probs = c(0.0+trim, 1.0-trim))
  }
  tq
}
