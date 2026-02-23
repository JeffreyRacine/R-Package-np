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

.npRmpi_safe_int <- function(expr) {
  tryCatch(as.integer(expr), error = function(e) NA_integer_)
}

.npRmpi_safe <- function(expr, fallback = NULL) {
  tryCatch(expr, error = function(e) fallback)
}

.npRmpi_session_ensure_comm <- function(comm = 1L) {
  size <- .npRmpi_safe_int(mpi.comm.size(comm))
  if (!is.na(size))
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

.npRmpi_session_attach_worker_loop <- function(comm = 1L,
                                               nonblock = TRUE,
                                               sleep = 0.1) {
  options(echo = FALSE)
  repeat {
    msg <- mpi.bcast.cmd(rank = 0, comm = comm, nonblock = nonblock, sleep = sleep)
    if (is.character(msg) && identical(msg, "kaerb"))
      break
    tryCatch(.npRmpi_eval_scmd(msg, envir = .GlobalEnv), error = function(e) invisible(e))
  }
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
                         np.messages = FALSE,
                         nonblock = TRUE,
                         sleep = 0.1,
                         quiet = FALSE) {
  if ("package:Rmpi" %in% search() &&
      !isTRUE(getOption("npRmpi.allow.attached.Rmpi", FALSE))) {
    stop(
      "package 'Rmpi' is attached. Detach it and use npRmpi APIs directly, ",
      "or set options(npRmpi.allow.attached.Rmpi=TRUE) to bypass this guard."
    )
  }

  mode <- match.arg(mode)
  autodispatch.option.sync <- match.arg(autodispatch.option.sync)
  world.size <- .npRmpi_safe_int(mpi.comm.size(0))
  world.size <- if (is.na(world.size)) 1L else as.integer(world.size)

  if (identical(mode, "auto")) {
    mode <- if (world.size > 1L) "attach" else "spawn"
  }

  .npRmpi_session_apply_options(
    autodispatch = autodispatch,
    np.messages = np.messages,
    autodispatch.verify.options = autodispatch.verify.options,
    autodispatch.option.sync = autodispatch.option.sync
  )

  if (identical(mode, "spawn")) {
    if (.npRmpi_session_has_active_pool(comm = comm)) {
      if (!quiet)
        slave.hostinfo(comm)
      return(invisible(TRUE))
    }
    mpi.spawn.Rslaves(..., nslaves = nslaves, comm = comm, quiet = quiet, nonblock = nonblock, sleep = sleep)
    mpi.bcast.cmd(np.mpi.initialize(), caller.execute = TRUE, comm = comm)
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

  rank <- .npRmpi_safe_int(mpi.comm.rank(comm))
  rank <- if (is.na(rank)) 0L else as.integer(rank)

  if (rank == 0L) {
    if (!quiet) slave.hostinfo(comm)
    return(invisible(TRUE))
  }

  .npRmpi_session_attach_worker_loop(comm = comm, nonblock = nonblock, sleep = sleep)
}

npRmpi.quit <- function(force = FALSE,
                        dellog = TRUE,
                        comm = 1,
                        mode = c("auto", "spawn", "attach")) {
  mode <- match.arg(mode)
  size.comm <- .npRmpi_safe_int(mpi.comm.size(comm))
  size.comm <- if (is.na(size.comm)) 0L else as.integer(size.comm)
  size.world <- .npRmpi_safe_int(mpi.comm.size(0))
  size.world <- if (is.na(size.world)) 1L else as.integer(size.world)

  if (identical(mode, "auto")) {
    mode <- if (size.world > 1L && size.comm > 1L) "attach" else "spawn"
  }

  if (size.comm < 2L)
    return(invisible(FALSE))

  if (identical(mode, "attach")) {
    comm <- 1L
    mpi.bcast.cmd(cmd = "kaerb", rank = 0, comm = comm, caller.execute = FALSE)
    if (comm != 0L) {
      .npRmpi_safe(mpi.comm.free(comm), fallback = NULL)
    }
    return(invisible(TRUE))
  }

  mpi.close.Rslaves(dellog = dellog, comm = comm, force = force)
  invisible(TRUE)
}

npRmpi.session.info <- function(comm=1){
  np_ver <- .npRmpi_safe(utils::packageVersion("npRmpi"), fallback = NA)
  rmpi_ver <- .npRmpi_safe(utils::packageVersion("Rmpi"), fallback = NA)
  reuse <- isTRUE(getOption("npRmpi.reuse.slaves", FALSE))
  env_no_reuse <- Sys.getenv("NP_RMPI_NO_REUSE_SLAVES", unset="")

  mpi_ver <- .npRmpi_safe({
    if (requireNamespace("Rmpi", quietly = TRUE) &&
        exists("mpi.get.version", envir = asNamespace("Rmpi"), inherits = FALSE)) {
      get("mpi.get.version", envir = asNamespace("Rmpi"))()
    } else {
      NA
    }
  }, fallback = NA)
  comm_size <- .npRmpi_safe_int(mpi.comm.size(comm))
  comm_rank <- .npRmpi_safe_int(mpi.comm.rank(comm))
  proc <- .npRmpi_safe(mpi.get.processor.name(), fallback = NA_character_)

  info <- list(
    npRmpi = if (is.na(np_ver)[1L]) NA_character_ else as.character(np_ver),
    Rmpi = if (is.na(rmpi_ver)[1L]) NA_character_ else as.character(rmpi_ver),
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
    processor = proc
  )

  cat("npRmpi session info\n")
  cat("  npRmpi:", info$npRmpi, "\n")
  cat("  Rmpi  :", info$Rmpi, "\n")
  cat("  OS    :", info$sysname, info$release, "|", info$platform, "\n")
  cat("  reuse :", info$reuse_slaves, "(NP_RMPI_NO_REUSE_SLAVES=", info$NP_RMPI_NO_REUSE_SLAVES, ")\n", sep="")
  cat("  comm  :", info$comm, "rank", info$comm_rank, "of", info$comm_size, "(nslaves=", info$nslaves, ")\n", sep=" ")
  if (!is.na(info$processor))
    cat("  host  :", info$processor, "\n")
  if (!is.na(info$mpi_version))
    cat("  mpi   :", paste(info$mpi_version, collapse=" "), "\n")

  invisible(info)
}
