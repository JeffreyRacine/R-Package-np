.npRmpi_session_apply_options <- function(autodispatch = NULL,
                                          np.messages = NULL) {
  if (!is.null(autodispatch))
    options(npRmpi.autodispatch = isTRUE(autodispatch))
  if (!is.null(np.messages))
    options(np.messages = isTRUE(np.messages))
  invisible(TRUE)
}

.npRmpi_session_ensure_comm <- function(comm = 1L) {
  size <- try(mpi.comm.size(comm), silent = TRUE)
  if (!inherits(size, "try-error") && !is.na(size))
    return(invisible(as.integer(size)))
  invisible(mpi.comm.dup(0, comm))
}

.npRmpi_session_attach_worker_loop <- function(comm = 1L,
                                               nonblock = TRUE,
                                               sleep = 0.1) {
  options(echo = FALSE)
  repeat {
    msg <- mpi.bcast.cmd(rank = 0, comm = comm, nonblock = nonblock, sleep = sleep)
    if (is.character(msg) && identical(msg, "kaerb"))
      break
    try(eval(msg, envir = .GlobalEnv), silent = TRUE)
  }
  try(if (comm != 0L) mpi.comm.free(comm), silent = TRUE)
  mpi.quit()
  invisible(FALSE)
}

npRmpi.init <- function(...,
                         nslaves = 1,
                         comm = 1,
                         mode = c("auto", "spawn", "attach"),
                         autodispatch = TRUE,
                         np.messages = FALSE,
                         nonblock = TRUE,
                         sleep = 0.1,
                         quiet = FALSE) {
  mode <- match.arg(mode)
  world.size <- try(mpi.comm.size(0), silent = TRUE)
  world.size <- if (inherits(world.size, "try-error") || is.na(world.size)) 1L else as.integer(world.size)

  if (identical(mode, "auto")) {
    mode <- if (world.size > 1L) "attach" else "spawn"
  }

  .npRmpi_session_apply_options(autodispatch = autodispatch, np.messages = np.messages)

  if (identical(mode, "spawn")) {
    mpi.spawn.Rslaves(..., nslaves = nslaves, comm = comm, quiet = quiet, nonblock = nonblock, sleep = sleep)
    mpi.bcast.cmd(np.mpi.initialize(), caller.execute = TRUE, comm = comm)
    return(invisible(TRUE))
  }

  if (world.size < 2L)
    stop("attach mode requires a pre-launched MPI world (e.g. mpiexec -n <master+slaves>)")

  .npRmpi_session_ensure_comm(comm = comm)
  np.mpi.initialize()
  mpi.barrier(0)

  rank <- try(mpi.comm.rank(comm), silent = TRUE)
  rank <- if (inherits(rank, "try-error") || is.na(rank)) 0L else as.integer(rank)

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
  size.comm <- try(mpi.comm.size(comm), silent = TRUE)
  size.comm <- if (inherits(size.comm, "try-error") || is.na(size.comm)) 0L else as.integer(size.comm)
  size.world <- try(mpi.comm.size(0), silent = TRUE)
  size.world <- if (inherits(size.world, "try-error") || is.na(size.world)) 1L else as.integer(size.world)

  if (identical(mode, "auto")) {
    mode <- if (size.world > 1L && size.comm > 1L) "attach" else "spawn"
  }

  if (size.comm < 2L)
    return(invisible(FALSE))

  if (identical(mode, "attach")) {
    mpi.bcast.cmd(cmd = "kaerb", rank = 0, comm = comm, caller.execute = FALSE)
    if (comm != 0L) {
      try(mpi.comm.free(comm), silent = TRUE)
    }
    return(invisible(TRUE))
  }

  mpi.close.Rslaves(dellog = dellog, comm = comm, force = force)
  invisible(TRUE)
}

npRmpi.session.info <- function(comm=1){
  np_ver <- try(utils::packageVersion("npRmpi"), silent=TRUE)
  rmpi_ver <- try(utils::packageVersion("Rmpi"), silent=TRUE)
  reuse <- isTRUE(getOption("npRmpi.reuse.slaves", FALSE))
  env_no_reuse <- Sys.getenv("NP_RMPI_NO_REUSE_SLAVES", unset="")

  mpi_ver <- try({
    if (requireNamespace("Rmpi", quietly = TRUE) &&
        exists("mpi.get.version", envir = asNamespace("Rmpi"), inherits = FALSE)) {
      get("mpi.get.version", envir = asNamespace("Rmpi"))()
    } else {
      NA
    }
  }, silent = TRUE)
  comm_size <- try(mpi.comm.size(comm), silent=TRUE)
  comm_rank <- try(mpi.comm.rank(comm), silent=TRUE)
  proc <- try(mpi.get.processor.name(), silent=TRUE)

  if (inherits(comm_size, "try-error")) comm_size <- NA_integer_
  if (inherits(comm_rank, "try-error")) comm_rank <- NA_integer_

  info <- list(
    npRmpi = if (inherits(np_ver, "try-error")) NA_character_ else as.character(np_ver),
    Rmpi = if (inherits(rmpi_ver, "try-error")) NA_character_ else as.character(rmpi_ver),
    platform = R.version$platform,
    sysname = Sys.info()[["sysname"]],
    release = Sys.info()[["release"]],
    reuse_slaves = reuse,
    NP_RMPI_NO_REUSE_SLAVES = env_no_reuse,
    mpi_version = if (inherits(mpi_ver, "try-error")) NA else mpi_ver,
    comm = comm,
    comm_rank = comm_rank,
    comm_size = comm_size,
    nslaves = if (is.na(comm_size)) NA_integer_ else max(as.integer(comm_size) - 1L, 0L),
    processor = if (inherits(proc, "try-error")) NA_character_ else proc
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
