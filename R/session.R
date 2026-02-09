npRmpi.start <- function(..., nslaves=1, comm=1){
  # Ensure a slave pool is available, then initialize npRmpi on all ranks.
  mpi.spawn.Rslaves(..., nslaves=nslaves, comm=comm)
  mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE, comm=comm)
  invisible(TRUE)
}

npRmpi.stop <- function(force=FALSE, dellog=TRUE, comm=1){
  # Idempotent stop: if no slaves are running, return silently.
  size <- try(mpi.comm.size(comm), silent=TRUE)
  if (inherits(size, "try-error") || is.na(size) || size < 2)
    return(invisible(FALSE))
  mpi.close.Rslaves(dellog=dellog, comm=comm, force=force)
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
