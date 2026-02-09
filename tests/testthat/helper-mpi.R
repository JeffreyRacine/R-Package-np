# Helper for MPI tests
spawn_mpi_slaves <- function(n=1) {
  if (!requireNamespace("Rmpi", quietly = TRUE)) {
    return(FALSE)
  }

  # Prefer the higher-level helper if available.
  if (exists("npRmpi.start", mode="function")) {
    ok <- try({
      npRmpi.start(nslaves=n, needlog=FALSE)
      TRUE
    }, silent=TRUE)
    return(isTRUE(ok) && isTRUE(try(mpi.comm.size(1) > 1, silent=TRUE)))
  }

  # Fallback to the legacy pattern.
  if (isTRUE(try(mpi.comm.size(1) > 1, silent=TRUE))) {
    return(TRUE)
  }

  spawn_status <- try(mpi.spawn.Rslaves(nslaves=n, quiet=FALSE), silent=TRUE)
  if (inherits(spawn_status, "try-error") || !isTRUE(try(mpi.comm.size(1) > 1, silent=TRUE))) {
    return(FALSE)
  }

  ok <- try(mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE), silent=TRUE)
  !inherits(ok, "try-error")
}

close_mpi_slaves <- function(force=TRUE) {
  if (requireNamespace("Rmpi", quietly = TRUE)) {
    if (isTRUE(try(mpi.comm.size(1) > 1, silent=TRUE))) {
      if (exists("npRmpi.stop", mode="function")) {
        try(npRmpi.stop(force=force), silent=TRUE)
      } else {
        try(mpi.close.Rslaves(force=force), silent=TRUE)
      }
    }
  }
}
