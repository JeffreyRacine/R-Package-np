# Helper for MPI tests
spawn_mpi_slaves <- function(n=1) {
  if (!requireNamespace("Rmpi", quietly = TRUE)) {
    return(FALSE)
  }
  
  # Check if slaves already spawned
  if (mpi.comm.size(1) > 1) {
    return(TRUE)
  }

  spawn_status <- try(mpi.spawn.Rslaves(nslaves=n, quiet=FALSE), silent=FALSE)
  
  if (inherits(spawn_status, "try-error") || mpi.comm.size(1) < 2) {
    return(FALSE)
  }
  
  # Initialize
  mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)
  
  return(TRUE)
}

close_mpi_slaves <- function() {
  if (requireNamespace("Rmpi", quietly = TRUE)) {
    if (mpi.comm.size(0) > 1) {
      try(mpi.close.Rslaves(), silent=TRUE)
    }
  }
}
