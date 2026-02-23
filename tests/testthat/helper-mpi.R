# Helper for MPI tests
.mpi_pool_active <- function() {
  if (!requireNamespace("Rmpi", quietly = TRUE))
    return(FALSE)
  if (!isTRUE(getOption("npRmpi.mpi.initialized", FALSE)))
    return(FALSE)
  isTRUE(try(mpi.comm.size(1) > 1, silent = TRUE))
}

spawn_mpi_slaves <- function(n=1) {
  if (!requireNamespace("Rmpi", quietly = TRUE)) {
    return(FALSE)
  }

  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)

  # Prefer the higher-level helper if available.
  if (exists("npRmpi.init", mode="function")) {
    ok <- try({
      npRmpi.init(nslaves=n, needlog=FALSE)
      TRUE
    }, silent=TRUE)
    return(isTRUE(ok) && .mpi_pool_active())
  }

  # Fallback to the legacy pattern.
  if (.mpi_pool_active()) {
    return(TRUE)
  }

  spawn_status <- try(mpi.spawn.Rslaves(nslaves=n, quiet=FALSE), silent=TRUE)
  if (inherits(spawn_status, "try-error") || !isTRUE(try(mpi.comm.size(1) > 1, silent=TRUE))) {
    return(FALSE)
  }

  ok <- try(np.mpi.initialize(), silent=TRUE)
  !inherits(ok, "try-error")
}

close_mpi_slaves <- function(force=TRUE) {
  if (!.mpi_pool_active())
    return(invisible())

  if (exists("npRmpi.quit", mode="function")) {
    try(npRmpi.quit(force=force), silent=TRUE)
  } else {
    try(mpi.close.Rslaves(force=force), silent=TRUE)
  }
}
