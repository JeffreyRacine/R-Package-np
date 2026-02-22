# npRmpi attach-mode startup profile template
#
# Usage (batch/cluster):
#   1) copy this file to your working directory as .Rprofile, or set
#      R_PROFILE to this file path.
#   2) launch with mpiexec/mpirun, e.g.
#      mpiexec -n 6 R CMD BATCH --no-save script.R
#
# Non-master ranks enter a receive/eval loop and execute commands sent by
# mpi.bcast.cmd(...). Rank 0 continues script execution.

.nonblock <- TRUE
.sleep <- 0.1
quiet <- FALSE

if (!invisible(library(npRmpi, logical.return = TRUE))) {
  warning("npRmpi cannot be loaded")
  q(save = "no")
}

options(error = quote(assign(".mpi.err", FALSE, envir = .GlobalEnv)))

if (mpi.comm.size(0) > 1)
  invisible(mpi.comm.dup(0, 1))

if (mpi.comm.rank(0) > 0) {
  options(echo = FALSE)
  .comm <- 1
  mpi.barrier(0)
  repeat {
    tmp.message <- mpi.bcast.cmd(rank = 0, comm = .comm,
                                 nonblock = .nonblock, sleep = .sleep)
    if (is.character(tmp.message) && identical(tmp.message, "kaerb"))
      break
    try(eval(tmp.message, envir = .GlobalEnv), TRUE)
  }
  mpi.comm.free(.comm)
  mpi.quit()
}

if (mpi.comm.rank(0) == 0) {
  mpi.barrier(0)
  if (mpi.comm.size(0) > 1 && !quiet)
    slave.hostinfo(1)
}

.Last <- function() {
  if (is.loaded("mpi_initialize")) {
    if (mpi.comm.size(1) > 1) {
      print("Please use mpi.close.Rslaves() to close slaves")
      mpi.close.Rslaves(comm = 1)
    }
  }
  if (is.loaded("mpi_initialize"))
    .Call("mpi_finalize", PACKAGE = "npRmpi")
}
