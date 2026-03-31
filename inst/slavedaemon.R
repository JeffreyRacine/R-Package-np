if (!library(npRmpi,logical.return = TRUE)){
    warning("npRmpi cannot be loaded")
    q(save = "no")
}
.npRmpi_ns <- asNamespace("npRmpi")
.mpi.comm.get.parent <- get("mpi.comm.get.parent", envir = .npRmpi_ns, inherits = FALSE)
.mpi.intercomm.merge <- get("mpi.intercomm.merge", envir = .npRmpi_ns, inherits = FALSE)
.mpi.comm.set.errhandler <- get("mpi.comm.set.errhandler", envir = .npRmpi_ns, inherits = FALSE)
.mpi.hostinfo <- get("mpi.hostinfo", envir = .npRmpi_ns, inherits = FALSE)
.mpi.comm.disconnect <- get("mpi.comm.disconnect", envir = .npRmpi_ns, inherits = FALSE)
.mpi.bcast <- get("mpi.bcast", envir = .npRmpi_ns, inherits = FALSE)
.npRmpi_worker_loop <- get(".npRmpi_worker_loop", envir = .npRmpi_ns, inherits = FALSE)
.mpi.comm.size <- get("mpi.comm.size", envir = .npRmpi_ns, inherits = FALSE)
.mpi.comm.free <- get("mpi.comm.free", envir = .npRmpi_ns, inherits = FALSE)
.mpi.quit <- get("mpi.quit", envir = .npRmpi_ns, inherits = FALSE)
options(error=quote(assign(".mpi.err", TRUE, envir = .GlobalEnv)))
.comm <- 1
.intercomm <- 2
invisible(.mpi.comm.get.parent(.intercomm))
invisible(.mpi.intercomm.merge(.intercomm,1,.comm))
invisible(.mpi.comm.set.errhandler(.comm))
.mpi.hostinfo(.comm)
invisible(.mpi.comm.disconnect(.intercomm))
.nonblock <- as.logical(.mpi.bcast(integer(1),type=1,rank=0,comm=.comm))
.sleep <- .mpi.bcast(double(1),type=2,rank=0,comm=.comm)
.recv.timeout <- suppressWarnings(as.numeric(Sys.getenv(
  "NP_RMPI_SESSION_RECV_TIMEOUT_SEC", "0"
)))
if (!is.finite(.recv.timeout) || .recv.timeout <= 0)
  .recv.timeout <- 0
.npRmpi_worker_loop(
  comm = .comm,
  nonblock = .nonblock,
  sleep = .sleep,
  recv.timeout = .recv.timeout,
  loop.label = "spawn slave",
  timeout.remediation = "verify FI_TCP_IFACE and relaunch with a fresh spawn session."
)
print("Done")
if (.mpi.comm.size(.comm) > 0) {
    # `.comm` is an intracommunicator (created via `MPI_Intercomm_merge()`).
    # Prefer `MPI_Comm_free()`; `MPI_Comm_disconnect()` can destabilize
    # subsequent spawn/merge cycles with MPICH.
    if (is.loaded("mpi_comm_free")) {
        invisible(.mpi.comm.free(.comm))
    } else if (is.loaded("mpi_comm_disconnect")) {
        invisible(.mpi.comm.disconnect(.comm))
    } else {
        invisible(.mpi.comm.free(.comm))
    }
}
invisible(.mpi.comm.set.errhandler(0))
.mpi.quit()
