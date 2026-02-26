if (!library(npRmpi,logical.return = TRUE)){
    warning("npRmpi cannot be loaded")
    q(save = "no")
}
options(error=quote(assign(".mpi.err", TRUE, envir = .GlobalEnv)))
.comm <- 1
.intercomm <- 2
invisible(npRmpi:::mpi.comm.get.parent(.intercomm))
invisible(npRmpi:::mpi.intercomm.merge(.intercomm,1,.comm))
invisible(npRmpi:::mpi.comm.set.errhandler(.comm))
npRmpi:::mpi.hostinfo(.comm)
invisible(npRmpi:::mpi.comm.disconnect(.intercomm))
.nonblock <- as.logical(npRmpi:::mpi.bcast(integer(1),type=1,rank=0,comm=.comm))
.sleep <- npRmpi:::mpi.bcast(double(1),type=2,rank=0,comm=.comm)
repeat {
	tmp.message=npRmpi:::mpi.bcast.cmd(rank=0,comm=.comm, nonblock=.nonblock, sleep=.sleep)
	if (is.character(tmp.message) && tmp.message =="kaerb")
		break
    res <- try(eval(tmp.message,envir=.GlobalEnv), silent=TRUE)
    if (inherits(res, "try-error")) {
        cmd <- paste(utils::capture.output(print(tmp.message)), collapse = " ")
        msg <- as.character(res)
        base::cat(sprintf("\n[spawn slave rank %d] CMD: %s\n[spawn slave rank %d] ERROR: %s\n",
                          npRmpi:::mpi.comm.rank(.comm), cmd, npRmpi:::mpi.comm.rank(.comm),
                          paste(msg, collapse = " ")),
                  file = stderr())
        flush(stderr())
        base::stop(msg)
    }
}
print("Done")
if (npRmpi:::mpi.comm.size(.comm) > 0) {
    # `.comm` is an intracommunicator (created via `MPI_Intercomm_merge()`).
    # Prefer `MPI_Comm_free()`; `MPI_Comm_disconnect()` can destabilize
    # subsequent spawn/merge cycles with MPICH.
    if (is.loaded("mpi_comm_free")) {
        invisible(npRmpi:::mpi.comm.free(.comm))
    } else if (is.loaded("mpi_comm_disconnect")) {
        invisible(npRmpi:::mpi.comm.disconnect(.comm))
    } else {
        invisible(npRmpi:::mpi.comm.free(.comm))
    }
}
invisible(npRmpi:::mpi.comm.set.errhandler(0))
npRmpi:::mpi.quit()
