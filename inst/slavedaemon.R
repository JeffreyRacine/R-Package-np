if (!library(npRmpi,logical.return = TRUE)){
    warning("npRmpi cannot be loaded")
    q(save = "no")
}
options(error=quote(assign(".mpi.err", TRUE, envir = .GlobalEnv)))
.comm <- 1
.intercomm <- 2
invisible(mpi.comm.get.parent(.intercomm))
invisible(mpi.intercomm.merge(.intercomm,1,.comm))
invisible(mpi.comm.set.errhandler(.comm))
mpi.hostinfo(.comm)
invisible(mpi.comm.disconnect(.intercomm))
.nonblock <- as.logical(mpi.bcast(integer(1),type=1,rank=0,comm=.comm))
.sleep <- mpi.bcast(double(1),type=2,rank=0,comm=.comm)
repeat 
    try(eval(mpi.bcast.cmd(rank=0,comm=.comm, nonblock=.nonblock, sleep=.sleep),envir=.GlobalEnv),TRUE)
print("Done")
invisible(mpi.comm.disconnect(.comm))
invisible(mpi.comm.set.errhandler(0))
mpi.quit()
