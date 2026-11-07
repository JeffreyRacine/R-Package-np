# Run in a supervised subprocess; no worker errors are caught inside the command.
# This transport test uses a synthetic package-internal agreed-refusal condition.
suppressPackageStartupMessages(library(npRmpi))
a <- commandArgs(TRUE); mode <- a[1L]; workers <- as.integer(a[2L])
if (mode == "session") {
  npRmpi.init(nslaves=workers, quiet=TRUE); master <- TRUE
} else if (mode == "attach") {
  master <- isTRUE(npRmpi.init(mode="attach", quiet=TRUE))
} else master <- mpi.comm.rank(1L) == 0L
if (master) {
  if (mode == "profile") mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)
  mpi.bcast.cmd({
    options(np.messages=FALSE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(73L)
    before <- .Random.seed; cleaned <- FALSE; after <- FALSE
    refusal <- structure(list(message="agreed transport fixture",call=NULL),
                         class=c("np_rng_scope_refusal","error","condition"))
    reject <- function() { on.exit(cleaned <<- TRUE); stop(refusal) }
  }, caller.execute=TRUE)
  got <- tryCatch(mpi.bcast.cmd({ reject(); after <- TRUE }, caller.execute=TRUE), error=identity)
  stopifnot(identical(got,refusal), cleaned, !after, identical(before,.Random.seed))
  remote <- get("mpi.remote.exec",asNamespace("npRmpi"))(
    list(rank=mpi.comm.rank(1L),cleaned=cleaned,after=after,seed=identical(before,.Random.seed)))
  stopifnot(length(remote)==workers,
            all(vapply(remote,`[[`,logical(1),"cleaned")),
            !any(vapply(remote,`[[`,logical(1),"after")),
            all(vapply(remote,`[[`,logical(1),"seed")))
  # Receiver do.call construction and subsequent numerical/collector reuse.
  got <- tryCatch(mpi.bcast.cmd(stop,refusal,caller.execute=TRUE),error=identity)
  stopifnot(identical(got,refusal))
  value <- mpi.bcast.cmd(2+3,caller.execute=TRUE); stopifnot(identical(value,5))
  returned <- mpi.bcast.cmd(simpleError("returned"),caller.execute=TRUE)
  stopifnot(inherits(returned,"simpleError"))
  stopifnot(length(get("mpi.remote.exec",asNamespace("npRmpi"))(list(rank=mpi.comm.rank(1L))))==workers)
  cat("R27_MANUAL_TRANSPORT_PASS", mode, workers, "\n")
  if(mode=="session")npRmpi.quit()
  else if(mode=="attach"){npRmpi.quit(mode="attach");mpi.quit()}
  else mpi.bcast.cmd(mpi.quit(),caller.execute=TRUE)
}
