\name{se}
\alias{se}
\title{Extract Standard Errors}
\usage{
se(x)
}
\arguments{
  \item{x}{an object for which the extraction of standard errors is
    meaningful.}
}
\value{
  Standard errors extracted from the model object \code{x}.
}
\description{
  \code{se} is a generic function which extracts standard errors 
  from objects. 
}
\details{
  This function provides a generic interface for extraction of standard
  errors from objects.
}
\seealso{
  \code{\link{fitted}}, \code{\link{residuals}}, \code{\link{coef}},
  and \code{\link{gradients}}, for related methods;
  \code{\link{npRmpi}} for supported objects. 
}

\author{
  Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\section{Note}{
  This method currently only supports objects from the \code{\link{npRmpi}} library.
}
\examples{
\dontrun{
## Not run in checks: excluded to keep MPI examples stable and check times short.
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

## Start npRmpi for interactive execution. If slaves are already running and
## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
## reuse the existing pool instead of respawning. To change the number of
## slaves, call `npRmpi.stop(force=TRUE)` then restart.
npRmpi.start(nslaves=1)

set.seed(42)

x <- rnorm(10)
mpi.bcast.Robj2slave(x)

mpi.bcast.cmd(bw <- npudensbw(~x),
              caller.execute=TRUE)

mpi.bcast.cmd(fhat <- npudens(bw),
              caller.execute=TRUE)

se(fhat)

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

## Note: on some systems (notably macOS+MPICH), repeatedly spawning and
## tearing down slaves in the same R session can lead to hangs/crashes.
## npRmpi may therefore keep slave daemons alive by default and
## `npRmpi.stop()` performs a "soft close". Use `force=TRUE` to
## actually shut down the slaves.
##
## You can disable reuse via `options(npRmpi.reuse.slaves=FALSE)` or by
## setting the environment variable `NP_RMPI_NO_REUSE_SLAVES=1` before
## loading the package.

npRmpi.stop()               ## soft close (may keep slaves alive)
## npRmpi.stop(force=TRUE)  ## hard close

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)
} % enddontrun
}
\keyword{nonparametric}
