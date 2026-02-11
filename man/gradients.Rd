\name{gradients}
\title{Extract Gradients}
\alias{gradients}
\alias{gradients.condensity}
\alias{gradients.condistribution}
\alias{gradients.npregression}
\alias{gradients.qregression}  
\alias{gradients.singleindex}
\usage{
gradients(x, \dots)

\method{gradients}{condensity}(x, errors = FALSE, \dots)

\method{gradients}{condistribution}(x, errors = FALSE, \dots)

\method{gradients}{npregression}(x, errors = FALSE, \dots)

\method{gradients}{qregression}(x, errors = FALSE, \dots)

\method{gradients}{singleindex}(x, errors = FALSE, \dots)
}

\arguments{
  \item{x}{an object for which the extraction of gradients is
    meaningful.}
  \item{\dots}{other arguments.}
  \item{errors}{ a logical value specifying whether or not standard
    errors of gradients are desired. Defaults to \code{FALSE}.}
}
\value{
  Gradients extracted from the model object \code{x}.
}
\description{
  \code{gradients} is a generic function which extracts gradients 
  from objects. 
}
\details{
  This function provides a generic interface for extraction of gradients
  from objects.
}
\seealso{
  \code{\link{fitted}}, \code{\link{residuals}}, \code{\link{coef}},
  and \code{\link{se}}, for related methods;
  \code{\link{npRmpi}} for supported objects. 
}
\references{
  See the references for the method being interrogated via
  \code{\link{gradients}} in the appropriate help file. For example, for
  the particulars of the gradients for nonparametric regression see the
  references in \code{\link{npreg}}
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

x <- runif(10)
y <- x + rnorm(10, sd = 0.1)
mydat <- data.frame(x,y)
rm(x,y)

mpi.bcast.Robj2slave(mydat)

mpi.bcast.cmd(model <- npreg(y~x, data=mydat, gradients=TRUE),
              caller.execute=TRUE)

gradients(model)

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
