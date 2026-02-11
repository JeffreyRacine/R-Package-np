\name{npseed}
\title{Set Random Seed}
\alias{npseed}

\usage{
npseed(seed)
}

\arguments{
  \item{seed}{an integer seed for the random number generator.}
}
\value{
  None.
}
\description{
  \code{npseed} is a function which sets the random seed in the
  \code{\link{npRmpi}} C backend, resetting the random number generator.
}

\details{
  \code{npseed} provides an interface for setting the random seed (and
  resetting the random number generator) used
  by \code{\link{npRmpi}}. The random number generator is used during the
  bandwidth search procedure to set the search starting point, and in
  subsequent searches when using multistarting, to avoid being trapped
  in local minima if the objective function is not globally concave.

  Calling \code{npseed} will only affect the numerical search if it is
  performed by the C backend. The affected functions include:
  \code{\link{npudensbw}}, \code{\link{npcdensbw}},
  \code{\link{npregbw}}, \code{\link{npplregbw}}, \code{\link{npqreg}},
  \code{\link{npcmstest}} (via \code{\link{npregbw}}),
  \code{\link{npqcmstest}} (via \code{\link{npregbw}}),
  \code{\link{npsigtest}} (via \code{\link{npregbw}}).
}

\seealso{
  \code{\link{set.seed}}
}

\references{
    Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
      Theory and Practice,} Princeton University Press.
}

\author{ Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey
  S. Racine \email{racinej@mcmaster.ca}
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

mpi.bcast.cmd(npseed(712),
              caller.execute=TRUE)

x <- runif(10)
y <- x + rnorm(10, sd = 0.1)
mydat <- data.frame(x,y)
rm(x,y)

mpi.bcast.Robj2slave(mydat)

mpi.bcast.cmd(bw <- npregbw(y~x, data=mydat),
              caller.execute=TRUE)

summary(bw)

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
