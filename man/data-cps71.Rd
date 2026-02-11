% $ID$
\name{cps71}
\docType{data}
\alias{cps71}
\title{ Canadian High School Graduate Earnings  }
\description{
	Canadian cross-section wage data consisting of a random sample taken
  from the 1971 Canadian Census Public Use Tapes for male individuals
  having common education (grade 13). There are 205 observations in total.
}
\usage{data("cps71")}
\format{ A data frame with 2 columns, and 205 rows.
  \describe{
    \item{logwage}{ the first column, of type \code{numeric}}
    \item{age}{ the second column, of type \code{integer}}
  }
}
\source{ Aman Ullah }
\examples{
\dontrun{
## Not run in checks: excluded to keep MPI examples stable and check times short.
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave). Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

## Start npRmpi for interactive execution. If slaves are already running and
## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
## reuse the existing pool instead of respawning. To change the number of
## slaves, call `npRmpi.stop(force=TRUE)` then restart.
npRmpi.start(nslaves=1)

data("cps71")
mpi.bcast.Robj2slave(cps71)

attach(cps71)

plot(age, logwage, xlab="Age", ylab="log(wage)")

detach(cps71)

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
\references{
  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge	University Press.
}
\keyword{datasets}
