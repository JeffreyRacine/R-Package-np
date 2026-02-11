\name{np.pairs}
\alias{np.pairs}
\alias{np.pairs.plot}
\title{Cross-Validated Pairs Plot (Helper Functions)}
\description{
Compute pairwise nonparametric regressions and densities for a set of
variables, then plot a pairs-style display with fitted smoothers.
}
\usage{
np.pairs(y_vars, y_dat, ...)
np.pairs.plot(pair_list)
}
\arguments{
  \item{y_vars}{
    character vector of column names in \code{y_dat}. If \code{y_vars}
    is named, the names are used as plot labels.
  }
  \item{y_dat}{
    data frame containing the variables listed in \code{y_vars}.
  }
  \item{...}{
    additional arguments passed to \code{\link{npudens}} and
    \code{\link{npreg}}.
  }
  \item{pair_list}{
    list returned by \code{np.pairs}.
  }
}
\value{
  \code{np.pairs} returns a list with components
  \code{y_vars}, \code{pair_names}, and \code{pair_kerns}.
  \code{np.pairs.plot} returns \code{NULL} (invisibly).
}
\details{
On the diagonal, \code{npudens} is used to compute kernel density
estimates. Off-diagonal panels use \code{npreg} with residuals to draw
scatterplots and smoothers.
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

data("USArrests")
y_vars <- c("Murder", "UrbanPop")
names(y_vars) <- c("Murder Arrests per 100K", "Pop. Percent Urban")

mpi.bcast.Robj2slave(USArrests)
mpi.bcast.Robj2slave(y_vars)

mpi.bcast.cmd(pair_list <- np.pairs(y_vars = y_vars, y_dat = USArrests,
                                    ckertype = "epanechnikov", 
                                    bwscaling = TRUE),
              caller.execute=TRUE)

np.pairs.plot(pair_list)

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
\seealso{\code{\link{npudens}}, \code{\link{npreg}}}
\keyword{hplot}
