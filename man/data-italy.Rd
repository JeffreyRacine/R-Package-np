% $ID$
\name{Italy}
\docType{data}
\alias{Italy}
\title{ Italian GDP Panel }
\description{
  Italian GDP growth panel for 21 regions covering the period 1951-1998
  (millions of Lire, 1990=base). There are 1008 observations in total.
}
\usage{data("Italy")}
\format{ A data frame with 2 columns, and 1008 rows.
  \describe{
    \item{year}{ the first column, of type \code{ordered}}

    \item{gdp}{ the second column, of type \code{numeric}: millions of
    Lire, 1990=base}
  }
}
\source{Giovanni Baiocchi}
\examples{
\donttest{
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave). Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

mpi.spawn.Rslaves(nslaves=1)
mpi.bcast.cmd(np.mpi.initialize(),caller.execute=TRUE)

data("Italy")
mpi.bcast.Robj2slave(Italy)

attach(Italy)

plot(ordered(year), gdp, xlab="Year (ordered factor)",
     ylab="GDP (millions of Lire, 1990=base)")

detach(Italy)

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

mpi.close.Rslaves()

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)
} % enddonttest
}
\references{
  Baiocchi, G. (2006), \dQuote{Economic Applications of Nonparametric
    Methods,} Ph.D. Thesis, University of York.
}
\keyword{datasets}
