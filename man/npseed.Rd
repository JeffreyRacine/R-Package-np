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
  \code{\link{np}} C backend, resetting the random number generator.
}

\details{
  \code{npseed} provides an interface for setting the random seed (and
  resetting the random number generator) used
  by \code{\link{np}}. The random number generator is used during the
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
  This method currently only supports objects from the \code{\link{np}} library.
}
\examples{
\dontrun{
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

mpi.spawn.Rslaves(nslaves=1)
mpi.bcast.cmd(np.mpi.initialize(),caller.execute=TRUE)

mpi.bcast.cmd(npseed(42),
              caller.execute=TRUE)

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

mpi.close.Rslaves()

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)
} % enddontrun     
}
\keyword{nonparametric}
