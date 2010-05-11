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

\author{ Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey
  S. Racine \email{racinej@mcmaster.ca}
}

\section{Note}{
  This method currently only supports objects from the \code{\link{np}} library.
}
\examples{
npseed(712)
x <- runif(10)
y <- x + rnorm(10, sd = 0.1)
bw <- npregbw(y~x)
}
\keyword{nonparametric}
