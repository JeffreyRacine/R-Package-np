\name{uocquantile}
\title{Compute Quantiles}
\usage{
uocquantile(x, prob)
}
\alias{uocquantile}
\arguments{
  \item{x}{an ordered, unordered or continuous variable.}
  \item{prob}{quantile to compute.}
}
\value{
  A quantile computed from \code{x}.
}
\description{
  \code{uocquantile} is a function which computes quantiles of an
  unordered, ordered, or continuous variable  \code{x}.
}
\details{
  \code{uocquantile} is a function which computes quantiles of
  unordered, ordered, or a continuous variable  \code{x}. If \code{x}
  is unordered, the mode is returned. If \code{x} is ordered, the level
  for which the cumulative distribution is >= prob is returned. If
  \code{x} is continuous, \code{\link{quantile}} is invoked and the
  result returned. 
}
\seealso{
  \code{\link{quantile}}
}
\author{
  Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\examples{
x <- rbinom(n = 100, size = 10, prob = 0.5)
uocquantile(x, 0.5)
}

\keyword{univar}
\keyword{nonparametric}

