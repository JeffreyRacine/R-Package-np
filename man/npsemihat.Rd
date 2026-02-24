\name{npsemihat}
\alias{npindexhat}
\alias{npplreghat}
\alias{npscoefhat}
\title{Experimental Hat Operators for Semiparametric Estimators}
\description{
Constructs hat operators for semiparametric estimators so that fitted values
or bootstrap draws can be computed by matrix application in one step.
These interfaces are currently experimental.
}
\usage{
npindexhat(bws,
  txdat = stop("training data 'txdat' missing"),
  exdat = txdat,
  y = NULL,
  output = c("matrix", "apply"),
  s = 0L,
  fd.step = NULL,
  ...)

npplreghat(bws,
  txdat = stop("training data 'txdat' missing"),
  tzdat = stop("training data 'tzdat' missing"),
  exdat = txdat,
  ezdat = tzdat,
  y = NULL,
  output = c("apply", "matrix"),
  ...)

npscoefhat(bws,
  txdat = stop("training data 'txdat' missing"),
  tzdat = NULL,
  exdat = txdat,
  ezdat = tzdat,
  y = NULL,
  output = c("matrix", "apply"),
  ridge = 1.0e-12,
  iterate = FALSE,
  leave.one.out = FALSE,
  ...)
}
\arguments{
  \item{bws}{A fitted bandwidth object. \code{npindexhat()} requires class \code{sibandwidth}, \code{npplreghat()} requires class \code{plbandwidth}, and \code{npscoefhat()} requires class \code{scbandwidth}.}
  \item{txdat}{Training \code{x} data.}
  \item{tzdat}{Training \code{z} data for partially linear and smooth-coefficient models.}
  \item{exdat}{Evaluation \code{x} data.}
  \item{ezdat}{Evaluation \code{z} data.}
  \item{y}{Optional response vector or matrix for apply mode.}
  \item{output}{Either \code{"matrix"} or \code{"apply"}.}
  \item{s}{For \code{npindexhat}, \code{s=0} returns fit operator and \code{s=1} returns index-derivative operator.}
  \item{fd.step}{Finite-difference step used for \code{npindexhat(..., s=1)}; must be a positive finite scalar.}
  \item{ridge}{Base ridge term for local linear solves in \code{npscoefhat}; must be a positive finite scalar.}
  \item{iterate}{Logical; \code{npscoefhat()} currently supports \code{iterate = FALSE} only.}
  \item{leave.one.out}{Logical; leave-one-out kernel weights for \code{npscoefhat}.}
  \item{...}{Reserved for future extensions.}
}
\value{
If \code{output = "matrix"}, returns a hat matrix \eqn{H}. If
\code{output = "apply"}, returns \eqn{H y} (or \eqn{H Y} for matrix
right-hand-side input).
}
\details{
These operators are intended for fixed-\eqn{X} workflows such as one-shot wild
bootstrap calculations where many response draws are projected through the same
operator. The implementation is intentionally conservative: class and scalar
argument contracts are validated explicitly, and unsupported iterative
\code{npscoefhat()} paths fail fast.
}
\examples{
\dontrun{
set.seed(42)
n <- 100
x <- runif(n)
y <- sin(2*pi*x) + rnorm(n, sd = 0.1)

bw <- npindexbw(xdat = data.frame(x, x2 = x^2), ydat = y, nmulti = 1)
H <- npindexhat(bws = bw, txdat = data.frame(x, x2 = x^2), output = "matrix")
fit.fast <- H %*% y
}
}
\keyword{nonparametric}
\keyword{semiparametric}
\keyword{regression}
