\name{nplsqreg}
\alias{nplsqreg}
\alias{nplsqreg.formula}
\alias{nplsqreg.default}
\alias{nplsqreg.lsqregressionbandwidth}
\alias{fitted.lsqregression}
\alias{quantile.lsqregression}
\alias{predict.lsqregression}
\alias{se.lsqregression}
\alias{gradients.lsqregression}
\alias{summary.lsqregression}
\alias{print.lsqregression}
\alias{plot.lsqregression}
\title{Location-scale conditional quantile regression}
\description{
Fits a scalar-\code{tau} location-scale conditional quantile regression
estimator. The final conditional quantile surface is obtained by applying
\code{\link{npreg}} to the transformed response
\code{Y + scale_hat * qnorm(delta)} using bandwidths and \code{delta} selected
by check-loss cross-validation.
}
\usage{
nplsqreg(bws, ...)
}
\arguments{
  \item{bws}{a formula, an \code{lsqregressionbandwidth} object, an
    \code{rbandwidth} object, a numeric bandwidth vector, or omitted for
    automatic bandwidth selection.}
  \item{...}{additional arguments passed to \code{nplsqregbw} and
    \code{npreg}.}
}
\details{
This scratch implementation intentionally accepts only a scalar
\code{tau}. Vector \code{tau} is not supported because the location-scale
bandwidths and shift parameter are optimized for a specific quantile index.

The estimator reuses \code{\link{npreg}} for mixed-data fitting, prediction,
standard errors of the transformed-response regression, gradients, and
categorical effects.
}
\value{
An object of class \code{lsqregression}.
}
\seealso{
\code{\link{nplsqregbw}}, \code{\link{npreg}}, \code{\link{npqreg}}
}
\keyword{nonparametric}
\keyword{regression}
