\name{nplsqregbw}
\alias{nplsqregbw}
\alias{nplsqregbw.formula}
\alias{nplsqregbw.default}
\alias{nplsqregbw.lsqregressionbandwidth}
\alias{lsqregressionbandwidth}
\alias{print.lsqregressionbandwidth}
\alias{summary.lsqregressionbandwidth}
\title{Bandwidth selection for location-scale conditional quantile regression}
\description{
Selects fixed bandwidths and the scalar location-scale shift parameter
\code{delta} for \code{\link{nplsqreg}} by minimizing leave-one-out check
loss.
}
\usage{
nplsqregbw(bws, ...)
}
\arguments{
  \item{bws}{a formula, an \code{lsqregressionbandwidth} object, an
    \code{rbandwidth} object, a numeric bandwidth vector, or omitted for
    automatic selection.}
  \item{...}{additional arguments, including \code{xdat}, \code{ydat},
    \code{tau}, \code{delta}, \code{scale}, \code{regtype}, \code{degree},
    and kernel/bandwidth arguments forwarded to \code{\link{npregbw}}.}
}
\details{
The default pilot scale is a Fan-Yao-style residual scale: first
\code{\link{npreg}} estimates the conditional mean, then \code{\link{npreg}}
smooths squared residuals, and the square root of that fitted variance is used
as \code{scale_hat}. A positive user-supplied scale vector may be provided via
\code{scale}.

This scratch implementation supports fixed bandwidths and scalar \code{tau}.
}
\value{
An object of class \code{lsqregressionbandwidth}.
}
\seealso{
\code{\link{nplsqreg}}, \code{\link{npregbw}}
}
\keyword{nonparametric}
\keyword{regression}
