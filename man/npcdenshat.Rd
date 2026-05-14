\name{npcdenshat}
\alias{npcdenshat}
\title{Conditional Density Hat Operator}
\description{
Constructs the conditional density hat operator associated with
\code{\link{npcdens}} bandwidth objects. The returned operator maps a
right-hand side \code{y} to \eqn{H y}; with \code{y = 1} this reproduces the
fitted conditional density.
}
\usage{
npcdenshat(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           exdat,
           eydat,
           y = NULL,
           output = c("matrix", "apply"),
           deriv = NULL,
           s = NULL)
}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the fitted bandwidth object, training data, and evaluation data.
}

  \item{bws}{A fitted conditional density bandwidth object of class \code{"conbandwidth"}.}


  \item{exdat}{Optional evaluation conditioning data. If omitted, the operator is built on the training conditioning data.}


  \item{eydat}{Optional evaluation response data. If omitted, the operator is built on the training response data.}


  \item{txdat}{Training conditioning data used to construct the operator.}


  \item{tydat}{Training response data used to construct the operator.}


\subsection{Operator Output}{
These arguments control whether the operator is returned as a matrix or applied directly.
}

  \item{deriv}{Convenience alias for \code{s}.}

  \item{output}{Either \code{"matrix"} for the hat matrix or \code{"apply"} for direct application to \code{y}.}

  \item{s}{Optional first-derivative multi-index over continuous conditioning variables. By default the conditional-density operator itself is returned. When \code{s} requests one first derivative, the returned operator maps \code{y} to the corresponding \code{npcdens(..., gradients = TRUE)} conditioning-variable gradient component.}

  \item{y}{Optional right-hand side vector or matrix with one row per training observation.}


}
\details{
For \code{output = "matrix"}, the return value is a matrix with class
\code{c("npcdenshat", "matrix")} and attributes storing the bandwidth object,
training data, evaluation data, and call metadata.

For \code{output = "apply"}, the function returns \code{H y} directly. Matrix
right-hand sides are applied column-wise.

The optional \code{s} argument follows the derivative multi-index convention
used by \code{\link{npreghat}} for continuous conditioning variables. The
default \code{s = NULL} returns the ordinary conditional-density hat operator.
Supplying one first derivative, for example \code{s = 1L} in a univariate
continuous-\eqn{x} problem or \code{s = c(x2 = 1L)} in a multivariate problem,
returns the hat operator for that component of the conditional-density
gradient. This derivative route is intended to match
\code{npcdens(..., gradients = TRUE)$congrad}.

This helper is intended for object-fed repeated evaluation once a bandwidth
object has already been constructed. It does not perform bandwidth selection.
}
\value{
Either a hat matrix of class \code{"npcdenshat"} or the applied result
\code{H y}, depending on \code{output}.
}
\examples{
\dontrun{
data(cps71)
tx <- data.frame(age = cps71$age)
ty <- data.frame(logwage = cps71$logwage)

bw <- npcdensbw(xdat = tx, ydat = ty, bwtype = "fixed",
                bandwidth.compute = FALSE, bws = c(1.0, 1.0))

H <- npcdenshat(bws = bw, txdat = tx, tydat = ty)
dens.hat <- npcdenshat(bws = bw, txdat = tx, tydat = ty,
                       y = rep(1, nrow(tx)),
                       output = "apply")
dens.core <- fitted(npcdens(bws = bw, txdat = tx, tydat = ty))

head(cbind(dens.core, dens.hat), n = 2L)

grad.hat <- npcdenshat(bws = bw, txdat = tx, tydat = ty,
                       y = rep(1, nrow(tx)),
                       output = "apply", s = 1L)
grad.core <- gradients(npcdens(bws = bw, txdat = tx, tydat = ty,
                               gradients = TRUE))[, 1L]
head(cbind(grad.core, grad.hat), n = 2L)
}
}
\keyword{nonparametric}
\keyword{density}
