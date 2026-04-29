\name{npcdisthat}
\alias{npcdisthat}
\title{Conditional Distribution Hat Operator}
\description{
Constructs the conditional distribution hat operator associated with
\code{\link{npcdist}} bandwidth objects. The returned operator maps a
right-hand side \code{y} to \eqn{H y}; with \code{y = 1} this reproduces the
fitted conditional distribution function.
}
\usage{
npcdisthat(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           exdat,
           eydat,
           y = NULL,
           output = c("matrix", "apply"))
}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the fitted bandwidth object, training data, and evaluation data.
}

  \item{bws}{A fitted conditional distribution bandwidth object of class \code{"condbandwidth"}.}


  \item{exdat}{Optional evaluation conditioning data. If omitted, the operator is built on the training conditioning data.}


  \item{eydat}{Optional evaluation response data. If omitted, the operator is built on the training response data.}


  \item{txdat}{Training conditioning data used to construct the operator.}


  \item{tydat}{Training response data used to construct the operator.}


\subsection{Operator Output}{
These arguments control whether the operator is returned as a matrix or applied directly.
}

  \item{output}{Either \code{"matrix"} for the hat matrix or \code{"apply"} for direct application to \code{y}.}

  \item{y}{Optional right-hand side vector or matrix with one row per training observation.}


}
\details{
For \code{output = "matrix"}, the return value is a matrix with class
\code{c("npcdisthat", "matrix")} and attributes storing the bandwidth object,
training data, evaluation data, and call metadata.

For \code{output = "apply"}, the function returns \code{H y} directly. Matrix
right-hand sides are applied column-wise.

This helper is intended for object-fed repeated evaluation once a bandwidth
object has already been constructed. It does not perform bandwidth selection.
}
\value{
Either a hat matrix of class \code{"npcdisthat"} or the applied result
\code{H y}, depending on \code{output}.
}
\examples{
\dontrun{
npRmpi.init(nslaves = 1)
data(cps71)
tx <- data.frame(age = cps71$age)
ty <- data.frame(logwage = cps71$logwage)

bw <- npcdistbw(xdat = tx, ydat = ty, bwtype = "fixed",
                bandwidth.compute = FALSE, bws = c(1.0, 1.0))

H <- npcdisthat(bws = bw, txdat = tx, tydat = ty)
dist.hat <- npcdisthat(bws = bw, txdat = tx, tydat = ty,
                       y = rep(1, nrow(tx)),
                       output = "apply")
dist.core <- fitted(npcdist(bws = bw, txdat = tx, tydat = ty))

head(cbind(dist.core, dist.hat), n = 2L)

npRmpi.quit()
}
}
\keyword{nonparametric}
\keyword{distribution}
