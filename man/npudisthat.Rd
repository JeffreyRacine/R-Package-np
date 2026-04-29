\name{npudisthat}
\alias{npudisthat}
\title{Unconditional Distribution Hat Operator}
\description{
Constructs the unconditional distribution hat operator associated with
\code{\link{npudist}} bandwidth objects. The returned operator maps a
right-hand side \code{y} to \eqn{H y}; with \code{y = 1} this reproduces the
fitted unconditional distribution function.
}
\usage{
npudisthat(bws,
           tdat = stop("training data 'tdat' missing"),
           edat,
           y = NULL,
           output = c("matrix", "apply"))
}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the fitted bandwidth object, training data, and evaluation data.
}

  \item{bws}{A fitted unconditional distribution bandwidth object of class \code{"dbandwidth"}.}


  \item{edat}{Optional evaluation data. If omitted, the operator is built on the training data.}


  \item{tdat}{Training data used to construct the operator.}


\subsection{Operator Output}{
These arguments control whether the operator is returned as a matrix or applied directly.
}

  \item{output}{Either \code{"matrix"} for the hat matrix or \code{"apply"} for direct application to \code{y}.}

  \item{y}{Optional right-hand side vector or matrix with one row per training observation.}


}
\details{
For \code{output = "matrix"}, the return value is a matrix with class
\code{c("npudisthat", "matrix")} and attributes storing the bandwidth object,
training data, evaluation data, and call metadata.

For \code{output = "apply"}, the function returns \code{H y} directly. Matrix
right-hand sides are applied column-wise.

This helper is intended for object-fed repeated evaluation once a bandwidth
object has already been constructed. It does not perform bandwidth selection.
}
\value{
Either a hat matrix of class \code{"npudisthat"} or the applied result
\code{H y}, depending on \code{output}.
}
\examples{
\dontrun{
npRmpi.init(nslaves = 1)
data(cps71)
tx <- data.frame(age = cps71$age)

bw <- npudistbw(dat = tx, bwtype = "fixed",
                bandwidth.compute = FALSE, bws = 1.0)

H <- npudisthat(bws = bw, tdat = tx)
dist.hat <- npudisthat(bws = bw, tdat = tx,
                       y = rep(1, nrow(tx)),
                       output = "apply")
dist.core <- fitted(npudist(bws = bw, tdat = tx))

head(cbind(dist.core, dist.hat), n = 2L)

npRmpi.quit()
}
}
\keyword{nonparametric}
\keyword{distribution}
