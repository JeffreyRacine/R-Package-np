\name{npreghat}
\alias{npreghat}
\alias{npreghat.formula}
\alias{npreghat.rbandwidth}
\alias{npreghat.npregression}
\alias{predict.npreghat}
\alias{print.npreghat}
\title{Nonparametric Regression Hat Operator}
\description{
Constructs nonparametric regression hat operators for \code{npreg}-compatible
bandwidth objects. The returned operator \eqn{H^{(s)}} maps responses to fitted
values or derivative estimates via \eqn{H^{(s)} y}.
}
\usage{
npreghat(bws, ...)

\method{npreghat}{formula}(bws, 
         data = NULL, 
         newdata = NULL, 
         ...)


\method{npreghat}{rbandwidth}(bws,
         txdat = stop("training data 'txdat' missing"),
         exdat, y = NULL,
         output = c("matrix", "apply"),
         basis = NULL,
         bernstein.basis = NULL,
         degree = NULL,
         deriv = NULL,
         leave.one.out = FALSE,
         ridge = 0,
         s = NULL, 
         ...)

\method{npreghat}{npregression}(bws, txdat, y, ...)

\method{predict}{npreghat}(object, 
        newdata = NULL, 
        y = NULL,
        output = c("matrix", "apply"),
        s = attr(object, "s"),
        leave.one.out = attr(object, "leave.one.out"),
        deriv = NULL, 
        ...)
}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the fitted bandwidth object, formula/data interface, training data, and evaluation data.
}

  \item{bws}{An object of class \code{rbandwidth} or \code{npregression}.}


  \item{data}{A data frame used with the formula interface.}


  \item{exdat}{Optional evaluation predictors.}


  \item{newdata}{Optional evaluation data for formula and predict methods.}


  \item{txdat}{Training predictors.}


\subsection{Local-Polynomial Controls}{
These arguments control local-polynomial basis, degree, derivatives, leave-one-out behavior, and ridge stabilization.
}

  \item{basis}{Local polynomial basis: \code{"glp"}, \code{"additive"}, or \code{"tensor"}.}


  \item{bernstein.basis}{Logical; use Bernstein basis for LP terms.}


  \item{degree}{Optional local polynomial degree vector override (LP path).}


  \item{deriv}{Convenience alias for \code{s}.}


  \item{leave.one.out}{Logical; if \code{TRUE}, compute in-sample leave-one-out hat weights. This cannot be combined with explicit \code{exdat}/\code{newdata}.}


  \item{ridge}{Base diagonal regularization used when local systems are ill-conditioned. The ridge sequence starts at \code{0} (no regularization) and then increments by \code{1/n.train} as needed for stable solves.}


  \item{s}{Derivative multi-index over continuous predictors.}


\subsection{Method Objects}{
This argument identifies a fitted hat-operator object supplied to an S3 method.
}

  \item{object}{An object returned by \code{npreghat}.}


\subsection{Operator Output}{
These arguments control whether the operator is returned as a matrix or applied directly.
}

  \item{output}{Either \code{"matrix"} for the hat matrix or \code{"apply"} for direct application to \code{y}.}


  \item{y}{Optional response vector or matrix for apply mode.}


\subsection{Additional Arguments}{
Further arguments are passed to methods.
}

  \item{...}{Additional arguments passed to methods.}

}
\details{
For \code{output = "matrix"}, the return value is a matrix with class
\code{c("npreghat", "matrix")} so it can be used directly in matrix products,
e.g. \code{H \%*\% y}. Attributes on the matrix store metadata used by
\code{predict.npreghat}.

For \code{output = "apply"}, the function returns \code{H^{(s)} y} directly and
accepts matrix right-hand sides for one-shot bootstrap-style calculations.
}
\value{
Either a hat matrix (class \code{"npreghat"}) or the applied result
\code{H^{(s)} y}, depending on \code{output}.
}
\examples{
\dontrun{
npRmpi.init(nslaves = 1)
data(cps71)
bw <- npregbw(xdat = cps71$age, ydat = cps71$logwage,
              regtype = "ll", bandwidth.compute = FALSE, bws = 1.0)
H <- npreghat(bws = bw, txdat = data.frame(age = cps71$age))
H.fitted <- H %*% cps71$logwage
ghat <- npreg(bws = bw)
head(cbind(fitted(ghat), H.fitted), n = 2L)

npRmpi.quit()
}
}
\keyword{nonparametric}
\keyword{regression}
