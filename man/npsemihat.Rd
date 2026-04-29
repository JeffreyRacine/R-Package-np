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
           ridge = 0,
           iterate = FALSE,
           leave.one.out = FALSE,
           ...)
}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the fitted bandwidth object, training data, and evaluation data.
}

  \item{bws}{A fitted bandwidth object. \code{npindexhat()} requires class \code{sibandwidth}, \code{npplreghat()} requires class \code{plbandwidth}, and \code{npscoefhat()} requires class \code{scbandwidth}.}


  \item{exdat}{Evaluation \code{x} data.}


  \item{txdat}{Training \code{x} data.}


\subsection{Operator Output And Derivatives}{
These arguments control operator output, derivative selection, finite-difference compatibility, and apply-mode right-hand sides.
}

  \item{fd.step}{Compatibility argument for \code{npindexhat(..., s=1)}. If supplied it must be a positive finite scalar, but the current \code{s=1} route uses the canonical exact derivative operator rather than finite-differencing the fit operator.}


  \item{output}{Either \code{"matrix"} or \code{"apply"}.}


  \item{s}{For \code{npindexhat}, \code{s=0} returns fit operator and \code{s=1} returns index-derivative operator.}


  \item{y}{Optional response vector or matrix for apply mode.}


\subsection{Partially Linear And Smooth-Coefficient Data}{
These arguments supply the additional \code{z} data used by partially linear and smooth-coefficient models.
}

  \item{ezdat}{Evaluation \code{z} data.}


  \item{tzdat}{Training \code{z} data for partially linear and smooth-coefficient models.}


\subsection{Smooth-Coefficient Controls}{
These arguments control smooth-coefficient iteration, leave-one-out behavior, and ridge stabilization.
}

  \item{iterate}{Logical; \code{npscoefhat()} currently supports \code{iterate = FALSE} only.}


  \item{leave.one.out}{Logical; leave-one-out kernel weights for \code{npscoefhat}. This currently requires evaluation \code{z} data to match training \code{z} data.}


  \item{ridge}{Base ridge term for local linear solves in \code{npscoefhat}; must be a non-negative finite scalar. The ridge sequence starts at \code{0} (no regularization) and then increments by \code{1/n.train} as needed for stable solves.}


\subsection{Additional Arguments}{
Reserved for future extensions.
}

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
\code{npscoefhat()} paths fail fast. \code{npscoefhat()} inherits
\code{regtype}/LP-basis controls from the supplied \code{scbandwidth}
object. For non-fixed \code{npscoef} bootstrap plotting, these operators can
support a frozen approximation, but they do not remove the need to recompute
the local smooth-coefficient vector itself for each resample: the local
weighted systems depend on the resample weights/counts at each evaluation
point, so unlike \code{npplreg} there is no single global coefficient vector
to update once per draw.

Method-specific argument map:
\code{npindexhat()} uses \code{s}; \code{fd.step} is accepted for compatibility but the current \code{s=1} route uses the canonical exact derivative operator;
\code{npplreghat()} and \code{npscoefhat()} use \code{tzdat}/\code{ezdat};
\code{npscoefhat()} additionally uses \code{ridge}, \code{iterate}, and
\code{leave.one.out}.
}
\examples{
\dontrun{
npRmpi.init(nslaves = 1)
set.seed(42)
n <- 100
x <- runif(n)
z <- runif(n)
y <- sin(2*pi*x) + 0.5 * z + rnorm(n, sd = 0.1)

tx <- data.frame(x = x)
tz <- data.frame(z = z)

ibw <- npindexbw(xdat = data.frame(x, x2 = x^2), ydat = y,
                 bws = c(0.5, 1.0, 1.0), bandwidth.compute = FALSE)
iH <- npindexhat(bws = ibw, txdat = data.frame(x, x2 = x^2), output = "matrix")
iH.fitted <- iH %*% y
ifit <- npindex(bws = ibw, txdat = data.frame(x, x2 = x^2), tydat = y)
head(cbind(fitted(ifit), iH.fitted), n = 2L)

pbw <- npplregbw(xdat = tx, zdat = tz, ydat = y,
                 bws = matrix(c(0.2, 0.2), nrow = 2L, ncol = 1L),
                 bandwidth.compute = FALSE)
pH <- npplreghat(bws = pbw, txdat = tx, tzdat = tz, output = "matrix")
pH.fitted <- pH %*% y
pfit <- npplreg(bws = pbw, txdat = tx, tydat = y, tzdat = tz)
head(cbind(fitted(pfit), pH.fitted), n = 2L)

sbw <- npscoefbw(xdat = tx, zdat = tz, ydat = y,
                 bws = 0.2, bandwidth.compute = FALSE)
sH <- npscoefhat(bws = sbw, txdat = tx, tzdat = tz,
                 output = "matrix", iterate = FALSE)
sH.fitted <- sH %*% y
sfit <- npscoef(bws = sbw, txdat = tx, tydat = y, tzdat = tz,
                iterate = FALSE)
head(cbind(fitted(sfit), sH.fitted), n = 2L)

npRmpi.quit()
}
}
\keyword{nonparametric}
\keyword{semiparametric}
\keyword{regression}
