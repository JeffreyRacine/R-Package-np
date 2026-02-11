\name{npuniden.sc}
\alias{npuniden.sc}

\title{
  Kernel Shape Constrained Bounded Univariate Density Estimation
}
\description{

  \code{npuniden.sc} computes shape constrained kernel univariate
  unconditional density estimates given a vector of continuously
  distributed training data and a bandwidth. Lower and upper bounds
  [\code{a},\code{b}] can be supplied (default is [0,1]) and if \code{a}
  is set to \code{-Inf} there is only one bound on the right, while if
  \code{b} is set to \code{Inf} there is only one bound on the left.

}
\usage{
npuniden.sc(X = NULL,
            Y = NULL,
            h = NULL,
            a = 0,
            b = 1,
            lb = NULL,
            ub = NULL,
            extend.range = 0,
            num.grid = 0,
            function.distance = TRUE,
            integral.equal = FALSE,
            constraint = c("density",
                           "mono.incr",
                           "mono.decr",
                           "concave",
                           "convex",
                           "log-concave",
                           "log-convex"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a required numeric vector of training data lying in \eqn{[a,b]}
}
  \item{Y}{
an optional numeric vector of evaluation data lying in \eqn{[a,b]}
}
  \item{h}{
a bandwidth (\eqn{>0})
}
  \item{a}{
an optional lower bound on the support of \code{X} or \code{Y} (defaults to 0)
}
  \item{b}{
an optional upper bound on the support of \code{X} or \code{Y} (defaults to 1)
}
\item{lb}{
  a scalar lower bound (\eqn{\ge 0}) to be used in conjunction with
  \code{constraint="density"}
}
\item{ub}{
  a scalar upper bound (\eqn{\ge 0} and \eqn{\ge} \code{lb}) to be used
  in conjunction with \code{constraint="density"}
}
\item{extend.range}{
  number specifying the fraction by which the range of the training data
  should be extended for the additional grid points (passed to the
  function \code{extendrange})
}
\item{num.grid}{
  number of additional grid points (in addition to \code{X} and
  \code{Y}) placed on an equi-spaced grid spanning
  \code{extendrange(c(X,Y),f=extend.range)} (if \code{num.grid=0} no
  additional grid points will be used regardless of the value of
  \code{extend.range})
}
\item{function.distance}{
  a logical value that, if \code{TRUE}, minimizes the squared deviation between
  the constrained and unconstrained estimates, otherwise, minimizes the
  squared deviation between the constrained and unconstrained weights
}
\item{integral.equal}{
  a logical value, that, if \code{TRUE}, adjusts the constrained estimate
  to have the same probability mass over the range \code{X}, \code{Y},
  and the additional grid points
}
\item{constraint}{
  a character string indicating whether the estimate is to be
  constrained to be monotonically increasing
  (\code{constraint="mono.incr"}), decreasing
  (\code{constraint="mono.decr"}), convex (\code{constraint="convex"}),
  concave (\code{constraint="concave"}), log-convex
  (\code{constraint="log-convex"}), or log-concave
  (\code{constraint="log-concave"})
}
}

\details{
  Typical usages are (see below for a complete list of options and also
  the examples at the end of this help file)
  \preformatted{
    model <- npuniden.sc(X,a=-2,b=3)
  }

  \code{npuniden.sc} implements a methods for estimating a univariate
  density function defined over a continuous random variable in the
  presence of bounds subject to a variety of shape constraints.  The
  bounded estimates use the truncated Gaussian kernel function.

  Note that for the log-constrained estimates, the derivative estimate
  returned is that for the log-constrained estimate not the non-log
  value of the estimate returned by the function. See Example 5 below
  hat manually plots the log-density and returned derivative (no
  transformation is needed when plotting the density estimate itself).

  If the quadratic program solver fails to find a solution, the
  unconstrained estimate is returned with an immediate warning. Possible
  causes to be investigated are undersmoothing, sparsity, and the
  presence of non-sample grid points. To investigate the possibility of
  undersmoothing try using a larger bandwidth, to investigate sparsity
  try decreasing \code{extend.range}, and to investigate non-sample grid
  points try setting \code{num.grid} to \code{0}.

  Mean square error performance seems to improve generally when using
  additional grid points in the empirical support of \code{X} and
  \code{Y} (i.e., in the observed range of the data sample) but appears
  to deteriorate when imposing constraints beyond the empirical support
  (i.e., when \code{extend.range} is positive). Increasing the number of
  additional points beyond a hundred or so appears to have a limited
  impact.

  The option \code{function.distance=TRUE} appears to perform better for
  imposing convexity, concavity, log-convexity and log-concavity, while
  \code{function.distance=FALSE} appears to perform better for imposing
  monotonicity, whether increasing or decreasing (based on simulations
  for the Beta(s1,s2) distribution with sample size \eqn{n=100}).

}
\value{
A list with the following elements:
\item{f }{unconstrained density estimate}
\item{f.sc }{shape constrained density estimate}
\item{se.f }{asymptotic standard error of the unconstrained density estimate}
\item{se.f.sc }{asymptotic standard error of the shape constrained density estimate}
\item{f.deriv }{unconstrained derivative estimate (of order 1 or 2 or log thereof)}
\item{f.sc.deriv }{shape constrained derivative estimate (of order 1 or 2 or log thereof)}
\item{F }{unconstrained distribution estimate}
\item{F.sc }{shape constrained distribution estimate}
\item{integral.f}{the integral of the unconstrained estimate over \code{X}, \code{Y}, and the additional grid points}
\item{integral.f.sc}{the integral of the constrained estimate over \code{X}, \code{Y}, and the additional grid points}
\item{solve.QP }{logical, if \code{TRUE} \code{solve.QP} succeeded, otherwise failed}
\item{attempts }{number of attempts when \code{solve.QP} fails (max = 9)}
}
\references{
  Du, P. and C. Parmeter and J. Racine (2024), \dQuote{Shape Constrained Kernel PDF and PMF Estimation}, Statistica Sinica, 34 (1), 257-289, \doi{10.5705/ss.202021.0112}
}
\author{
Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\seealso{
  The \pkg{logcondens}, \pkg{LogConDEAD}, and \pkg{scdensity} packages,
  and the function \code{\link{npuniden.boundary}}.
}
\examples{
\dontrun{
## Not run in checks: excluded to keep MPI examples stable and check times short.
n <- 100
set.seed(42)

## Example 1: N(0,1), constrain the density to lie within lb=.1 and ub=.2

X <- sort(rnorm(n))
h <- npuniden.boundary(X,a=-Inf,b=Inf)$h
foo <- npuniden.sc(X,h=h,constraint="density",a=-Inf,b=Inf,lb=.1,ub=.2)
ylim <- range(c(foo$f.sc,foo$f))
plot(X,foo$f.sc,type="l",ylim=ylim,xlab="X",ylab="Density")
lines(X,foo$f,col=2,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

## Example 2: Beta(5,1), DGP is monotone increasing, impose valid
## restriction

X <- sort(rbeta(n,5,1))
h <- npuniden.boundary(X)$h

foo <- npuniden.sc(X=X,h=h,constraint=c("mono.incr"))

oldpar <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
ylim <- range(c(foo$f.sc,foo$f))
plot(X,foo$f.sc,type="l",ylim=ylim,xlab="X",ylab="Density")
lines(X,foo$f,col=2,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

ylim <- range(c(foo$f.sc.deriv,foo$f.deriv))
plot(X,foo$f.sc.deriv,type="l",ylim=ylim,xlab="X",ylab="First Derivative")
lines(X,foo$f.deriv,col=2,lty=2)
abline(h=0,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

## Example 3: Beta(1,5), DGP is monotone decreasing, impose valid
## restriction

X <- sort(rbeta(n,1,5))
h <- npuniden.boundary(X)$h

foo <- npuniden.sc(X=X,h=h,constraint=c("mono.decr"))

par(mfrow=c(1,2))
ylim <- range(c(foo$f.sc,foo$f))
plot(X,foo$f.sc,type="l",ylim=ylim,xlab="X",ylab="Density")
lines(X,foo$f,col=2,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

ylim <- range(c(foo$f.sc.deriv,foo$f.deriv))
plot(X,foo$f.sc.deriv,type="l",ylim=ylim,xlab="X",ylab="First Derivative")
lines(X,foo$f.deriv,col=2,lty=2)
abline(h=0,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")
 
## Example 4: N(0,1), DGP is log-concave, impose invalid concavity
## restriction

X <- sort(rnorm(n))
h <- npuniden.boundary(X,a=-Inf,b=Inf)$h

foo <- npuniden.sc(X=X,h=h,a=-Inf,b=Inf,constraint=c("concave"))

par(mfrow=c(1,2))
ylim <- range(c(foo$f.sc,foo$f))
plot(X,foo$f.sc,type="l",ylim=ylim,xlab="X",ylab="Density")
lines(X,foo$f,col=2,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")
ylim <- range(c(foo$f.sc.deriv,foo$f.deriv))

plot(X,foo$f.sc.deriv,type="l",ylim=ylim,xlab="X",ylab="Second Derivative")
lines(X,foo$f.deriv,col=2,lty=2)
abline(h=0,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

## Example 45: Beta(3/4,3/4), DGP is convex, impose valid restriction

X <- sort(rbeta(n,3/4,3/4))
h <- npuniden.boundary(X)$h

foo <- npuniden.sc(X=X,h=h,constraint=c("convex"))

par(mfrow=c(1,2))
ylim <- range(c(foo$f.sc,foo$f))
plot(X,foo$f.sc,type="l",ylim=ylim,xlab="X",ylab="Density")
lines(X,foo$f,col=2,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

ylim <- range(c(foo$f.sc.deriv,foo$f.deriv))
plot(X,foo$f.sc.deriv,type="l",ylim=ylim,xlab="X",ylab="Second Derivative")
lines(X,foo$f.deriv,col=2,lty=2)
abline(h=0,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

## Example 6: N(0,1), DGP is log-concave, impose log-concavity
## restriction

X <- sort(rnorm(n))
h <- npuniden.boundary(X,a=-Inf,b=Inf)$h

foo <- npuniden.sc(X=X,h=h,a=-Inf,b=Inf,constraint=c("log-concave"))

par(mfrow=c(1,2))

ylim <- range(c(log(foo$f.sc),log(foo$f)))
plot(X,log(foo$f.sc),type="l",ylim=ylim,xlab="X",ylab="Log-Density")
lines(X,log(foo$f),col=2,lty=2)
rug(X)
legend("topleft",c("Constrained-log","Unconstrained-log"),lty=1:2,col=1:2,bty="n")

ylim <- range(c(foo$f.sc.deriv,foo$f.deriv))
plot(X,foo$f.sc.deriv,type="l",ylim=ylim,xlab="X",ylab="Second Derivative of Log-Density")
lines(X,foo$f.deriv,col=2,lty=2)
abline(h=0,lty=2)
rug(X)
legend("topleft",c("Constrained-log","Unconstrained-log"),lty=1:2,col=1:2,bty="n")
par(oldpar)
} % enddontrun
}
\keyword{ nonparametric }
\keyword{ smooth }
