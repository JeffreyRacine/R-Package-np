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
            constraint = c("mono.incr",
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
an bandwidth (>0)
}
  \item{a}{
an optional lower bound (defaults to 0)
}
  \item{b}{
an optional upper bound (defaults to 1)
}
  
  \item{constraint}{
a character string indicating whether the estimate is to be constrained
to be monotonically increasing (\code{constraint="mono.incr"}),
decreasing (\code{constraint="mono.incr"}), convex
(\code{constraint="convex"}), concave (\code{constraint="concave"}), log-convex
(\code{constraint="log-convex"}), or log-concave (\code{constraint="log-concave"})
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

}
\value{
A list with the following elements:
\item{f }{unconstrained density estimate}
\item{f.deriv }{unconstrained derivative estimate (or order 1 or 2)}
\item{f.sc }{shape constrained density estimate}
\item{f.sc.deriv }{shape constrained derivative estimate (or order 1 or 2)}
}
\references{
TBD
}
\author{
Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\seealso{
  The \pkg{logcondens} and \pkg{scdensity} packages, and the function
  \code{\link{npuniden.boundary}}.
}
\examples{
\dontrun{
n <- 250
set.seed(42)

## Example 1: Beta(5,1), DGP is monotone increasing, impose valid
## restriction

X <- sort(rbeta(n,5,1))
Y <- seq(max(c(0,extendrange(X,f=.25)[1])),min(c(1,extendrange(X,f=.25)[2])),,1000)
h <- bw.nrd0(X)

foo <- npuniden.sc(X=X,Y=Y,h=h,constraint=c("mono.incr"))

par(mfrow=c(1,2))
ylim <- range(c(foo$f.sc,foo$f))
plot(Y,foo$f.sc,type="l",ylim=ylim,xlab="X",ylab="Density")
lines(Y,foo$f,col=2,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

ylim <- range(c(foo$f.sc.deriv,foo$f.deriv))
plot(Y,foo$f.sc.deriv,type="l",ylim=ylim,xlab="X",ylab="First Derivative")
lines(Y,foo$f.deriv,col=2,lty=2)
abline(h=0,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

## Example 2: Beta(1,5), DGP is monotone decreasing, impose valid
## restriction

X <- sort(rbeta(n,1,5))
Y <- seq(max(c(0,extendrange(X,f=.25)[1])),min(c(1,extendrange(X,f=.25)[2])),,1000)
h <- bw.nrd0(X)

foo <- npuniden.sc(X=X,Y=Y,h=h,constraint=c("mono.decr"))

par(mfrow=c(1,2))
ylim <- range(c(foo$f.sc,foo$f))
plot(Y,foo$f.sc,type="l",ylim=ylim,xlab="X",ylab="Density")
lines(Y,foo$f,col=2,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

ylim <- range(c(foo$f.sc.deriv,foo$f.deriv))
plot(Y,foo$f.sc.deriv,type="l",ylim=ylim,xlab="X",ylab="First Derivative")
lines(Y,foo$f.deriv,col=2,lty=2)
abline(h=0,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")
 
## Example 3: N(0,1), DGP is log-concave, impose invalid concavity
## restriction

X <- sort(rnorm(n))
Y <- seq(min(X),max(X),,1000)
h <- bw.nrd0(X)

foo <- npuniden.sc(X=X,Y=Y,h=h,a=-Inf,b=Inf,constraint=c("concave"))

par(mfrow=c(1,2))
ylim <- range(c(foo$f.sc,foo$f))
plot(Y,foo$f.sc,type="l",ylim=ylim,xlab="X",ylab="Density")
lines(Y,foo$f,col=2,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")
ylim <- range(c(foo$f.sc.deriv,foo$f.deriv))

plot(Y,foo$f.sc.deriv,type="l",ylim=ylim,xlab="X",ylab="Second Derivative")
lines(Y,foo$f.deriv,col=2,lty=2)
abline(h=0,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

## Example 4: Beta(3/4,3/4), DGP is convex, impose valid restriction

X <- sort(rbeta(n,3/4,3/4))
Y <- seq(max(c(0,extendrange(X,f=.25)[1])),min(c(1,extendrange(X,f=.25)[2])),,1000)
h <- bw.nrd0(X)

foo <- npuniden.sc(X=X,Y=Y,h=h,constraint=c("convex"))

par(mfrow=c(1,2))
ylim <- range(c(foo$f.sc,foo$f))
plot(Y,foo$f.sc,type="l",ylim=ylim,xlab="X",ylab="Density")
lines(Y,foo$f,col=2,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

ylim <- range(c(foo$f.sc.deriv,foo$f.deriv))
plot(Y,foo$f.sc.deriv,type="l",ylim=ylim,xlab="X",ylab="Second Derivative")
lines(Y,foo$f.deriv,col=2,lty=2)
abline(h=0,lty=2)
rug(X)
legend("topleft",c("Constrained","Unconstrained"),lty=1:2,col=1:2,bty="n")

## Example 5: N(0,1), DGP is log-concave, impose log-concavity
## restriction

X <- sort(rnorm(n))
Y <- seq(min(X),max(X),,1000)
h <- bw.nrd0(X)

foo <- npuniden.sc(X=X,Y=Y,h=h,a=-Inf,b=Inf,constraint=c("log-concave"))

par(mfrow=c(1,2))

ylim <- range(c(log(foo$f.sc),log(foo$f)))
plot(Y,log(foo$f.sc),type="l",ylim=ylim,xlab="X",ylab="Log-Density")
lines(Y,log(foo$f),col=2,lty=2)
rug(X)
legend("topleft",c("Constrained-log","Unconstrained-log"),lty=1:2,col=1:2,bty="n")

ylim <- range(c(foo$f.sc.deriv,foo$f.deriv))
plot(Y,foo$f.sc.deriv,type="l",ylim=ylim,xlab="X",ylab="Second Derivative of Log-Density")
lines(Y,foo$f.deriv,col=2,lty=2)
abline(h=0,lty=2)
rug(X)
legend("topleft",c("Constrained-log","Unconstrained-log"),lty=1:2,col=1:2,bty="n")
} % enddontrun
}
\keyword{ nonparametric }
\keyword{ smooth }
