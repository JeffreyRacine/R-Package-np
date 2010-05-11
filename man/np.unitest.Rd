\name{npunitest}
\alias{npunitest}
\title{ Kernel Consistent Univariate Density Equality Test with Mixed Data Types}
\description{
  \code{npunitest} implements the consistent metric entropy test of
  Maasoumi and Racine (2002) for two arbitrary, stationary
  univariate nonparametric densities on common support.
}
\usage{
npunitest(data.x = NULL,
          data.y = NULL,
          method = c("integration","summation"),
          bootstrap = TRUE,
          boot.num = 399,
          bw.x = NULL,
          bw.y = NULL,
          random.seed = 42,
          \dots)
}

\arguments{
  \item{data.x, data.y}{
    common support univariate vectors containing the variables.
  }

  \item{method}{
    a character string used to specify whether to compute
    the integral version or the summation version of the statistic. Can
    be set as \code{integration} or \code{summation}. Defaults to
    \code{integration}. See `Details' below for important information
    regarding the use of \code{summation} when \code{data.x} and
    \code{data.y} lack common support and/or are sparse.
  }

  \item{bootstrap}{
    a logical value which specifies whether to conduct the bootstrap
    test or not. If set to \code{FALSE}, only the statistic will be
    computed. Defaults to \code{TRUE}.
  }
  \item{boot.num}{
    an integer value specifying the number of bootstrap
    replications to use. Defaults to \code{399}.
  }

  \item{bw.x, bw.y}{
    numeric (scalar) bandwidths. Defaults to plug-in (see details below).
  }

  \item{random.seed}{
    an integer used to seed R's random number generator. This is to
    ensure replicability. Defaults to 42.
  }

  \item{\dots}{ additional arguments supplied to specify the bandwidth
    type, kernel types, and so on.  This is used since we specify bw as
    a numeric scalar and not a \code{bandwidth} object, and is of
    interest if you do not desire the default behaviours. To change the
    defaults, you may specify any of \code{bwscaling}, \code{bwtype},
    \code{ckertype}, \code{ckerorder}, \code{ukertype},
    \code{okertype}.}

}

\value{
  \code{npunitest} returns an object of type \code{unitest} with the
  following components

  \item{Srho}{ the statistic \code{Srho} }
  \item{Srho.bootstrap}{ contains the bootstrap replications of \code{Srho} }
  \item{P}{ the P-value of the statistic }
  \item{boot.num}{ number of bootstrap replications }
  \item{bw.x, bw.y}{ scalar bandwidths for \code{data.x, data.y} }

  \code{\link{summary}} supports object of type \code{unitest}.
}
\references{

  Granger, C.W. and E. Maasoumi and J.S. Racine (2004), \dQuote{A
  dependence metric for possibly nonlinear processes}, Journal of Time
  Series Analysis, 25, 649-669.

  Maasoumi, E. and J.S. Racine (2002), \dQuote{Entropy and
    Predictability of Stock Market Returns,} Journal of Econometrics,
    107, 2, pp 291-312.

}
\author{
    Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\details{

  \code{npunitest} computes the nonparametric metric entropy (normalized
   Hellinger of Granger, Maasoumi and Racine (2004)) for testing
   equality of two univariate density/probability functions,
   \eqn{D[f(x), f(y)]}{D[f(x), f(y)]}. See Maasoumi and Racine (2002)
   for details. Default bandwidths are of the plug-in variety
   (\code{\link{bw.SJ}} for continuous variables and direct plug-in for
   discrete variables). The bootstrap is conducted via simple resampling
   with replacement from the pooled \code{data.x} and \code{data.y}
   (\code{data.x} only for \code{summation}).

   The summation version of this statistic can be numerically unstable
   when \code{data.x} and \code{data.y} lack common support or when the
   overlap is sparse (the summation version involves division of
   densities while the integration version involves differences, and the
   statistic in such cases can be reported as exactly 0.5 or 0). Warning
   messages are produced when this occurs (`integration recommended')
   and should be heeded.

   Numerical integration can occasionally fail when the \code{data.x}
   and \code{data.y} distributions lack common support and/or lie an
   extremely large distance from one another (the statistic in such
   cases will be reported as exactly 0.5 or 0). However, in these
   extreme cases, simple tests will reveal the obvious differences in
   the distributions and entropy-based tests for equality will be
   clearly unnecessary.

}

\section{Usage Issues}{

  See the example below for proper usage.

}

\seealso{
  \code{\link{npdeneqtest},\link{npdeptest},\link{npsdeptest},\link{npsymtest}}
}

\examples{
set.seed(1234)
n <- 1000

## Compute the statistic only for data drawn from same distribution

x <- rnorm(n)
y <- rnorm(n)

npunitest(x,y,bootstrap=FALSE)

Sys.sleep(5)

## Conduct the test for this data

npunitest(x,y,boot.num=99)

Sys.sleep(5)

## Conduct the test for data drawn from different distributions having
## the same mean and variance

x <- rchisq(n,df=5)
y <- rnorm(n,mean=5,sd=sqrt(10))
mean(x)
mean(y)
sd(x)
sd(y)

npunitest(x,y,boot.num=99)

Sys.sleep(5)

\dontrun{
## Two sample t-test for equality of means
t.test(x,y)
## F test for equality of variances and asymptotic
## critical values
F <- var(x)/var(y)
qf(c(0.025,0.975),df1=n-1,df2=n-1)
} % enddontrun

## Plot the nonparametric density estimates on the same axes

fx <- density(x)
fy <- density(y)
xlim <- c(min(fx$x,fy$x),max(fx$x,fy$x))
ylim <- c(min(fx$y,fy$y),max(fx$y,fy$y))
plot(fx,xlim=xlim,ylim=ylim,xlab="Data",main="f(x), f(y)")
lines(fy$x,fy$y,col="red")

Sys.sleep(5)

## Test for equality of log(wage) distributions

data(wage1)
attach(wage1)
lwage.male <- lwage[female=="Male"]
lwage.female <- lwage[female=="Female"]

npunitest(lwage.male,lwage.female,boot.num=99)

Sys.sleep(5)

## Plot the nonparametric density estimates on the same axes

f.m <- density(lwage.male)
f.f <- density(lwage.female)
xlim <- c(min(f.m$x,f.f$x),max(f.m$x,f.f$x))
ylim <- c(min(f.m$y,f.f$y),max(f.m$y,f.f$y))
plot(f.m,xlim=xlim,ylim=ylim,
     xlab="log(wage)",
     main="Male/Female log(wage) Distributions")
lines(f.f$x,f.f$y,col="red",lty=2)
rug(lwage.male)
legend(-1,1.2,c("Male","Female"),lty=c(1,2),col=c("black","red"))

detach(wage1)

Sys.sleep(5)

## Conduct the test for data drawn from different discrete probability
## distributions

x <- factor(rbinom(n,2,.5))
y <- factor(rbinom(n,2,.1))

npunitest(x,y,boot.num=99)
}
\keyword{ univar }
\keyword{ nonparametric }
