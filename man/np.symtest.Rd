\name{npsymtest}
\alias{npsymtest}
\title{ Kernel Consistent Density Asymmetry Test with Mixed Data Types }
\description{
  \code{npsymtest} implements the consistent metric entropy test of
  asymmetry as described in Maasoumi and Racine (2009).
}
\usage{
npsymtest(data = NULL,
          method = c("integration","summation"),
          boot.num = 399,
          bw = NULL,
          boot.method = c("iid", "geom"),
          random.seed = 42,
          \dots)
}

\arguments{
  \item{data}{
    a vector containing the variable.
  }
  
  \item{method}{
    a character string used to specify whether to compute the integral
    version or the summation version of the statistic. Can be set as
    \code{integration} or \code{summation} (see below for
    details). Defaults to \code{integration}.
  }
  
  \item{boot.num}{
    an integer value specifying the number of bootstrap
    replications to use. Defaults to \code{399}.
  }
  
  \item{bw}{
    a numeric (scalar) bandwidth. Defaults to plug-in (see details below).
  }
  
  \item{boot.method}{ a character string used to specify the
    bootstrap method. Can be set as \code{iid} or \code{geom} (see below
    for details). Defaults to \code{iid}.
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
  \code{npsymtest} returns an object of type \code{symtest} with the
  following components
  
  \item{Srho}{ the statistic \code{Srho} }
  \item{Srho.bootstrap}{ contains the bootstrap replications of \code{Srho} }
  \item{P}{ the P-value of the statistic }
  \item{boot.num}{ number of bootstrap replications }
  \item{data.rotate}{ the rotated data series }
  \item{bw}{ the numeric (scalar) bandwidth }
    
  \code{\link{summary}} supports object of type \code{symtest}.
}
\references{

  Granger, C.W. and E. Maasoumi and J.S. Racine (2004), \dQuote{A
  dependence metric for possibly nonlinear processes}, Journal of Time
  Series Analysis, 25, 649-669.

  Maasoumi, E. and J.S. Racine (2009), \dQuote{A robust entropy-based
  test of asymmetry for discrete and continuous processes,} Econometric
  Reviews, 28, 246-261.

  Politis, D.N. and J.P. Romano (1994), \dQuote{The stationary
    bootstrap,} Journal of the American Statistical Association, 89,
    1303-1313.

}
\author{
    Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\details{

  \code{npsymtest} computes the nonparametric metric entropy (normalized
   Hellinger of Granger, Maasoumi and Racine (2004)) for testing
   symmetry using the densities/probabilities of the data and the
   rotated data, \eqn{D[f(y), f(\tilde y)]}{D[f(y), f(y-tilde)]}. See
   Maasoumi and Racine (2009) for details. Default bandwidths are of the
   plug-in variety (\code{\link{bw.SJ}} for continuous variables and
   direct plug-in for discrete variables).
 
  For bootstrapping the null distribution of the statistic, \code{iid}
  conducts simple random resampling, while \code{geom} conducts Politis
  and Romano's (1994) stationary bootstrap using automatic block length
  selection via the \code{\link{b.star}} function in the
  \code{\link{np}} package. See the \code{\link{boot}} package for
  details.

  The summation version of this statistic may be numerically unstable
  when \code{y} is sparse (the summation version involves division of
  densities while the integration version involves differences). Warning
  messages are produced should this occur (`integration recommended')
  and should be heeded.
 
}

\section{Usage Issues}{
  
  When using data of type \code{\link{factor}} it is crucial that the
  variable not be an alphabetic character string (i.e. the factor must
  be integer-valued). The rotation is conducted about the median after
  conversion to type \code{\link{numeric}} which is then converted back
  to type \code{\link{factor}}. Failure to do so will have unpredictable
  results. See the example below for proper usage.
  
}

\seealso{
  \code{\link{npdeneqtest},\link{npdeptest},\link{npsdeptest},\link{npunitest}}
}

\examples{
set.seed(1234)

n <- 100

## Asymmetric discrete probability distribution function

x <- factor(rbinom(n,2,.8))
npsymtest(x,boot.num=99)

Sys.sleep(5)

## Symmetric discrete probability distribution function

x <- factor(rbinom(n,2,.5))
npsymtest(x,boot.num=99)

Sys.sleep(5)

## Asymmetric continuous distribution function

y <- rchisq(n,df=2)
npsymtest(y,boot.num=99)

Sys.sleep(5)

## Symmetric continuous distribution function

y <- rnorm(n)
npsymtest(y,boot.num=99)

\dontrun{

## Time-series bootstrap

ar.series <- function(phi,epsilon) {
  n <- length(epsilon)
  series <- numeric(n)
  series[1] <- epsilon[1]/(1-phi)
  for(i in 2:n) {
    series[i] <- phi*series[i-1] + epsilon[i]
  }
  return(series)
}

## Asymmetric time-series

yt <- ar.series(0.5,rchisq(n,df=3))
npsymtest(yt,boot.num=99,boot.method="geom")

## Symmetric time-series

yt <- ar.series(0.5,rnorm(n))
npsymtest(yt,boot.num=99,boot.method="geom")

} % enddontrun

}
\keyword{ univar }
\keyword{ nonparametric }
