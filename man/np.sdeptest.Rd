% $Id: np.cmstest.Rd,v 1.58 2006/11/03 21:17:20 tristen Exp $
\name{npsdeptest}
\alias{npsdeptest}
\title{ Kernel Consistent Serial Dependence Test for Univariate Nonlinear Processes }
\description{
  \code{npsdeptest} implements the consistent metric entropy test of
  nonlinear serial dependence as described in Granger, Maasoumi and
  Racine (2004).
}
\usage{

npsdeptest(data = NULL,
           lag.num = 1,
           method = c("integration","summation"),
           bootstrap = TRUE,
           boot.num = 399,
           random.seed = 42)

}

\arguments{
  \item{data}{
    a vector containing the variable that can be of type
    \code{\link{numeric}} or \code{\link{ts}}.
  }
  
  \item{lag.num}{
    an integer value specifying the maximum number of lags to
    use. Defaults to \code{1}.
  }

  \item{method}{
    a character string used to specify whether to compute the integral
    version or the summation version of the statistic. Can be set as
    \code{integration} or \code{summation} (see below for
    details). Defaults to \code{integration}.
  }
  
  \item{bootstrap}{
    a logical value which specifies whether to conduct
    the bootstrap test or not. If set to \code{FALSE}, only the
    statistic will be computed. Defaults to \code{TRUE}.
  }
  
  \item{boot.num}{
    an integer value specifying the number of bootstrap
    replications to use. Defaults to \code{399}.
  }
  
  \item{random.seed}{
    an integer used to seed R's random number generator. This is to
    ensure replicability. Defaults to 42.
  }

}
            
\value{
  \code{npsdeptest} returns an object of type \code{deptest} with the
  following components
  
  \item{Srho}{ the statistic vector \code{Srho} }
  \item{Srho.cumulant}{ the cumulant statistic vector \code{Srho.cumulant} }            
  \item{Srho.bootstrap.mat}{ contains the bootstrap replications of
    \code{Srho} }
  \item{Srho.cumulant.bootstrap.mat}{ contains the bootstrap
    replications of \code{Srho.cumulant} }
  \item{P}{ the P-value vector of the Srho statistic vector }
  \item{P.cumulant}{ the P-value vector of the cumulant Srho statistic vector }  
  \item{bootstrap}{ a logical value indicating whether bootstrapping was
  performed }
  \item{boot.num}{ number of bootstrap replications }
  \item{lag.num}{ the number of lags }
  \item{bw.y}{ the numeric vector of bandwidths for \code{data}
  marginal density at lag \code{num.lag}}
  \item{bw.y.lag}{ the numeric vector of bandwidths for lagged
    \code{data} marginal density  at lag \code{num.lag}}  
  \item{bw.joint}{ the numeric matrix of bandwidths for \code{data}
  and lagged \code{data} joint density  at lag \code{num.lag}}      

\code{\link{summary}} supports object of type \code{deptest}.
}
\references{
  
  Granger, C.W. and E. Maasoumi and J.S. Racine (2004), \dQuote{A
  dependence metric for possibly nonlinear processes}, Journal of Time
  Series Analysis, 25, 649-669.

}
\author{
    Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\details{
  
  \code{npsdeptest} computes the nonparametric metric entropy
   (normalized Hellinger of Granger, Maasoumi and Racine (2004)) for
   testing for nonlinear serial dependence, \eqn{D[f(y_t, \hat y_{t-k}),
   f(y_t)\times f(\hat y_{t-k})]}{D[f(y[t], y[t-k]),
   f(y[t])*f(y[t-k])]}.  Default bandwidths are of the Kullback-Leibler
   variety obtained via likelihood cross-validation.

  The test may be applied to a raw data series or to residuals of user
  estimated models.

  The summation version of this statistic may be numerically unstable
  when \code{data} is sparse (the summation version involves division of
  densities while the integration version involves differences). Warning
  messages are produced should this occur (`integration recommended')
  and should be heeded.
 
}

\section{Usage Issues}{
  
  The \code{integration} version of the statistic uses multdimensional
  numerical methods from the \code{\link{cubature}} package. See
  \code{\link{adaptIntegrate}} for details. The \code{integration}
  version of the statistic will be substantially slower than the
  \code{summation} version, however, it will likely be both more
  accurate and powerful.
  
}

\seealso{
  \code{\link{npdeptest},\link{npdeneqtest},\link{npsymtest},\link{npunitest}}
}

\examples{
set.seed(1234)

## A function to create a time series

ar.series <- function(phi,epsilon) {
  n <- length(epsilon)
  series <- numeric(n)
  series[1] <- epsilon[1]/(1-phi)
  for(i in 2:n) {
    series[i] <- phi*series[i-1] + epsilon[i]
  }
  return(series)
}

n <- 100

## Stationary persistent time-series

yt <- ar.series(0.95,rnorm(n))
npsdeptest(yt,lag.num=2,boot.num=99,method="summation")

Sys.sleep(5)

## Stationary independent time-series

yt <- ar.series(0.0,rnorm(n))
npsdeptest(yt,lag.num=2,boot.num=99,method="summation")

\dontrun{

## Stationary persistent time-series

yt <- ar.series(0.95,rnorm(n))
npsdeptest(yt,lag.num=2,boot.num=99,method="integration")

Sys.sleep(5)

## Stationary independent time-series

yt <- ar.series(0.0,rnorm(n))
npsdeptest(yt,lag.num=2,boot.num=99,method="integration")

} % enddontrun

}
\keyword{ univar }
\keyword{ nonparametric }
