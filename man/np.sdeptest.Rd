%
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
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the data series, lag count, and statistic variant.
}

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


\subsection{Bootstrap Controls}{
These arguments control bootstrap execution and reproducibility settings.
}

  \item{boot.num}{
    an integer value specifying the number of bootstrap
    replications to use. Defaults to \code{399}.
  }


  \item{bootstrap}{
    a logical value which specifies whether to conduct
    the bootstrap test or not. If set to \code{FALSE}, only the
    statistic will be computed. Defaults to \code{TRUE}.
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
\section{Book And Method Pointers}{
  \code{npsdeptest} applies the dependence-test idea to a time series by
  comparing the joint law of current and lagged values with the product
  of their marginal laws. The lag and bootstrap choices determine the
  serial-dependence alternatives being probed and the null calibration.

  For book-length background, see Li and Racine (2007), Chapter 13
  \emph{Nonsmoothing Tests}, especially Sections 13.4 and 13.5, and
  Chapter 12 \emph{Model Specification Tests} for related distributional
  tests. For entropy and density context, see Racine (2019), Chapter 2.
}

\references{

  Granger, C.W. and E. Maasoumi and J.S. Racine (2004), \dQuote{A
  dependence metric for possibly nonlinear processes}, Journal of Time
  Series Analysis, 25, 649-669.

}
\author{
    Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, and \code{\link{plot}} for plotting options.


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
  messages are produced should this occur (\sQuote{integration recommended})
  and should be heeded.

}

\section{Usage Issues}{

  The \code{integration} version of the statistic uses multidimensional
  numerical methods from the \pkg{cubature} package. See
  \pkg{adaptIntegrate} for details. The \code{integration}
  version of the statistic will be substantially slower than the
  \code{summation} version, however, it will likely be both more
  accurate and powerful.

}

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}
  \code{\link{npdeptest},\link{npdeneqtest},\link{npsymtest},\link{npunitest}}
}

\examples{
\dontrun{
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

if (interactive()) Sys.sleep(5)

## Stationary independent time-series

yt <- ar.series(0.0,rnorm(n))
npsdeptest(yt,lag.num=2,boot.num=99,method="summation")

## Stationary persistent time-series

yt <- ar.series(0.95,rnorm(n))
npsdeptest(yt,lag.num=2,boot.num=99,method="integration")

if (interactive()) Sys.sleep(5)

## Stationary independent time-series

yt <- ar.series(0.0,rnorm(n))
npsdeptest(yt,lag.num=2,boot.num=99,method="integration")

} % enddontrun

}
\keyword{ univar }
\keyword{ nonparametric }
