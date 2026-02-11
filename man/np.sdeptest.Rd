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
    Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine
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
  \code{\link{npdeptest},\link{npdeneqtest},\link{npsymtest},\link{npunitest}}
}

\examples{
\dontrun{
## Not run in checks: excluded to keep MPI examples stable and check times short.
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

## Start npRmpi for interactive execution. If slaves are already running and
## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
## reuse the existing pool instead of respawning. To change the number of
## slaves, call `npRmpi.stop(force=TRUE)` then restart.
npRmpi.start(nslaves=1)

mpi.bcast.cmd(set.seed(42),
              caller.execute=TRUE)

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

yt <- ar.series(0.95,rnorm(n))

mpi.bcast.Robj2slave(yt)

mpi.bcast.cmd(output <- npsdeptest(yt,
                                   lag.num=2,
                                   boot.num=29,
                                   method="summation"),
              caller.execute=TRUE)

summary(output)

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

## Note: on some systems (notably macOS+MPICH), repeatedly spawning and
## tearing down slaves in the same R session can lead to hangs/crashes.
## npRmpi may therefore keep slave daemons alive by default and
## `npRmpi.stop()` performs a "soft close". Use `force=TRUE` to
## actually shut down the slaves.
##
## You can disable reuse via `options(npRmpi.reuse.slaves=FALSE)` or by
## setting the environment variable `NP_RMPI_NO_REUSE_SLAVES=1` before
## loading the package.

npRmpi.stop()               ## soft close (may keep slaves alive)
## npRmpi.stop(force=TRUE)  ## hard close

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)
} % enddontrun     
}
\keyword{ univar }
\keyword{ nonparametric }
