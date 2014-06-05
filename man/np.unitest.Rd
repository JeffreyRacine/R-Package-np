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
    predictability of stock market returns,} Journal of Econometrics,
    107, 2, pp 291-312.

}
\author{
    Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine
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
\dontrun{
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

mpi.spawn.Rslaves(nslaves=1)
mpi.bcast.cmd(np.mpi.initialize(),caller.execute=TRUE)

mpi.bcast.cmd(set.seed(42),
              caller.execute=TRUE)

n <- 250

x <- rnorm(n)
y <- rnorm(n)

mpi.bcast.Robj2slave(x)
mpi.bcast.Robj2slave(y)

mpi.bcast.cmd(output <- npunitest(x,y,
                                  method="summation",
                                  bootstrap=TRUE,
                                  boot.num=99),
              caller.execute=TRUE)

summary(output)

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

mpi.close.Rslaves()

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)
} % enddontrun     
}
\keyword{ univar }
\keyword{ nonparametric }
