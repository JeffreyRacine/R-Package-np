\name{npquantile}
\alias{npquantile}

\title{ Kernel Univariate Quantile Estimation }
\description{

  \code{npquantile} computes smooth quantiles from a univariate
  unconditional kernel cumulative distribution estimate given data and,
  optionally, a bandwidth specification i.e. a \code{dbandwidth} object
  using the bandwidth selection method of Li, Li and Racine (2017).

}
\usage{
npquantile(x = NULL,
           tau = c(0.01,0.05,0.25,0.50,0.75,0.95,0.99),
           num.eval = 10000,
           bws = NULL,
           f = 1,
           \dots)
}

\arguments{
  \item{x}{
    a univariate vector of type \code{\link{numeric}} containing sample
    realizations (training data) used to estimate the cumulative
    distribution (must be the same training data used to compute the
    bandwidth object \code{bws} passed in).
  }
  \item{tau}{
    an optional vector containing the probabilities for quantile(s) to
    be estimated (must contain numbers in \eqn{[0,1]}). Defaults to
    \code{c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)}.
  }
  \item{num.eval}{
    an optional integer specifying the length of the grid on which the
    quasi-inverse is computed. Defaults to \code{10000}.
  }
  \item{bws}{
    an optional \code{dbandwidth} specification (if already computed
    avoid unnecessary computation inside \code{npquantile}). This must
    be set as a \code{dbandwidth} object returned from an invocation of
    \code{\link{npudistbw}}. If not provided \code{\link{npudistbw}} is
    invoked with optional arguments passed via \code{\dots}.
  }
  \item{f}{
    an optional argument fed to \code{\link{extendrange}}. Defaults to
    \code{1}. See \code{?\link{extendrange}} for details.
  }
  \item{\dots}{
    additional arguments supplied to specify the bandwidth type, kernel
    types, bandwidth selection methods, and so on. See
    \code{?\link{npudistbw}} for details.
  }
  
}

\details{

  Typical usage is
  \preformatted{
    ## Start npRmpi for interactive execution. If slaves are already running and
    ## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
    ## reuse the existing pool instead of respawning. To change the number of
    ## slaves, call `npRmpi.stop(force=TRUE)` then restart.
    npRmpi.start(nslaves=1)
    set.seed(42)
    x <- rchisq(100,df=10)
    mpi.bcast.Robj2slave(x)
    mpi.bcast.cmd(npquantile(x),caller.execute=TRUE)
    npRmpi.stop()
  }

  The quantile function \eqn{q_\tau} is defined to be the
  left-continuous inverse of the distribution function \eqn{F(x)},
  i.e. \eqn{q_\tau = \inf\{x: F(x) \ge \tau\}}.

  A traditional estimator of \eqn{q_\tau} is the \eqn{\tau}th sample
  quantile. However, these estimates suffer from lack of efficiency
  arising from variability of individual order statistics; see Sheather
  and Marron (1990) and Hyndman and Fan (1996) for methods that
  interpolate/smooth the order statistics, each of which discussed in
  the latter can be invoked through \code{\link{quantile}} via
  \code{type=j}, \code{j=1,\dots,9}.
 
  The function \code{npquantile} implements a method for estimating
  smooth quantiles based on the quasi-inverse of a \code{\link{npudist}}
  object where \eqn{F(x)} is replaced with its kernel estimator and
  bandwidth selection is that appropriate for such objects; see
  Definition 2.3.6, page 21, Nelsen 2006 for a definition of the
  quasi-inverse of \eqn{F(x)}.

  For construction of the quasi-inverse we create a grid of evaluation
  points based on the function \code{\link{extendrange}} along with the
  sample quantiles themselves computed from invocation of
  \code{\link{quantile}}. The coarseness of the grid defined by
  \code{\link{extendrange}} (which has been passed the option
  \code{f=1}) is controlled by \code{num.eval}.
  
  Note that for any value of \eqn{\tau} less/greater than the
  smallest/largest value of \eqn{F(x)} computed for the evaluation data
  (i.e. that outlined in the paragraph above), the quantile returned for
  such values is that associated with the smallest/largest value of
  \eqn{F(x)}, respectively.

}
\value{
  
  \code{\link{npquantile}} returns a vector of quantiles corresponding
  to \code{tau}.
  
}
\references{


  Cheng, M.-Y. and Sun, S. (2006), \dQuote{Bandwidth selection for
  kernel quantile estimation,} Journal of the Chinese Statistical
  Association, \bold{44}, 271-295.
  
  Hyndman, R.J. and Fan, Y. (1996), \dQuote{Sample quantiles in
  statistical packages,} American Statistician, \bold{50}, 361-365.

  Li, Q. and J.S. Racine (2017), \dQuote{Smooth Unconditional Quantile
  Estimation,} Manuscript.

  Li, C. and H. Li and J.S. Racine (2017), \dQuote{Cross-Validated Mixed
  Datatype Bandwidth Selection for Nonparametric Cumulative
  Distribution/Survivor Functions,} Econometric Reviews, \bold{36},
  970-987.

  Nelsen, R.B. (2006), \emph{An Introduction to Copulas,} Second
  Edition, Springer-Verlag.

  Sheather, S. and J.S. Marron (1990), \dQuote{Kernel quantile
  estimators,} Journal of the American Statistical Association, Vol. 85,
  No. 410, 410-416.

  Yang, S.-S. (1985), \dQuote{A Smooth Nonparametric Estimator of a
  Quantile Function,} Journal of the American Statistical Association,
  \bold{80}, 1004-1011.
  
}
\section{Usage Issues}{
  Cross-validated bandwidth selection is used by default
  (\code{\link{npudistbw}}). For large datasets this can be
  computationally demanding. In such cases one might instead consider a
  rule-of-thumb bandwidth (\code{bwmethod="normal-reference"}) or,
  alternatively, use kd-trees (\code{options(np.tree=TRUE)} along with a
  bounded kernel (\code{ckertype="epanechnikov"})), both of which will
  reduce the computational burden appreciably.
}
\author{
    Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\seealso{

  \code{\link{quantile}} for various types of sample quantiles;
  \code{\link{ecdf}} for empirical distributions of which
  \code{\link{quantile}} is an inverse; \code{\link{boxplot.stats}} and
  \code{\link{fivenum}} for computing other versions of quartiles;
  \code{\link[logspline]{qlogspline}} for logspline density quantiles;
  \code{\link[ks]{qkde}} for alternative kernel quantiles, etc.

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

set.seed(42)

## Simulate data from a chi-square distribution
df <- 50
x <- rchisq(100,df=df)

## Vector of quantiles desired
tau <- c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)

mpi.bcast.Robj2slave(x)
mpi.bcast.Robj2slave(tau)

## Compute kernel smoothed sample quantiles
mpi.bcast.cmd(q <- npquantile(x,tau),
              caller.execute=TRUE)

q

## Compute sample quantiles using the default method in R (Type 7)
quantile(x,tau)

## True quantiles based on known distribution
qchisq(tau,df=df)

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
\keyword{ nonparametric }% at least one, from doc/KEYWORDS
