% $ID$
\name{npscoef}
\alias{npscoef}
\alias{npscoef.call}
\alias{npscoef.default}
\alias{npscoef.formula}
\alias{npscoef.scbandwidth}
\title{Smooth Coefficient Kernel Regression}

\description{
  \code{npscoef} computes a kernel regression estimate of a one (1)
  dimensional dependent variable on \eqn{p}-variate explanatory data,
  using the model \eqn{Y_i = W_{i}^{\prime} \gamma (Z_i) + u_i}{Y_i =
  t(W_i) * gamma(Z_i) + u_i} where
  \eqn{W_i'=(1,X_i')}{t(W_i)=(1,t(X_i))}, given a set of evaluation
  points, training points (consisting of explanatory data and dependent
  data), and a bandwidth specification. A bandwidth specification can be
  a \code{scbandwidth} object, or a bandwidth vector, bandwidth type and
  kernel type.
}

\usage{
npscoef(bws, \dots)

\method{npscoef}{formula}(bws, data = NULL, newdata = NULL, y.eval =
FALSE, \dots)

\method{npscoef}{call}(bws, \dots)

\method{npscoef}{default}(bws, txdat, tydat, tzdat, \dots)

\method{npscoef}{scbandwidth}(bws,
        txdat = stop("training data 'txdat' missing"),
        tydat = stop("training data 'tydat' missing"),
        tzdat = NULL,
        exdat,
        eydat,
        ezdat,
        residuals = FALSE,
        errors = TRUE,
        iterate = TRUE,
        maxiter = 100,
        tol = .Machine$double.eps,
        leave.one.out = FALSE,
        betas = FALSE,
        \dots)
}

\arguments{
  \item{bws}{
  a bandwidth specification. This can be set as a \code{scbandwidth}
  object returned from an invocation of \code{\link{npscoefbw}}, or
    as a vector of bandwidths, with each element \eqn{i} corresponding
    to the bandwidth for column \eqn{i} in \code{tzdat}. If specified as
    a vector
    additional arguments will need to be supplied as necessary to
    specify the bandwidth type, kernel types, training data, and so on.

  }

  \item{\dots}{
    additional arguments supplied to specify the regression type,
    bandwidth type, kernel types, selection methods, and so on.
    To do this, you may specify any of \code{bwscaling}, \code{bwtype}
    (one of \code{fixed}, \code{generalized_nn}, \code{adaptive_nn}),
    \code{ckertype}, \code{ckerorder}, as described in
    \code{\link{npscoefbw}}. 
    
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing
    the variables 
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
   \code{\link{npscoefbw}} was called.
  }

  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.  
  }

  \item{y.eval}{
    If \code{newdata} contains dependent data and \code{y.eval = TRUE},
    \code{\link{npRmpi}} will compute goodness of fit statistics on these
    data and return them. Defaults to \code{FALSE}.
  }

  \item{txdat}{
    a \eqn{p}-variate data frame of explanatory data (training data),
    which, by default, populates the columns \eqn{2} through \eqn{p+1}
    of \eqn{W} in the model equation, and in the 
    absence of \code{zdat}, will also correspond to 
    \eqn{Z} from the model equation. Defaults to
    the training data used to 
    compute the bandwidth object.
  }


  \item{tydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{txdat}. Defaults to
    the training data used to 
    compute the bandwidth object.
  }

  \item{tzdat}{
    an optionally specified \eqn{q}-variate data frame of explanatory data
    (training data), which corresponds to \eqn{Z} in the model
    equation. Defaults to 
    the training data used to 
    compute the bandwidth object.
  }

  \item{exdat}{
    a \eqn{p}-variate data frame of points on which the regression will be
    estimated (evaluation data).By default,
    evaluation takes place on the data provided by \code{txdat}.
  }

  \item{eydat}{
    a one (1) dimensional numeric or integer vector of the true values
    of the dependent variable. Optional, and used only to calculate the
    true errors. 
  }

  \item{ezdat}{
    an optionally specified \eqn{q}-variate data frame of points on
    which
    the regression will be estimated
    (evaluation data), which corresponds to \eqn{Z}
    in the model equation. Defaults to be the same as \code{txdat}.
  }

  \item{errors}{
    a logical value indicating whether or not asymptotic standard errors
    should be computed and returned in the resulting
    \code{smoothcoefficient} object. Defaults to \code{TRUE}.
  }

  \item{residuals}{
    a logical value indicating that you want residuals computed and
    returned in the resulting \code{smoothcoefficient} object. Defaults to
    \code{FALSE}.
  }


  \item{iterate}{
    a logical value indicating whether or not backfitted estimates
    should be iterated for self-consistency. Defaults to \code{TRUE}.
  }

  \item{maxiter}{
    integer specifying the maximum number of times to iterate the
    backfitted estimates while attempting to make the backfitted estimates
    converge to the desired tolerance. Defaults to \code{100}.
  }

  \item{tol}{
    desired tolerance on the relative convergence of backfit
    estimates. Defaults to \code{.Machine$double.eps}.
  }

  \item{leave.one.out}{
    a logical value to specify whether or not to compute the leave one
    out estimates. Will not work if \code{e[xyz]dat} is specified. Defaults to
    \code{FALSE}.
  }

  \item{betas}{
    a logical value indicating whether or not estimates of the
    components of \eqn{\gamma}{gamma} should be returned in the
    \code{smoothcoefficient} object along with the regression estimates.
    Defaults to \code{FALSE}.
  }
  
}

\value{
  \code{npscoef} returns a \code{smoothcoefficient} object.   The generic
  functions \code{\link{fitted}}, \code{\link{residuals}}, \code{\link{coef}},
  \code{\link{se}}, and \code{\link{predict}},
  extract (or generate) estimated values,
  residuals, coefficients, bootstrapped standard
  errors on estimates, and predictions, respectively, from
  the returned object. Furthermore, the functions \code{\link{summary}}
  and \code{\link{plot}} support objects of this type. The returned object
  has the following components:

  \item{eval}{ evaluation points }
  \item{mean}{ estimation of the regression function (conditional mean) at the
    evaluation points } 
  \item{merr}{ if \code{errors = TRUE}, standard errors of the
    regression estimates } 
  \item{beta}{ if \code{betas = TRUE}, estimates of the coefficients
    \eqn{\gamma}{gamma} at the 
    evaluation points }
  \item{resid}{ if \code{residuals = TRUE}, in-sample or out-of-sample
    residuals where appropriate (or possible) }
  \item{R2}{ coefficient of determination (Doksum and Samarov (1995))}
  \item{MSE}{ mean squared error }
  \item{MAE}{ mean absolute error }
  \item{MAPE}{ mean absolute percentage error }
  \item{CORR}{ absolute value of Pearson's correlation coefficient }
  \item{SIGN}{ fraction of observations where fitted and observed values
  agree in sign }
  
}
\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Cai Z. (2007), \dQuote{Trending time-varying coefficient time series
    models with serially correlated errors,} Journal of Econometrics,
  136, 163-188.

  Doksum, K. and A. Samarov (1995), \dQuote{Nonparametric estimation of
  global functionals and a measure of the explanatory power of
  covariates in regression,} The Annals of Statistics, 23 1443-1473.

  Hastie, T. and R. Tibshirani (1993), \dQuote{Varying-coefficient
    models,} Journal of the Royal Statistical Society, B 55, 757-796.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
    Theory and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2010), \dQuote{Smooth varying-coefficient
    estimation and inference for qualitative and quantitative data,}
    Econometric Theory, 26, 1-31.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.

 Li, Q. and D. Ouyang and J.S. Racine (2013), \dQuote{Categorical
 semiparametric varying-coefficient models,} Journal of Applied
 Econometrics, 28, 551-589.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth
    estimators for discrete distributions,}  Biometrika, 68, 301-309.

}

\author{
  Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.

  Support for backfitted bandwidths is experimental and is limited in
  functionality. The code does not support asymptotic standard errors
  or out of sample estimates with backfitting.
}

\seealso{
  \code{\link{bw.nrd}}, \code{\link{bw.SJ}}, \code{\link{hist}},
  \code{\link{npudens}}, \code{\link{npudist}},
  \code{\link{npudensbw}}, \code{\link{npscoefbw}} }

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

n <- 500

x <- runif(n)
z <- runif(n, min=-2, max=2)
y <- x*exp(z)*(1.0+rnorm(n,sd = 0.2))
mydat <- data.frame(x,y,z)
rm(x,y,z)

mpi.bcast.Robj2slave(mydat)

## A smooth coefficient model example

mpi.bcast.cmd(bw <- npscoefbw(y~x|z,data=mydat),
              caller.execute=TRUE)

summary(bw)

mpi.bcast.cmd(model <- npscoef(bws=bw, gradients=TRUE),
              caller.execute=TRUE)

summary(model)

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
\keyword{nonparametric}
