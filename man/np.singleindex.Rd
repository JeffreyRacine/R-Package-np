% $ID$
\name{npindex}
\alias{npindex}
\alias{npindex.call}
\alias{npindex.default}
\alias{npindex.formula}
\alias{npindex.sibandwidth}

\title{Semiparametric Single Index Model}

\description{
  \code{npindex} computes a semiparametric single index model
  for a dependent variable and \eqn{p}-variate explanatory data using
  the model \eqn{Y = G(X\beta) + \epsilon}{Y = G(XB) + epsilon}, given a
  set of evaluation points, training points (consisting of explanatory
  data and dependent data), and a \code{npindexbw} bandwidth
  specification. Note that for this semiparametric estimator, the
  bandwidth object contains parameters for the single index model and
  the (scalar) bandwidth for the index function.
}

\usage{
npindex(bws, \dots)

\method{npindex}{formula}(bws,
        data = NULL,
        newdata = NULL,
        y.eval = FALSE,
        \dots)

\method{npindex}{call}(bws,
        \dots)

\method{npindex}{default}(bws,
        txdat,
        tydat,
        \dots) 

\method{npindex}{sibandwidth}(bws,
        txdat = stop("training data 'txdat' missing"),
        tydat = stop("training data 'tydat' missing"),
        exdat,
        eydat,
        gradients = FALSE,
        residuals = FALSE,
        errors = FALSE,
        boot.num = 399,
        \dots)
}

\arguments{

  \item{bws}{
    a bandwidth specification. This can be set as a
    \code{sibandwidth}
    object returned from an invocation of \code{npindexbw}, or
    as a vector of parameters (beta) with each element \eqn{i}
    corresponding to the coefficient for column \eqn{i} in \code{txdat}
    where the first element is normalized to 1, and a scalar bandwidth
    (h). 
  }

  \item{gradients}{
    a logical value indicating that you want gradients and the
    asymptotic covariance matrix for beta computed and returned in the
    resulting \code{singleindex} object. Defaults to \code{FALSE}.
  }

  \item{residuals}{
    a logical value indicating that you want residuals computed and
    returned in the resulting \code{singleindex} object. Defaults to
    \code{FALSE}.
  }

  \item{errors}{
    a logical value indicating that you want (bootstrapped)
    standard errors for the conditional mean, gradients (when
    \code{gradients=TRUE} is set), and average gradients (when
    \code{gradients=TRUE} is set), computed and returned in the
    resulting \code{singleindex} object. Defaults to \code{FALSE}.
  }

  \item{boot.num}{
    an integer specifying the number of bootstrap replications to use
    when performing standard error calculations. Defaults to
    \code{399}. 
  }

  \item{\dots}{
    additional arguments supplied to specify the parameters to the
    \code{sibandwidth} S3 method, which is called during estimation. 
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npindexbw}} was called.
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
    a \eqn{p}-variate data frame of explanatory data (training data) used to
    calculate the regression estimators. Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{tydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{txdat}. Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{exdat}{
    a \eqn{p}-variate data frame of points on which the regression will be
    estimated (evaluation data). By default,
    evaluation takes place on the data provided by \code{txdat}.
  }

  \item{eydat}{
    a one (1) dimensional numeric or integer vector of the true values
    of the dependent variable. Optional, and used only to calculate the
    true errors. 
  }
  
}

\details{

  A matrix of gradients along with average derivatives are computed and
  returned if \code{gradients=TRUE} is used.

}

\value{

  \code{npindex} returns a \code{npsingleindex} object.  The generic
  functions \code{\link{fitted}}, \code{\link{residuals}},
  \code{\link{coef}}, \code{\link{vcov}}, \code{\link{se}},
  \code{\link{predict}}, and \code{\link{gradients}}, extract (or
  generate) estimated values, residuals, coefficients,
  variance-covariance matrix, bootstrapped standard errors on estimates,
  predictions, and gradients, respectively, from the returned
  object. Furthermore, the functions \code{\link{summary}} and
  \code{\link{plot}} support objects of this type. The returned object
  has the following components:

  \item{eval}{ evaluation points }
  \item{mean}{ estimates of the regression function (conditional mean) at the
  evaluation points } 
  \item{beta}{ the model coefficients }
  \item{betavcov}{ the asymptotic covariance matrix for the model coefficients}
  \item{merr}{ standard errors of the regression function estimates } 
  \item{grad}{ estimates of the gradients at each evaluation point }
  \item{gerr}{ standard errors of the gradient estimates }
  \item{mean.grad}{ mean (average) gradient over the evaluation points}
  \item{mean.gerr}{ bootstrapped standard error of the mean gradient estimates }
  \item{R2}{ if \code{method="ichimura"}, coefficient of determination
  (Doksum and Samarov (1995))}
  \item{MSE}{if \code{method="ichimura"}, mean squared error }
  \item{MAE}{if \code{method="ichimura"}, mean absolute error }
  \item{MAPE}{if \code{method="ichimura"}, mean absolute percentage error }
  \item{CORR}{if \code{method="ichimura"}, absolute value of Pearson's correlation coefficient }
  \item{SIGN}{if \code{method="ichimura"}, fraction of observations where fitted and observed values
    agree in sign }
  \item{confusion.matrix}{if \code{method="kleinspady"}, the confusion matrix or \code{NA} if outcomes
    are not available }
  \item{CCR.overall}{if \code{method="kleinspady"},  the overall correct
    classification ratio, or \code{NA} if outcomes are not available }
  \item{CCR.byoutcome}{if \code{method="kleinspady"},  a numeric vector containing the correct
    classification ratio by outcome, or \code{NA} if outcomes are not
    available }
  \item{fit.mcfadden}{if \code{method="kleinspady"},  the McFadden-Puig-Kerschner performance measure
    or \code{NA} if outcomes are not available }
}



\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Doksum, K. and A. Samarov (1995), \dQuote{Nonparametric estimation of
  global functionals and a measure of the explanatory power of
  covariates  regression,} The Annals of Statistics, 23 1443-1473.

  Ichimura, H., (1993), \dQuote{Semiparametric least squares (SLS) and
  weighted SLS estimation of single-index models,} Journal of
  Econometrics, 58, 71-120.

  Klein, R. W. and R. H. Spady (1993), \dQuote{An efficient semiparametric
  estimator for binary response models,} Econometrica, 61, 387-421.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

  McFadden, D. and C. Puig and D. Kerschner (1977), \dQuote{Determinants
    of the long-run demand for electricity,} Proceedings of the
  American Statistical Association (Business and Economics Section),
  109-117. 

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth
    estimators for discrete distributions,} Biometrika, 68, 301-309.
}

    

\author{
  Tristen Hayfield  \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.

  \code{\link{vcov}} requires that \code{gradients=TRUE} be set.
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

n <- 500

x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)
y <- x1 - x2 + rnorm(n)
mydat <- data.frame(x1,x2,y)
rm(y,x1,x2)

mpi.bcast.Robj2slave(mydat)

## Ichimura, continuous y

mpi.bcast.cmd(bw <- npindexbw(formula=y~x1+x2,
                              data=mydat),
              caller.execute=TRUE)

summary(bw)

mpi.bcast.cmd(model <- npindex(bws=bw,
                               gradients=TRUE),
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
