%
\name{npindex}
\alias{npindex}
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


\method{npindex}{default}(bws,
        txdat,
        tydat,
        nomad = FALSE,
        \dots) 

\method{npindex}{sibandwidth}(bws,
        txdat = stop("training data 'txdat' missing"),
        tydat = stop("training data 'tydat' missing"),
        exdat,
        eydat,
        boot.num = 399,
        errors = FALSE,
        gradients = FALSE,
        residuals = FALSE,
        \dots)
}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification, formula/data interface, and training data.
}

  \item{bws}{
    a bandwidth specification. This can be set as a
    \code{sibandwidth}
    object returned from an invocation of \code{npindexbw}, or
    as a vector of parameters (beta) with each element \eqn{i}
    corresponding to the coefficient for column \eqn{i} in \code{txdat}
    where the first element is normalized to 1, and a scalar bandwidth
    (h).
  }


  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npindexbw}} was called.
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


\subsection{Bandwidth Search Shortcut}{
This argument passes the recommended automatic local-polynomial NOMAD preset to \code{\link{npindexbw}} when bandwidths are computed inside \code{npindex}.
}

  \item{nomad}{
    logical shortcut passed through to \code{\link{npindexbw}} when
    bandwidths are computed inside \code{npindex}. When \code{TRUE},
    the single-index bandwidth route fills any missing values among
    \code{regtype}, \code{search.engine}, \code{degree.select},
    \code{bernstein.basis}, \code{degree.min}, \code{degree.max},
    \code{degree.verify}, and \code{bwtype} with the recommended
    automatic LP NOMAD preset documented in \code{\link{npindexbw}}.
    Additional NOMAD tuning arguments such as \code{nomad.nmulti} may
    also be supplied through \code{\dots}; \code{nmulti} remains the
    outer restart count while \code{nomad.nmulti} controls inner
    \code{crs::snomadr()} multistarts within each outer restart.
    After fitting, inspect \code{fit$bws$nomad.shortcut} on the
    returned object \code{fit} to see the normalized shortcut metadata.
  }


\subsection{Evaluation Data And Returned Quantities}{
These arguments control where the single-index fit is evaluated and which evaluation quantities are returned.
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


  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.
  }


  \item{y.eval}{
    If \code{newdata} contains dependent data and \code{y.eval = TRUE},
    \code{\link{npRmpi}} will compute goodness of fit statistics on these
    data and return them. Defaults to \code{FALSE}.
  }


\subsection{Fitted Quantities And Inference}{
These arguments control residuals, gradients, and bootstrap standard errors.
}

  \item{boot.num}{
    an integer specifying the number of bootstrap replications to use
    when performing standard error calculations. Defaults to
    \code{399}.
  }


  \item{errors}{
    a logical value indicating that you want (bootstrapped)
    standard errors for the conditional mean, gradients (when
    \code{gradients=TRUE} is set), and average gradients (when
    \code{gradients=TRUE} is set), computed and returned in the
    resulting \code{singleindex} object. Defaults to \code{FALSE}.
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


\subsection{Additional Arguments}{
Further arguments are passed to the bandwidth-selection counterpart when bandwidths are not supplied.
}

  \item{\dots}{
    additional arguments supplied to specify the parameters to the
    \code{sibandwidth} S3 method, which is called during estimation.
  }
}
\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, \code{\link{plot}}, \code{\link{plot.np}} for plotting options, and \code{\link{npRmpi.init}} for interactive/cluster MPI startup. See \code{\link{npRmpi.init}} details for performance tradeoffs (message passing/startup mode) and the \code{inst/Rprofile} manual-broadcast template.

  For S3 plotting help, use \code{methods("plot")} and query
  class-specific help topics such as \code{?plot.npregression} and
  \code{?plot.rbandwidth}. You can inspect implementations with
  \code{getS3method("plot","npregression")}.


  A matrix of gradients along with average derivatives are computed and
  returned if \code{gradients=TRUE} is used.

  For practitioners who want the recommended automatic LP NOMAD route
  without spelling out all LP tuning arguments,
  \code{npindex(..., nomad=TRUE)} and \code{npindexbw(..., nomad=TRUE)}
  expand missing settings to the same documented preset. Explicit
  incompatible settings fail fast rather than being silently rewritten.

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

\seealso{
  \code{\link{npRmpi.init}} for MPI startup and workflow guidance.
  \code{\link{np.kernels}}, \code{\link{np.options}},
  \code{\link{plot}}, \code{\link{plot.np}}, \code{\link{npindexbw}}
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
## vignette("npRmpi_getting_started", package = "npRmpi").

## Start npRmpi for interactive execution. If slaves are already running and
## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
## reuse the existing pool instead of respawning. To change the number of
## slaves, call `npRmpi.quit(force=TRUE)` then restart.

force.run <- nzchar(Sys.getenv("NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK"))

in.check <- !force.run && (
            nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
            nzchar(Sys.getenv("_R_CHECK_INTERNALS2_")) ||
            nzchar(Sys.getenv("_R_CHECK_CRAN_INCOMING_")) ||
            nzchar(Sys.getenv("_R_CHECK_EXAMPLE_TIMING_THRESHOLD_"))
          )

if (!in.check) {
npRmpi.init(nslaves=1)
set.seed(42)

n <- 500

x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)
y <- x1 - x2 + rnorm(n)

## Ichimura, continuous y

bw <- npindexbw(formula=y~x1+x2)

summary(bw)

model <- npindex(bws=bw,
                 gradients=TRUE)

summary(model)

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

## Note: on some systems (notably macOS+MPICH), repeatedly spawning and
## tearing down slaves in the same R session can lead to hangs/crashes.
## npRmpi may therefore keep slave daemons alive by default and
## `npRmpi.quit()` performs a "soft close". Use `force=TRUE` to
## actually shut down the slaves.
##
## You can disable reuse via `options(npRmpi.reuse.slaves=FALSE)` or by
## setting the environment variable `NP_RMPI_NO_REUSE_SLAVES=1` before
## loading the package.

npRmpi.quit()               ## soft close (may keep slaves alive)
## npRmpi.quit(force=TRUE)  ## hard close
} else {
  message("Skipping MPI spawn in check context (set NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK=1).")
}
} % enddontrun
}
\keyword{nonparametric}
