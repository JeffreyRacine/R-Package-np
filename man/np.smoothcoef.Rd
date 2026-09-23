%
\name{npscoef}
\alias{npscoef}
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

\method{npscoef}{formula}(bws, 
        data = NULL, 
        newdata = NULL, 
        y.eval = FALSE, 
        se = FALSE,
        \dots)


\method{npscoef}{default}(bws, 
        txdat, 
        tydat, 
        tzdat, 
        nomad = FALSE, 
        se = FALSE,
        \dots)

\method{npscoef}{scbandwidth}(bws,
        txdat = stop("training data 'txdat' missing"),
        tydat = stop("training data 'tydat' missing"),
        tzdat = NULL,
        exdat,
        eydat,
        ezdat,
        betas = FALSE,
        se = FALSE,
        iterate = TRUE,
        leave.one.out = FALSE,
        maxiter = 100,
        residuals = FALSE,
        tol = .Machine$double.eps,
        \dots)
}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification, formula/data interface, and smooth-coefficient training data.
}

  \item{bws}{
  a bandwidth specification. This can be set as a \code{scbandwidth}
  object returned from an invocation of \code{\link{npscoefbw}}, or
    as a vector of bandwidths, with each element \eqn{i} corresponding
    to the bandwidth for column \eqn{i} in \code{tzdat}. If specified as
    a vector
    additional arguments will need to be supplied as necessary to
    specify the bandwidth type, kernel types, training data, and so on.

  }



  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing
    the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
   \code{\link{npscoefbw}} was called.
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


\subsection{Local-Polynomial Degree And Bandwidth Search}{
This argument controls the recommended automatic local-polynomial NOMAD route, which jointly selects continuous polynomial degree and bandwidths when these are computed inside \code{npscoef}.
}

  \item{nomad}{
    logical shortcut passed through to \code{\link{npscoefbw}} when
    bandwidths are computed inside \code{npscoef}. When \code{TRUE},
    the smooth-coefficient bandwidth route fills any missing values
    among \code{regtype}, \code{search.engine},
    \code{degree.select}, \code{bernstein.basis},
    \code{degree.min}, \code{degree.max}, \code{degree.verify}, and
    \code{bwtype} with the recommended automatic local-polynomial degree-and-bandwidth NOMAD preset
    documented in \code{\link{npscoefbw}}. Additional NOMAD tuning
    arguments such as \code{nomad.nmulti} may also be supplied through
    \code{\dots}; \code{nmulti} remains the outer restart count while
    \code{nomad.nmulti} controls inner
    native \pkg{crs} NOMAD API multistarts within each outer restart. After fitting, inspect
    \code{fit$bws$nomad.shortcut} on the returned object \code{fit} to
    see the normalized shortcut metadata.
  }


\subsection{Evaluation Data And Returned Quantities}{
These arguments control where the smooth-coefficient fit is evaluated and which evaluation quantities are returned.
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
    in the model equation. With separate training \code{tzdat}, this defaults
    to \code{tzdat}; otherwise it follows \code{exdat}. Supplying only
    \code{ezdat} evaluates at the original training \code{txdat}. The two
    evaluation roles must have the same row count before omitting missing rows.
    Native \code{exdat} and \code{ezdat} columns correspond positionally to
    their training roles; different column labels do not change this mapping.
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


\subsection{Fitted Quantities And Backfitting}{
These arguments control returned coefficient estimates, errors, residuals, and iterative backfitting.
}

  \item{betas}{
    a logical value indicating whether or not estimates of the
    components of \eqn{\gamma}{gamma} should be returned in the
    \code{smoothcoefficient} object along with the regression estimates.
    Defaults to \code{FALSE}.
  }


  \item{se}{
    a logical value indicating whether or not asymptotic standard errors
    should be computed and returned in the resulting
    \code{smoothcoefficient} object. Defaults to \code{FALSE}.
  }


  \item{iterate}{
    a logical value indicating whether or not backfitted estimates
    should be iterated for self-consistency. Defaults to \code{TRUE}.
    Fitted objects retain this control, \code{maxiter}, \code{tol}, and
    \code{leave.one.out} for prediction and lazy residual extraction.
    Explicit controls supplied to \code{predict} override retained values.
  }


  \item{leave.one.out}{
    a logical value to specify whether or not to compute the leave one
    out estimates. Will not work if \code{e[xyz]dat} is specified. Defaults to
    \code{FALSE}.
  }


  \item{maxiter}{
    integer specifying the maximum number of times to iterate the
    backfitted estimates while attempting to make the backfitted estimates
    converge to the desired tolerance. Defaults to \code{100}.
  }


  \item{residuals}{
    a logical value indicating that you want residuals computed and
    returned in the resulting \code{smoothcoefficient} object. Defaults to
    \code{FALSE}. Residuals refer to the training sample, including when
    evaluation data are supplied, and do not require \code{se = TRUE}.
  }


  \item{tol}{
    desired tolerance on the relative convergence of backfit
    estimates. Defaults to \code{.Machine$double.eps}.
  }


\subsection{Additional Arguments}{
Further arguments are passed to the bandwidth-selection counterpart when bandwidths are not supplied.
}

  \item{\dots}{
    additional arguments supplied to specify the regression type,
    bandwidth type, kernel types, selection methods, and so on.
    To do this, you may specify any of \code{bwscaling}, \code{bwtype}
    (one of \code{fixed}, \code{generalized_nn}, \code{adaptive_nn}),
    \code{ckertype}, \code{ckerorder}, as described in
    \code{\link{npscoefbw}}.

  }
}
\details{
  In native calls, three unnamed arguments follow the fitter's bandwidth-first
  convention: \code{npscoef(bws, txdat, tydat)}; the smoothing data default
  to the training predictors. Numeric argument types cannot distinguish this
  from an intended data-first call. To select bandwidths with separate smoothing
  data, write \code{npscoef(xdat, ydat, zdat = zdat)} or name all data roles.
  The selector \code{npscoefbw(xdat, ydat, zdat)} retains its data-first
  convention. No argument roles are inferred from vector lengths or values.

  With a formula-derived bandwidth object, supplying only \code{tydat}
  replaces the evaluated response in retained-training-sample order: one row
  per observation retained after the original subset and missing-value selection.
  Original excluded rows are not restored; row names or equal lengths do not
  imply original-data order. The effective \code{na.action} is then applied
  jointly to the retained predictors and replacement response; any new omission
  indices refer to this retained sample. The selected bandwidths remain fixed,
  and the returned object retains the replacement sample for subsequent refits.
  Use \code{data} for full formula-data replacement, not together with a partial
  response replacement. Explicit native training arguments retain their separate
  positional contract.

  For native (non-formula) objects, named \code{newdata} columns are matched
  to the retained variable names and reordered as needed. Extra columns are
  ignored; missing, duplicate or ambiguous names are errors. An unnamed
  matrix or vector remains positional for a single evaluation role.
  Explicit native evaluation arguments retain their positional contract and
  take precedence over \code{newdata}.

  For a native bandwidth object, \code{newdata} is split into
  \code{exdat} and \code{ezdat} using the retained predictor and smoothing
  variable names. Names must identify both roles unambiguously;
  explicitly supplied native evaluation arguments take precedence. This
  evaluates at the requested rows without repeating bandwidth selection.

}
\value{
  When an external evaluation row has no usable kernel support or its
  kernel-weighted system cannot be solved, the corresponding fitted value,
  coefficients and requested standard errors are \code{NA}. One warning
  identifies affected rows; valid rows are retained in their original order.
  A larger bandwidth or another kernel may provide usable support.
  Required training and bootstrap fits instead signal an error on this
  failure; they do not silently discard rows. Being outside the training
  range alone does not imply that a local fit is undefined.

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
  \item{merr}{ if \code{se = TRUE}, standard errors of the
    regression estimates }
  \item{beta}{ if \code{betas = TRUE}, estimates of the coefficients
    \eqn{\gamma}{gamma} at the
    evaluation points }
  \item{grad}{ estimated derivatives of the conditional mean with
    respect to the regressors in \code{xdat}; these correspond to the
    non-intercept smooth coefficient estimates at each evaluation point}
  \item{gerr}{ if \code{se = TRUE}, asymptotic standard errors for
    \code{grad} }
  \item{resid}{ if \code{residuals = TRUE}, training responses minus their
    fitted training values, restored to training rows under \code{na.exclude}.
    Evaluation responses, when supplied, are used for evaluation goodness-of-fit
    summaries, not to redefine these residuals. }
  \item{R2}{ coefficient of determination (Doksum and Samarov (1995))}
  \item{MSE}{ mean squared error }
  \item{MAE}{ mean absolute error }
  \item{MAPE}{ mean absolute percentage error }
  \item{CORR}{ absolute value of Pearson's correlation coefficient }
  \item{SIGN}{ fraction of observations where fitted and observed values
  agree in sign }

}
\section{Book And Method Pointers}{
  Smooth-coefficient moment solves use their own bounded numerical ridging;
  see \code{\link{npreg}} for its Seifert--Gasser motivation.
  The smooth-coefficient covariance is a separate moment calculation, not
  the ordinary-regression residual-normalized sandwich.

  The smooth-coefficient model lets slopes vary with conditioning
  variables, typically written
  \eqn{Y=X^\prime\beta(Z)+\epsilon}. The functions estimate the
  coefficient functions \eqn{\beta_j(z)} using mixed-data kernel
  smoothing over \eqn{Z}, with the reported fitted values obtained by
  combining the estimated coefficient functions with the corresponding
  columns of \eqn{X}.

  For book-length derivations, see Li and Racine (2007), Chapter 9
  \emph{Additive and Smooth (Varying) Coefficient Semiparametric Models},
  especially Sections 9.3-9.3.4, and Racine (2019), Chapter 8
  \emph{Semiparametric Conditional Mean Function Estimation}, especially
  the varying-coefficient material.
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

  For formula fits, character-valued missing-data policies are resolved
  through the \code{stats} model-frame owner; pass a caller-local policy as
  a function object rather than its name. With replacement training data,
  an explicit \code{na.action}, including \code{NULL}, takes precedence
  over the saved policy. Native fitting still requires complete rows, so
  a permissive policy does not guarantee the same retained sample as
  \code{stats::model.frame}.

  When reusing a bandwidth object without explicit training or evaluation
  data, newly constructed bandwidth objects use retained training values.
  A fit with explicit replacement training data retains that sample in its
  returned bandwidth object; evaluation grids do not replace training data.
  Older objects without retained samples use their stored call environment
  when available and may require explicit training inputs after
  \code{saveRDS()} and \code{readRDS()}. Custom formula transformations used
  on new evaluation data must still be available.

  For practitioners who want the recommended automatic local-polynomial degree-and-bandwidth NOMAD route
  without spelling out all LP tuning arguments,
  \code{npscoef(..., nomad=TRUE)} and \code{npscoefbw(..., nomad=TRUE)}
  expand missing settings to the same documented preset. Explicit
  incompatible settings fail fast rather than being silently rewritten.

  For plotting options for fitted smooth-coefficient objects and their
  bandwidth objects, see \code{\link{plot.np}}.

  Support for backfitted bandwidths is experimental and is limited in
  functionality. The code does not support asymptotic standard errors
  or out of sample estimates with backfitting.
}

\seealso{
  \code{\link{npRmpi.init}},
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}, \code{\link{plot.np}},
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

x <- runif(n)
z <- runif(n, min=-2, max=2)
y <- x*exp(z)*(1.0+rnorm(n,sd = 0.2))

## A smooth coefficient model example

bw <- npscoefbw(y~x|z)

summary(bw)

model <- npscoef(bws=bw)

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
