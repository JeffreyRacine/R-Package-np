%
\name{npindexbw}
\alias{npindexbw}
\alias{npindexbw.default}
\alias{npindexbw.formula}
\alias{npindexbw.sibandwidth}
\title{Semiparametric Single Index Model Parameter and Bandwidth Selection}

\description{
  \code{npindexbw} computes a \code{npindexbw} bandwidth specification
  using the model \eqn{Y = G(X\beta) + \epsilon}{Y = G(XB) +
  epsilon}. For continuous \eqn{Y}, the approach is that of Hardle, Hall
  and Ichimura (1993) which jointly minimizes a least-squares
  cross-validation function with respect to the parameters and
  bandwidth. For binary \eqn{Y}, a likelihood-based cross-validation
  approach is employed which jointly maximizes a likelihood
  cross-validation function with respect to the parameters and
  bandwidth. The bandwidth object contains parameters for the single
  index model and the (scalar) bandwidth for the index function.
}

\usage{
npindexbw(\dots)

\method{npindexbw}{formula}(formula, 
          data, 
          subset, 
          na.action, 
          call, 
          \dots)

\method{npindexbw}{default}(xdat = stop("training data xdat missing"),
          ydat = stop("training data ydat missing"),
          bws,
          bandwidth.compute = TRUE,
          basis = c("glp", "additive", "tensor"),
          bernstein.basis = FALSE,
          degree = NULL,
          degree.select = c("manual", "coordinate", "exhaustive"),
          search.engine = c("nomad+powell", "cell", "nomad"),
          nomad = FALSE,
          nomad.nmulti = 0L,
          degree.min = NULL,
          degree.max = NULL,
          degree.start = NULL,
          degree.restarts = 0L,
          degree.max.cycles = 20L,
          degree.verify = FALSE,
          nmulti,
          only.optimize.beta,
          optim.abstol,
          optim.maxattempts,
          optim.maxit,
          optim.method,
          optim.reltol,
          random.seed,
          regtype = c("lc", "ll", "lp"),
          scale.factor.init.lower = 0.1,
          scale.factor.init.upper = 2.0,
          scale.factor.init = 0.5,
          scale.factor.search.lower = NULL,
          \dots)

\method{npindexbw}{sibandwidth}(xdat = stop("training data xdat missing"),
          ydat = stop("training data ydat missing"),
          bws,
          bandwidth.compute = TRUE,
          nmulti,
          only.optimize.beta = FALSE,
          optim.abstol = .Machine$double.eps,
          optim.maxattempts = 10,
          optim.maxit = 500,
          optim.method = c("Nelder-Mead", "BFGS", "CG"),
          optim.reltol = sqrt(.Machine$double.eps),
          random.seed = 42,
          scale.factor.init.lower = 0.1,
          scale.factor.init.upper = 2.0,
          scale.factor.init = 0.5,
          scale.factor.search.lower = NULL,
          \dots)

}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the data, formula interface, method label, and whether bandwidths are supplied or computed.
}

  \item{bandwidth.compute}{
    a logical value which specifies whether to do a numerical search for
    bandwidths or not. If set to \code{FALSE}, a \code{bandwidth} object
    will be returned with bandwidths set to those specified
    in \code{bws}. Defaults to \code{TRUE}.
  }

  \item{bws}{
    a bandwidth specification. This can be set as a
    \code{singleindexbandwidth}
    object returned from an invocation of \code{npindexbw}, or
    as a vector of parameters (beta) with each element \eqn{i}
    corresponding to the coefficient for column \eqn{i} in \code{xdat}
    where the first element is normalized to 1, and a scalar bandwidth
    (h). If specified as a vector, then additional arguments will need
    to be supplied as necessary to specify the bandwidth type, kernel
    types, and so on.
  }

  \item{call}{
    the original function call. This is passed internally by
    \code{\link{npRmpi}} when a bandwidth search has been implied by a call to
    another function. It is not recommended that the user set this.
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing
    the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which the
    function is called.
  }

  \item{formula}{
    a symbolic description of variables on which bandwidth selection is
    to be performed. The details of constructing a formula are
    described below.
  }

  \item{method}{
    the single index model method, one of either \dQuote{ichimura}
    (Ichimura (1993)) or \dQuote{kleinspady} (Klein and Spady
    (1993)). Defaults to
    \code{ichimura}.
  }

  \item{na.action}{
    a function which indicates what should happen when the data contain
    \code{NA}s. The default is set by the
    \code{\link{na.action}} setting of options, and is
    \code{\link{na.fail}} if that is unset. The (recommended) default is
    \code{\link{na.omit}}.
  }

  \item{subset}{
    an optional vector specifying a subset of observations to be used in
    the fitting process.
  }

  \item{xdat}{
    a \eqn{p}-variate data frame of explanatory data (training data) used to
    calculate the regression estimators.
  }

  \item{ydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{xdat}.
  }


\subsection{Automatic Degree Search Controls}{
These arguments control automatic local-polynomial degree search.
}

  \item{degree.max}{
    optional scalar or integer vector giving upper bounds for automatic
    degree search when \code{degree.select != "manual"}.
  }

  \item{degree.max.cycles}{
    positive integer giving the maximum number of coordinate-search
    sweeps over the degree vector. Ignored for \code{"manual"} and
    \code{"exhaustive"} degree selection.
  }

  \item{degree.min}{
    optional scalar or integer vector giving lower bounds for automatic
    degree search when \code{degree.select != "manual"}.
  }

  \item{degree.restarts}{
    non-negative integer giving the number of additional deterministic
    coordinate-search restarts. Ignored for \code{"manual"} and
    \code{"exhaustive"} degree selection.
  }

  \item{degree.select}{
    character string controlling local-polynomial degree handling when
    \code{regtype="lp"}. \code{"manual"} (default) treats
    \code{degree} as fixed. \code{"coordinate"} performs cached
    coordinate-wise search over admissible degree values for the index
    smoother. \code{"exhaustive"} evaluates the full admissible degree
    grid when \code{search.engine="cell"}. For NOMAD-based search
    engines, any non-\code{"manual"} value requests direct joint search
    over the degree and bandwidth coordinates.
  }

  \item{degree.start}{
    optional starting degree vector for automatic coordinate search. If
    omitted, the search starts from the degree-zero local-constant
    baseline for the index smoother.
  }

  \item{degree.verify}{
    logical value indicating whether a coordinate-search solution should
    be exhaustively verified over the admissible degree grid after the
    heuristic phase completes. Available only for
    \code{search.engine="cell"}.
  }


\subsection{Continuous Scale-Factor Search Initialization}{
These controls define deterministic and random continuous scale-factor starts and the lower admissibility floor for fixed-bandwidth search.
}

  \item{scale.factor.init}{
    deterministic initial scale factor for continuous fixed-bandwidth
    search. Defaults to \code{0.5}. The value supplied by the user is
    not rewritten, but the effective first start passed to the optimizer
    is \code{max(scale.factor.init, scale.factor.search.lower)}. See
    Details.
  }
  \item{scale.factor.init.lower}{
    lower endpoint for random continuous scale-factor starts. Defaults
    to \code{0.1}. The value supplied by the user is not rewritten, but
    the effective random-start lower endpoint is
    \code{max(scale.factor.init.lower, scale.factor.search.lower)}. See
    Details.
  }
  \item{scale.factor.init.upper}{
    upper endpoint for random continuous scale-factor starts. Defaults
    to \code{2.0}. It must be greater than or equal to the effective
    lower endpoint,
    \code{max(scale.factor.init.lower, scale.factor.search.lower)};
    otherwise bandwidth search errors rather than silently expanding the
    interval. See Details.
  }
  \item{scale.factor.search.lower}{
    optional nonnegative scalar giving the hard lower admissibility
    bound for continuous fixed-bandwidth search candidates. Defaults to
    \code{NULL}. If \code{NULL}, an existing bandwidth object's stored
    value is inherited when available; otherwise the package default
    \code{0.1} is used. This floor applies to computed/search bandwidth
    candidates and to effective search starts only. It does not rewrite
    explicit bandwidths supplied for storage with
    \code{bandwidth.compute = FALSE}. Final fixed-bandwidth search
    candidates must also have a finite valid raw objective value.
  }

\subsection{Local-Polynomial Model Specification}{
These arguments control the index smoother, local-polynomial basis, and fixed degree specification.
}

  \item{basis}{
    local polynomial basis selector used when
    \code{regtype="lp"}: one of \code{"glp"}, \code{"additive"}, or
    \code{"tensor"}. Ignored for \code{"lc"} and \code{"ll"}.
  }

  \item{bernstein.basis}{
    logical flag used when \code{regtype="lp"}; if \code{TRUE}, use a
    Bernstein/B-spline basis for local polynomial terms. When automatic
    degree search is requested and \code{bernstein.basis} is not
    explicitly supplied, the search route defaults to \code{TRUE} for
    numerical stability. Explicit \code{bernstein.basis=FALSE} is
    honored, but raw-polynomial search can be poorly conditioned at
    higher degrees.
  }

  \item{degree}{
    integer degree vector for continuous predictors when
    \code{regtype="lp"}. When \code{degree.select="manual"}, this must
    be supplied explicitly. For single-index regression this is
    typically length one.
  }

  \item{regtype}{
    a character string specifying local smoothing type for the
    nonparametric index regression fit used downstream in
    \code{npindex}. Supported values are \code{"lc"}, \code{"ll"}, and
    \code{"lp"}.
  }


\subsection{NOMAD Search Controls}{
These arguments control the optional NOMAD direct-search route for local-polynomial degree and bandwidth search.
}

  \item{nomad}{
    logical shortcut for the recommended automatic local-polynomial
    NOMAD route. When \code{TRUE}, any missing values among
    \code{regtype}, \code{search.engine}, \code{degree.select},
    \code{bernstein.basis}, \code{degree.min}, \code{degree.max},
    \code{degree.verify}, and \code{bwtype} are filled with
    \code{regtype="lp"}, \code{search.engine="nomad+powell"},
    \code{degree.select="coordinate"},
    \code{bernstein.basis=TRUE}, \code{degree.min=0L},
    \code{degree.max=10L}, \code{degree.verify=FALSE}, and
    \code{bwtype="fixed"}. Explicit incompatible settings error
    immediately; in particular, \code{nomad=TRUE} currently requires
    \code{regtype="lp"}, \code{bwtype="fixed"}, automatic degree
    search, \code{bernstein.basis=TRUE}, no explicit \code{degree},
    and \code{search.engine \%in\% c("nomad", "nomad+powell")}.
    This shortcut does not change the meaning of \code{nmulti} or
    \code{nomad.nmulti}: \code{nmulti} remains the outer restart count,
    while \code{nomad.nmulti} controls inner
    \code{crs::snomadr()} multistarts within each outer restart.
    Returned bandwidth objects retain this normalized preset metadata in
    \code{bw$nomad.shortcut} for a returned object \code{bw}; when
    available, \code{nomad.time} and \code{powell.time} record the
    direct-search and Powell-polish timing components.
  }

  \item{nomad.nmulti}{
    non-negative integer controlling the inner
    \code{crs::snomadr()} multistart count used within each outer
    NOMAD restart when \code{regtype="lp"} and automatic degree search
    uses \code{search.engine="nomad"} or \code{"nomad+powell"}.
    Defaults to \code{0L}, which preserves the current one-start-per-
    restart behavior. This does not replace \code{nmulti}:
    \code{nmulti} controls outer restarts, while \code{nomad.nmulti}
    controls inner NOMAD multistarts within each outer restart.
  }

  \item{search.engine}{
    character string controlling the automatic local-polynomial search
    backend when \code{regtype="lp"} and
    \code{degree.select != "manual"}. \code{"nomad+powell"} (default)
    performs direct joint search over the index-smoother bandwidth and
    degree using \code{crs::snomadr()}, then applies one Powell hot
    start from the NOMAD solution. \code{"nomad"} omits the Powell
    refinement. \code{"cell"} profiles the criterion over the
    admissible degree grid using repeated fixed-degree solves.
    NOMAD-based search currently requires \code{degree.verify=FALSE}
    and the suggested package \pkg{crs} to be installed.
  }


\subsection{Numerical Search Controls}{
These arguments control outer restart behavior for bandwidth and index-parameter search.
}

  \item{nmulti}{
    integer number of times to restart the process of finding extrema
    of the cross-validation function from different (random) initial
    points. Defaults to \code{min(2,ncol(xdat))}.
  }


\subsection{Optimization Controls}{
These arguments control outer optimization behavior for the semiparametric search.
}

   \item{only.optimize.beta}{
     signals the routine to only minimize the objective function with
     respect to beta
   }

  \item{optim.abstol}{
    the absolute convergence tolerance used by \code{\link{optim}}. Only useful
    for non-negative functions, as a tolerance for reaching
    zero. Defaults to \code{.Machine$double.eps}.
  }

  \item{optim.maxattempts}{
    maximum number of attempts taken trying to achieve successful
    convergence in \code{\link{optim}}. Defaults to \code{100}.
  }

  \item{optim.maxit}{
    maximum number of iterations used by \code{\link{optim}}. Defaults
     to \code{500}.
   }

  \item{optim.method}{ method used by \code{\link{optim}} for minimization of
    the objective function. See \code{?optim} for references. Defaults
    to \code{"Nelder-Mead"}.

    the default method is an implementation of that of Nelder and Mead
    (1965), that uses only function values and is robust but relatively
    slow.  It will work reasonably well for non-differentiable
    functions.

    method \code{"BFGS"} is a quasi-Newton method (also known as a
    variable metric algorithm), specifically that published
    simultaneously in 1970 by Broyden, Fletcher, Goldfarb and Shanno.
    This uses function values and gradients to build up a picture of the
    surface to be optimized.

    method \code{"CG"} is a conjugate gradients method based
    on that by Fletcher and Reeves (1964) (but with the option of
    Polak-Ribiere or Beale-Sorenson updates).  Conjugate gradient
    methods will generally be more fragile than the BFGS method, but as
    they do not store a matrix they may be successful in much larger
    optimization problems.
  }

  \item{optim.reltol}{
    relative convergence tolerance used by \code{\link{optim}}.  The algorithm
    stops if it is unable to reduce the value by a factor of 'reltol *
    (abs(val) + reltol)' at a step.  Defaults to
    \code{sqrt(.Machine$double.eps)}, typically about \code{1e-8}.
  }

  \item{random.seed}{
    an integer used to seed R's random number generator. This ensures
    replicability of the numerical search. Defaults to 42.
  }


\subsection{Additional Arguments}{
These arguments collect remaining controls passed through S3 methods.
}

  \item{\dots}{
    additional arguments supplied to specify the parameters to the
    \code{sibandwidth} S3 method, which is called during the numerical
    search. In particular, \code{bwtype} may be supplied here to request
    \code{"fixed"}, \code{"generalized_nn"}, or \code{"adaptive_nn"} bandwidth
    types.
  }




}

\details{

The \code{scale.factor.*} controls are dimensionless search
controls. The package converts scale factors to bandwidths using the
estimator-specific scaling encoded in the bandwidth object, including
kernel order and the number of continuous variables relevant for the
estimator. Users should not pre-multiply these controls by sample-size
or standard-deviation factors.

  \code{scale.factor.init} controls the deterministic first search
  start when that control is exposed. \code{scale.factor.init.lower}
  and \code{scale.factor.init.upper} define the random multistart
  interval when exposed. \code{scale.factor.search.lower} is the lower
  admissibility bound for continuous fixed-bandwidth search candidates.
  The effective first start is \code{max(scale.factor.init,
  scale.factor.search.lower)} when both controls are present, and the
  effective random-start lower endpoint is
  \code{max(scale.factor.init.lower, scale.factor.search.lower)}.
  \code{scale.factor.init.upper} must be at least that effective lower
  endpoint; the package errors rather than silently expanding the user's
  interval.

  When \code{scale.factor.search.lower} is \code{NULL}, an existing
  bandwidth object's stored floor is inherited when available;
  otherwise the package default \code{0.1} is used. Explicit bandwidths
  supplied for storage with \code{bandwidth.compute = FALSE} are not
  rewritten by the search floor.

  Categorical search-start controls such as \code{dfac.init},
  \code{lbd.init}, and \code{hbd.init} have separate semantics and are
  not affected by \code{scale.factor.search.lower}.

  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, \code{\link{plot}} for plotting options, and \code{\link{npRmpi.init}} for interactive/cluster MPI startup. See \code{\link{npRmpi.init}} details for performance tradeoffs (message passing/startup mode) and the \code{inst/Rprofile} manual-broadcast template.

  For S3 plotting help, use \code{methods("plot")} and query
  class-specific help topics such as \code{?plot.npregression} and
  \code{?plot.rbandwidth}. You can inspect implementations with
  \code{getS3method("plot","npregression")}.

  We implement Ichimura's (1993) method via
  joint estimation of the bandwidth and coefficient vector using
  leave-one-out nonlinear least squares. We implement Klein and Spady's
  (1993) method maximizing the leave-one-out log likelihood function
  jointly with respect to the bandwidth and coefficient vector. Note
  that Klein and Spady's (1993) method is for \emph{binary outcomes
    only}, while Ichimura's (1993) method can be applied for any outcome
  data type (i.e., continuous or discrete).

  We impose the identification condition that the first element of the
  coefficient vector beta is equal to one, while identification also
  requires that the explanatory variables contain \emph{at least one}
  continuous variable.

  \code{npindexbw} may be invoked \emph{either} with a formula-like
  symbolic description of variables on which bandwidth selection is to
  be performed \emph{or} through a simpler interface whereby data is
  passed directly to the function via the \code{xdat} and \code{ydat}
  parameters. Use of these two interfaces is \bold{mutually exclusive}.

  Note that, unlike most other bandwidth methods in the \code{npRmpi}
  package, this implementation uses the R \code{\link{optim}} nonlinear
  minimization routines and \code{\link{npksum}}. We have implemented
  multistarting and strongly encourage its use in practice. For
  exploratory purposes, you may wish to override the default search
  tolerances, say, setting \code{optim.reltol=.1} and conduct
  multistarting (the default is to restart min(2, ncol(xdat)) times) as is done
  for a number of examples.

  Data for which bandwidths are to be estimated may be specified
  symbolically. A typical description has the form \code{dependent data
  ~ explanatory data}, where \code{dependent data} is a univariate
  response, and \code{explanatory data} is a series of variables
  specified by name, separated by the separation character '+'. For
  example \code{ y1 ~ x1 + x2 } specifies that the bandwidth object for
  the regression of response \code{y1} and semiparametric regressors
  \code{x1} and \code{x2} are to be estimated.  See below for further
  examples.

  When \code{regtype="lp"} and \code{degree.select != "manual"},
  \code{npindexbw} can jointly determine the local-polynomial degree for
  the index smoother together with its bandwidth coordinate. With
  \code{search.engine="cell"}, the criterion is profiled over the
  admissible degree grid using cached coordinate-wise or exhaustive
  search. With \code{search.engine="nomad"} or
  \code{"nomad+powell"}, the criterion is optimized directly over the
  joint degree/bandwidth space using \code{crs::snomadr()};
  \code{"nomad+powell"} then performs one Powell hot start and retains
  the better of the direct NOMAD and polished solutions. For the
  index-smoother local-polynomial component, this polynomial-adaptive
  joint-search route follows Hall and Racine (2015).

  Setting \code{nomad=TRUE} is a convenience preset for this automatic
  LP route, not a generic optimizer alias. For single-index bandwidth
  selection it expands any missing values to the equivalent long-form
  call
  \preformatted{
    npindexbw(...,
              regtype = "lp",
              search.engine = "nomad+powell",
              degree.select = "coordinate",
              bernstein.basis = TRUE,
              degree.min = 0L,
              degree.max = 10L,
              degree.verify = FALSE,
              bwtype = "fixed")
  }
  Compatible explicit tuning arguments are respected. Incompatible
  explicit settings fail fast so the shortcut never silently changes
  user-selected semantics.

}

\value{
  \code{npindexbw} returns a \code{sibandwidth} object, with the
  following components:

  \item{bw}{ bandwidth(s), scale factor(s) or nearest neighbours for the
    data, \code{xdat} }
  \item{beta}{ coefficients of the model }
  \item{fval}{ objective function value at minimum }

  If \code{bwtype} is set to \code{fixed}, an object containing a scalar
  bandwidth for the function \eqn{G(X\beta)}{G(XB)} and an estimate of
  the parameter vector \eqn{\beta}{B} is returned.

  If \code{bwtype} is set to \code{generalized_nn} or
  \code{adaptive_nn}, then instead the scalar \eqn{k}th nearest neighbor
  is returned.

  The functions \code{\link{coef}}, \code{\link{predict}},
  \code{\link{summary}}, and  \code{\link{plot}} support
  objects of this class.
}

\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Hardle, W. and P. Hall and H. Ichimura (1993), \dQuote{Optimal
  Smoothing in Single-Index Models,} The Annals of Statistics, 21,
  157-178.

  Ichimura, H., (1993), \dQuote{Semiparametric least squares (SLS) and
  weighted SLS estimation of single-index models,} Journal of
  Econometrics, 58, 71-120.

  Klein, R. W. and R. H. Spady (1993), \dQuote{An efficient semiparametric
  estimator for binary response models,} Econometrica, 61, 387-421.

  Hall, P. and J.S. Racine (2015), \dQuote{Infinite Order Cross-Validated
  Local Polynomial Regression,} Journal of Econometrics, 185, 510-525.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth
    estimators for discrete distributions,} Biometrika, 68, 301-309.
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

  Caution: multivariate data-driven bandwidth selection methods are, by
  their nature, \emph{computationally intensive}. Virtually all methods
  require dropping the \eqn{i}th observation from the data set,
  computing an object, repeating this for all observations in the
  sample, then averaging each of these leave-one-out estimates for a
  \emph{given} value of the bandwidth vector, and only then repeating
  this a large number of times in order to conduct multivariate
  numerical minimization/maximization. Furthermore, due to the potential
  for local minima/maxima, \emph{restarting this procedure a large
  number of times may often be necessary}. This can be frustrating for
  users possessing large datasets. For exploratory purposes, you may
  wish to override the default search tolerances, say, setting
  \code{optim.reltol=.1} and conduct multistarting (the default is to
  restart min(2, ncol(xdat)) times). Once the procedure terminates, you can
  restart search with default tolerances using those bandwidths obtained
  from the less rigorous search (i.e., set \code{bws=bw} on subsequent
  calls to this routine where \code{bw} is the initial bandwidth
  object).  A version of this package using the \code{Rmpi} wrapper is
  under development that allows one to deploy this software in a
  clustered computing environment to facilitate computation involving
  large datasets.  }

\seealso{\code{\link{npRmpi.init}} for MPI startup and workflow guidance.}
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
