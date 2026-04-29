%
\name{npplregbw}
\alias{npplregbw}
\alias{npplregbw.default}
\alias{npplregbw.formula}
\alias{npplregbw.plbandwidth}

\title{Partially Linear Kernel Regression Bandwidth Selection with Mixed Data Types}

\description{
  \code{npplregbw} computes a bandwidth object for a partially linear
  kernel regression estimate of a one (1) dimensional dependent variable
  on \eqn{p+q}-variate explanatory data, using the model \eqn{Y = X\beta
  + \Theta (Z) + \epsilon}{Y = XB + theta(Z) + epsilon} given a set of
  estimation points, training points (consisting of explanatory data and
  dependent data), and a bandwidth specification, which can be a
  \code{rbandwidth} object, or a bandwidth vector, bandwidth type and
  kernel type.
}

\usage{
npplregbw(\dots)

\method{npplregbw}{formula}(formula, 
          data, 
          subset, 
          na.action, 
          call, 
          \dots)

\method{npplregbw}{default}(xdat = stop("invoked without data `xdat'"),
          ydat = stop("invoked without data `ydat'"),
          zdat = stop("invoked without data `zdat'"),
          bandwidth.compute = TRUE,
          bws,
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
          scale.factor.search.lower = NULL,
          ftol,
          itmax,
          nmulti,
          remin,
          small,
          tol,
          \dots)

\method{npplregbw}{plbandwidth}(xdat = stop("invoked without data `xdat'"),
          ydat = stop("invoked without data `ydat'"),
          zdat = stop("invoked without data `zdat'"),
          bws,
          nmulti,
          \dots)

}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the linear, nonparametric, formula, and bandwidth inputs.
}

  \item{bandwidth.compute}{
    a logical value which specifies whether to do a numerical search for
    bandwidths or not. If set to \code{FALSE}, a \code{plbandwidth} object
    will be returned with bandwidths set to those specified
    in \code{bws}. Defaults to \code{TRUE}.
  }

  \item{bws}{
    a bandwidth specification. This can be set as a \code{plbandwidth}
    object returned from an invocation of \code{npplregbw}, or as
    a matrix of bandwidths, where each row is a set of bandwidths for \eqn{Z},
    with a column for each variable \eqn{Z_i}{Z[i]}. In the first row
    are the bandwidths for the regression of \eqn{Y} on \eqn{Z}. The
    following rows contain the bandwidths for the regressions of the
    columns of \eqn{X} on \eqn{Z}.  If specified as a matrix,
    additional arguments will need to be supplied as necessary to
    specify the bandwidth type, kernel types, and so on.

    If left unspecified, \code{npplregbw} will search for optimal
    bandwidths using \code{\link{npregbw}} in the course of
    calculations. If specified, \code{npplregbw} will use the given
    bandwidths as the starting point for the numerical search for optimal
    bandwidths, unless you specify bandwidth.compute = FALSE.
  }

  \item{call}{
    the original function call. This is passed internally by
    \code{\link{npRmpi}} when a bandwidth search has been implied by a call to
    another function. It is not recommended that the user set this.
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which the
    function is called.
  }

  \item{formula}{
    a symbolic description of variables on which bandwidth selection is
    to be performed. The details of constructing a formula are
    described below.
  }

  \item{na.action}{
    a function which indicates what should happen when the data contain
   \code{NA}s. The default is set by the \code{\link{na.action}} setting of options, and is
    \code{\link{na.fail}} if that is unset. The (recommended) default is
    \code{\link{na.omit}}.
  }

  \item{subset}{
    an optional vector specifying a subset of observations to be used in
    the fitting process.
  }

  \item{xdat}{
  a \eqn{p}-variate data frame of explanatory data (training data),
  corresponding to \eqn{X} in the model equation, whose linear
  relationship with the dependent data \eqn{Y} is posited.
  }

  \item{ydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
  \code{xdat}.
  }

  \item{zdat}{
    a \eqn{q}-variate data frame of explanatory data (training data),
    corresponding to \eqn{Z} in the model equation, whose relationship
    to the dependent variable is unspecified (nonparametric)
  }


\subsection{Automatic Degree Search Controls}{
These arguments control automatic local-polynomial degree search.
}

  \item{degree.max}{
    optional scalar or integer vector giving upper bounds for automatic
    degree search over continuous \code{zdat} predictors when
    \code{degree.select != "manual"}.
  }

  \item{degree.max.cycles}{
    positive integer giving the maximum number of coordinate-search
    sweeps over the degree vector. Ignored for \code{"manual"} and
    \code{"exhaustive"} degree selection.
  }

  \item{degree.min}{
    optional scalar or integer vector giving lower bounds for automatic
    degree search over continuous \code{zdat} predictors when
    \code{degree.select != "manual"}.
  }

  \item{degree.restarts}{
    non-negative integer giving the number of additional deterministic
    coordinate-search restarts. Ignored for \code{"manual"} and
    \code{"exhaustive"} degree selection.
  }

  \item{degree.select}{
    character string controlling local-polynomial degree handling for
    the nonparametric \code{zdat} component. \code{"manual"} (default)
    treats \code{degree} as fixed. \code{"coordinate"} performs cached
    coordinate-wise search over admissible degree vectors.
    \code{"exhaustive"} evaluates the full admissible degree grid when
    \code{search.engine="cell"}. For NOMAD-based search engines, any
    non-\code{"manual"} value requests direct joint search over degree
    and bandwidth coordinates for the nonparametric component.
  }

  \item{degree.start}{
    optional starting degree vector for automatic coordinate search. If
    omitted, the search starts from the degree-zero local-constant
    baseline on the continuous \code{zdat} predictors.
  }

  \item{degree.verify}{
    logical value indicating whether a coordinate-search solution should
    be exhaustively verified over the admissible degree grid after the
    heuristic phase completes. Available only for
    \code{search.engine="cell"}.
  }


\subsection{Continuous Scale-Factor Search Controls}{
These controls define lower admissibility bounds for continuous fixed-bandwidth search.
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
These arguments control fixed local-polynomial specification for the nonparametric component.
}

  \item{degree}{
    for local-polynomial partially linear fits, polynomial degree
    specification for each continuous nonparametric regressor in
    \code{zdat}. When supplied with \code{degree.select="manual"},
    bandwidth optimization treats this vector as fixed input.
  }


\subsection{NOMAD Search Controls}{
These arguments control the optional NOMAD direct-search route for local-polynomial degree and bandwidth search.
}

  \item{nomad}{
    logical shortcut for the recommended automatic local-polynomial
    NOMAD route for the nonparametric \code{zdat} component. When
    \code{TRUE}, any missing values among \code{regtype},
    \code{search.engine}, \code{degree.select},
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
    backend for the nonparametric \code{zdat} component when
    \code{degree.select != "manual"}. \code{"nomad+powell"} (default)
    performs direct joint search over fixed bandwidths and the degree
    vector using \code{crs::snomadr()}, then applies one Powell hot
    start from the NOMAD solution. \code{"nomad"} omits the Powell
    refinement. \code{"cell"} profiles the criterion over the
    admissible degree grid using repeated fixed-degree bandwidth solves.
    NOMAD-based search currently requires fixed-bandwidth child
    templates, \code{degree.verify=FALSE}, and the suggested package
    \pkg{crs} to be installed.
  }


\subsection{Numerical Search And Tolerance Controls}{
These controls set optimizer tolerances and restart behavior.
}

  \item{ftol}{
    tolerance on the value of the cross-validation function
    evaluated at located minima. Defaults to \code{1.19e-07
      (FLT_EPSILON)}
  }

  \item{itmax}{
    integer number of iterations before failure in the numerical
    optimization routine. Defaults to \code{10000}
  }

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points. Defaults to \code{min(2,ncol(zdat))}.
  }

  \item{remin}{
    a logical value which when set as \code{TRUE} the search routine
    restarts from located minima for a minor gain in accuracy. Defaults
    to \code{TRUE}
  }

  \item{small}{ a small number, at about the precision of the data type
    used. Defaults to \code{2.22e-16 (DBL_EPSILON)}

  }

  \item{tol}{
    tolerance on the position of located minima of the
    cross-validation function. Defaults to \code{1.49e-08
      (sqrt(DBL_EPSILON))}
  }


\subsection{Additional Arguments}{
These arguments collect remaining controls passed through S3 methods.
}

  \item{\dots}{
    additional arguments supplied to specify the regression type,
    bandwidth type, kernel types, selection methods, and so on. To do
    this, you may specify any of \code{regtype},
    \code{bwmethod}, \code{bwscaling}, \code{bwtype} (one of
    \code{fixed}, \code{generalized_nn}, \code{adaptive_nn}), \code{ckertype},
    \code{ckerorder}, \code{ukertype}, \code{okertype},
    \code{bernstein.basis}, and \code{basis}, as described in
    \code{\link{npregbw}}.
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

  Documentation guide: see \code{\link{npregbw}} for component
  nonparametric regression bandwidth controls, \code{\link{np.kernels}}
  for kernels, \code{\link{np.options}} for global options,
  \code{\link{plot}} for plotting options, and
  \code{\link{npRmpi.init}} for interactive/cluster MPI startup. See
  \code{\link{npRmpi.init}} details for performance tradeoffs (message
  passing/startup mode) and the \code{inst/Rprofile} manual-broadcast
  template.

  The partially linear bandwidth-selection argument surface is easiest
  to read by decision group: linear \code{xdat} inputs,
  nonparametric \code{zdat} inputs, and existing bandwidth inputs;
  local-polynomial/NOMAD controls for the nonparametric component;
  numerical search and feasibility controls; formula-interface
  controls; and additional bandwidth, kernel, and support controls that
  are passed to the component \code{\link{npregbw}} searches.

  For S3 plotting help, use \code{methods("plot")} and query
  class-specific help topics such as \code{?plot.npregression} and
  \code{?plot.rbandwidth}. You can inspect implementations with
  \code{getS3method("plot","npregression")}.

  \code{npplregbw} implements a variety of methods for nonparametric
 regression on multivariate (\eqn{q}-variate) explanatory data defined
 over a set of possibly continuous and/or discrete (unordered, ordered)
 data. The approach is based on Li and Racine (2003), who employ
 \sQuote{generalized product kernels} that admit a mix of continuous and
 discrete data types.

 Three classes of kernel estimators for the continuous data types are
 available: fixed, adaptive nearest-neighbor, and generalized
 nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
 each sample realization in the set, \eqn{x_i}{x[i]}, when estimating the
 density at the point \eqn{x}. Generalized nearest-neighbor bandwidths change
 with the point at which the density is estimated, \eqn{x}. Fixed bandwidths
 are constant over the support of \eqn{x}.

 \code{npplregbw} may be invoked \emph{either} with a formula-like
 symbolic
 description of variables on which bandwidth selection is to be
 performed \emph{or} through a simpler interface whereby data is passed
 directly to the function via the \code{xdat}, \code{ydat}, and
 \code{zdat}
 parameters. Use of these two interfaces is \bold{mutually exclusive}.

 Data contained in the data frame \code{zdat} may be a mix of continuous
 (default), unordered discrete (to be specified in the data frame
 \code{zdat} using \code{\link{factor}}), and ordered discrete (to be
 specified in the data frame \code{zdat} using
 \code{\link{ordered}}). Data can be entered in an arbitrary order and
 data types will be detected automatically by the routine (see
 \code{\link{npRmpi}} for details).

 Data for which bandwidths are to be estimated may be specified
 symbolically. A typical description has the form \code{dependent}
 \code{data} \code{~} \code{parametric} \code{explanatory} \code{data}
 \code{|} \code{nonparametric} \code{explanatory} \code{data},
 where \code{dependent} \code{data} is a univariate response, and
 \code{parametric} \code{explanatory} \code{data} and
 \code{nonparametric} \code{explanatory}
 \code{data} are both series of variables specified by name, separated by
 the separation character '+'. For example, \code{ y1 ~ x1 + x2 | z1 }
 specifies that the bandwidth object for the partially linear model with
 response \code{y1}, linear parametric regressors \code{x1} and
 \code{x2}, and
 nonparametric regressor \code{z1} is to be estimated. See below for
 further examples.

 A variety of kernels may be specified by the user. Kernels implemented
 for continuous data types include the second, fourth, sixth, and eighth
 order Gaussian and Epanechnikov kernels, and the uniform
 kernel. Unordered discrete data types use a variation on Aitchison and
 Aitken's (1976) kernel, while ordered data types use a variation of the
 Wang and van Ryzin (1981) kernel.

 When the nonparametric component is estimated with
 \code{regtype="lp"} and \code{degree.select != "manual"},
 \code{npplregbw} can jointly determine the \code{zdat}-side degree
 vector and the associated bandwidth coordinates. With
 \code{search.engine="cell"}, the criterion is profiled over the degree
 grid using cached coordinate-wise or exhaustive search together with
 repeated fixed-degree bandwidth solves. With
 \code{search.engine="nomad"} or \code{"nomad+powell"}, the criterion
 is optimized directly over the joint degree/bandwidth space using
 \code{crs::snomadr()}; \code{"nomad+powell"} then performs one Powell
 hot start and keeps the better of the direct NOMAD and polished
 solutions. For the nonparametric regression component, this
 polynomial-adaptive joint-search route follows Hall and Racine (2015).

 Setting \code{nomad=TRUE} is a convenience preset for this automatic
 LP route, not a generic optimizer alias. For partially linear
 regression it expands any missing values to the equivalent long-form
 call
 \preformatted{
   npplregbw(...,
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
  if \code{bwtype} is set to \code{fixed}, an object containing bandwidths
  (or scale factors if \code{bwscaling = TRUE}) is returned. If it is set to
  \code{generalized_nn} or \code{adaptive_nn}, then instead the \eqn{k}th nearest
  neighbors are returned for the continuous variables while the discrete
  kernel bandwidths are returned for the discrete variables. Bandwidths
  are stored in a list under the component name \code{bw}. Each element
  is an \code{rbandwidth} object. The first
  element of the list corresponds to the regression of \eqn{Y} on \eqn{Z}.
  Each subsequent element is the bandwidth object corresponding to the
  regression of the \eqn{i}th column of \eqn{X} on \eqn{Z}. See examples
  for more information.
}

\references{

 Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
 discrimination by the kernel method,} Biometrika, 63, 413-420.

 Gao, Q. and L. Liu and J.S. Racine (2015), \dQuote{A partially linear
 kernel estimator for categorical data,} Econometric Reviews, 34 (6-10),
958-977.

 Hall, P. and J.S. Racine (2015), \dQuote{Infinite Order Cross-Validated
 Local Polynomial Regression,} Journal of Econometrics, 185, 510-525.

 Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
 and Practice,} Princeton University Press.

 Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated local linear
   nonparametric regression,} Statistica Sinica, 14, 485-512.

 Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
 Cambridge University Press.

 Racine, J.S. and Q. Li (2004), \dQuote{Nonparametric estimation of regression
 functions with both categorical and continuous data,} Journal of
 Econometrics, 119, 99-130.

 Robinson, P.M. (1988), \dQuote{Root-n-consistent semiparametric regression,}
 Econometrica, 56, 931-954.

 Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth estimators
 for discrete distributions,}  Biometrika, 68, 301-309.

}
\author{
    Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.

  Caution: multivariate data-driven bandwidth selection methods are, by
  their nature, \emph{computationally intensive}. Virtually all methods
  require dropping the \eqn{i}th observation from the data set, computing an
  object, repeating this for all observations in the sample, then
  averaging each of these leave-one-out estimates for a \emph{given}
  value of the bandwidth vector, and only then repeating this a large
  number of times in order to conduct multivariate numerical
  minimization/maximization. Furthermore, due to the potential for local
  minima/maxima, \emph{restarting this procedure a large number of times may
  often be necessary}. This can be frustrating for users possessing
  large datasets. For exploratory purposes, you may wish to override the
  default search tolerances, say, setting ftol=.01 and tol=.01 and
  conduct multistarting (the default is to restart min(2, ncol(zdat))
  times) as is done for a number of examples. Once the procedure
  terminates, you can restart search with default tolerances using those
  bandwidths obtained from the less rigorous search (i.e., set
  \code{bws=bw} on subsequent calls to this routine where \code{bw} is
  the initial bandwidth object).  A version of this package using the
  \code{Rmpi} wrapper is under development that allows one to deploy
  this software in a clustered computing environment to facilitate
  computation involving large datasets.
}

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}
  \code{\link{npregbw}}, \code{\link{npreg}} }
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

n <- 250
x1 <- rnorm(n)
x2 <- rbinom(n, 1, .5)

z1 <- rbinom(n, 1, .5)
z2 <- rnorm(n)

y <- 1 + x1 + x2 + z1 + sin(z2) + rnorm(n)

x2 <- factor(x2)
z1 <- factor(z1)

bw <- npplregbw(formula=y~x1+x2|z1+z2)

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
