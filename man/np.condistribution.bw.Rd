
\name{npcdistbw}
\alias{npcdistbw}
\alias{npcdistbw.condbandwidth}
\alias{npcdistbw.default}
\alias{npcdistbw.formula}
\title{Kernel Conditional Distribution Bandwidth Selection with Mixed
Data Types}

\description{
  \code{npcdistbw} computes a \code{condbandwidth} object for estimating
  a \eqn{p+q}-variate kernel conditional cumulative distribution
  estimator defined over mixed continuous and discrete (unordered
  \code{xdat}, ordered \code{xdat} and \code{ydat}) data using either
  the normal-reference rule-of-thumb or least-squares cross validation
  method of Li and Racine (2008) and Li, Lin and Racine
  (2013).
}

\usage{
npcdistbw(\dots)

\method{npcdistbw}{formula}(formula, 
          data, 
          subset, 
          na.action, 
          call, 
          gdata = NULL,
          \dots)


\method{npcdistbw}{condbandwidth}(xdat = stop("data 'xdat' missing"),
          ydat = stop("data 'ydat' missing"),
          gydat = NULL,
          bws,
          bandwidth.compute = TRUE,
          cfac.dir = 2.5*(3.0-sqrt(5)),
          scale.factor.init = 0.5,
          dfac.dir = 0.25*(3.0-sqrt(5)),
          dfac.init = 0.375,
          dfc.dir = 3,
          do.full.integral = FALSE,
          ftol = 1.490116e-07,
          scale.factor.init.upper = 2.0,
          hbd.dir = 1,
          hbd.init = 0.9,
          initc.dir = 1.0,
          initd.dir = 1.0,
          invalid.penalty = c("baseline","dbmax"),
          itmax = 10000,
          lbc.dir = 0.5,
          scale.factor.init.lower = 0.1,
          lbd.dir = 0.1,
          lbd.init = 0.1,
          memfac = 500.0,
          ngrid = 100,
          nmulti,
          penalty.multiplier = 10,
          remin = TRUE,
          scale.init.categorical.sample = FALSE,
          scale.factor.search.lower = NULL,
          small = 1.490116e-05,
          tol = 1.490116e-04,
          transform.bounds = FALSE,
          \dots)

\method{npcdistbw}{default}(xdat = stop("data 'xdat' missing"),
          ydat = stop("data 'ydat' missing"),
          gydat,
          bws,
          bandwidth.compute = TRUE,
          bwmethod,
          bwscaling,
          bwtype,
          cfac.dir,
          scale.factor.init,
          cxkerbound,
          cxkerlb,
          cxkerorder,
          cxkertype,
          cxkerub,
          cykerbound,
          cykerlb,
          cykerorder,
          cykertype,
          cykerub,
          dfac.dir,
          dfac.init,
          dfc.dir,
          do.full.integral,
          ftol,
          scale.factor.init.upper,
          hbd.dir,
          hbd.init,
          initc.dir,
          initd.dir,
          invalid.penalty,
          itmax,
          lbc.dir,
          scale.factor.init.lower,
          lbd.dir,
          lbd.init,
          memfac,
          ngrid,
          nmulti,
          oxkertype,
          oykertype,
          penalty.multiplier,
          remin,
          scale.init.categorical.sample,
          scale.factor.search.lower = NULL,
          small,
          tol,
          transform.bounds,
          uxkertype,
          regtype = c("lc", "ll", "lp"),
          basis = c("glp", "additive", "tensor"),
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
          bernstein.basis = FALSE,
          \dots)
}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the data, formula interface, optional distribution grid, and whether bandwidths are supplied or computed.
}

  \item{bandwidth.compute}{
    a logical value which specifies whether to do a numerical search for
    bandwidths or not. If set to \code{FALSE}, a \code{condbandwidth} object
    will be returned with bandwidths set to those specified
    in \code{bws}. Defaults to \code{TRUE}.
  }

  \item{bws}{
    a bandwidth specification. This can be set as a \code{condbandwidth}
    object returned from a previous invocation, or as a \eqn{p+q}-vector
    of bandwidths, with each element \eqn{i} up to \eqn{i=q}
    corresponding to the bandwidth for column \eqn{i} in \code{ydat},
    and each element \eqn{i} from \eqn{i=q+1} to \eqn{i=p+q}
    corresponding to the bandwidth for column \eqn{i-q} in
    \code{xdat}. In either case, the bandwidth supplied will serve as a
    starting point in the numerical search for optimal bandwidths. If
    specified as a vector, then additional arguments will need to be
    supplied as necessary to specify the bandwidth type, kernel types,
    selection methods, and so on. This can be left unset.
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

  \item{gdata}{
    a grid of data on which the indicator function for
    least-squares cross-validation is to be computed (can be the sample
    or a grid of quantiles).
  }

  \item{gydat}{
    a grid of data on which the indicator function for
    least-squares cross-validation is to be computed (can be the sample
    or a grid of quantiles for \code{ydat}).
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
    a \eqn{p}-variate data frame of explanatory data on which bandwidth selection will
    be performed. The data types may be continuous, discrete (unordered
    and ordered factors), or some combination thereof.
  }

  \item{ydat}{
    a \eqn{q}-variate data frame of dependent data on which bandwidth
    selection will be performed. The data types may be continuous,
    discrete (ordered factors), or some combination thereof.
  }


\subsection{Automatic Degree Search Controls}{
These arguments control automatic local-polynomial degree search when \code{regtype="lp"}.
}

  \item{degree.max}{
    optional scalar or integer vector giving upper bounds for automatic
    degree search over continuous \code{xdat} predictors when
    \code{degree.select != "manual"}.
  }

  \item{degree.max.cycles}{
    positive integer giving the maximum number of coordinate-search
    sweeps over the degree vector. Ignored for \code{"manual"} and
    \code{"exhaustive"} degree selection.
  }

  \item{degree.min}{
    optional scalar or integer vector giving lower bounds for automatic
    degree search over continuous \code{xdat} predictors when
    \code{degree.select != "manual"}.
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
    coordinate-wise search over admissible degree vectors for the
    continuous \code{xdat} predictors. \code{"exhaustive"} evaluates
    the full admissible degree grid when \code{search.engine="cell"}.
    For NOMAD-based search engines, any non-\code{"manual"} value
    requests direct joint search over degree and bandwidth
    coordinates.
  }

  \item{degree.start}{
    optional starting degree vector for automatic coordinate search. If
    omitted, the search starts from the degree-zero local-constant
    baseline on the continuous \code{xdat} predictors.
  }

  \item{degree.verify}{
    logical value indicating whether a coordinate-search solution should
    be exhaustively verified over the admissible degree grid after the
    heuristic phase completes. Available only for
    \code{search.engine="cell"}.
  }


\subsection{Bandwidth Criterion And Representation}{
These arguments choose the selection criterion and the way continuous bandwidths are represented.
}

  \item{bwmethod}{ which method to use to select bandwidths.
    \code{cv.ls} specifies least-squares cross-validation (Li, Lin and
    Racine (2013), and \code{normal-reference} just computes the
    \sQuote{rule-of-thumb} bandwidth \eqn{h_j}{h[j]} using the standard
    formula \eqn{h_j = 1.06 \sigma_j n^{-1/(2P+l)}}{h[j] =
    1.06*sigma[j]*n^(-1.0/(2.0*P+l))}, where \eqn{\sigma_j}{sigma[j]} is
    an adaptive measure of spread of the \eqn{j}th continuous variable
    defined as min(standard deviation, mean absolute deviation/1.4826,
    interquartile range/1.349), \eqn{n} the number of observations,
    \eqn{P} the order of the kernel, and \eqn{l} the number of
    continuous variables. Note that when there exist factors and the
    normal-reference rule is used, there is zero smoothing of the
    factors. Defaults to \code{cv.ls}.  }

  \item{bwscaling}{ a logical value that when set to \code{TRUE} the
    supplied bandwidths are interpreted as \sQuote{scale factors}
    (\eqn{c_j}{c[j]}), otherwise when the value is \code{FALSE} they are
    interpreted as \sQuote{raw bandwidths} (\eqn{h_j}{h[j]} for continuous data
    types, \eqn{\lambda_j}{lambda[j]} for discrete data types). For
    continuous data types, \eqn{c_j}{c[j]} and \eqn{h_j}{h[j]} are
    related by the formula \eqn{h_j = c_j \sigma_j n^{-1/(2P+l)}}{h[j] =
    c[j]*sigma[j]*n^(-1/(2*P+l))}, where \eqn{\sigma_j}{sigma[j]} is an
    adaptive measure of spread of continuous variable \eqn{j} defined as
    min(standard deviation, mean absolute deviation/1.4826, interquartile
    range/1.349), \eqn{n} the number of observations, \eqn{P}{P} the
    order of the kernel, and \eqn{l}{l} the number of continuous
    variables. For discrete data types, \eqn{c_j}{c[j]} and
    \eqn{h_j}{h[j]} are related by the formula \eqn{h_j =
    c_jn^{-2/(2P+l)}}{h[j] = c[j]*n^(-2/(2*P+l))}, where here
    \eqn{j}{[j]} denotes discrete variable \eqn{j}{j}.  Defaults to
    \code{FALSE}.  }

  \item{bwtype}{
    character string used for the continuous variable bandwidth type,
    specifying the type of bandwidth to compute and return in the
    \code{condbandwidth} object. Defaults to \code{fixed}. Option
    summary:\cr
    \code{fixed}: compute fixed bandwidths \cr
    \code{generalized_nn}: compute generalized nearest neighbors \cr
    \code{adaptive_nn}: compute adaptive nearest neighbors
  }


\subsection{Categorical Search Initialization}{
These controls set categorical search starts and categorical direction-set initialization.
}

  \item{dfac.dir}{stretch factor for direction set search for Powell's algorithm for categorical variables. See Details}

  \item{dfac.init}{non-random initial values for scale factors for categorical variables for Powell's algorithm. See Details}

  \item{hbd.dir}{upper bound for direction set search for Powell's algorithm for categorical variables. See Details}

  \item{hbd.init}{upper bound for scale factors for categorical variables for Powell's algorithm. See Details}

  \item{initd.dir}{initial non-random values for direction set search for Powell's algorithm for categorical variables. See Details}

  \item{lbd.dir}{lower bound for direction set search for Powell's algorithm for categorical variables. See Details}

  \item{lbd.init}{lower bound for scale factors for categorical variables for Powell's algorithm. See Details}

  \item{scale.init.categorical.sample}{ a logical value that when set
    to \code{TRUE} scales \code{lbd.dir}, \code{hbd.dir},
    \code{dfac.dir}, and \code{initd.dir} by \eqn{n^{-2/(2P+l)}},
    \eqn{n} the number of observations, \eqn{P}{P} the order of the
    kernel, and \eqn{l}{l} the number of \code{numeric} variables. See
    Details}


\subsection{Continuous Direction-Set Search Controls}{
These controls set Powell direction-set initialization for continuous variables.
}

  \item{cfac.dir}{stretch factor for direction set search for Powell's algorithm for \code{numeric} variables. See Details}

  \item{dfc.dir}{chi-square degrees of freedom for direction set search for Powell's algorithm for \code{numeric} variables. See Details}

  \item{initc.dir}{initial non-random values for direction set search for Powell's algorithm for \code{numeric} variables. See Details}

  \item{lbc.dir}{lower bound for direction set search for Powell's algorithm for \code{numeric} variables. See Details}


\subsection{Continuous Kernel Support Controls}{
These controls choose and parameterize bounded support for continuous kernels.
}

  \item{cxkerbound}{
    character string controlling continuous-kernel support handling for
    \code{xdat}. Can be set as \code{none} (default kernel on full support),
    \code{range} (use sample min/max), or \code{fixed} (use
    \code{cxkerlb}/\code{cxkerub}). The bounded-kernel route reuses the selected continuous kernel and renormalizes it on the chosen support; see \code{\link{np.kernels}}.
  }

  \item{cxkerlb}{
    numeric scalar/vector of lower bounds for continuous \code{xdat}
    variables used when \code{cxkerbound="fixed"}. Must satisfy
    lower-bound validity for each variable (e.g.,
    \code{<= min(variable)}). Use \code{-Inf} for unbounded below. See \code{\link{np.kernels}} for bounded-kernel normalization details.
  }

  \item{cxkerub}{
    numeric scalar/vector of upper bounds for continuous \code{xdat}
    variables used when \code{cxkerbound="fixed"}. Must satisfy
    upper-bound validity for each variable (e.g.,
    \code{>= max(variable)}). Use \code{Inf} for unbounded above. See \code{\link{np.kernels}} for bounded-kernel normalization details.
  }

  \item{cykerbound}{
    character string controlling continuous-kernel support handling for
    \code{ydat}. Can be set as \code{none} (default kernel on full support),
    \code{range} (use sample min/max), or \code{fixed} (use
    \code{cykerlb}/\code{cykerub}). The bounded-kernel route reuses the selected continuous kernel and renormalizes it on the chosen support; see \code{\link{np.kernels}}.
  }

  \item{cykerlb}{
    numeric scalar/vector of lower bounds for continuous \code{ydat}
    variables used when \code{cykerbound="fixed"}. Must satisfy
    lower-bound validity for each variable (e.g.,
    \code{<= min(variable)}). Use \code{-Inf} for unbounded below. See \code{\link{np.kernels}} for bounded-kernel normalization details.
  }

  \item{cykerub}{
    numeric scalar/vector of upper bounds for continuous \code{ydat}
    variables used when \code{cykerbound="fixed"}. Must satisfy
    upper-bound validity for each variable (e.g.,
    \code{>= max(variable)}). Use \code{Inf} for unbounded above. See \code{\link{np.kernels}} for bounded-kernel normalization details.
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

\subsection{Distribution Integral And Grid Controls}{
These controls tune the conditional distribution-function integral and grid calculations.
}

  \item{do.full.integral}{
    a logical value which when set as \code{TRUE} evaluates the
    moment-based integral on the entire sample.
  }

  \item{memfac}{
    The algorithm to compute the least-squares objective function uses
    a block-based algorithm to eliminate or minimize redundant kernel
    evaluations. Due to memory, hardware and software constraints, a
    maximum block size must be imposed by the algorithm. This block size
    is roughly equal to memfac*10^5 elements. Empirical tests on
    modern hardware find that a memfac of around 500 performs well. If
    you experience out of memory errors, or strange behaviour for
    large data sets (>100k elements) setting memfac to a lower value may
    fix the problem.
  }

  \item{ngrid}{
    integer number of grid points to use when computing the moment-based
    integral. Defaults to \code{100}.
  }


\subsection{Kernel Type Controls}{
These controls choose continuous, unordered, and ordered kernels for \code{xdat} and \code{ydat}.
}

  \item{cxkerorder}{
    numeric value specifying kernel order for
    \code{xdat} (one of
    \code{(2,4,6,8)}). Kernel order specified along with a
  \code{uniform} continuous kernel type will be ignored. Defaults to
  \code{2}.
  }

  \item{cxkertype}{
    character string used to specify the continuous kernel type for
    \code{xdat}.  Can be set as \code{gaussian},
    \code{epanechnikov}, or \code{uniform}. Defaults to \code{gaussian}.
  }

  \item{cykerorder}{
    numeric value specifying kernel order for
    \code{ydat} (one of
    \code{(2,4,6,8)}). Kernel order specified along with a
  \code{uniform} continuous kernel type will be ignored. Defaults to
  \code{2}.
  }

  \item{cykertype}{
    character string used to specify the continuous kernel type for
    \code{ydat}.
    Can be set as \code{gaussian}, \code{epanechnikov}, or
    \code{uniform}. Defaults to \code{gaussian}.
  }

  \item{oxkertype}{
    character string used to specify the ordered categorical kernel
    type for \code{xdat}. Can be set as \code{wangvanryzin},
    \code{liracine}, or \code{racineliyan}. Defaults to
    \code{liracine}.
  }
  \item{oykertype}{
    character string used to specify the ordered categorical kernel
    type for \code{ydat}. Can be set as \code{wangvanryzin},
    \code{liracine}, or \code{racineliyan}. Defaults to
    \code{liracine}.
  }
  \item{uxkertype}{
    character string used to specify the unordered categorical kernel
    type for \code{xdat}. Can be set as \code{aitchisonaitken} or
    \code{liracine}. Defaults to \code{aitchisonaitken}.
  }

\subsection{Local-Polynomial Model Specification}{
These arguments control the local-polynomial estimator, basis, and fixed degree specification.
}

  \item{basis}{
    character string specifying the polynomial basis used when
    \code{regtype="lp"}. Options are \code{"glp"},
    \code{"additive"}, and \code{"tensor"}.
  }

  \item{bernstein.basis}{
    logical value controlling Bernstein basis evaluation for
    \code{regtype="lp"}. When automatic degree search is requested and
    \code{bernstein.basis} is not explicitly supplied, the search route
    defaults to \code{TRUE} for numerical stability. Explicit
    \code{bernstein.basis=FALSE} is honored, but raw-polynomial search
    can be poorly conditioned at higher degrees.
  }

  \item{degree}{
    integer scalar or integer vector of polynomial degrees for
    continuous \code{xdat} variables when \code{regtype="lp"}.
    If scalar, the value is recycled to all continuous \code{xdat}
    variables.
  }

  \item{regtype}{
    character string specifying the conditional local method used for
    the \code{xdat} regression weight operator. Options are
    \code{"lc"}, \code{"ll"}, and \code{"lp"}.
    For \code{npc*} methods, \code{"ll"} is implemented via the
    canonical local polynomial engine with \code{degree = 1} and
    \code{basis = "glp"}. If local-linear \code{cv.ls} search fails
    while using this canonical raw basis, retrying explicitly with
    \code{regtype="lp"}, \code{degree=1}, and
    \code{bernstein.basis=TRUE}, or centering/scaling the continuous
    regressors, can improve numerical conditioning without changing package
    defaults or invoking an automatic fallback.
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
    performs direct joint search over the \code{xdat}-side continuous
    bandwidth coordinates and degree vector using
    \code{crs::snomadr()}, then applies one Powell hot start from the
    NOMAD solution. \code{"nomad"} omits the Powell refinement.
    \code{"cell"} profiles the criterion over the admissible degree
    grid using repeated fixed-degree bandwidth solves. NOMAD-based
    search currently requires \code{bwtype="fixed"},
    \code{degree.verify=FALSE}, and the suggested package \pkg{crs} to
    be installed.
  }


\subsection{Numerical Search And Tolerance Controls}{
These controls set optimizer tolerances, restart behavior, invalid-candidate penalties, memory blocking, and bounded search transformations.
}

 \item{ftol}{
    fractional tolerance on the value of the cross-validation function
    evaluated at located minima (of order the machine precision or
    perhaps slightly larger so as not to be diddled by
    roundoff). Defaults to \code{1.490116e-07}
    (1.0e+01*sqrt(.Machine$double.eps)).
  }

  \item{invalid.penalty}{ a character string specifying the penalty
    used when the optimizer encounters invalid bandwidths.
    \code{"baseline"} returns a finite penalty based on a baseline
    objective; \code{"dbmax"} returns \code{DBL\_MAX}. Defaults to
    \code{"baseline"}.}

  \item{itmax}{
    integer number of iterations before failure in the numerical
    optimization routine. Defaults to \code{10000}.
  }

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial points
  }

  \item{penalty.multiplier}{ a numeric multiplier applied to the
    baseline penalty when \code{invalid.penalty="baseline"}.
    Defaults to \code{10}.}

  \item{remin}{
    a logical value which when set as \code{TRUE} the search routine
    restarts from located minima for a minor gain in accuracy. Defaults
    to \code{TRUE}.
  }

  \item{small}{
    a small number used to bracket a minimum (it is hopeless to ask for
    a bracketing interval of width less than sqrt(epsilon) times its
    central value, a fractional width of only about 10-04 (single
    precision) or 3x10-8 (double precision)). Defaults to \code{small
      = 1.490116e-05 (1.0e+03*sqrt(.Machine$double.eps))}.
  }

  \item{tol}{
    tolerance on the position of located minima of the cross-validation
    function (tol should generally be no smaller than the square root of
    your machine's floating point precision). Defaults to \code{
      1.490116e-04 (1.0e+04*sqrt(.Machine$double.eps))}.
  }

  \item{transform.bounds}{ a logical value that when set to \code{TRUE}
    applies an internal transformation that maps the unconstrained
    search to the feasible bandwidth domain. Defaults to
    \code{FALSE}.}


\subsection{Additional Arguments}{
These arguments collect remaining controls passed through S3 methods.
}

  \item{\dots}{
    additional arguments supplied to specify the bandwidth type,
    kernel types, selection methods, and so on, detailed below.
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
  start. \code{scale.factor.init.lower} and
  \code{scale.factor.init.upper} define the random multistart interval.
  \code{scale.factor.search.lower} is the lower admissibility bound for
  continuous fixed-bandwidth search candidates. The effective first
  start is \code{max(scale.factor.init, scale.factor.search.lower)},
  and the effective random-start lower endpoint is
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

  The bandwidth-selection argument surface is easiest to read by
  decision group. Start by initializing MPI execution if needed, then
  choose the data and bandwidth inputs (\code{xdat}, \code{ydat},
  \code{gydat}, \code{bws}, and \code{bandwidth.compute}), then choose
  the bandwidth criterion and representation (\code{bwmethod},
  \code{bwscaling}, and \code{bwtype}). Next choose continuous kernel
  and support controls (\code{cxker*} and \code{cyker*}), categorical
  kernel controls (\code{uxkertype}, \code{oxkertype}, and
  \code{oykertype}), and numerical search controls including
  \code{nmulti}, tolerances, penalties, and the \code{scale.factor.*}
  search-start and admissibility controls. Local-polynomial and NOMAD
  controls (\code{regtype}, \code{basis}, \code{degree*},
  \code{search.engine}, \code{nomad}, \code{nomad.nmulti}, and
  \code{bernstein.basis}) are relevant when using the explicit
  local-polynomial route.

  For S3 plotting help, use \code{methods("plot")} and query
  class-specific help topics such as \code{?plot.npregression} and
  \code{?plot.rbandwidth}. You can inspect implementations with
  \code{getS3method("plot","npregression")}.

  \code{npcdistbw} implements a variety of methods for choosing
  bandwidths for multivariate distributions (\eqn{p+q}-variate) defined
  over a set of possibly continuous and/or discrete (unordered
  \code{xdat}, ordered \code{xdat} and \code{ydat}) data. The approach
  is based on Li and Racine (2004) who employ \sQuote{generalized
  product kernels} that admit a mix of continuous and discrete data
  types.

  The cross-validation methods employ multivariate numerical search
  algorithms. For fixed local-constant/local-linear fits, and for
  local-polynomial fits with \code{degree.select="manual"}, bandwidth
  search uses multidimensional Powell direction-set optimization.

  Bandwidths can (and will) differ for each variable which is, of
  course, desirable.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating
  the cumulative distribution at the point \eqn{x}. Generalized nearest-neighbor
  bandwidths change with the point at which the cumulative distribution is estimated,
  \eqn{x}. Fixed bandwidths are constant over the support of \eqn{x}.

  \code{npcdistbw} may be invoked \emph{either} with a formula-like
  symbolic
  description of variables on which bandwidth selection is to be
  performed \emph{or} through a simpler interface whereby data is passed
  directly to the function via the \code{xdat} and \code{ydat}
  parameters. Use of these two interfaces is \bold{mutually exclusive}.

  Data contained in the data frame \code{xdat} may be a mix of
  continuous (default), unordered discrete (to be specified in the data
  frames using \code{\link{factor}}), and ordered discrete (to be
  specified in the data frames using \code{\link{ordered}}).  Data
  contained in the data frame \code{ydat} may be a mix of continuous
  (default) and ordered discrete (to be specified in the data frames
  using \code{\link{ordered}}). Data can be entered in an arbitrary
  order and data types will be detected automatically by the routine
  (see \code{\link{npRmpi}} for details).

  Data for which bandwidths are to be estimated may be specified
  symbolically. A typical description has the form \code{dependent data
  ~ explanatory data},
  where \code{dependent data} and \code{explanatory data} are both
  series of variables specified by name, separated by
  the separation character '+'. For example, \code{ y1 + y2 ~ x1 + x2 }
  specifies that the bandwidths for the joint distribution of variables
  \code{y1} and \code{y2} conditioned on \code{x1} and \code{x2} are to
  be estimated. See below for further examples.

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and eighth
  order Gaussian and Epanechnikov kernels, and the uniform
  kernel. Unordered discrete data types use a variation on Aitchison and
  Aitken's (1976) kernel, while ordered data types use a variation of the
  Wang and van Ryzin (1981) kernel.

  When \code{regtype="lp"} and \code{degree.select != "manual"},
  \code{npcdistbw} can jointly determine the \code{xdat}-side local
  polynomial degree vector and the fixed bandwidth coordinates entering
  the conditional distribution criterion. With
  \code{search.engine="cell"}, the criterion is profiled over the
  admissible degree grid using cached coordinate-wise or exhaustive
  search. With \code{search.engine="nomad"} or
  \code{"nomad+powell"}, the criterion is optimized directly over the
  joint degree/bandwidth space using \code{crs::snomadr()};
  \code{"nomad+powell"} then performs one Powell hot start from the
  NOMAD solution and keeps the better of the direct NOMAD and polished
  answers. This polynomial-adaptive joint-search route is motivated by
  Hall and Racine (2015) together with Li, Li, and Racine (under
  revision). When \code{bernstein.basis} is not explicitly supplied,
  the automatic search route defaults to \code{bernstein.basis=TRUE}
  for numerical stability.

  Setting \code{nomad=TRUE} is a convenience preset for this automatic
  LP route, not a generic optimizer alias. For conditional distribution
  bandwidth selection it expands any missing values to the equivalent
  long-form call
  \preformatted{
    npcdistbw(...,
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

  The optimizer invoked for search is Powell's conjugate direction
  method which requires the setting of (non-random) initial values and
  search directions for bandwidths, and, when restarting, random values
  for successive invocations. Bandwidths for \code{numeric} variables
  are scaled by robust measures of spread, the sample size, and the
  number of \code{numeric} variables where appropriate. Two sets of
  parameters for bandwidths for \code{numeric} can be modified, those
  for initial values for the parameters themselves, and those for the
  directions taken (Powell's algorithm does not involve explicit
  computation of the function's gradient). The default values are set by
  considering search performance for a variety of difficult test cases
  and simulated cases. We highly recommend restarting search a large
  number of times to avoid the presence of local minima (achieved by
  modifying \code{nmulti}). Further refinement for difficult cases can
  be achieved by modifying these sets of parameters. However, these
  parameters are intended more for the authors of the package to enable
  \sQuote{tuning} for various methods rather than for the user themselves.

}
\value{

  \code{npcdistbw} returns a \code{condbandwidth} object, with the
  following components:

  \item{xbw}{ bandwidth(s), scale factor(s) or nearest neighbours for the
    explanatory data, \code{xdat} }
  \item{ybw}{ bandwidth(s), scale factor(s) or nearest neighbours for the
    dependent data, \code{ydat} }
  \item{fval}{ objective function value at minimum }

  if \code{bwtype} is set to \code{fixed}, an object containing
  bandwidths (or scale factors if \code{bwscaling = TRUE}) is
  returned. If it is set to \code{generalized_nn} or \code{adaptive_nn},
  then instead the \eqn{k}th nearest neighbors are returned for the
  continuous variables while the discrete kernel bandwidths are returned
  for the discrete variables.

  The functions \code{\link{predict}}, \code{\link{summary}} and \code{\link{plot}} support
  objects of type \code{condbandwidth}.
}

\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
  discrimination by the kernel method,} Biometrika, 63, 413-420.

  Hall, P. and J.S. Racine and Q. Li (2004), \dQuote{Cross-validation and
  the estimation of conditional probability densities,} Journal of the
  American Statistical Association, 99, 1015-1026.

  Hall, P. and J.S. Racine (2015), \dQuote{Infinite Order Cross-Validated
  Local Polynomial Regression,} Journal of Econometrics, 185, 510-525.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2008), \dQuote{Nonparametric estimation of
  conditional CDF and quantile functions with mixed categorical and
  continuous data,} Journal of Business and Economic Statistics, 26,
  423-434.

  Li, Q. and J. Lin and J.S. Racine (2013), \dQuote{Optimal bandwidth
  selection for nonparametric conditional distribution and quantile
  functions}, Journal of Business and Economic Statistics, 31, 57-65.

  Li, A. and Q. Li and J.S. Racine (under revision), \dQuote{Boundary
  Adjusted, Polynomial Adaptive, Nonparametric Kernel Conditional
  Density Estimation,} Econometric Reviews.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.

  Scott, D.W. (1992), \emph{Multivariate Density Estimation. Theory,
  Practice and Visualization,} New York: Wiley.

  Silverman, B.W. (1986), \emph{Density Estimation,} London: Chapman and
  Hall.

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
  conduct multistarting (the default is to restart min(2, ncol(xdat,ydat))
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
  \code{\link{bw.nrd}}, \code{\link{bw.SJ}}, \code{\link{hist}},
  \code{\link{npudens}}, \code{\link{npudist}}
}

\examples{
\dontrun{
## Not run in checks: data-driven conditional CDF bandwidth selection is
## computationally intensive and may exceed check limits under MPI.
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
npRmpi.init(nslaves=1)

data("Italy")

bw <- npcdistbw(formula=gdp~ordered(year),
                data=Italy)

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
} % enddontrun
}
\keyword{nonparametric}
