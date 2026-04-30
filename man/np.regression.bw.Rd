% 
\name{npregbw}
\alias{npregbw}
\alias{npregbw.default}
\alias{npregbw.formula}
\alias{npregbw.rbandwidth}

\title{Kernel Regression Bandwidth Selection with Mixed Data Types}

\description{
  \code{npregbw} computes a bandwidth object for a
  \eqn{p}-variate kernel regression estimator defined over mixed
  continuous and discrete (unordered, ordered) data using expected
  Kullback-Leibler cross-validation, or least-squares cross validation
  using the method of Racine and Li (2004) and Li and Racine (2004).
}

\usage{
npregbw(\dots)

\method{npregbw}{formula}(formula, 
        data, 
        subset, 
        na.action, 
        call, 
        \dots)

\method{npregbw}{default}(xdat = stop("invoked without data 'xdat'"),
        ydat = stop("invoked without data 'ydat'"),
        bws,
        bandwidth.compute = TRUE,
        basis,
        bernstein.basis,
        bwmethod,
        bwscaling,
        bwtype,
        cfac.dir,
        scale.factor.init,
        ckerbound,
        ckerlb,
        ckerorder,
        ckertype,
        ckerub,
        degree,
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
        dfac.dir,
        dfac.init,
        dfc.dir,
        ftol,
        scale.factor.init.upper,
        hbd.dir,
        hbd.init,
        initc.dir,
        initd.dir,
        invalid.penalty = c("baseline","dbmax"),
        itmax,
        lbc.dir,
        scale.factor.init.lower,
        lbd.dir,
        lbd.init,
        nmulti,
        okertype,
        penalty.multiplier = 10,
        regtype,
        remin,
        scale.init.categorical.sample,
        scale.factor.search.lower = NULL,
        small,
        tol,
        transform.bounds = FALSE,
        ukertype,
        \dots)

\method{npregbw}{rbandwidth}(xdat = stop("invoked without data 'xdat'"),
        ydat = stop("invoked without data 'ydat'"),
        bws,
        bandwidth.compute = TRUE,
        cfac.dir = 2.5*(3.0-sqrt(5)),
        scale.factor.init = 0.5,
        dfac.dir = 0.25*(3.0-sqrt(5)),
        dfac.init = 0.375,
        dfc.dir = 3,
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
        nmulti,
        penalty.multiplier = 10,
        remin = TRUE,
        scale.init.categorical.sample = FALSE,
        scale.factor.search.lower = NULL,
        small = 1.490116e-05,
        tol = 1.490116e-04,
        transform.bounds = FALSE,
        \dots)
}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the data, formula interface, and whether bandwidths are supplied or computed.
}

  \item{bandwidth.compute}{
    a logical value which specifies whether to do a numerical search for
    bandwidths or not. If set to \code{FALSE}, a \code{rbandwidth} object
    will be returned with bandwidths set to those specified
    in \code{bws}. Defaults to \code{TRUE}.
  }

  \item{bws}{
    a bandwidth specification. This can be set as a \code{rbandwidth}
    object returned from a previous invocation, or as a vector of
    bandwidths, with each element \eqn{i} corresponding to the bandwidth
    for column \eqn{i} in \code{xdat}. In either case, the bandwidth
    supplied will serve as a starting point in the numerical search for
    optimal bandwidths. If specified as a vector, then additional
    arguments will need to be supplied as necessary to specify the
    bandwidth type, kernel types, selection methods, and so on. This can
    be left unset.
  }

  \item{call}{
    the original function call. This is passed internally by
    \code{\link{np}} when a bandwidth search has been implied by a call to
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
    a \eqn{p}-variate data frame of regressors on which bandwidth selection will
    be performed. The data types may be continuous, discrete (unordered
    and ordered factors), or some combination thereof.
  }

  \item{ydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{xdat}. 
  }

\subsection{Automatic Degree Search Controls}{
These arguments control automatic local-polynomial degree search when \code{regtype="lp"}.
}

  \item{degree.max}{
    optional scalar or integer vector giving upper bounds for automatic
    degree search when \code{degree.select != "manual"}. If scalar, the
    value is recycled over continuous predictors.
  }

  \item{degree.max.cycles}{
    positive integer giving the maximum number of coordinate-search
    sweeps over the continuous-predictor degree vector. Ignored for
    \code{degree.select="manual"} and \code{"exhaustive"}.
  }

  \item{degree.min}{
    optional scalar or integer vector giving lower bounds for automatic
    degree search when \code{degree.select != "manual"}. If scalar, the
    value is recycled over continuous predictors.
  }

  \item{degree.restarts}{
    non-negative integer giving the number of additional deterministic
    restarts used by coordinate search. Ignored for
    \code{degree.select="manual"} and \code{"exhaustive"}.
  }

  \item{degree.select}{
    character string controlling local-polynomial degree handling when
    \code{regtype="lp"}. \code{"manual"} (default) treats
    \code{degree} as fixed. \code{"coordinate"} performs cached
    coordinate-wise search over admissible degree vectors.
    \code{"exhaustive"} evaluates the full admissible degree grid when
    \code{search.engine="cell"}. For NOMAD-based search engines, any
    non-\code{"manual"} value requests direct joint search over degree
    and bandwidth coordinates.
  }

  \item{degree.start}{
    optional starting degree vector for automatic degree search when
    \code{degree.select="coordinate"}. If omitted, cell-based search starts
    from the degree-zero local-constant baseline on the continuous
    predictors, while NOMAD-based search starts from a clipped
    degree-one vector on the searchable continuous predictors. For
    NOMAD multistarts, later restart starts are generated
    reproducibly from a conservative proposal box and screened using
    \code{dim_basis()} so that the initial basis dimension remains well
    below the training-sample limit. This avoids wasting starts on flat
    penalty or heavily ridged designs while leaving the full user
    requested degree search region unchanged.
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

  \item{bwmethod}{
    which method to use to select bandwidths. \code{cv.aic} specifies
    expected Kullback-Leibler cross-validation (Hurvich, Simonoff, and
    Tsai (1998)), and \code{cv.ls} specifies least-squares
    cross-validation. Defaults to \code{cv.ls}.
  }

  \item{bwscaling}{ a logical value that when set to \code{TRUE} the
    supplied bandwidths are interpreted as \sQuote{scale factors}
    (\eqn{c_j}{c[j]}), otherwise when the value is \code{FALSE} they are
    interpreted as \sQuote{raw bandwidths} (\eqn{h_j}{h[j]} for
    continuous data types, \eqn{\lambda_j}{lambda[j]} for discrete data
    types). For continuous data types, \eqn{c_j}{c[j]} and
    \eqn{h_j}{h[j]} are related by the formula \eqn{h_j = c_j \sigma_j
    n^{-1/(2P+l)}}{h[j] = c[j]*sigma[j]*n^(-1/(2*P+l))}, where
    \eqn{\sigma_j}{sigma[j]} is an adaptive measure of spread of
    continuous variable \eqn{j} defined as min(standard deviation, mean
    absolute deviation/1.4826, interquartile range/1.349), \eqn{n} the
    number of observations, \eqn{P}{P} the order of the kernel, and
    \eqn{l}{l} the number of continuous variables. For discrete data
    types, \eqn{c_j}{c[j]} and \eqn{h_j}{h[j]} are related by the
    formula \eqn{h_j = c_jn^{-2/(2P+l)}}{h[j] = c[j]*n^(-2/(2*P+l))},
    where here \eqn{j}{[j]} denotes discrete variable \eqn{j}{j}.
    Defaults to \code{FALSE}.  }

  \item{bwtype}{
    character string used for the continuous variable bandwidth type,
    specifying the type of bandwidth to compute and return in the
    \code{bandwidth} object. Defaults to \code{fixed}. Option
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

  \item{ckerbound}{
    character string controlling continuous-kernel support handling.
    Can be set as \code{none} (default kernel on full support),
    \code{range} (use sample min/max), or \code{fixed} (use
    \code{ckerlb}/\code{ckerub}). The bounded-kernel route reuses the selected continuous kernel and renormalizes it on the chosen support; see \code{\link{np.kernels}}.
  }

  \item{ckerlb}{
    numeric scalar/vector of lower bounds for continuous variables used
    when \code{ckerbound="fixed"}. Must satisfy lower-bound validity for
    each continuous variable (e.g., \code{<= min(variable)}). Use \code{-Inf} for unbounded below. See \code{\link{np.kernels}} for bounded-kernel normalization details.
  }

  \item{ckerub}{
    numeric scalar/vector of upper bounds for continuous variables used
    when \code{ckerbound="fixed"}. Must satisfy upper-bound validity for
    each continuous variable (e.g., \code{>= max(variable)}). Use \code{Inf} for unbounded above. See \code{\link{np.kernels}} for bounded-kernel normalization details.
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

\subsection{Kernel Type Controls}{
These controls choose continuous, unordered, and ordered kernels.
}

  \item{ckerorder}{
    numeric value specifying kernel order (one of
    \code{(2,4,6,8)}). Kernel order specified along with a
  \code{uniform} continuous kernel type will be ignored. Defaults to
  \code{2}.
  }

  \item{ckertype}{
    character string used to specify the continuous kernel type.
    Can be set as \code{gaussian}, \code{epanechnikov}, or
    \code{uniform}. Defaults to \code{gaussian}.
  }

  \item{okertype}{
    character string used to specify the ordered categorical kernel type.
    Can be set as \code{wangvanryzin}, \code{liracine}, or \code{racineliyan}. Defaults to
    \code{liracine}.
  }

  \item{ukertype}{
    character string used to specify the unordered categorical kernel type.
    Can be set as \code{aitchisonaitken} or \code{liracine}. Defaults to
    \code{aitchisonaitken}.
  }

\subsection{Local-Polynomial Model Specification}{
These arguments control the local-polynomial estimator, basis, and fixed degree specification.
}

  \item{basis}{
    basis selector relevant only when \code{regtype="lp"}.
    Supported values are \code{"glp"}, \code{"additive"}, and
    \code{"tensor"}.
    Let \eqn{d_j}{d[j]} denote the degree for continuous predictor
    \eqn{j}{j} and \eqn{q}{q} the number of continuous predictors.
    With one segment per predictor (no internal knots), basis
    dimensions are:
    \eqn{1+\sum_{j=1}^q d_j}{1+sum(d[j])} for additive,
    \eqn{\prod_{j=1}^q (d_j+1)}{prod(d[j]+1)} for tensor, and
    \eqn{1 + |\{\alpha:\alpha_j \le d_j,\ 0<\sum_j \alpha_j \le \max_j d_j\}|}{1 + number of admissible multi-indices}
    for generalized local-polynomial (GLP) basis construction.
  }

  \item{bernstein.basis}{
    logical flag relevant only when \code{regtype="lp"}. If
    \code{FALSE} (default), the GLP basis uses raw local-polynomial
    powers (stable for extrapolation). If \code{TRUE}, a Bernstein
    (B-spline) basis is used for continuous predictors. For
    \code{bernstein.basis=TRUE}, prediction/evaluation points must lie
    within the training support of each continuous predictor. For
    automatic degree search, if \code{bernstein.basis} is not
    explicitly supplied, the search route defaults to \code{TRUE} for
    numerical stability. Explicit \code{bernstein.basis=FALSE} is
    honored, but raw-polynomial search can be poorly conditioned at
    higher degrees. For
    \code{regtype="ll"} and \code{regtype="lp"}, a pre-optimization
    design-conditioning check is performed on the training continuous
    design: rank deficiency triggers an error, and large condition
    number (\eqn{\kappa(B)}{kappa(B)}) triggers warning/error thresholds
    to avoid unstable optimization dominated by ridging.
  }

  \item{degree}{
    a user-supplied vector of fixed polynomial degrees for the
    continuous predictors (exactly one degree per continuous predictor),
    relevant only when \code{regtype="lp"}. When
    \code{degree.select="manual"}, this must be supplied explicitly.
    Entries must be non-negative integers in \code{[0,12]}.
    Bandwidth optimization treats this vector as fixed input and
    optimizes only bandwidths.
  }

  \item{regtype}{
    a character string specifying which type of kernel regression
    estimator to use. \code{lc} specifies a local-constant estimator
    (Nadaraya-Watson) and \code{ll} specifies a local-linear
    estimator. \code{lp} specifies a local polynomial
    estimator with polynomial degree(s) given by \code{degree} for
    continuous predictors, or selected automatically when
    \code{degree.select != "manual"}. Defaults to \code{lc}.
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
    performs direct joint mixed discrete/continuous search over fixed
    bandwidths and the degree vector using \code{crs::snomadr()},
    followed by one Powell hot start from the NOMAD solution.
    \code{"nomad"} performs the direct joint NOMAD search without the
    Powell refinement. \code{"cell"} uses the legacy profiled
    degree-grid search built from repeated fixed-degree bandwidth
    solves. NOMAD-based search currently requires \code{bwtype="fixed"},
    \code{degree.verify=FALSE}, and the suggested package
    \pkg{crs} to be installed.
  }

\subsection{Numerical Search And Tolerance Controls}{
These controls set optimizer tolerances, restart behavior, invalid-candidate penalties, and bounded search transformations.
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
    the cross-validation function from different (random) initial
    points. Defaults to \code{min(2,ncol(xdat))}.
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

  Documentation guide: see \code{\link{np.kernels}} for kernels,
  \code{\link{np.options}} for global options, and \code{\link{plot}}
  for plotting options.

  The bandwidth-selection argument surface is easiest to read by
  decision group: data and existing bandwidth inputs;
  local-polynomial/NOMAD controls when polynomial-adaptive regression
  is requested; bandwidth criterion and representation; continuous
  kernel and support controls beginning with \code{cker*};
  categorical kernel controls \code{ukertype} and \code{okertype}; and
  numerical search initialization, tolerances, and feasibility
  controls. Users who call \code{\link{npreg}} without a bandwidth
  object can pass these same bandwidth-selection controls through that
  function's \code{\dots}.

  For S3 plotting help, use \code{methods("plot")} and query
  class-specific help topics such as \code{?plot.npregression} and
  \code{?plot.rbandwidth}. You can inspect implementations with
  \code{getS3method("plot","npregression")}.

  \code{npregbw} implements a variety of methods for choosing
  bandwidths for multivariate (\eqn{p}-variate) regression data defined
  over a set of possibly continuous and/or discrete (unordered, ordered)
  data. The approach is based on Li and Racine (2003) who employ
  \sQuote{generalized product kernels} that admit a mix of continuous
  and discrete data types.

  The cross-validation methods employ multivariate numerical search
  algorithms. For fixed-degree local-constant/local-linear regression,
  and for local-polynomial regression with \code{degree.select="manual"},
  the bandwidth search uses multidimensional Powell direction-set
  optimization.

  Bandwidths can (and will) differ for each variable which is, of
  course, desirable.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating the
  density at the point \eqn{x}. Generalized nearest-neighbor bandwidths change
  with the point at which the density is estimated, \eqn{x}. Fixed bandwidths
  are constant over the support of \eqn{x}.

  \code{npregbw} may be invoked \emph{either} with a formula-like
  symbolic 
  description of variables on which bandwidth selection is to be
  performed \emph{or} through a simpler interface whereby data is passed
  directly to the function via the \code{xdat} and \code{ydat}
  parameters. Use of these two interfaces is \bold{mutually exclusive}.

  Data contained in the data frame \code{xdat} may be a mix of
  continuous (default), unordered discrete (to be specified in the data
  frame \code{xdat} using \code{\link{factor}}), and ordered discrete
  (to be specified in the data frame \code{xdat} using
  \code{\link{ordered}}). Data can be entered in an arbitrary order and
  data types will be detected automatically by the routine (see
  \code{\link{np}} for details).

  Data for which bandwidths are to be estimated may be specified
  symbolically. A typical description has the form \code{dependent data
    ~ explanatory data},
  where \code{dependent data} is a univariate response, and
  \code{explanatory data} is a
  series of variables specified by name, separated by 
  the separation character '+'. For example, \code{ y1 ~ x1 + x2 }
  specifies that the bandwidths for the regression of response \code{y1}
  and 
  nonparametric regressors \code{x1} and \code{x2} are to be estimated.
  See below for further examples. 

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and eighth
  order Gaussian and Epanechnikov kernels, and the uniform
  kernel. Unordered discrete data types use a variation on Aitchison and
  Aitken's (1976) kernel, while ordered data types use a variation of the
  Wang and van Ryzin (1981) kernel.

  When \code{regtype="lp"} and \code{degree.select != "manual"},
  \code{npregbw} can jointly determine the continuous-predictor degree
  vector and bandwidth coordinates. With \code{search.engine="cell"},
  the objective is profiled over the degree grid using cached
  coordinate-wise or exhaustive search together with the existing
  fixed-degree bandwidth optimizer. With
  \code{search.engine="nomad"} or \code{"nomad+powell"}, the package
  instead evaluates the cross-validation criterion directly over the
  joint space of fixed bandwidths and polynomial degrees using
  \code{crs::snomadr()}. \code{"nomad+powell"} then performs one Powell
  hot start from the NOMAD solution and retains the better of the
  direct NOMAD and polished solutions. This direct joint-search route
  follows the polynomial-adaptive cross-validation rationale of Hall
  and Racine (2015). When \code{bernstein.basis} is not explicitly
  supplied, the automatic search route defaults to
  \code{bernstein.basis=TRUE} for numerical stability; explicit
  \code{bernstein.basis=FALSE} is honored but can be poorly conditioned
  at higher degrees. NOMAD multistarts are initialized more
  conservatively than the full degree search box: start 1 is the
  user-supplied degree/bandwidth vector when provided and otherwise a
  clipped degree-one vector, while later starts are reproducible random
  draws from a reduced degree proposal box whose candidates are screened
  using \code{dim_basis()}. This heuristic is used only to obtain
  feasible, numerically safer, and quicker initial evaluations; it does
  not restrict the admissible degree region searched by NOMAD. The
  direct NOMAD backend is provided by the suggested package
  \code{crs}, so install \code{crs} before using
  \code{search.engine="nomad"}, \code{"nomad+powell"}, or
  \code{nomad=TRUE}.

  Setting \code{nomad=TRUE} is a convenience preset for this automatic
  LP route, not a generic optimizer alias. For regression it expands any
  missing values to the equivalent long-form call
  \preformatted{
    npregbw(...,
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
  When the direct NOMAD route is active, \code{nmulti} controls the
  package-level outer restart count while \code{nomad.nmulti}
  controls the inner \code{crs::snomadr()} multistart count used within
  each outer restart. The default \code{nomad.nmulti=0L} preserves the
  current single-start inner NOMAD behavior.

  The use of compactly supported kernels or the occurrence of small
  bandwidths during cross-validation can lead to numerical problems for
  the local linear estimator when computing the locally weighted least
  squares solution. To overcome this problem we rely on a form or
  \sQuote{ridging} proposed by Cheng, Hall, and Titterington (1997),
  modified so that we solve the problem pointwise rather than globally
  (i.e. only when it is needed).

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
  \sQuote{tuning} for various methods rather than for the user
  themselves.

}
\value{
  \code{npregbw} returns a \code{rbandwidth} object, with the
  following components:

  \item{bw}{ bandwidth(s), scale factor(s) or nearest neighbours for the
    data, \code{xdat} }
  \item{fval}{ objective function value at minimum }

  if \code{bwtype} is set to \code{fixed}, an object containing bandwidths
  (or scale factors if \code{bwscaling = TRUE}) is returned. If it is set to
  \code{generalized_nn} or \code{adaptive_nn}, then instead the \eqn{k}th nearest
  neighbors are returned for the continuous variables while the discrete
  kernel bandwidths are returned for the discrete variables. Bandwidths
  are stored under the component name \code{bw}, with each
  element \eqn{i} corresponding to column \eqn{i} of input data
  \code{xdat}.

  The functions \code{\link{predict}}, \code{\link{summary}}, and \code{\link{plot}} support
  objects of this class.
}

\references{
  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Cheng, M.-Y. and P. Hall and D.M. Titterington (1997), \dQuote{On the
  shrinkage of local linear curve estimators,} Statistics and Computing,
  7, 11-17.

  Fan, J. and I. Gijbels (1996), \emph{Local Polynomial Modelling and
  Its Applications,} Chapman and Hall.

  Hall, P. and J.S. Racine (2015), \dQuote{Infinite Order Cross-Validated
  Local Polynomial Regression,} Journal of Econometrics, 185, 510-525.

  Hall, P. and Q. Li and J.S. Racine (2007), \dQuote{Nonparametric
  estimation of regression functions in the presence of irrelevant
  regressors,} The Review of Economics and Statistics, 89, 784-789.

  Hurvich, C.M. and J.S. Simonoff and C.L. Tsai (1998),
  \dQuote{Smoothing parameter selection in nonparametric regression
    using an improved Akaike information criterion,} Journal of the
  Royal Statistical Society B, 60, 271-293.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated local linear
    nonparametric regression,} Statistica Sinica, 14, 485-512.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.
  
  Racine, J.S. and Q. Li (2004), \dQuote{Nonparametric estimation of regression
    functions with both categorical and continuous data,} Journal of
  Econometrics, 119, 99-130.
  
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
  conduct multistarting (the default is to restart min(2, ncol(xdat))
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
  \code{\link{npreg}}
}
\examples{
\dontrun{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we compute a
# Bivariate nonparametric regression estimate for Giovanni Baiocchi's
# Italian income panel (see Italy for details)

data("Italy")
attach(Italy)

# Compute the least-squares cross-validated bandwidths for the local
# constant estimator (default)

bw <- npregbw(formula=gdp~ordered(year))

summary(bw)

# Sleep for 5 seconds so that we can examine the output...

if (interactive()) Sys.sleep(5)

# Supply your own bandwidth...

bw <- npregbw(formula=gdp~ordered(year), bws=c(0.75),
              bandwidth.compute=FALSE)

summary(bw)

# Sleep for 5 seconds so that we can examine the output...

if (interactive()) Sys.sleep(5)

# Treat year as continuous and supply your own scaling factor c in
# c sigma n^{-1/(2p+q)}

bw <- npregbw(formula=gdp~year, bws=c(1.06),
              bandwidth.compute=FALSE, 
              bwscaling=TRUE)

summary(bw)

# Note - see also the example for npudensbw() for more extensive
# multiple illustrations of how to change the kernel function, kernel
# order, bandwidth type and so forth.

detach(Italy)

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we compute a
# Bivariate nonparametric regression estimate for Giovanni Baiocchi's
# Italian income panel (see Italy for details)

data("Italy")
attach(Italy)

# Compute the least-squares cross-validated bandwidths for the local
# constant estimator (default)

bw <- npregbw(xdat=ordered(year), ydat=gdp)

summary(bw)

# Sleep for 5 seconds so that we can examine the output...

if (interactive()) Sys.sleep(5)

# Supply your own bandwidth...

bw <- npregbw(xdat=ordered(year), ydat=gdp, bws=c(0.75),
              bandwidth.compute=FALSE)

summary(bw)

# Sleep for 5 seconds so that we can examine the output...

if (interactive()) Sys.sleep(5)

# Treat year as continuous and supply your own scaling factor c in
# c sigma n^{-1/(2p+q)}

bw <- npregbw(xdat=year, ydat=gdp, bws=c(1.06),
              bandwidth.compute=FALSE, 
              bwscaling=TRUE)

summary(bw)

# Note - see also the example for npudensbw() for more extensive
# multiple illustrations of how to change the kernel function, kernel
# order, bandwidth type and so forth.

detach(Italy)
} % enddontrun
}
\keyword{nonparametric}
