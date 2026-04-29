%
\name{npscoefbw}
\alias{npscoefbw}
\alias{npscoefbw.default}
\alias{npscoefbw.formula}
\alias{npscoefbw.scbandwidth}
\title{Smooth Coefficient Kernel Regression Bandwidth Selection}

\description{
  \code{npscoefbw} computes a bandwidth object for a smooth
  coefficient kernel regression estimate of a one (1) dimensional
  dependent variable on
  \eqn{p+q}-variate explanatory data, using the model
  \eqn{Y_i = W_{i}^{\prime} \gamma (Z_i) + u_i}{Y_i = t(W_i) * gamma(Z_i) +
    u_i} where \eqn{W_i'=(1,X_i')}{t(W_i) = (1,t(X_i))}
  given training points (consisting of explanatory data and dependent
  data), and a bandwidth specification, which can be a \code{rbandwidth}
  object, or a bandwidth vector, bandwidth type and kernel type.
}

\usage{
npscoefbw(\dots)

\method{npscoefbw}{formula}(formula, 
          data, 
          subset, 
          na.action, 
          call, 
          \dots)

\method{npscoefbw}{default}(xdat = stop("invoked without data 'xdat'"),
          ydat = stop("invoked without data 'ydat'"),
          zdat = NULL,
          bws,
          backfit.iterate,
          backfit.maxiter,
          backfit.tol,
          bandwidth.compute = TRUE,
          basis,
          bernstein.basis,
          bwmethod,
          bwscaling,
          bwtype,
          ckerbound,
          ckerlb,
          ckerorder,
          ckertype,
          ckerub,
          cv.iterate,
          cv.num.iterations,
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
          nmulti,
          okertype,
          optim.abstol,
          optim.maxattempts,
          optim.maxit,
          optim.method,
          optim.reltol,
          random.seed,
          regtype,
          ukertype,
          scale.factor.init.lower = 0.1,
          scale.factor.init.upper = 2.0,
          scale.factor.init = 0.5,
          lbd.init = 0.5,
          hbd.init = 1.5,
          dfac.init = 1.0,
          scale.factor.search.lower = NULL,
          \dots)

\method{npscoefbw}{scbandwidth}(xdat = stop("invoked without data 'xdat'"),
          ydat = stop("invoked without data 'ydat'"),
          zdat = NULL,
          bws,
          backfit.iterate = FALSE,
          backfit.maxiter = 100,
          backfit.tol = .Machine$double.eps,
          bandwidth.compute = TRUE,
          cv.iterate = FALSE,
          cv.num.iterations = 1,
          nmulti,
          optim.abstol = .Machine$double.eps,
          optim.maxattempts = 10,
          optim.maxit = 500,
          optim.method = c("Nelder-Mead", "BFGS", "CG"),
          optim.reltol = sqrt(.Machine$double.eps),
          random.seed = 42,
          scale.factor.init.lower = 0.1,
          scale.factor.init.upper = 2.0,
          scale.factor.init = 0.5,
          lbd.init = 0.5,
          hbd.init = 1.5,
          dfac.init = 1.0,
          scale.factor.search.lower = NULL,
          \dots)
}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the smooth-coefficient data, formula interface, and whether bandwidths are supplied or computed.
}

  \item{bandwidth.compute}{
    a logical value which specifies whether to do a numerical search for
    bandwidths or not. If set to \code{FALSE}, a \code{scbandwidth} object
    will be returned with bandwidths set to those specified
    in \code{bws}. Defaults to \code{TRUE}.
  }

  \item{bws}{
    a bandwidth specification. This can be set as a \code{scbandwidth}
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
    which, by default, populates the columns \eqn{2} through \eqn{p+1}
    of \eqn{W} in the model equation, and in the
    absence of \code{zdat}, will also correspond to
    \eqn{Z} from the model equation.
  }

  \item{ydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{xdat}.
  }

  \item{zdat}{
    an optionally specified \eqn{q}-variate data frame of explanatory
    data (training data), which corresponds to \eqn{Z}
    in the model equation. Defaults to be the same as \code{xdat}.
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
    coordinate-wise search over admissible degree vectors for the
    continuous \code{z} variables. \code{"exhaustive"} evaluates the
    full admissible degree grid when \code{search.engine="cell"}.
    For NOMAD-based search engines, any non-\code{"manual"} value
    requests direct joint search over degree and bandwidth
    coordinates.
  }

  \item{degree.start}{
    optional starting degree vector for automatic coordinate search. If
    omitted, the search starts from the degree-zero local-constant
    baseline on the continuous \code{z} variables.
  }

  \item{degree.verify}{
    logical value indicating whether a coordinate-search solution should
    be exhaustively verified over the admissible degree grid after the
    heuristic phase completes. Available only for
    \code{search.engine="cell"}.
  }


\subsection{Backfitting Controls}{
These controls tune the optional smooth-coefficient backfitting iterations.
}

  \item{backfit.iterate}{
    boolean value specifying whether or not to iterate evaluations of
    the smooth coefficient estimator, for extra accuracy, during the
    cross-validated backfitting procedure. Defaults to \code{FALSE}.
  }

  \item{backfit.maxiter}{
    integer specifying the maximum number of times to iterate the
    evaluation of the smooth coefficient estimator in the attempt to
    obtain the desired accuracy. Defaults to \code{100}.
  }

  \item{backfit.tol}{
    tolerance to determine convergence of iterated evaluations of the
    smooth coefficient estimator. Defaults to \code{.Machine$double.eps}.
  }


\subsection{Bandwidth Criterion And Representation}{
These arguments choose the selection criterion and the way continuous bandwidths are represented.
}

  \item{bwmethod}{
    which method was used to select bandwidths. \code{cv.ls}
    specifies least-squares cross-validation, which is all that is
    currently supported. Defaults to \code{cv.ls}.
  }

  \item{bwscaling}{ a logical value that when set to \code{TRUE} the
    supplied bandwidths are interpreted as \sQuote{scale factors}
    (\eqn{c_j}{c[j]}), otherwise when the value is \code{FALSE} they are
    interpreted as \sQuote{raw bandwidths} (\eqn{h_j}{h[j]} for continuous data
    types, \eqn{\lambda_j}{lambda[j]} for discrete data types). For
    continuous data types, \eqn{c_j}{c[j]} and \eqn{h_j}{h[j]} are
    related by the formula \eqn{h_j = c_j \sigma_j n^{-1/(2P+l)}}{h[j] =
    c[j]*sigma[j]*n^(-1/(2*P+l))}, where \eqn{\sigma_j}{sigma[j]} is an
    adaptive measure of spread of continuous variable \eqn{j} defined as
    min(standard deviation, mean absolute deviation, interquartile
    range/1.349), \eqn{n} the number of observations, \eqn{P}{P} the
    order of the kernel, and \eqn{l}{l} the number of continuous
    variables. For discrete data types, \eqn{c_j}{c[j]} and
    \eqn{h_j}{h[j]} are related by the formula \eqn{h_j =
    c_jn^{-2/(2P+l)}}{h[j] = c[j]*n^(-2/(2*P+l))}, where here
    \eqn{j}{[j]} denotes discrete variable \eqn{j}{j}.  Defaults to
    \code{FALSE}.  }

  \item{bwtype}{
    character string used for the continuous variable bandwidth type,
    specifying the type of bandwidth provided. Defaults to
    \code{fixed}. Option summary:\cr
    \code{fixed}: fixed bandwidths or scale factors \cr
    \code{generalized_nn}: generalized nearest neighbors \cr
    \code{adaptive_nn}:  adaptive nearest neighbors
  }


\subsection{Categorical Search Initialization}{
These controls set categorical search starts.
}

  \item{dfac.init}{
    deterministic fixed-bandwidth start factor for ordered and
    unordered categorical coordinates. Used only when
    \code{bwtype="fixed"}. Defaults to \code{1.0}. Values must not
    exceed \code{2}.
  }

  \item{hbd.init}{
    upper bound for random fixed-bandwidth start factors for ordered
    and unordered categorical coordinates. Used only when
    \code{bwtype="fixed"}. Defaults to \code{1.5}. Must be greater than
    or equal to \code{lbd.init} and must not exceed \code{2}.
  }

  \item{lbd.init}{
    lower bound for random fixed-bandwidth start factors for ordered
    and unordered categorical coordinates. Used only when
    \code{bwtype="fixed"}. Defaults to \code{0.5}. Values must not
    exceed \code{2} so that categorical fixed starts remain within
    lawful bandwidth bounds.
  }


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

\subsection{Cross-Validation Iteration Controls}{
These controls tune iterative cross-validation behavior.
}

  \item{cv.iterate}{
    boolean value specifying whether or not to perform iterative,
    cross-validated backfitting on the data. See details for limitations
    of the backfitting procedure. Defaults to \code{FALSE}.
  }

  \item{cv.num.iterations}{
    integer specifying the number of times to iterate the backfitting
    process over all covariates. Defaults to \code{1}.
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
    Can be set as \code{wangvanryzin}, \code{liracine}, or
    \code{racineliyan}. Defaults to \code{liracine}.
  }
  \item{ukertype}{
    character string used to specify the unordered categorical kernel type.
    Can be set as \code{aitchisonaitken} or \code{liracine}. Defaults
    to \code{aitchisonaitken}.
  }

\subsection{Local-Polynomial Model Specification}{
These arguments control the local-polynomial estimator, basis, and fixed degree specification.
}

  \item{basis}{
    for \code{regtype="lp"}, the polynomial basis family used for local
    polynomial fitting. Options are \code{"glp"} (default),
    \code{"additive"}, and \code{"tensor"}.
  }

  \item{bernstein.basis}{
    for \code{regtype="lp"}, logical flag selecting Bernstein-basis
    representation where supported. When automatic degree search is
    requested and \code{bernstein.basis} is not explicitly supplied,
    the search route defaults to \code{TRUE} for numerical stability.
    Explicit \code{bernstein.basis=FALSE} is honored, but raw-polynomial
    search can be poorly conditioned at higher degrees.
  }

  \item{degree}{
    for \code{regtype="lp"}, polynomial degree specification for each
    continuous \code{z} variable.
  }

  \item{regtype}{
    a character string specifying local smoothing type for the \code{z}
    surface. Options are \code{"lc"} (default), \code{"ll"}, and
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
    performs direct joint search over the \code{zdat}-side continuous
    bandwidth coordinates and degree vector using
    \code{crs::snomadr()}, then applies one Powell hot start from the
    NOMAD solution. \code{"nomad"} omits the Powell refinement.
    \code{"cell"} profiles the criterion over the admissible degree
    grid using repeated fixed-degree bandwidth solves. NOMAD-based
    search currently requires \code{bwtype="fixed"},
    \code{degree.verify=FALSE}, and the suggested package \pkg{crs} to
    be installed.
  }


\subsection{Numerical Search Controls}{
These controls set search restart behavior.
}

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points. Defaults to \code{min(2,ncol(xdat))}.
  }


\subsection{Optimization Controls}{
These arguments control outer optimization behavior for the semiparametric search.
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
  additional arguments supplied to specify the regression type,
  bandwidth type, kernel types, selection methods, and so on, detailed
  below.
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

  \code{npscoefbw} implements a variety of methods for semiparametric
 regression on multivariate (\eqn{p+q}-variate) explanatory data defined
 over a set of possibly continuous data. The approach is based on Li and
 Racine (2003) who employ \sQuote{generalized product kernels} that
 admit a mix of continuous and discrete data types.

 Three classes of kernel estimators for the continuous data types are
 available: fixed, adaptive nearest-neighbor, and generalized
 nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
 each sample realization in the set, \eqn{x_i}{x[i]}, when estimating the
 density at the point \eqn{x}. Generalized nearest-neighbor bandwidths change
 with the point at which the density is estimated, \eqn{x}. Fixed bandwidths
 are constant over the support of \eqn{x}.

 \code{npscoefbw} may be invoked \emph{either} with a formula-like
 symbolic description of variables on which bandwidth selection is to be
 performed \emph{or} through a simpler interface whereby data is passed
 directly to the function via the \code{xdat}, \code{ydat}, and
 \code{zdat} parameters. Use of these two interfaces is \bold{mutually
 exclusive}.

 Data contained in the data frame \code{xdat} may be continuous and in
 \code{zdat} may be of mixed type. Data can be entered in an arbitrary
 order and data types will be detected automatically by the routine (see
 \code{\link{npRmpi}} for details).

 Data for which bandwidths are to be estimated may be specified
 symbolically. A typical description has the form \code{dependent}
 \code{data} \code{~} \code{parametric} \code{explanatory} \code{data}
 \code{|} \code{nonparametric} \code{explanatory} \code{data}, where
 \code{dependent} \code{data} is a univariate response, and
 \code{parametric} \code{explanatory} \code{data} and
 \code{nonparametric} \code{explanatory} \code{data} are both series of
 variables specified by name, separated by the separation character
 '+'. For example, \code{ y1 ~ x1 + x2 | z1 } specifies that the
 bandwidth object for the smooth coefficient model with response
 \code{y1}, linear parametric regressors \code{x1} and \code{x2}, and
 nonparametric regressor (that is, the slope-changing variable)
 \code{z1} is to be estimated. See below for further examples.  In the
 case where the nonparametric (slope-changing) variable is not
 specified, it is assumed to be the same as the parametric variable.

 A variety of kernels may be specified by the user. Kernels implemented
 for continuous data types include the second, fourth, sixth, and eighth
 order Gaussian and Epanechnikov kernels, and the uniform
 kernel. Unordered discrete data types use a variation on Aitchison and
 Aitken's (1976) kernel, while ordered data types use a variation of the
 Wang and van Ryzin (1981) kernel.

 Setting \code{nomad=TRUE} is a convenience preset for this automatic
 LP route, not a generic optimizer alias. For smooth coefficient
 regression it expands any missing values to the equivalent long-form
 call
 \preformatted{
   npscoefbw(...,
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

  When \code{regtype="lp"} and \code{degree.select != "manual"},
  \code{npscoefbw} can jointly determine the \code{zdat}-side local
  polynomial degree vector together with the associated bandwidth
  coordinates. With \code{search.engine="cell"}, the criterion is
  profiled over the admissible degree grid using cached
  coordinate-wise or exhaustive search together with repeated
  fixed-degree bandwidth solves. With \code{search.engine="nomad"} or
  \code{"nomad+powell"}, the criterion is optimized directly over the
  joint degree/bandwidth space using \code{crs::snomadr()};
  \code{"nomad+powell"} then performs one Powell hot start from the
  NOMAD solution and keeps the better of the direct NOMAD and polished
  answers. This polynomial-adaptive joint-search route is motivated by
  Hall and Racine (2015). When \code{bernstein.basis} is not explicitly
  supplied, the automatic search route defaults to
  \code{bernstein.basis=TRUE} for numerical stability.

}

\value{

 if \code{bwtype} is set to \code{fixed}, an object containing
 bandwidths (or scale factors if \code{bwscaling = TRUE}) is
 returned. If it is set to \code{generalized_nn} or \code{adaptive_nn},
 then instead the \eqn{k}th nearest neighbors are returned for the
 continuous variables while the discrete kernel bandwidths are returned
 for the discrete variables. Bandwidths are stored in a vector under the
 component name \code{bw}. Backfitted bandwidths are stored under the
 component name \code{bw.fitted}.

  The functions \code{\link{predict}}, \code{\link{summary}}, and
  \code{\link{plot}} support
  objects of this class.

}

\references{

 Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
   discrimination by the kernel method,} Biometrika, 63, 413-420.

 Cai Z. (2007), \dQuote{Trending time-varying coefficient time series
   models with serially correlated errors,} Journal of Econometrics,
 136, 163-188.

 Hastie, T. and R. Tibshirani (1993), \dQuote{Varying-coefficient
  models,} Journal of the Royal Statistical Society, B 55, 757-796.

 Hall, P. and J.S. Racine (2015), \dQuote{Infinite Order Cross-Validated
 Local Polynomial Regression,} Journal of Econometrics, 185, 510-525.

 Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
 and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2010), \dQuote{Smooth varying-coefficient
    estimation and inference for qualitative and quantitative data,}
    Econometric Theory, 26, 1-31.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
 Cambridge University Press.

 Li, Q. and D. Ouyang and J.S. Racine (2013), \dQuote{Categorical
 semiparametric varying-coefficient models,} Journal of Applied
 Econometrics, 28, 551-589.

 Li, A. and Q. Li and J.S. Racine (under revision), \dQuote{Boundary
 Adjusted, Polynomial Adaptive, Nonparametric Kernel Conditional
 Density Estimation,} Econometric Reviews.

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
  optim.reltol=.1 and conduct multistarting (the default is to restart
  min(2,ncol(zdat)) times). Once the procedure terminates, you can restart
  search with default tolerances using those bandwidths obtained from
  the less rigorous search (i.e., set \code{bws=bw} on subsequent calls
  to this routine where \code{bw} is the initial bandwidth object).  A
  version of this package using the \code{Rmpi} wrapper is under
  development that allows one to deploy this software in a clustered
  computing environment to facilitate computation involving large
  datasets.

  Support for backfitted bandwidths is experimental and is limited in
  functionality. The code does not support asymptotic standard errors
  or out of sample estimates with backfitting.

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

n <- 500

x <- runif(n)
z <- runif(n, min=-2, max=2)
y <- x*exp(z)*(1.0+rnorm(n,sd = 0.2))

## A smooth coefficient model example

bw <- npscoefbw(y~x|z)

summary(bw)

model <- npscoef(bws=bw, gradients=TRUE)

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
