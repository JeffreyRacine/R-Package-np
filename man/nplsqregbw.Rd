\name{nplsqregbw}
\alias{nplsqregbw}
\alias{nplsqregbw.formula}
\alias{nplsqregbw.default}
\alias{nplsqregbw.lsqregressionbandwidth}
\alias{lsqregressionbandwidth}
\alias{print.lsqregressionbandwidth}
\alias{summary.lsqregressionbandwidth}
\title{Bandwidth Selection for Location-Scale Kernel Quantile Regression}
\description{
  \code{nplsqregbw} selects bandwidths and the scalar location-scale shift
  parameter \code{delta} for \code{\link{nplsqreg}} by minimizing
  leave-one-out check loss.
}
\usage{
nplsqregbw(bws, \dots)

\method{nplsqregbw}{formula}(bws, data = NULL, tau = 0.5,
       subset, na.action, \dots)

\method{nplsqregbw}{lsqregressionbandwidth}(bws, tau = bws$tau, \dots)

\method{nplsqregbw}{default}(bws,
       xdat = stop("invoked without data 'xdat'"),
       ydat = stop("invoked without data 'ydat'"),
       tau = 0.5,
       tau.search = c("full", "refined"),
       delta = NULL,
       scale = NULL,
       regtype = c("lc", "ll", "lp"),
       regtype.pilot = c("auto", "ll", "lc", "lp"),
       nomad = FALSE,
       nomad.pilot = FALSE,
       pilot.args = list(),
       bandwidth.compute = TRUE,
       delta.bounds = c(1e-4, 1 - 1e-4),
       optim.control = list(maxit = 50L),
       \dots,
       nomad.opts = list())
}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification, formula/data interface,
and training data.
}

  \item{bws}{
    a formula, an \code{lsqregressionbandwidth} object, an
    \code{rbandwidth} object, a numeric bandwidth vector, or omitted for
    automatic selection. Ordinary \code{\link{npregbw}} \code{rbandwidth}
    objects may be used as starting values. Exact nplsqreg reuse is through
    an \code{lsqregressionbandwidth} object, typically \code{fit$bws};
    \code{fit$reg.bws} is internal regression state.
  }

  \item{data}{
    an optional data frame, list or environment containing the variables
    in the model. If not found in \code{data}, the variables are taken
    from \code{environment(bws)}.
  }

  \item{subset}{
    an optional vector specifying a subset of observations to be used by
    the formula method.
  }

  \item{na.action}{
    a function specifying the action to take when missing values are found
    by the formula method.
  }

  \item{xdat}{
    a \eqn{p}-variate data frame of explanatory data used as training data.
  }

  \item{ydat}{
    a one dimensional numeric vector of dependent data.
  }

\subsection{Quantile Index And Vector-Tau Search}{
These arguments identify the quantile probabilities and how vector-\code{tau}
fits are orchestrated.
}

  \item{tau}{
    a numeric scalar or vector specifying the quantile probability or
    probabilities \eqn{\tau}{tau}. Values must lie strictly in \eqn{(0,1)}.
  }

  \item{tau.search}{
    either \code{"full"} or \code{"refined"}. The default
    \code{"full"} route performs a full bandwidth/\code{delta} search
    for each \code{tau}, sharing only the pilot scale. The explicit
    \code{"refined"} route fits the central \code{tau} first and uses
    the resulting bandwidths, \code{delta}, and when applicable selected
    degree vector as warm starts for the remaining quantiles.
  }

\subsection{Scale Pilot}{
These arguments control or supply the conditional standard deviation used in
the location-scale transformation.
}

  \item{scale}{
    an optional strictly positive scale vector interpreted as the
    conditional standard deviation at the training observations. If
    omitted, a pilot scale is computed automatically.
  }

  \item{regtype.pilot}{
    regression type used for the pilot mean and residual-variance
    smooths. The default \code{"auto"} uses local linear smoothing when
    continuous predictors are present and local constant smoothing for
    categorical-only designs.
  }

  \item{nomad.pilot}{
    logical value indicating whether the pilot regressions should use the
    automatic NOMAD local-polynomial route when compatible with
    \code{regtype.pilot}. Defaults to \code{FALSE}.
  }

  \item{pilot.args}{
    optional named list of additional arguments supplied to the pilot
    \code{\link{npreg}} calls.
  }

\subsection{Bandwidth, Delta And Degree Search}{
These arguments control the main check-loss bandwidth and
\code{delta} optimization.
}

  \item{regtype}{
    regression type for the main transformed-response fit and check-loss
    bandwidth search. This is independent of \code{regtype.pilot} and
    follows \code{\link{npregbw}} semantics.
  }

  \item{delta}{
    optional starting value for the shift parameter. If omitted, the
    first start uses \code{0.5}; subsequent multistarts draw starting
    values from \eqn{(0,1)}.
  }

  \item{delta.bounds}{
    two numeric values giving the lower and upper bounds for
    \code{delta}. The interval must lie strictly inside \eqn{(0,1)}.
  }

  \item{nomad}{
    logical value indicating whether to use the local-polynomial NOMAD
    route for joint bandwidth, \code{delta}, and continuous degree-vector
    search. Defaults to \code{FALSE}.
  }

  \item{nomad.opts}{
    optional named list of native NOMAD parameters forwarded to the
    \pkg{crs} NOMAD C API for NOMAD-backed bandwidth searches.
    \code{NULL} is treated as no extra NOMAD options. Non-list values and
    invalid or unsupported parameters are rejected before solver entry.
  }

  \item{bandwidth.compute}{
    logical value indicating whether automatic bandwidth selection should
    be performed. Defaults to \code{TRUE}.
  }

  \item{optim.control}{
    a list of controls passed to the Powell optimizer used for fixed
    degree searches and NOMAD hot-start refinement.
  }

  \item{\dots}{
    additional bandwidth, kernel, local-polynomial, and search controls
    forwarded to the package bandwidth machinery. Common examples include
    \code{regtype}, \code{bwtype}, \code{bwmethod}, \code{nmulti},
    \code{degree}, \code{basis}, \code{bernstein.basis},
    \code{search.engine}, \code{degree.min}, \code{degree.max},
    \code{degree.start}, \code{degree.restarts}, \code{powell.remin},
    \code{nomad.remin}, and \code{nomad.nmulti}.
    The explanatory-variable kernel aliases used by \code{\link{npqreg}}
    and \code{\link{npcdistbw}}, such as \code{cxkertype},
    \code{cxkerorder}, \code{uxkertype}, and \code{oxkertype}, are
    accepted and mapped to the corresponding regression-kernel controls.
    Response-side conditional-distribution kernel controls such as
    \code{cykertype} are not meaningful for the location-scale
    transformed-response estimator and fail clearly rather than being
    ignored.
  }
}
\details{
  For a requested quantile probability \eqn{\tau}{tau},
  \code{nplsqregbw} chooses bandwidths and \code{delta} for the
  transformed response
  \deqn{
    Y_i^\delta = Y_i + \hat\sigma(X_i)\Phi^{-1}(\delta)
  }
  by minimizing the leave-one-out check loss. The selected bandwidths
  are then used by \code{\link{nplsqreg}} for the final mixed-data
  \code{\link{npreg}} fit of the transformed response.

  The default pilot scale is a residual scale: first \code{\link{npreg}}
  estimates the conditional mean, then \code{\link{npreg}} smooths
  squared residuals using the same bandwidth object, and the square root
  of the floored fitted variance is used as the scale. The local-linear
  pilot option follows Fan and Yao (1998). A user-supplied \code{scale}
  must be strictly positive and is interpreted as a conditional standard
  deviation, not a variance.

  When \code{nomad=TRUE}, the search can include the continuous
  local-polynomial degree vector in addition to bandwidths and
  \code{delta}. The \code{search.engine="nomad+powell"} route uses NOMAD
  for degree search and Powell for hot-start refinement at the selected
  degree.
}
\value{
  \code{nplsqregbw} returns an object of class
  \code{lsqregressionbandwidth}. For vector \code{tau}, the object stores
  the per-\code{tau} bandwidth objects, selected \code{delta} values,
  objectives, fit order, warm-start provenance, and shared pilot scale.
}
\references{
  Fan, J. and Q. Yao (1998), \dQuote{Efficient Estimation of Conditional
  Variance Functions in Stochastic Regression,} Biometrika, 85, 645-660.
  \doi{10.1093/biomet/85.3.645}

  Racine, J.S. and K. Li (2017), \dQuote{Nonparametric conditional quantile
  estimation: A locally weighted quantile kernel approach,} Journal of
  Econometrics, 201, 72-94. \doi{10.1016/j.jeconom.2017.06.020}

  Racine, J.S. and I. Van Keilegom (2020), \dQuote{A smooth
  nonparametric, multivariate, mixed-data location-scale test,} Journal
  of Business & Economic Statistics, 38, 784-795.
  \doi{10.1080/07350015.2019.1574227}
}
\author{
  Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\seealso{
  \code{\link{nplsqreg}}, \code{\link{npregbw}}, \code{\link{npreg}},
  \code{\link{np.kernels}}, \code{\link{np.options}}
}
\keyword{nonparametric}
\keyword{regression}
