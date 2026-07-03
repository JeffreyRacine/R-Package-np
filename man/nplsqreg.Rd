\name{nplsqreg}
\alias{nplsqreg}
\alias{nplsqreg.formula}
\alias{nplsqreg.default}
\alias{nplsqreg.lsqregressionbandwidth}
\alias{fitted.lsqregression}
\alias{quantile.lsqregression}
\alias{predict.lsqregression}
\alias{residuals.lsqregression}
\alias{se.lsqregression}
\alias{gradients.lsqregression}
\alias{summary.lsqregression}
\alias{print.lsqregression}
\alias{plot.lsqregression}
\title{Location-Scale Kernel Quantile Regression with Mixed Data Types}
\description{
  \code{nplsqreg} computes a location-scale kernel estimate of the
  conditional quantile function for a one dimensional dependent variable
  and mixed continuous, unordered factor, and ordered factor explanatory
  data. Unlike \code{\link{npqreg}}, which obtains conditional quantiles
  by inverting an estimated conditional distribution, \code{nplsqreg}
  estimates the requested conditional quantile surface directly using a
  locally weighted quantile-kernel construction.
}
\usage{
nplsqreg(bws, \dots)

\method{nplsqreg}{formula}(bws, data = NULL, newdata = NULL, tau = 0.5,
       gradients = FALSE, residuals = FALSE, subset, na.action, \dots)

\method{nplsqreg}{lsqregressionbandwidth}(bws, txdat = NULL, tydat = NULL,
       tau = bws$tau, \dots)

\method{nplsqreg}{default}(bws,
       txdat = stop("training data 'txdat' missing"),
       tydat = stop("training data 'tydat' missing"),
       tau = 0.5,
       exdat,
       gradients = FALSE,
       residuals = FALSE,
       \dots)
}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification, formula/data interface,
and training data.
}

  \item{bws}{
    a formula, an \code{lsqregressionbandwidth} object returned by
    \code{\link{nplsqregbw}}, an \code{rbandwidth} object, a numeric
    bandwidth vector, or omitted for automatic bandwidth selection.
    Exact nplsqreg reuse is through the fitted object's \code{$bws}
    component; \code{$reg.bws} is internal regression state.
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

  \item{txdat}{
    a \eqn{p}-variate data frame of explanatory data used as training data.
    Defaults to the training data stored in \code{bws}.
  }

  \item{tydat}{
    a one dimensional numeric vector of dependent data. Defaults to the
    training response stored in \code{bws}.
  }

\subsection{Evaluation Data And Returned Quantities}{
These arguments control where the quantile regression is evaluated and which
fitted quantities are returned.
}

  \item{newdata}{
    an optional data frame in which to look for evaluation covariates for
    formula fits. If omitted, the training data are used.
  }

  \item{exdat}{
    a \eqn{p}-variate data frame of evaluation points. By default,
    evaluation takes place on \code{txdat}. The native \code{exdat}
    argument takes precedence over \code{newdata} when both are supplied.
  }

  \item{gradients}{
    a logical value indicating whether gradients and categorical effects
    of the conditional quantile with respect to the conditioning variables
    should be computed and returned. Defaults to \code{FALSE}.
  }

  \item{residuals}{
    a logical value indicating whether residuals should be returned for
    training-data fits. Defaults to \code{FALSE}.
  }

\subsection{Quantile Index And Additional Controls}{
These arguments control the requested quantile probability and the bandwidth
selection, prediction, or plotting route.
}

  \item{tau}{
    a numeric scalar or vector specifying the quantile probability or
    probabilities \eqn{\tau}{tau}. Values must lie strictly in \eqn{(0,1)}.
  }

  \item{\dots}{
    additional arguments supplied to \code{\link{nplsqregbw}} when
    bandwidths are computed internally, to \code{\link{npreg}} for the
    final transformed-response fit, or to plotting and prediction methods
    as appropriate. Common bandwidth-selection controls include
    \code{regtype}, \code{bwtype}, \code{nmulti}, \code{degree},
    \code{nomad}, \code{search.engine}, \code{tau.search},
    \code{delta}, \code{scale}, \code{regtype.pilot},
    \code{nomad.pilot}, and \code{pilot.args}.
  }
}
\details{
  The estimator follows the locally weighted quantile-kernel approach of
  Racine and Li (2017). Given a conditional scale pilot
  \eqn{\hat\sigma(X_i)}, define
  \deqn{
    Y_i^\delta = Y_i + \hat\sigma(X_i)\Phi^{-1}(\delta),
  }
  where \eqn{0 < \delta < 1} and \eqn{\Phi^{-1}} is the standard normal
  quantile function. For a requested quantile probability
  \eqn{\tau}{tau}, \code{nplsqregbw} selects the bandwidths and
  \eqn{\delta}{delta} by leave-one-out check-loss cross-validation. With
  the selected bandwidths and \eqn{\delta}{delta}, \code{nplsqreg} then
  fits a kernel regression of \eqn{Y_i^\delta} on \eqn{X_i} using the
  ordinary mixed-data machinery in \code{\link{npreg}}. The fitted mean
  of the transformed response is the estimated conditional quantile.

  The scale pilot is interpreted as a conditional standard deviation. The
  default pilot estimates the conditional mean, smooths squared residuals,
  floors the fitted variance before taking square roots, and then uses the
  resulting positive scale vector in the quantile-kernel transformation.
  The local-linear residual-variance pilot follows the idea of Fan and
  Yao (1998); \code{regtype.pilot} can be used to select the pilot
  regression type independently of the final quantile fit.

  If \code{tau} has length greater than one, \code{tau.search="full"}
  performs a separate bandwidth/\eqn{\delta}{delta} search for each
  quantile while sharing the same pilot scale. The explicit
  \code{tau.search="refined"} route fits the central quantile first and
  warm-starts the remaining quantiles, recording the search order and
  warm-start provenance in the returned object. The conservative default
  is \code{tau.search="full"}.

  Gradients and categorical effects are those returned by the final
  \code{\link{npreg}} fit on the transformed response. Thus ordered-factor
  effects are finite-difference contrasts and unordered-factor effects
  follow the corresponding mixed-data regression semantics.

  Standard errors are the asymptotic standard errors from that final
  transformed-response regression, conditional on the selected
  bandwidths, pilot scale, and \eqn{\delta}{delta}. For vector
  \code{tau}, per-quantile fits are independent, so estimated quantile
  curves are not monotonized and may cross in finite samples.
}
\value{
  \code{nplsqreg} returns an object of class \code{lsqregression}. The
  generic functions \code{\link{fitted}}, \code{\link{quantile}},
  \code{\link{se}}, \code{\link{predict}}, \code{\link{residuals}}, and
  \code{\link{gradients}} extract estimated conditional quantiles,
  asymptotic standard errors from the transformed-response regression,
  predictions, residuals when requested, and gradients or categorical
  effects. The functions \code{\link{summary}},
  \code{\link{print}}, and \code{\link{plot}} support objects of this type.
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
\section{Usage Issues}{
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same type.
}
\seealso{
  \code{\link{nplsqregbw}}, \code{\link{npreg}}, \code{\link{npqreg}},
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot.np}}
}
\examples{
\dontrun{
data("Italy")

model.q <- nplsqreg(gdp ~ ordered(year), data = Italy,
                    tau = c(0.25, 0.50, 0.75))
plot(model.q)

model.med <- nplsqreg(gdp ~ ordered(year), data = Italy, tau = 0.50)
model.q2 <- nplsqreg(bws = model.med$bws, tau = 0.50)
plot(model.med, gradient = TRUE)
}
}
\keyword{nonparametric}
\keyword{regression}
