%
\name{npqreg}
\alias{npqreg}
\alias{npqreg.condbandwidth}
\alias{npqreg.default}
\alias{npqreg.formula}
\alias{predict.qregression}
\alias{plot.qregression}
\title{Kernel Quantile Regression with Mixed Data Types}

\description{
  \code{npqreg} computes a kernel quantile regression estimate of a one
  (1) dimensional dependent variable on \eqn{p}-variate explanatory
  data, given a set of evaluation points, training points (consisting of
  explanatory data and dependent data), and a bandwidth specification
  using the methods of Li and Racine (2008) and Li, Lin and Racine
  (2013). A bandwidth specification can be a \code{condbandwidth} object,
  or a bandwidth vector, bandwidth type and kernel type.
}

\usage{
npqreg(bws, \dots)

\method{npqreg}{formula}(bws, data = NULL, newdata = NULL, \dots)


\method{npqreg}{condbandwidth}(bws,
       txdat = stop("training data 'txdat' missing"),
       tydat = stop("training data 'tydat' missing"),
       exdat,
       tau = 0.5,
       gradients = FALSE,
       tol = 1.490116e-04,
       small = 1.490116e-05,
       itmax = 10000,
       \dots)

\method{npqreg}{default}(bws, txdat, tydat, nomad = FALSE, \dots)

\method{predict}{qregression}(object, se.fit = FALSE, \dots)

\method{plot}{qregression}(x, \dots)

}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification, formula/data interface, and training data.
}

  \item{bws}{
    a bandwidth specification. This can be set as a \code{condbandwidth}
    object returned from an invocation of \code{\link{npcdistbw}}, or
    as a vector of bandwidths, with each element \eqn{i} corresponding
    to the bandwidth for column \eqn{i} in \code{txdat}. If specified as
    a vector, then additional arguments will need to be supplied as
    necessary to specify the bandwidth type, kernel types, and so on.
  }


  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npcdistbw}} was called.
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

  \item{object}{
    an object of class \code{"qregression"} returned by \code{npqreg}.
  }

  \item{x}{
    an object of class \code{"qregression"} returned by \code{npqreg}.
  }


\subsection{Local-Polynomial Degree And Bandwidth Search}{
This argument controls the recommended automatic local-polynomial NOMAD route, which jointly selects continuous polynomial degree and bandwidths when conditional-distribution bandwidths are computed inside \code{npqreg}.
}

  \item{nomad}{
    logical shortcut passed through to \code{\link{npcdistbw}} when
    bandwidths are computed inside \code{npqreg}. When \code{TRUE}, the
    conditional-distribution bandwidth route fills any missing values
    among \code{regtype}, \code{search.engine}, \code{degree.select},
    \code{bernstein.basis}, \code{degree.min}, \code{degree.max},
    \code{degree.verify}, and \code{bwtype} with the recommended
    automatic local-polynomial degree-and-bandwidth NOMAD preset documented in \code{\link{npcdistbw}}.
    Additional NOMAD tuning arguments such as \code{nomad.nmulti} may
    also be supplied through \code{\dots}; \code{nmulti} remains the
    outer restart count while \code{nomad.nmulti} controls inner
    \code{crs::snomadr()} multistarts within each outer restart.
    After fitting, inspect \code{fit$bws$nomad.shortcut} on the
    returned object \code{fit} to see the normalized shortcut metadata.
  }


\subsection{Evaluation Data And Returned Quantities}{
These arguments control where the quantile regression is evaluated and which fitted quantities are returned.
}

  \item{exdat}{
    a \eqn{p}-variate data frame of points on which the regression will be
    estimated (evaluation data). By default,
    evaluation takes place on the data provided by \code{txdat}.
  }


  \item{gradients}{
    a logical value indicating that you want gradients of the conditional
    quantile with respect to the conditioning variables computed and returned
    in the resulting \code{qregression} object. Defaults to \code{FALSE}.
  }


  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.
  }

  \item{se.fit}{
    logical value. If \code{TRUE}, \code{predict.qregression} returns a
    list with components \code{fit} and \code{se.fit}; otherwise it
    returns fitted conditional quantiles.
  }


  \item{tau}{
    a numeric scalar or vector specifying the quantile probability or
    probabilities \eqn{\tau}{tau}. Defaults to \code{0.5}.
  }



\subsection{Quantile Solver Controls}{
These arguments control the one-dimensional numerical quantile extraction step.
}

  \item{itmax}{
    integer maximum number of iterations allowed in the one-dimensional
    quantile refinement. Defaults to \code{10000}.
  }


  \item{small}{
    minimum interval width used by the one-dimensional quantile
    refinement. Defaults to \code{1.490116e-05}
    (approximately \code{1000*sqrt(.Machine$double.eps)}).
  }


  \item{tol}{
    tolerance on the one-dimensional quantile location refinement.
    Defaults to \code{1.490116e-04}
    (approximately \code{10000*sqrt(.Machine$double.eps)}).
  }


\subsection{Additional Arguments}{
Further arguments are passed to the bandwidth-selection counterpart,
prediction/evaluation route, or plot route as appropriate.
}

  \item{\dots}{
    additional arguments supplied to \code{\link{npcdistbw}} when
    \code{npqreg} computes bandwidths internally, or arguments needed
    to interpret a numeric \code{bws} vector. This is where bandwidth
    selection controls such as \code{bwmethod}, \code{bwtype}, and
    \code{bwscaling}, kernel/support controls such as
    \code{cxkertype}, \code{cykertype}, \code{cxkerorder},
    \code{cykerorder}, \code{cxkerbound}, and \code{cykerbound},
    categorical kernel controls such as \code{uxkertype},
    \code{uykertype}, \code{oxkertype}, and \code{oykertype}, search
    controls such as \code{nmulti} and
    \code{scale.factor.search.lower}, and local-polynomial/NOMAD
    controls such as \code{regtype}, \code{degree},
    \code{bernstein.basis}, \code{degree.select}, and
    \code{nomad.nmulti} are supplied. In \code{predict.qregression},
    additional arguments are passed to \code{npqreg} for evaluation
    with the stored bandwidth object; common examples are
    \code{newdata}, native \code{exdat}, and \code{tau}. In
    \code{plot.qregression}, additional arguments are passed through
    the package plot route; common controls include \code{tau},
    \code{gradients}, \code{output}, \code{legend}, and graphics
    arguments. See \code{\link{npcdistbw}} and \code{\link{plot.np}}
    for the complete bandwidth-selection and plotting argument
    surfaces.
  }
}
\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, and \code{\link{plot}}, \code{\link{plot.np}} for plotting options.


  Given a conditional distribution bandwidth object, \code{npqreg}
  estimates the conditional distribution function
  \eqn{F(y|x)} and extracts the requested conditional quantile. For
  \eqn{0 < \tau < 1}, the conditional quantile at probability
  \eqn{\tau}{tau} is
  \deqn{
    q_\tau(x) = \inf\{y : F(y|x) \ge \tau\}.
  }
  Equivalently, \eqn{q_\tau(x)} is a quasi-inverse of the conditional
  distribution in the sense of Nelsen (2006): an inverse agrees with
  \eqn{F} on the range of \eqn{F}, while outside that range the
  generalized inverse is defined by the lower endpoint at which
  \eqn{F} reaches or exceeds the requested probability. Numerically,
  \code{npqreg} inverts the selected conditional distribution estimator
  represented by \code{bws}. This includes the selected bandwidth type,
  kernels, local-polynomial regression type, selected polynomial degree,
  basis, and Bernstein-basis setting inherited from \code{\link{npcdistbw}}.
  If the bandwidth object was selected with \code{nomad=TRUE}, the returned
  conditional-distribution bandwidth object is an LP object: its
  \code{regtype}/\code{regtype.engine} metadata identify the selected
  local-polynomial route and its \code{degree}/\code{degree.engine} metadata
  record the selected continuous-coordinate polynomial degree. \code{npqreg},
  \code{\link{predict}}, and \code{\link{plot}} reuse this stored LP metadata;
  plotting additional \code{tau} values does not recompute or downgrade the
  selected degree.
  The one-dimensional inversion is carried out over the observed support
  of the dependent variable using the same selected conditional CDF
  estimator that is later used for quantile standard errors and gradients.
  The arguments \code{tol}, \code{small}, and \code{itmax} control this
  one-dimensional refinement.

  Let \eqn{f(y|x) = \partial F(y|x)/\partial y} denote the conditional
  density. The asymptotic standard error of the conditional quantile is
  computed by the first-order delta method,
  \deqn{
    se\{\hat q_\tau(x)\}
      =
    \frac{se\{\hat F(\hat q_\tau(x)|x)\}}
         {\hat f(\hat q_\tau(x)|x)} ,
  }
  using the selected conditional distribution standard-error machinery
  and the selected conditional density evaluated at the fitted quantile.
  This corresponds to the quantile variance expression in Li, Lin and
  Racine (2013).

  If \code{gradients=TRUE}, \code{npqreg} also computes gradients of
  the conditional quantile with respect to the conditioning variables
  for which gradients are defined. Differentiating
  \eqn{F(q_\tau(x)|x) = \tau} gives
  \deqn{
    \nabla_x q_\tau(x)
      =
    -\frac{F_x(q_\tau(x)|x)}{f(q_\tau(x)|x)},
  }
  where \eqn{F_x(y|x)} is the derivative of the same selected conditional
  distribution estimator with respect to \eqn{x}. For \code{regtype="lc"},
  this uses the local-constant conditional-gradient machinery; for
  \code{regtype="ll"} it uses the canonical local-polynomial degree-one
  route; and for \code{regtype="lp"} it uses the selected or supplied
  degree vector. The corresponding first-order gradient standard errors
  are computed componentwise as
  \deqn{
    se\{\nabla_x \hat q_\tau(x)\}
      =
    \frac{se\{\hat F_x(\hat q_\tau(x)|x)\}}
         {\hat f(\hat q_\tau(x)|x)} .
  }

  When \code{npqreg} is called without an explicit \code{bws} object, it
  first computes conditional distribution bandwidths using
  \code{\link{npcdistbw}} and stores them in the returned object's
  \code{bws} component. If a scalar \code{tau} was used initially and
  additional quantiles are later desired as fitted objects, reuse those
  selected bandwidths directly, for example
  \code{npqreg(bws = fit$bws, tau = c(0.25, 0.5, 0.75))}. If the goal
  is only to inspect additional quantiles graphically, use
  \code{plot(fit, tau = c(0.25, 0.5, 0.75))}; this reuses the stored
  bandwidths and recomputes only the one-dimensional quantile extraction
  step for the requested \code{tau} values. Vector-\code{tau} plots are
  overlaid and include a legend; use \code{legend=FALSE},
  \code{legend=NULL}, or a \code{legend=list(...)} control to suppress
  or customize it.

  The \code{\link{predict}} method follows the usual S3
  \code{newdata} convention. For formula fits, supply a data frame of
  evaluation covariates via \code{predict(fit, newdata=...)}. For
  non-formula fits, \code{newdata} is translated to the native
  evaluator argument \code{exdat} when \code{exdat} is not supplied.
  The native \code{exdat} argument remains available for advanced
  workflows and takes precedence if both \code{newdata} and
  \code{exdat} are supplied. If \code{tau} is omitted in
  \code{predict}, the fitted object's stored \code{tau} value is used.

}

\value{

  \code{npqreg} returns a \code{qregression} object.  The generic
  functions \code{\link{fitted}} (or \code{\link{quantile}}),
  \code{\link{se}}, \code{\link{predict}}, and \code{\link{gradients}}
  extract (or generate) estimated values, asymptotic standard errors on
  estimates, predictions, and gradients, respectively, from the returned
  object. \code{\link{predict}} uses the object's stored \code{tau}
  value by default; supply \code{tau=} to override it. Furthermore, the functions
  \code{\link{summary}} and \code{\link{plot}} support objects of this
  type. The returned object has the following components:

  \item{eval}{ evaluation points }
  \item{quantile}{ estimation of the quantile regression function
    (conditional quantile) at the evaluation points. If \code{tau} has
    length greater than one this is an evaluation-by-\code{tau} matrix. }
  \item{quanterr}{ asymptotic standard errors of the quantile
    regression estimates, obtained from the conditional distribution
    standard error and the estimated conditional density at the fitted
    quantile. If \code{tau} has length greater than one this is an
    evaluation-by-\code{tau} matrix. }
  \item{quantgrad}{ gradients of the conditional quantile with respect
    to the conditioning variables at each evaluation point, when
    \code{gradients=TRUE}. If \code{tau} has length greater than one
    this is an evaluation-by-gradient-by-\code{tau} array. }
  \item{quantgerr}{ asymptotic standard errors for gradients, when
    \code{gradients=TRUE}. If \code{tau} has length greater than one
    this is an evaluation-by-gradient-by-\code{tau} array. }
  \item{tau}{ the quantile probability or probabilities computed }


}

\section{Book And Method Pointers}{
  The conditional quantile target is the generalized inverse
  \eqn{q_\tau(x)=\inf\{y:F(y\mid x)\ge \tau\}} of the conditional
  distribution. The standard errors and gradients described above are
  first-order delta-method quantities evaluated using the same selected
  conditional CDF, conditional density, bandwidths, kernels, and
  local-polynomial degree inherited from the supplied
  \code{\link{npcdistbw}} object.

  For book-length derivations, see Li and Racine (2007), Chapter 6
  \emph{Conditional CDF and Quantile Estimation}, especially Sections
  6.3-6.5, and Racine (2019), Chapter 4 \emph{Conditional Probability
  Density and Cumulative Distribution Functions}. The quasi-inverse
  terminology follows Nelsen (2006).
}

\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Hall, P. and J.S. Racine and Q. Li (2004), \dQuote{Cross-validation
    and the estimation of conditional probability densities,} Journal of
    the American Statistical Association, 99, 1015-1026.

  Koenker, R. W. and G.W. Bassett (1978), \dQuote{Regression
  quantiles,} Econometrica, 46, 33-50.

  Koenker, R. (2005), \emph{Quantile Regression,} Econometric Society
  Monograph Series, Cambridge University Press.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2008), \dQuote{Nonparametric estimation of
  conditional CDF and quantile functions with mixed categorical and
  continuous data,} Journal of Business and Economic Statistics, 26,
  423-434.

  Li, Q. and J. Lin and J.S. Racine (2013), \dQuote{Optimal Bandwidth
  Selection for Nonparametric Conditional Distribution and Quantile
  Functions}, Journal of Business and Economic Statistics, 31, 57-65.

  Nelsen, R.B. (2006), \emph{An Introduction to Copulas}, Second
  Edition, Springer.

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
}

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}, \code{\link{plot.np}}, \pkg{quantreg} }

\examples{
\dontrun{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we compute a
# bivariate nonparametric quantile regression estimate for Giovanni
# Baiocchi's Italian income panel (see Italy for details)

data("Italy")
with(Italy, {

# Compute conditional distribution bandwidths and extract three
# conditional quantiles using the same selected bandwidths.

model.q <- npqreg(gdp~ordered(year), tau=c(0.25, 0.50, 0.75))

# Plot the overlaid quantiles.

plot(model.q)

# If a scalar tau was used first, additional quantiles can reuse the
# selected bandwidths without recomputing cross-validation. Use npqreg()
# when the additional fitted values are needed as an object, or plot()
# when graphical inspection is all that is desired.

model.med <- npqreg(gdp~ordered(year), tau=0.50)
model.q <- npqreg(bws=model.med$bws, tau=c(0.25, 0.50, 0.75))
plot(model.med, tau=c(0.25, 0.50, 0.75))

})

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we compute a
# bivariate nonparametric quantile regression estimate for Giovanni
# Baiocchi's Italian income panel (see Italy for details)

data("Italy")
with(Italy, {
data <- data.frame(ordered(year), gdp)

# First, compute the likelihood cross-validation bandwidths (default).
# Note - this may take a few minutes depending on the speed of your
# computer...

bw <- npcdistbw(xdat=ordered(year), ydat=gdp)

# Note - numerical search for computing the quantiles will take a
# minute or so...

model.q <- npqreg(bws=bw, tau=c(0.25, 0.50, 0.75))

plot(model.q)

})
} % enddontrun
}
\keyword{nonparametric}
