\name{npregivderiv}
\alias{npregivderiv}
\alias{npregivderiv.default}
\alias{npregivderiv.formula}
\alias{summary.npregivderiv}
\alias{print.summary.npregivderiv}
\alias{fitted.npregivderiv}
\alias{gradients.npregivderiv}
\alias{residuals.npregivderiv}
\alias{plot.npregivderiv}

\title{
Nonparametric Instrumental Derivatives
}
\description{
\code{npregivderiv} uses the approach of Florens, Racine and Centorrino
(2018) to compute the partial derivative of a nonparametric
estimation of an instrumental regression function \eqn{\varphi}{phi}
defined by conditional moment restrictions stemming from a structural
econometric model: \eqn{E [Y - \varphi (Z,X) | W ] = 0}{E [Y - phi (Z,X)
| W ] = 0}, and involving endogenous variables \eqn{Y} and \eqn{Z} and
exogenous variables \eqn{X} and instruments \eqn{W}. The derivative
function \eqn{\varphi'}{phi'} is the solution of an ill-posed inverse
problem, and is computed using Landweber-Fridman regularization.
}
\usage{
npregivderiv(y, \dots)

\method{npregivderiv}{default}(y,
             z,
             w,
             x = NULL,
             zeval = NULL,
             weval = NULL,
             xeval = NULL,
             constant = 0.5,
             iterate.break = TRUE,
             iterate.max = 1000,
             nmulti = NULL,
             random.seed = 42,
             smooth.residuals = TRUE,
             start.from = c("Eyz","EEywz"),
             starting.values = NULL,
             stop.on.increase = TRUE,
             regtype = c("ll","lc","lp"),
             degree = NULL,
             nomad = FALSE,
             ...)

\method{npregivderiv}{formula}(y,
             data = NULL,
             subset,
             na.action,
             newdata = NULL,
             zeval = NULL,
             weval = NULL,
             xeval = NULL,
             ...)

\method{plot}{npregivderiv}(x, plot.data = FALSE, phi = FALSE, ...)
}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the response, endogenous variables, instruments, exogenous covariates, and evaluation data.
}

  \item{data}{
    an optional data frame used with the formula interface. This estimator
    accepts \code{y ~ z | w}. The shared parser recognizes the optional
    \code{| x} partition but this univariate derivative operator rejects it
    explicitly; see Details.
  }

  \item{subset}{
    an optional expression selecting training rows in the formula interface.
  }

  \item{na.action}{
    a function specifying missing-value handling for the formula model frame.
  }

  \item{newdata}{
    an optional data frame supplying fit-time values of formula variable
    \code{z} at which to evaluate the structural derivative. It may contain
    any positive number of rows. Explicit \code{zeval} takes precedence.
  }

  \item{w}{
    a \eqn{q}-variate data frame of instruments. The data types may be
    continuous, discrete (unordered and ordered factors), or some
    combination thereof.
  }

  \item{weval}{
    retained for compatibility. The inverse problem is defined on the
    training instruments, so a supplied value must equal \code{w}; use
    \code{zeval} to request derivative evaluation.
  }

  \item{x}{
    for \code{npregivderiv.default}, reserved for interface compatibility. The present derivative estimator
    supports one univariate endogenous \code{z} and does not yet support a
    separate exogenous \code{x}; a non-\code{NULL} value fails explicitly. For
    \code{plot.npregivderiv}, an \code{npregivderiv} object.
  }

  \item{xeval}{
    reserved with \code{x}; it is not available until the derivative operator
    supports separate exogenous regressors.
  }


  \item{y}{
    a one (1) dimensional numeric or integer vector of dependent data, or a
    formula using the explicit IV grammar described under \code{data}. Each
    response element \eqn{i} corresponds to row \eqn{i} of \code{z}.
  }


  \item{z}{
    one finite continuous numeric endogenous regressor with at least two
    distinct support points.
  }

  \item{zeval}{
    finite values of the endogenous regressor at which to evaluate the
    structural derivative. The training inverse problem, residuals, centering,
    stopping rule, and selected state remain based on \code{z} and \code{w}.
    By default no separate evaluation surface is stored.
  }

  \item{plot.data}{
    logical; when \code{phi=TRUE}, add the selected training structural curve
    to a plot of the training response. It is ignored for derivative plots.
  }

  \item{phi}{
    logical; plot the reconstructed training structural function instead of
    the structural derivative. Derivative plots use \code{zeval} and
    \code{phi.prime.eval} when a separate evaluation grid was supplied.
  }

\subsection{Iteration Controls}{
These arguments control derivative iteration and reproducibility settings.
}

\item{constant}{
  the constant to use for Landweber-Fridman iteration.
}


  \item{iterate.break}{
    a logical value indicating whether to compute all objects up to
    \code{iterate.max} or to break when a potential optimum arises
    (useful for inspecting full stopping rule profile up to
    \code{iterate.max})
  }


  \item{iterate.max}{
    an integer indicating the maximum number of complete positive
    Landweber-Fridman states to evaluate. The initialization state is
    numbered zero and does not consume this budget.
  }


  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points.
  }


  \item{random.seed}{
    an integer used to seed R's random number generator. This ensures
    replicability of the numerical search. Defaults to 42.
  }


\subsection{Residual Smoothing And Starting Values}{
These arguments control residual smoothing and the initial derivative path.
}

\item{smooth.residuals}{
  a logical value (defaults to \code{TRUE}) indicating whether to
  optimize bandwidths for the regression of \eqn{y-\varphi(z)}{y-phi(z)}
  on \eqn{w} or for the regression of \eqn{\varphi(z)}{phi(z)} on
  \eqn{w} during Landweber-Fridman iteration.
 }


\item{start.from}{ a character string indicating whether to start from
  \eqn{E(Y|z)} (default, \code{"Eyz"}) or from \eqn{E(E(Y|w)|z)} (this can
  be overridden by providing \code{starting.values} below)
}  

\item{starting.values}{ an optional vector supplying the derivative
  immediately before the first reported Landweber-Fridman update. When
  omitted, the initialization derivative is obtained from \code{start.from}
  (see details below). A separate \code{zeval} grid is not accepted with a
  supplied initialization because no interpolation rule for the externally
  supplied derivative has been declared.
}

  \item{stop.on.increase}{
    a logical value (defaults to \code{TRUE}) indicating whether to halt
    iteration if the stopping criterion (see below) increases over the
    course of one iteration (i.e. it may be above the iteration tolerance
    but increased).
  }

\subsection{Regression Smoothing Controls}{
These arguments control all internal \code{npreg} stages consistently.
}

  \item{regtype}{
    a fixed regression type: \code{"ll"} (the default), \code{"lc"}, or
    \code{"lp"}. The default applies local-linear, degree-one smoothing to
    every internal stage containing continuous predictors, matching the
    default polynomial order of \code{npregiv}. A categorical-only internal
    stage uses local constant because no continuous polynomial term exists.
    Supply \code{regtype="lc"} to recover the former omitted-default
    computation.
  }

  \item{degree}{
    one nonnegative integer required with \code{regtype="lp"}. It is expanded
    only across the continuous coordinates of each internal regression stage;
    vectors are not recycled across heterogeneous stages, and a
    categorical-only stage uses local constant.
  }

  \item{nomad}{
    logical or \code{"auto"}. Only \code{FALSE} is currently supported for
    the multi-stage derivative estimator. Other values fail explicitly until
    the stage-wise search/state contract is validated.
  }


\subsection{Additional Arguments}{
Further arguments are passed to \code{\link{npreg}} and \code{\link{npksum}}.
}

  \item{...}{ 
    additional arguments supplied to \code{\link{npreg}} and
    \code{\link{npksum}}. The private integral kernel sum in the empirical
    adjoint always uses \code{bandwidth.divide = TRUE}, as required for the
    ordinary CDF normalization; a \code{bandwidth.divide} value in
    \code{...} continues to be supplied to \code{npreg} but does not override
    this internal normalization. Likewise, \code{ukertype} and
    \code{okertype} in \code{...} control applicable internal regressions,
    while the private adjoint retains its required
    \code{ukertype = "liracine"} and \code{okertype = "liracine"}. For the
    plot method, graphical arguments are passed to \code{plot} and
    \code{lines} as applicable.
  }

}
\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, \code{\link{plot}} for plotting options, and \code{\link{npRmpi.init}} for interactive/cluster MPI startup. See \code{\link{npRmpi.init}} details for performance tradeoffs (message passing/startup mode) and the \code{inst/Rprofile} manual-broadcast template.

  The formula interface uses explicit role order:
  \code{y ~ z | w} and \code{y ~ z | w | x}. Additive transformations and
  factors are retained through \code{model.frame}. Dot expansion and formula
  interaction operators are rejected because they do not identify a unique
  mixed-data kernel role; create an explicit transformed variable (for example
  with \code{I(a * b)}) instead. Native vector/data-frame calls remain
  supported unchanged.

  The current derivative operator is univariate in \code{z}. Consequently,
  although the shared parser recognizes the package-wide optional \code{| x}
  partition, \code{npregivderiv} rejects it explicitly rather than silently
  changing the derivative target. Post-fit \code{predict()} is not available;
  use fit-time \code{newdata} or the native evaluation arguments.

  Fit-time evaluation is output-only. Conditional regressions, the empirical
  adjoint, fitted-residual centering, structural integration, stopping values,
  and state selection are always computed on aligned training rows. A separate
  \code{zeval} accumulates the derivative corrections from those same training
  states into \code{phi.prime.eval} and \code{phi.prime.eval.mat}; it cannot
  change \code{phi}, \code{phi.prime}, their training paths, or
  \code{norm.stop}. Structural levels at arbitrary new points are not returned
  because this derivative estimator's integration and centering normalization
  is defined on the training support.


  Let \eqn{\varphi'_0}{phi'[0]} denote the initialization derivative,
  obtained from \code{start.from} or supplied directly through
  \code{starting.values}. Its integrated and centered curve is
  \eqn{\varphi_0}{phi[0]}. This state is retained in
  \code{starting.values.phi.prime} and \code{starting.values.phi}; it is
  not a column of the positive-iteration matrices. For \eqn{N=1,2,\ldots},
  iteration \eqn{N} means exactly \eqn{N} completed residual-correction
  updates after this initialization. Column \eqn{N} of
  \code{phi.prime.mat} and \code{phi.mat} therefore contains the coherent
  pair \eqn{(\varphi'_N,\varphi_N)}{(phi'[N],phi[N])}, and
  \code{norm.stop[N]} is evaluated at that same state using the multiplier
  \eqn{N}. Check the stopping vector for anomalies such as an immediate
  increase.
  
  Landweber-Fridman iteration uses an optimal stopping rule based upon
  \eqn{||E(y|w)-E(\varphi_k(z,x)|w)||^2 }{||E(y|w)-E(phi(z,x)|w)||^2
  }. However, if local rather than global optima are encountered the
  resulting estimates can be overly noisy. To best guard against this
  eventuality set \code{nmulti} to a larger number than the default
  \code{nmulti=min(2,p)} for the first iteration, where \code{p} is
  the dimension of the current smoothing problem.

  Note that for subsequent Landweber-Fridman iterations, a \dQuote{warm
  start} strategy is employed. The optimal bandwidths from the previous
  iteration are used as starting values for the current iteration. The
  user-supplied \code{nmulti} is respected for all iterations. For
  iterations after the first successful one, these optimal bandwidths
  serve as the first of the multiple initial points (a warm start),
  while any remaining restarts are cold starts. If \code{nmulti} is not
  explicitly supplied by the user, it defaults to \code{min(2,p)} for the first
  iteration and to 1 for all subsequent iterations. This strategy
  provides a balance between computational efficiency and robustness,
  allowing the numerical optimizer to refine the structural bandwidths
  as the residuals evolve incrementally while still guarding against
  local optima.

  With \code{iterate.break=TRUE}, iteration can stop after the stopping rule
  turns upward once the initial guard interval has passed; otherwise all
  states through \code{iterate.max} are evaluated. The selected state is then
  obtained from the computed stopping profile.  }

\value{

  \code{npregivderiv} returns a \code{npregivderiv} object.  The
    generic functions \code{\link{print}}, \code{\link{summary}},
    \code{\link{plot}}, \code{\link{fitted}}, \code{\link{residuals}}, and
    \code{\link{gradients}} support objects of this type. \code{fitted}
    returns the reconstructed training structural function while \code{gradients}
    returns the training derivative estimate. Residuals use the aligned
    training structural curve even when a separate evaluation grid was
    supplied. For formula fits using an exclusion-style
    \code{na.action}, these accessors restore omitted rows. \code{summary}
    returns a structured
    \code{summary.npregivderiv} object with a dedicated print method. The
    printed smoothing description reports the effective continuous-stage
    policy using the same \code{Kernel Regression Estimator:} labels as
    \code{npreg}: \code{Local-Constant}, \code{Local-Linear}, or, for
    \code{regtype="lp"}, the local-polynomial degree and basis details.

  \code{npregivderiv} returns a list with components \code{phi.prime},
  \code{phi}, \code{phi.prime.mat}, \code{phi.mat}, and, when a separate
  evaluation grid was supplied, \code{phi.prime.eval} and
  \code{phi.prime.eval.mat},
  \code{num.iterations}, \code{norm.stop}, \code{convergence}, and the
  initialization components \code{starting.values.phi.prime} and
  \code{starting.values.phi}. \code{num.iterations} is the selected state
  returned in \code{phi.prime} and \code{phi}; the number of states evaluated,
  including any stopping-rule overshoot, is \code{length(norm.stop)}. The
  compact \code{bws}, \code{smoothing.spec}, and \code{stage.specs}
  components record stage bandwidth payloads and requested/effective
  fixed-degree smoothing state without replacing the legacy bandwidth fields.
  
}

\references{
  
  Carrasco, M. and J.P. Florens and E. Renault (2007), \dQuote{Linear
  Inverse Problems in Structural Econometrics Estimation Based on
  Spectral Decomposition and Regularization,} In: James J. Heckman and
  Edward E. Leamer, Editor(s), Handbook of Econometrics, Elsevier, 2007,
  Volume 6, Part 2, Chapter 77, Pages 5633-5751

  Darolles, S. and Y. Fan and J.P. Florens and E. Renault (2011),
  \dQuote{Nonparametric instrumental regression,} Econometrica, 79,
  1541-1565.
  
  Feve, F. and J.P. Florens (2010), \dQuote{The practice of
  non-parametric estimation by solving inverse problems: the example of
  transformation models,} Econometrics Journal, 13, S1-S27.

  Florens, J.P. and J.S. Racine and S. Centorrino (2018),
  \dQuote{Nonparametric instrumental derivatives,} Journal of
  Nonparametric Statistics, 30 (2), 368-391.

  Fridman, V. M. (1956), \dQuote{A method of successive approximations
  for Fredholm integral equations of the first kind,} Uspeskhi,
  Math. Nauk., 11, 233-334, in Russian.

  Horowitz, J.L. (2011), \dQuote{Applied nonparametric instrumental
  variables estimation,} Econometrica, 79, 347-394.
  
  Landweber, L. (1951), \dQuote{An iterative formula for Fredholm
  integral equations of the first kind,} American Journal of
  Mathematics, 73, 615-24.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated Local Linear
  Nonparametric Regression,} Statistica Sinica, 14, 485-512.

}
\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\note{
  This function currently supports univariate \code{z} only.  This
  function should be considered to be in \sQuote{beta test} status until
  further notice.
}

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}
\code{\link{npregiv},\link{npreg}}
}
\examples{
\dontrun{
## This illustration was made possible by Samuele Centorrino
## <samuele.centorrino@univ-tlse1.fr>

npRmpi.init(nslaves = 1)

set.seed(42)
n <- 500
trim <- 0.005

## Endogenous-regressor DGP:
##
## y = phi(z) + u
## E(u | z) != 0
## z = 0.2 w + v
## E(u | w) = 0

v <- rnorm(n, mean = 0, sd = 0.27)
eps <- rnorm(n, mean = 0, sd = 0.05)
u <- -0.5 * v + eps
w <- rnorm(n)

phi <- function(z) z^2
phi.prime <- function(z) 2 * z

z <- 0.2 * w + v
y <- phi(z) + u

ivdata <- data.frame(y, z, w)
ivdata <- ivdata[order(ivdata$z), ]

model.ivderiv <- npregivderiv(y ~ z | w, data = ivdata)
deriv.iv <- gradients(model.ivderiv)

## Restrict the display to the central 99 percent of the z observations.

zlim <- quantile(ivdata$z, c(trim, 1 - trim))
keep <- ivdata$z >= zlim[1L] & ivdata$z <= zlim[2L]

plot(ivdata$z[keep], deriv.iv[keep],
     main = "",
     xlab = "Z",
     ylab = "Derivative",
     type = "l",
     lwd = 2)

curve(phi.prime(x), add = TRUE, col = 2, lty = 2, lwd = 2)
rug(ivdata$z[keep])

npRmpi.quit()
} % enddontrun
}
\keyword{ instrument }
