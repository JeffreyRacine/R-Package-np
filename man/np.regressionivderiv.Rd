\name{npregivderiv}
\alias{npregivderiv}

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
npregivderiv(y,
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
             ...)
}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the response, endogenous variables, instruments, exogenous covariates, and evaluation data.
}

  \item{w}{
    a \eqn{q}-variate data frame of instruments. The data types may be
    continuous, discrete (unordered and ordered factors), or some
    combination thereof.
  }

  \item{weval}{
    a \eqn{q}-variate data frame of instruments on which the regression
    will be estimated (evaluation data). By default, evaluation
    takes place on the data provided by \code{w}.
  }

  \item{x}{
    an \eqn{r}-variate data frame of exogenous regressors. The data
    types may be continuous, discrete (unordered and ordered factors),
    or some combination thereof.
  }

  \item{xeval}{
    an \eqn{r}-variate data frame of exogenous regressors on which the
    regression will be estimated (evaluation data). By default,
    evaluation takes place on the data provided by \code{x}.
  }


  \item{y}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{z}.
  }


  \item{z}{
    a \eqn{p}-variate data frame of endogenous regressors. The data
    types may be continuous, discrete (unordered and ordered factors),
    or some combination thereof.
  }

  \item{zeval}{
    a \eqn{p}-variate data frame of endogenous regressors on which the
    regression will be estimated (evaluation data). By default, evaluation
    takes place on the data provided by \code{z}.
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
  \eqn{E(Y|z)} (default, \code{"Eyz"}) or from \eqn{E(E(Y|z)|z)} (this can
  be overridden by providing \code{starting.values} below)
}  

\item{starting.values}{ an optional vector supplying the derivative
  immediately before the first reported Landweber-Fridman update. When
  omitted, the initialization derivative is obtained from \code{start.from}
  (see details below).
}

  \item{stop.on.increase}{
    a logical value (defaults to \code{TRUE}) indicating whether to halt
    iteration if the stopping criterion (see below) increases over the
    course of one iteration (i.e. it may be above the iteration tolerance
    but increased).
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
    this internal normalization.
  }

}
\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, and \code{\link{plot}} for plotting options.


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

  Iteration will terminate when either the change in the value of
  \eqn{||(E(y|w)-E(\varphi_k(z,x)|w))/E(y|w)||^2
  }{||(E(y|w)-E(phi(z,x)|w))/E(y|w)||^2 } from iteration to iteration is
  less than \code{iterate.diff.tol} or we hit \code{iterate.max} or
  \eqn{||(E(y|w)-E(\varphi_k(z,x)|w))/E(y|w)||^2
  }{||(E(y|w)-E(phi(z,x)|w))/E(y|w)||^2 } stops falling in value and
  starts rising.  }

\value{

  \code{npregivderiv} returns a \code{npregivderiv} object.  The
    generic functions \code{\link{print}}, \code{\link{summary}}, and
    \code{\link{plot}} support objects of this type.

  \code{npregivderiv} returns a list with components \code{phi.prime},
  \code{phi}, \code{phi.prime.mat}, \code{phi.mat},
  \code{num.iterations}, \code{norm.stop}, \code{convergence}, and the
  initialization components \code{starting.values.phi.prime} and
  \code{starting.values.phi}. \code{num.iterations} is the selected state
  returned in \code{phi.prime} and \code{phi}; the number of states evaluated,
  including any stopping-rule overshoot, is \code{length(norm.stop)}.
  
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

set.seed(42)
n <- 500

## For trimming the plot (trim .5\% from each tail)

trim <- 0.005

## The DGP is as follows:

## 1) y = phi(z) + u

## 2) E(u|z) != 0 (endogeneity present)

## 3) Suppose there exists an instrument w such that z = f(w) + v and
## E(u|w) = 0

## 4) We generate v, w, and generate u such that u and z are
## correlated. To achieve this we express u as a function of v (i.e. u =
## gamma v + eps)

v <- rnorm(n,mean=0,sd=0.27)
eps <- rnorm(n,mean=0,sd=0.05)
u <- -0.5*v + eps
w <- rnorm(n,mean=0,sd=1)

## In Darolles et al (2011) there exist two DGPs. The first is
## phi(z)=z^2 and the second is phi(z)=exp(-abs(z)) (which is
## discontinuous and has a kink at zero).

fun1 <- function(z) { z^2 }
fun2 <- function(z) { exp(-abs(z)) }

z <- 0.2*w + v

## Generate two y vectors for each function.

y1 <- fun1(z) + u
y2 <- fun2(z) + u

## You set y to be either y1 or y2 (ditto for phi) depending on which
## DGP you are considering:

y <- y1
phi <- fun1

## Sort on z (for plotting)

ivdata <- data.frame(y,z,w,u,v)
ivdata <- ivdata[order(ivdata$z),]
rm(y,z,w,u,v)
with(ivdata, {

model.ivderiv <- npregivderiv(y=y,z=z,w=w)

ylim <-c(quantile(model.ivderiv$phi.prime,trim),
         quantile(model.ivderiv$phi.prime,1-trim))

plot(z,model.ivderiv$phi.prime,
     xlim=quantile(z,c(trim,1-trim)),
     main="",
     ylim=ylim,
     xlab="Z",
     ylab="Derivative",
     type="l",
     lwd=2)
rug(z)     
})
} % enddontrun
}
\keyword{ instrument }
