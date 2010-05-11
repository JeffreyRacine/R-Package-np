% $ID$
\name{npindexbw}
\alias{npindexbw}
\alias{npindexbw.NULL}
\alias{npindexbw.default}
\alias{npindexbw.formula}
\alias{npindexbw.sibandwidth}
\title{Semiparametric Single Index Model Parameter and Bandwidth Selection}

\description{
  \code{npindexbw} computes a \code{npindexbw} bandwidth
  specification using the model \eqn{Y = G(X\beta) + \epsilon}{Y = G(XB)
    + epsilon}. The bandwidth object contains parameters for the single
  index model and the (scalar) bandwidth for the index function.
}

\usage{
npindexbw(\dots)

\method{npindexbw}{formula}(formula, data, subset, na.action, call, \dots)

\method{npindexbw}{NULL}(xdat = stop("training data xdat missing"),
          ydat = stop("training data ydat missing"),
          bws,
          \dots)

\method{npindexbw}{default}(xdat = stop("training data xdat missing"),
          ydat = stop("training data ydat missing"),
          bws,
          bandwidth.compute = TRUE,
          nmulti,
          random.seed,
          optim.method,
          optim.maxattempts,
          optim.reltol,
          optim.abstol,
          optim.maxit,
          \dots)

\method{npindexbw}{sibandwidth}(xdat = stop("training data xdat missing"),
          ydat = stop("training data ydat missing"),
          bws,
          bandwidth.compute = TRUE,
          nmulti,
          random.seed = 42,
          optim.method = c("Nelder-Mead", "BFGS", "CG"),
          optim.maxattempts = 10,
          optim.reltol = sqrt(.Machine$double.eps),
          optim.abstol = .Machine$double.eps,
          optim.maxit = 500,
          \dots)

}

\arguments{
  \item{formula}{
    a symbolic description of variables on which bandwidth selection is
    to be performed. The details of constructing a formula are
    described below.
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing
    the variables 
    in the model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which the
    function is called.
  }

  \item{subset}{
    an optional vector specifying a subset of observations to be used in
    the fitting process. 
  }

  \item{na.action}{
    a function which indicates what should happen when the data contain
    \code{NA}s. The default is set by the
    \code{\link{na.action}} setting of options, and is
    \code{\link{na.fail}} if that is unset. The (recommended) default is
    \code{\link{na.omit}}.  
  }
  
  \item{call}{
    the original function call. This is passed internally by
    \code{\link{np}} when a bandwidth search has been implied by a call to
    another function. It is not recommended that the user set this.
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

  \item{method}{
    the single index model method, one of either \dQuote{ichimura}
    (default) (Ichimura (1993)) or \dQuote{kleinspady} (Klein and Spady
    (1993)). Defaults to
    \code{ichimura}.
  }

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points. Defaults to \code{min(5,ncol(xdat))}.
  }

  \item{random.seed}{
    an integer used to seed R's random number generator. This ensures
    replicability of the numerical search. Defaults to 42.
  }

  \item{bandwidth.compute}{
    a logical value which specifies whether to do a numerical search for
    bandwidths or not. If set to \code{FALSE}, a \code{bandwidth} object
    will be returned with bandwidths set to those specified
    in \code{bws}. Defaults to \code{TRUE}.
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

  \item{optim.maxattempts}{
    maximum number of attempts taken trying to achieve successful
    convergence in \code{\link{optim}}. Defaults to \code{100}.
  }

  \item{optim.abstol}{
    the absolute convergence tolerance used by \code{\link{optim}}. Only useful
    for non-negative functions, as a tolerance for reaching
    zero. Defaults to \code{.Machine$double.eps}.
  }

  \item{optim.reltol}{
    relative convergence tolerance used by \code{\link{optim}}.  The algorithm
    stops if it is unable to reduce the value by a factor of 'reltol *
    (abs(val) + reltol)' at a step.  Defaults to
    \code{sqrt(.Machine$double.eps)}, typically about \code{1e-8}.
  }

  \item{optim.maxit}{
    maximum number of iterations used by \code{\link{optim}}. Defaults
     to \code{500}.
   }

  \item{\dots}{
    additional arguments supplied to specify the parameters to the
    \code{sibandwidth} S3 method, which is called during the numerical
    search. 
  }

  

}

\details{
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
  
  Note that, unlike most other bandwidth methods in the \code{np}
  package, this implementation uses the R \code{\link{optim}} nonlinear
  minimization routines and \code{\link{npksum}}. We have implemented
  multistarting and strongly encourage its use in practice. For
  exploratory purposes, you may wish to override the default search
  tolerances, say, setting \code{optim.reltol=.1} and conduct
  multistarting (the default is to restart min(5, ncol(xdat)) times) as is done
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

  Ichimura, H., (1993), \dQuote{Semiparametric least squares (SLS) and
  weighted SLS estimation of single-index models,} Journal of
  Econometrics, 58, 71-120.

  Klein, R. W. and R. H. Spady (1993), \dQuote{An efficient semiparametric
  estimator for binary response models,} Econometrica, 61, 387-421.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth
    estimators for discrete distributions,} Biometrika, 68, 301-309.
}

    

\author{
  Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
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
  restart min(5, ncol(xdat)) times). Once the procedure terminates, you can
  restart search with default tolerances using those bandwidths obtained
  from the less rigorous search (i.e., set \code{bws=bw} on subsequent
  calls to this routine where \code{bw} is the initial bandwidth
  object).  A version of this package using the \code{Rmpi} wrapper is
  under development that allows one to deploy this software in a
  clustered computing environment to facilitate computation involving
  large datasets.  }

\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): Generate a simple linear model then
# compute coefficients and the bandwidth using Ichimura's nonlinear
# least squares approach.

set.seed(12345)

n <- 100

x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)

y <- x1 - x2 + rnorm(n)

# Note - this may take a minute or two depending on the speed of your
# computer. Note also that the first element of the vector beta is
# normalized to one for identification purposes, and that X must contain
# at least one continuous variable.

bw <- npindexbw(formula=y~x1+x2, method="ichimura")

summary(bw)

\dontrun{

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# EXAMPLE 1 (INTERFACE=DATA FRAME): Generate a simple linear model then
# compute coefficients and the bandwidth using Ichimura's nonlinear
# least squares approach.

set.seed(12345)

n <- 100

x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)

y <- x1 - x2 + rnorm(n)

X <- cbind(x1, x2)

# Note - this may take a minute or two depending on the speed of your
# computer. Note also that the first element of the vector beta is
# normalized to one for identification purposes, and that X must contain
# at least one continuous variable.

bw <- npindexbw(xdat=X, ydat=y, method="ichimura")

summary(bw)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# EXAMPLE 2 (INTERFACE=DATA FRAME): Generate a simple binary outcome
# model then compute coefficients and the bandwidth using Klein and
# Spady's likelihood-based approach.

n <- 100

x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)

y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)

# Note that the first element of the vector beta is normalized to one
# for identification purposes, and that X must contain at least one
# continuous variable.

bw <- npindexbw(formula=y~x1+x2, method="kleinspady")

summary(bw)

# EXAMPLE 2 (INTERFACE=DATA FRAME): Generate a simple binary outcome
# model then compute coefficients and the bandwidth using Klein and
# Spady's likelihood-based approach.

n <- 100

x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)

y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)

X <- cbind(x1, x2)

# Note that the first element of the vector beta is normalized to one
# for identification purposes, and that X must contain at least one
# continuous variable.

bw <- npindexbw(xdat=X, ydat=y, method="kleinspady")

summary(bw)
} % enddontrun
}
\keyword{nonparametric}
