% $ID$
\name{npscoef}
\alias{npscoef}
\alias{npscoef.call}
\alias{npscoef.default}
\alias{npscoef.formula}
\alias{npscoef.scbandwidth}
\title{Smooth Coefficient Kernel Regression}

\description{
  \code{npscoef} computes a kernel regression estimate of a one (1)
  dimensional dependent variable on \eqn{p}-variate explanatory data,
  using the model \eqn{Y_i = W_{i}^{\prime} \gamma (Z_i) + u_i}{Y_i =
  t(W_i) * gamma(Z_i) + u_i} where
  \eqn{W_i'=(1,X_i')}{t(W_i)=(1,t(X_i))}, given a set of evaluation
  points, training points (consisting of explanatory data and dependent
  data), and a bandwidth specification. A bandwidth specification can be
  a \code{scbandwidth} object, or a bandwidth vector, bandwidth type and
  kernel type.
}

\usage{
npscoef(bws, \dots)

\method{npscoef}{formula}(bws, data = NULL, newdata = NULL, \dots)

\method{npscoef}{call}(bws, \dots)

\method{npscoef}{default}(bws, txdat, tydat, tzdat, \dots)

\method{npscoef}{scbandwidth}(bws,
        txdat = stop("training data 'txdat' missing"),
        tydat = stop("training data 'tydat' missing"),
        tzdat = NULL,
        exdat,
        eydat,
        ezdat,
        residuals = FALSE,
        errors = TRUE,
        iterate = TRUE,
        maxiter = 100,
        tol = .Machine$double.eps,
        leave.one.out = FALSE,
        betas = FALSE,
        \dots)
}

\arguments{
  \item{bws}{
  a bandwidth specification. This can be set as a \code{scbandwidth}
  object returned from an invocation of \code{\link{npscoefbw}}, or
    as a vector of bandwidths, with each element \eqn{i} corresponding
    to the bandwidth for column \eqn{i} in \code{tzdat}. If specified as
    a vector
    additional arguments will need to be supplied as necessary to
    specify the bandwidth type, kernel types, training data, and so on.

  }

  \item{\dots}{
    additional arguments supplied to specify the regression type,
    bandwidth type, kernel types, selection methods, and so on.
    To do this, you may specify any of \code{bwscaling}, \code{bwtype},
    \code{ckertype}, \code{ckerorder}, as described in
    \code{\link{npscoefbw}}. 
    
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing
    the variables 
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
   \code{\link{npscoefbw}} was called.
  }

  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.  
  }


  \item{txdat}{
    a \eqn{p}-variate data frame of explanatory data (training data),
    which, by default, populates the columns \eqn{2} through \eqn{p+1}
    of \eqn{W} in the model equation, and in the 
    absence of \code{zdat}, will also correspond to 
    \eqn{Z} from the model equation. Defaults to
    the training data used to 
    compute the bandwidth object.
  }


  \item{tydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{txdat}. Defaults to
    the training data used to 
    compute the bandwidth object.
  }

  \item{tzdat}{
    a optionally specified \eqn{q}-variate data frame of explanatory data
    (training data), which corresponds to \eqn{Z} in the model
    equation. Defaults to 
    the training data used to 
    compute the bandwidth object.
  }

  \item{exdat}{
    a \eqn{p}-variate data frame of points on which the regression will be
    estimated (evaluation data).By default,
    evaluation takes place on the data provided by \code{txdat}.
  }

  \item{eydat}{
    a one (1) dimensional numeric or integer vector of the true values
    of the dependent variable. Optional, and used only to calculate the
    true errors. 
  }

  \item{ezdat}{
    an optionally specified \eqn{q}-variate data frame of points on
    which
    the regression will be estimated
    (evaluation data), which corresponds to \eqn{Z}
    in the model equation. Defaults to be the same as \code{txdat}.
  }

  \item{errors}{
    a logical value indicating whether or not asymptotic standard errors
    should be computed and returned in the resulting
    \code{smoothcoefficient} object. Defaults to \code{TRUE}.
  }

  \item{residuals}{
    a logical value indicating that you want residuals computed and
    returned in the resulting \code{smoothcoefficient} object. Defaults to
    \code{FALSE}.
  }


  \item{iterate}{
    a logical value indicating whether or not backfitted estimates
    should be iterated for self-consistency. Defaults to \code{TRUE}.
  }

  \item{maxiter}{
    integer specifying the maximum number of times to iterate the
    backfitted estimates while attempting make the backfitted estimates
    converge to the desired tolerance. Defaults to \code{100}.
  }

  \item{tol}{
    desired tolerance on the relative convergence of backfit
    estimates. Defaults to \code{.Machine$double.eps}.
  }

  \item{leave.one.out}{
    a logical value to specify whether or not to compute the leave one
    out estimates. Will not work if \code{e[xyz]dat} is specified. Defaults to
    \code{FALSE}.
  }

  \item{betas}{
    a logical value indicating whether or not estimates of the
    components of \eqn{\gamma}{gamma} should be returned in the
    \code{smoothcoefficient} object along with the regression estimates.
    Defaults to \code{FALSE}.
  }
  
}

\value{
  \code{npscoef} returns a \code{smoothcoefficient} object.   The generic
  functions \code{\link{fitted}}, \code{\link{residuals}}, \code{\link{coef}},
  \code{\link{se}}, and \code{\link{predict}},
  extract (or generate) estimated values,
  residuals, coefficients, bootstrapped standard
  errors on estimates, and predictions, respectively, from
  the returned object. Furthermore, the functions \code{\link{summary}}
  and \code{\link{plot}} support objects of this type. The returned object
  has the following components:

  \item{eval}{ evaluation points }
  \item{mean}{ estimation of the regression function (conditional mean) at the
    evaluation points } 
  \item{merr}{ if \code{errors = TRUE}, standard errors of the
    regression estimates } 
  \item{beta}{ if \code{betas = TRUE}, estimates of the coefficients
    \eqn{\gamma}{gamma} at the 
    evaluation points }
  \item{resid}{ if \code{residuals = TRUE}, in-sample or out-of-sample
    residuals where appropriate (or possible) }
  \item{R2}{ coefficient of determination }
  \item{MSE}{ mean squared error }
  \item{MAE}{ mean absolute error }
  \item{MAPE}{ mean absolute percentage error }
  \item{CORR}{ absolute value of Pearson's correlation coefficient }
  \item{SIGN}{ fraction of observations where fitted and observed values
  agree in sign }
  
}
\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Cai Z. (2007), \dQuote{Trending time-varying coefficient time series
    models with serially correlated errors,} Journal of Econometrics,
  136, 163-188.

  Hastie, T. and R. Tibshirani (1993), \dQuote{Varying-coefficient
    models,} Journal of the Royal Statistical Society, B 55, 757-796.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
    Theory and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2010), \dQuote{Smooth varying-coefficient
    estimation and inference for qualitative and quantitative data,}
    Econometric Theory, 26, 1-31.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.

 Racine, J.S. and D. Ouyang and Q. Li (2010), \dQuote{Semiparametric
    Hierarchical Models,} manuscript.
 
  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth
    estimators for discrete distributions,}  Biometrika, 68, 301-309.

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

  Support for backfitted bandwidths is experimental and is limited in
  functionality. The code does not support asymptotic standard errors
  or out of sample estimates with backfitting.
}

\seealso{
  \code{\link{bw.nrd}}, \code{\link{bw.SJ}}, \code{\link{hist}},
  \code{\link{npudens}}, \code{\link{npudist}},
  \code{\link{npudensbw}}, \code{\link{npscoefbw}} }

\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): 

n <- 250
x <- runif(n)
z <- runif(n, min=-2, max=2)
y <- x*exp(z)*(1.0+rnorm(n,sd = 0.2))
bw <- npscoefbw(y~x|z)
model <- npscoef(bw)
plot(model)

\dontrun{
# EXAMPLE 1 (INTERFACE=DATA FRAME): 

n <- 250
x <- runif(n)
z <- runif(n, min=-2, max=2)
y <- x*exp(z)*(1.0+rnorm(n,sd = 0.2))
bw <- npscoefbw(xdat=x, ydat=y, zdat=z)
model <- npscoef(bw)
plot(model)
} % enddontrun
}
\keyword{nonparametric}

