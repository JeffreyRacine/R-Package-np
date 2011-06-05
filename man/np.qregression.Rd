% $ID$
\name{npqreg}
\alias{npqreg}
\alias{npqreg.call}
\alias{npqreg.conbandwidth}
\alias{npqreg.default}
\alias{npqreg.formula}
\title{Kernel Quantile Regression with Mixed Data Types}

\description{
  \code{npqreg} computes a kernel quantile regression estimate of a one
  (1) dimensional dependent variable on \eqn{p}-variate explanatory
  data, given a set of evaluation points, training points (consisting of
  explanatory data and dependent data), and a bandwidth specification
  using the method of Li and Racine (2008). A bandwidth specification
  can be a \code{conbandwidth} object, or a bandwidth vector, bandwidth
  type and kernel type.
}

\usage{
npqreg(bws, \dots)

\method{npqreg}{formula}(bws, data = NULL, newdata = NULL, \dots)

\method{npqreg}{call}(bws, \dots)

\method{npqreg}{conbandwidth}(bws,
       txdat = stop("training data 'txdat' missing"),
       tydat = stop("training data 'tydat' missing"),
       exdat,
       tau = 0.5,
       gradients = FALSE,
       ftol = 1.19209e-07,
       tol = 1.49012e-08,
       small = 2.22045e-16,
       itmax = 10000,
       \dots)

\method{npqreg}{default}(bws, txdat, tydat, \dots)

}

\arguments{
  \item{bws}{
    a bandwidth specification. This can be set as a \code{conbandwidth}
    object returned from an invocation of \code{\link{npcdensbw}}, or
    as a vector of bandwidths, with each element \eqn{i} corresponding
    to the bandwidth for column \eqn{i} in \code{txdat}. If specified as
    a vector, then additional arguments will need to be supplied as
    necessary to specify the bandwidth type, kernel types, and so on.
  }

  \item{tau}{
    a numeric value specifying the \eqn{\tau}{tau}th quantile is
    desired. Defaults to \code{0.5}.
  }

  \item{\dots}{
    additional arguments supplied to specify the regression type,
    bandwidth type, kernel types, training data, and so on.
    To do this,
    you may specify any of \code{bwmethod}, \code{bwscaling},
    \code{bwtype}, \code{cxkertype}, \code{cxkerorder},
    \code{cykertype}, \code{cykerorder}, \code{uxkertype},
    \code{uykertype}, \code{oxkertype}, \code{oykertype}, as described
    in \code{\link{npcdensbw}}. 
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npcdensbw}} was called.
  }

  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.  
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

  \item{exdat}{
    a \eqn{p}-variate data frame of points on which the regression will be
    estimated (evaluation data). By default,
    evaluation takes place on the data provided by \code{txdat}.
  }

  \item{gradients}{
    [currently not supported] a logical value indicating that you want
    gradients computed and returned in the resulting \code{npregression}
    object. Defaults to \code{FALSE}.
  }

  \item{itmax}{
    integer number of iterations before failure in the numerical
    optimization routine. Defaults to \code{10000}.
  }

  \item{ftol}{
    tolerance on the value of the objective function 
    evaluated at located minima. Defaults to \code{1.19e-07 (FLT_EPSILON)}.
  }

  \item{tol}{
    tolerance on the position of located minima of the 
    objective function. Defaults to \code{1.49e-08 (sqrt(DBL_EPSILON))}.
  }

  \item{small}{ a small number, at about the precision of the data type
    used. Defaults to \code{2.22e-16 (DBL_EPSILON)}.
  }

}

\value{
  \code{npqreg} returns a \code{npqregression} object.
  The generic
  functions \code{\link{fitted}} (or \code{\link{quantile}}),
  \code{\link{se}}, \code{\link{predict}}, and
  \code{\link{gradients}}, extract (or generate) estimated values,
  asymptotic standard 
  errors on estimates, predictions, and gradients, respectively, from
  the returned object. Furthermore, the functions \code{\link{summary}}
  and \code{\link{plot}} support objects of this type. The returned object
  has the following components:

  \item{eval}{ evaluation points }
  \item{quantile}{ estimation of the quantile regression function
    (conditional quantile) at the  evaluation points }
  \item{quanterr}{ standard errors of the quantile regression estimates } 
  \item{quantgrad}{ gradients at each evaluation point }
  \item{tau}{ the \eqn{\tau}{tau}th quantile computed }

  
}

\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Hall, P. and J.S. Racine and Q. Li (2004), \dQuote{Cross-validation and the
    estimation of conditional probability densities,} Journal of the
  American Statistical Association, 99, 1015-1026.

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
}

\seealso{ \pkg{quantreg} }

\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we compute a
# bivariate nonparametric quantile regression estimate for Giovanni
# Baiocchi's Italian income panel (see Italy for details)

data("Italy")
attach(Italy)

# First, compute the likelihood cross-validation bandwidths (default). We
# override the default tolerances for the search method as the objective
# function is well-behaved (don't of course do this in general).  Note -
# this may take a few minutes depending on the speed of your computer...

bw <- npcdensbw(formula=gdp~ordered(year), tol=.1, ftol=.1)

# Note - numerical search for computing the quantiles will take a
# minute or so...

model.q0.25 <- npqreg(bws=bw, tau=0.25)
model.q0.50 <- npqreg(bws=bw, tau=0.50)
model.q0.75 <- npqreg(bws=bw, tau=0.75)

# Plot the resulting quantiles manually...

plot(ordered(year), gdp, 
     main="CDF Quantile Estimates for the Italian Income Panel", 
     xlab="Year", 
     ylab="GDP Quantiles")

lines(ordered(year), model.q0.25$quantile, col="red", lty=2)
lines(ordered(year), model.q0.50$quantile, col="blue", lty=3)
lines(ordered(year), model.q0.75$quantile, col="red", lty=2)

legend(ordered(1951), 32, c("tau = 0.25", "tau = 0.50", "tau = 0.75"), 
       lty=c(2, 3, 2), col=c("red", "blue", "red"))

detach(Italy)

\dontrun{

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we compute a
# bivariate nonparametric quantile regression estimate for Giovanni
# Baiocchi's Italian income panel (see Italy for details)

data("Italy")
attach(Italy)
data <- data.frame(ordered(year), gdp)

# First, compute the likelihood cross-validation bandwidths (default). We
# override the default tolerances for the search method as the objective
# function is well-behaved (don't of course do this in general).  Note -
# this may take a few minutes depending on the speed of your computer...

bw <- npcdensbw(xdat=ordered(year), ydat=gdp, tol=.1, ftol=.1)

# Note - numerical search for computing the quantiles will take a
# minute or so...

model.q0.25 <- npqreg(bws=bw, tau=0.25)
model.q0.50 <- npqreg(bws=bw, tau=0.50)
model.q0.75 <- npqreg(bws=bw, tau=0.75)

# Plot the resulting quantiles manually...

plot(ordered(year), gdp, 
     main="CDF Quantile Estimates for the Italian Income Panel", 
     xlab="Year", 
     ylab="GDP Quantiles")

lines(ordered(year), model.q0.25$quantile, col="red", lty=2)
lines(ordered(year), model.q0.50$quantile, col="blue", lty=3)
lines(ordered(year), model.q0.75$quantile, col="red", lty=2)

legend(ordered(1951), 32, c("tau = 0.25", "tau = 0.50", "tau = 0.75"), 
       lty=c(2, 3, 2), col=c("red", "blue", "red"))

detach(Italy)
} % enddontrun
}
\keyword{nonparametric}
