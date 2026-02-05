% $ID$
\name{npqreg}
\alias{npqreg}
\alias{npqreg.call}
\alias{npqreg.condbandwidth}
\alias{npqreg.default}
\alias{npqreg.formula}
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

\method{npqreg}{call}(bws, \dots)

\method{npqreg}{condbandwidth}(bws,
       txdat = stop("training data 'txdat' missing"),
       tydat = stop("training data 'tydat' missing"),
       exdat,
       tau = 0.5,
       gradients = FALSE,
       ftol = 1.490116e-07,
       tol = 1.490116e-04,
       small = 1.490116e-05,
       itmax = 10000,
       lbc.dir = 0.5,
       dfc.dir = 3,
       cfac.dir = 2.5*(3.0-sqrt(5)),
       initc.dir = 1.0,
       lbd.dir = 0.1,
       hbd.dir = 1,
       dfac.dir = 0.25*(3.0-sqrt(5)),
       initd.dir = 1.0,
       \dots)

\method{npqreg}{default}(bws, txdat, tydat, \dots)

}

\arguments{
  \item{bws}{
    a bandwidth specification. This can be set as a \code{condbandwidth}
    object returned from an invocation of \code{\link{npcdistbw}}, or
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
    in \code{\link{npcdistbw}}. 
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npcdistbw}} was called.
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
    fractional tolerance on the value of the cross-validation function
    evaluated at located minima (of order the machine precision or
    perhaps slightly larger so as not to be diddled by
    roundoff). Defaults to \code{1.490116e-07}
    (1.0e+01*sqrt(.Machine$double.eps)).
  }

  \item{tol}{
    tolerance on the position of located minima of the cross-validation
    function (tol should generally be no smaller than the square root of
    your machine's floating point precision). Defaults to \code{
      1.490116e-04 (1.0e+04*sqrt(.Machine$double.eps))}.
  }

  \item{small}{
    a small number used to bracket a minimum (it is hopeless to ask for
    a bracketing interval of width less than sqrt(epsilon) times its
    central value, a fractional width of only about 10-04 (single
    precision) or 3x10-8 (double precision)). Defaults to \code{small
      = 1.490116e-05 (1.0e+03*sqrt(.Machine$double.eps))}.
  }

  \item{lbc.dir,dfc.dir,cfac.dir,initc.dir}{ lower bound, chi-square
    degrees of freedom, stretch factor, and initial non-random values
    for direction set search for Powell's algorithm for \code{numeric}
    variables. See Details}

  \item{lbd.dir,hbd.dir,dfac.dir,initd.dir}{ lower bound, upper bound,
    stretch factor, and initial non-random values for direction set
    search for Powell's algorithm for categorical variables. See
    Details}

}

\details{

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
  \sQuote{tuning} for various methods rather than for the user themselves.

}

\value{

  \code{npqreg} returns a \code{npqregression} object.  The generic
  functions \code{\link{fitted}} (or \code{\link{quantile}}),
  \code{\link{se}}, \code{\link{predict}} (when using
  \code{\link{predict}} you must add the argument \code{tau=} to
  generate predictions other than the median), and
  \code{\link{gradients}}, extract (or generate) estimated values,
  asymptotic standard errors on estimates, predictions, and gradients,
  respectively, from the returned object. Furthermore, the functions
  \code{\link{summary}} and \code{\link{plot}} support objects of this
  type. The returned object has the following components:

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

\seealso{ \pkg{quantreg} }

\examples{
\dontrun{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we compute a
# bivariate nonparametric quantile regression estimate for Giovanni
# Baiocchi's Italian income panel (see Italy for details)

data("Italy")
attach(Italy)

# First, compute the cross-validated bandwidths.  Note - this may take a
# few minutes depending on the speed of your computer...

bw <- npcdistbw(formula=gdp~ordered(year))

# Note - numerical search for computing the quantiles may take a minute
# or so...

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

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we compute a
# bivariate nonparametric quantile regression estimate for Giovanni
# Baiocchi's Italian income panel (see Italy for details)

data("Italy")
attach(Italy)
data <- data.frame(ordered(year), gdp)

# First, compute the likelihood cross-validation bandwidths (default).
# Note - this may take a few minutes depending on the speed of your
# computer...

bw <- npcdistbw(xdat=ordered(year), ydat=gdp)

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
