\name{npquantile}
\alias{npquantile}

\title{ Kernel Univariate Quantile Estimation }
\description{

  \code{npquantile} computes smooth quantiles from a univariate kernel
  unconditional cumulative distribution estimate given data and,
  optionally, a bandwidth specification i.e. a \code{dbandwidth} object
  using the bandwidth selection method of Li and Racine (2013).

}
\usage{
npquantile(x = NULL,
           tau = c(0.01,0.05,0.25,0.50,0.75,0.95,0.99),
           num.eval = 10000,
           bws = NULL,
           f = 1,
           \dots)
}

\arguments{
  \item{x}{
    a univariate vector of type \code{\link{numeric}} containing sample
    realizations (training data) used to estimate the cumulative
    distribution (must be the same training data used to compute the
    bandwidth object \code{bws} passed in).
  }
  \item{tau}{
    an optional vector containing the probabilities for quantile(s) to
    be estimated (must contain numbers in \eqn{[0,1]}). Defaults to
    \code{c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)}.
  }
  \item{num.eval}{
    an optional integer specifying the length of the grid on which the
    quasi-inverse is computed. Defaults to \code{10000}.
  }
  \item{bws}{
    an optional \code{dbandwidth} specification (if already computed
    avoid unnecessary computation inside \code{npquantile}). This must
    be set as a \code{dbandwidth} object returned from an invocation of
    \code{\link{npudistbw}}. If not provided \code{\link{npudistbw}} is
    invoked with optional arguments passed via \code{\dots}.
  }
  \item{f}{
    an optional argument fed to \code{\link{extendrange}}. Defaults to
    \code{1}. See \code{?\link{extendrange}} for details.
  }
  \item{\dots}{
    additional arguments supplied to specify the bandwidth type, kernel
    types, bandwidth selection methods, and so on. See
    \code{?\link{npudistbw}} for details.
  }
  
}

\details{

  Typical usage is
  \preformatted{
    x <- rchisq(100,df=10)
    npquantile(x)
  }

  The quantile function \eqn{q_\tau} is defined to be the
  left-continuous inverse of the distribution function \eqn{F(x)},
  i.e. \eqn{q_\tau = \inf\{x: F(x) \ge \tau\}}.

  A traditional estimator of \eqn{q_\tau} is the \eqn{\tau}th sample
  quantile. However, these estimates suffer from lack of efficiency
  arising from variability of individual order statistics; see Sheather
  and Marron (1990) and Hyndman and Fan (1996) for methods that
  interpolate/smooth the order statistics, each of which discussed in
  the latter can be invoked through \code{\link{quantile}} via
  \code{type=j}, \code{j=1,\dots,9}.
 
  The function \code{npquantile} implements a method for estimating
  smooth quantiles based on the quasi-inverse of a \code{\link{npudist}}
  object where \eqn{F(x)} is replaced with its kernel estimator and
  bandwidth selection is that appropriate for such objects; see
  Definition 2.3.6, page 21, Nelsen 2006 for a definition of the
  quasi-inverse of \eqn{F(x)}.

  For construction of the quasi-inverse we create a grid of evaluation
  points based on the function \code{\link{extendrange}} along with the
  sample quantiles themselves computed from invocation of
  \code{\link{quantile}}. The coarseness of the grid defined by
  \code{\link{extendrange}} (which has been passed the option
  \code{f=1}) is controlled by \code{num.eval}.
  
  Note that for any value of \eqn{\tau} less/greater than the
  smallest/largest value of \eqn{F(x)} computed for the evaluation data
  (i.e. that outlined in the paragraph above), the quantile returned for
  such values is that associated with the smallest/largest value of
  \eqn{F(x)}, respectively.

}
\value{
  
  \code{\link{npquantile}} returns a vector of quantiles corresponding
  to \code{tau}.
  
}
\references{

  Hyndman, R. J. and Fan, Y. (1996), \dQuote{Sample quantiles in
  statistical packages,} American Statistician, \bold{50}, 361-365.

  Li, H. and J.S. Racine (2013), \dQuote{Cross-validated estimation of
  cumulative distribution functions with categorical and continuous
  data,} manuscript.

  Nelsen, R. B. (2006), \emph{An Introduction to Copulas,} Second
  Edition, Springer-Verlag.

  Sheather, S. and J. S. Marron (1990), \dQuote{Kernel Quantile
  Estimators,} Journal of the American Statistical Association,
  Vol. 85, No. 410, 410-416.
  
}
\section{Usage Issues}{
  Cross-validated bandwidth selection is used by default
  (\code{\link{npudistbw}}). For large datasets this can be
  computationally demanding. In such cases one might instead consider a
  rule-of-thumb bandwidth (\code{bwmethod="normal-reference"}) or,
  alternatively, use kd-trees (\code{options(np.tree=TRUE)} along with a
  bounded kernel (\code{ckertype="epanechnikov"})), both of which will
  reduce the computational burden appreciably.
}
\author{
    Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\seealso{

  \code{\link{quantile}} for various types of sample quantiles;
  \code{\link{ecdf}} for empirical distributions of which
  \code{\link{quantile}} is an inverse; \code{\link{boxplot.stats}} and
  \code{\link{fivenum}} for computing other versions of quartiles;
  \code{\link[logspline]{qlogspline}} for logspline density quantiles;
  \code{\link[ks]{qkde}} for alternative kernel quantiles, etc.

}

\examples{
\dontrun{
## Simulate data from a chi-square distribution
df <- 50
x <- rchisq(100,df=df)

## Vector of quantiles desired
tau <- c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)

## Compute kernel smoothed sample quantiles
npquantile(x,tau)

## Compute sample quantiles using the default method in R (Type 7)
quantile(x,tau)

## True quantiles based on known distribution
qchisq(tau,df=df)
} % enddontrun
}
\keyword{ nonparametric }% at least one, from doc/KEYWORDS

