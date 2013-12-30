\name{npquantile}
\alias{npquantile}

\title{ Kernel Univariate Quantile Estimation }
\description{

  \code{npquantile} computes smooth quantiles from a univariate kernel
  unconditional cumulative distribution estimate given data and a
  bandwidth specification (a \code{dbandwidth} object) using the
  bandwidth selection method of Li and Racine (2013).

}
\usage{
npquantile(bws = NULL,
           dat = NULL,
           tau = c(0.01,0.05,0.25,0.50,0.75,0.95,0.99),
           num.eval = 5000)
}

\arguments{
  \item{bws}{
    a \code{dbandwidth} specification. This must be set as a \code{dbandwidth}
    object returned from an invocation of \code{\link{npudistbw}}.
  }
  \item{dat}{
    a univariate vector of type \code{\link{numeric}} containing sample
    realizations (training data) used to estimate the cumulative
    distribution (must be the same training data used to compute the
    bandwidth object \code{bws} passed in).
  }
  \item{tau}{
    an optional vector containing the quantile(s) to be estimated (must
    contain numbers strictly between 0 and 1). Defaults to
    \code{c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)}.
  }
  \item{num.eval}{
    an optional integer specifying the length of the grid on which the
    quasi-inverse is computed. Defaults to \code{5000}.
  }
  
}

\details{

  The quantile function \eqn{q_\tau} is defined as the quasi-inverse of
  the distribution function \eqn{F(x)}, i.e. \eqn{q_\tau = \inf\{x: F(x)
  \ge \tau\}}.
 
  The function \code{npquantile} implements a method for estimating
  smooth quantiles based on the quasi-inverse of a \code{\link{npudist}}
  object where \eqn{F(x)} is replaced with its kernel estimator and
  bandwidth selection is that appropriate for such objects. For
  construction of the quasi-inverse we create a grid of points based on
  the function \code{\link{extendrange}} along with the sample quantiles
  themselves computed from invocation of \code{\link{quantile}}. The
  coarseness of the grid defined by \code{\link{extendrange}} (which has
  been passed the option \code{f=1}) is controlled by \code{num.eval}.
  
  See Definition 2.3.6, page 21, Nelsen 2006 for further details on
  computation of the quasi-inverse of \eqn{F(x)}.


}
\value{
  
  \code{\link{npquantile}} returns a vector of quantiles corresponding
  to \code{tau}.
  
}
\references{
  Nelsen, R. B. (2006), \emph{An Introduction to Copulas,} Second
  Edition, Springer-Verlag.

  Li, H. and J.S. Racine (2013), \dQuote{Cross-validated estimation of
  cumulative distribution functions with categorical and continuous
  data,} manuscript.
}
\author{
    Tristen Hayfield \email{hayfield@mpia.de}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\seealso{
  \code{\link{quantile}}, \code{\link[logspline]{qlogspline}},
  \code{\link[ks]{qkde}}, \code{\link{ecdf}}
  }

\examples{
\dontrun{
## Simulate data
set.seed(42)
df <- 50
x <- rchisq(100,df=df)

## Compute npudist bandwidth object
bw <- npudistbw(~x)

## Pass the bandwidth object and data to npquantile()
npquantile(bw,x)

## Compute the nonsmooth empirical quantiles
quantile(x,c(0.01,0.05,0.25,0.50,0.75,0.95,0.99))

## Compute the true quantiles
qchisq(c(0.01,0.05,0.25,0.50,0.75,0.95,0.99),df=df)
} % enddontrun
}
\keyword{ nonparametric }% at least one, from doc/KEYWORDS

