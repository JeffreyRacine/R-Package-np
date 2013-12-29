\name{npquantile}
\alias{npquantile}

\title{ Kernel Quantile Estimation }
\description{

  \code{npquantile} computes smooth quantiles from a kernel
  unconditional cumulative distribution estimate given data and
  a bandwidth specification (a \code{dbandwidth} object) using the
  method of Li and Racine (2013).

}
\usage{
npquantile(bws, dat, \dots)
}

\arguments{
  \item{bws}{
    a \code{dbandwidth} specification. This must be set as a \code{dbandwidth}
    object returned from an invocation of \code{\link{npudistbw}}.
  }
  \item{dat}{
    a univariate vector of type \code{numeric} containing sample
    realizations (training data) used to estimate the cumulative
    distribution (must be the same training data used to compute the
    bandwidth object \code{bws} passed in).
  }
  \item{tau}{
    an optional vector containing the quantile(s) to be estimated (must
    contain numbers strictly between 0 and 1). Defaults to
    \code{c(0.01,0.025,0.50,0.75,0.95,0.99)}.
  }
  \item{num.eval}{
    an optional integer specifying the length of the grid on which the
    quasi-inverse is computed. Defaults to \code{5000}.
  }
  
}

\details{

  \code{npquantile} implements a method for estimating smooth quantiles
  based on the quasi-inverse of a \code{\link{npudist}} object. See
  Definition 2.3.6, page 21, Nelsen 2006 for details.

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

\seealso{ \code{\link{ecdf}} }

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
quantile(x,c(0.01,0.025,0.50,0.75,0.95,0.99))

## Compute the true quantiles
qchisq(c(0.01,0.025,0.50,0.75,0.95,0.99),df=df)
} % enddontrun
}
\keyword{ nonparametric }% at least one, from doc/KEYWORDS

