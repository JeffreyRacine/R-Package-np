\name{b.star}
\alias{b.star}
\title{Compute Optimal Block Length for Stationary and Circular Bootstrap}
\usage{
b.star(data,
       Kn = NULL,
       mmax= NULL,
       Bmax = NULL,
       c = NULL,
       round = FALSE)
}
\arguments{
  \item{data}{data, an n x k matrix, each column being a data series.}
  
  \item{Kn}{ See footnote c, page 59, Politis and White (2004). Defaults
    to \code{ceiling(log10(n))}.}

  \item{mmax}{See Politis and White (2004). Defaults to
  \code{ceiling(sqrt(n))+Kn}.}

  \item{Bmax}{See Politis and White (2004). Defaults to
    \code{ceiling(min(3*sqrt(n),n/3))}.}

  \item{c}{See Politis and White (2004). Defaults to
  \code{qnorm(0.975)}.}

  \item{round}{whether to round the result or not. Defaults to FALSE.}
 
}
\value{
  
  A kx2 matrix of optimal bootstrap block lengths computed from
 \code{data} for the stationary bootstrap and circular bootstrap (column
 1 is for the stationary bootstrap, column 2 the circular).

}
\description{
  \code{b.star} is a function which computes the optimal block length
  for the continuous variable \code{data} using the method described in
  Patton, Politis and White (2009).
}
\details{
  \code{b.star} is a function which computes optimal block lengths for
  the stationary and circular bootstraps. This allows the use of
  \code{tsboot} from the \code{\link{boot}} package to be fully
  automatic by using the output from \code{b.star} as an input to the
  argument \code{l = } in \code{tsboot}. See below for an example.
}

\references{

  Patton, A. and D.N. Politis and H. White (2009), \dQuote{CORRECTION TO
    "Automatic block-length selection for the dependent bootstrap" by
    D. Politis and H. White}, Econometric Reviews 28(4), 372-375.

  Politis, D.N. and J.P. Romano (1994), \dQuote{Limit theorems for
    weakly dependent Hilbert space valued random variables with
    applications to the stationary bootstrap}, Statistica Sinica 4,
    461-476.

  Politis, D.N. and H. White (2004), \dQuote{Automatic block-length
    selection for the dependent bootstrap}, Econometric Reviews 23(1),
    53-70.

}

\author{
  Tristen Hayfield \email{hayfield@mpia.de}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\examples{
\dontrun{
set.seed(42)

# Function to generate an AR(1) series

ar.series <- function(phi,epsilon) {
  n <- length(epsilon)
  series <- numeric(n)
  series[1] <- epsilon[1]/(1-phi)
  for(i in 2:n) {
    series[i] <- phi*series[i-1] + epsilon[i]
  }
  return(series)
}

yt <- ar.series(0.1,rnorm(10000))
b.star(yt,round=TRUE)

yt <- ar.series(0.9,rnorm(10000))
b.star(yt,round=TRUE)
} % enddontrun
}

\keyword{univar}
\keyword{nonparametric}
