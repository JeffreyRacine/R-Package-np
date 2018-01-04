\name{npuniden.reflect}
\alias{npuniden.reflect}

\title{ Kernel Bounded Univariate Density Estimation Via Data-Reflection }

\description{

  \code{npuniden.reflect} computes kernel univariate unconditional
  density estimates given a vector of continuously distributed training
  data and, optionally, a bandwidth (otherwise likelihood
  cross-validation is used for its selection). Lower and upper bounds
  [\code{a},\code{b}] can be supplied (default is [0,1]) and if \code{a}
  is set to \code{-Inf} there is only one bound on the right, while if
  \code{b} is set to \code{Inf} there is only one bound on the left.

}

\usage{
npuniden.reflect(X = NULL,
                 h = NULL,
                 a = 0,
                 b = 1,
                 ...)
}
\arguments{
  \item{X}{
a numeric vector of data
}
  \item{h}{
an optional bandwidth (>0)
}
  \item{a}{
an optional lower bound (defaults to 0)
}
  \item{b}{
an optional upper bound (defaults to 1)
}
 \item{...}{
optional arguments passed to \code{npudensbw} and \code{npudens}
}
}
\details{
  Typical usages are (see below for a complete list of options and also
  the examples at the end of this help file)
  \preformatted{
    model <- npuniden.reflect(X,a=-2,b=3)
  }

  \code{npuniden.reflect} implements the data-reflection method for
  estimating a univariate density function defined over a continuous
  random variable in the presence of bounds.
  
}
\value{

  \code{npuniden.reflect} returns the following components:
  
\item{f}{ estimated density at the points X }
\item{h}{ bandwidth used }
\item{nmulti}{ number of multistarts used }
\item{int}{numeric integral of the estimated density }
}
\references{

  Boneva, L. I., Kendall, D., and Stefanov, I. (1971). \dQuote{Spline
  transformations: Three new diagnostic aids for the statistical data-
  analyst}. Journal of the Royal Statistical Society. Series B (Method-
  ological), 33(1):1-71.

  Cline, D. B. H. and Hart, J. D. (1991). \dQuote{Kernel estimation of
  densities with discontinuities or discontinuous
  derivatives}. Statistics, 22(1):69-84.

  Hall, P. and Wehrly, T. E. (1991). \dQuote{A geometrical method for
  removing edge effects from kernel- type nonparametric regression
  estimators}. Journal of the American Statistical Association,
  86(415):665-672.
  
}
\author{
Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\seealso{
The \pkg{bde} package.
}
\examples{
\dontrun{
set.seed(42)
n <- 100
X <- sort(rbeta(n,5,1))
dgp <- dbeta(X,5,1)
model <- npuniden.reflect(X)
model.unadjusted <- npuniden.boundary(X,a=-Inf,b=Inf)
plot(X,model$f,type="l")
lines(X,model.unadjusted$f,lty=2,col=2)
lines(X,dgp,lty=3,col=3)
legend("topleft",c("Data-Reflection","Unadjusted","DGP"),col=1:3,lty=1:3,bty="n")
} % enddontrun
}
\keyword{ nonparametric }
\keyword{ smooth }
