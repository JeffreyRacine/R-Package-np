\name{npuniden.boundary}
\alias{npuniden.boundary}

\title{ Kernel Bounded Univariate Density Estimation Via Boundary Kernel
  Functions }

\description{

  \code{npuniden.boundary} computes kernel univariate unconditional
  density estimates given a vector of continuously distributed training
  data and, optionally, a bandwidth (otherwise likelihood
  cross-validation is used for its selection). Lower and upper bounds
  [\code{a},\code{b}] can be supplied (default is [0,1]) and if \code{a}
  is set to \code{-Inf} there is only one bound on the right, while if
  \code{b} is set to \code{Inf} there is only one bound on the left. If
  \code{a} is set to \code{-Inf} and \code{b} to \code{Inf} and either
  the Gaussian or Epanechnikov kernel functions are used, this will
  deliver the standard unadjusted kernel density estimate.

}

\usage{
npuniden.boundary(X = NULL,
                  h = NULL,
                  a = 0,
                  b = 1,
                  kertype = c("gaussian", "epanechnikov", "beta1", "beta2", "gamma"),
                  cv = c("grid-hybrid", "numeric"),
                  grid = NULL,
                  nmulti = 5)
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
  \item{kertype}{
an optional kernel specification (defaults to "gaussian")
}
  \item{cv}{
an optional argument for search (default is likely more reliable in the
presence of local maxima)
}
  \item{grid}{
an optional grid used for the initial grid search when \code{cv="grid-hybrid"}
}
  \item{nmulti}{
number of multistarts used when \code{cv="numeric"} (defaults to 5)
}
}
\details{
  Typical usages are (see below for a complete list of options and also
  the examples at the end of this help file)
  \preformatted{
    model <- npuniden.boundary(X,a=-2,b=3)
  }

  \code{npuniden.boundary} implements a variety of methods for
  estimating a univariate density functions defined over a continuous
  random variable in the presence of bounds via the use of so-called
  boundary or edge kernel functions.

  The kernel functions \code{"beta1"} and \code{"beta2"} are Chen's
  (1999) type 1 and 2 kernel functions, the \code{"gamma"} is from Chen
  (2000), and \code{"gaussian"}, and \code{"epanechnikov"} use Diggle's
  (1985) end-correction.

}
\value{

  \code{npuniden.boundary} returns the following components:
  
\item{f}{ estimated density at the points X }
\item{h}{ bandwidth used }
\item{nmulti}{ number of multistarts used }
\item{int}{numeric integral of the estimated density }
}
\references{

  Bouezmarni, T. and Rolin, J.-M. (2003). \dQuote{Consistency of the
  beta kernel density function estimator}.  The Canadian Journal of
  Statistics / La Revue Canadienne de Statistique, 31(1):89-98.

  Chen, S. X. (1999). \dQuote{Beta kernel estimators for density
  functions}. Computational Statistics & Data Analysis, 31(2):131-145.

  Chen, S. X. (2000). \dQuote{Probability density function estimation using
  gamma kernels}. Annals of the Institute of Statistical Mathematics,
  52(3):471-480.

  Diggle, P. (1985). \dQuote{A kernel method for smoothing point process
  data}. Journal of the Royal Statistical Society. Series C (Applied
  Statistics), 34(2):138-147.
  
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
model <- npuniden.boundary(X)
model.unadjusted <- npuniden.boundary(X,a=-Inf,b=Inf)
plot(X,model$f,type="l")
lines(X,model.unadjusted$f,lty=2,col=2)
lines(X,dgp,lty=3,col=3)
legend("topleft",c("Boundary Kernel","Unadjusted","DGP"),col=1:3,lty=1:3,bty="n")
} % enddontrun
}
\keyword{ nonparametric }
\keyword{ smooth }
