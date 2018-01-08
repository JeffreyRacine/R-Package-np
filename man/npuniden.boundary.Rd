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
  \code{a} is set to \code{-Inf} and \code{b} to \code{Inf} and the
  Gaussian type 1 kernel function is used, this will deliver the
  standard unadjusted kernel density estimate.

}

\usage{
npuniden.boundary(X = NULL,
                  h = NULL,
                  a = 0,
                  b = 1,
                  kertype = c("gaussian1","gaussian2","beta1","beta2","gamma"),
                  cv = c("grid-hybrid","numeric"),
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
an optional kernel specification (defaults to "gaussian1")
}
  \item{cv}{
an optional argument for search (default is likely more reliable in the
presence of local maxima)
}
  \item{grid}{
an optional grid used for the initial grid search when \code{cv="grid-hybrid"}
}
  \item{nmulti}{
number of multi-starts used when \code{cv="numeric"} (defaults to 5)
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
  (1999) type 1 and 2 kernel functions with biases of \eqn{O(h)}, the
  \code{"gamma"} is from Chen (2000), and \code{"gaussian1"} and
  \code{"gaussian2"} are renormalized Gaussian kernel functions with
  biases of \eqn{O(h)} and \eqn{O(h^2)}.

}
\value{
  \code{npuniden.boundary} returns the following components:
  
  \item{f}{ estimated density at the points X }
  \item{F}{ estimated distribution at the points X (numeric integral of f) }
  \item{sd.f}{ asymptotic standard error of the estimated density at the points X }
  \item{sd.F}{ asymptotic standard error of the estimated distribution at the points X }
  \item{h}{ bandwidth used }
  \item{nmulti}{ number of multi-starts used }
}
\references{

  Bouezmarni, T. and Rolin, J.-M. (2003). \dQuote{Consistency of the
  beta kernel density function estimator,}  The Canadian Journal of
  Statistics / La Revue Canadienne de Statistique, 31(1):89-98.

  Chen, S. X. (1999). \dQuote{Beta kernel estimators for density
  functions,} Computational Statistics & Data Analysis, 31(2):131-145.

  Chen, S. X. (2000). \dQuote{Probability density function estimation
  using gamma kernels,} Annals of the Institute of Statistical
  Mathematics, 52(3):471-480.

  Diggle, P. (1985). \dQuote{A kernel method for smoothing point process
  data,} Journal of the Royal Statistical Society. Series C (Applied
  Statistics), 34(2):138-147.

  Igarashi, G. (2016). \dQuote{Bias Reductions for Beta Kernel
  Estimation,} Journal of Nonparametric Statistics, 28(1):1-30.
  
  Zhang, S. and R. J. Karunamuni (2010). \dQuote{Boundary Performance of
  the Beta Kernel Estimators,} Journal of Nonparametric Statistics,
  22(1):81-104.
  
}
\author{
Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\seealso{
The \pkg{bde} package and the function \code{\link{npuniden.reflect}}.
}
\examples{
\dontrun{
## Example 1: f(0)=0, f(1)=1, plot boundary corrected density,
## unadjusted density, and DGP
set.seed(42)
n <- 100
X <- sort(rbeta(n,5,1))
dgp <- dbeta(X,5,1)
model.g1 <- npuniden.boundary(X)
model.g2 <- npuniden.boundary(X,kertype="gaussian2")
model.b1 <- npuniden.boundary(X,kertype="beta1")
model.b2 <- npuniden.boundary(X,kertype="beta2")
model.unadjusted <- npuniden.boundary(X,a=-Inf,b=Inf)
ylim <- c(0,max(c(dgp,model.g1$f,model.g2$f)))
plot(X,dgp,ylab="Density",ylim=ylim,type="l")
lines(X,model.g1$f,lty=2,col=2)
lines(X,model.g2$f,lty=3,col=3)
lines(X,model.b1$f,lty=4,col=4)
lines(X,model.b2$f,lty=5,col=5)
lines(X,model.unadjusted$f,lty=6,col=6)
rug(X)
legend("topleft",c("DGP",
                   "Boundary Kernel (gaussian1)",
                   "Boundary Kernel (gaussian2)",
                   "Boundary Kernel (beta1)",
                   "Boundary Kernel (beta2)",
                   "Unadjusted"),col=1:6,lty=1:6,bty="n")

## Example 2: f(0)=0, f(1)=0, plot density, distribution, DGP, and
## asymptotic point-wise confidence intervals
set.seed(42)
X <- sort(rbeta(100,5,3))
model <- npuniden.boundary(X)
par(mfrow=c(1,2))
ylim=range(c(model$f,model$f+1.96*model$sd.f,model$f-1.96*model$sd.f,dbeta(X,5,3)))
plot(X,model$f,ylim=ylim,ylab="Density",type="l",)
lines(X,model$f+1.96*model$sd.f,lty=2)
lines(X,model$f-1.96*model$sd.f,lty=2)
lines(X,dbeta(X,5,3),col=2)
rug(X)
legend("topleft",c("Density","DGP"),lty=c(1,1),col=1:2,bty="n")

plot(X,model$F,ylab="Distribution",type="l")
lines(X,model$F+1.96*model$sd.F,lty=2)
lines(X,model$F-1.96*model$sd.F,lty=2)
lines(X,pbeta(X,5,3),col=2)
rug(X)
legend("topleft",c("Distribution","DGP"),lty=c(1,1),col=1:2,bty="n")

## Example 3: Age for working age males in the cps71 data set bounded
## below by 21 and above by 65
data(cps71)
attach(cps71)
model <- npuniden.boundary(age,a=21,b=65)
par(mfrow=c(1,1))
hist(age,prob=T,main="")
lines(age,model$f)
lines(density(age,bw=model$h),col=2)
legend("topright",c("Boundary Kernel","Unadjusted"),lty=c(1,1),col=1:2,bty="n")
detach(cps71)
} % enddontrun
}
\keyword{ nonparametric }
\keyword{ smooth }
