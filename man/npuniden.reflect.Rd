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

  Note that data-reflection imposes a zero derivative at the boundary,
  i.e., \eqn{f'(a)=f'(b)=0}.
  
}
\value{
  \code{npuniden.reflect} returns the following components:
  
  \item{f}{ estimated density at the points X }
  \item{F}{ estimated distribution at the points X (numeric integral of f) }
  \item{sd.f}{ asymptotic standard error of the estimated density at the points X }
  \item{sd.F}{ asymptotic standard error of the estimated distribution at the points X }
  \item{h}{ bandwidth used }
  \item{nmulti}{ number of multi-starts used }
}
\references{

  Boneva, L. I., Kendall, D., and Stefanov, I. (1971). \dQuote{Spline
  transformations: Three new diagnostic aids for the statistical data-
  analyst,} Journal of the Royal Statistical Society. Series B
  (Methodological), 33(1):1-71.

  Cline, D. B. H. and Hart, J. D. (1991). \dQuote{Kernel estimation of
  densities with discontinuities or discontinuous
  derivatives,} Statistics, 22(1):69-84.

  Hall, P. and Wehrly, T. E. (1991). \dQuote{A geometrical method for
  removing edge effects from kernel- type nonparametric regression
  estimators,} Journal of the American Statistical Association,
  86(415):665-672.
  
}
\author{
Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\seealso{
The \pkg{Ake}, \pkg{bde}, and \pkg{Conake} packages and the function \code{\link{npuniden.boundary}}.
}
\examples{
\dontrun{
## Example 1: f(0)=0, f(1)=1, plot boundary corrected density,
## unadjusted density, and DGP
set.seed(42)
n <- 100
X <- sort(rbeta(n,5,1))
dgp <- dbeta(X,5,1)
model <- npuniden.reflect(X)
model.unadjusted <- npuniden.boundary(X,a=-Inf,b=Inf)
ylim <- c(0,max(c(dgp,model$f,model.unadjusted$f)))
plot(X,model$f,ylab="Density",ylim=ylim,type="l")
lines(X,model.unadjusted$f,lty=2,col=2)
lines(X,dgp,lty=3,col=3)
rug(X)
legend("topleft",c("Data-Reflection","Unadjusted","DGP"),col=1:3,lty=1:3,bty="n")

## Example 2: f(0)=0, f(1)=0, plot density, distribution, DGP, and
## asymptotic point-wise confidence intervals
set.seed(42)
X <- sort(rbeta(100,5,3))
model <- npuniden.reflect(X)
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
model <- npuniden.reflect(age,a=21,b=65)
par(mfrow=c(1,1))
hist(age,prob=T,main="",ylim=c(0,max(model$f)))
lines(age,model$f)
lines(density(age,bw=model$h),col=2)
legend("topright",c("Data-Reflection","Unadjusted"),lty=c(1,1),col=1:2,bty="n")
detach(cps71)
} % enddontrun
}
\keyword{ nonparametric }
\keyword{ smooth }
