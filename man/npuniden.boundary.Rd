\name{npuniden.boundary}
\alias{npuniden.boundary}

\title{ Kernel Bounded Univariate Density Estimation Via Boundary Kernel
  Functions }

\description{

  \code{npuniden.boundary} computes kernel univariate unconditional
  density estimates given a vector of continuously distributed training
  data and, optionally, a bandwidth (otherwise least squares
  cross-validation is used for its selection). Lower and upper bounds
  [\code{a},\code{b}] can be supplied (default is the empirical support 
  \eqn{[\min(X),\max(X)]}) and if \code{a}
  is set to \code{-Inf} there is only one bound on the right, while if
  \code{b} is set to \code{Inf} there is only one bound on the left. If
  \code{a} is set to \code{-Inf} and \code{b} to \code{Inf} and the
  Gaussian type 1 kernel function is used, this will deliver the
  standard unadjusted kernel density estimate.

}

\usage{
npuniden.boundary(X = NULL,
                  Y = NULL,
                  h = NULL,
                  a = min(X),
                  b = max(X),
                  bwmethod = c("cv.ls","cv.ml"),
                  cv = c("grid-hybrid","numeric"),
                  grid = NULL,
                  kertype = c("gaussian1","gaussian2",
                              "beta1","beta2",
                              "fb","fbl","fbu",
                              "rigaussian","gamma"),
                  nmulti = 5,
                  proper = FALSE)
}
\arguments{
  \item{X}{
a required numeric vector of training data lying in \eqn{[a,b]}
}
  \item{Y}{
an optional numeric vector of evaluation data lying in \eqn{[a,b]}
}
  \item{h}{
an optional bandwidth (>0)
}
  \item{a}{
an optional lower bound (defaults to lower bound of empirical support \eqn{\min(X)})
}
  \item{b}{
an optional upper bound (defaults to upper bound of empirical support \eqn{\max(X)})
}
  \item{bwmethod}{
whether to conduct bandwidth search via least squares cross-validation
(\code{"cv.ls"}) or likelihood cross-validation (\code{"cv.ml"})
}
  \item{cv}{
an optional argument for search (default is likely more reliable in the
presence of local maxima)
}
  \item{grid}{
an optional grid used for the initial grid search when \code{cv="grid-hybrid"}
}
  \item{kertype}{
an optional kernel specification (defaults to "gaussian1")
}
  \item{nmulti}{
number of multi-starts used when \code{cv="numeric"} (defaults to 5)
}
  \item{proper}{
an optional logical value indicating whether to enforce proper density
and distribution function estimates over the range \eqn{[a,b]}
}
}
\details{
  Typical usages are (see below for a complete list of options and also
  the examples at the end of this help file)
  \preformatted{
    model <- npuniden.boundary(X,a=-2,b=3)
  }

  \code{npuniden.boundary} implements a variety of methods for
  estimating a univariate density function defined over a continuous
  random variable in the presence of bounds via the use of so-called
  boundary or edge kernel functions.

  The kernel functions \code{"beta1"} and \code{"beta2"} are Chen's
  (1999) type 1 and 2 kernel functions with biases of \eqn{O(h)}, the
  \code{"gamma"} kernel function is from Chen (2000) with a bias of
  \eqn{O(h)}, \code{"rigaussian"} is the reciprocal inverse Gaussian
  kernel function (Scaillet (2004), Igarashi & Kakizawa (2014)) with
  bias of \eqn{O(h)}, and \code{"gaussian1"} and \code{"gaussian2"} are
  truncated Gaussian kernel functions with biases of \eqn{O(h)} and
  \eqn{O(h^2)}, respectively. The kernel functions \code{"fb"},
  \code{"fbl"} and \code{"fbu"} are floating boundary polynomial
  biweight kernels with biases of \eqn{O(h^2)} (Scott (1992), Page
  146). Without exception, these kernel functions are asymmetric in
  general with shape that changes depending on where the density is
  being estimated (i.e., how close the estimation point \eqn{x} in
  \eqn{\hat f(x)} is to a boundary). This function is written purely in
  R, so to see the exact form for each of these kernel functions, simply
  enter the name of this function in R (i.e., enter
  \code{npuniden.boundary} after loading this package) and scroll up for
  their definitions.

  The kernel functions \code{"gamma"}, \code{"rigaussian"}, and
  \code{"fbl"} have support \eqn{[a,\infty]}. The kernel function
  \code{"fbu"} has support \eqn{[-\infty,b]}. The rest have support on
  \eqn{[a,b]}.  Note that the two sided support default values are
  \code{a=min(X)} and \code{b=max(X)}. 

  Note that data-driven bandwidth selection is more nuanced in bounded
  settings, therefore it would be prudent to manually select a bandwidth
  that is, say, 1/25th of the range of the data and manually inspect the
  estimate (say \code{h=0.05} when \eqn{X\in [0,1]}). Also, it may be
  wise to compare the density estimate with that from a histogram with
  the option \code{breaks=25}. Note also that the kernel functions
  \code{"gaussian2"}, \code{"fb"}, \code{"fbl"} and \code{"fbu"} can
  assume negative values leading to potentially negative density
  estimates, and must be trimmed when conducting likelihood
  cross-validation which can lead to oversmoothing. Least squares
  cross-validation is unaffected and appears to be more reliable in such
  instances hence is the default here.

  Scott (1992, Page 149) writes \dQuote{While boundary kernels can be
  very useful, there are potentially serious problems with real
  data. There are an infinite number of boundary kernels reflecting the
  spectrum of possible design constraints, and these kernels are not
  interchangeable. Severe artifacts can be introduced by any one of them
  in inappropriate situations. Very careful examination is required to
  avoid being victimized by the particular boundary kernel
  chosen. Artifacts can unfortunately be introduced by the choice of the
  support interval for the boundary kernel.}

  Note that since some kernel functions can assume negative values, this
  can lead to improper density estimates. The estimated distribution
  function is obtained via numerical integration of the estimated
  density function and may itself not be proper even when evaluated on
  the full range of the data \eqn{[a,b]}. Setting the option
  \code{proper=TRUE} will render the density and distribution estimates
  proper over the full range of the data, though this may not in
  general be a mean square error optimal strategy.

  Finally, note that this function is pretty bare-bones relative to
  other functions in this package. For one, at this time there is no
  automatic print support so kindly see the examples for illustrations
  of its use, among other differences.

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

  Igarashi, G. and Y. Kakizawa (2014). \dQuote{Re-formulation of inverse
  Gaussian, reciprocal inverse Gaussian, and Birnbaum-Saunders kernel
  estimators,} Statistics & Probability Letters, 84:235-246.

  Igarashi, G. and Y. Kakizawa (2015). \dQuote{Bias corrections for some
  asymmetric kernel estimators,} Journal of Statistical Planning and
  Inference, 159:37-63.
  
  Igarashi, G. (2016). \dQuote{Bias reductions for beta kernel
  estimation,} Journal of Nonparametric Statistics, 28(1):1-30.
  
  Racine, J. S. and Q. Li and Q. Wang, \dQuote{Boundary-adaptive 
  kernel density estimation: the case of (near) uniform density}, 
  Journal of Nonparametric Statistics, 2024, 36 (1), 146-164, 
  https://doi.org/10.1080/10485252.2023.2250011.

  Scaillet, O. (2004). \dQuote{Density estimation using inverse and
  reciprocal inverse Gaussian kernels,} Journal of Nonparametric
  Statistics, 16(1-2):217-226.

  Scott, D. W. (1992). \dQuote{Multivariate density estimation: Theory,
  practice, and visualization,} New York: Wiley.

  Zhang, S. and R. J. Karunamuni (2010). \dQuote{Boundary performance of
  the beta kernel estimators,} Journal of Nonparametric Statistics,
  22(1):81-104.
  
}
\author{
Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\seealso{
The \pkg{Ake}, \pkg{bde}, and \pkg{Conake} packages and the function \code{\link{npuniden.reflect}}.
}
\examples{
\dontrun{
## Not run in checks: excluded to keep MPI examples stable and check times short.
## Example 1: f(0)=0, f(1)=1, plot boundary corrected density,
## unadjusted density, and DGP
set.seed(42)
n <- 100
X <- sort(rbeta(n,5,1))
dgp <- dbeta(X,5,1)
model.g1 <- npuniden.boundary(X,kertype="gaussian1")
model.g2 <- npuniden.boundary(X,kertype="gaussian2")
model.b1 <- npuniden.boundary(X,kertype="beta1")
model.b2 <- npuniden.boundary(X,kertype="beta2")
model.fb <- npuniden.boundary(X,kertype="fb")
model.unadjusted <- npuniden.boundary(X,a=-Inf,b=Inf)
ylim <- c(0,max(c(dgp,model.g1$f,model.g2$f,model.b1$f,model.b2$f,model.fb$f)))
plot(X,dgp,ylab="Density",ylim=ylim,type="l")
lines(X,model.g1$f,lty=2,col=2)
lines(X,model.g2$f,lty=3,col=3)
lines(X,model.b1$f,lty=4,col=4)
lines(X,model.b2$f,lty=5,col=5)
lines(X,model.fb$f,lty=6,col=6)
lines(X,model.unadjusted$f,lty=7,col=7)
rug(X)
legend("topleft",c("DGP",
                   "Boundary Kernel (gaussian1)",
                   "Boundary Kernel (gaussian2)",
                   "Boundary Kernel (beta1)",
                   "Boundary Kernel (beta2)",
                   "Boundary Kernel (floating boundary)",
                   "Unadjusted"),col=1:7,lty=1:7,bty="n")

## Example 2: f(0)=0, f(1)=0, plot density, distribution, DGP, and
## asymptotic point-wise confidence intervals
set.seed(42)
X <- sort(rbeta(100,5,3))
model <- npuniden.boundary(X)
oldpar <- par(no.readonly = TRUE)
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
hist(age,prob=TRUE,main="")
lines(age,model$f)
lines(density(age,bw=model$h),col=2)
legend("topright",c("Boundary Kernel","Unadjusted"),lty=c(1,1),col=1:2,bty="n")
detach(cps71)
par(oldpar)
} % enddontrun
}
\keyword{ nonparametric }
\keyword{ smooth }
