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
                 Y = NULL,
                 h = NULL,
                 a = 0,
                 b = 1,
                 ...)
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
## Not run in checks: excluded to keep MPI examples stable and check times short.
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

## Start npRmpi for interactive execution. If slaves are already running and
## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
## reuse the existing pool instead of respawning. To change the number of
## slaves, call `npRmpi.stop(force=TRUE)` then restart.
npRmpi.start(nslaves=1)

## Example 1: f(0)=0, f(1)=1, plot boundary corrected density,
## unadjusted density, and DGP

mpi.bcast.cmd(set.seed(42),
              caller.execute=TRUE)

n <- 100
X <- sort(rbeta(n,5,1))
dgp <- dbeta(X,5,1)

mpi.bcast.Robj2slave(X)

mpi.bcast.cmd(model <- npuniden.reflect(X),
              caller.execute=TRUE)

mpi.bcast.cmd(model.unadjusted <- npuniden.boundary(X,a=-Inf,b=Inf),
              caller.execute=TRUE)

ylim <- c(0,max(c(dgp,model$f,model.unadjusted$f)))
plot(X,model$f,ylab="Density",ylim=ylim,type="l")
lines(X,model.unadjusted$f,lty=2,col=2)
lines(X,dgp,lty=3,col=3)
rug(X)
legend("topleft",c("Data-Reflection","Unadjusted","DGP"),col=1:3,lty=1:3,bty="n")

## Example 2: f(0)=0, f(1)=0, plot density, distribution, DGP, and
## asymptotic point-wise confidence intervals

mpi.bcast.cmd(set.seed(42),
              caller.execute=TRUE)

X <- sort(rbeta(100,5,3))

mpi.bcast.Robj2slave(X)

mpi.bcast.cmd(model <- npuniden.reflect(X),
              caller.execute=TRUE)

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

mpi.bcast.cmd(data(cps71),
              caller.execute=TRUE)

mpi.bcast.cmd(model <- npuniden.reflect(cps71$age,a=21,b=65),
              caller.execute=TRUE)

par(mfrow=c(1,1))
hist(cps71$age,prob=TRUE,main="",ylim=c(0,max(model$f)))
lines(cps71$age,model$f)
lines(density(cps71$age,bw=model$h),col=2)
legend("topright",c("Data-Reflection","Unadjusted"),lty=c(1,1),col=1:2,bty="n")

par(oldpar)

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

## Note: on some systems (notably macOS+MPICH), repeatedly spawning and
## tearing down slaves in the same R session can lead to hangs/crashes.
## npRmpi may therefore keep slave daemons alive by default and
## `npRmpi.stop()` performs a "soft close". Use `force=TRUE` to
## actually shut down the slaves.
##
## You can disable reuse via `options(npRmpi.reuse.slaves=FALSE)` or by
## setting the environment variable `NP_RMPI_NO_REUSE_SLAVES=1` before
## loading the package.

npRmpi.stop()               ## soft close (may keep slaves alive)
## npRmpi.stop(force=TRUE)  ## hard close

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)
} % enddontrun
}
\keyword{ nonparametric }
\keyword{ smooth }
