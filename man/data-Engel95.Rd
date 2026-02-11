\name{Engel95}
\docType{data}
\alias{Engel95}
\title{ 1995 British Family Expenditure Survey  }
\description{
British cross-section data consisting of a random sample taken from
the British Family Expenditure Survey for 1995. The households consist
of married couples with an employed head-of-household between the ages
of 25 and 55 years. There are 1655 household-level observations in
total.
}
\usage{data("Engel95")}
\format{ A data frame with 10 columns, and 1655 rows.
  \describe{
 \item{food}{ expenditure share on food, of type \code{numeric}}
 \item{catering}{ expenditure share on catering, of type \code{numeric}}
 \item{alcohol}{ expenditure share on alcohol, of type \code{numeric}}
 \item{fuel}{ expenditure share on fuel, of type \code{numeric}}
 \item{motor}{ expenditure share on motor, of type \code{numeric}}
 \item{fares}{ expenditure share on fares, of type \code{numeric}}
 \item{leisure}{ expenditure share on leisure, of type \code{numeric}}
 \item{logexp}{ logarithm of total expenditure, of type \code{numeric}}
 \item{logwages}{ logarithm of total earnings, of type \code{numeric}}
 \item{nkids}{ number of children, of type \code{numeric}}  
  }
}
\source{ Richard Blundell and Dennis Kristensen }
\references{

  Blundell, R. and X. Chen and D. Kristensen (2007),
  \dQuote{Semi-Nonparametric IV Estimation of Shape-Invariant Engel
  Curves,} Econometrica, 75, 1613-1669.
    
  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.
  
}

\keyword{datasets}

\examples{
\dontrun{
## Not run in checks: this IV example is computationally expensive and can
## exceed check time limits in MPI environments.
## Example - compute nonparametric instrumental regression using
## Landweber-Fridman iteration of Fredholm integral equations of the
## first kind.

## We consider an equation with an endogenous regressor (`z') and an
## instrument (`w'). Let y = phi(z) + u where phi(z) is the function of
## interest. Here E(u|z) is not zero hence the conditional mean E(y|z)
## does not coincide with the function of interest, but if there exists
## an instrument w such that E(u|w) = 0, then we can recover the
## function of interest by solving an ill-posed inverse problem.

## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave). Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

## Start npRmpi for interactive execution. If slaves are already running and
## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
## reuse the existing pool instead of respawning. To change the number of
## slaves, call `npRmpi.stop(force=TRUE)` then restart.
npRmpi.start(nslaves=1)

data(Engel95)

## Sort on logexp (the endogenous regressor) for plotting purposes

Engel95 <- Engel95[order(Engel95$logexp),] 
mpi.bcast.Robj2slave(Engel95)

mpi.bcast.cmd(attach(Engel95),
              caller.execute=TRUE)

mpi.bcast.cmd(model.iv <- npregiv(y=food,z=logexp,w=logwages,method="Landweber-Fridman"),
              caller.execute=TRUE)
phi <- model.iv$phi

## Compute the non-IV regression (i.e. regress y on z)

mpi.bcast.cmd(ghat <- npreg(food~logexp,regtype="ll"),
              caller.execute=TRUE)

## For the plots, restrict focal attention to the bulk of the data
## (i.e. for the plotting area trim out 1/4 of one percent from each
## tail of y and z)

trim <- 0.0025

plot(logexp,food,
     ylab="Food Budget Share",
     xlab="log(Total Expenditure)",
     xlim=quantile(logexp,c(trim,1-trim)),
     ylim=quantile(food,c(trim,1-trim)),
     main="Nonparametric Instrumental Kernel Regression",
     type="p",
     cex=.5,
     col="lightgrey")

lines(logexp,phi,col="blue",lwd=2,lty=2)

lines(logexp,fitted(ghat),col="red",lwd=2,lty=4)

legend(quantile(logexp,trim),quantile(food,1-trim),
       c(expression(paste("Nonparametric IV: ",hat(varphi)(logexp))),
         "Nonparametric Regression: E(food | logexp)"),
       lty=c(2,4),
       col=c("blue","red"),
       lwd=c(2,2))

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
