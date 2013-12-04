## $Id: lc_k1_mean.R,v 1.11 2013/12/01 16:07:16 jracine Exp jracine $

## We illustrate constrained kernel estimation using the local
## constant estimator where the constraints are l(x) <= \hat g(x) <=
## u(x) (see Du, P. and C. Parmeter and J.S. Racine (2013),
## "Nonparametric Kernel Regression with Multiple Predictors and
## Multiple Shape Constraints," Statistica Sinica, Volume 23, Number
## 3, 1343-1372).

rm(list=ls())

## Load required packages, set options.

library(np)
library(quadprog)
options(np.tree=TRUE,np.messages=FALSE)

## Set the kernel function.

ckertype <- "epanechnikov"

## Simulate a sample of data. We can control signal/noise ratio by
## multiplying epsilon by sd(dgp), so rnorm(n,sd=...) can be set to
## sd=(.25,.5,1,2) which would yield an R-squared for the Oracle model
## of (.95,.8,.5,and .2).

n <- 1000
x <- sort(runif(n))
dgp <- cos(2*pi*x)
y <- dgp + sd(dgp)*rnorm(n,sd=1)

## X (data frame of regressors) and y are passed below, so if you add
## extra regressors simply add them to X here and be done.

X <- data.frame(x)

## Set up bounds for the quadratic program. We are going to require
## the lower and upper constraints l(x) and u(x) for g(x). Here they
## are constant, but in general can depend on x.

lower <- rep(-0.5,n)
upper <- rep(0.5,n)

## Generate the matrix of kernel weights using data-driven bandwidths
## that are optimal for the unconstrained model.

K <- npksum(txdat=X,
            bws=npregbw(xdat=X,ydat=y,ckertype=ckertype)$bw,
            ckertype=ckertype,
            return.kernel.weights=TRUE)$kw

## Compute the matrix for which t(A)%*%p is the local constant
## estimator

A <- K*y/colSums(K)

## Solve the quadratic program. The function solve.QP in the quadprog
## package solves the problem min (p-p.u)'(p-p.u) subject to the
## constraints Amat^T p >= bvec. Note that we construct Amat to
## contain the constraints a) ) A^T p >= lower, and b) -A^Tp >= -upper
## (here A^T = \hat g(x|p)).

p.hat <- solve.QP(Dmat=diag(n),
                  dvec=rep(1,n),
                  Amat=cbind(A,-A),
                  bvec=c(lower,-upper))$solution

## Compute the unconstrained and constrained estimator.

p.u <- rep(1,n)
fit.unres <- t(A)%*%p.u
fit.res <- t(A)%*%p.hat

## Plot the DGP, unrestricted, and restricted fits along with the
## constraints and data.

ylim <- c(min(dgp,fit.unres,fit.res,y),
          max(dgp,fit.unres,fit.res,y))

plot(x,y,cex=.25,ylab="g(x)",ylim=ylim,col="grey")
lines(x,dgp,col=1,lty=1)
lines(x,fit.unres,col=2,lty=1,lwd=2)
lines(x,fit.res,col=3,lty=1,lwd=2)
lines(x,lower,lty=1,col="lightgrey")
lines(x,upper,lty=1,col="lightgrey")

legend("bottomleft",
       c("DGP","Unconstrained","Constrained"),
       lty=c(1,1,1),
       col=1:3,
       bty="n",
       cex=0.75)
