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
library(crs)
options(np.tree=TRUE,np.messages=FALSE)

## Set the kernel function.

ckertype <- "epanechnikov"

## Set the order of the local polynomial `p' (p=0 is the
## Nadaraya-Watson estimator, p=1 local linear, p=2 local quadratic
## etc.). If p < 0 both p and the bandwidth (bw) are determined by
## cross-validation, otherwise only bw is determined.

p <- 1

## Simulate a sample of data. We can control signal/noise ratio by
## multiplying epsilon by sd(dgp), so rnorm(n,sd=...) can be set to
## sd=(.25,.5,1,2) which would yield an R-squared for the Oracle model
## of (.95,.8,.5,and .2). Changing n changes the sample size,
## modifying cos() the data generating process (dgp).

n <- 1000
x <- sort(runif(n))
dgp <- cos(2*pi*x)
y <- dgp + sd(dgp)*rnorm(n,sd=0.25)

## Set up bounds for the quadratic program. We are going to require
## the lower and upper constraints l(x) and u(x) for g(x). Here the
## upper bound is a constant, the lower bound depends on x.

lower <- -1.5+2*x
upper <- rep(0.5,n)
if(any(lower>upper)) stop(" constraints are inconsistent")

## X (data frame of regressors) and y are passed below, so if you add
## extra regressors simply add them to X here and to the formula (used
## to compute cross-validated smoothing parameters) and be done.

X <- data.frame(x)

formula.glp <- formula(y~x)

model.glp <- crs:::npglpreg(formula=formula.glp,
                            cv=ifelse(p>=0,"bandwidth","degree-bandwidth"),
                            degree=ifelse(p>=0,rep(p,NCOL(X)),rep(0,NCOL(X))),
                            ckertype=ckertype,
                            nmulti=min(NCOL(X),5))

bw <- model.glp$bws
p <- model.glp$degree

## W is the local polynomial model matrix (includes an intercept).

W <- crs:::W.glp(xdat=X,
                 degree=rep(p,NCOL(X)))

## Generate the matrix of kernel weights using data-driven bandwidths
## that are optimal for the unconstrained model.

K <- npksum(txdat=X,
            bws=bw,
            ckertype=ckertype,
            return.kernel.weights=TRUE)$kw

## Compute the matrix for which \hat g(x|p) = t(A)%*%p (i.e. the local
## polynomial fit).

A <- sapply(1:n,function(i){W[i,,drop=FALSE]%*%chol2inv(chol(t(W)%*%(K[,i]*W)))%*%t(W)*y*K[,i]})

## Compute the unconstrained estimator.

p.u <- rep(1,n)
fit.unres <- t(A)%*%p.u

## Solve the quadratic program. The function solve.QP in the quadprog
## package solves the problem min (p-p.u)'(p-p.u) subject to the
## constraints Amat^T p >= bvec. Note that we construct Amat to
## contain the constraints a) ) A^T p >= lower, and b) -A^Tp >= -upper
## (here A^T = \hat g(x|p)).

output.QP <- solve.QP(Dmat=diag(n),
                      dvec=rep(1,n),
                      Amat=cbind(A,-A),
                      bvec=c(lower,-upper))

if(is.nan(output.QP$value)) stop(" solve.QP was unable to find a solution: adjust constraints and restart")

p.hat <- output.QP$solution

## Compute the constrained estimator.

fit.res <- t(A)%*%p.hat

## Plot the DGP, unrestricted, and restricted fits along with the
## constraints and data.

ylim <- c(min(dgp,lower,fit.unres,fit.res,y),
          max(dgp,upper,fit.unres,fit.res,y))
subtext <- paste("Local polynomial order = ", p,", bandwidth = ",formatC(bw,digits=3,format="f"),", n = ",n,sep="")
plot(x,y,cex=.25,ylab="g(x)",ylim=ylim,col="grey",sub=subtext)
lines(x,dgp,col=1,lty=1)
lines(x,fit.unres,col=2,lty=1,lwd=2)
lines(x,fit.res,col=3,lty=1,lwd=2)
lines(x,lower,lty=1,col="lightgrey")
lines(x,upper,lty=1,col="lightgrey")

legend("topleft",
       c("DGP","Unconstrained [g(x)]","Constrained [g(x|p)]","Constraints [l(x) <= g(x|p) <= u(x)]"),
       lty=c(1,1,1,1),
       col=c(1:3,"lightgrey"),
       bty="n",
       cex=0.75)
