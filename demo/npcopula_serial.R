## This is the serial version of npudistml_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
library(MASS)
options(np.messages=FALSE)

## Generate some data

set.seed(42)

require(MASS)

set.seed(42)

n <- 2500
n.eval <- 25
rho <- 0.95
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
mydat <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=mydat[,1],y=mydat[,2])

q.minm <- 0.0
q.max <- 1.0
grid.seq <- seq(q.minm,q.max,length=n.eval)
grid.dat <- cbind(grid.seq,grid.seq)

## Estimate the copula

t.0 <- system.time(bw <- npudistbw(~x+y,data=mydat))
t.1 <- system.time(copula <- npcopula(bws=bw,data=mydat,u=grid.dat))

t <- t.0+t.1

summary(bw)

cat("Elapsed time =", t[3], "\n")
