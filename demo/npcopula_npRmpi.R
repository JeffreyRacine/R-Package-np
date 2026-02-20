## Make sure you have the .Rprofile file from npRmpi/inst/ in your
## current directory or home directory. It is necessary.

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npudistml_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

## Initialize master and slaves.

npRmpi.start(nslaves=1)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)

options(npRmpi.autodispatch=TRUE, np.messages=FALSE)

## Generate some data and broadcast it to all slaves (it will be known
## to the master node)

set.seed(42)

library(MASS)

set.seed(42)

n <- 5000
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

## Clean up properly then quit()

npRmpi.stop(force=TRUE)
