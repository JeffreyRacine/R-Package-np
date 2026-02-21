## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npscoef_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

## Initialize master and slaves.

npRmpi.init(nslaves=1)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)

options(npRmpi.autodispatch=TRUE, np.messages=FALSE)

## Generate some data and broadcast it to all slaves (it will be known
## to the master node so no need to broadcast it)

set.seed(42)

n <- 10000

x <- runif(n)
z <- runif(n, min=-2, max=2)
y <- x*exp(z)*(1.0+rnorm(n,sd = 0.2))
mydat <- data.frame(x,y,z)
rm(x,y,z)
## A smooth coefficient model example

t <- system.time(bw <- npscoefbw(y~x|z,data=mydat))

summary(bw)

t <- t + system.time(model <- npscoef(bws=bw, gradients=TRUE))

summary(model)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.quit(force=TRUE)
