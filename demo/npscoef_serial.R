## This is the serial version of npscoef_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

## Generate some data

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

