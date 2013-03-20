## This is the serial version of npindexich_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

## Generate some data

set.seed(42)

n <- 5000

x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)
y <- x1 - x2 + rnorm(n)
mydat <- data.frame(x1,x2,y)
rm(y,x1,x2)

## A single index model example (Ichimura, continuous y)

t <- system.time(bw <- npindexbw(formula=y~x1+x2,
                                 data=mydat))

summary(bw)

t <- t + system.time(model <- npindex(bws=bw,
                                      gradients=TRUE))

summary(model)

cat("Elapsed time =", t[3], "\n")

