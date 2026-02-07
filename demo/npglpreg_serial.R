## This is the serial version of npglpreg_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
library(crs)
options(np.messages=FALSE,crs.messages=FALSE)

set.seed(42)
n <-  1500
x1 <- runif(n)
x2 <- runif(n)
dgp <- cos(8*pi*x1)
y <- dgp+rnorm(n,sd=0.1)
  
t <- system.time(model.glp <- npglpreg(y~x1+x2))

summary(model.glp)

cat("Elapsed time =", t[3], "\n")

