## This is the serial version of npdeptest_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

## Generate some data

set.seed(42)

n <- 2500

x <- rnorm(n)
y <- 1 + x + rnorm(n)
model <- lm(y~x)
y.fit <- fitted(model)
     

## A simple example for the consistent dependence metric test 

t <- system.time(output <- npdeptest(y,y.fit,
                                     boot.num=99,
                                     method="summation"))
                 
output

cat("Elapsed time =", t[3], "\n")
