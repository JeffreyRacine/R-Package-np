## This is the serial version of npunitest_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

## Generate some data

set.seed(42)

n <- 5000

x <- rnorm(n)
y <- rnorm(n)

## A simple example of the test for equality of univariate densities

t <- system.time(output <- npunitest(x,y,
                                     method="summation",
                                     bootstrap=TRUE))

output

cat("Elapsed time =", t[3], "\n")

