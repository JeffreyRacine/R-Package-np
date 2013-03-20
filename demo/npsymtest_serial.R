## This is the serial version of npsymtest_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

## Generate some data

set.seed(42)

## A function to create a time series

ar.series <- function(phi,epsilon) {
  m <- length(epsilon)
  series <- numeric(m)
  series[1] <- epsilon[1]/(1-phi)
  for(i in 2:m) {
    series[i] <- phi*series[i-1] + epsilon[i]
  }
  return(series)
}

n <- 2500

## Stationary persistent symmetric time-series

yt <- ar.series(0.5,rnorm(n))

## A simple example of the test for symmetry

t <- system.time(output <- npsymtest(yt,
                                     boot.num=399,
                                     boot.method="geom",
                                     method="summation"))

output

cat("Elapsed time =", t[3], "\n")
