## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npsymtest_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

## Initialize master and slaves.

npRmpi.init(nslaves=1)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)

options(npRmpi.autodispatch=TRUE, np.messages=FALSE)

## Generate some data and broadcast it to all slaves (it will be known
## to the master node)

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

## Clean up properly then quit()

npRmpi.quit(force=TRUE)
