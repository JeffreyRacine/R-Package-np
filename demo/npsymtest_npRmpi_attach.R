library(npRmpi)

## Attach mode demo (session routing under mpiexec).
## Run with two ranks (master + one worker), e.g.
##   mpiexec -n 2 Rscript --no-save <script>.R
## or
##   mpiexec -n 2 R CMD BATCH --no-save <script>.R
##
## If running under mpiexec, keep profile env vars cleared for attach mode:
##   -env R_PROFILE_USER "" -env R_PROFILE ""
##
## Initialize master and slaves.
npRmpi.init(mode="attach", comm=1, autodispatch=TRUE)
options(np.messages=FALSE)

if (mpi.comm.rank(0L) == 0L) {

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)
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

n <- as.integer(Sys.getenv("NP_DEMO_N", "2500"))
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

npRmpi.quit(mode="attach", comm=1)
}
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
