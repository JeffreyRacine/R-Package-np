## Profile/manual-broadcast demo (mpiexec + .Rprofile + mpi.bcast.*).
## Run with two ranks (master + one worker), e.g.
##   mpiexec -env R_PROFILE_USER ../.Rprofile -env R_PROFILE "" \\
##           -n 2 R CMD BATCH --no-save <script>.R
## Do not use R CMD BATCH --vanilla for profile mode.
##
## Initialize master and slaves.

mpi.bcast.cmd(np.mpi.initialize(),
              caller.execute=TRUE)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)

mpi.bcast.cmd(options(np.messages=FALSE),
              caller.execute=TRUE)

## Generate some data and broadcast it to all slaves (it will be known
## to the master node)

mpi.bcast.cmd(set.seed(42),
              caller.execute=TRUE)

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

mpi.bcast.Robj2slave(yt)

## A simple example of the test for symmetry

t <- system.time(mpi.bcast.cmd(output <- npsymtest(yt,
                                                   boot.num=399,
                                                   boot.method="geom",
                                                   method="summation"),
                                                   caller.execute=TRUE))

output

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
