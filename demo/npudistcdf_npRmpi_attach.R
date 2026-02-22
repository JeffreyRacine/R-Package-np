library(npRmpi)

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npudistml_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

## Initialize master and slaves.

## Batch/cluster usage (attach mode under mpiexec):
##   mpiexec -n <master+slaves> R CMD BATCH --vanilla <script>.R
## Inside the script, use attach mode instead of spawning:
##   try(mpi.comm.dup(0, 1), silent = TRUE)
##   npRmpi.init(mode="attach", comm=1, autodispatch=TRUE, np.messages=FALSE)
##
npRmpi.init(mode="attach", comm=1, autodispatch=TRUE)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)
## Generate some data and broadcast it to all slaves (it will be known
## to the master node)

set.seed(42)

n <- as.integer(Sys.getenv("NP_DEMO_N", "10000"))
mydat <- data.frame(x=rnorm(n))
## A simple example with cross-validation

t <- system.time(bw <- npudistbw(~x,
                                               bwmethod="cv.cdf",
                                               do.full.integral = TRUE,
                                               data=mydat))

summary(bw)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.quit(mode="attach", comm=1)
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
