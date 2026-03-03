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

n <- as.integer(Sys.getenv("NP_DEMO_N", "5000"))
x <- rnorm(n)
y <- rnorm(n)

mpi.bcast.Robj2slave(x)
mpi.bcast.Robj2slave(y)

## A simple example of the test for equality of univariate densities

t <- system.time(mpi.bcast.cmd(output <- npunitest(x,y,
                                                   method="summation",
                                                   bootstrap=TRUE),
                               caller.execute=TRUE))

output

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
