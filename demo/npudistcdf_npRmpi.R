## Make sure you have the .Rprofile file from npRmpi/inst/ in your
## current directory or home directory. It is necessary.

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npudistml_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

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

n <- 10000

mydat <- data.frame(x=rnorm(n))
mpi.bcast.Robj2slave(mydat)

## A simple example with cross-validation

t <- system.time(mpi.bcast.cmd(bw <- npudistbw(~x,
                                               bwmethod="cv.cdf",
                                               data=mydat),
                               caller.execute=TRUE))

summary(bw)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
