## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npunitest_npRmpi. Check the time in the
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

n <- 5000

x <- rnorm(n)
y <- rnorm(n)
## A simple example of the test for equality of univariate densities

t <- system.time(output <- npunitest(x,y,
                                                   method="summation",
                                                   bootstrap=TRUE))

output

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.quit(force=TRUE)
