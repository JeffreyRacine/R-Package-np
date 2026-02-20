## Make sure you have the .Rprofile file from npRmpi/inst/ in your
## current directory or home directory. It is necessary.

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npudensml_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

## Initialize master and slaves.

npRmpi.start(nslaves=1)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)

options(npRmpi.autodispatch=TRUE, np.messages=FALSE)

## Generate some data and broadcast it to all slaves (it will be known
## to the master node)

set.seed(42)

n <- 10000

mydat <- data.frame(x=rnorm(n))
## A simple example with likelihood cross-validation

t <- system.time(bw <- npudensbw(~x,
                                               bwmethod="cv.ml",
                                               data=mydat))

summary(bw)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.stop(force=TRUE)
