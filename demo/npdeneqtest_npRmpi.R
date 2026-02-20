## Make sure you have the .Rprofile file from npRmpi/inst/ in your
## current directory or home directory. It is necessary.

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npdeneqtest_npRmpi. Check the time in the
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

n <- 2500

sample.A <- data.frame(x=rnorm(n))
sample.B <- data.frame(x=rnorm(n))
## A consistent density equality test example

t <- system.time(output <- npdeneqtest(sample.A,sample.B,boot.num=99))

output

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.stop(force=TRUE)
