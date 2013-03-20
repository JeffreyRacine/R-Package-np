## Make sure you have the .Rprofile file from npRmpi/inst/ in your
## current directory or home directory. It is necessary.

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npglpreg_npRmpi. Check the time in the
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

set.seed(42)
n <-  1000

degree.max <- 20

mpi.bcast.cmd(library(crs),
              caller.execute=TRUE)

x1 <- runif(n)
x2 <- runif(n)
dgp <- cos(8*pi*x1)
y <- dgp+rnorm(n,sd=0.1)
  
mpi.bcast.Robj2slave(y)
mpi.bcast.Robj2slave(x1)
mpi.bcast.Robj2slave(x2)
mpi.bcast.Robj2slave(degree.max)

t <- system.time(mpi.bcast.cmd(model.glp <-
              npglpreg(y~x1+x2,degree.max=degree.max,mpi=TRUE),
              caller.execute=TRUE))

summary(model.glp)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
