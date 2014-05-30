## Make sure you have the .Rprofile file from npRmpi/inst/ in your
## current directory or home directory. It is necessary.

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npreglcaic_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

## Initialize master and slaves.

mpi.bcast.cmd(np.mpi.initialize(),
              caller.execute=TRUE)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)

mpi.bcast.cmd(options(np.messages=FALSE),
              caller.execute=TRUE)

## Generate data and broadcast it to all slave nodes in the form of a
## data frame

set.seed(42)

n <- 5000

x <- runif(n)
z1 <- rbinom(n,1,.5)
z2 <- rbinom(n,1,.5)
y <- cos(2*pi*x) + z1 + rnorm(n,sd=.25)
z1 <- factor(z1)
z2 <- factor(z2)
mydat <- data.frame(y,x,z1,z2)
rm(x,y,z1,z2)

mpi.bcast.Robj2slave(mydat)

## A regression example (local constant, aic cross-validation)

t <- system.time(mpi.bcast.cmd(bw <- npregbw(y~x+z1+z2,
                                             regtype="lc",
                                             bwmethod="cv.aic",
                                             data=mydat),
                               caller.execute=TRUE))

summary(bw)

t <- t + system.time(mpi.bcast.cmd(model <- npreg(bws=bw,
                                                  data=mydat),
                                   caller.execute=TRUE))

summary(model)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
