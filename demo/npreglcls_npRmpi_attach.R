library(npRmpi)

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npreglcls_npRmpi. Check the time in the
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
## Generate data and broadcast it to all slave nodes in the form of a
## data frame

set.seed(42)

n <- as.integer(Sys.getenv("NP_DEMO_N", "5000"))
x <- runif(n)
z1 <- rbinom(n,1,.5)
z2 <- rbinom(n,1,.5)
y <- cos(2*pi*x) + z1 + rnorm(n,sd=.25)
z1 <- factor(z1)
z2 <- factor(z2)
mydat <- data.frame(y,x,z1,z2)
rm(x,y,z1,z2)
## A regression example (local constant, least-squares cross-validation)

t <- system.time(bw <- npregbw(y~x+z1+z2,
                                             regtype="lc",
                                             bwmethod="cv.ls",
                                             data=mydat))

summary(bw)

t <- t + system.time(model <- npreg(bws=bw,
                                                  data=mydat))

summary(model)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.quit(mode="attach", comm=1)
mpi.quit()
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
