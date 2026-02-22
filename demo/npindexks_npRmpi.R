## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npindexks_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

## Initialize master and slaves.

## Batch/cluster usage (attach mode under mpiexec):
##   mpiexec -n <master+slaves> R CMD BATCH --vanilla <script>.R
## Inside the script, use attach mode instead of spawning:
##   try(mpi.comm.dup(0, 1), silent = TRUE)
##   npRmpi.init(mode="attach", comm=1, autodispatch=TRUE, np.messages=FALSE)
##
npRmpi.init(nslaves=1)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)

options(npRmpi.autodispatch=TRUE, np.messages=FALSE)

## Generate some data and broadcast it to all slaves (it will be known
## to the master node so no need to broadcast it)

set.seed(42)

n <- 5000

x <- rchisq(n, df=3)
x1 <- (ifelse(x < 6, x, 6) - 2.348)/1.511
x <- rnorm(n)
x2 <- ifelse(abs(x) < 2 , x, 2) / 0.8796
y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)
mydat <- data.frame(x1,x2,y)
rm(x,x1,x2,y)
## A single index model example (Klein & Spady, binary y)

t <- system.time(bw <- npindexbw(formula=y~x1+x2,
                                               method="kleinspady",
                                               data=mydat))

summary(bw)

t <- t + system.time(model <- npindex(bws=bw, gradients=TRUE))

summary(model)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.quit(force=TRUE)
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
