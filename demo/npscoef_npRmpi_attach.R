library(npRmpi)

## Attach mode demo (session routing under mpiexec).
## Run with two ranks (master + one worker), e.g.
##   mpiexec -n 2 Rscript --no-save <script>.R
## or
##   mpiexec -n 2 R CMD BATCH --no-save <script>.R
##
## If running under mpiexec, keep profile env vars cleared for attach mode:
##   -env R_PROFILE_USER "" -env R_PROFILE ""
##
## Initialize master and slaves.
npRmpi.init(mode="attach", comm=1, autodispatch=TRUE)
options(np.messages=FALSE)

if (mpi.comm.rank(0L) == 0L) {

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)
## Generate some data and broadcast it to all slaves (it will be known
## to the master node so no need to broadcast it)

set.seed(42)

n <- as.integer(Sys.getenv("NP_DEMO_N", "10000"))
x <- runif(n)
z <- runif(n, min=-2, max=2)
y <- x*exp(z)*(1.0+rnorm(n,sd = 0.2))
mydat <- data.frame(x,y,z)
rm(x,y,z)
## A smooth coefficient model example

t <- system.time(bw <- npscoefbw(y~x|z,data=mydat))

summary(bw)

t <- t + system.time(model <- npscoef(bws=bw, gradients=TRUE))

summary(model)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.quit(mode="attach", comm=1)
}
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
