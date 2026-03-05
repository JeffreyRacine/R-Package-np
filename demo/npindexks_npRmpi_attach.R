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

n <- as.integer(Sys.getenv("NP_DEMO_N", "5000"))
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

npRmpi.quit(mode="attach", comm=1)
}
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
