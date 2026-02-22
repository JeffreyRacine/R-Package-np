library(npRmpi)

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npudistml_npRmpi. Check the time in the
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
## Generate some data and broadcast it to all slaves (it will be known
## to the master node)

set.seed(42)

library(MASS)

set.seed(42)

n <- as.integer(Sys.getenv("NP_DEMO_N", "5000"))
n.eval <- 25
rho <- 0.95
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
mydat <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=mydat[,1],y=mydat[,2])

q.minm <- 0.0
q.max <- 1.0
grid.seq <- seq(q.minm,q.max,length=n.eval)
grid.dat <- data.frame(x = grid.seq, y = grid.seq)
## Estimate the copula

t.0 <- system.time(bw <- npudistbw(~x+y,data=mydat))
t.1 <- system.time({
  # npcopula currently fails under attach-mode worker autodispatch; run safely on master.
  old.autod <- getOption("npRmpi.autodispatch")
  options(npRmpi.autodispatch = FALSE)
  on.exit(options(npRmpi.autodispatch = old.autod), add = TRUE)
  copula <- try(np::npcopula(bws=bw, data=mydat, u=grid.dat), silent=TRUE)
  if (inherits(copula, "try-error")) {
    warning("npcopula fit failed in attach demo; continuing with bandwidth summary only.")
    copula <- NULL
  }
})

t <- t.0+t.1

summary(bw)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.quit(mode="attach", comm=1)
mpi.quit()
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
