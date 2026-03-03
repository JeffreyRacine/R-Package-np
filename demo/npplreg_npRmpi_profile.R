## Profile/manual-broadcast demo (mpiexec + .Rprofile + mpi.bcast.*).
## Run with two ranks (master + one worker), e.g.
##   mpiexec -env R_PROFILE_USER ../.Rprofile -env R_PROFILE "" \\
##           -n 2 R CMD BATCH --no-save <script>.R
## Do not use R CMD BATCH --vanilla for profile mode.
##
## Initialize master and slaves.

mpi.bcast.cmd(np.mpi.initialize(),
              caller.execute=TRUE)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)

mpi.bcast.cmd(options(np.messages=FALSE),
              caller.execute=TRUE)

## Generate some data and broadcast it to all slaves (it will be known
## to the master node)

mpi.bcast.cmd(set.seed(42),
              caller.execute=TRUE)

n <- as.integer(Sys.getenv("NP_DEMO_N", "1000"))
x1 <- rnorm(n)
x2 <- rbinom(n, 5, .3)

z1 <- rbinom(n, 2, .3)
z2 <- rnorm(n)

y <- 1 + x1 + x2 + z1 + sin(z2) + rnorm(n)

x2 <- ordered(x2)
z1 <- ordered(z1)

mpi.bcast.Robj2slave(x1)
mpi.bcast.Robj2slave(x2)
mpi.bcast.Robj2slave(z1)
mpi.bcast.Robj2slave(z2)
mpi.bcast.Robj2slave(y)

## Partially linear model

t <- system.time(mpi.bcast.cmd(bw <- npplregbw(formula=y~x1+x2|z1+z2),
                               caller.execute=TRUE))

summary(bw)

t <- t + system.time(mpi.bcast.cmd(pl <- npplreg(bws=bw),
                                   caller.execute=TRUE))

summary(pl)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
