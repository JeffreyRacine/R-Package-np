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

mpi.bcast.cmd(set.seed(42),
              caller.execute=TRUE)

## Significance testing with z irrelevant

n <- as.integer(Sys.getenv("NP_DEMO_N", "1000"))
z <- factor(rbinom(n,1,.5))
x1 <- rnorm(n)
x2 <- runif(n,-2,2)
y <- x1 + x2 + rnorm(n)
mydat <- data.frame(z,x1,x2,y)
rm(z,x1,x2,y)

mpi.bcast.Robj2slave(mydat)

t <- system.time(mpi.bcast.cmd(model <- npreg(y~z+x1+x2,
                                              regtype="ll",
                                              bwmethod="cv.aic",
                                              data=mydat),
                               caller.execute=TRUE))

## An example of the consistent nonparametric significance test

t <- t + system.time(mpi.bcast.cmd(output <- npsigtest(model),
                               caller.execute=TRUE ))

output

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
