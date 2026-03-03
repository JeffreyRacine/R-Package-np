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

data("Italy")



n <- 1008 # n <- nrow(Italy)


mpi.bcast.Robj2slave(Italy)

## A quantile regression example

t <- system.time(mpi.bcast.cmd(bw <- npcdistbw(gdp~ordered(year),data=Italy),
                               caller.execute=TRUE))

summary(bw)

t <- t + system.time(mpi.bcast.cmd(model.q0.25 <- npqreg(bws=bw, tau=0.25),
                                   caller.execute=TRUE))
t <- t + system.time(mpi.bcast.cmd(model.q0.50 <- npqreg(bws=bw, tau=0.50),
                                   caller.execute=TRUE))
t <- t + system.time(mpi.bcast.cmd(model.q0.75 <- npqreg(bws=bw, tau=0.75),
                                   caller.execute=TRUE))

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
