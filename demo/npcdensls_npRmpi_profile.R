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

## Load your data and broadcast it to all slave nodes

mpi.bcast.cmd(library(MASS),
              caller.execute=TRUE)

mpi.bcast.cmd(set.seed(42),
              caller.execute=TRUE)

n <- as.integer(Sys.getenv("NP_DEMO_N", "1000"))
rho <- 0.25
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
data <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=data[,2],y=data[,1])

mpi.bcast.Robj2slave(mydat)

## A conditional density estimation example. 

t <- system.time(mpi.bcast.cmd(bw <- npcdensbw(y~x,
                                               bwmethod="cv.ls",
                                               data=mydat),
                               caller.execute=TRUE))

summary(bw)

t <- t + system.time(mpi.bcast.cmd(model <- npcdens(bws=bw),
                                   caller.execute=TRUE))

summary(model)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
