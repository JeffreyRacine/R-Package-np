## Profile/manual-broadcast demo (mpiexec + .Rprofile + mpi.bcast.*).
## Run with two ranks (master + one worker), e.g.
##   mpiexec -env R_PROFILE_USER ../.Rprofile -env R_PROFILE "" \\
##           -n 2 R CMD BATCH --no-save <script>.R
## Do not use R CMD BATCH --vanilla for profile mode.
##
## Initialize master and slaves.

mpi.bcast.cmd(np.mpi.initialize(),
              caller.execute=TRUE)

mpi.bcast.cmd(options(np.messages=FALSE),
              caller.execute=TRUE)


## Generate data and broadcast it to all slave nodes in the form of a
## data frame

set.seed(42)

.np_demo_src <- Sys.getenv("NP_DEMO_SRC", "")
.np_demo_utils <- c(if (nzchar(.np_demo_src)) file.path(.np_demo_src, "..", "inst", "demo_utils.R"),
                    if (nzchar(.np_demo_src)) file.path(.np_demo_src, "demo_utils.R"),
                    "demo_utils", "demo_utils.R", "../demo_utils", "../demo_utils.R",
                    file.path("demo", "demo_utils"), file.path("demo", "demo_utils.R"),
                    system.file("demo_utils.R", package = "npRmpi"),
                    system.file("demo", "demo_utils", package = "npRmpi"),
                    system.file("demo", "demo_utils.R", package = "npRmpi"))
.np_demo_utils <- .np_demo_utils[nzchar(.np_demo_utils) & file.exists(.np_demo_utils)]
source(.np_demo_utils[[1L]])

default_n <- 5000L
n <- np_demo_n(default_n)

x <- runif(n)
z1 <- rbinom(n,1,.5)
z2 <- rbinom(n,1,.5)
y <- cos(2*pi*x) + z1 + rnorm(n,sd=.25)
z1 <- factor(z1)
z2 <- factor(z2)
mydat <- data.frame(y,x,z1,z2)
rm(x,y,z1,z2)

mpi.bcast.Robj2slave(mydat)

## A regression example (local linear, least squares cross-validation)

t <- system.time(mpi.bcast.cmd(bw <- npregbw(y~x+z1+z2,
                                             regtype="ll",
                                             bwmethod="cv.ls",
                                             data=mydat),
                               caller.execute=TRUE))

summary(bw)

t <- t + system.time(mpi.bcast.cmd(model <- npreg(bws=bw,
                                                  data=mydat),
                                   caller.execute=TRUE))

summary(model)

cat("Elapsed time =", t[3], "\n")
np_demo_result("npregllls", "profile", n, default_n, t[3],
               bwmethod = "cv.ls", regtype = "ll")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
