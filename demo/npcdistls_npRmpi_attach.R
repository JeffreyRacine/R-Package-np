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
## Load your data and broadcast it to all slave nodes

library(MASS)

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

set.seed(42)

default_n <- 2000L
n <- np_demo_n(default_n)
rho <- 0.25
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
data <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=data[,2],y=data[,1])
## A conditional distribution estimation example. 

t <- system.time(bw <- npcdistbw(y~x,
                                               bwmethod="cv.ls",
                                               data=mydat))

summary(bw)

t <- t + system.time(model <- npcdist(bws=bw))

summary(model)

cat("Elapsed time =", t[3], "\n")
np_demo_result("npcdistls", "attach", n, default_n, t[3], bwmethod = "cv.ls")

## Clean up properly then quit()

npRmpi.quit(mode="attach", comm=1)
}
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
