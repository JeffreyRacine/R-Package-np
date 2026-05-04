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
## A regression example (local constant, least-squares cross-validation)

t <- system.time(bw <- npregbw(y~x+z1+z2,
                                             regtype="lc",
                                             bwmethod="cv.ls",
                                             data=mydat))

summary(bw)

t <- t + system.time(model <- npreg(bws=bw,
                                                  data=mydat))

summary(model)

cat("Elapsed time =", t[3], "\n")
np_demo_result("npreglcls", "attach", n, default_n, t[3],
               bwmethod = "cv.ls", regtype = "lc")

## Clean up properly then quit()

npRmpi.quit(mode="attach", comm=1)
}
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
