## This is the serial version of npcdensls_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

library(MASS)

.np_demo_src <- Sys.getenv("NP_DEMO_SRC", "")
.np_demo_utils <- c(if (nzchar(.np_demo_src)) file.path(.np_demo_src, "..", "inst", "demo_utils.R"),
                    if (nzchar(.np_demo_src)) file.path(.np_demo_src, "demo_utils.R"),
                    "demo_utils", "demo_utils.R", "../demo_utils", "../demo_utils.R",
                    file.path("demo", "demo_utils"), file.path("demo", "demo_utils.R"),
                    system.file("demo_utils.R", package = "npRmpi"))
.np_demo_utils <- .np_demo_utils[file.exists(.np_demo_utils)]
source(.np_demo_utils[[1L]])

set.seed(42)

default_n <- 1000L
n <- np_demo_n(default_n)
rho <- 0.25
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
data <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=data[,2],y=data[,1])

## A simple example with least-squares cross-validation

t <- system.time(bw <- npcdensbw(y~x,
                                 bwmethod="cv.ls",
                                 data=mydat))

summary(bw)

t <- t + system.time(model <- npcdens(bws=bw))

summary(model)

cat("Elapsed time =", t[3], "\n")
np_demo_result("npcdensls", "serial", n, default_n, t[3], bwmethod = "cv.ls")
