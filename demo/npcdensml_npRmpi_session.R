library(npRmpi)
options(np.messages=FALSE)

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

library(MASS)

nslaves <- np_demo_n(default = 1L, floor = 1L,
                     exact_env = "NP_DEMO_NSLAVES",
                     frac_env = "NP_DEMO_NSLAVES_FRAC")
npRmpi.init(nslaves = nslaves)

set.seed(42)

default_n <- 2500L
n <- np_demo_n(default_n)
rho <- 0.25
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
data <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=data[,2],y=data[,1])

t <- system.time(bw <- npcdensbw(y~x,
                                 bwmethod="cv.ml",
                                 data=mydat))

summary(bw)

t <- t + system.time(model <- npcdens(bws=bw))

summary(model)

cat("Elapsed time =", t[3], "\n")
np_demo_result("npcdensml", "session", n, default_n, t[3],
               bwmethod = "cv.ml", slaves = nslaves)

## The demo harness cleans spawned slave daemons from the launcher layer.
