## Profile/manual-broadcast route smoke script.
## Launch with mpiexec and npRmpi profile loaded via R_PROFILE_USER.

mpi.bcast.cmd(np.mpi.initialize(), caller.execute = TRUE)
mpi.bcast.cmd(options(np.messages = FALSE), caller.execute = TRUE)

set.seed(42)
n <- as.integer(Sys.getenv("NP_RMPI_PROFILE_N", "120"))
x <- runif(n)
y <- rnorm(n)
d <- data.frame(y = y, x = x)

mpi.bcast.Robj2slave(d)

t <- system.time(
  mpi.bcast.cmd(
    bw <- npregbw(y ~ x, data = d, regtype = "lc", bwmethod = "cv.ls", nmulti = 1),
    caller.execute = TRUE
  )
)
t <- t + system.time(
  mpi.bcast.cmd(
    fit <- npreg(bws = bw, data = d, gradients = FALSE),
    caller.execute = TRUE
  )
)

stopifnot(inherits(fit, "npregression"))
cat("PROFILE_ROUTE_OK\n")
cat("Elapsed time =", t[3], "\n")

mpi.bcast.cmd(mpi.quit(), caller.execute = TRUE)
