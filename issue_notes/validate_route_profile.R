## Profile/manual-broadcast route smoke script.
## Launch with mpiexec and npRmpi profile loaded via R_PROFILE_USER.

mpi.bcast.cmd(np.mpi.initialize(),
              caller.execute = TRUE)

mpi.bcast.cmd(options(np.messages = FALSE),
              caller.execute = TRUE)

set.seed(42)

n <- as.integer(Sys.getenv("NP_RMPI_PROFILE_N", "120"))
x <- runif(n)
z1 <- rbinom(n, 1, .5)
z2 <- rbinom(n, 1, .5)
y <- cos(2 * pi * x) + z1 + rnorm(n, sd = .25)
z1 <- factor(z1)
z2 <- factor(z2)
mydat <- data.frame(y, x, z1, z2)
rm(x, y, z1, z2)

mpi.bcast.Robj2slave(mydat)

t <- system.time(mpi.bcast.cmd(
  bw <- npregbw(y ~ x + z1 + z2,
                regtype = "lc",
                bwmethod = "cv.ls",
                data = mydat),
  caller.execute = TRUE
))

t <- t + system.time(mpi.bcast.cmd(
  fit <- npreg(bws = bw, data = mydat),
  caller.execute = TRUE
))

stopifnot(inherits(fit, "npregression"))
cat("PROFILE_ROUTE_OK\n")
cat("Elapsed time =", t[3], "\n")

mpi.bcast.cmd(mpi.quit(),
              caller.execute = TRUE)
