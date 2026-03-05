## Profile/manual-broadcast route smoke script.
## Launch with mpiexec and npRmpi profile loaded via R_PROFILE_USER.

mpi.bcast.cmd(np.mpi.initialize(),
              caller.execute = TRUE)

mpi.bcast.cmd(options(np.messages = FALSE, npRmpi.autodispatch = FALSE),
              caller.execute = TRUE)

suppressPackageStartupMessages(library(MASS))
set.seed(42)

n <- as.integer(Sys.getenv("NP_RMPI_PROFILE_N", "120"))
x <- runif(n)
zc <- runif(n)
z1 <- rbinom(n, 1, .5)
z2 <- rbinom(n, 1, .5)
y <- cos(2 * pi * x) + z1 + rnorm(n, sd = .25)
z1 <- factor(z1)
z2 <- factor(z2)
mydat <- data.frame(y, x, zc, z1, z2)
d.cop <- data.frame(x = mydat$x, y = mydat$y)
u.cop <- data.frame(x = c(0.25, 0.5, 0.75), y = c(0.25, 0.5, 0.75))
data(birthwt)
bdat <- birthwt
bdat$low <- factor(bdat$low)
bdat$smoke <- factor(bdat$smoke)
bdat$race <- factor(bdat$race)
bdat$ht <- factor(bdat$ht)
bdat$ui <- factor(bdat$ui)
bdat$ftv <- ordered(bdat$ftv)
rm(x, y, zc, z1, z2)

mpi.bcast.Robj2slave(mydat)
mpi.bcast.Robj2slave(d.cop)
mpi.bcast.Robj2slave(u.cop)
mpi.bcast.Robj2slave(bdat)

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

mpi.bcast.cmd(
  bw.sc <- npscoefbw(y ~ x | zc, regtype = "lc", nmulti = 1, data = mydat),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  fit.sc <- npscoef(bws = bw.sc, data = mydat, gradients = FALSE),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  bw.pl <- npplregbw(y ~ x | zc, regtype = "lc", nmulti = 1, data = mydat),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  fit.pl <- npplreg(bws = bw.pl, data = mydat, gradients = FALSE),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  bw.si <- npindexbw(y ~ x + zc, regtype = "lc", nmulti = 1, data = mydat),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  fit.si <- npindex(bws = bw.si, data = mydat, gradients = FALSE),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  bw.cop <- npudistbw(~x + y, data = d.cop),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  cop <- npcopula(bws = bw.cop, data = d.cop, u = u.cop, n.quasi.inv = 60),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  bw.cm <- npcdensbw(low ~ smoke + race + ht + ui + ftv + age + lwt, data = bdat, nmulti = 1),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  fit.cm <- npconmode(bws = bw.cm),
  caller.execute = TRUE
)

stopifnot(inherits(fit, "npregression"))
stopifnot(inherits(fit.sc, "smoothcoefficient"))
stopifnot(inherits(fit.pl, "plregression"))
stopifnot(inherits(fit.si, "singleindex"))
stopifnot(inherits(cop, "data.frame"))
stopifnot(nrow(cop) == 9L)
stopifnot(all(is.finite(cop$copula)))
stopifnot(inherits(fit.cm, "conmode"))
cat("PROFILE_NPCONMODE_ROUTE_OK\n")
cat("PROFILE_NPCOPULA_ROUTE_OK\n")
cat("PROFILE_ROUTE_OK\n")
cat("Elapsed time =", t[3], "\n")

mpi.bcast.cmd(mpi.quit(),
              caller.execute = TRUE)
