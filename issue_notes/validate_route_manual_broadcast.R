#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(npRmpi))

npRmpi.init(nslaves = 1, quiet = TRUE)
on.exit(try(npRmpi.quit(), silent = TRUE), add = TRUE)

options(npRmpi.autodispatch = FALSE, np.messages = FALSE)

set.seed(42)
n <- 120
x <- runif(n)
z <- runif(n)
y <- sin(2 * pi * x) + 0.5 * z + rnorm(n, sd = 0.1)
d <- data.frame(x = x, y = y)
d1 <- data.frame(x = x)
d2 <- data.frame(x = x, z = z)
dz <- data.frame(z = z)

mpi.bcast.Robj2slave(d)
mpi.bcast.Robj2slave(d1)
mpi.bcast.Robj2slave(d2)
mpi.bcast.Robj2slave(dz)
mpi.bcast.Robj2slave(y)
mpi.bcast.cmd(
  bw <- npregbw(y ~ x, data = d, regtype = "lc", bwmethod = "cv.ls", nmulti = 1),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  fit <- npreg(bws = bw, data = d, gradients = FALSE),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  bw.sc <- npscoefbw(xdat = d1, ydat = y, zdat = dz, regtype = "lc", nmulti = 1),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  fit.sc <- npscoef(bws = bw.sc, gradients = FALSE),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  bw.pl <- npplregbw(xdat = d1, ydat = y, zdat = dz, regtype = "lc", nmulti = 1),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  fit.pl <- npplreg(bws = bw.pl, gradients = FALSE),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  bw.si <- npindexbw(xdat = d2, ydat = y, regtype = "lc", nmulti = 1),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  fit.si <- npindex(bws = bw.si, gradients = FALSE),
  caller.execute = TRUE
)

stopifnot(inherits(fit, "npregression"))
stopifnot(inherits(fit.sc, "smoothcoefficient"))
stopifnot(inherits(fit.pl, "plregression"))
stopifnot(inherits(fit.si, "singleindex"))
cat("MANUAL_BCAST_ROUTE_OK\n")
