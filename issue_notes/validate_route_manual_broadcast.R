#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(npRmpi))

npRmpi.init(nslaves = 1, quiet = TRUE)
on.exit(try(npRmpi.quit(), silent = TRUE), add = TRUE)

options(npRmpi.autodispatch = FALSE, np.messages = FALSE)

set.seed(42)
n <- 120
x <- runif(n)
y <- rnorm(n)
d <- data.frame(x = x, y = y)

mpi.bcast.Robj2slave(d)
mpi.bcast.cmd(
  bw <- npregbw(y ~ x, data = d, regtype = "lc", bwmethod = "cv.ls", nmulti = 1),
  caller.execute = TRUE
)
mpi.bcast.cmd(
  fit <- npreg(bws = bw, data = d, gradients = FALSE),
  caller.execute = TRUE
)

stopifnot(inherits(fit, "npregression"))
cat("MANUAL_BCAST_ROUTE_OK\n")
