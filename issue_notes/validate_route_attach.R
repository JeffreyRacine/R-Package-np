#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(npRmpi))

npRmpi.init(mode = "attach", quiet = TRUE)

if (mpi.comm.rank(1L) == 0L) {
  set.seed(42)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * x) + 0.5 * z + rnorm(n, sd = 0.1)
  d1 <- data.frame(x = x)
  d2 <- data.frame(x = x, z = z)
  dz <- data.frame(z = z)
  bw <- npregbw(y ~ x, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  fit <- npreg(bws = bw, gradients = FALSE)
  bw.sc <- npscoefbw(xdat = d1, ydat = y, zdat = dz, regtype = "lc", nmulti = 1)
  fit.sc <- npscoef(bws = bw.sc, gradients = FALSE)
  bw.pl <- npplregbw(xdat = d1, ydat = y, zdat = dz, regtype = "lc", nmulti = 1)
  fit.pl <- npplreg(bws = bw.pl, gradients = FALSE)
  bw.si <- npindexbw(xdat = d2, ydat = y, regtype = "lc", nmulti = 1)
  fit.si <- npindex(bws = bw.si, gradients = FALSE)
  stopifnot(inherits(fit, "npregression"))
  stopifnot(inherits(fit.sc, "smoothcoefficient"))
  stopifnot(inherits(fit.pl, "plregression"))
  stopifnot(inherits(fit.si, "singleindex"))
  cat("ATTACH_ROUTE_OK\n")
  npRmpi.quit(mode = "attach")
}
