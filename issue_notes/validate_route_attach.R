#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(npRmpi))

npRmpi.init(mode = "attach", quiet = TRUE)

if (mpi.comm.rank(1L) == 0L) {
  set.seed(42)
  n <- 120
  x <- runif(n)
  y <- rnorm(n)
  bw <- npregbw(y ~ x, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  fit <- npreg(bws = bw, gradients = FALSE)
  stopifnot(inherits(fit, "npregression"))
  cat("ATTACH_ROUTE_OK\n")
  npRmpi.quit(mode = "attach")
}
