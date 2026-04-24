library(npRmpi)

test_that("public npcdensbw cv.ml keeps lc adjacency live during containment", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260306)
  n <- 30L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + rnorm(n, sd = 0.15))

  bw.lc <- npcdensbw(xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ml", nmulti = 1)
  fit.lc <- npcdens(bws = bw.lc)

  expect_true(is.finite(bw.lc$fval))
  expect_true(all(is.finite(fitted(fit.lc))))
})

test_that("public npcdensbw cv.ml fixed LP/LL route activates with ll == lp parity", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260307)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 - x$x2 + rnorm(n, sd = 0.2))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwmethod = "cv.ml",
    nmulti = 1
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ml",
    nmulti = 1
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdensbw cv.ml generalized-nn LP route activates with ll == lp parity", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260308)
  n <- 24L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 1,
    itmax = 1
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 1,
    itmax = 1
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdensbw cv.ml adaptive-nn LP route activates with ll == lp parity", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260309)
  n <- 24L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ml",
    nmulti = 1,
    itmax = 1
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ml",
    nmulti = 1,
    itmax = 1
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})
