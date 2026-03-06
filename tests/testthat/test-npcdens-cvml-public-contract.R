library(npRmpi)

test_that("public npcdensbw cv.ml enforces ll == lp(glp, degree=1)", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260306)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 - x$x2 + rnorm(n, sd = 0.2))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(xdat = x, ydat = y, regtype = "ll", bwmethod = "cv.ml", nmulti = 0)
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ml",
    nmulti = 0
  )

  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-10)
  expect_equal(bw.ll$xbw, bw.lp$xbw, tolerance = 1e-10)
  expect_equal(bw.ll$ybw, bw.lp$ybw, tolerance = 1e-10)
})

test_that("public npcdens fixed LP fit preserves tree parity at fixed bandwidths", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260307)
  n <- 34L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  old <- options(np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ml",
    nmulti = 0
  )

  options(np.tree = FALSE)
  fit.nt <- npcdens(bws = bw.lp)
  options(np.tree = TRUE)
  fit.tr <- npcdens(bws = bw.lp)

  expect_equal(fitted(fit.nt), fitted(fit.tr), tolerance = 1e-10)
})

test_that("public npcdens cv.ml LP route does not collapse gaussian bandwidths", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  n <- 100L
  x <- data.frame(x = rnorm(n))
  y <- data.frame(y = x$x + rnorm(n))

  bw.lc <- npcdensbw(xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ml", nmulti = 0)
  bw.ll <- npcdensbw(xdat = x, ydat = y, regtype = "ll", bwmethod = "cv.ml", nmulti = 0)
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 3L,
    bwmethod = "cv.ml",
    nmulti = 0
  )

  fit.ll <- npcdens(bws = bw.ll)
  fit.lp <- npcdens(bws = bw.lp)

  expect_true(all(is.finite(fitted(fit.ll))))
  expect_true(all(is.finite(fitted(fit.lp))))
  expect_gte(min(fitted(fit.ll)), 0)
  expect_gte(min(fitted(fit.lp)), 0)
  expect_lt(max(fitted(fit.ll)), 10)
  expect_lt(max(fitted(fit.lp)), 10)
  expect_gt(unname(bw.ll$bandwidth$y), 1e-4)
  expect_gt(unname(bw.lp$bandwidth$y), 1e-4)
  expect_lt(max(fitted(fit.ll)) / max(fitted(npcdens(bws = bw.lc))), 20)
  expect_lt(max(fitted(fit.lp)) / max(fitted(npcdens(bws = bw.lc))), 20)
})
