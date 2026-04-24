library(npRmpi)

test_that("adaptive-nn conditional density lp matches ll and stays off the search boundary", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(42)
  n <- 120L
  x <- runif(n)
  y <- x^2 + rnorm(n, sd = 0.1)

  bw.ll <- npcdensbw(
    y ~ x,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ml",
    nmulti = 1,
    itmax = 1
  )
  bw.lp <- npcdensbw(
    y ~ x,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ml",
    nmulti = 1,
    itmax = 1
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_gt(bw.lp$xbw, 0)
  expect_gt(bw.lp$ybw, 0)
  expect_lt(bw.lp$xbw, 50)
  expect_lt(bw.lp$ybw, 50)
  expect_equal(c(bw.ll$xbw, bw.ll$ybw), c(bw.lp$xbw, bw.lp$ybw), tolerance = 1e-8)
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)

  fit.ll <- npcdens(bws = bw.ll, txdat = data.frame(x = x), tydat = data.frame(y = y))
  fit.lp <- npcdens(bws = bw.lp, txdat = data.frame(x = x), tydat = data.frame(y = y))

  expect_true(all(is.finite(fitted(fit.ll))))
  expect_true(all(is.finite(fitted(fit.lp))))
  expect_equal(fitted(fit.ll), fitted(fit.lp), tolerance = 1e-8)
})

test_that("adaptive-nn conditional distribution lp matches ll and stays off the search boundary", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(42)
  n <- 120L
  x <- runif(n)
  y <- x^2 + rnorm(n, sd = 0.1)

  bw.ll <- npcdistbw(
    y ~ x,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1
  )
  bw.lp <- npcdistbw(
    y ~ x,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_gt(bw.lp$xbw, 0)
  expect_gt(bw.lp$ybw, 0)
  expect_lt(bw.lp$xbw, 50)
  expect_lt(bw.lp$ybw, 50)
  expect_equal(c(bw.ll$xbw, bw.ll$ybw), c(bw.lp$xbw, bw.lp$ybw), tolerance = 1e-8)
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)

  fit.ll <- npcdist(bws = bw.ll, txdat = data.frame(x = x), tydat = data.frame(y = y))
  fit.lp <- npcdist(bws = bw.lp, txdat = data.frame(x = x), tydat = data.frame(y = y))

  expect_true(all(is.finite(fitted(fit.ll))))
  expect_true(all(is.finite(fitted(fit.lp))))
  expect_equal(fitted(fit.ll), fitted(fit.lp), tolerance = 1e-8)
})
