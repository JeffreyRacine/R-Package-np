library(np)

test_that("public npcdensbw cv.ml lc matches canonical one-shot evaluation", {
  set.seed(202)
  n <- 32L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = rnorm(n))

  bw.lc <- npcdensbw(xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ml", nmulti = 1)
  objective <- np:::.npcdensbw_eval_only(x, y, bw.lc)$objective

  expect_equal(bw.lc$fval, objective, tolerance = 1e-12)
})

test_that("public npcdensbw cv.ml fixed LP/LL route activates with ll == lp parity", {
  set.seed(101)
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
 
  expect_equal(-bw.ll$fval, -bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdensbw cv.ml fixed LP route preserves tree parity", {
  set.seed(102)
  n <- 34L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + x$x2 + rnorm(n, sd = 0.15))
  degree <- rep.int(1L, ncol(x))

  bw.serial <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ml",
    nmulti = 1
  )

  old_opt <- getOption("np.tree")
  old_compress <- getOption("np.categorical.compress")
  on.exit({
    options(np.tree = old_opt)
    options(np.categorical.compress = old_compress)
  }, add = TRUE)
  options(np.tree = FALSE, np.categorical.compress = TRUE)

  bw.tree <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ml",
    nmulti = 1
  )

  expect_equal(-bw.tree$fval, -bw.serial$fval, tolerance = 1e-8)
})

test_that("public npcdensbw cv.ml generalized-nn LP route activates with ll == lp parity", {
  set.seed(103)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 - x$x2 + rnorm(n, sd = 0.2))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 1
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 1
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdensbw cv.ml adaptive-nn LP route activates with ll == lp parity", {
  set.seed(104)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 - x$x2 + rnorm(n, sd = 0.2))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ml",
    nmulti = 1
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ml",
    nmulti = 1
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})
