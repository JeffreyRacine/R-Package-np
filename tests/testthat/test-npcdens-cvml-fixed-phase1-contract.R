library(npRmpi)

phase1_npcdens_cvml_fixed_fixture <- function() {
  set.seed(20260309)
  n <- 36L
  x <- data.frame(
    x1 = runif(n),
    x2 = runif(n)
  )
  y <- data.frame(
    y1 = x$x1^2 - 0.35 * x$x2 + 0.2 * sin(2 * pi * x$x1) + rnorm(n, sd = 0.08)
  )
  list(x = x, y = y)
}

test_that("phase1 npcdensbw cv.ml fixed lc matches the frozen public baseline", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  dat <- phase1_npcdens_cvml_fixed_fixture()

  bw.lc <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    regtype = "lc",
    bwtype = "fixed",
    bwmethod = "cv.ml",
    nmulti = 1L
  )

  expect_true(is.finite(bw.lc$fval))
  expect_equal(bw.lc$fval, 23.4336801011909, tolerance = 1e-10)
})

test_that("phase1 npcdensbw cv.ml fixed keeps ll on canonical lp degree-1 glp", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  dat <- phase1_npcdens_cvml_fixed_fixture()
  degree <- rep.int(1L, ncol(dat$x))

  bw.ll <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    regtype = "ll",
    bwtype = "fixed",
    bwmethod = "cv.ml",
    nmulti = 1L
  )
  bw.lp <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "fixed",
    bwmethod = "cv.ml",
    nmulti = 1L
  )

  expect_identical(bw.ll$regtype.engine, "lp")
  expect_identical(bw.ll$basis.engine, "glp")
  expect_identical(as.integer(bw.ll$degree.engine), degree)
  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, 15.7534838033289, tolerance = 1e-10)
  expect_equal(bw.lp$fval, 15.7534838033289, tolerance = 1e-10)
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-10)
})

test_that("phase1 npcdensbw cv.ml fixed lp degree-2 succeeds on a higher-order fixture", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  dat <- phase1_npcdens_cvml_fixed_fixture()
  degree1 <- rep.int(1L, ncol(dat$x))
  degree2 <- rep.int(2L, ncol(dat$x))

  bw.d1 <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    regtype = "lp",
    basis = "glp",
    degree = degree1,
    bwtype = "fixed",
    bwmethod = "cv.ml",
    nmulti = 1L
  )
  bw.d2 <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    regtype = "lp",
    basis = "glp",
    degree = degree2,
    bwtype = "fixed",
    bwmethod = "cv.ml",
    nmulti = 1L
  )

  expect_identical(as.integer(bw.d2$degree.engine), degree2)
  expect_true(is.finite(bw.d2$fval))
  expect_gt(abs(bw.d2$fval - bw.d1$fval), 1e-6)
})

test_that("phase1 npcdensbw cv.ml fixed avoids boundary collapse on bounded support", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260309)
  n <- 60L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = rbeta(n, shape1 = 1 + 2 * x$x1, shape2 = 2))

  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ml",
    nmulti = 1L,
    cykerbound = "fixed",
    cykerlb = 0,
    cykerub = 1
  )
  fit <- npcdens(bws = bw)
  grid <- data.frame(y1 = c(0, 0.02, 0.98, 1), x1 = rep(0.5, 4L))
  pred <- predict(fit, newdata = grid)

  expect_true(is.finite(bw$fval))
  expect_true(all(is.finite(pred)))
  expect_true(all(pred > 0))
})
