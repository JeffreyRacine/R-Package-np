library(npRmpi)

phase1_npcdens_cvml_gnn_fixture <- function() {
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

test_that("phase1 npcdensbw cv.ml generalized-nn lc matches the frozen public baseline", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  dat <- phase1_npcdens_cvml_gnn_fixture()

  bw.lc <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    regtype = "lc",
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 1L,
    itmax = 1L
  )

  expect_true(is.finite(bw.lc$fval))
  expect_equal(bw.lc$fval, 23.289998754972199, tolerance = 1e-10)
})

test_that("phase1 npcdensbw cv.ml generalized-nn keeps ll on canonical lp degree-1 glp", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  dat <- phase1_npcdens_cvml_gnn_fixture()
  degree <- rep.int(1L, ncol(dat$x))

  bw.ll <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 1L,
    itmax = 1L
  )
  bw.lp <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 1L,
    itmax = 1L
  )

  expect_identical(bw.ll$regtype.engine, "lp")
  expect_identical(bw.ll$basis.engine, "glp")
  expect_identical(as.integer(bw.ll$degree.engine), degree)
  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, 19.7753475006634, tolerance = 1e-10)
  expect_equal(bw.lp$fval, 19.7753475006634, tolerance = 1e-10)
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-10)
})

test_that("phase1 npcdensbw cv.ml generalized-nn lp degree-2 succeeds on a higher-order fixture", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  dat <- phase1_npcdens_cvml_gnn_fixture()
  degree1 <- rep.int(1L, ncol(dat$x))
  degree2 <- rep.int(2L, ncol(dat$x))

  bw.d1 <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    regtype = "lp",
    basis = "glp",
    degree = degree1,
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 1L,
    itmax = 1L
  )
  bw.d2 <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    regtype = "lp",
    basis = "glp",
    degree = degree2,
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 1L,
    itmax = 1L
  )

  expect_identical(as.integer(bw.d2$degree.engine), degree2)
  expect_true(is.finite(bw.d2$fval))
  expect_equal(bw.d2$fval, 25.9256868232211, tolerance = 1e-10)
  expect_gt(abs(bw.d2$fval - bw.d1$fval), 1e-6)
})

test_that("phase1 npcdensbw cv.ml generalized-nn avoids search-boundary collapse", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260309)
  n <- 60L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = x$x1^2 + rnorm(n, sd = 0.06))

  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 1L
  )
  fit <- npcdens(bws = bw)

  expect_true(is.finite(bw$fval))
  expect_true(all(is.finite(fitted(fit))))
  expect_true(all(bw$xbw > 1))
  expect_true(all(bw$ybw > 1))
  expect_true(all(bw$xbw < n))
  expect_true(all(bw$ybw < n))
})
