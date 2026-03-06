library(npRmpi)

test_that("public npcdistbw cv.ls enforces ll == lp(glp, degree=1)", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260310)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdistbw(xdat = x, ydat = y, regtype = "ll", bwmethod = "cv.ls", nmulti = 0)
  bw.lp <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 0
  )

  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-10)
  expect_equal(bw.ll$xbw, bw.lp$xbw, tolerance = 1e-10)
  expect_equal(bw.ll$ybw, bw.lp$ybw, tolerance = 1e-10)
})

test_that("public npcdistbw cv.ls preserves tree parity on the LP route", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old <- options(np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20260311)
  n <- 34L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  options(np.tree = FALSE)
  bw.nt <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 0
  )

  options(np.tree = TRUE)
  bw.tr <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 0
  )

  expect_equal(bw.nt$fval, bw.tr$fval, tolerance = 1e-10)
})
