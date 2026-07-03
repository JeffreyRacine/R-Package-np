library(npRmpi)

test_that("public npcdistbw cv.ls keeps lc adjacency live during containment", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260310)
  n <- 30L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + rnorm(n, sd = 0.15))

  bw.lc <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lc",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )
  fit.lc <- npcdist(bws = bw.lc)

  expect_true(is.finite(bw.lc$fval))
  expect_true(all(is.finite(fitted(fit.lc))))
})

test_that("public npcdistbw cv.ls fixed LP/LL route activates with ll == lp parity", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260311)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )
  bw.lp <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("npcdistbw cv.ls fixed continuous stream follows the tree option", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old_dispatch <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_dispatch), add = TRUE)

  data("wage1")
  x <- wage1[, c("married", "female", "nonwhite", "educ", "exper", "tenure")]
  y <- wage1["lwage"]

  bw <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    degree = rep(1L, 3L),
    bwmethod = "cv.ls",
    cxkertype = "epanechnikov",
    cykertype = "epanechnikov",
    bws = c(1, 0.1, 0.1, 0.1, 5, 10, 5),
    bandwidth.compute = FALSE
  )

  old_opts <- options(np.tree = TRUE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  expect_equal(
    npRmpi:::.npcdistbw_tree_code(
      bw,
      ncon = bw$yncon + bw$xncon,
      ncat = bw$ynuno + bw$ynord + bw$xnuno + bw$xnord
    ),
    npRmpi:::DO_TREE_YES
  )

  options(np.tree = FALSE)
  expect_equal(
    npRmpi:::.npcdistbw_tree_code(
      bw,
      ncon = bw$yncon + bw$xncon,
      ncat = bw$ynuno + bw$ynord + bw$xnuno + bw$xnord
    ),
    npRmpi:::DO_TREE_NO
  )
})

test_that("public npcdistbw cv.ls generalized-nn LP route activates with ll == lp parity", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260312)
  n <- 24L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )
  bw.lp <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "generalized_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdistbw cv.ls adaptive-nn LP route activates with ll == lp parity", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260313)
  n <- 24L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )
  bw.lp <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})
