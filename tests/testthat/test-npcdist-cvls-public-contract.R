library(npRmpi)

test_that("public npcdistbw cv.ls keeps lc adjacency live during containment", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260310)
  n <- 30L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + rnorm(n, sd = 0.15))

  bw.lc <- npcdistbw(xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ls", nmulti = 0)
  fit.lc <- npcdist(bws = bw.lc)

  expect_true(is.finite(bw.lc$fval))
  expect_true(all(is.finite(fitted(fit.lc))))
})

test_that("public npcdistbw cv.ls LP/LL routes fail fast during containment", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260311)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  expect_error(
    npcdistbw(xdat = x, ydat = y, regtype = "ll", bwmethod = "cv.ls", nmulti = 0),
    "temporarily disabled pending low-memory shadow CV remediation"
  )
  expect_error(
    npcdistbw(
      xdat = x,
      ydat = y,
      regtype = "lp",
      basis = "glp",
      degree = degree,
      bwmethod = "cv.ls",
      nmulti = 0
    ),
    "temporarily disabled pending low-memory shadow CV remediation"
  )
  expect_error(
    npcdistbw(
      xdat = x,
      ydat = y,
      regtype = "lp",
      basis = "glp",
      degree = degree,
      bwtype = "generalized_nn",
      bwmethod = "cv.ls",
      nmulti = 0
    ),
    "temporarily disabled pending low-memory shadow CV remediation"
  )
})
