test_that("npplregbw ll stores canonical engine metadata and child regression bandwidths", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260306)
  n <- 40
  x <- data.frame(x = runif(n))
  z <- data.frame(z = runif(n))
  y <- sin(2 * pi * z$z) + x$x + rnorm(n, sd = 0.05)

  bw <- npplregbw(
    xdat = x,
    ydat = y,
    zdat = z,
    bws = matrix(0.25, nrow = 2L, ncol = 1L),
    regtype = "ll",
    bandwidth.compute = FALSE
  )

  expect_identical(bw$regtype, "ll")
  expect_identical(bw$basis, "glp")
  expect_identical(as.integer(bw$degree), 1L)
  expect_false(isTRUE(bw$bernstein.basis))
  expect_identical(bw$regtype.engine, "lp")
  expect_identical(bw$basis.engine, "glp")
  expect_identical(as.integer(bw$degree.engine), 1L)
  expect_false(isTRUE(bw$bernstein.basis.engine))

  expect_identical(bw$bw$yzbw$regtype, "lp")
  expect_identical(bw$bw$yzbw$basis, "glp")
  expect_identical(as.integer(bw$bw$yzbw$degree), 1L)
  expect_false(isTRUE(bw$bw$yzbw$bernstein.basis))

  expect_identical(bw$bw[[2]]$regtype, "lp")
  expect_identical(bw$bw[[2]]$basis, "glp")
  expect_identical(as.integer(bw$bw[[2]]$degree), 1L)
  expect_false(isTRUE(bw$bw[[2]]$bernstein.basis))
})

test_that("npplregbw ll rejects non-canonical LP controls", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260306)
  n <- 30
  x <- data.frame(x = runif(n))
  z <- data.frame(z = runif(n))
  y <- x$x + z$z + rnorm(n, sd = 0.05)

  expect_error(
    npplregbw(
      xdat = x,
      ydat = y,
      zdat = z,
      regtype = "ll",
      degree = 2L,
      nmulti = 1L
    ),
    "canonical LP\\(degree=1, basis='glp'\\)"
  )

  expect_error(
    npplregbw(
      xdat = x,
      ydat = y,
      zdat = z,
      regtype = "ll",
      basis = "tensor",
      nmulti = 1L
    ),
    "canonical basis='glp'"
  )

  expect_error(
    npplregbw(
      xdat = x,
      ydat = y,
      zdat = z,
      regtype = "ll",
      bernstein.basis = TRUE,
      nmulti = 1L
    ),
    "canonical bernstein.basis=FALSE"
  )
})

test_that("npindexbw ll stores canonical engine metadata", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260306)
  n <- 50
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(2 * pi * x$x1) + 0.5 * x$x2 + rnorm(n, sd = 0.05)

  bw <- npindexbw(
    xdat = x,
    ydat = y,
    bws = c(1, 0, 0.25),
    regtype = "ll",
    bandwidth.compute = FALSE
  )

  expect_identical(bw$regtype, "ll")
  expect_identical(bw$basis, "glp")
  expect_identical(as.integer(bw$degree), 1L)
  expect_false(isTRUE(bw$bernstein.basis))
  expect_identical(bw$regtype.engine, "lp")
  expect_identical(bw$basis.engine, "glp")
  expect_identical(as.integer(bw$degree.engine), 1L)
  expect_false(isTRUE(bw$bernstein.basis.engine))
})

test_that("npindexbw ll rejects non-canonical LP controls", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260306)
  n <- 40
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- x$x1 + x$x2 + rnorm(n, sd = 0.05)

  expect_error(
    npindexbw(
      xdat = x,
      ydat = y,
      regtype = "ll",
      degree = 2L,
      method = "ichimura",
      nmulti = 1L
    ),
    "canonical LP\\(degree=1, basis='glp'\\)"
  )

  expect_error(
    npindexbw(
      xdat = x,
      ydat = y,
      regtype = "ll",
      basis = "tensor",
      method = "ichimura",
      nmulti = 1L
    ),
    "canonical basis='glp'"
  )

  expect_error(
    npindexbw(
      xdat = x,
      ydat = y,
      regtype = "ll",
      bernstein.basis = TRUE,
      method = "ichimura",
      nmulti = 1L
    ),
    "canonical bernstein.basis=FALSE"
  )
})
