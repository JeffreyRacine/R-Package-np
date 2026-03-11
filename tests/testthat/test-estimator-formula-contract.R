test_that("npindex formula newdata path matches explicit data path", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = TRUE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260223)
  dat <- data.frame(
    y = rnorm(36),
    x1 = runif(36),
    x2 = runif(36)
  )
  ex <- data.frame(
    x1 = seq(0.1, 0.9, length.out = 10),
    x2 = seq(0.9, 0.1, length.out = 10)
  )

  bw <- npRmpi::npindexbw(
    y ~ x1 + x2,
    data = dat,
    bws = c(1, 0.35, 0.45),
    bandwidth.compute = FALSE,
    method = "ichimura"
  )

  fit_formula <- npRmpi::npindex(bws = bw, data = dat, newdata = ex)
  fit_default <- npRmpi::npindex(
    bws = bw,
    txdat = dat[c("x1", "x2")],
    tydat = dat$y,
    exdat = ex[c("x1", "x2")]
  )

  expect_equal(fit_formula$mean, fit_default$mean)
  expect_equal(fit_formula$merr, fit_default$merr)
})

test_that("npscoef formula newdata path matches explicit data path", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = TRUE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260223)
  dat <- data.frame(
    y = rnorm(34),
    x1 = runif(34),
    z1 = runif(34)
  )
  ex <- data.frame(
    x1 = seq(0.05, 0.95, length.out = 11),
    z1 = seq(0.9, 0.2, length.out = 11)
  )

  bw <- npRmpi::npscoefbw(
    y ~ x1 | z1,
    data = dat,
    bws = 0.4,
    bandwidth.compute = FALSE
  )

  fit_formula <- npRmpi::npscoef(bws = bw, data = dat, newdata = ex)
  fit_default <- npRmpi::npscoef(
    bws = bw,
    txdat = dat["x1"],
    tydat = dat$y,
    tzdat = dat["z1"],
    exdat = ex["x1"],
    ezdat = ex["z1"]
  )

  expect_equal(fit_formula$mean, fit_default$mean)
  expect_equal(fit_formula$merr, fit_default$merr)
})

test_that("npplreg formula newdata path matches explicit data path", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = TRUE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260223)
  dat <- data.frame(
    y = rnorm(30),
    x1 = runif(30),
    z1 = runif(30)
  )
  ex <- data.frame(
    x1 = seq(0.1, 0.8, length.out = 9),
    z1 = seq(0.8, 0.1, length.out = 9)
  )

  bw <- npRmpi::npplregbw(
    y ~ x1 | z1,
    data = dat,
    bws = matrix(c(0.45, 0.35), nrow = 2, ncol = 1),
    bandwidth.compute = FALSE
  )

  fit_formula <- npRmpi::npplreg(bws = bw, data = dat, newdata = ex)
  fit_default <- npRmpi::npplreg(
    bws = bw,
    txdat = dat["x1"],
    tydat = dat$y,
    tzdat = dat["z1"],
    exdat = ex["x1"],
    ezdat = ex["z1"]
  )

  expect_equal(fit_formula$mean, fit_default$mean)
  expect_equal(fit_formula$merr, fit_default$merr)
})
