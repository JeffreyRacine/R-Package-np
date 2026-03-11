test_that("npreg formula newdata path matches explicit data path", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(32),
    x = runif(32)
  )
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 9))

  bw <- npRmpi::npregbw(
    y ~ x,
    data = dat,
    bws = 0.45,
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  fit_formula <- npRmpi::npreg(bws = bw, data = dat, newdata = ex)
  fit_default <- npRmpi::npreg(
    bws = bw,
    txdat = dat["x"],
    tydat = dat$y,
    exdat = ex["x"]
  )

  expect_equal(fit_formula$mean, fit_default$mean)
  expect_equal(fit_formula$merr, fit_default$merr)
})

test_that("npreg formula y.eval path matches explicit data path", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(36),
    x = runif(36)
  )
  ex <- data.frame(
    y = rnorm(12),
    x = seq(0.1, 0.9, length.out = 12)
  )

  bw <- npRmpi::npregbw(
    y ~ x,
    data = dat,
    bws = 0.4,
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  fit_formula <- npRmpi::npreg(bws = bw, data = dat, newdata = ex, y.eval = TRUE)
  fit_default <- npRmpi::npreg(
    bws = bw,
    txdat = dat["x"],
    tydat = dat$y,
    exdat = ex["x"],
    eydat = ex$y
  )

  expect_equal(fit_formula$mean, fit_default$mean)
  expect_equal(fit_formula$merr, fit_default$merr)
})
