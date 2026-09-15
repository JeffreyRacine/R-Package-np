test_that("one-call regression aligns lagged time series before fitting", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = TRUE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260915)
  y <- stats::ts(rnorm(28), start = c(2001, 1), frequency = 4)
  form <- y ~ stats::lag(y, -1) + stats::lag(y, -2)
  set.seed(42)
  bw <- npRmpi::npregbw(form, nmulti = 1L, regtype = "ll")
  set.seed(42)
  direct <- npRmpi::npreg(form, nmulti = 1L, regtype = "ll")
  two_step <- npRmpi::npreg(bws = bw)
  manual <- npRmpi::npreg(bws = bw,
    txdat = data.frame(lag1 = as.numeric(y)[2:27], lag2 = as.numeric(y)[1:26]),
    tydat = as.numeric(y)[3:28])

  expect_equal(direct$nobs, 26L)
  expect_equal(direct$bws$bw, bw$bw, tolerance = 0)
  expect_equal(fitted(direct), fitted(two_step), tolerance = 0)
  expect_equal(fitted(direct), fitted(manual), tolerance = 0)
})

test_that("npreg formula newdata path matches explicit data path", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = TRUE)
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
  old.opts <- options(npRmpi.autodispatch = TRUE)
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

test_that("npreg formula fixed-bandwidth predict works with autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = TRUE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260512)
  dat <- data.frame(
    y = rnorm(48),
    x = runif(48)
  )
  ex <- data.frame(x = seq(0.05, 0.95, length.out = 7))

  bw <- npRmpi::npregbw(
    y ~ x,
    data = dat,
    regtype = "ll",
    bws = 0.35,
    bandwidth.compute = FALSE
  )

  fit_formula <- npRmpi::npreg(bws = bw)
  pred_formula <- predict(fit_formula, newdata = ex)

  fit_explicit <- npRmpi::npreg(
    bws = bw,
    txdat = dat["x"],
    tydat = dat$y,
    exdat = ex["x"]
  )

  expect_equal(pred_formula, fitted(fit_explicit))
})

test_that("npreg direct formula scalar bandwidth works with autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = TRUE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260512)
  dat <- data.frame(
    y = rnorm(48),
    x = runif(48)
  )

  fit_formula <- npRmpi::npreg(
    y ~ x,
    data = dat,
    regtype = "ll",
    bws = 0.35
  )
  bw <- npRmpi::npregbw(
    y ~ x,
    data = dat,
    regtype = "ll",
    bws = 0.35,
    bandwidth.compute = FALSE
  )
  fit_explicit <- npRmpi::npreg(
    bws = bw,
    txdat = dat["x"],
    tydat = dat$y
  )

  expect_equal(fitted(fit_formula), fitted(fit_explicit))
})
