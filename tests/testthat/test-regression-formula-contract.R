test_that("npreg formula newdata path matches explicit data path", {
  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(32),
    x = runif(32)
  )
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 9))

  bw <- np::npregbw(
    y ~ x,
    data = dat,
    bws = 0.45,
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  fit_formula <- np::npreg(bws = bw, data = dat, newdata = ex)
  fit_default <- np::npreg(
    bws = bw,
    txdat = dat["x"],
    tydat = dat$y,
    exdat = ex["x"]
  )

  expect_equal(fit_formula$mean, fit_default$mean)
  expect_equal(fit_formula$merr, fit_default$merr)
})

test_that("npreg formula y.eval path matches explicit data path", {
  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(36),
    x = runif(36)
  )
  ex <- data.frame(
    y = rnorm(12),
    x = seq(0.1, 0.9, length.out = 12)
  )

  bw <- np::npregbw(
    y ~ x,
    data = dat,
    bws = 0.4,
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  fit_formula <- np::npreg(bws = bw, data = dat, newdata = ex, y.eval = TRUE)
  fit_default <- np::npreg(
    bws = bw,
    txdat = dat["x"],
    tydat = dat$y,
    exdat = ex["x"],
    eydat = ex$y
  )

  expect_equal(fit_formula$mean, fit_default$mean)
  expect_equal(fit_formula$merr, fit_default$merr)
})
