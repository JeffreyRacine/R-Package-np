library(npRmpi)

test_that("session-route qreg plots preserve quantreg mode and fitted tau", {
  npRmpi.init(nslaves = 1, quiet = TRUE)
  on.exit(npRmpi.quit(), add = TRUE)

  set.seed(108)
  n <- 24
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.05)
  xdat <- data.frame(x = x)
  ydat <- data.frame(y = y)

  bw <- npcdistbw(xdat = xdat, ydat = ydat, bws = c(0.25, 0.25), bandwidth.compute = FALSE)
  fit <- npqreg(bws = bw, txdat = xdat, tydat = ydat, tau = 0.4)

  bw.fixed <- plot(
    bw,
    xdat = xdat,
    ydat = ydat,
    plot.behavior = "data",
    perspective = TRUE,
    view = "fixed",
    quantreg = TRUE,
    tau = 0.4,
    neval = 8
  )
  fit.fixed <- plot(
    fit,
    plot.behavior = "data",
    perspective = TRUE,
    view = "fixed",
    neval = 8
  )
  bw.slice <- plot(
    bw,
    xdat = xdat,
    ydat = ydat,
    plot.behavior = "data",
    perspective = FALSE,
    quantreg = TRUE,
    tau = 0.4,
    neval = 8
  )
  fit.slice <- plot(
    fit,
    plot.behavior = "data",
    perspective = FALSE,
    neval = 8
  )

  expect_true(all(vapply(fit.fixed, inherits, logical(1), "qregression")))
  expect_true(all(vapply(fit.slice, inherits, logical(1), "qregression")))
  expect_true(all(vapply(fit.fixed, function(xi) identical(xi$tau, 0.4), logical(1))))
  expect_true(all(vapply(fit.slice, function(xi) identical(xi$tau, 0.4), logical(1))))
  expect_equal(fit.fixed[[1L]]$quantile, bw.fixed[[1L]]$quantile, tolerance = 1e-5)
  expect_equal(fit.fixed[[1L]]$xeval, bw.fixed[[1L]]$xeval)
  expect_equal(fit.slice[[1L]]$quantile, bw.slice[[1L]]$quantile, tolerance = 1e-5)
  expect_equal(fit.slice[[1L]]$xeval, bw.slice[[1L]]$xeval)
})
