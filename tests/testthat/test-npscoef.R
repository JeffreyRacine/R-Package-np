test_that("npscoef basic functionality works", {
  set.seed(42)
  n <- 100
  x1 <- runif(n) # smoothing variable
  z1 <- runif(n) # parametric variable
  # Smooth coefficient model: y = beta(x1) * z1 + e
  y <- (x1^2) * z1 + rnorm(n, sd=0.1)
  
  bw <- npscoefbw(xdat=x1, zdat=z1, ydat=y, bws=0.1, bandwidth.compute=FALSE)
  
  model <- npscoef(bws=bw)
  
  expect_s3_class(model, "smoothcoefficient")
  expect_type(predict(model), "double")
  expect_output(summary(model))
})

test_that("npscoefbw records ll/lp controls", {
  set.seed(43)
  n <- 80
  x1 <- runif(n)
  z1 <- runif(n)
  y <- (0.2 + x1) * cos(2 * pi * z1) + rnorm(n, sd = 0.08)

  bw.ll <- npscoefbw(
    xdat = x1,
    zdat = z1,
    ydat = y,
    regtype = "ll",
    bws = 0.15,
    bandwidth.compute = FALSE
  )
  expect_identical(bw.ll$regtype, "ll")

  bw.lp <- npscoefbw(
    xdat = x1,
    zdat = z1,
    ydat = y,
    regtype = "lp",
    basis = "tensor",
    degree = 2L,
    bws = 0.15,
    bandwidth.compute = FALSE
  )
  expect_identical(bw.lp$regtype, "lp")
  expect_identical(bw.lp$basis, "tensor")
})

test_that("npscoef direct route accepts omitted tzdat with explicit bandwidths", {
  x <- data.frame(x = c(0.1, 0.3, 0.5, 0.8, 0.9, 1.1))
  y <- c(0.12, 0.25, 0.44, 0.61, 0.73, 0.95)

  bw <- npscoefbw(xdat = x, ydat = y, bws = 0.35, bandwidth.compute = FALSE)
  fit <- npscoef(bws = bw, txdat = x, tydat = y, errors = FALSE, iterate = FALSE)

  expect_s3_class(fit, "smoothcoefficient")
  expect_equal(length(fit$mean), nrow(x))
})

test_that("npscoef explicit tzdat direct route matches stored-data route", {
  x <- data.frame(x = c(0.1, 0.3, 0.5, 0.8, 0.9, 1.1))
  z <- data.frame(z = c(1.0, 0.5, 1.5, 0.8, 1.2, 0.7))
  y <- c(0.12, 0.25, 0.44, 0.61, 0.73, 0.95)

  bw <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = 0.35,
    bandwidth.compute = FALSE
  )

  fit.stored <- npscoef(bws = bw, errors = FALSE, iterate = FALSE)
  fit.direct <- npscoef(
    bws = bw,
    txdat = x,
    tydat = y,
    tzdat = z,
    errors = FALSE,
    iterate = FALSE
  )

  expect_equal(fit.direct$mean, fit.stored$mean, tolerance = 0)
  expect_equal(fit.direct$bws$bw, fit.stored$bws$bw, tolerance = 0)
})
