test_that("npudens core smoke stays alive", {
  set.seed(1)
  x <- data.frame(x = seq(0.1, 0.9, length.out = 24L))

  bw <- npudensbw(dat = x, bws = 0.2, bandwidth.compute = FALSE)
  fit <- npudens(bws = bw)

  expect_s3_class(fit, "npdensity")
  expect_equal(length(predict(fit)), nrow(x))
  expect_true(all(is.finite(predict(fit))))
})

test_that("npudist core smoke stays alive", {
  set.seed(2)
  x <- data.frame(x = seq(0.1, 0.9, length.out = 24L))

  bw <- npudistbw(dat = x, bws = 0.2, bandwidth.compute = FALSE)
  fit <- npudist(bws = bw)

  expect_s3_class(fit, "npdistribution")
  expect_equal(length(predict(fit)), nrow(x))
  expect_true(all(predict(fit) >= 0 & predict(fit) <= 1))
})

test_that("npreg core smoke aligns lagged time series before fitting", {
  set.seed(20260915)
  y <- stats::ts(rnorm(28), start = c(2001, 1), frequency = 4)
  form <- y ~ stats::lag(y, -1) + stats::lag(y, -2)
  set.seed(42)
  bw <- npregbw(form, nmulti = 1L, regtype = "ll")
  set.seed(42)
  direct <- npreg(form, nmulti = 1L, regtype = "ll")
  two_step <- npreg(bws = bw)
  manual <- npreg(bws = bw,
    txdat = data.frame(lag1 = as.numeric(y)[2:27], lag2 = as.numeric(y)[1:26]),
    tydat = as.numeric(y)[3:28])

  expect_s3_class(direct, "npregression")
  expect_equal(direct$nobs, 26L)
  expect_true(all(is.finite(predict(direct))))
  expect_equal(direct$bws$bw, bw$bw, tolerance = 0)
  expect_equal(fitted(direct), fitted(two_step), tolerance = 0)
  expect_equal(fitted(direct), fitted(manual), tolerance = 0)
})

test_that("npcdens core smoke stays alive", {
  set.seed(4)
  x <- data.frame(x = seq(0.1, 1.0, length.out = 24L))
  y <- data.frame(y = x$x^2)

  bw <- npcdensbw(xdat = x, ydat = y, bws = c(0.25, 0.25), bandwidth.compute = FALSE)
  fit <- npcdens(bws = bw)

  expect_s3_class(fit, "condensity")
  expect_equal(length(predict(fit)), nrow(x))
  expect_true(all(is.finite(predict(fit))))
})

test_that("npcdist core smoke stays alive", {
  set.seed(5)
  x <- data.frame(x = seq(0.1, 1.0, length.out = 24L))
  y <- data.frame(y = x$x^2)

  bw <- npcdistbw(xdat = x, ydat = y, bws = c(0.25, 0.25), bandwidth.compute = FALSE)
  fit <- npcdist(bws = bw)

  expect_s3_class(fit, "condistribution")
  expect_equal(length(predict(fit)), nrow(x))
  expect_true(all(predict(fit) >= 0 & predict(fit) <= 1))
})

test_that("npplreg core smoke stays alive", {
  set.seed(6)
  n <- 24L
  z <- data.frame(z = seq(0.1, 1.0, length.out = n))
  x <- data.frame(x = seq(0.2, 1.1, length.out = n))
  y <- z$z^2 + 2 * x$x

  bw <- npplregbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = matrix(c(0.25, 0.25), nrow = 2L),
    bandwidth.compute = FALSE
  )
  fit <- npplreg(bws = bw)

  expect_s3_class(fit, "plregression")
  expect_equal(length(predict(fit)), n)
  expect_true(all(is.finite(predict(fit))))
})

test_that("npindex raw formulas align lagged time series", {
  set.seed(42)
  y <- ts(rnorm(28), frequency = 4)
  f <- y ~ lag(y, -1) + lag(y, -2)
  h <- c(1, .5, .8)
  direct <- npindex(f, bws = h, bandwidth.compute = FALSE)
  bw <- npindexbw(f, bws = h, bandwidth.compute = FALSE)
  control <- npindex(bws = bw)
  expect_identical(direct$nobs, 26L)
  expect_equal(fitted(direct), fitted(control), tolerance = 1e-14)
  expect_equal(vcov(direct), vcov(control), tolerance = 1e-14)
})

test_that("npindex core smoke stays alive", {
  set.seed(7)
  n <- 24L
  x <- data.frame(
    x1 = seq(0.1, 1.0, length.out = n),
    x2 = seq(1.0, 0.1, length.out = n)
  )
  y <- x$x1 - x$x2

  bw <- npindexbw(
    xdat = x,
    ydat = y,
    method = "ichimura",
    bws = c(1, 0.25, 0.25),
    bandwidth.compute = FALSE
  )
  fit <- npindex(bws = bw)

  expect_s3_class(fit, "singleindex")
  expect_equal(length(predict(fit)), n)
  expect_true(all(is.finite(predict(fit))))
})

test_that("npscoef core smoke stays alive", {
  set.seed(8)
  n <- 24L
  x <- data.frame(x = seq(0.1, 1.0, length.out = n))
  z <- data.frame(z = seq(1.0, 0.1, length.out = n))
  y <- x$x * z$z

  bw <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  fit <- npscoef(bws = bw, iterate = FALSE, se = FALSE)

  expect_s3_class(fit, "smoothcoefficient")
  expect_equal(length(predict(fit)), n)
  expect_true(all(is.finite(predict(fit))))
})
test_that("copula grid coordinates cannot fall back to equal-length training data", {
  ns <- asNamespace(getNamespaceName(environment(npcopula)))
  coordinates <- get(".npcopula_eval_xgrid", ns)
  x <- structure(list(evaluation = "grid", xnames = c("x", "log(y)"),
    copula = c(.2, .4, .6, .8), grid.dim = c(2L, 2L),
    data = data.frame(x = 1:4, check.names = FALSE, "log(y)" = 5:8),
    eval = data.frame(copula = c(.2, .4, .6, .8), u1 = c(.2, .8, .2, .8),
      u2 = c(.2, .2, .8, .8), x = 9:12, log.y. = 13:16)),
    class = "npcopula")
  expect_identical(coordinates(x)[[1L]], 9:12)
  expect_identical(names(coordinates(x)), x$xnames)
  x$evaluation <- "sample"
  expect_identical(coordinates(x)[[1L]], 1:4)
})
