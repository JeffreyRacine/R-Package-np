test_that("npplreg basic functionality works", {
  set.seed(42)
  n <- 100
  x1 <- runif(n) # nonparametric part
  z1 <- runif(n) # parametric part
  y <- x1^2 + 2*z1 + rnorm(n, sd=0.1)
  
  # Partially linear model: y = g(x1) + z1*beta + e
  # bws needs to be a matrix. Row 1: y on x1, Row 2: z1 on x1
  bw_mat <- matrix(c(0.1, 0.1), nrow=2, ncol=1)
  bw <- npplregbw(xdat=z1, zdat=x1, ydat=y, bws=bw_mat, bandwidth.compute=FALSE)
  
  model <- npplreg(bws=bw)
  
  expect_s3_class(model, "plregression")
  expect_type(predict(model), "double")
  expect_output(summary(model))
  expect_length(coef(model), 1)
})

test_that("npplreg direct-fit objects expose canonical and compatibility bandwidth slots", {
  set.seed(7)
  n <- 80
  x1 <- runif(n)
  z1 <- runif(n)
  y <- x1^2 + 2 * z1 + rnorm(n, sd = 0.1)

  bw_mat <- matrix(c(0.1, 0.1), nrow = 2, ncol = 1)
  bw <- npplregbw(xdat = z1, zdat = x1, ydat = y, bws = bw_mat, bandwidth.compute = FALSE)
  model <- npplreg(bws = bw)

  expect_true("bws" %in% names(model))
  expect_true("bw" %in% names(model))
  expect_s3_class(model$bws, "plbandwidth")
  expect_s3_class(model$bw, "plbandwidth")
  expect_equal(model$bw$bw, model$bws$bw)
  expect_equal(model$bw$fval, model$bws$fval)
})

test_that("plregression methods remain compatible with legacy objects lacking bws", {
  set.seed(11)
  n <- 80
  x1 <- runif(n)
  z1 <- runif(n)
  y <- x1^2 + 2 * z1 + rnorm(n, sd = 0.1)

  bw_mat <- matrix(c(0.1, 0.1), nrow = 2, ncol = 1)
  bw <- npplregbw(xdat = z1, zdat = x1, ydat = y, bws = bw_mat, bandwidth.compute = FALSE)
  model <- npplreg(bws = bw)
  legacy <- model
  legacy$bws <- NULL

  expect_silent(capture.output(print(legacy)))
  expect_silent(capture.output(summary(legacy)))
  expect_type(predict(legacy), "double")
  expect_type(residuals(legacy), "double")
})

test_that("predict.plregression se.fit returns prediction standard errors", {
  set.seed(13)
  n <- 80
  x1 <- runif(n)
  z1 <- runif(n)
  y <- x1^2 + 2 * z1 + rnorm(n, sd = 0.1)

  bw_mat <- matrix(c(0.1, 0.1), nrow = 2, ncol = 1)
  bw <- npplregbw(xdat = z1, zdat = x1, ydat = y, bws = bw_mat, bandwidth.compute = FALSE)
  model <- npplreg(bws = bw)

  expect_equal(as.numeric(predict(model)), as.numeric(fitted(model)), tolerance = 0)
  pred.se <- predict(model, se.fit = TRUE)
  expect_named(pred.se, c("fit", "se.fit"))
  expect_equal(as.numeric(pred.se$fit), as.numeric(fitted(model)), tolerance = 0)
  expect_length(pred.se$se.fit, length(fitted(model)))
  expect_true(all(is.finite(pred.se$se.fit)))
  expect_true(all(pred.se$se.fit >= 0))
})

test_that("predict.plregression se.fit supports evaluation data and mixed x/z", {
  set.seed(20260506)
  n <- 72L
  z1 <- runif(n)
  z2 <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  x1 <- rnorm(n)
  x2 <- factor(sample(c("low", "high"), n, replace = TRUE))
  x2num <- as.numeric(x2 == "high")
  y <- 0.8 * x1 - 0.4 * x2num + cos(2 * pi * z1) +
    0.25 * as.numeric(z2 == "b") + rnorm(n, sd = 0.08)

  tx <- data.frame(x1 = x1, x2 = x2)
  tz <- data.frame(z1 = z1, z2 = z2)
  ex <- data.frame(
    x1 = seq(-1, 1, length.out = 11L),
    x2 = factor(rep(c("low", "high", "low"), length.out = 11L), levels = levels(x2))
  )
  ez <- data.frame(
    z1 = seq(0.05, 0.95, length.out = 11L),
    z2 = factor(rep(c("a", "b", "c"), length.out = 11L), levels = levels(z2))
  )
  bw <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    bws = matrix(c(0.25, 0.25, 0.35, 0.35, 0.3, 0.3), nrow = 3L, ncol = 2L),
    bandwidth.compute = FALSE
  )
  model <- npplreg(bws = bw)
  pred <- predict(model, exdat = ex, ezdat = ez)
  pred.se <- predict(model, exdat = ex, ezdat = ez, se.fit = TRUE)

  expect_equal(as.numeric(pred.se$fit), as.numeric(pred), tolerance = 0)
  expect_length(pred.se$se.fit, length(pred))
  expect_true(all(is.finite(pred.se$se.fit)))
  expect_true(all(pred.se$se.fit >= 0))
})

test_that("predict.plregression se.fit follows formula newdata evaluation rows", {
  set.seed(202605061)
  n <- 65L
  dat <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    z = runif(n)
  )
  dat$y <- 0.7 * dat$x + sin(2 * pi * dat$z) + rnorm(n, sd = 0.12)
  bw <- npplregbw(
    y ~ x | z,
    data = dat,
    bws = matrix(c(0.3, 0.3), nrow = 2L),
    bandwidth.compute = FALSE
  )
  fit <- npplreg(bws = bw)
  nd <- data.frame(
    x = seq(-1, 1, length.out = 9L),
    z = seq(0.1, 0.9, length.out = 9L)
  )

  pred <- predict(fit, newdata = nd)
  pred.se <- predict(fit, newdata = nd, se.fit = TRUE)

  expect_equal(as.numeric(pred.se$fit), as.numeric(pred), tolerance = 0)
  expect_length(pred.se$se.fit, nrow(nd))
  expect_true(all(is.finite(pred.se$se.fit)))
  expect_true(all(pred.se$se.fit >= 0))
})

test_that("npplreg plot-data asymptotic and bootstrap intervals are comparable", {
  set.seed(202605062)
  n <- 80L
  x <- rnorm(n)
  z <- runif(n)
  y <- 0.5 + x + sin(2 * pi * z) + rnorm(n, sd = 0.2)
  bw <- npplregbw(
    xdat = data.frame(x = x),
    zdat = data.frame(z = z),
    ydat = y,
    bws = matrix(c(0.3, 0.3), nrow = 2L),
    bandwidth.compute = FALSE
  )
  fit <- npplreg(bws = bw)
  extract_err <- function(obj)
    unlist(lapply(obj, function(x) rowMeans(abs(x$merr[, 1:2, drop = FALSE]))),
           use.names = FALSE)

  asym <- plot(
    fit,
    output = "data",
    perspective = FALSE,
    errors = "asymptotic",
    band = "pointwise",
    neval = 8L
  )
  boot <- plot(
    fit,
    output = "data",
    perspective = FALSE,
    errors = "bootstrap",
    band = "pointwise",
    B = 39L,
    random.seed = 202605062,
    neval = 8L
  )

  asym.err <- extract_err(asym)
  boot.err <- extract_err(boot)
  ratio <- median(asym.err, na.rm = TRUE) / median(boot.err, na.rm = TRUE)

  expect_true(is.finite(ratio))
  expect_true(ratio > 0.2)
  expect_true(ratio < 10)
})

test_that("npplreg residuals stay training-length when evaluation x/z are supplied", {
  set.seed(20260323)
  n <- 24L
  x1 <- seq(0.1, 1.2, length.out = n)
  z1 <- seq(0.2, 1.4, length.out = n)
  y <- x1^2 + 2 * z1

  bw <- npplregbw(
    xdat = z1,
    zdat = x1,
    ydat = y,
    bws = matrix(c(0.35, 0.35), nrow = 2),
    bandwidth.compute = FALSE
  )
  ex <- data.frame(x1 = z1[c(2, 5, 8)])
  ez <- data.frame(z1 = x1[c(3, 9, 15)])
  ey <- y[c(3, 9, 15)]

  train.fit <- npplreg(bws = bw, residuals = TRUE)
  eval.fit.noey <- npplreg(bws = bw, residuals = TRUE, exdat = ex, ezdat = ez)
  eval.fit.ey <- npplreg(bws = bw, residuals = TRUE, exdat = ex, ezdat = ez, eydat = ey)

  expect_equal(length(eval.fit.noey$mean), nrow(ex), tolerance = 0)
  expect_equal(length(eval.fit.ey$mean), nrow(ex), tolerance = 0)
  expect_equal(length(eval.fit.noey$resid), length(y), tolerance = 0)
  expect_equal(length(eval.fit.ey$resid), length(y), tolerance = 0)
  expect_equal(as.numeric(eval.fit.noey$resid), as.numeric(train.fit$resid), tolerance = 0)
  expect_equal(as.numeric(eval.fit.ey$resid), as.numeric(train.fit$resid), tolerance = 0)
})
