test_that("npplreg basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 100
  x1 <- runif(n) # nonparametric part
  z1 <- runif(n) # parametric part
  y <- x1^2 + 2*z1 + rnorm(n, sd=0.1)
  
  mydat <- data.frame(y, x1, z1)
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
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

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
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

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

test_that("predict.plregression fail-fast is explicit and plain predict is unchanged", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(13)
  n <- 80
  x1 <- runif(n)
  z1 <- runif(n)
  y <- x1^2 + 2 * z1 + rnorm(n, sd = 0.1)

  bw_mat <- matrix(c(0.1, 0.1), nrow = 2, ncol = 1)
  bw <- npplregbw(xdat = z1, zdat = x1, ydat = y, bws = bw_mat, bandwidth.compute = FALSE)
  model <- npplreg(bws = bw)

  expect_equal(as.numeric(predict(model)), as.numeric(fitted(model)), tolerance = 0)
  expect_error(
    predict(model, se.fit = TRUE),
    "se.fit = TRUE is not implemented"
  )
})

test_that("npplreg residuals stay training-length when evaluation x/z are supplied", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

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
