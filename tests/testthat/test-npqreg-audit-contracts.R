test_that("npqreg selected-CDF inversion honors iteration, midpoint, and clamp contracts", {
  inv <- getFromNamespace(".npqreg_invert_selected_cdf", "npRmpi")
  validate_itmax <- getFromNamespace(".npqreg_validate_itmax", "npRmpi")
  quantile_clamp <- getFromNamespace(".npqreg_quantile_clamp", "npRmpi")
  mark_clamped_delta <- getFromNamespace(".npqreg_mark_clamped_delta", "npRmpi")

  bws <- list(regtype.engine = "lc")
  xdat <- data.frame(x = 0)
  ydat <- data.frame(y = c(-1, 1))
  exdat <- data.frame(x = 0)
  linear_cdf <- function(bws, xdat, ydat, exdat, ycand) {
    (as.double(ycand) + 1) / 2
  }

  expect_identical(validate_itmax(1500), 1500L)
  expect_error(validate_itmax(.Machine$integer.max + 1), "positive integer")
  expect_error(
    inv(bws, xdat, ydat, exdat, tau = 0.6, tol = .Machine$double.eps,
        small = .Machine$double.eps, itmax = 1L, cdf.values = linear_cdf),
    "failed to converge"
  )
  expect_error(
    inv(list(regtype.engine = "lc", cykerorder = 4L),
        xdat, ydat, exdat, tau = 0.6, tol = 0.5,
        small = .Machine$double.eps, itmax = 1500L,
        cdf.values = linear_cdf),
    "cykerorder = 2"
  )

  out <- inv(bws, xdat, ydat, exdat, tau = 0.6, tol = 0.5,
             small = .Machine$double.eps, itmax = 1500L,
             cdf.values = linear_cdf)
  expect_equal(as.numeric(out), 0.25, tolerance = 0)
  expect_identical(quantile_clamp(out), "none")
})

test_that("npqreg MPI internal local inversion carries clamp metadata", {
  inv <- getFromNamespace(".npqreg_invert_selected_cdf", "npRmpi")
  quantile_clamp <- getFromNamespace(".npqreg_quantile_clamp", "npRmpi")
  mark_clamped_delta <- getFromNamespace(".npqreg_mark_clamped_delta", "npRmpi")

  bws <- list(regtype.engine = "lc")
  xdat <- data.frame(x = 0)
  ydat <- data.frame(y = c(-1, 1))
  exdat3 <- data.frame(x = 1:3)

  # npRmpi selects its CDF evaluator internally; use constant support to test
  # the shared clamp contract without requiring an active MPI communicator.
  out <- inv(bws, xdat, data.frame(y = c(0, 0)), exdat3, tau = 0.5,
             tol = 0.5, small = .Machine$double.eps, itmax = 10L)
  expect_equal(as.numeric(out), c(0, 0, 0), tolerance = 0)
  expect_identical(quantile_clamp(out), c("constant", "constant", "constant"))

  delta <- list(
    quanterr = c(1, 2, 3),
    quantgrad = matrix(1, nrow = 3, ncol = 2),
    quantgerr = matrix(2, nrow = 3, ncol = 2)
  )
  masked <- mark_clamped_delta(delta, quantile_clamp(out))
  expect_true(all(is.na(masked$quanterr)))
  expect_true(all(is.na(masked$quantgrad)))
  expect_true(all(is.na(masked$quantgerr)))
})

test_that("npRmpi npqreg direct evaluation preserves NA rows for vector tau gradients", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260703)
  xdat <- data.frame(x = seq(0.05, 0.95, length.out = 40))
  ydat <- data.frame(y = sin(2 * pi * xdat$x) + rnorm(nrow(xdat), sd = 0.03))
  exdat <- data.frame(x = c(0.2, NA, 0.8))

  bw <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.25, 0.25),
    bandwidth.compute = FALSE
  )
  fit <- npqreg(
    bws = bw,
    txdat = xdat,
    tydat = ydat,
    exdat = exdat,
    tau = c(0.25, 0.5),
    gradients = TRUE
  )

  expect_equal(dim(fit$quantile), c(3L, 2L))
  expect_equal(dim(fit$quanterr), c(3L, 2L))
  expect_equal(dim(fit$quantgrad), c(3L, 1L, 2L))
  expect_equal(dim(fit$quantgerr), c(3L, 1L, 2L))
  expect_true(all(is.na(fit$quantile[2, ])))
  expect_true(all(is.na(fit$quanterr[2, ])))
  expect_true(all(is.na(fit$quantgrad[2, , ])))
  expect_true(all(is.na(fit$quantgerr[2, , ])))

  train.fit <- npqreg(
    bws = bw,
    txdat = xdat,
    tydat = ydat,
    tau = c(0.25, 0.5),
    gradients = TRUE
  )
  pred <- predict(train.fit, exdat = exdat)
  expect_equal(dim(pred), c(3L, 2L))
  expect_true(all(is.na(pred[2, ])))
})

test_that("npRmpi npqreg direct evaluation follows conditional distribution eval guards", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260703)
  xdat <- data.frame(
    g = factor(rep(c("a", "b"), each = 12L)),
    x = seq(0.05, 0.95, length.out = 24)
  )
  ydat <- data.frame(y = rnorm(nrow(xdat)))
  bw.cat <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.3, 0.5, 0.3),
    bandwidth.compute = FALSE
  )
  ex.new <- data.frame(
    g = factor(c("a", "c"), levels = c("a", "b", "c")),
    x = c(0.2, 0.4)
  )
  expect_s3_class(
    npqreg(bws = bw.cat, txdat = xdat, tydat = ydat, exdat = ex.new, tau = 0.5),
    "qregression"
  )

  xcon <- data.frame(x = seq(0, 1, length.out = 25))
  ycon <- data.frame(y = rnorm(nrow(xcon)))
  bw.bound <- npcdistbw(
    xdat = xcon,
    ydat = ycon,
    bws = c(0.25, 0.25),
    bandwidth.compute = FALSE,
    cxkerbound = "range"
  )
  expect_error(
    npqreg(bws = bw.bound, txdat = xcon, tydat = ycon,
           exdat = data.frame(x = 2), tau = 0.5),
    "cxkerbound"
  )
})
