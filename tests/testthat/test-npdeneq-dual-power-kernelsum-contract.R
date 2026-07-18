test_that("internal dual-power kernel sums equal separate scalar calls", {
  skip_if_not(isTRUE(getOption("npRmpi.pool.active", FALSE)))
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  local_eval <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")

  set.seed(20260717)
  x <- data.frame(
    c = rnorm(48),
    u = factor(sample(letters[1:3], 48, TRUE)),
    o = ordered(sample(1:4, 48, TRUE))
  )
  bw <- c(0.6, 0.2, 0.3)

  scalar1 <- local_eval(npksum(
    bws = bw, txdat = x, leave.one.out = TRUE,
    bandwidth.divide = TRUE))
  scalar2 <- local_eval(npksum(
    bws = bw, txdat = x, kernel.pow = 2,
    leave.one.out = TRUE, bandwidth.divide = TRUE))
  dual_sum <- getFromNamespace(".npksum_power12", "npRmpi")
  dual <- local_eval(dual_sum(
    bws = bw, txdat = x, leave.one.out = TRUE,
    bandwidth.divide = TRUE))

  expect_identical(dual$ksum, scalar1$ksum)
  expect_identical(dual$ksum.power2, scalar2$ksum)
  expect_false("ksum.power2" %in% names(scalar1))
})

test_that("dual-power route is exact for nearest-neighbor bandwidths", {
  skip_if_not(isTRUE(getOption("npRmpi.pool.active", FALSE)))
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)
  dual_sum <- getFromNamespace(".npksum_power12", "npRmpi")
  local_eval <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")

  set.seed(20260718)
  x <- data.frame(x = rnorm(52))
  ex <- data.frame(x = rnorm(31))

  for (type in c("generalized_nn", "adaptive_nn")) {
    args <- list(bws = 8, txdat = x, exdat = ex, bwtype = type,
                 bandwidth.divide = TRUE)
    scalar1 <- local_eval(do.call(npksum, args))
    scalar2 <- local_eval(do.call(npksum, c(args, list(kernel.pow = 2))))
    dual <- local_eval(do.call(dual_sum, args))

    expect_identical(dual$ksum, scalar1$ksum)
    expect_identical(dual$ksum.power2, scalar2$ksum)
  }
})

test_that("dual-power entry rejects unsupported kernel-sum configurations", {
  skip_if_not(isTRUE(getOption("npRmpi.pool.active", FALSE)))
  default_sum <- getFromNamespace("npksum.default", "npRmpi")
  local_eval <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")
  set.seed(20260719)
  x <- data.frame(x = rnorm(20))

  expect_error(
    local_eval(default_sum(
      bws = 0.5, txdat = x, tydat = rnorm(20),
      .np.internal.power12 = TRUE)),
    "invalid use"
  )
  expect_error(
    local_eval(default_sum(
      bws = 0.5, txdat = x, operator = "integral",
      .np.internal.power12 = TRUE)),
    "invalid use"
  )
})
