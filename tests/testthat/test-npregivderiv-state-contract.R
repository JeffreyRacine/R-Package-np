integrate_npregivderiv_state <- function(z, derivative, y_mean) {
  z <- as.numeric(z)
  derivative <- as.numeric(derivative)
  order_z <- order(z)
  z_sorted <- z[order_z]
  derivative_sorted <- derivative[order_z]
  dz <- diff(z_sorted)
  cz <- dz[[1L]]
  derivative_diff <- diff(derivative_sorted)
  correction <- cz^2 / 12 *
    (derivative_diff[[length(derivative_diff)]] / cz -
       derivative_diff[[1L]] / cz)
  if (!is.finite(correction)) correction <- 0
  integrated <- c(
    0,
    cumsum(dz * (derivative_sorted[-length(z_sorted)] +
                   derivative_sorted[-1L]) / 2)
  ) - correction
  integrated <- integrated[order(order_z)]
  integrated - mean(integrated) + y_mean
}

test_that("npregivderiv stop selection is boundary safe and deterministic", {
  select <- np:::.npregivderiv_select_stop_index

  expect_identical(select(1)$index, 1L)
  expect_false(select(1)$monotone.failure)
  expect_identical(select(c(1, 2, 3))$index, 1L)
  expect_true(select(c(1, 2, 3))$monotone.failure)
  expect_identical(select(c(3, 2, 1))$index, 3L)
  expect_identical(select(c(1, 3, 2, 1, 2))$index, 4L)
  expect_identical(select(c(1, 3, 2, 2, 3))$index, 3L)
})

test_that("npregivderiv columns and returned curves identify one state", {
  set.seed(20260720)
  n <- 48L
  w <- runif(n, -2, 2)
  v <- rnorm(n, sd = 0.15)
  z <- 0.4 * w + v
  y <- z^2 - 0.3 * v + rnorm(n, sd = 0.05)

  fit <- suppressWarnings(
    npregivderiv(
      y = y,
      z = z,
      w = w,
      starting.values = rep(0, n),
      iterate.max = 3L,
      iterate.break = FALSE,
      stop.on.increase = FALSE,
      nmulti = 1L,
      random.seed = 20260720
    )
  )

  expect_equal(length(fit$norm.stop), 3L)
  expect_equal(ncol(fit$phi.mat), length(fit$norm.stop))
  expect_equal(ncol(fit$phi.prime.mat), length(fit$norm.stop))
  expect_true(fit$num.iterations %in% seq_along(fit$norm.stop))
  expect_identical(as.numeric(fit$phi),
                   as.numeric(fit$phi.mat[, fit$num.iterations]))
  expect_identical(as.numeric(fit$phi.prime),
                   as.numeric(fit$phi.prime.mat[, fit$num.iterations]))
  expect_equal(mean(fit$starting.values.phi), mean(y), tolerance = 1e-14)

  for (N in seq_along(fit$norm.stop)) {
    expected_phi <- integrate_npregivderiv_state(
      z,
      fit$phi.prime.mat[, N],
      mean(y)
    )
    expect_equal(as.numeric(fit$phi.mat[, N]),
                 expected_phi,
                 tolerance = 1e-14)
  }

  summary_text <- capture.output(summary(fit))
  selected_value <- format(fit$norm.stop[fit$num.iterations], digits = 8)
  expect_true(any(grepl(selected_value, summary_text, fixed = TRUE)))
})
