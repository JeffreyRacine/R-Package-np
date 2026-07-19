beta_adaptive_nn_cdf_radius <- function(train, k) {
  vapply(seq_along(train), function(i) {
    distance <- abs(train - train[i])
    duplicate_count <- sum(distance == 0) - 1L
    positive <- sort(distance[distance > 0])
    if (duplicate_count >= k) positive[1L] else
      positive[max(1L, k - duplicate_count)]
  }, numeric(1L))
}

beta_adaptive_nn_cdf_oracle <- function(train, evaluation, k, lower, upper) {
  support_length <- upper - lower
  target_unit <- (evaluation - lower) / support_length
  observation_unit <- (train - lower) / support_length
  bandwidth <- beta_adaptive_nn_cdf_radius(train, k)
  concentration <- (support_length / bandwidth)^2

  vapply(target_unit, function(target) {
    pbeta(target,
          1 + observation_unit * concentration,
          1 + (1 - observation_unit) * concentration)
  }, numeric(length(train)))
}

test_that("adaptive-NN beta CDF matches observation-specific incomplete betas", {
  train <- data.frame(x = c(0.02, 0.09, 0.21, 0.43, 0.68, 0.87, 0.98))
  evaluation <- data.frame(x = c(0, 0.01, 0.13, 0.37, 0.72, 0.95, 1))
  oracle <- beta_adaptive_nn_cdf_oracle(
    train$x, evaluation$x, k = 3L, lower = 0, upper = 1
  )

  fit <- npksum(
    bws = 3,
    txdat = train,
    exdat = evaluation,
    bwtype = "adaptive_nn",
    ckertype = "beta",
    ckerorder = 2,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    operator = "integral",
    return.kernel.weights = TRUE
  )

  expect_equal(fit$kw, oracle, tolerance = 3e-13)
  expect_equal(as.double(fit$ksum), colSums(oracle), tolerance = 3e-13)
  expect_identical(fit$kw[, 1L], rep(0, nrow(train)))
  expect_identical(fit$kw[, ncol(fit$kw)], rep(1, nrow(train)))
  expect_true(all(apply(fit$kw, 1L, diff) >= 0))
})

test_that("adaptive-NN beta CDF preserves affine support and mixed operators", {
  train <- data.frame(
    x = c(0.03, 0.15, 0.38, 0.62, 0.89, 0.97),
    y = c(0.06, 0.22, 0.41, 0.71, 0.84, 0.99)
  )
  evaluation <- data.frame(
    x = c(0, 0.17, 0.55, 0.93, 1),
    y = c(0.01, 0.29, 0.64, 0.91, 0.995)
  )
  lower <- c(-4, 20)
  support_length <- c(9, 0.05)

  unit <- npksum(
    bws = c(2, 3), txdat = train, exdat = evaluation,
    bwtype = "adaptive_nn", ckertype = "beta", ckerbound = "fixed",
    ckerlb = c(0, 0), ckerub = c(1, 1),
    operator = c("integral", "normal"), return.kernel.weights = TRUE
  )
  transformed <- npksum(
    bws = c(2, 3),
    txdat = data.frame(x = lower[1] + support_length[1] * train$x,
                       y = lower[2] + support_length[2] * train$y),
    exdat = data.frame(x = lower[1] + support_length[1] * evaluation$x,
                       y = lower[2] + support_length[2] * evaluation$y),
    bwtype = "adaptive_nn", ckertype = "beta", ckerbound = "fixed",
    ckerlb = lower, ckerub = lower + support_length,
    operator = c("integral", "normal"), return.kernel.weights = TRUE
  )

  expect_equal(transformed$kw * support_length[2], unit$kw,
               tolerance = 3e-12)
  expect_equal(as.double(transformed$ksum) * support_length[2],
               as.double(unit$ksum), tolerance = 3e-12)
})

test_that("adaptive-NN order-2 beta CDF remains monotone on boundary designs", {
  designs <- list(
    lower_heavy = qchisq(seq(0.01, 0.99, length.out = 41), df = 1) / 8,
    upper_heavy = 1 - qchisq(seq(0.01, 0.99, length.out = 41), df = 1) / 8
  )
  grid <- data.frame(x = seq(0, 1, length.out = 501))

  for (train in designs) {
    train <- pmin(1, pmax(0, train))
    fit <- npksum(
      bws = 5,
      txdat = data.frame(x = train),
      exdat = grid,
      bwtype = "adaptive_nn",
      ckertype = "beta",
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1,
      operator = "integral"
    )
    estimate <- as.double(fit$ksum) / length(train)
    expect_true(all(diff(estimate) >= -1e-15))
    expect_true(all(estimate >= 0 & estimate <= 1))
    expect_identical(estimate[c(1L, length(estimate))], c(0, 1))
  }
})
