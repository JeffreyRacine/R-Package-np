beta_generalized_nn_cdf_radius <- function(train, evaluation, k) {
  vapply(evaluation, function(target) {
    distance <- abs(train - target)
    exact_count <- sum(distance == 0)
    positive <- sort(distance[distance > 0])
    if (exact_count >= k) positive[1L] else positive[k - exact_count]
  }, numeric(1L))
}

beta_generalized_nn_cdf_oracle <- function(train, evaluation, k,
                                           lower, upper) {
  support_length <- upper - lower
  observation_unit <- (train - lower) / support_length
  target_unit <- (evaluation - lower) / support_length
  bandwidth <- beta_generalized_nn_cdf_radius(train, evaluation, k)

  vapply(seq_along(evaluation), function(j) {
    concentration <- (support_length / bandwidth[j])^2
    pbeta(target_unit[j],
          1 + observation_unit * concentration,
          1 + (1 - observation_unit) * concentration)
  }, numeric(length(train)))
}

test_that("generalized-NN beta CDF matches target-specific incomplete betas", {
  train <- data.frame(x = c(0.02, 0.09, 0.21, 0.43, 0.68, 0.87, 0.98))
  evaluation <- data.frame(x = c(0, 0.01, 0.13, 0.37, 0.72, 0.95, 1))
  oracle <- beta_generalized_nn_cdf_oracle(
    train$x, evaluation$x, k = 3L, lower = 0, upper = 1
  )

  fit <- npksum(
    bws = 3,
    txdat = train,
    exdat = evaluation,
    bwtype = "generalized_nn",
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
  expect_true(all(fit$kw >= 0 & fit$kw <= 1))
})

test_that("generalized-NN beta CDF exposes rather than repairs nonmonotonicity", {
  train <- data.frame(x = seq(0.086, 0.146, length.out = 12L))
  evaluation <- data.frame(x = seq(0, 1, length.out = 2001L))
  oracle <- beta_generalized_nn_cdf_oracle(
    train$x, evaluation$x, k = 4L, lower = 0, upper = 1
  )
  fit <- npksum(
    bws = 4,
    txdat = train,
    exdat = evaluation,
    bwtype = "generalized_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    operator = "integral"
  )
  estimate <- as.double(fit$ksum) / nrow(train)
  oracle_estimate <- colMeans(oracle)

  expect_equal(estimate, oracle_estimate, tolerance = 3e-13)
  expect_lt(min(diff(estimate)), -1e-5)
  expect_true(all(estimate >= 0 & estimate <= 1))
  expect_identical(estimate[c(1L, length(estimate))], c(0, 1))
})

test_that("generalized-NN beta CDF preserves affine support units", {
  train <- c(0.03, 0.14, 0.31, 0.57, 0.76, 0.94)
  evaluation <- c(0, 0.08, 0.39, 0.81, 1)
  lower <- -7
  support_length <- 13
  unit <- npksum(
    bws = 2, txdat = data.frame(x = train),
    exdat = data.frame(x = evaluation), bwtype = "generalized_nn",
    ckertype = "beta", ckerbound = "fixed", ckerlb = 0, ckerub = 1,
    operator = "integral", return.kernel.weights = TRUE
  )
  transformed <- npksum(
    bws = 2, txdat = data.frame(x = lower + support_length * train),
    exdat = data.frame(x = lower + support_length * evaluation),
    bwtype = "generalized_nn", ckertype = "beta", ckerbound = "fixed",
    ckerlb = lower, ckerub = lower + support_length,
    operator = "integral", return.kernel.weights = TRUE
  )
  expect_equal(transformed$kw, unit$kw, tolerance = 3e-13)
  expect_equal(as.double(transformed$ksum), as.double(unit$ksum),
               tolerance = 3e-13)
})
