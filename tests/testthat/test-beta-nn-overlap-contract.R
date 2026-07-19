beta_nn_overlap_oracle <- function(center_one, bandwidth_one,
                                   center_two, bandwidth_two,
                                   lower, upper) {
  support_length <- upper - lower
  unit_one <- (center_one - lower) / support_length
  unit_two <- (center_two - lower) / support_length
  concentration_one <- (support_length / bandwidth_one)^2
  concentration_two <- (support_length / bandwidth_two)^2
  alpha_one <- 1 + unit_one * concentration_one
  beta_one <- 1 + (1 - unit_one) * concentration_one
  alpha_two <- 1 + unit_two * concentration_two
  beta_two <- 1 + (1 - unit_two) * concentration_two

  exp(lbeta(alpha_one + alpha_two - 1, beta_one + beta_two - 1) -
        lbeta(alpha_one, beta_one) - lbeta(alpha_two, beta_two)) /
    support_length
}

beta_nn_adaptive_radius <- function(train, k) {
  vapply(seq_along(train), function(i) {
    distance <- abs(train - train[i])
    duplicate_count <- sum(distance == 0) - 1L
    positive <- sort(distance[distance > 0])
    if (duplicate_count >= k) positive[1L] else
      positive[max(1L, k - duplicate_count)]
  }, numeric(1L))
}

beta_nn_generalized_radius <- function(train, evaluation, k) {
  vapply(evaluation, function(target) {
    distance <- abs(train - target)
    exact_count <- sum(distance == 0)
    positive <- sort(distance[distance > 0])
    if (exact_count >= k) positive[1L] else
      positive[k - exact_count]
  }, numeric(1L))
}

test_that("nearest-neighbor beta overlap uses both center-specific bandwidths", {
  train <- c(0.04, 0.16, 0.34, 0.59, 0.81, 0.96)
  evaluation <- c(0.02, 0.24, 0.53, 0.88)
  k <- 2L
  eval_bandwidth <- beta_nn_generalized_radius(train, evaluation, k)
  generalized_train_bandwidth <-
    beta_nn_generalized_radius(train, train, k)
  adaptive_train_bandwidth <- beta_nn_adaptive_radius(train, k)
  generalized_oracle <- vapply(seq_along(evaluation), function(j) {
    beta_nn_overlap_oracle(evaluation[j], eval_bandwidth[j],
                           train, generalized_train_bandwidth, 0, 1)
  }, numeric(length(train)))
  adaptive_oracle <- vapply(seq_along(evaluation), function(j) {
    beta_nn_overlap_oracle(evaluation[j], eval_bandwidth[j],
                           train, adaptive_train_bandwidth, 0, 1)
  }, numeric(length(train)))

  generalized <- npksum(
    bws = k,
    txdat = data.frame(x = train),
    exdat = data.frame(x = evaluation),
    bwtype = "generalized_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    operator = "convolution",
    return.kernel.weights = TRUE
  )
  adaptive <- npksum(
    bws = k,
    txdat = data.frame(x = train),
    exdat = data.frame(x = evaluation),
    bwtype = "adaptive_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    operator = "convolution",
    return.kernel.weights = TRUE
  )

  expect_equal(generalized$kw, generalized_oracle, tolerance = 2e-12)
  expect_equal(adaptive$kw, adaptive_oracle, tolerance = 2e-12)
  expect_false(isTRUE(all.equal(generalized$kw, adaptive$kw)))
  expect_equal(as.double(generalized$ksum), colSums(generalized_oracle),
               tolerance = 2e-12)
})

test_that("nearest-neighbor beta overlap is symmetric on training centers", {
  train <- data.frame(x = c(0.03, 0.17, 0.36, 0.62, 0.83, 0.98))
  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    fit <- npksum(
      bws = 3,
      txdat = train,
      bwtype = bwtype,
      ckertype = "beta",
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1,
      operator = "convolution",
      return.kernel.weights = TRUE
    )
    expect_equal(fit$kw, t(fit$kw), tolerance = 2e-13)
    expect_true(all(diag(fit$kw) > 0))
  }
})

test_that("nearest-neighbor beta overlap preserves affine units and mixed products", {
  train <- data.frame(
    x = c(0.05, 0.18, 0.39, 0.67, 0.91),
    y = c(0.08, 0.29, 0.44, 0.73, 0.97)
  )
  evaluation <- data.frame(
    x = c(0.02, 0.31, 0.84),
    y = c(0.13, 0.58, 0.93)
  )
  lower <- c(-2, 10)
  length <- c(5, 0.04)

  unit <- npksum(
    bws = c(2, 3),
    txdat = train,
    exdat = evaluation,
    bwtype = "adaptive_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = c(0, 0),
    ckerub = c(1, 1),
    operator = c("normal", "convolution"),
    return.kernel.weights = TRUE
  )
  transformed <- npksum(
    bws = c(2, 3),
    txdat = data.frame(x = lower[1] + length[1] * train$x,
                       y = lower[2] + length[2] * train$y),
    exdat = data.frame(x = lower[1] + length[1] * evaluation$x,
                       y = lower[2] + length[2] * evaluation$y),
    bwtype = "adaptive_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = lower,
    ckerub = lower + length,
    operator = c("normal", "convolution"),
    return.kernel.weights = TRUE
  )

  expect_equal(transformed$kw * prod(length), unit$kw,
               tolerance = 2e-12)
  expect_equal(as.double(transformed$ksum) * prod(length),
               as.double(unit$ksum), tolerance = 2e-12)
})
