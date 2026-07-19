beta_generalized_nn_oracle <- function(train, evaluation, k, lower, upper) {
  support_length <- upper - lower
  bandwidth <- vapply(evaluation, function(target) {
    sort(abs(train - target), partial = k)[k]
  }, numeric(1L))

  weights <- vapply(seq_along(evaluation), function(j) {
    target_unit <- (evaluation[j] - lower) / support_length
    observation_unit <- (train - lower) / support_length
    concentration <- (support_length / bandwidth[j])^2
    dbeta(observation_unit,
          1 + target_unit * concentration,
          1 + (1 - target_unit) * concentration) / support_length
  }, numeric(length(train)))

  list(bandwidth = bandwidth, weights = weights)
}

test_that("generalized-NN beta PDF matches an independent distance oracle", {
  train <- data.frame(x = c(0.05, 0.15, 0.30, 0.55, 0.80, 0.95))
  evaluation <- data.frame(x = c(0.10, 0.40, 0.90))
  oracle <- beta_generalized_nn_oracle(
    train$x, evaluation$x, k = 2L, lower = 0, upper = 1
  )

  fit <- npksum(
    bws = 2,
    txdat = train,
    exdat = evaluation,
    bwtype = "generalized_nn",
    ckertype = "beta",
    ckerorder = 2,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    return.kernel.weights = TRUE
  )

  expect_equal(fit$kw, oracle$weights, tolerance = 2e-13)
  expect_equal(as.double(fit$ksum), colSums(oracle$weights),
               tolerance = 2e-13)

  for (j in seq_along(evaluation$x)) {
    fixed <- npksum(
      bws = oracle$bandwidth[j],
      txdat = train,
      exdat = evaluation[j, , drop = FALSE],
      bwtype = "fixed",
      ckertype = "beta",
      ckerorder = 2,
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1,
      return.kernel.weights = TRUE
    )
    expect_identical(fit$kw[, j], fixed$kw[, 1L])
  }
})

test_that("generalized-NN beta PDF preserves affine support units", {
  train_unit <- c(0.04, 0.18, 0.33, 0.58, 0.77, 0.96)
  eval_unit <- c(0.08, 0.42, 0.88)
  lower <- -3
  upper <- 5
  support_length <- upper - lower

  unit <- npksum(
    bws = 3,
    txdat = data.frame(x = train_unit),
    exdat = data.frame(x = eval_unit),
    bwtype = "generalized_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    return.kernel.weights = TRUE
  )
  transformed <- npksum(
    bws = 3,
    txdat = data.frame(x = lower + support_length * train_unit),
    exdat = data.frame(x = lower + support_length * eval_unit),
    bwtype = "generalized_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = lower,
    ckerub = upper,
    return.kernel.weights = TRUE
  )

  expect_equal(transformed$kw * support_length, unit$kw,
               tolerance = 2e-13)
  expect_equal(as.double(transformed$ksum) * support_length,
               as.double(unit$ksum), tolerance = 2e-13)
})

test_that("generalized-NN beta PDF preserves observation-support tie handling", {
  train <- data.frame(x = c(0.10, 0.10, 0.40, 0.80, 0.95))
  evaluation <- data.frame(x = c(0.10, 0.40, 0.50))
  expected_bandwidth <- c(0.30, 0.30, 0.30)

  fit <- npksum(
    bws = 2,
    txdat = train,
    exdat = evaluation,
    bwtype = "generalized_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    return.kernel.weights = TRUE
  )

  for (j in seq_along(expected_bandwidth)) {
    fixed <- npksum(
      bws = expected_bandwidth[j],
      txdat = train,
      exdat = evaluation[j, , drop = FALSE],
      bwtype = "fixed",
      ckertype = "beta",
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1,
      return.kernel.weights = TRUE
    )
    expect_equal(fit$kw[, j], fixed$kw[, 1L], tolerance = 2e-14)
  }

  constant <- npksum(
    bws = 2,
    txdat = data.frame(x = rep(0.4, 5)),
    exdat = data.frame(x = 0.4),
    bwtype = "generalized_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    return.kernel.weights = TRUE
  )
  expect_identical(constant$kw[, 1L], rep(1, 5))
  expect_identical(as.double(constant$ksum), 5)
})

test_that("generalized-NN beta PDF retains the extended-NN policy", {
  train <- data.frame(x = c(0.05, 0.20, 0.45, 0.70, 0.95))
  evaluation <- data.frame(x = c(0.10, 0.60))

  old <- options(np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  expect_error(
    npksum(
      bws = 6,
      txdat = train,
      exdat = evaluation,
      bwtype = "generalized_nn",
      ckertype = "beta",
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1
    ),
    "invalid beta nearest-neighbor bandwidth"
  )

  options(np.extendednn = TRUE)
  extended <- npksum(
    bws = 6,
    txdat = train,
    exdat = evaluation,
    bwtype = "generalized_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    return.kernel.weights = TRUE
  )

  base_k <- nrow(train) - 1L
  base_radius <- vapply(evaluation$x, function(target) {
    sort(abs(train$x - target), partial = base_k)[base_k]
  }, numeric(1L))
  scaled_bandwidth <- base_radius * 6 / base_k
  oracle <- vapply(seq_along(evaluation$x), function(j) {
    target_unit <- evaluation$x[j]
    concentration <- 1 / scaled_bandwidth[j]^2
    dbeta(train$x,
          1 + target_unit * concentration,
          1 + (1 - target_unit) * concentration)
  }, numeric(nrow(train)))

  expect_equal(extended$kw, oracle, tolerance = 2e-13)
})
