beta_adaptive_nn_radius <- function(train, k) {
  vapply(seq_along(train), function(i) {
    distance <- abs(train - train[i])
    positive <- sort(distance[distance > 0])
    duplicate_count <- sum(distance == 0) - 1L
    if (duplicate_count >= k) {
      positive[1L]
    } else {
      positive[max(1L, k - duplicate_count)]
    }
  }, numeric(1L))
}

beta_adaptive_nn_oracle <- function(train, evaluation, k, lower, upper) {
  support_length <- upper - lower
  bandwidth <- beta_adaptive_nn_radius(train, k)
  weights <- vapply(evaluation, function(target) {
    target_unit <- (target - lower) / support_length
    observation_unit <- (train - lower) / support_length
    concentration <- (support_length / bandwidth)^2
    dbeta(observation_unit,
          1 + target_unit * concentration,
          1 + (1 - target_unit) * concentration) / support_length
  }, numeric(length(train)))
  list(bandwidth = bandwidth, weights = weights)
}

test_that("adaptive-NN beta PDF matches observation-specific bandwidths", {
  train <- data.frame(x = c(0.03, 0.12, 0.28, 0.51, 0.79, 0.96))
  evaluation <- data.frame(x = c(0.02, 0.18, 0.46, 0.73, 0.98))
  oracle <- beta_adaptive_nn_oracle(
    train$x, evaluation$x, k = 2L, lower = 0, upper = 1
  )

  fit <- npksum(
    bws = 2,
    txdat = train,
    exdat = evaluation,
    bwtype = "adaptive_nn",
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

  for (i in seq_along(train$x)) {
    singleton <- npksum(
      bws = oracle$bandwidth[i],
      txdat = train[i, , drop = FALSE],
      exdat = evaluation,
      bwtype = "fixed",
      ckertype = "beta",
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1,
      return.kernel.weights = TRUE
    )
    expect_identical(fit$kw[i, ], as.double(singleton$kw[1L, ]))
  }
})

test_that("adaptive-NN beta PDF preserves ties and affine units", {
  train_unit <- c(0.08, 0.08, 0.31, 0.57, 0.84, 0.97)
  eval_unit <- c(0.02, 0.22, 0.63, 0.99)
  lower <- -2
  upper <- 4
  support_length <- upper - lower

  unit <- npksum(
    bws = 2,
    txdat = data.frame(x = train_unit),
    exdat = data.frame(x = eval_unit),
    bwtype = "adaptive_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    return.kernel.weights = TRUE
  )
  transformed <- npksum(
    bws = 2,
    txdat = data.frame(x = lower + support_length * train_unit),
    exdat = data.frame(x = lower + support_length * eval_unit),
    bwtype = "adaptive_nn",
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

  oracle <- beta_adaptive_nn_oracle(
    train_unit, eval_unit, k = 2L, lower = 0, upper = 1
  )
  expect_equal(unit$kw, oracle$weights, tolerance = 2e-13)
})

test_that("adaptive-NN beta PDF preserves leave-one-out and extended-NN policy", {
  train <- data.frame(x = c(0.05, 0.19, 0.38, 0.64, 0.88))
  fit <- npksum(
    bws = 2,
    txdat = train,
    bwtype = "adaptive_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    leave.one.out = TRUE,
    return.kernel.weights = TRUE
  )
  expect_identical(diag(fit$kw), rep(0, nrow(train)))
  expect_equal(as.double(fit$ksum), colSums(fit$kw), tolerance = 0)

  old <- options(np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  expect_error(
    npksum(
      bws = 7,
      txdat = train,
      exdat = data.frame(x = c(0.2, 0.7)),
      bwtype = "adaptive_nn",
      ckertype = "beta",
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1
    ),
    "invalid beta nearest-neighbor bandwidth"
  )

  options(np.extendednn = TRUE)
  extended <- npksum(
    bws = 7,
    txdat = train,
    exdat = data.frame(x = c(0.2, 0.7)),
    bwtype = "adaptive_nn",
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    return.kernel.weights = TRUE
  )
  base <- beta_adaptive_nn_radius(train$x, nrow(train) - 1L)
  bandwidth <- base * 7 / (nrow(train) - 1L)
  oracle <- vapply(c(0.2, 0.7), function(target) {
    concentration <- 1 / bandwidth^2
    dbeta(train$x,
          1 + target * concentration,
          1 + (1 - target) * concentration)
  }, numeric(nrow(train)))
  expect_equal(extended$kw, oracle, tolerance = 2e-13)
})
