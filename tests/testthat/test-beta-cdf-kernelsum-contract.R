test_that("order-2 beta integral operator matches observation-centred CDF", {
  training <- data.frame(x = c(0, 0.01, 0.12, 0.4, 0.73, 0.96, 1))
  evaluation <- data.frame(x = c(0, 1e-12, 0.04, 0.5, 0.98, 1 - 1e-12, 1))
  h <- 0.13
  tau <- (1 / h)^2

  fit <- npksum(
    bws = h, txdat = training, exdat = evaluation,
    operator = "integral", return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  expected <- vapply(evaluation$x, function(target) {
    pbeta(target,
          shape1 = 1 + training$x * tau,
          shape2 = 1 + (1 - training$x) * tau)
  }, numeric(nrow(training)))

  expect_equal(fit$kw, expected, tolerance = 2e-12)
  expect_equal(as.double(fit$ksum), colSums(expected), tolerance = 2e-12)
  expect_identical(fit$kw[, 1L], rep(0, nrow(training)))
  expect_identical(fit$kw[, nrow(evaluation)], rep(1, nrow(training)))
  expect_true(all(fit$kw >= 0 & fit$kw <= 1))
  expect_true(all(apply(fit$kw, 1L, function(value) all(diff(value) >= 0))))
})

test_that("beta normal and integral operators compose dimension by dimension", {
  training <- data.frame(
    x = c(0, 0.08, 0.3, 0.64, 0.91, 1),
    y = c(-3, -2.6, -0.2, 2.8, 6.4, 7)
  )
  evaluation <- data.frame(
    x = c(0, 0.2, 0.7, 1),
    y = c(-3, -1, 4, 7)
  )
  bandwidth <- c(0.14, 1.1)
  tau <- c((1 / bandwidth[1L])^2, (10 / bandwidth[2L])^2)

  fit <- npksum(
    bws = bandwidth, txdat = training, exdat = evaluation,
    operator = c("normal", "integral"), return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = c(0, -3), ckerub = c(1, 7)
  )
  expected <- vapply(seq_len(nrow(evaluation)), function(index) {
    pdf <- dbeta(training$x,
                 1 + evaluation$x[index] * tau[1L],
                 1 + (1 - evaluation$x[index]) * tau[1L])
    cdf <- pbeta((evaluation$y[index] + 3) / 10,
                 1 + ((training$y + 3) / 10) * tau[2L],
                 1 + ((7 - training$y) / 10) * tau[2L])
    pdf * cdf
  }, numeric(nrow(training)))

  expect_equal(fit$kw, expected, tolerance = 2e-12)
  expect_equal(as.double(fit$ksum), colSums(expected), tolerance = 2e-12)
})

test_that("beta CDF is affine invariant and has the uniform large-h limit", {
  training <- data.frame(x = c(0, 0.05, 0.2, 0.55, 0.9, 1))
  evaluation <- data.frame(x = c(0, 0.03, 0.3, 0.8, 1))
  unit <- npksum(
    bws = 0.12, txdat = training, exdat = evaluation,
    operator = "integral", return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  shifted <- npksum(
    bws = 1.2, txdat = data.frame(x = -3 + 10 * training$x),
    exdat = data.frame(x = -3 + 10 * evaluation$x),
    operator = "integral", return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = -3, ckerub = 7
  )
  uniform_limit <- npksum(
    bws = 1e200, txdat = training, exdat = evaluation,
    operator = "integral", return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )

  expect_equal(shifted$kw, unit$kw, tolerance = 2e-12)
  expect_equal(as.double(shifted$ksum), as.double(unit$ksum),
               tolerance = 2e-12)
  expect_equal(uniform_limit$kw,
               matrix(rep(evaluation$x, each = nrow(training)),
                      nrow = nrow(training)),
               tolerance = 2e-15)
})

test_that("beta CDF supports responses, case weights, and leave-one-out", {
  training <- data.frame(x = seq(0, 1, length.out = 7))
  response <- cbind(y1 = seq_len(7), y2 = cos(seq_len(7)))
  case_weights <- cbind(w1 = 1, w2 = seq(0.4, 1, length.out = 7))
  fit <- npksum(
    bws = 0.15, txdat = training, tydat = response,
    weights = case_weights, leave.one.out = TRUE,
    operator = "integral", return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )

  expect_true(all(diag(fit$kw) == 0))
  for (index in seq_len(nrow(training)))
    for (response_column in seq_len(ncol(response)))
      for (weight_column in seq_len(ncol(case_weights)))
        expect_equal(
          fit$ksum[weight_column, response_column, index],
          sum(fit$kw[, index] * response[, response_column] *
                case_weights[, weight_column]),
          tolerance = 2e-12
        )
})
