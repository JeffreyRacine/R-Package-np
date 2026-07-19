beta_overlap_expected <- function(center_one, bandwidth_one,
                                  center_two, bandwidth_two,
                                  lower, upper) {
  support_length <- upper - lower
  unit_one <- (center_one - lower) / support_length
  unit_two <- (center_two - lower) / support_length
  tau_one <- (support_length / bandwidth_one)^2
  tau_two <- (support_length / bandwidth_two)^2
  alpha_one <- 1 + unit_one * tau_one
  beta_one <- 1 + (1 - unit_one) * tau_one
  alpha_two <- 1 + unit_two * tau_two
  beta_two <- 1 + (1 - unit_two) * tau_two
  exp(
    -log(support_length) +
      lbeta(alpha_one + alpha_two - 1, beta_one + beta_two - 1) -
      lbeta(alpha_one, beta_one) - lbeta(alpha_two, beta_two)
  )
}

test_that("order-2 beta convolution is the analytic beta-weight overlap", {
  training <- data.frame(x = c(0, 0.01, 0.17, 0.5, 0.83, 0.99, 1))
  evaluation <- data.frame(x = c(0, 0.04, 0.35, 0.7, 0.98, 1))
  bandwidth <- 0.14
  fit <- npksum(
    bws = bandwidth, txdat = training, exdat = evaluation,
    operator = "convolution", return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  expected <- vapply(evaluation$x, function(center) {
    beta_overlap_expected(training$x, bandwidth,
                          center, bandwidth, 0, 1)
  }, numeric(nrow(training)))

  expect_equal(fit$kw, expected, tolerance = 2e-12)
  expect_equal(as.double(fit$ksum), colSums(expected), tolerance = 2e-12)
  expect_true(all(is.finite(fit$kw) & fit$kw >= 0))
})

test_that("beta overlap is symmetric and its diagonal is exact roughness", {
  centers <- c(0, 0.03, 0.27, 0.5, 0.91, 1)
  bandwidth <- 0.17
  forward <- npksum(
    bws = bandwidth, txdat = data.frame(x = centers),
    exdat = data.frame(x = rev(centers)),
    operator = "convolution", return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )$kw
  reverse <- npksum(
    bws = bandwidth, txdat = data.frame(x = rev(centers)),
    exdat = data.frame(x = centers),
    operator = "convolution", return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )$kw
  expect_equal(forward, t(reverse), tolerance = 2e-12)

  diagonal <- vapply(centers, function(center) {
    npksum(
      bws = bandwidth, txdat = data.frame(x = center),
      exdat = data.frame(x = center), operator = "convolution",
      ckertype = "beta", ckerorder = 2,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )$ksum
  }, numeric(1L))
  integrated <- vapply(centers, function(center) {
    tau <- (1 / bandwidth)^2
    stats::integrate(
      function(value) dbeta(value, 1 + center * tau,
                            1 + (1 - center) * tau)^2,
      lower = 0, upper = 1, rel.tol = 1e-11
    )$value
  }, numeric(1L))
  expect_equal(diagonal, integrated, tolerance = 2e-10)
})

test_that("beta overlap has original-scale units and composes with other operators", {
  training <- data.frame(
    x = c(0, 0.08, 0.3, 0.64, 0.91, 1),
    y = c(-3, -2.6, -0.2, 2.8, 6.4, 7)
  )
  evaluation <- data.frame(
    x = c(0, 0.2, 0.7, 1),
    y = c(-3, -1, 4, 7)
  )
  operators <- c("integral", "convolution")
  unit <- npksum(
    bws = c(0.14, 0.11),
    txdat = data.frame(x = training$x, y = (training$y + 3) / 10),
    exdat = data.frame(x = evaluation$x, y = (evaluation$y + 3) / 10),
    operator = operators, return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
  )
  original <- npksum(
    bws = c(0.14, 1.1), txdat = training, exdat = evaluation,
    operator = operators, return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = c(0, -3), ckerub = c(1, 7)
  )

  expect_equal(unit$kw, 10 * original$kw, tolerance = 2e-12)
  expect_equal(as.double(unit$ksum), 10 * as.double(original$ksum),
               tolerance = 2e-12)
})

test_that("beta overlap supports responses, weights, and leave-one-out", {
  training <- data.frame(x = seq(0, 1, length.out = 7))
  response <- cbind(y1 = seq_len(7), y2 = cos(seq_len(7)))
  case_weights <- cbind(w1 = 1, w2 = seq(0.4, 1, length.out = 7))
  fit <- npksum(
    bws = 0.15, txdat = training, tydat = response,
    weights = case_weights, leave.one.out = TRUE,
    operator = "convolution", return.kernel.weights = TRUE,
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
