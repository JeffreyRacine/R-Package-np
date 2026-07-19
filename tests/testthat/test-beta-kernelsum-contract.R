test_that("second-order beta npksum matches its associated-kernel formula", {
  training <- data.frame(x = c(0, 0.03, 0.2, 0.55, 0.91, 1))
  evaluation <- data.frame(x = c(0, 0.08, 0.5, 0.96, 1))
  h <- 0.14
  tau <- (1 / h)^2

  fit <- npksum(
    bws = h,
    txdat = training,
    exdat = evaluation,
    ckertype = "beta",
    ckerorder = 2,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    return.kernel.weights = TRUE
  )
  expected <- vapply(evaluation$x, function(target) {
    dbeta(training$x,
          shape1 = 1 + target * tau,
          shape2 = 1 + (1 - target) * tau)
  }, numeric(nrow(training)))

  expect_equal(fit$kw, expected, tolerance = 2e-12)
  expect_equal(as.double(fit$ksum), colSums(expected), tolerance = 2e-12)
  expect_true(all(fit$kw >= 0))
})

test_that("beta npksum carries full support scaling and composes in products", {
  training <- data.frame(
    x1 = c(0, 0.11, 0.43, 0.82, 1),
    x2 = c(-3, -1, 2, 5, 7)
  )
  evaluation <- data.frame(x1 = c(0, 0.5, 1), x2 = c(-3, 2, 7))

  fit <- npksum(
    bws = c(0.12, 1.2), txdat = training, exdat = evaluation,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = c(0, -3), ckerub = c(1, 7),
    return.kernel.weights = TRUE
  )
  unit_shifted <- data.frame(x1 = training$x1, x2 = (training$x2 + 3) / 10)
  eval_shifted <- data.frame(x1 = evaluation$x1,
                             x2 = (evaluation$x2 + 3) / 10)
  fit_unit <- npksum(
    bws = c(0.12, 0.12), txdat = unit_shifted, exdat = eval_shifted,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1),
    return.kernel.weights = TRUE
  )

  expect_equal(fit_unit$kw, 10 * fit$kw, tolerance = 2e-12)
  expect_equal(as.double(fit_unit$ksum),
               10 * as.double(fit$ksum), tolerance = 2e-12)
})

test_that("beta npksum supports weighting and leave-one-out", {
  training <- data.frame(x = seq(0, 1, length.out = 9))
  response <- cbind(y1 = seq_len(9), y2 = cos(seq_len(9)))
  case_weights <- cbind(w1 = 1, w2 = seq(0.5, 1.3, length.out = 9))
  args <- list(
    bws = 0.11, txdat = training,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )

  loo <- do.call(npksum, c(args, list(
    tydat = response,
    weights = case_weights,
    leave.one.out = TRUE,
    return.kernel.weights = TRUE
  )))
  expect_true(all(diag(loo$kw) == 0))
  for (j in seq_len(nrow(training)))
    for (y in seq_len(ncol(response)))
      for (w in seq_len(ncol(case_weights)))
        expect_equal(
          loo$ksum[w, y, j],
          sum(loo$kw[, j] * response[, y] * case_weights[, w]),
          tolerance = 2e-12
        )
})

test_that("unsupported beta npksum surfaces fail explicitly", {
  args <- list(
    bws = 0.1,
    txdat = data.frame(x = c(0, 0.5, 1)),
    ckertype = "beta",
    ckerorder = 2,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  expect_error(do.call(npksum, c(args, list(bandwidth.divide = TRUE))),
               "bandwidth.divide = TRUE", fixed = TRUE)
  expect_error(do.call(npksum, c(args, list(operator = "derivative"))),
               "only operator = \"normal\", \"convolution\", or \"integral\"",
               fixed = TRUE)
  expect_error(do.call(npksum, c(args, list(kernel.pow = 2))),
               "require kernel.pow = 1", fixed = TRUE)
  expect_error(npksum(0.1, args$txdat, ckertype = "beta",
                     ckerbound = "range"),
               "require ckerbound = \"fixed\"", fixed = TRUE)
})
