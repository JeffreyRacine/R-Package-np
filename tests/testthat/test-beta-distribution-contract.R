test_that("manual order-2 beta distribution matches exact CDF contributions", {
  training <- data.frame(x = c(0, 0.01, 0.12, 0.4, 0.73, 0.96, 1))
  evaluation <- data.frame(x = c(0, 1e-12, 0.04, 0.5, 0.98, 1 - 1e-12, 1))
  h <- 0.13
  tau <- (1 / h)^2

  fit <- npudist(
    bws = h, tdat = training, edat = evaluation,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  contributions <- vapply(evaluation$x, function(target) {
    pbeta(target,
          1 + training$x * tau,
          1 + (1 - training$x) * tau)
  }, numeric(nrow(training)))
  expected <- colMeans(contributions)
  expected_se <- sqrt(colSums(sweep(contributions, 2L, expected, "-")^2) /
                        (nrow(training) * (nrow(training) - 1L)))

  expect_equal(fitted(fit), expected, tolerance = 2e-12)
  expect_equal(se(fit), expected_se, tolerance = 2e-12)
  expect_identical(fitted(fit)[1L], 0)
  expect_identical(fitted(fit)[nrow(evaluation)], 1)
  expect_true(all(diff(fitted(fit)) >= 0))
  expect_true(all(fitted(fit) >= 0 & fitted(fit) <= 1))
})

test_that("beta distribution supports dimensions, formulas, objects, and prediction", {
  training <- data.frame(
    x = c(0, 0.08, 0.3, 0.64, 0.91, 1),
    y = c(-3, -2.6, -0.2, 2.8, 6.4, 7)
  )
  evaluation <- data.frame(
    x = c(0, 0.2, 0.7, 1),
    y = c(-3, -1, 4, 7)
  )
  args <- list(
    bws = c(0.14, 1.1), ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = c(0, -3), ckerub = c(1, 7)
  )
  direct <- do.call(npudist, c(list(tdat = training, edat = evaluation), args))
  sums <- do.call(npksum, c(list(
    txdat = training, exdat = evaluation, operator = "integral",
    return.kernel.weights = TRUE
  ), args))
  bw <- do.call(npudistbw, c(list(
    dat = training, bandwidth.compute = FALSE
  ), args))
  object_fit <- npudist(bws = bw, tdat = training, edat = evaluation)
  formula_fit <- do.call(npudist, c(list(
    ~ x + y, data = training, newdata = evaluation
  ), args))
  predicted <- predict(formula_fit, newdata = evaluation, se.fit = TRUE)

  expected <- colMeans(sums$kw)
  expected_se <- sqrt(colSums(sweep(sums$kw, 2L, expected, "-")^2) /
                        (nrow(training) * (nrow(training) - 1L)))
  expect_s3_class(bw, "dbandwidth")
  expect_identical(bw$pmethod, "Manual")
  expect_identical(bw$ckertype, "beta")
  expect_equal(fitted(direct), expected, tolerance = 2e-12)
  expect_equal(se(direct), expected_se, tolerance = 2e-12)
  expect_equal(fitted(object_fit), fitted(direct), tolerance = 2e-12)
  expect_equal(se(object_fit), se(direct), tolerance = 2e-12)
  expect_equal(fitted(formula_fit), fitted(direct), tolerance = 2e-12)
  expect_equal(as.numeric(predicted$fit), fitted(direct), tolerance = 2e-12)
  expect_equal(as.numeric(predicted$se.fit), se(direct), tolerance = 2e-12)
})

test_that("unsupported beta distribution surfaces fail explicitly", {
  training <- data.frame(x = c(0, 0.03, 0.2, 0.6, 1))

  expect_error(
    suppressWarnings(npudist(
      bws = 0.1, tdat = training, ckertype = "beta",
      ckerbound = "range"
    )),
    "require ckerbound = \"fixed\"",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(npudist(
      bws = c(0.1, 0.2),
      tdat = data.frame(x = training$x,
                        group = ordered(c(1, 1, 1, 2, 2))),
      ckertype = "beta", ckerbound = "fixed",
      ckerlb = 0, ckerub = 1
    )),
    "continuous variables only",
    fixed = TRUE
  )
})
