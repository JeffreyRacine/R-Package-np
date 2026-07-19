test_that("manual order-2 beta density matches npksum and exact contribution SE", {
  training <- data.frame(x = c(0, 0.01, 0.04, 0.2, 0.55, 0.9, 1))
  evaluation <- data.frame(x = seq(0, 1, length.out = 21))
  args <- list(
    bws = 0.12,
    ckertype = "beta",
    ckerorder = 2,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  fit <- do.call(npudens, c(list(tdat = training, edat = evaluation), args))
  sums <- do.call(npksum, c(list(
    txdat = training,
    exdat = evaluation,
    return.kernel.weights = TRUE
  ), args))
  expected_density <- colMeans(sums$kw)
  expected_se <- sqrt(apply(sums$kw, 2L, var) / nrow(training))

  expect_equal(fitted(fit), expected_density, tolerance = 2e-12)
  expect_equal(se(fit), expected_se, tolerance = 2e-12)
  expect_equal(
    fit$log_likelihood,
    sum(log(pmax(expected_density, .Machine$double.xmin))),
    tolerance = 2e-12
  )
  expect_true(all(fitted(fit) >= 0))
})

test_that("manual beta density supports formulas, objects, and dimensions", {
  training <- data.frame(
    x1 = c(0, 0.03, 0.2, 0.5, 0.8, 0.97, 1),
    x2 = c(-3, -2.2, -0.4, 1.7, 4.2, 6.1, 7)
  )
  evaluation <- data.frame(
    x1 = c(0, 0.1, 0.5, 0.9, 1),
    x2 = c(-3, -1, 2, 5.7, 7)
  )
  args <- list(
    bws = c(0.13, 1.3),
    ckertype = "beta",
    ckerorder = 2,
    ckerbound = "fixed",
    ckerlb = c(0, -3),
    ckerub = c(1, 7)
  )

  direct <- do.call(npudens, c(list(tdat = training, edat = evaluation), args))
  formula <- do.call(npudens, c(list(
    ~ x1 + x2,
    data = training,
    newdata = evaluation
  ), args))
  bw <- do.call(npudensbw, c(list(
    dat = training,
    bandwidth.compute = FALSE
  ), args))
  object_fit <- npudens(bws = bw, tdat = training, edat = evaluation)

  expect_s3_class(bw, "bandwidth")
  expect_identical(bw$pmethod, "Manual")
  expect_identical(bw$ckertype, "beta")
  expect_equal(fitted(formula), fitted(direct), tolerance = 2e-12)
  expect_equal(fitted(object_fit), fitted(direct), tolerance = 2e-12)
  expect_equal(se(object_fit), se(direct), tolerance = 2e-12)
})

test_that("beta density keeps unsupported selection and operator routes closed", {
  training <- data.frame(x = c(0, 0.05, 0.25, 0.6, 1))

  expect_error(
    suppressWarnings(npudens(tdat = training, ckertype = "beta", ckerorder = 2,
                             ckerbound = "fixed", ckerlb = 0, ckerub = 1)),
    "does not yet support automatic bandwidth selection",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(npudens(tdat = training, bws = 0.1, ckertype = "beta",
                             ckerbound = "range")),
    "require ckerbound = \"fixed\"",
    fixed = TRUE
  )
  expect_error(
    npudens(tdat = training, edat = data.frame(x = 1.1), bws = 0.1,
            ckertype = "beta", ckerbound = "fixed",
            ckerlb = 0, ckerub = 1),
    "Evaluation data violate",
    fixed = TRUE
  )
})
