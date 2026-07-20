test_that("manual order-2 beta local-constant regression matches exact weights", {
  training <- data.frame(
    x1 = c(0, 0.002, 0.02, 0.1, 0.32, 0.61, 0.86, 0.985, 1),
    x2 = c(-3, -2.95, -2.3, -0.5, 1.4, 3.8, 5.9, 6.9, 7)
  )
  response <- sin(3 * training$x1) + 0.2 * training$x2 + training$x1^2
  evaluation <- data.frame(
    x1 = c(0, 1e-12, 0.01, 0.5, 0.99, 1 - 1e-12, 1),
    x2 = c(-3, -3 + 1e-11, -2, 2, 6, 7 - 1e-11, 7)
  )
  args <- list(
    bws = c(0.11, 1.1),
    regtype = "lc",
    ckertype = "beta",
    ckerorder = 2,
    ckerbound = "fixed",
    ckerlb = c(0, -3),
    ckerub = c(1, 7)
  )

  fit <- do.call(npreg, c(list(
    txdat = training,
    tydat = response,
    exdat = evaluation
  ), args))
  sums <- do.call(npksum, c(list(
    txdat = training,
    exdat = evaluation,
    return.kernel.weights = TRUE
  ), args[setdiff(names(args), "regtype")]))
  normalized <- sweep(sums$kw, 2L, colSums(sums$kw), "/")
  expected_mean <- colSums(normalized * response)
  centered <- response - matrix(expected_mean,
                                nrow = nrow(training),
                                ncol = nrow(evaluation),
                                byrow = TRUE)
  expected_variance <- colSums(normalized * centered^2)
  expected_se <- sqrt(expected_variance * colSums(normalized^2))

  expect_equal(fitted(fit), expected_mean, tolerance = 2e-12)
  expect_lt(max(abs(se(fit) - expected_se)), 2e-12)
  expect_true(all(is.finite(fitted(fit))))
  expect_true(all(is.finite(se(fit))))
  expect_true(all(se(fit) >= 0))
})

test_that("beta regression supports formula, object, prediction, and residual routes", {
  training <- data.frame(
    x = c(0, 0.01, 0.04, 0.15, 0.4, 0.72, 0.91, 0.99, 1)
  )
  training$y <- cos(2 * training$x) + training$x
  evaluation <- data.frame(x = c(0, 0.02, 0.3, 0.8, 1))
  args <- list(
    bws = 0.12,
    regtype = "lc",
    ckertype = "beta",
    ckerorder = 2,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  direct <- do.call(npreg, c(list(
    txdat = training["x"],
    tydat = training$y,
    exdat = evaluation
  ), args))
  formula <- do.call(npreg, c(list(
    y ~ x,
    data = training,
    newdata = evaluation
  ), args))
  bw <- do.call(npregbw, c(list(
    xdat = training["x"],
    ydat = training$y,
    bandwidth.compute = FALSE
  ), args))
  object_fit <- npreg(bws = bw, txdat = training["x"],
                      tydat = training$y, exdat = evaluation)
  direct_helper <- np:::.np_regression_direct(
    bws = bw, txdat = training["x"], tydat = training$y,
    exdat = evaluation
  )
  in_sample <- npreg(bws = bw, txdat = training["x"],
                     tydat = training$y, residuals = TRUE)
  predicted <- predict(formula, newdata = evaluation, se.fit = TRUE)

  expect_s3_class(bw, "rbandwidth")
  expect_identical(bw$pmethod, "Manual")
  expect_identical(bw$ckertype, "beta")
  expect_equal(fitted(formula), fitted(direct), tolerance = 2e-12)
  expect_equal(fitted(object_fit), fitted(direct), tolerance = 2e-12)
  expect_equal(se(object_fit), se(direct), tolerance = 2e-12)
  expect_equal(direct_helper$mean, fitted(direct), tolerance = 2e-12)
  expect_equal(as.numeric(predicted$fit), fitted(direct), tolerance = 2e-12)
  expect_equal(as.numeric(predicted$se.fit), se(direct), tolerance = 2e-12)
  expect_equal(residuals(in_sample), training$y - fitted(in_sample),
               tolerance = 2e-12)
})

test_that("beta regression log-sum-exp survives complete raw-weight underflow", {
  training <- data.frame(x = 0.9 + c(0, 1e-7, 2e-7))
  response <- c(1, 4, 9)
  evaluation <- data.frame(x = 0)
  h <- 0.001
  tau <- (1 / h)^2
  log_weights <- dbeta(training$x, shape1 = 1, shape2 = 1 + tau,
                       log = TRUE)
  normalized <- exp(log_weights - max(log_weights))
  normalized <- normalized / sum(normalized)
  expected_mean <- sum(normalized * response)
  expected_variance <- sum(normalized * (response - expected_mean)^2)
  expected_se <- sqrt(expected_variance * sum(normalized^2))

  raw <- npksum(
    bws = h, txdat = training, exdat = evaluation,
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1,
    return.kernel.weights = TRUE
  )
  fit <- npreg(
    bws = h, txdat = training, tydat = response, exdat = evaluation,
    regtype = "lc", ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )

  expect_true(all(raw$kw == 0))
  expect_equal(fitted(fit), expected_mean, tolerance = 2e-12)
  expect_equal(se(fit), expected_se, tolerance = 2e-12)
})

test_that("unsupported beta regression surfaces fail explicitly", {
  training <- data.frame(x = c(0, 0.03, 0.2, 0.6, 1))
  response <- c(0, 1, 0.5, 2, 1.5)

  expect_error(
    suppressWarnings(npreg(
      bws = 0.1, txdat = training, tydat = response,
      ckertype = "beta", ckerbound = "none"
    )),
    "require ckerbound = \"fixed\" or \"range\"",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(npreg(
      bws = 0.1, txdat = training, tydat = response,
      regtype = "ll", ckertype = "beta",
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )),
    "only regtype = \"lc\"",
    fixed = TRUE
  )

  bw <- npregbw(
    xdat = training, ydat = response, bws = 0.1,
    bandwidth.compute = FALSE, regtype = "lc",
    ckertype = "beta", ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  expect_error(
    npreg(bws = bw, txdat = training, tydat = response, gradients = TRUE),
    "gradients are not yet available",
    fixed = TRUE
  )
  expect_error(
    npreg(bws = bw, txdat = training, tydat = response,
          exdat = data.frame(x = 1.01)),
    "Evaluation data violate",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(npreg(
      bws = c(0.1, 0.2),
      txdat = data.frame(x = training$x, group = factor(c(1, 1, 1, 2, 2))),
      tydat = response, regtype = "lc", ckertype = "beta",
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )),
    "continuous variables only",
    fixed = TRUE
  )
})
