beta_conditional_side_weights <- function(train, evaluation, bandwidth,
                                          bwtype, kertype, order,
                                          lower, upper,
                                          operator = "normal") {
  args <- list(
    bws = bandwidth, txdat = data.frame(value = train),
    exdat = data.frame(value = evaluation), bwtype = bwtype,
    ckertype = kertype, ckerorder = order,
    operator = operator, return.kernel.weights = TRUE
  )
  if (identical(kertype, "beta")) {
    args$ckerbound <- "fixed"
    args$ckerlb <- lower
    args$ckerub <- upper
  }
  weights <- do.call(npksum, args)$kw
  if (!identical(kertype, "beta") && identical(operator, "normal") &&
      identical(bwtype, "fixed"))
    weights <- weights / bandwidth
  weights
}

beta_conditional_oracle <- function(x_weights, y_weights) {
  denominator <- colSums(x_weights)
  estimate <- colSums(x_weights * y_weights) / denominator
  influence <- x_weights * sweep(y_weights, 2L, estimate, "-")
  stderr <- sqrt(colSums(influence^2) / (nrow(x_weights) - 1L)) /
    abs(denominator)
  list(estimate = estimate, stderr = stderr)
}

test_that("conditional beta X and Y kernels match side-weight ratios", {
  training_x <- data.frame(x = c(0.02, 0.08, 0.2, 0.39, 0.61, 0.8, 0.94, 0.99))
  training_y <- data.frame(y = c(0.03, 0.11, 0.18, 0.42, 0.55, 0.76, 0.88, 0.97))
  evaluation_x <- data.frame(x = c(0.01, 0.15, 0.47, 0.83, 0.98))
  evaluation_y <- data.frame(y = c(0.02, 0.21, 0.5, 0.79, 0.99))

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bandwidth_y <- if (identical(bwtype, "fixed")) 0.16 else 3
    bandwidth_x <- if (identical(bwtype, "fixed")) 0.14 else 3
    common <- list(
      xdat = training_x, ydat = training_y,
      bws = c(bandwidth_y, bandwidth_x),
      bandwidth.compute = FALSE, bwtype = bwtype,
      cxkertype = "beta", cxkerorder = 4,
      cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
      cykertype = "beta", cykerorder = 6,
      cykerbound = "fixed", cykerlb = 0, cykerub = 1
    )
    density_bw <- do.call(npcdensbw, common)
    distribution_bw <- do.call(npcdistbw, common)
    density <- npcdens(
      bws = density_bw, txdat = training_x, tydat = training_y,
      exdat = evaluation_x, eydat = evaluation_y
    )
    distribution <- npcdist(
      bws = distribution_bw, txdat = training_x, tydat = training_y,
      exdat = evaluation_x, eydat = evaluation_y
    )
    x_weights <- beta_conditional_side_weights(
      training_x$x, evaluation_x$x, bandwidth_x, bwtype,
      "beta", 4, 0, 1
    )
    y_pdf <- beta_conditional_side_weights(
      training_y$y, evaluation_y$y, bandwidth_y, bwtype,
      "beta", 6, 0, 1
    )
    y_cdf <- beta_conditional_side_weights(
      training_y$y, evaluation_y$y, bandwidth_y, bwtype,
      "beta", 6, 0, 1, operator = "integral"
    )
    density_expected <- beta_conditional_oracle(x_weights, y_pdf)
    distribution_expected <- beta_conditional_oracle(x_weights, y_cdf)

    expect_equal(fitted(density), density_expected$estimate, tolerance = 3e-10)
    expect_equal(se(density), density_expected$stderr, tolerance = 3e-10)
    expect_equal(fitted(distribution), distribution_expected$estimate,
                 tolerance = 3e-10)
    expect_equal(se(distribution), distribution_expected$stderr,
                 tolerance = 3e-10)
  }
})

test_that("beta and legacy conditional kernels can be mixed by side", {
  training_x <- data.frame(x = c(0.02, 0.08, 0.2, 0.39, 0.61, 0.8, 0.94, 0.99))
  training_y <- data.frame(y = c(-1.2, -0.8, -0.3, 0.1, 0.35, 0.7, 1.0, 1.3))
  evaluation_x <- data.frame(x = c(0.01, 0.15, 0.47, 0.83, 0.98))
  evaluation_y <- data.frame(y = c(-1.1, -0.4, 0.2, 0.8, 1.2))
  combinations <- list(
    list(xkernel = "beta", xorder = 8L, ykernel = "gaussian", yorder = 4L),
    list(xkernel = "gaussian", xorder = 4L, ykernel = "beta", yorder = 8L)
  )

  for (specification in combinations) {
    y_train <- if (identical(specification$ykernel, "beta")) {
      data.frame(y = (training_y$y + 1.5) / 3)
    } else training_y
    y_eval <- if (identical(specification$ykernel, "beta")) {
      data.frame(y = (evaluation_y$y + 1.5) / 3)
    } else evaluation_y
    x_train <- if (identical(specification$xkernel, "beta")) training_x else
      data.frame(x = qnorm(pmin(pmax(training_x$x, 0.01), 0.99)))
    x_eval <- if (identical(specification$xkernel, "beta")) evaluation_x else
      data.frame(x = qnorm(pmin(pmax(evaluation_x$x, 0.01), 0.99)))
    x_bandwidth <- 0.16
    y_bandwidth <- if (identical(specification$ykernel, "beta")) 0.15 else 0.35
    common <- list(
      xdat = x_train, ydat = y_train,
      bws = c(y_bandwidth, x_bandwidth),
      bandwidth.compute = FALSE,
      cxkertype = specification$xkernel,
      cxkerorder = specification$xorder,
      cykertype = specification$ykernel,
      cykerorder = specification$yorder
    )
    if (identical(specification$xkernel, "beta")) {
      common$cxkerbound <- "fixed"
      common$cxkerlb <- 0
      common$cxkerub <- 1
    }
    if (identical(specification$ykernel, "beta")) {
      common$cykerbound <- "fixed"
      common$cykerlb <- 0
      common$cykerub <- 1
    }
    density_bw <- do.call(npcdensbw, common)
    distribution_bw <- do.call(npcdistbw, common)
    density <- npcdens(
      bws = density_bw, txdat = x_train, tydat = y_train,
      exdat = x_eval, eydat = y_eval
    )
    distribution <- npcdist(
      bws = distribution_bw, txdat = x_train, tydat = y_train,
      exdat = x_eval, eydat = y_eval
    )
    x_weights <- beta_conditional_side_weights(
      x_train$x, x_eval$x, x_bandwidth, "fixed",
      specification$xkernel, specification$xorder, 0, 1
    )
    y_pdf <- beta_conditional_side_weights(
      y_train$y, y_eval$y, y_bandwidth, "fixed",
      specification$ykernel, specification$yorder, 0, 1
    )
    y_cdf <- beta_conditional_side_weights(
      y_train$y, y_eval$y, y_bandwidth, "fixed",
      specification$ykernel, specification$yorder, 0, 1,
      operator = "integral"
    )
    density_expected <- beta_conditional_oracle(x_weights, y_pdf)
    distribution_expected <- beta_conditional_oracle(x_weights, y_cdf)

    expect_equal(fitted(density), density_expected$estimate, tolerance = 4e-9)
    expect_equal(se(density), density_expected$stderr, tolerance = 4e-9)
    expect_equal(fitted(distribution), distribution_expected$estimate,
                 tolerance = 4e-9)
    expect_equal(se(distribution), distribution_expected$stderr,
                 tolerance = 4e-9)
  }
})

test_that("beta conditional distribution is exact at dependent support endpoints", {
  training_x <- data.frame(x = c(0.03, 0.12, 0.28, 0.53, 0.77, 0.96))
  training_y <- data.frame(y = c(0.02, 0.1, 0.31, 0.59, 0.82, 0.98))
  evaluation_x <- data.frame(x = rep(0.45, 2L))
  evaluation_y <- data.frame(y = c(0, 1))
  bw <- npcdistbw(
    xdat = training_x, ydat = training_y, bws = c(0.15, 0.17),
    bandwidth.compute = FALSE,
    cxkertype = "beta", cxkerorder = 4,
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykertype = "beta", cykerorder = 8,
    cykerbound = "fixed", cykerlb = 0, cykerub = 1
  )
  fit <- npcdist(
    bws = bw, txdat = training_x, tydat = training_y,
    exdat = evaluation_x, eydat = evaluation_y
  )

  expect_identical(fitted(fit), c(0, 1))
})

test_that("conditional beta supports formula objects and prediction", {
  training <- data.frame(
    y = c(0.03, 0.11, 0.18, 0.42, 0.55, 0.76, 0.88, 0.97),
    x = c(0.02, 0.08, 0.2, 0.39, 0.61, 0.8, 0.94, 0.99)
  )
  evaluation <- data.frame(
    y = c(0.02, 0.21, 0.5, 0.79, 0.99),
    x = c(0.01, 0.15, 0.47, 0.83, 0.98)
  )
  common <- list(
    formula = y ~ x, data = training,
    bws = c(0.16, 0.14), bandwidth.compute = FALSE,
    cxkertype = "beta", cxkerorder = 4,
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykertype = "beta", cykerorder = 8,
    cykerbound = "fixed", cykerlb = 0, cykerub = 1
  )
  density_bw <- do.call(npcdensbw, common)
  distribution_bw <- do.call(npcdistbw, common)
  density_formula <- npcdens(bws = density_bw, data = training,
                             newdata = evaluation)
  distribution_formula <- npcdist(bws = distribution_bw, data = training,
                                   newdata = evaluation)
  density_native <- npcdens(
    bws = density_bw,
    txdat = training["x"], tydat = training["y"],
    exdat = evaluation["x"], eydat = evaluation["y"]
  )
  distribution_native <- npcdist(
    bws = distribution_bw,
    txdat = training["x"], tydat = training["y"],
    exdat = evaluation["x"], eydat = evaluation["y"]
  )

  expect_equal(fitted(density_formula), fitted(density_native),
               tolerance = 3e-10)
  expect_equal(se(density_formula), se(density_native), tolerance = 3e-10)
  expect_equal(fitted(distribution_formula), fitted(distribution_native),
               tolerance = 3e-10)
  expect_equal(se(distribution_formula), se(distribution_native),
               tolerance = 3e-10)
  expect_equal(predict(density_formula, newdata = evaluation),
               fitted(density_native), tolerance = 3e-10)
  expect_equal(predict(distribution_formula, newdata = evaluation),
               fitted(distribution_native), tolerance = 3e-10)
  density_se <- predict(density_formula, newdata = evaluation, se.fit = TRUE)
  distribution_se <- predict(distribution_formula, newdata = evaluation,
                             se.fit = TRUE)
  expect_equal(density_se$fit, fitted(density_native), tolerance = 3e-10)
  expect_equal(density_se$se.fit, se(density_native), tolerance = 3e-10)
  expect_equal(distribution_se$fit, fitted(distribution_native),
               tolerance = 3e-10)
  expect_equal(distribution_se$se.fit, se(distribution_native),
               tolerance = 3e-10)
})

test_that("conditional beta guards unsupported automatic and gradient routes", {
  training_x <- data.frame(x = c(0.02, 0.12, 0.34, 0.68, 0.94))
  training_y <- data.frame(y = c(0.03, 0.2, 0.4, 0.75, 0.97))

  expect_error(
    suppressWarnings(npcdensbw(
      xdat = training_x, ydat = training_y,
      cxkertype = "beta", cxkerbound = "fixed",
      cxkerlb = 0, cxkerub = 1
    )),
    "does not yet support automatic bandwidth selection"
  )
  bw <- npcdensbw(
    xdat = training_x, ydat = training_y, bws = c(0.15, 0.15),
    bandwidth.compute = FALSE,
    cxkertype = "beta", cxkerbound = "fixed",
    cxkerlb = 0, cxkerub = 1
  )
  expect_error(
    npcdens(bws = bw, txdat = training_x, tydat = training_y,
            gradients = TRUE),
    "gradients are not yet available",
    fixed = TRUE
  )
})
