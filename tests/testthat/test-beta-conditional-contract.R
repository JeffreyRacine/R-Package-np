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
  if (!identical(kertype, "beta") && identical(operator, "normal") &&
      identical(bwtype, "generalized_nn")) {
    realized <- vapply(evaluation, function(value) {
      sort(abs(train - value), method = "quick")[[bandwidth]]
    }, numeric(1L))
    weights <- sweep(weights, 2L, realized, "/")
  }
  if (!identical(kertype, "beta") && identical(operator, "normal") &&
      identical(bwtype, "adaptive_nn")) {
    realized <- vapply(train, function(value) {
      sort(abs(train - value), method = "quick")[[bandwidth + 1L]]
    }, numeric(1L))
    weights <- sweep(weights, 1L, realized, "/")
  }
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

beta_conditional_lp_oracle <- function(training_x, evaluation_x, degree,
                                       bernstein, x_weights, y_weights) {
  basis <- np:::W.lp(
    training_x, degree = degree, basis = "glp",
    Bernstein = bernstein
  )
  evaluation_basis <- np:::W.lp(
    training_x, exdat = evaluation_x, degree = degree, basis = "glp",
    Bernstein = bernstein
  )
  evaluation_derivative <- np:::W.lp(
    training_x, exdat = evaluation_x, degree = degree,
    gradient.vec = rep(1L, length(degree)), basis = "glp",
    Bernstein = bernstein
  )
  estimate <- stderr <- gradient <- gradient_stderr <-
    numeric(nrow(evaluation_x))
  for (evaluation in seq_len(nrow(evaluation_x))) {
    weight <- x_weights[, evaluation]
    response <- y_weights[, evaluation]
    gram <- crossprod(basis, weight * basis)
    coefficient <- solve(gram, crossprod(basis, weight * response))
    variance <- sum(weight * response * response) / sum(weight) -
      (sum(weight * response) / sum(weight))^2
    gram_power_two <- crossprod(basis, (weight * weight) * basis)
    fit_vector <- solve(gram, evaluation_basis[evaluation, ])
    derivative_vector <- solve(gram, evaluation_derivative[evaluation, ])
    estimate[[evaluation]] <-
      drop(crossprod(evaluation_basis[evaluation, ], coefficient))
    stderr[[evaluation]] <- sqrt(max(
      0, variance * drop(crossprod(fit_vector,
                                  gram_power_two %*% fit_vector))
    ))
    gradient[[evaluation]] <-
      drop(crossprod(evaluation_derivative[evaluation, ], coefficient))
    gradient_stderr[[evaluation]] <- sqrt(max(
      0, variance * drop(crossprod(derivative_vector,
                                  gram_power_two %*% derivative_vector))
    ))
  }
  list(
    estimate = estimate, stderr = stderr,
    gradient = gradient, gradient_stderr = gradient_stderr
  )
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

test_that("conditional beta supports local-constant gradient routes", {
  training_x <- data.frame(x = c(0.02, 0.12, 0.34, 0.68, 0.94))
  training_y <- data.frame(y = c(0.03, 0.2, 0.4, 0.75, 0.97))

  bw <- npcdensbw(
    xdat = training_x, ydat = training_y, bws = c(0.15, 0.15),
    bandwidth.compute = FALSE,
    cxkertype = "beta", cxkerbound = "fixed",
    cxkerlb = 0, cxkerub = 1
  )
  fit <- npcdens(bws = bw, txdat = training_x, tydat = training_y,
                 gradients = TRUE)
  expect_equal(dim(gradients(fit)), c(nrow(training_x), 1L))
  expect_true(all(is.finite(gradients(fit))))
})

test_that("beta X uses the canonical conditional LP engine", {
  training_x <- data.frame(
    x = c(0.025, 0.07, 0.13, 0.21, 0.34, 0.46, 0.58, 0.67,
          0.76, 0.84, 0.91, 0.975)
  )
  training_y <- data.frame(
    y = c(-1.1, -0.82, -0.55, -0.31, -0.08, 0.13,
          0.29, 0.48, 0.66, 0.81, 1.02, 1.21)
  )
  evaluation_x <- data.frame(x = c(0.09, 0.27, 0.49, 0.72, 0.9))
  evaluation_y <- data.frame(y = c(-0.75, -0.25, 0.2, 0.6, 0.95))
  x_weights <- beta_conditional_side_weights(
    training_x$x, evaluation_x$x, 0.17, "fixed", "beta", 8L, 0, 1
  )
  y_weights <- beta_conditional_side_weights(
    training_y$y, evaluation_y$y, 0.3, "fixed", "gaussian", 4L, 0, 1
  )

  for (specification in list(
    list(degree = 2L, bernstein = FALSE),
    list(degree = 3L, bernstein = TRUE)
  )) {
    bw <- npcdensbw(
      xdat = training_x, ydat = training_y, bws = c(0.3, 0.17),
      bandwidth.compute = FALSE, regtype = "lp",
      degree = specification$degree, basis = "glp",
      bernstein.basis = specification$bernstein,
      cxkertype = "beta", cxkerorder = 8L,
      cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
      cykertype = "gaussian", cykerorder = 4L
    )
    fit <- npcdens(
      bws = bw, txdat = training_x, tydat = training_y,
      exdat = evaluation_x, eydat = evaluation_y, gradients = TRUE
    )
    expected <- beta_conditional_lp_oracle(
      training_x, evaluation_x, specification$degree,
      specification$bernstein, x_weights, y_weights
    )
    expect_equal(fitted(fit), expected$estimate, tolerance = 2e-8)
    expect_equal(se(fit), expected$stderr, tolerance = 2e-8)
    expect_equal(drop(fit$congrad), expected$gradient, tolerance = 3e-7)
    expect_equal(drop(fit$congerr), expected$gradient_stderr,
                 tolerance = 3e-7)
  }

  ll_bw <- npcdensbw(
    xdat = training_x, ydat = training_y, bws = c(0.3, 0.17),
    bandwidth.compute = FALSE, regtype = "ll",
    cxkertype = "beta", cxkerorder = 8L,
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykertype = "gaussian", cykerorder = 4L
  )
  lp1_bw <- npcdensbw(
    xdat = training_x, ydat = training_y, bws = c(0.3, 0.17),
    bandwidth.compute = FALSE, regtype = "lp", degree = 1L,
    basis = "glp", bernstein.basis = FALSE,
    cxkertype = "beta", cxkerorder = 8L,
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykertype = "gaussian", cykerorder = 4L
  )
  ll_fit <- npcdens(
    bws = ll_bw, txdat = training_x, tydat = training_y,
    exdat = evaluation_x, eydat = evaluation_y, gradients = TRUE
  )
  lp1_fit <- npcdens(
    bws = lp1_bw, txdat = training_x, tydat = training_y,
    exdat = evaluation_x, eydat = evaluation_y, gradients = TRUE
  )
  expect_identical(fitted(ll_fit), fitted(lp1_fit))
  expect_identical(se(ll_fit), se(lp1_fit))
  expect_identical(ll_fit$congrad, lp1_fit$congrad)
  expect_identical(ll_fit$congerr, lp1_fit$congerr)
})

test_that("beta X prepared NN bandwidths preserve scalar and LP algebra", {
  training_x <- data.frame(
    x = c(0.025, 0.07, 0.13, 0.21, 0.34, 0.46, 0.58, 0.67,
          0.76, 0.84, 0.91, 0.975)
  )
  training_y <- data.frame(
    y = c(-1.1, -0.82, -0.55, -0.31, -0.08, 0.13,
          0.29, 0.48, 0.66, 0.81, 1.02, 1.21)
  )
  evaluation_x <- data.frame(x = c(0.09, 0.27, 0.49, 0.72, 0.9))
  evaluation_y <- data.frame(y = c(-0.75, -0.25, 0.2, 0.6, 0.95))
  neighbours <- 4L

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    x_weights <- beta_conditional_side_weights(
      training_x$x, evaluation_x$x, neighbours, bwtype,
      "beta", 8L, 0, 1
    )
    for (kind in c("density", "distribution")) {
      y_operator <- if (identical(kind, "density")) "normal" else "integral"
      y_weights <- beta_conditional_side_weights(
        training_y$y, evaluation_y$y, neighbours, bwtype,
        "gaussian", 4L, 0, 1, operator = y_operator
      )
      constructor <- if (identical(kind, "density")) npcdensbw else npcdistbw
      estimator <- if (identical(kind, "density")) npcdens else npcdist
      common <- list(
        xdat = training_x, ydat = training_y,
        bws = c(neighbours, neighbours), bandwidth.compute = FALSE,
        bwtype = bwtype,
        cxkertype = "beta", cxkerorder = 8L,
        cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
        cykertype = "gaussian", cykerorder = 4L
      )

      scalar_bw <- do.call(constructor, c(common, list(regtype = "lc")))
      scalar <- estimator(
        bws = scalar_bw, txdat = training_x, tydat = training_y,
        exdat = evaluation_x, eydat = evaluation_y
      )
      scalar_expected <- beta_conditional_oracle(x_weights, y_weights)
      expect_equal(fitted(scalar), scalar_expected$estimate,
                   tolerance = 2e-8)
      expect_equal(se(scalar), scalar_expected$stderr, tolerance = 2e-8)

      lp_bw <- do.call(constructor, c(
        common,
        list(regtype = "lp", degree = 2L, basis = "glp",
             bernstein.basis = FALSE)
      ))
      lp <- estimator(
        bws = lp_bw, txdat = training_x, tydat = training_y,
        exdat = evaluation_x, eydat = evaluation_y, gradients = TRUE
      )
      lp_expected <- beta_conditional_lp_oracle(
        training_x, evaluation_x, 2L, FALSE, x_weights, y_weights
      )
      expect_equal(fitted(lp), lp_expected$estimate, tolerance = 2e-8)
      expect_equal(se(lp), lp_expected$stderr, tolerance = 2e-8)
      expect_equal(drop(lp$congrad), lp_expected$gradient, tolerance = 3e-7)
      expect_equal(drop(lp$congerr), lp_expected$gradient_stderr,
                   tolerance = 3e-7)
    }
  }
})

test_that("beta-X mixed conditional fits honor categorical compression", {
  training_x <- data.frame(
    x = seq(0.035, 0.965, length.out = 24L),
    group = factor(rep(c("a", "b", "c"), length.out = 24L))
  )
  training_y <- data.frame(
    y = sin(2 * pi * training_x$x) +
      as.integer(training_x$group) / 9
  )
  evaluation_x <- data.frame(
    x = seq(0.09, 0.91, length.out = 8L),
    group = factor(rep(c("a", "c"), length.out = 8L),
                   levels = levels(training_x$group))
  )
  evaluation_y <- data.frame(y = seq(-0.6, 0.7, length.out = 8L))
  results <- list()

  for (compress in c(FALSE, TRUE)) {
    old_options <- options(np.categorical.compress = compress)
    on.exit(options(old_options), add = TRUE)
    for (regtype in c("lc", "lp")) {
      arguments <- list(
        xdat = training_x, ydat = training_y,
        bws = c(0.28, 0.17, 0.22), bandwidth.compute = FALSE,
        regtype = regtype,
        cxkertype = "beta", cxkerorder = 6L,
        cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
        cykertype = "gaussian", cykerorder = 4L
      )
      if (identical(regtype, "lp")) {
        arguments$degree <- 2L
        arguments$basis <- "glp"
        arguments$bernstein.basis <- FALSE
      }
      bw <- do.call(npcdensbw, arguments)
      results[[paste(regtype, compress)]] <- npcdens(
        bws = bw, txdat = training_x, tydat = training_y,
        exdat = evaluation_x, eydat = evaluation_y, gradients = TRUE
      )
    }
    options(old_options)
  }

  for (regtype in c("lc", "lp")) {
    dense <- results[[paste(regtype, FALSE)]]
    compressed <- results[[paste(regtype, TRUE)]]
    expect_identical(fitted(dense), fitted(compressed))
    expect_identical(se(dense), se(compressed))
    expect_identical(dense$congrad, compressed$congrad)
    expect_identical(dense$congerr, compressed$congerr)
  }
})
