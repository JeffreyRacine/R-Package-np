beta_overlap_order_coefficients <- function(order) {
  switch(as.character(order),
         `2` = 1,
         `4` = c(2, -1),
         `6` = c(3, -3, 1),
         `8` = c(4, -6, 4, -1))
}

beta_overlap_component <- function(center_one, bandwidth_one, scale_one,
                                   center_two, bandwidth_two, scale_two,
                                   lower, upper) {
  support_length <- upper - lower
  unit_one <- (center_one - lower) / support_length
  unit_two <- (center_two - lower) / support_length
  tau_one <- (support_length / bandwidth_one)^2 / scale_one
  tau_two <- (support_length / bandwidth_two)^2 / scale_two
  alpha_one <- 1 + unit_one * tau_one
  beta_one <- 1 + (1 - unit_one) * tau_one
  alpha_two <- 1 + unit_two * tau_two
  beta_two <- 1 + (1 - unit_two) * tau_two

  exp(-log(support_length) +
        lbeta(alpha_one + alpha_two - 1, beta_one + beta_two - 1) -
        lbeta(alpha_one, beta_one) - lbeta(alpha_two, beta_two))
}

beta_higher_overlap_expected <- function(center_one, bandwidth_one,
                                         center_two, bandwidth_two,
                                         lower, upper, order) {
  coefficients <- beta_overlap_order_coefficients(order)
  sum(vapply(seq_along(coefficients), function(scale_one) {
    sum(vapply(seq_along(coefficients), function(scale_two) {
      coefficients[scale_one] * coefficients[scale_two] *
        beta_overlap_component(
          center_one, bandwidth_one, scale_one,
          center_two, bandwidth_two, scale_two,
          lower, upper
        )
    }, numeric(1L)))
  }, numeric(1L)))
}

beta_higher_overlap_pdf <- function(observation, center, bandwidth,
                                    lower, upper, order) {
  support_length <- upper - lower
  unit_observation <- (observation - lower) / support_length
  unit_center <- (center - lower) / support_length
  tau <- (support_length / bandwidth)^2
  coefficients <- beta_overlap_order_coefficients(order)
  Reduce(`+`, Map(function(scale, coefficient) {
    coefficient * dbeta(
      unit_observation,
      1 + unit_center * tau / scale,
      1 + (1 - unit_center) * tau / scale
    ) / support_length
  }, seq_along(coefficients), coefficients))
}

beta_overlap_adaptive_radius <- function(train, k) {
  vapply(seq_along(train), function(index) {
    distance <- abs(train - train[index])
    duplicate_count <- sum(distance == 0) - 1L
    positive <- sort(distance[distance > 0])
    if (duplicate_count >= k) positive[1L] else
      positive[max(1L, k - duplicate_count)]
  }, numeric(1L))
}

beta_overlap_generalized_radius <- function(train, evaluation, k) {
  vapply(evaluation, function(target) {
    distance <- abs(train - target)
    exact_count <- sum(distance == 0)
    positive <- sort(distance[distance > 0])
    if (exact_count >= k) positive[1L] else
      positive[k - exact_count]
  }, numeric(1L))
}

test_that("higher-order beta convolution matches the analytic double sum", {
  training <- data.frame(x = c(-2, -1.8, -0.6, 1.1, 2.7, 3))
  evaluation <- data.frame(x = c(-2, -1.2, 0.4, 2.2, 3))
  bandwidth <- 0.72

  for (order in c(2L, 4L, 6L, 8L)) {
    fit <- npksum(
      bws = bandwidth, txdat = training, exdat = evaluation,
      operator = "convolution", return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = -2, ckerub = 3
    )
    expected <- vapply(evaluation$x, function(center) {
      vapply(training$x, function(observation_center) {
        beta_higher_overlap_expected(
          center, bandwidth, observation_center, bandwidth,
          -2, 3, order
        )
      }, numeric(1L))
    }, numeric(nrow(training)))

    expect_equal(fit$kw, expected, tolerance = 8e-11)
    expect_equal(as.double(fit$ksum), colSums(expected), tolerance = 8e-11)
  }
})

test_that("analytic higher-order overlap agrees with direct quadrature", {
  pairs <- list(c(0, 0), c(0.03, 0.42), c(0.5, 0.87), c(1, 1))
  bandwidth <- 0.17

  for (order in c(2L, 4L, 6L, 8L)) {
    for (pair in pairs) {
      analytic <- npksum(
        bws = bandwidth, txdat = data.frame(x = pair[2L]),
        exdat = data.frame(x = pair[1L]), operator = "convolution",
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1
      )$ksum
      quadrature <- integrate(
        function(observation) {
          beta_higher_overlap_pdf(observation, pair[1L], bandwidth,
                                  0, 1, order) *
            beta_higher_overlap_pdf(observation, pair[2L], bandwidth,
                                    0, 1, order)
        },
        lower = 0, upper = 1, rel.tol = 2e-10, subdivisions = 600L
      )$value

      expect_equal(as.double(analytic), quadrature, tolerance = 2e-8)
    }
  }
})

test_that("higher-order nearest-neighbor overlap uses both local bandwidths", {
  train <- c(0.04, 0.16, 0.34, 0.59, 0.81, 0.96)
  evaluation <- c(0.02, 0.24, 0.53, 0.88)
  k <- 2L
  eval_bandwidth <- beta_overlap_generalized_radius(train, evaluation, k)
  generalized_train_bandwidth <-
    beta_overlap_generalized_radius(train, train, k)
  adaptive_train_bandwidth <- beta_overlap_adaptive_radius(train, k)

  for (order in c(2L, 4L, 6L, 8L)) {
    generalized_expected <- vapply(seq_along(evaluation), function(index) {
      vapply(seq_along(train), function(train_index) {
        beta_higher_overlap_expected(
          evaluation[index], eval_bandwidth[index],
          train[train_index], generalized_train_bandwidth[train_index],
          0, 1, order
        )
      }, numeric(1L))
    }, numeric(length(train)))
    adaptive_expected <- vapply(seq_along(evaluation), function(index) {
      vapply(seq_along(train), function(train_index) {
        beta_higher_overlap_expected(
          evaluation[index], eval_bandwidth[index],
          train[train_index], adaptive_train_bandwidth[train_index],
          0, 1, order
        )
      }, numeric(1L))
    }, numeric(length(train)))
    generalized <- npksum(
      bws = k, txdat = data.frame(x = train),
      exdat = data.frame(x = evaluation), bwtype = "generalized_nn",
      operator = "convolution", return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )
    adaptive <- npksum(
      bws = k, txdat = data.frame(x = train),
      exdat = data.frame(x = evaluation), bwtype = "adaptive_nn",
      operator = "convolution", return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )

    expect_equal(generalized$kw, generalized_expected, tolerance = 2e-10)
    expect_equal(adaptive$kw, adaptive_expected, tolerance = 2e-10)
  }
})

test_that("higher-order overlap is symmetric with positive roughness diagonal", {
  centers <- data.frame(x = c(0, 0.04, 0.31, 0.66, 0.95, 1))

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      fit <- npksum(
        bws = if (identical(bwtype, "fixed")) 0.16 else 3,
        txdat = centers, bwtype = bwtype,
        operator = "convolution", return.kernel.weights = TRUE,
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1
      )

      expect_equal(fit$kw, t(fit$kw), tolerance = 3e-11)
      expect_true(all(diag(fit$kw) > 0))
    }
  }
})

test_that("higher-order overlap preserves inverse support units", {
  training <- data.frame(x = c(0.04, 0.2, 0.49, 0.78, 0.97))
  evaluation <- data.frame(x = c(0.02, 0.36, 0.71, 0.99))

  for (order in c(2L, 4L, 6L, 8L)) {
    unit <- npksum(
      bws = 0.15, txdat = training, exdat = evaluation,
      operator = "convolution", return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )
    shifted <- npksum(
      bws = 1.5, txdat = data.frame(x = -3 + 10 * training$x),
      exdat = data.frame(x = -3 + 10 * evaluation$x),
      operator = "convolution", return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = -3, ckerub = 7
    )

    expect_equal(unit$kw, 10 * shifted$kw, tolerance = 3e-10)
  }
})
