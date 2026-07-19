beta_order_coefficients <- function(order) {
  switch(as.character(order),
         `2` = 1,
         `4` = c(2, -1),
         `6` = c(3, -3, 1),
         `8` = c(4, -6, 4, -1))
}

beta_higher_pdf_expected <- function(target, observation, bandwidth,
                                     lower, upper, order) {
  support_length <- upper - lower
  target_unit <- (target - lower) / support_length
  observation_unit <- (observation - lower) / support_length
  tau <- (support_length / bandwidth)^2
  coefficients <- beta_order_coefficients(order)
  Reduce(`+`, Map(function(scale, coefficient) {
    coefficient * dbeta(
      observation_unit,
      1 + target_unit * tau / scale,
      1 + (1 - target_unit) * tau / scale
    ) / support_length
  }, seq_along(coefficients), coefficients))
}

beta_higher_cdf_expected <- function(target, observation, bandwidth,
                                     lower, upper, order) {
  support_length <- upper - lower
  target_unit <- (target - lower) / support_length
  observation_unit <- (observation - lower) / support_length
  tau <- (support_length / bandwidth)^2
  coefficients <- beta_order_coefficients(order)
  Reduce(`+`, Map(function(scale, coefficient) {
    coefficient * pbeta(
      target_unit,
      1 + observation_unit * tau / scale,
      1 + (1 - observation_unit) * tau / scale
    )
  }, seq_along(coefficients), coefficients))
}

test_that("higher-order beta PDFs match signed multiscale beta combinations", {
  training <- data.frame(x = c(-3, -2.7, -0.4, 2.2, 6.1, 7))
  evaluation <- data.frame(x = c(-3, -2.999999, -1.8, 1.7, 5.4, 6.999999, 7))
  bandwidth <- 1.25

  for (order in c(2L, 4L, 6L, 8L)) {
    fit <- npksum(
      bws = bandwidth, txdat = training, exdat = evaluation,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = -3, ckerub = 7,
      return.kernel.weights = TRUE
    )
    expected <- vapply(evaluation$x, function(target) {
      beta_higher_pdf_expected(target, training$x, bandwidth,
                               -3, 7, order)
    }, numeric(nrow(training)))

    expect_equal(fit$kw, expected, tolerance = 3e-12)
    expect_equal(as.double(fit$ksum), colSums(expected), tolerance = 3e-12)
  }
})

test_that("higher-order beta CDFs match signed multiscale beta combinations", {
  training <- data.frame(x = c(0, 0.03, 0.24, 0.57, 0.91, 1))
  evaluation <- data.frame(x = c(0, 1e-10, 0.07, 0.41, 0.82, 1 - 1e-10, 1))
  bandwidth <- 0.13

  for (order in c(2L, 4L, 6L, 8L)) {
    fit <- npksum(
      bws = bandwidth, txdat = training, exdat = evaluation,
      operator = "integral", return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )
    expected <- vapply(evaluation$x, function(target) {
      beta_higher_cdf_expected(target, training$x, bandwidth,
                               0, 1, order)
    }, numeric(nrow(training)))

    expect_equal(fit$kw, expected, tolerance = 4e-12)
    expect_equal(as.double(fit$ksum), colSums(expected), tolerance = 4e-12)
    expect_identical(fit$kw[, 1L], rep(0, nrow(training)))
    expect_identical(fit$kw[, nrow(evaluation)], rep(1, nrow(training)))
  }
})

test_that("each higher-order beta PDF has unit mass on its declared support", {
  bandwidth <- 0.16
  centers <- c(0, 0.04, 0.5, 0.93, 1)

  for (order in c(2L, 4L, 6L, 8L)) {
    for (center in centers) {
      integrated <- integrate(
        function(observation) {
          beta_higher_pdf_expected(center, observation, bandwidth,
                                   0, 1, order)
        },
        lower = 0, upper = 1, rel.tol = 2e-11, subdivisions = 400L
      )$value
      expect_equal(integrated, 1, tolerance = 3e-9)
    }
  }
})

test_that("signed higher-order beta values are returned without clipping", {
  grid <- data.frame(x = seq(0, 1, length.out = 2001L))

  for (order in c(4L, 6L, 8L)) {
    pdf <- npksum(
      bws = 0.12, txdat = grid, exdat = data.frame(x = 0.5),
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      return.kernel.weights = TRUE
    )$kw[, 1L]
    cdf <- npksum(
      bws = 0.12, txdat = data.frame(x = 0.08), exdat = grid,
      operator = "integral", ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      return.kernel.weights = TRUE
    )$kw[1L, ]

    expect_lt(min(pdf), 0)
    expect_true(min(cdf) < 0 || max(cdf) > 1)
    expect_identical(cdf[c(1L, length(cdf))], c(0, 1))
  }
})

test_that("higher-order beta kernels retain affine units and uniform limits", {
  observations <- data.frame(x = c(0.03, 0.27, 0.72, 0.98))
  targets <- data.frame(x = c(0.02, 0.31, 0.65, 0.97))

  for (order in c(2L, 4L, 6L, 8L)) {
    unit_pdf <- npksum(
      bws = 0.14, txdat = observations, exdat = targets,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      return.kernel.weights = TRUE
    )$kw
    shifted_pdf <- npksum(
      bws = 1.4, txdat = data.frame(x = -3 + 10 * observations$x),
      exdat = data.frame(x = -3 + 10 * targets$x),
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = -3, ckerub = 7,
      return.kernel.weights = TRUE
    )$kw
    uniform_pdf <- npksum(
      bws = 1e200, txdat = observations, exdat = targets,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      return.kernel.weights = TRUE
    )$kw
    uniform_cdf <- npksum(
      bws = 1e200, txdat = observations, exdat = targets,
      operator = "integral", ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      return.kernel.weights = TRUE
    )$kw

    expect_equal(unit_pdf, 10 * shifted_pdf, tolerance = 5e-12)
    expect_equal(uniform_pdf, matrix(1, nrow(observations), nrow(targets)),
                 tolerance = 3e-14)
    expect_equal(uniform_cdf,
                 matrix(rep(targets$x, each = nrow(observations)),
                        nrow = nrow(observations)),
                 tolerance = 3e-14)
  }
})
