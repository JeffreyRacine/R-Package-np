beta_regression_order_coefficients <- function(order) {
  switch(as.character(order),
         `2` = 1,
         `4` = c(2, -1),
         `6` = c(3, -3, 1),
         `8` = c(4, -6, 4, -1))
}

test_that("higher-order beta local-constant regression matches signed weights", {
  training <- data.frame(
    x1 = c(0.01, 0.05, 0.14, 0.31, 0.55, 0.73, 0.9, 0.99),
    x2 = c(-2, -1.8, -1.1, 0.2, 1.8, 3.1, 4.5, 5)
  )
  response <- sin(2 * training$x1) + 0.15 * training$x2^2
  evaluation <- data.frame(
    x1 = c(0.03, 0.12, 0.4, 0.78, 0.97),
    x2 = c(-1.9, -1.4, 1, 4, 4.9)
  )

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      common <- list(
        bws = if (identical(bwtype, "fixed")) c(0.16, 1.15) else c(3, 3),
        bwtype = bwtype, ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = c(0, -2), ckerub = c(1, 5)
      )
      fit <- do.call(npreg, c(list(
        txdat = training, tydat = response, exdat = evaluation,
        regtype = "lc"
      ), common))
      sums <- do.call(npksum, c(list(
        txdat = training, exdat = evaluation,
        return.kernel.weights = TRUE
      ), common))
      denominator <- colSums(sums$kw)
      expect_true(all(abs(denominator) > 1e-8))
      normalized <- sweep(sums$kw, 2L, denominator, "/")
      expected_mean <- colSums(normalized * response)
      centered <- response - matrix(expected_mean,
                                    nrow = nrow(training),
                                    ncol = nrow(evaluation),
                                    byrow = TRUE)
      expected_variance <- pmax(colSums(normalized * centered^2), 0)
      expected_se <- sqrt(expected_variance * colSums(normalized^2))

      expect_equal(fitted(fit), expected_mean, tolerance = 2e-10)
      expect_equal(se(fit), expected_se, tolerance = 2e-10)
      expect_true(all(is.finite(fitted(fit))))
      expect_true(all(is.finite(se(fit))))
    }
  }
})

test_that("higher-order beta regression preserves formula, object, and prediction routes", {
  training <- data.frame(x = c(0.01, 0.06, 0.17, 0.35, 0.58, 0.79, 0.93, 0.99))
  training$y <- cos(3 * training$x) + training$x
  evaluation <- data.frame(x = c(0.02, 0.11, 0.48, 0.87, 0.98))

  for (order in c(4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      common <- list(
        bws = if (identical(bwtype, "fixed")) 0.16 else 3,
        bwtype = bwtype, regtype = "lc",
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1
      )
      direct <- do.call(npreg, c(list(
        txdat = training["x"], tydat = training$y, exdat = evaluation
      ), common))
      formula <- do.call(npreg, c(list(
        y ~ x, data = training, newdata = evaluation
      ), common))
      bw <- do.call(npregbw, c(list(
        xdat = training["x"], ydat = training$y,
        bandwidth.compute = FALSE
      ), common))
      object_fit <- npreg(
        bws = bw, txdat = training["x"], tydat = training$y,
        exdat = evaluation
      )
      prediction <- predict(formula, newdata = evaluation, se.fit = TRUE)

      expect_identical(bw$ckerorder, order)
      expect_identical(bw$type, bwtype)
      expect_equal(fitted(formula), fitted(direct), tolerance = 2e-10)
      expect_equal(fitted(object_fit), fitted(direct), tolerance = 2e-10)
      expect_equal(se(object_fit), se(direct), tolerance = 2e-10)
      expect_equal(as.numeric(prediction$fit), fitted(direct), tolerance = 2e-10)
      expect_equal(as.numeric(prediction$se.fit), se(direct), tolerance = 2e-10)
    }
  }
})

test_that("higher-order beta regression survives complete raw-weight underflow", {
  training <- data.frame(x = 0.9 + c(0, 1e-7, 2e-7))
  response <- c(1, 4, 9)
  evaluation <- data.frame(x = 0)
  bandwidth <- 0.001

  for (order in c(4L, 6L, 8L)) {
    coefficients <- beta_regression_order_coefficients(order)
    log_absolute_weights <- vapply(training$x, function(observation) {
      component_logs <- vapply(seq_along(coefficients), function(scale) {
        dbeta(observation,
              1,
              1 + (1 / bandwidth)^2 / scale,
              log = TRUE) + log(abs(coefficients[scale]))
      }, numeric(1L))
      maximum <- max(component_logs)
      signed_scaled <- sum(sign(coefficients) *
                             exp(component_logs - maximum))
      maximum + log(abs(signed_scaled))
    }, numeric(1L))
    normalized <- exp(log_absolute_weights - max(log_absolute_weights))
    normalized <- normalized / sum(normalized)
    expected_mean <- sum(normalized * response)
    expected_variance <- sum(normalized * (response - expected_mean)^2)
    expected_se <- sqrt(expected_variance * sum(normalized^2))
    raw <- npksum(
      bws = bandwidth, txdat = training, exdat = evaluation,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      return.kernel.weights = TRUE
    )
    fit <- npreg(
      bws = bandwidth, txdat = training, tydat = response,
      exdat = evaluation, regtype = "lc",
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )

    expect_true(all(raw$kw == 0))
    expect_equal(fitted(fit), expected_mean, tolerance = 3e-10)
    expect_equal(se(fit), expected_se, tolerance = 3e-10)
  }
})

test_that("higher-order beta regression retains the uniform large-bandwidth limit", {
  training <- data.frame(x = c(0.1, 0.9))
  response <- c(1, 2)

  for (order in c(4L, 6L, 8L)) {
    fit <- npreg(
      bws = 1e200, txdat = training, tydat = response,
      exdat = data.frame(x = 0.5), regtype = "lc",
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )
    expect_equal(fitted(fit), mean(response), tolerance = 2e-14)
  }
})
