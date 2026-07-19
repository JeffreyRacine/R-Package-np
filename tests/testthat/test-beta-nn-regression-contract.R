test_that("nearest-neighbor beta local-constant regression matches kernel weights", {
  training <- data.frame(
    x1 = c(0.01, 0.05, 0.14, 0.31, 0.55, 0.73, 0.9, 0.99),
    x2 = c(-2, -1.8, -1.1, 0.2, 1.8, 3.1, 4.5, 5)
  )
  response <- sin(2 * training$x1) + 0.15 * training$x2^2
  evaluation <- data.frame(
    x1 = c(0, 0.08, 0.4, 0.81, 1),
    x2 = c(-2, -1.4, 1, 4, 5)
  )

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    common <- list(
      bws = c(3, 3), bwtype = bwtype,
      ckertype = "beta", ckerorder = 2,
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
    normalized <- sweep(sums$kw, 2L, colSums(sums$kw), "/")
    expected_mean <- colSums(normalized * response)
    centered <- response - matrix(expected_mean,
                                  nrow = nrow(training),
                                  ncol = nrow(evaluation),
                                  byrow = TRUE)
    expected_variance <- colSums(normalized * centered^2)
    expected_se <- sqrt(expected_variance * colSums(normalized^2))

    expect_equal(fitted(fit), expected_mean, tolerance = 3e-12)
    expect_equal(se(fit), expected_se, tolerance = 3e-12)
    expect_true(all(is.finite(fitted(fit))))
    expect_true(all(is.finite(se(fit))))
  }
})

test_that("nearest-neighbor beta regression preserves object and prediction routes", {
  training <- data.frame(x = c(0.01, 0.06, 0.17, 0.35, 0.58, 0.79, 0.93, 0.99))
  training$y <- cos(3 * training$x) + training$x
  evaluation <- data.frame(x = c(0, 0.11, 0.48, 0.87, 1))

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    common <- list(
      bws = 3, bwtype = bwtype, regtype = "lc",
      ckertype = "beta", ckerorder = 2,
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

    expect_s3_class(bw, "rbandwidth")
    expect_identical(bw$type, bwtype)
    expect_equal(fitted(formula), fitted(direct), tolerance = 3e-12)
    expect_equal(fitted(object_fit), fitted(direct), tolerance = 3e-12)
    expect_equal(se(object_fit), se(direct), tolerance = 3e-12)
    expect_equal(as.numeric(prediction$fit), fitted(direct), tolerance = 3e-12)
    expect_equal(as.numeric(prediction$se.fit), se(direct), tolerance = 3e-12)
  }
})

test_that("adaptive beta regression train-is-eval uses self-excluding NN distances", {
  training <- data.frame(x = c(0.02, 0.08, 0.21, 0.43, 0.68, 0.84, 0.97))
  response <- c(1, 1.5, 0.7, 2.1, 1.2, 2.8, 2.4)
  common <- list(
    bws = 2, bwtype = "adaptive_nn",
    ckertype = "beta", ckerorder = 2,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  fit <- do.call(npreg, c(list(
    txdat = training, tydat = response, regtype = "lc"
  ), common))
  sums <- do.call(npksum, c(list(
    txdat = training, return.kernel.weights = TRUE
  ), common))
  normalized <- sweep(sums$kw, 2L, colSums(sums$kw), "/")

  expect_equal(fitted(fit), colSums(normalized * response), tolerance = 3e-12)
})
