test_that("higher-order beta densities match kernel-sum contributions", {
  training <- data.frame(x = c(0.02, 0.09, 0.22, 0.44, 0.69, 0.86, 0.98))
  evaluation <- data.frame(x = c(0.01, 0.13, 0.38, 0.74, 0.99))

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bandwidth <- if (identical(bwtype, "fixed")) 0.14 else 3
      common <- list(
        bws = bandwidth, bwtype = bwtype,
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1
      )
      fit <- do.call(npudens, c(list(
        tdat = training, edat = evaluation
      ), common))
      sums <- do.call(npksum, c(list(
        txdat = training, exdat = evaluation,
        return.kernel.weights = TRUE
      ), common))
      expected <- colMeans(sums$kw)
      expected_se <- sqrt(apply(sums$kw, 2L, var) / nrow(training))

      expect_equal(fitted(fit), expected, tolerance = 4e-12)
      expect_equal(se(fit), expected_se, tolerance = 4e-12)
    }
  }
})

test_that("higher-order beta distributions match kernel-sum contributions", {
  training <- data.frame(x = c(0.02, 0.09, 0.22, 0.44, 0.69, 0.86, 0.98))
  evaluation <- data.frame(x = c(0, 0.01, 0.13, 0.38, 0.74, 0.99, 1))

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bandwidth <- if (identical(bwtype, "fixed")) 0.14 else 3
      common <- list(
        bws = bandwidth, bwtype = bwtype,
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1
      )
      fit <- do.call(npudist, c(list(
        tdat = training, edat = evaluation
      ), common))
      sums <- do.call(npksum, c(list(
        txdat = training, exdat = evaluation,
        operator = "integral", return.kernel.weights = TRUE
      ), common))
      expected <- colMeans(sums$kw)
      expected_se <- sqrt(colSums(sweep(sums$kw, 2L, expected, "-")^2) /
                            (nrow(training) * (nrow(training) - 1L)))

      expect_equal(fitted(fit), expected, tolerance = 4e-12)
      expect_equal(se(fit), expected_se, tolerance = 4e-12)
      expect_identical(fitted(fit)[c(1L, nrow(evaluation))], c(0, 1))
    }
  }
})

test_that("higher-order beta bandwidth objects preserve order and NN type", {
  training <- data.frame(x = c(0.03, 0.12, 0.28, 0.53, 0.77, 0.96))
  evaluation <- data.frame(x = c(0.01, 0.21, 0.49, 0.82, 0.99))

  for (order in c(4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bandwidth <- if (identical(bwtype, "fixed")) 0.15 else 2
      common <- list(
        bws = bandwidth, bandwidth.compute = FALSE,
        bwtype = bwtype, ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1
      )
      density_bw <- do.call(npudensbw, c(list(dat = training), common))
      distribution_bw <- do.call(npudistbw, c(list(dat = training), common))
      density <- npudens(
        bws = density_bw, tdat = training, edat = evaluation
      )
      distribution <- npudist(
        bws = distribution_bw, tdat = training, edat = evaluation
      )

      expect_identical(density_bw$ckerorder, order)
      expect_identical(distribution_bw$ckerorder, order)
      expect_identical(density_bw$type, bwtype)
      expect_identical(distribution_bw$type, bwtype)
      expect_true(all(is.finite(fitted(density))))
      expect_true(all(is.finite(fitted(distribution))))
    }
  }
})

test_that("higher-order beta estimators do not clamp signed results", {
  density_training <- data.frame(x = c(0.498, 0.499, 0.5, 0.501, 0.502))
  distribution_training <- data.frame(x = c(0.078, 0.079, 0.08, 0.081, 0.082))
  grid <- data.frame(x = seq(0, 1, length.out = 2001L))

  for (order in c(4L, 6L, 8L)) {
    density <- npudens(
      bws = 0.12, tdat = density_training, edat = grid,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )
    distribution <- npudist(
      bws = 0.12, tdat = distribution_training, edat = grid,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )

    expect_lt(min(fitted(density)), 0)
    expect_true(min(fitted(distribution)) < 0 ||
                  max(fitted(distribution)) > 1)
    expect_identical(fitted(distribution)[c(1L, nrow(grid))], c(0, 1))
  }
})
