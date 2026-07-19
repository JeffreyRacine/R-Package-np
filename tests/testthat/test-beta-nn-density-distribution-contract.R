test_that("nearest-neighbor beta density matches its kernel-sum contributions", {
  training <- data.frame(x = c(0.02, 0.09, 0.22, 0.44, 0.69, 0.86, 0.98))
  evaluation <- data.frame(x = c(0.01, 0.13, 0.38, 0.74, 0.99))

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    fit <- npudens(
      bws = 3, tdat = training, edat = evaluation,
      bwtype = bwtype, ckertype = "beta", ckerorder = 2,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )
    sums <- npksum(
      bws = 3, txdat = training, exdat = evaluation,
      bwtype = bwtype, ckertype = "beta", ckerorder = 2,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      return.kernel.weights = TRUE
    )
    expected <- colMeans(sums$kw)
    expected_se <- sqrt(apply(sums$kw, 2L, var) / nrow(training))

    expect_equal(fitted(fit), expected, tolerance = 3e-13)
    expect_equal(se(fit), expected_se, tolerance = 3e-13)
    expect_equal(fit$log_likelihood,
                 sum(log(pmax(expected, .Machine$double.xmin))),
                 tolerance = 3e-13)
  }
})

test_that("nearest-neighbor beta distribution matches its CDF contributions", {
  training <- data.frame(x = c(0.02, 0.09, 0.22, 0.44, 0.69, 0.86, 0.98))
  evaluation <- data.frame(x = c(0, 0.01, 0.13, 0.38, 0.74, 0.99, 1))

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    fit <- npudist(
      bws = 3, tdat = training, edat = evaluation,
      bwtype = bwtype, ckertype = "beta", ckerorder = 2,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )
    sums <- npksum(
      bws = 3, txdat = training, exdat = evaluation,
      bwtype = bwtype, ckertype = "beta", ckerorder = 2,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      operator = "integral", return.kernel.weights = TRUE
    )
    expected <- colMeans(sums$kw)
    expected_se <- sqrt(colSums(sweep(sums$kw, 2L, expected, "-")^2) /
                          (nrow(training) * (nrow(training) - 1L)))

    expect_equal(fitted(fit), expected, tolerance = 3e-13)
    expect_equal(se(fit), expected_se, tolerance = 3e-13)
    expect_identical(fitted(fit)[c(1L, nrow(evaluation))], c(0, 1))
  }
})

test_that("nearest-neighbor beta density and distribution retain bandwidth objects", {
  training <- data.frame(x = c(0.03, 0.12, 0.28, 0.53, 0.77, 0.96))
  evaluation <- data.frame(x = c(0.01, 0.21, 0.49, 0.82, 0.99))

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    density_bw <- npudensbw(
      dat = training, bws = 2, bandwidth.compute = FALSE,
      bwtype = bwtype, ckertype = "beta", ckerorder = 2,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )
    distribution_bw <- npudistbw(
      dat = training, bws = 2, bandwidth.compute = FALSE,
      bwtype = bwtype, ckertype = "beta", ckerorder = 2,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )
    density <- npudens(bws = density_bw, tdat = training, edat = evaluation)
    distribution <- npudist(
      bws = distribution_bw, tdat = training, edat = evaluation
    )

    expect_s3_class(density_bw, "bandwidth")
    expect_s3_class(distribution_bw, "dbandwidth")
    expect_identical(density_bw$type, bwtype)
    expect_identical(distribution_bw$type, bwtype)
    expect_true(all(is.finite(fitted(density))))
    expect_true(all(is.finite(fitted(distribution))))
  }
})

test_that("generalized-NN beta distribution keeps documented nonmonotonicity", {
  training <- data.frame(x = seq(0.086, 0.146, length.out = 12L))
  evaluation <- data.frame(x = seq(0, 1, length.out = 2001L))
  fit <- npudist(
    bws = 4, tdat = training, edat = evaluation,
    bwtype = "generalized_nn", ckertype = "beta",
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  expect_lt(min(diff(fitted(fit))), -1e-5)
  expect_true(all(fitted(fit) >= 0 & fitted(fit) <= 1))
  expect_identical(fitted(fit)[c(1L, nrow(evaluation))], c(0, 1))
})
