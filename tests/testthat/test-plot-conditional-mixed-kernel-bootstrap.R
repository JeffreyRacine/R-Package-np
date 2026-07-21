test_that("beta conditional fixed bootstrap matches the signed-log estimator", {
  xdat <- data.frame(x = c(
    0.02340804, 0.19247293, 0.34607979, 0.43130701,
    0.50471434, 0.52169417, 0.63018924, 0.66711685,
    0.74800750, 0.74852143, 0.79852423, 0.88447980
  ))
  ydat <- data.frame(y = c(
    0.4140633, 0.2824111, 0.4724171, 0.1185962,
    0.2713280, 0.3527147, 0.8634327, 0.4683063,
    0.1998146, 0.1090273, 0.9564857, 0.4932104
  ))
  exdat <- data.frame(x = c(
    1e-8, 1e-5, 1e-3, 0.64429388, 0.85635732,
    0.49983058, 0.42668611, 0.99465242, 0.999, 0.99999999
  ))
  eydat <- data.frame(y = c(
    0.47083275, 0.56110883, 0.03229574, 0.19628559, 0.09090910,
    0.43176428, 0.09647046, 0.80853377, 0.93223586, 0.93556179
  ))
  bw <- npcdensbw(
    xdat = xdat, ydat = ydat,
    bws = c(0.51419008, 0.01038156), bandwidth.compute = FALSE,
    cxkertype = "beta", cxkerorder = 8L,
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykertype = "beta", cykerorder = 8L,
    cykerbound = "fixed", cykerlb = 0, cykerub = 1
  )
  truth <- fitted(npcdens(
    bws = bw, txdat = xdat, tydat = ydat,
    exdat = exdat, eydat = eydat
  ))
  counts <- matrix(1, nrow = nrow(xdat), ncol = 1L)
  out <- np:::.np_inid_boot_from_ksum_conditional(
    xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
    bws = bw, B = 1L, cdf = FALSE, counts = counts
  )

  expect_equal(out$t0, truth, tolerance = 2e-11)
  expect_equal(as.numeric(out$t[1L, ]), truth, tolerance = 2e-11)
})

test_that("fixed mixed continuous families use exact conditional bootstrap levels", {
  xdat <- data.frame(x = c(0.04, 0.13, 0.29, 0.46, 0.68, 0.91))
  ydat <- data.frame(y = c(0.05, 0.19, 0.37, 0.55, 0.79, 0.96))
  exdat <- data.frame(x = c(0.09, 0.38, 0.83))
  eydat <- data.frame(y = c(0.08, 0.49, 0.88))
  counts <- cbind(
    c(2, 0, 1, 1, 1, 1),
    c(0, 2, 1, 1, 1, 1)
  )

  for (cdf in c(FALSE, TRUE)) {
    bw.fun <- if (cdf) npcdistbw else npcdensbw
    fit.fun <- if (cdf) npcdist else npcdens
    bw <- bw.fun(
      xdat = xdat, ydat = ydat,
      bws = c(0.18, 0.14), bandwidth.compute = FALSE,
      cxkertype = "epanechnikov", cxkerorder = 4L,
      cykertype = "beta", cykerorder = 6L,
      cykerbound = "fixed", cykerlb = 0, cykerub = 1
    )
    out <- np:::.np_inid_boot_from_ksum_conditional(
      xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
      bws = bw, B = ncol(counts), cdf = cdf, counts = counts
    )
    expect_type(out, "list")

    t0 <- fitted(fit.fun(
      bws = bw, txdat = xdat, tydat = ydat,
      exdat = exdat, eydat = eydat
    ))
    expect_equal(out$t0, t0, tolerance = 2e-11)

    for (jj in seq_len(ncol(counts))) {
      idx <- rep.int(seq_len(nrow(xdat)), counts[, jj])
      oracle <- fitted(fit.fun(
        bws = bw,
        txdat = xdat[idx, , drop = FALSE],
        tydat = ydat[idx, , drop = FALSE],
        exdat = exdat, eydat = eydat
      ))
      expect_equal(as.numeric(out$t[jj, ]), oracle, tolerance = 2e-11)
    }
  }
})

test_that("beta on X and mixed all-legacy fixed levels match public refits", {
  xdat <- data.frame(x = c(0.03, 0.14, 0.31, 0.48, 0.66, 0.93))
  ydat <- data.frame(y = c(0.06, 0.21, 0.39, 0.58, 0.76, 0.95))
  exdat <- data.frame(x = c(0.08, 0.42, 0.87))
  eydat <- data.frame(y = c(0.11, 0.51, 0.90))
  counts <- cbind(
    c(2, 0, 1, 1, 1, 1),
    c(1, 1, 0, 2, 1, 1)
  )
  specifications <- list(
    list(
      cxkertype = "beta", cxkerorder = 4L,
      cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
      cykertype = "gaussian", cykerorder = 8L
    ),
    list(
      cxkertype = "epanechnikov", cxkerorder = 6L,
      cykertype = "gaussian", cykerorder = 4L
    )
  )

  for (specification in specifications) {
    bw <- do.call(npcdensbw, c(list(
      xdat = xdat, ydat = ydat,
      bws = c(0.17, 0.15), bandwidth.compute = FALSE
    ), specification))
    out <- np:::.np_inid_boot_from_ksum_conditional(
      xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
      bws = bw, B = ncol(counts), cdf = FALSE, counts = counts
    )
    truth <- fitted(npcdens(
      bws = bw, txdat = xdat, tydat = ydat,
      exdat = exdat, eydat = eydat
    ))
    expect_equal(out$t0, truth, tolerance = 3e-11)

    for (jj in seq_len(ncol(counts))) {
      idx <- rep.int(seq_len(nrow(xdat)), counts[, jj])
      oracle <- fitted(npcdens(
        bws = bw,
        txdat = xdat[idx, , drop = FALSE],
        tydat = ydat[idx, , drop = FALSE],
        exdat = exdat, eydat = eydat
      ))
      expect_equal(as.numeric(out$t[jj, ]), oracle, tolerance = 3e-10)
    }
  }
})

test_that("beta empirical-range bounds replay exactly in count bootstrap", {
  xdat <- data.frame(x = c(0.08, 0.16, 0.34, 0.57, 0.73, 0.91))
  ydat <- data.frame(y = c(0.04, 0.18, 0.32, 0.61, 0.81, 0.96))
  exdat <- data.frame(x = c(0.10, 0.44, 0.86))
  eydat <- data.frame(y = c(0.07, 0.52, 0.90))
  counts <- cbind(c(2, 0, 1, 1, 1, 1), c(1, 1, 1, 1, 0, 2))
  bw <- npcdensbw(
    xdat = xdat, ydat = ydat,
    bws = c(0.14, 0.13), bandwidth.compute = FALSE,
    cxkertype = "beta", cxkerorder = 4L, cxkerbound = "range",
    cykertype = "beta", cykerorder = 6L, cykerbound = "range"
  )
  resolved <- list(
    xlb = bw$cxkerlb, xub = bw$cxkerub,
    ylb = bw$cykerlb, yub = bw$cykerub
  )
  out <- np:::.np_inid_boot_from_ksum_conditional(
    xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
    bws = bw, B = ncol(counts), cdf = FALSE, counts = counts
  )

  expect_identical(bw$cxkerlb, resolved$xlb)
  expect_identical(bw$cxkerub, resolved$xub)
  expect_identical(bw$cykerlb, resolved$ylb)
  expect_identical(bw$cykerub, resolved$yub)
  for (jj in seq_len(ncol(counts))) {
    idx <- rep.int(seq_len(nrow(xdat)), counts[, jj])
    oracle <- fitted(npcdens(
      bws = bw,
      txdat = xdat[idx, , drop = FALSE],
      tydat = ydat[idx, , drop = FALSE],
      exdat = exdat, eydat = eydat
    ))
    expect_equal(as.numeric(out$t[jj, ]), oracle, tolerance = 3e-10)
  }
})

test_that("mixed exact and frozen NN levels preserve their declared semantics", {
  xdat <- data.frame(x = c(0.04, 0.12, 0.24, 0.39, 0.55, 0.71, 0.84, 0.96))
  ydat <- data.frame(y = c(0.03, 0.16, 0.27, 0.44, 0.63, 0.74, 0.89, 0.97))
  exdat <- data.frame(x = c(0.09, 0.36, 0.68, 0.91))
  eydat <- data.frame(y = c(0.08, 0.40, 0.69, 0.92))
  counts <- cbind(
    rep.int(1, nrow(xdat)),
    c(2, 0, 1, 1, 1, 1, 1, 1),
    c(1, 1, 2, 0, 1, 1, 1, 1)
  )

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    bw <- npcdensbw(
      xdat = xdat, ydat = ydat,
      bws = c(3L, 3L), bwtype = bwtype, bandwidth.compute = FALSE,
      cxkertype = "epanechnikov", cxkerorder = 4L,
      cykertype = "beta", cykerorder = 6L,
      cykerbound = "fixed", cykerlb = 0, cykerub = 1
    )
    exact <- np:::.np_inid_boot_from_ksum_conditional_exact(
      xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
      bws = bw, B = ncol(counts), cdf = FALSE, counts = counts
    )
    frozen <- np:::.np_inid_boot_from_hat_conditional_frozen(
      xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
      bws = bw, B = ncol(counts), cdf = FALSE, counts = counts
    )
    truth <- fitted(npcdens(
      bws = bw, txdat = xdat, tydat = ydat,
      exdat = exdat, eydat = eydat
    ))
    expect_equal(exact$t0, truth, tolerance = 3e-10)
    expect_equal(frozen$t0, truth, tolerance = 3e-10)

    for (jj in seq_len(ncol(counts))) {
      idx <- rep.int(seq_len(nrow(xdat)), counts[, jj])
      oracle <- fitted(npcdens(
        bws = bw,
        txdat = xdat[idx, , drop = FALSE],
        tydat = ydat[idx, , drop = FALSE],
        exdat = exdat, eydat = eydat
      ))
      expect_equal(as.numeric(exact$t[jj, ]), oracle, tolerance = 3e-9)
    }
    expect_equal(as.numeric(frozen$t[1L, ]), truth, tolerance = 3e-10)
  }
})

test_that("all three conditional resampling methods accept mixed beta levels", {
  set.seed(20260721)
  xdat <- data.frame(x = sort(runif(30)))
  ydat <- data.frame(y = plogis(2 * xdat$x + rnorm(30, sd = 0.4)))
  exdat <- data.frame(x = seq(0.2, 0.8, length.out = 5L))
  eydat <- data.frame(y = seq(0.1, 0.9, length.out = 5L))
  compute <- np:::compute.bootstrap.errors.conbandwidth

  for (cdf in c(FALSE, TRUE)) {
    bw.fun <- if (cdf) npcdistbw else npcdensbw
    bw <- bw.fun(
      xdat = xdat, ydat = ydat,
      bws = c(0.2, 0.15), bandwidth.compute = FALSE,
      cxkertype = "gaussian", cxkerorder = 4L,
      cykertype = "beta", cykerorder = 6L,
      cykerbound = "fixed", cykerlb = 0, cykerub = 1
    )
    for (method in c("inid", "fixed", "geom")) {
      set.seed(100L + match(method, c("inid", "fixed", "geom")))
      out <- compute(
        xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
        cdf = cdf, quantreg = FALSE, tau = NULL,
        gradients = FALSE, gradient.index = 1L, slice.index = 1L,
        plot.errors.boot.method = method,
        plot.errors.boot.nonfixed = "exact",
        plot.errors.boot.blocklen = 3L,
        plot.errors.boot.num = 39L,
        plot.errors.center = "estimate",
        plot.errors.type = "pointwise",
        plot.errors.alpha = 0.05,
        bws = bw
      )
      expect_identical(dim(out$boot.err), c(nrow(exdat), 3L))
      expect_true(all(is.finite(out$boot.err)))
    }
  }
})

test_that("mixed all-legacy exact NN retains duplicate-row refit semantics", {
  xdat <- data.frame(x = c(0.04, 0.13, 0.25, 0.41, 0.58, 0.72, 0.86, 0.96))
  ydat <- data.frame(y = c(0.03, 0.17, 0.29, 0.43, 0.62, 0.77, 0.88, 0.97))
  exdat <- data.frame(x = c(0.11, 0.37, 0.69, 0.90))
  eydat <- data.frame(y = c(0.09, 0.39, 0.71, 0.91))
  counts <- cbind(
    rep.int(1, nrow(xdat)),
    c(2, 0, 1, 1, 1, 1, 1, 1),
    c(1, 1, 2, 0, 1, 1, 1, 1)
  )
  bw <- npcdistbw(
    xdat = xdat, ydat = ydat,
    bws = c(3L, 3L), bwtype = "generalized_nn",
    bandwidth.compute = FALSE,
    cxkertype = "epanechnikov", cxkerorder = 4L,
    cykertype = "gaussian", cykerorder = 8L
  )
  out <- np:::.np_inid_boot_from_ksum_conditional_exact(
    xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
    bws = bw, B = ncol(counts), cdf = TRUE, counts = counts
  )
  for (jj in seq_len(ncol(counts))) {
    idx <- rep.int(seq_len(nrow(xdat)), counts[, jj])
    oracle <- fitted(npcdist(
      bws = bw,
      txdat = xdat[idx, , drop = FALSE],
      tydat = ydat[idx, , drop = FALSE],
      exdat = exdat, eydat = eydat
    ))
    expect_equal(as.numeric(out$t[jj, ]), oracle, tolerance = 3e-9)
  }
})
