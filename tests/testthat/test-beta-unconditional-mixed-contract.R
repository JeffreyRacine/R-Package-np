test_that("mixed beta density uses the canonical categorical row", {
  old <- options(np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)
  training <- data.frame(
    x = c(0.01, 0.04, 0.11, 0.2, 0.31, 0.43, 0.56, 0.68, 0.79, 0.88,
          0.95, 0.99),
    u = factor(rep(letters[1:3], each = 4L)),
    o = ordered(rep(1:4, length.out = 12L), levels = 1:4)
  )
  evaluation <- training[c(2L, 5L, 9L, 12L), , drop = FALSE]

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bandwidth <- if (identical(bwtype, "fixed")) {
      c(0.13, 0.21, 0.29)
    } else {
      c(5, 0.21, 0.29)
    }
    for (order in c(2L, 8L)) {
      bw <- npudensbw(
        dat = training, bws = bandwidth, bandwidth.compute = FALSE,
        bwtype = bwtype, bwscaling = FALSE,
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1
      )
      dense <- npudens(bws = bw, tdat = training, edat = evaluation)
      ## Density/distribution routes use the normalized ordered Li--Racine
      ## kernel.  npksum() names that explicit low-level variant nliracine.
      oracle_bw <- bw
      oracle_bw[["okertype"]] <- "nliracine"
      sum_one <- npksum(
        bws = oracle_bw, txdat = training, exdat = evaluation
      )$ksum
      sum_two <- npksum(
        bws = oracle_bw, txdat = training, exdat = evaluation,
        kernel.pow = 2L
      )$ksum
      options(np.categorical.compress = TRUE)
      compressed <- npudens(bws = bw, tdat = training, edat = evaluation)
      options(np.categorical.compress = FALSE)

      expected <- as.vector(sum_one) / nrow(training)
      expected_variance <- as.vector(sum_two) / nrow(training) - expected^2
      expected_se <- sqrt(pmax(expected_variance, 0) /
                            (nrow(training) - 1L))
      expect_equal(fitted(dense), expected, tolerance = 3e-13)
      expect_equal(se(dense), expected_se, tolerance = 3e-13)
      expect_identical(fitted(dense), fitted(compressed))
      expect_identical(se(dense), se(compressed))
    }
  }
})

test_that("mixed beta distribution uses the canonical categorical row", {
  old <- options(np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)
  training <- data.frame(
    x = c(0.01, 0.04, 0.11, 0.2, 0.31, 0.43, 0.56, 0.68, 0.79, 0.88,
          0.95, 0.99),
    o = ordered(rep(1:4, length.out = 12L), levels = 1:4)
  )
  evaluation <- training[c(2L, 5L, 9L, 12L), , drop = FALSE]

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bandwidth <- if (identical(bwtype, "fixed")) {
      c(0.13, 0.29)
    } else {
      c(5, 0.29)
    }
    for (order in c(2L, 8L)) {
      bw <- npudistbw(
        dat = training, bws = bandwidth, bandwidth.compute = FALSE,
        bwtype = bwtype, bwscaling = FALSE,
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1
      )
      dense <- npudist(bws = bw, tdat = training, edat = evaluation)
      oracle_bw <- bw
      oracle_bw[["okertype"]] <- "nliracine"
      weights <- npksum(
        bws = oracle_bw, txdat = training, exdat = evaluation,
        operator = c("integral", "integral"),
        return.kernel.weights = TRUE
      )$kw
      options(np.categorical.compress = TRUE)
      compressed <- npudist(bws = bw, tdat = training, edat = evaluation)
      options(np.categorical.compress = FALSE)

      expected <- colMeans(weights)
      expected_se <- sqrt(
        colSums(sweep(weights, 2L, expected, "-")^2) /
          (nrow(training) * (nrow(training) - 1L))
      )
      expect_equal(fitted(dense), expected, tolerance = 3e-13)
      expect_equal(se(dense), expected_se, tolerance = 3e-13)
      expect_identical(fitted(dense), fitted(compressed))
      expect_identical(se(dense), se(compressed))
    }
  }
})

test_that("mixed beta density objectives are active while distribution stays closed", {
  training <- data.frame(
    x = seq(0.05, 0.95, length.out = 12L),
    u = factor(rep(letters[1:3], each = 4L)),
    o = ordered(rep(1:4, length.out = 12L), levels = 1:4)
  )

  bw <- npudensbw(
    dat = training, bws = c(0.15, 0.2, 0.25),
    bandwidth.compute = FALSE, bwmethod = "cv.ml", bwscaling = FALSE,
    ckertype = "beta", ckerorder = 4,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  objective <- np:::npudensbw.bandwidth(
    dat = training, bws = bw, bandwidth.compute = TRUE,
    eval.only = TRUE, nmulti = 1L, invalid.penalty = "dbmax"
  )
  expect_true(is.finite(as.numeric(objective$fval[[1L]])))

  expect_error(
    npudistbw(
      dat = training[c("x", "o")], bandwidth.compute = TRUE,
      ckertype = "beta", ckerbound = "fixed", ckerlb = 0, ckerub = 1
    ),
    "continuous variables only",
    fixed = TRUE
  )
})

test_that("native density objective ingress requires compression state", {
  training <- data.frame(
    x = seq(0.05, 0.95, length.out = 12L),
    u = factor(rep(letters[1:3], each = 4L)),
    o = ordered(rep(1:4, length.out = 12L), levels = 1:4)
  )
  bw <- npudensbw(
    dat = training, bws = c(0.15, 0.2, 0.25),
    bandwidth.compute = FALSE, bwmethod = "cv.ml", bwscaling = FALSE,
    ckertype = "beta", ckerorder = 4,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  prepare <- getFromNamespace(
    ".npudensbw_nomad_native_prepare_args", "np"
  )
  prep <- prepare(training, bw, invalid.penalty = "dbmax")
  native_eval <- function(myopti) {
    .Call(
      "C_np_density_nomad_native_fixed_eval",
      prep$duno, prep$dord, prep$dcon, prep$mysd,
      myopti, prep$myoptd, as.double(bw$bw),
      prep$penalty_mode, prep$penalty_multiplier,
      prep$ckerlb, prep$ckerub,
      PACKAGE = "np"
    )
  }

  expect_error(
    native_eval(prep$myopti[-length(prep$myopti)]),
    "categorical-compression state is missing",
    fixed = TRUE
  )
  valid <- native_eval(prep$myopti)
  expect_identical(valid$status, 0L)
  expect_true(is.finite(valid$objective))
})

test_that("native distribution objective ingress requires compression state", {
  training <- data.frame(x = seq(0.05, 0.95, length.out = 12L))
  bw <- npudistbw(
    dat = training, bws = 0.15,
    bandwidth.compute = FALSE, bwmethod = "cv.cdf",
    bwscaling = FALSE, ckertype = "beta", ckerorder = 4,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  prepare <- getFromNamespace(
    ".npudistbw_nomad_native_prepare_args", "np"
  )
  prep <- prepare(training, bw, invalid.penalty = "dbmax")
  native_eval <- function(myopti) {
    .Call(
      "C_np_distribution_bw_eval",
      prep$duno, prep$dord, prep$dcon,
      prep$guno, prep$gord, prep$gcon, prep$mysd,
      myopti, prep$myoptd, as.double(bw$bw), 1L,
      prep$penalty_mode, prep$penalty_multiplier,
      prep$ckerlb, prep$ckerub,
      PACKAGE = "np"
    )
  }

  expect_error(
    native_eval(prep$myopti[-length(prep$myopti)]),
    "categorical-compression state is missing",
    fixed = TRUE
  )
  valid <- native_eval(prep$myopti)
  expect_true(is.finite(valid$fval[[1L]]))
})

test_that("native MADS searches mixed beta density objectives", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  set.seed(9134)
  training <- data.frame(
    x = runif(24L),
    u = factor(rep(letters[1:3], length.out = 24L)),
    o = ordered(rep(1:4, length.out = 24L), levels = 1:4)
  )

  for (method in c("cv.ml", "cv.ls")) {
    bw <- npudensbw(
      dat = training, bwmethod = method,
      bwtype = "fixed", bwscaling = FALSE,
      ckertype = "beta", ckerorder = 4,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      bwsolver = "mads", nmulti = 1L,
      nomad.opts = list(MAX_BB_EVAL = 12)
    )
    fit <- npudens(
      bws = bw, tdat = training,
      edat = training[seq_len(4L), , drop = FALSE]
    )
    expect_true(all(is.finite(bw$bw)))
    expect_true(is.finite(as.numeric(bw$fval[[1L]])))
    expect_true(all(is.finite(fitted(fit))))
  }
})

test_that("mixed beta fits reject invalid categorical compression state", {
  training <- data.frame(
    x = seq(0.05, 0.95, length.out = 12L),
    u = factor(rep(letters[1:3], each = 4L))
  )
  bw <- npudensbw(
    dat = training, bws = c(0.15, 0.2), bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  old <- options(np.categorical.compress = NA)
  on.exit(options(old), add = TRUE)

  expect_error(
    npudens(bws = bw, tdat = training),
    "must be a single non-missing logical value",
    fixed = TRUE
  )
})
