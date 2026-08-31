hc0_hat_oracle <- function(bws, txdat, tydat, exdat = NULL) {
  H.train <- unclass(npreghat(
    bws = bws,
    txdat = txdat,
    output = "matrix"
  ))
  H.eval <- if (is.null(exdat)) {
    H.train
  } else {
    unclass(npreghat(
      bws = bws,
      txdat = txdat,
      exdat = exdat,
      output = "matrix"
    ))
  }
  residual <- as.double(tydat) - drop(H.train %*% as.double(tydat))

  sqrt(drop((H.eval^2) %*% (residual^2)))
}

hc0_explicit_lc_bw <- function(xdat, ydat, bws, bwtype = "fixed", ...) {
  npregbw(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = FALSE,
    bwscaling = FALSE,
    bwtype = bwtype,
    regtype = "lc",
    ...
  )
}

test_that("scalar fixed, GNN, and adaptive-NN mean SEs equal the HC0 oracle", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(x = seq(-1.2, 1.3, length.out = 31L))
  ydat <- sin(1.4 * xdat$x) + seq_len(nrow(xdat)) / 50
  exdat <- data.frame(x = c(-1.1, -0.35, 0.2, 1.15))

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- hc0_explicit_lc_bw(
      xdat,
      ydat,
      bws = if (identical(bwtype, "fixed")) 0.43 else 7,
      bwtype = bwtype
    )

    for (external in c(FALSE, TRUE)) {
      eval <- if (external) exdat else NULL
      args <- list(bws = bw, txdat = xdat, tydat = ydat)
      if (external)
        args$exdat <- eval
      fit.no.se <- do.call(npreg, c(args, list(se = FALSE)))
      fit.se <- do.call(npreg, c(args, list(se = TRUE)))
      oracle <- hc0_hat_oracle(bw, xdat, ydat, eval)

      expect_identical(fit.se$mean, fit.no.se$mean)
      expect_identical(fit.se$xtra, fit.no.se$xtra)
      expect_equal(fit.se$merr, oracle, tolerance = 5e-14)
      expect_true(all(is.finite(fit.se$merr)))
      expect_true(all(fit.se$merr >= 0))
    }
  }
})

test_that("beta scalar mean SEs equal the HC0 oracle across orders and NN modes", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(x = c(0.03, 0.10, 0.24, 0.43, 0.67, 0.82, 0.96))
  ydat <- sin(2 * xdat$x) + xdat$x +
    c(-0.03, 0.02, 0.01, -0.02, 0.04, -0.01, 0.02)
  exdat <- data.frame(x = c(0.08, 0.38, 0.74, 0.92))

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bw <- hc0_explicit_lc_bw(
        xdat,
        ydat,
        bws = if (identical(bwtype, "fixed")) 0.18 else 3,
        bwtype = bwtype,
        ckertype = "beta",
        ckerorder = order,
        ckerbound = "fixed",
        ckerlb = 0,
        ckerub = 1
      )

      for (external in c(FALSE, TRUE)) {
        eval <- if (external) exdat else NULL
        args <- list(bws = bw, txdat = xdat, tydat = ydat)
        if (external)
          args$exdat <- eval
        fit.no.se <- do.call(npreg, c(args, list(se = FALSE)))
        fit.se <- do.call(npreg, c(args, list(se = TRUE)))

        expect_identical(fit.se$mean, fit.no.se$mean)
        expect_equal(
          fit.se$merr,
          hc0_hat_oracle(bw, xdat, ydat, eval),
          tolerance = 2e-14
        )
        expect_true(all(is.finite(fit.se$merr)))
        expect_true(all(fit.se$merr >= 0))
      }
    }
  }
})

test_that("mixed beta scalar HC0 preserves dense and compressed point routes", {
  old <- options(np.messages = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  x <- seq(0.025, 0.975, length.out = 24L)
  xdat <- data.frame(
    x = x,
    u = factor(rep(letters[1:3], length.out = length(x))),
    o = ordered(rep(1:4, each = 6L))
  )
  ydat <- sin(4 * x) + 0.12 * as.integer(xdat$u) +
    0.07 * as.integer(xdat$o)
  exdat <- xdat[c(2L, 5L, 9L, 14L, 20L, 23L), , drop = FALSE]

  for (order in c(2L, 6L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bw <- hc0_explicit_lc_bw(
        xdat,
        ydat,
        bws = if (identical(bwtype, "fixed")) {
          c(0.15, 0.2, 0.25)
        } else {
          c(7, 0.2, 0.25)
        },
        bwtype = bwtype,
        ckertype = "beta",
        ckerorder = order,
        ckerbound = "fixed",
        ckerlb = 0,
        ckerub = 1
      )

      for (external in c(FALSE, TRUE)) {
        eval <- if (external) exdat else NULL
        args <- list(
          bws = bw,
          txdat = xdat,
          tydat = ydat,
          gradients = TRUE
        )
        if (external)
          args$exdat <- eval
        fit.no.se <- do.call(npreg, c(args, list(se = FALSE)))

        options(np.categorical.compress = FALSE)
        dense <- do.call(npreg, c(args, list(se = TRUE)))
        options(np.categorical.compress = TRUE)
        compressed <- do.call(npreg, c(args, list(se = TRUE)))

        expect_identical(dense$mean, fit.no.se$mean)
        expect_identical(dense$grad, fit.no.se$grad)
        expect_equal(
          dense$merr,
          hc0_hat_oracle(bw, xdat, ydat, eval),
          tolerance = 3e-14
        )
        expect_true(all(is.finite(dense$gerr[, bw$icon, drop = FALSE])))
        expect_true(all(is.finite(dense$gerr[, !bw$icon, drop = FALSE])))
        expect_identical(compressed$mean, dense$mean)
        expect_identical(compressed$merr, dense$merr)
        expect_identical(compressed$grad, dense$grad)
        expect_identical(compressed$gerr, dense$gerr)
      }
    }
  }
})

test_that("beta scalar HC0 survives complete raw-weight underflow", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(x = 0.9 + c(0, 1e-7, 2e-7))
  ydat <- c(1, 4, 9)
  exdat <- data.frame(x = 0)

  for (order in c(2L, 4L, 6L, 8L)) {
    bw <- hc0_explicit_lc_bw(
      xdat,
      ydat,
      bws = 0.001,
      ckertype = "beta",
      ckerorder = order,
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1
    )
    fit.no.se <- npreg(
      bws = bw, txdat = xdat, tydat = ydat, exdat = exdat, se = FALSE
    )
    fit.se <- npreg(
      bws = bw, txdat = xdat, tydat = ydat, exdat = exdat, se = TRUE
    )

    expect_identical(fit.se$mean, fit.no.se$mean)
    expect_equal(
      fit.se$merr,
      hc0_hat_oracle(bw, xdat, ydat, exdat),
      tolerance = 3e-12
    )
    expect_true(all(is.finite(fit.se$merr)))
  }
})

test_that("compact scalar kernels retain the HC0 row identity", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(
    x = seq(-1.25, 1.35, length.out = 37L),
    z = factor(rep(c("a", "b", "c"), length.out = 37L))
  )
  ydat <- sin(1.3 * xdat$x) + c(a = 0, b = 0.4, c = -0.25)[xdat$z] +
    seq_len(nrow(xdat)) / 80
  exdat <- xdat[c(2L, 9L, 16L, 25L, 34L), , drop = FALSE]

  for (kernel in c("epanechnikov", "uniform")) {
    bw <- suppressWarnings(hc0_explicit_lc_bw(
      xdat,
      ydat,
      c(0.44, 0.19),
      ckertype = kernel
    ))
    fit <- suppressWarnings(npreg(
      bws = bw,
      txdat = xdat,
      tydat = ydat,
      exdat = exdat,
      se = TRUE
    ))

    expect_equal(
      fit$merr,
      suppressWarnings(hc0_hat_oracle(bw, xdat, ydat, exdat)),
      tolerance = 5e-14
    )
  }
})

test_that("mixed and categorical-only scalar HC0 follows the same hat-row algebra", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  mixed <- data.frame(
    z = factor(rep(c("a", "b", "c"), c(11L, 13L, 12L))),
    x = seq(-1, 1, length.out = 36L)
  )
  y.mixed <- 0.4 + 1.1 * mixed$x +
    c(a = 0, b = 0.6, c = -0.3)[mixed$z] + sin(seq_len(36L)) / 8
  bw.mixed <- hc0_explicit_lc_bw(mixed, y.mixed, c(0.47, 0.18))
  fit.mixed <- npreg(
    bws = bw.mixed,
    txdat = mixed,
    tydat = y.mixed,
    se = TRUE
  )

  expect_equal(
    fit.mixed$merr,
    hc0_hat_oracle(bw.mixed, mixed, y.mixed),
    tolerance = 5e-14
  )

  categorical <- data.frame(
    z = factor(rep(c("a", "b", "c"), c(8L, 15L, 13L)))
  )
  y.categorical <- c(
    seq(-0.2, 0.5, length.out = 8L),
    seq(0.4, 1.4, length.out = 15L),
    seq(-1, 0.1, length.out = 13L)
  )
  bw.categorical <- hc0_explicit_lc_bw(categorical, y.categorical, 0.12)
  fit.categorical.no.se <- npreg(
    bws = bw.categorical,
    txdat = categorical,
    tydat = y.categorical,
    se = FALSE
  )
  fit.categorical <- npreg(
    bws = bw.categorical,
    txdat = categorical,
    tydat = y.categorical,
    se = TRUE
  )

  expect_identical(fit.categorical$mean, fit.categorical.no.se$mean)
  expect_equal(
    fit.categorical$merr,
    hc0_hat_oracle(bw.categorical, categorical, y.categorical),
    tolerance = 5e-14
  )
  expect_true(any(fit.categorical$merr > 0))

  cat.eval <- data.frame(
    z = factor(c("a", "b", "c"), levels = levels(categorical$z))
  )
  fit.cat.eval.no.se <- npreg(
    bws = bw.categorical,
    txdat = categorical,
    tydat = y.categorical,
    exdat = cat.eval,
    se = FALSE
  )
  fit.cat.eval <- npreg(
    bws = bw.categorical,
    txdat = categorical,
    tydat = y.categorical,
    exdat = cat.eval,
    se = TRUE
  )
  expect_identical(fit.cat.eval$mean, fit.cat.eval.no.se$mean)
  expect_equal(
    fit.cat.eval$merr,
    hc0_hat_oracle(bw.categorical, categorical, y.categorical, cat.eval),
    tolerance = 5e-14
  )
})

test_that("the scalar all-large point shortcut cedes without point drift", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(x = seq(-2, 2, length.out = 25L))
  ydat <- 1.1 + 0.65 * xdat$x +
    c(-0.12, 0.03, 0.08, -0.06, 0.04, 0.11, -0.02, -0.07, 0.05, 0.01,
      -0.09, 0.06, 0.02, -0.04, 0.07, -0.01, 0.09, -0.08, 0.03, 0.05,
      -0.03, 0.04, -0.05, 0.08, -0.02)
  exdat <- data.frame(x = c(-1.5, 0, 1.5))
  bw <- hc0_explicit_lc_bw(xdat, ydat, 1e16)

  fit.no.se <- npreg(
    bws = bw,
    txdat = xdat,
    tydat = ydat,
    exdat = exdat,
    se = FALSE
  )
  fit.se <- npreg(
    bws = bw,
    txdat = xdat,
    tydat = ydat,
    exdat = exdat,
    se = TRUE
  )

  expect_identical(fit.se$mean, fit.no.se$mean)
  expect_identical(fit.se$xtra, fit.no.se$xtra)
  expect_equal(
    fit.se$merr,
    hc0_hat_oracle(bw, xdat, ydat, exdat),
    tolerance = 5e-14
  )
})

test_that("scalar HC0 preserves point gradients and exposes every gerr", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(
    z = factor(rep(c("a", "b"), each = 10L)),
    x = seq(-1, 1, length.out = 20L)
  )
  ydat <- sin(xdat$x) + 0.5 * as.numeric(xdat$z == "b") +
    seq_len(nrow(xdat)) / 100
  exdat <- xdat[c(2L, 7L, 13L, 18L), , drop = FALSE]
  bw <- hc0_explicit_lc_bw(xdat, ydat, c(0.4, 0.2))

  for (external in c(FALSE, TRUE)) {
    args <- list(
      bws = bw,
      txdat = xdat,
      tydat = ydat,
      gradients = TRUE
    )
    if (external)
      args$exdat <- exdat
    fit.no.se <- do.call(npreg, c(args, list(se = FALSE)))
    fit.se <- do.call(npreg, c(args, list(se = TRUE)))

    expect_identical(fit.se$mean, fit.no.se$mean)
    expect_identical(fit.se$grad, fit.no.se$grad)
    expect_identical(fit.se$xtra, fit.no.se$xtra)
    expect_true(all(is.finite(fit.se$gerr[, bw$icon, drop = FALSE])))
    expect_true(all(is.finite(fit.se$gerr[, !bw$icon, drop = FALSE])))
  }
})

test_that("predict se.fit consumes the canonical scalar HC0 result", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(
    x = seq(-1, 1, length.out = 29L),
    z = factor(rep(c("a", "b", "c"), length.out = 29L))
  )
  ydat <- cos(1.2 * xdat$x) +
    c(a = -0.2, b = 0.35, c = 0.1)[xdat$z] +
    seq_len(nrow(xdat)) / 90
  exdat <- xdat[c(3L, 10L, 18L, 27L), , drop = FALSE]
  bw <- hc0_explicit_lc_bw(xdat, ydat, c(0.42, 0.16))
  fit <- npreg(bws = bw, txdat = xdat, tydat = ydat, se = FALSE)
  direct <- npreg(
    bws = bw,
    txdat = xdat,
    tydat = ydat,
    exdat = exdat,
    se = TRUE
  )
  predicted <- predict(fit, newdata = exdat, se.fit = TRUE)

  expect_identical(predicted$fit, direct$mean)
  expect_identical(predicted$se.fit, direct$merr)
  expect_equal(
    predicted$se.fit,
    hc0_hat_oracle(bw, xdat, ydat, exdat),
    tolerance = 5e-14
  )
})

test_that("scalar HC0 respects tree, donor-order, and label invariants", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(
    x = c(seq(-1.2, -0.1, length.out = 15L),
          seq(0.03, 1.4, length.out = 18L)),
    z = factor(rep(c("low", "mid", "high"), length.out = 33L),
               levels = c("low", "mid", "high"))
  )
  ydat <- sin(1.5 * xdat$x) +
    c(low = -0.25, mid = 0.15, high = 0.5)[xdat$z] +
    seq_len(nrow(xdat)) / 110
  bw <- hc0_explicit_lc_bw(xdat, ydat, c(0.39, 0.17))

  options(np.tree = FALSE)
  tree.off <- npreg(bws = bw, txdat = xdat, tydat = ydat, se = TRUE)
  options(np.tree = TRUE)
  tree.on <- npreg(bws = bw, txdat = xdat, tydat = ydat, se = TRUE)

  expect_identical(tree.on$mean, tree.off$mean)
  expect_identical(tree.on$merr, tree.off$merr)

  permutation <- c(seq(2L, 32L, by = 2L), seq(1L, 33L, by = 2L))
  permuted.x <- xdat[permutation, , drop = FALSE]
  permuted.y <- ydat[permutation]
  permuted.bw <- hc0_explicit_lc_bw(permuted.x, permuted.y, c(0.39, 0.17))
  permuted <- npreg(
    bws = permuted.bw,
    txdat = permuted.x,
    tydat = permuted.y,
    se = TRUE
  )
  inverse <- order(permutation)

  expect_equal(permuted$mean[inverse], tree.on$mean, tolerance = 2e-14)
  expect_equal(permuted$merr[inverse], tree.on$merr, tolerance = 2e-14)

  renamed.x <- xdat
  renamed.x$z <- factor(
    c(low = "A", mid = "B", high = "C")[xdat$z],
    levels = c("A", "B", "C")
  )
  renamed.bw <- hc0_explicit_lc_bw(renamed.x, ydat, c(0.39, 0.17))
  renamed <- npreg(
    bws = renamed.bw,
    txdat = renamed.x,
    tydat = ydat,
    se = TRUE
  )

  expect_identical(renamed$mean, tree.on$mean)
  expect_identical(renamed$merr, tree.on$merr)
})

test_that("scalar HC0 is response-equivariant and constant responses are exact zero", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(
    x = seq(-1, 1, length.out = 31L),
    z = factor(rep(c("a", "b", "c"), length.out = 31L))
  )
  ydat <- sin(xdat$x) + seq_len(nrow(xdat)) / 40
  bw <- hc0_explicit_lc_bw(xdat, ydat, c(0.4, 0.2))
  base <- npreg(bws = bw, txdat = xdat, tydat = ydat, se = TRUE)
  shifted <- npreg(bws = bw, txdat = xdat, tydat = ydat + 17, se = TRUE)
  scaled <- npreg(bws = bw, txdat = xdat, tydat = -3.5 * ydat, se = TRUE)
  near.limit <- npreg(
    bws = bw,
    txdat = xdat,
    tydat = 1e150 * ydat,
    se = TRUE
  )
  constant <- npreg(
    bws = bw,
    txdat = xdat,
    tydat = rep(5.25, nrow(xdat)),
    se = TRUE
  )

  expect_equal(shifted$merr, base$merr, tolerance = 2e-13)
  expect_equal(scaled$merr, 3.5 * base$merr, tolerance = 2e-13)
  expect_true(all(is.finite(near.limit$merr)))
  expect_equal(near.limit$merr / 1e150, base$merr, tolerance = 2e-13)
  expect_identical(constant$merr, rep(0, nrow(xdat)))
})

test_that("HC0 activation remains private and linear-memory", {
  source <- paste(
    readLines(test_path("..", "..", "src", "np.c"), warn = FALSE),
    collapse = "\n"
  )
  reducer <- paste(
    readLines(test_path("..", "..", "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  hc0.start <- regexpr(
    "memset(&ordinary_hc0_context, 0, sizeof(ordinary_hc0_context));",
    source,
    fixed = TRUE
  )
  hc0.tail <- substring(source, hc0.start)
  hc0.end <- regexpr(
    "np_progress_fit_set_offset(num_obs_train_extern);",
    hc0.tail,
    fixed = TRUE
  )
  hc0.source <- substring(
    hc0.tail,
    1L,
    hc0.end + nchar("np_progress_fit_set_offset(num_obs_train_extern);") - 1L
  )

  expect_match(source, "NPRegressionHC0Context ordinary_hc0_context", fixed = TRUE)
  expect_match(
    source,
    paste0(
      "ordinary_hc0_active = do_merr &&\n",
      "    (np_lp_engine_extern == NP_LP_ENGINE_SCALAR ||\n",
      "     np_lp_engine_extern == NP_LP_ENGINE_GENERAL);"
    ),
    fixed = TRUE
  )
  expect_gt(hc0.start, 0L)
  expect_gt(hc0.end, 0L)
  expect_false(grepl("npreghat", hc0.source, fixed = TRUE))
  expect_false(grepl("np_regression_lp_hat_matrix", hc0.source, fixed = TRUE))
  expect_match(reducer, "ordinary_hc0 ? &hc0_dual_power_ctx : NULL", fixed = TRUE)
  expect_match(
    reducer,
    "regression_moment_context.hc0_scaled_residual",
    fixed = TRUE
  )
  expect_match(reducer, "hc0_context == NULL", fixed = TRUE)
})

test_that("the compiled-fit total accounts for both HC0 phases only", {
  ordinary <- list(type = "fixed", ckertype = "gaussian")
  beta <- list(type = "fixed", ckertype = "beta")
  adaptive <- list(type = "adaptive_nn", ckertype = "gaussian")

  expect_identical(
    .np_reg_fit_total(ordinary, 24L, 3L, TRUE, REGTYPE_LP0),
    27L
  )
  expect_identical(
    .np_reg_fit_total(ordinary, 24L, 3L, FALSE, REGTYPE_LP0),
    3L
  )
  expect_identical(
    .np_reg_fit_total(ordinary, 24L, 3L, TRUE, REGTYPE_LP),
    3L
  )
  expect_identical(
    .np_reg_fit_total(beta, 24L, 3L, TRUE, REGTYPE_LP0),
    27L
  )
  expect_identical(
    .np_reg_fit_total(adaptive, 24L, 3L, TRUE, REGTYPE_LP0),
    48L
  )
})
