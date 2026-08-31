h3_explicit_lc_bw <- function(xdat, ydat, bws, bwtype = "fixed", ...) {
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

h3_hc0_derivative_oracle <- function(bws, txdat, tydat, exdat = NULL) {
  H0 <- unclass(suppressWarnings(npreghat(
    bws = bws,
    txdat = txdat,
    output = "matrix"
  )))
  residual <- as.double(tydat) - drop(H0 %*% as.double(tydat))
  eval <- if (is.null(exdat)) txdat else exdat
  H.eval <- if (is.null(exdat)) {
    H0
  } else {
    unclass(suppressWarnings(npreghat(
      bws = bws,
      txdat = txdat,
      exdat = exdat,
      output = "matrix"
    )))
  }
  gradient <- matrix(NA_real_, nrow = nrow(eval), ncol = bws$ncon)
  stderr <- matrix(NA_real_, nrow = nrow(eval), ncol = bws$ncon)

  for (coordinate in seq_len(bws$ncon)) {
    derivative <- integer(bws$ncon)
    derivative[coordinate] <- 1L
    hat.args <- list(
      bws = bws,
      txdat = txdat,
      s = derivative,
      output = "matrix"
    )
    if (!is.null(exdat))
      hat.args$exdat <- exdat
    Hd <- unclass(suppressWarnings(do.call(npreghat, hat.args)))

    gradient[, coordinate] <- drop(Hd %*% as.double(tydat))
    stderr[, coordinate] <- sqrt(drop((Hd^2) %*% (residual^2)))
  }

  list(
    mean.stderr = sqrt(drop((H.eval^2) %*% (residual^2))),
    gradient = gradient,
    stderr = stderr
  )
}

h3_expect_fit_matches_oracle <- function(bws, xdat, ydat, exdat = NULL,
                                         tolerance = 8e-13) {
  args <- list(
    bws = bws,
    txdat = xdat,
    tydat = ydat,
    gradients = TRUE
  )
  if (!is.null(exdat))
    args$exdat <- exdat

  without.se <- do.call(npreg, c(args, list(se = FALSE)))
  with.se <- do.call(npreg, c(args, list(se = TRUE)))
  oracle <- h3_hc0_derivative_oracle(bws, xdat, ydat, exdat)

  expect_identical(with.se$mean, without.se$mean)
  expect_identical(with.se$grad, without.se$grad)
  expect_identical(with.se$xtra, without.se$xtra)
  expect_equal(with.se$merr, oracle$mean.stderr, tolerance = tolerance)
  expect_equal(
    with.se$grad[, bws$icon, drop = FALSE],
    oracle$gradient,
    tolerance = tolerance
  )
  expect_equal(
    with.se$gerr[, bws$icon, drop = FALSE],
    oracle$stderr,
    tolerance = tolerance
  )
  expect_true(all(is.finite(with.se$gerr[, bws$icon, drop = FALSE])))
  expect_true(all(with.se$gerr[, bws$icon, drop = FALSE] >= 0))
  if (any(!bws$icon))
    expect_true(all(is.na(with.se$gerr[, !bws$icon, drop = FALSE])))

  invisible(with.se)
}

test_that("legacy scalar derivative SEs equal explicit HC0 derivative hats", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(
    x = seq(-1.3, 1.4, length.out = 37L),
    z = factor(rep(c("low", "mid", "high"), length.out = 37L))
  )
  ydat <- sin(1.35 * xdat$x) +
    c(low = -0.25, mid = 0.2, high = 0.55)[xdat$z] +
    seq_len(nrow(xdat)) / 95
  exdat <- xdat[c(2L, 8L, 15L, 23L, 31L), , drop = FALSE]
  cases <- list(
    list(type = "fixed", kernel = "gaussian", bw = c(0.42, 0.17)),
    list(type = "generalized_nn", kernel = "epanechnikov", bw = c(9, 0.17)),
    list(type = "adaptive_nn", kernel = "uniform", bw = c(9, 0.17))
  )

  for (case in cases) {
    bw <- suppressWarnings(h3_explicit_lc_bw(
      xdat,
      ydat,
      case$bw,
      bwtype = case$type,
      ckertype = case$kernel
    ))
    h3_expect_fit_matches_oracle(bw, xdat, ydat)
    h3_expect_fit_matches_oracle(bw, xdat, ydat, exdat)
  }
})

test_that("bounded adaptive-NN derivative hats share the estimator kernel row", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(x = seq(0.02, 0.98, length.out = 33L))
  ydat <- sin(4 * xdat$x) + seq_len(nrow(xdat)) / 100
  exdat <- data.frame(x = c(0.03, 0.11, 0.42, 0.79, 0.96))
  bw <- h3_explicit_lc_bw(
    xdat,
    ydat,
    8,
    bwtype = "adaptive_nn",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  h3_expect_fit_matches_oracle(bw, xdat, ydat, exdat)
})

test_that("mixed multivariate scalar derivative SEs retain one HC0 row contract", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(
    x = seq(-1.1, 1.2, length.out = 39L),
    w = cos(seq(-0.7, 1.4, length.out = 39L)),
    u = factor(rep(letters[1:3], length.out = 39L)),
    o = ordered(rep(1:3, each = 13L))
  )
  ydat <- sin(1.2 * xdat$x) + 0.45 * xdat$w +
    0.12 * as.integer(xdat$u) - 0.08 * as.integer(xdat$o) +
    seq_len(nrow(xdat)) / 130
  exdat <- xdat[c(3L, 10L, 18L, 27L, 36L), , drop = FALSE]

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- h3_explicit_lc_bw(
      xdat,
      ydat,
      if (identical(bwtype, "fixed")) {
        c(0.46, 0.39, 0.18, 0.22)
      } else {
        c(10, 10, 0.18, 0.22)
      },
      bwtype = bwtype
    )
    h3_expect_fit_matches_oracle(bw, xdat, ydat)
    h3_expect_fit_matches_oracle(bw, xdat, ydat, exdat)
  }
})

test_that("beta scalar derivative SEs equal HC0 hats across orders and modes", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(x = c(0.03, 0.10, 0.24, 0.43, 0.67, 0.82, 0.96))
  ydat <- sin(2 * xdat$x) + xdat$x +
    c(-0.03, 0.02, 0.01, -0.02, 0.04, -0.01, 0.02)
  exdat <- data.frame(x = c(0.08, 0.38, 0.74, 0.92))

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bw <- h3_explicit_lc_bw(
        xdat,
        ydat,
        if (identical(bwtype, "fixed")) 0.18 else 3,
        bwtype = bwtype,
        ckertype = "beta",
        ckerorder = order,
        ckerbound = "fixed",
        ckerlb = 0,
        ckerub = 1
      )
      h3_expect_fit_matches_oracle(bw, xdat, ydat)
      h3_expect_fit_matches_oracle(bw, xdat, ydat, exdat, tolerance = 8e-12)
    }
  }
})

test_that("scalar derivative HC0 is response-equivariant and exact for constants", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(x = seq(-1, 1, length.out = 31L))
  ydat <- sin(xdat$x) + seq_len(nrow(xdat)) / 40
  bw <- h3_explicit_lc_bw(
    xdat, ydat, 8, bwtype = "adaptive_nn"
  )
  base <- npreg(
    bws = bw, txdat = xdat, tydat = ydat, gradients = TRUE, se = TRUE
  )
  shifted <- npreg(
    bws = bw, txdat = xdat, tydat = ydat + 1e6,
    gradients = TRUE, se = TRUE
  )
  scaled <- npreg(
    bws = bw, txdat = xdat, tydat = -3.25 * ydat,
    gradients = TRUE, se = TRUE
  )
  constant <- npreg(
    bws = bw, txdat = xdat, tydat = rep(5.25, nrow(xdat)),
    gradients = TRUE, se = TRUE
  )

  expect_equal(shifted$gerr, base$gerr, tolerance = 2e-9)
  expect_equal(scaled$gerr, 3.25 * base$gerr, tolerance = 3e-13)
  expect_identical(constant$gerr, matrix(0, nrow(xdat), 1L))
})

test_that("H3 keeps derivative HC0 streamed and scalar-only", {
  source <- paste(
    readLines(test_path("..", "..", "src", "np.c"), warn = FALSE),
    collapse = "\n"
  )
  reducer <- paste(
    readLines(test_path("..", "..", "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(source, "eg_hc0_storage = alloc_tmatd", fixed = TRUE)
  expect_match(
    reducer,
    "np_regression_hc0_derivative_moments_accumulate",
    fixed = TRUE
  )
  expect_match(
    reducer,
    "hc0_dual_power_ctx.regression_derivative = &hc0_derivative_ctx",
    fixed = TRUE
  )
  expect_false(grepl("npreghat", substring(
    source,
    regexpr("void np_regression(", source, fixed = TRUE)
  ), fixed = TRUE))
})
