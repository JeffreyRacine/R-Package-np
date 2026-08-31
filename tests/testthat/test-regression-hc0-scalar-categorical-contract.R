h4_explicit_lc_bw <- function(xdat, ydat, bws, bwtype = "fixed", ...) {
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

h4_actual_categorical_influence <- function(bws, txdat, exdat = NULL) {
  eval <- if (is.null(exdat)) txdat else exdat
  categorical <- which(!bws$icon)
  influence <- array(
    NA_real_,
    dim = c(nrow(eval), nrow(txdat), length(categorical))
  )

  for (observation in seq_len(nrow(txdat))) {
    unit <- numeric(nrow(txdat))
    unit[observation] <- 1
    args <- list(
      bws = bws,
      txdat = txdat,
      tydat = unit,
      gradients = TRUE,
      se = FALSE
    )
    if (!is.null(exdat))
      args$exdat <- exdat
    fit <- suppressWarnings(do.call(npreg, args))
    influence[, observation, ] <-
      fit$grad[, categorical, drop = FALSE]
  }
  influence
}

h4_expect_categorical_hc0 <- function(bws, xdat, ydat, exdat = NULL,
                                       tolerance = 2e-11) {
  args <- list(
    bws = bws,
    txdat = xdat,
    tydat = ydat,
    gradients = TRUE
  )
  if (!is.null(exdat))
    args$exdat <- exdat
  without.se <- suppressWarnings(do.call(npreg, c(args, list(se = FALSE))))
  with.se <- suppressWarnings(do.call(npreg, c(args, list(se = TRUE))))
  training.hat <- unclass(suppressWarnings(npreghat(
    bws = bws,
    txdat = xdat,
    output = "matrix"
  )))
  residual <- as.double(ydat) - drop(training.hat %*% as.double(ydat))
  influence <- h4_actual_categorical_influence(bws, xdat, exdat)
  expected <- vapply(
    seq_len(dim(influence)[3L]),
    function(coordinate) {
      sqrt(drop((influence[, , coordinate]^2) %*% (residual^2)))
    },
    numeric(dim(influence)[1L])
  )
  if (is.null(dim(expected)))
    expected <- matrix(expected, ncol = 1L)
  categorical <- which(!bws$icon)

  expect_identical(with.se$mean, without.se$mean)
  expect_identical(with.se$grad, without.se$grad)
  expect_identical(with.se$xtra, without.se$xtra)
  expect_equal(
    with.se$gerr[, categorical, drop = FALSE],
    expected,
    tolerance = tolerance
  )
  expect_true(all(is.finite(with.se$gerr[, categorical, drop = FALSE])))
  expect_true(all(with.se$gerr[, categorical, drop = FALSE] >= 0))
  invisible(with.se)
}

test_that("legacy scalar categorical HC0 matches the actual contrast map", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(
    x = seq(-1.1, 1.2, length.out = 24L),
    u = factor(rep(c("a", "b", "c"), length.out = 24L)),
    o = ordered(rep(1:3, each = 8L))
  )
  ydat <- sin(1.3 * xdat$x) + 0.18 * as.integer(xdat$u) -
    0.11 * as.integer(xdat$o) + seq_len(nrow(xdat)) / 120
  exdat <- xdat[c(1L, 4L, 9L, 16L, 24L), , drop = FALSE]

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- h4_explicit_lc_bw(
      xdat,
      ydat,
      if (identical(bwtype, "fixed")) c(0.42, 0.23, 0.19) else
        c(8, 0.23, 0.19),
      bwtype = bwtype
    )
    h4_expect_categorical_hc0(bw, xdat, ydat)
    h4_expect_categorical_hc0(bw, xdat, ydat, exdat)
  }
})

test_that("categorical-only scalar HC0 is finite under every bandwidth mode", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(
    u = factor(rep(c("a", "b", "c"), length.out = 21L)),
    o = ordered(rep(1:3, each = 7L))
  )
  ydat <- 0.3 * as.integer(xdat$u) - 0.17 * as.integer(xdat$o) +
    sin(seq_len(nrow(xdat)) / 4)
  exdat <- xdat[c(1L, 7L, 8L, 14L, 15L, 21L), , drop = FALSE]

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- h4_explicit_lc_bw(
      xdat, ydat, c(0.22, 0.18), bwtype = bwtype
    )
    h4_expect_categorical_hc0(bw, xdat, ydat)
    h4_expect_categorical_hc0(bw, xdat, ydat, exdat)
  }
})

test_that("beta scalar categorical HC0 uses paired normalized rows", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(
    x = c(0.03, 0.10, 0.18, 0.27, 0.39, 0.51, 0.64, 0.76, 0.87, 0.96),
    u = factor(rep(c("a", "b"), 5L)),
    o = ordered(rep(1:3, length.out = 10L))
  )
  ydat <- sin(2 * xdat$x) + 0.2 * as.integer(xdat$u) -
    0.12 * as.integer(xdat$o) + seq_len(nrow(xdat)) / 100
  exdat <- xdat[c(1L, 4L, 7L, 10L), , drop = FALSE]

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- h4_explicit_lc_bw(
      xdat,
      ydat,
      if (identical(bwtype, "fixed")) c(0.2, 0.2, 0.2) else
        c(4, 0.2, 0.2),
      bwtype = bwtype,
      ckertype = "beta",
      ckerorder = 4,
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1
    )
    h4_expect_categorical_hc0(bw, xdat, ydat)
    h4_expect_categorical_hc0(bw, xdat, ydat, exdat, tolerance = 4e-11)
  }
})

test_that("paired endpoint covariance reduces correctly at zero overlap", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(u = factor(rep(c("a", "b"), each = 8L)))
  ydat <- c(seq(-0.4, 0.3, length.out = 8L),
            seq(0.6, 1.5, length.out = 8L)) +
    c(0.03, -0.02, 0.01, -0.01, 0.02, -0.03, 0.01, 0,
      -0.02, 0.03, -0.01, 0.02, -0.03, 0.01, 0, 0.02)
  exdat <- data.frame(u = factor("b", levels = levels(xdat$u)))

  endpoint_rows <- function(lambda) {
    bw <- h4_explicit_lc_bw(xdat, ydat, lambda)
    training.hat <- unclass(npreghat(bws = bw, txdat = xdat,
                                     output = "matrix"))
    residual <- ydat - drop(training.hat %*% ydat)
    upper <- unclass(npreghat(
      bws = bw, txdat = xdat, exdat = exdat, output = "matrix"
    ))
    lower.data <- data.frame(u = factor("a", levels = levels(xdat$u)))
    lower <- unclass(npreghat(
      bws = bw, txdat = xdat, exdat = lower.data, output = "matrix"
    ))
    fit <- npreg(
      bws = bw, txdat = xdat, tydat = ydat, exdat = exdat,
      gradients = TRUE, se = TRUE
    )
    list(
      actual = fit$gerr[1L, 1L]^2,
      paired = sum(((upper - lower) * residual)^2),
      marginal = sum((upper * residual)^2) + sum((lower * residual)^2),
      cross = sum(upper * lower * residual^2)
    )
  }

  disjoint <- endpoint_rows(0)
  overlapping <- endpoint_rows(0.35)

  expect_equal(disjoint$cross, 0, tolerance = 1e-15)
  expect_equal(disjoint$paired, disjoint$marginal, tolerance = 1e-14)
  expect_equal(disjoint$actual, disjoint$paired, tolerance = 1e-14)
  expect_gt(abs(overlapping$cross), 1e-8)
  expect_equal(overlapping$actual, overlapping$paired, tolerance = 1e-14)
  expect_gt(abs(overlapping$marginal - overlapping$paired), 1e-8)
})

test_that("ordered first interior last and singleton conventions remain intact", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  xdat <- data.frame(o = ordered(rep(1:3, each = 6L)))
  ydat <- cos(seq_len(nrow(xdat)) / 3) + 0.2 * as.integer(xdat$o)
  exdat <- data.frame(o = ordered(1:3, levels = levels(xdat$o)))
  bw <- h4_explicit_lc_bw(xdat, ydat, 0.2)
  fit <- h4_expect_categorical_hc0(bw, xdat, ydat, exdat)

  expect_true(all(is.finite(fit$grad)))
  expect_true(all(is.finite(fit$gerr)))

  singleton.x <- data.frame(o = ordered(rep("only", 9L)))
  singleton.y <- seq_len(9L) / 7
  singleton.bw <- h4_explicit_lc_bw(singleton.x, singleton.y, 0.2)
  singleton.fit <- npreg(
    bws = singleton.bw,
    txdat = singleton.x,
    tydat = singleton.y,
    gradients = TRUE,
    se = TRUE
  )
  expect_identical(singleton.fit$grad, matrix(0, 9L, 1L))
  expect_identical(singleton.fit$gerr, matrix(0, 9L, 1L))
})

test_that("H4 retains only streamed categorical HC0 moments", {
  reducer <- paste(
    readLines(test_path("..", "..", "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  beta <- paste(
    readLines(
      test_path("..", "..", "src", "continuous_kernel_row.c"),
      warn = FALSE
    ),
    collapse = "\n"
  )

  expect_match(reducer, "level - alternate * ratio", fixed = TRUE)
  expect_match(
    beta,
    "np_continuous_kernel_beta_regression_paired_hc0_rows_validated",
    fixed = TRUE
  )
  expect_false(grepl("npreghat(", reducer, fixed = TRUE))
})
