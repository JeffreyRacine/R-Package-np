h7a_explicit_lp_bw <- function(txdat,
                               tydat,
                               bws,
                               bwtype = "fixed",
                               basis = "glp",
                               bernstein = FALSE,
                               ckertype = "gaussian",
                               ckerorder = 2L) {
  beta <- identical(ckertype, "beta")
  npregbw(
    xdat = txdat,
    ydat = tydat,
    bws = bws,
    bandwidth.compute = FALSE,
    bwscaling = FALSE,
    bwtype = bwtype,
    regtype = "lp",
    degree = 2L,
    degree.select = "manual",
    basis = basis,
    bernstein.basis = bernstein,
    ckertype = ckertype,
    ckerorder = ckerorder,
    ckerbound = if (beta) "fixed" else "none",
    ckerlb = if (beta) 0 else NULL,
    ckerub = if (beta) 1 else NULL
  )
}

h7a_fit <- function(bw, txdat, tydat, exdat = NULL, se = FALSE) {
  args <- list(
    bws = bw,
    txdat = txdat,
    tydat = tydat,
    gradients = TRUE,
    se = se,
    warn.glp.gradient = FALSE
  )
  if (!is.null(exdat))
    args$exdat <- exdat
  do.call(npreg, args)
}

h7a_native_fit <- function(bw, txdat, tydat, exdat = NULL, se = FALSE) {
  h7a_fit(bw, txdat, tydat, exdat = exdat, se = se)
}

h7a_expect_native_point <- function(bw,
                                    txdat,
                                    tydat,
                                    exdat = NULL,
                                    se = FALSE,
                                    tolerance = 5e-12,
                                    exact = FALSE) {
  native <- h7a_native_fit(bw, txdat, tydat, exdat = exdat, se = se)
  oracle <- native
  oracle$grad <- getFromNamespace(
    ".npreg_glp_categorical_gradients_from_npreghat", "np"
  )(
    bws = bw,
    txdat = txdat,
    tydat = tydat,
    exdat = exdat,
    grad = oracle$grad,
    where = "H7A oracle"
  )
  cat.index <- which(bw$iuno | bw$iord)

  expect_identical(native$mean, oracle$mean)
  expect_identical(native$grad[, bw$icon, drop = FALSE],
                   oracle$grad[, bw$icon, drop = FALSE])
  if (exact) {
    expect_identical(native$grad[, cat.index, drop = FALSE],
                     oracle$grad[, cat.index, drop = FALSE])
  } else {
    expect_equal(native$grad[, cat.index, drop = FALSE],
                 oracle$grad[, cat.index, drop = FALSE],
                 tolerance = tolerance)
  }
  invisible(list(oracle = oracle, native = native))
}

test_that("native general-LP categorical points match fixed/GNN/ANN oracles", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 36L
  txdat <- data.frame(
    x = seq(-1, 1, length.out = n),
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(
      rep(c("low", "mid", "high"), each = 12L),
      levels = c("low", "mid", "high")
    )
  )
  tydat <- sin(2 * txdat$x) + 0.4 * (txdat$u == "b") -
    0.3 * (txdat$o == "high") + seq_len(n) / 200
  exdat <- txdat[c(1L, 2L, 13L, 25L, 36L), , drop = FALSE]

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bws <- if (identical(bwtype, "fixed")) c(0.42, 0.3, 0.3) else
      c(9, 0.3, 0.3)
    bw <- h7a_explicit_lp_bw(txdat, tydat, bws, bwtype = bwtype)

    exact <- !identical(bwtype, "adaptive_nn")
    h7a_expect_native_point(bw, txdat, tydat, exact = exact)
    h7a_expect_native_point(bw, txdat, tydat, se = TRUE, exact = exact)
    h7a_expect_native_point(
      bw, txdat, tydat, exdat = exdat, exact = exact
    )
    h7a_expect_native_point(
      bw, txdat, tydat, exdat = exdat, se = TRUE, exact = exact
    )
  }
})

test_that("native categorical frames cover bases kernels and support edges", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  n <- 30L
  txdat <- data.frame(
    x = seq(0.03, 0.97, length.out = n),
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(
      rep(c("low", "mid", "high"), length.out = n),
      levels = c("low", "mid", "high")
    ),
    singleton = factor(rep("only", n))
  )
  tydat <- sin(2 * txdat$x) + 0.4 * (txdat$u == "b") -
    0.3 * (txdat$o == "high") + seq_len(n) / 200
  exdat <- txdat[c(1L, 2L, 3L, 14L, 29L), , drop = FALSE]
  cases <- list(
    list(bws = c(0.28, 0.3, 0.3, 0), basis = "glp",
         bernstein = FALSE, ckertype = "gaussian", ckerorder = 2L),
    list(bws = c(0.4, 0.3, 0.3, 0), basis = "additive",
         bernstein = FALSE, ckertype = "epanechnikov", ckerorder = 4L),
    list(bws = c(0.16, 0.3, 0.3, 0), basis = "tensor",
         bernstein = TRUE, ckertype = "beta", ckerorder = 4L)
  )

  for (case in cases) {
    bw <- h7a_explicit_lp_bw(
      txdat,
      tydat,
      case$bws,
      basis = case$basis,
      bernstein = case$bernstein,
      ckertype = case$ckertype,
      ckerorder = case$ckerorder
    )
    result <- h7a_expect_native_point(
      bw, txdat, tydat, exdat = exdat, se = TRUE
    )
    expect_identical(result$native$grad[, 4L], rep(0, nrow(exdat)))
  }
})

test_that("all-large native categorical point retains the accepted point map", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  n <- 29L
  txdat <- data.frame(
    x = seq(-1, 1, length.out = n),
    u = factor(rep(c("a", "b", "c"), length.out = n))
  )
  tydat <- sin(2 * txdat$x) + 0.3 * (txdat$u == "b") + seq_len(n) / 100
  exdat <- txdat[c(1L, 5L, 13L, 21L, 29L), , drop = FALSE]
  bw <- h7a_explicit_lp_bw(txdat, tydat, c(1e8, 2 / 3))
  result <- h7a_expect_native_point(
    bw, txdat, tydat, exdat = exdat, tolerance = 2e-12
  )

  expect_identical(result$native$grad[, 2L], rep(0, nrow(exdat)))
})

test_that("H7A remains a point-only bounded-work adapter", {
  source <- paste(
    readLines(test_path("..", "..", "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  start <- regexpr(
    "static int np_regression_general_lp_point_at_frame(",
    source,
    fixed = TRUE
  )[[1L]]
  finish <- regexpr(
    "static SEXP np_regression_general_lp_fit_execute(void *data)",
    source,
    fixed = TRUE
  )[[1L]]
  adapter <- substr(source, start, finish - 1L)

  expect_gt(start, 0L)
  expect_gt(finish, start)
  expect_false(grepl("npreghat[[:space:]]*\\(", adapter))
  expect_false(grepl("alloc_matd", adapter, fixed = TRUE))
  expect_false(grepl("power2_moments", adapter, fixed = TRUE))
  expect_match(source,
               "categorical_matrix_bandwidth = matrix_bandwidth_deriv",
               fixed = TRUE)
})
