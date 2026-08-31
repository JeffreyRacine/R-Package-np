h7b_explicit_bw <- function(x, y, bws, bwtype = "fixed",
                            basis = "glp", bernstein = FALSE,
                            ckertype = "gaussian", ckerorder = 2L) {
  beta <- identical(ckertype, "beta")
  npregbw(
    xdat = x, ydat = y, bws = bws,
    bandwidth.compute = FALSE, bwscaling = FALSE,
    bwtype = bwtype, regtype = "lp", degree = 2L,
    degree.select = "manual", basis = basis,
    bernstein.basis = bernstein, ckertype = ckertype,
    ckerorder = ckerorder,
    ckerbound = if (beta) "fixed" else "none",
    ckerlb = if (beta) 0 else NULL,
    ckerub = if (beta) 1 else NULL
  )
}

h7b_fit <- function(bw, x, y, eval = NULL, se = FALSE) {
  call <- list(
    bws = bw, txdat = x, tydat = y,
    gradients = TRUE, se = se, warn.glp.gradient = FALSE
  )
  if (!is.null(eval))
    call$exdat <- eval
  suppressWarnings(do.call(npreg, call))
}

h7b_unit_response_categorical_se <- function(bw, x, y, eval = NULL) {
  training <- h7b_fit(bw, x, y, se = FALSE)
  residual <- y - training$mean
  n.eval <- if (is.null(eval)) nrow(x) else nrow(eval)
  cat.index <- which(bw$iuno | bw$iord)
  out <- matrix(NA_real_, nrow = n.eval, ncol = length(cat.index))

  for (column in seq_along(cat.index)) {
    influence <- vapply(seq_len(nrow(x)), function(donor) {
      unit <- numeric(nrow(x))
      unit[donor] <- 1
      h7b_fit(bw, x, unit, eval = eval, se = FALSE)$grad[, cat.index[column]]
    }, numeric(n.eval))
    out[, column] <- sqrt(rowSums(sweep(influence, 2L, residual, `*`)^2))
  }
  out
}

test_that("positive-degree categorical HC0 SE matches literal influence rows", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 13L
  x <- data.frame(
    z = seq(0.04, 0.96, length.out = n),
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(c("low", "mid", "high"), length.out = n),
                levels = c("low", "mid", "high"))
  )
  y <- sin(2 * x$z) + 0.35 * (x$u == "b") -
    0.25 * (x$o == "high") + seq_len(n) / 100
  external <- x[c(1L, 4L, 8L, 13L), , drop = FALSE]
  cases <- list(
    list(type = "fixed", bws = c(0.31, 0.25, 0.3), eval = NULL),
    list(type = "generalized_nn", bws = c(7, 0.25, 0.3), eval = NULL),
    list(type = "adaptive_nn", bws = c(7, 0.25, 0.3), eval = external),
    list(type = "fixed", bws = c(0.27, 0.25, 0.3), eval = external,
         basis = "tensor", bernstein = TRUE,
         ckertype = "beta", ckerorder = 4L)
  )

  for (case in cases) {
    bw <- do.call(h7b_explicit_bw, c(
      list(x = x, y = y, bws = case$bws, bwtype = case$type),
      case[intersect(names(case),
                     c("basis", "bernstein", "ckertype", "ckerorder"))]
    ))
    no.se <- h7b_fit(bw, x, y, eval = case$eval, se = FALSE)
    native <- h7b_fit(bw, x, y, eval = case$eval, se = TRUE)
    oracle <- h7b_unit_response_categorical_se(
      bw, x, y, eval = case$eval
    )
    cat.index <- which(bw$iuno | bw$iord)

    expect_identical(native$mean, no.se$mean)
    expect_identical(native$grad, no.se$grad)
    expect_true(all(is.finite(native$gerr[, cat.index, drop = FALSE])))
    expect_true(all(native$gerr[, cat.index, drop = FALSE] >= 0))
    expect_equal(native$gerr[, cat.index, drop = FALSE], oracle,
                 tolerance = 2e-11)
  }
})

test_that("ordered boundaries and all-large invariant map remain typed", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  n <- 15L
  x <- data.frame(
    z = seq(-1, 1, length.out = n),
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(c("low", "mid", "high"), length.out = n),
                levels = c("low", "mid", "high"))
  )
  y <- sin(x$z) + 0.3 * (x$u == "b") -
    0.2 * (x$o == "high") + seq_len(n) / 100
  external <- x[c(1L, 2L, 8L, 14L, 15L), , drop = FALSE]

  ordinary <- h7b_explicit_bw(x, y, c(0.44, 0.25, 0.3))
  fit <- h7b_fit(ordinary, x, y, eval = external, se = TRUE)
  expect_true(all(is.finite(fit$gerr[, 2:3, drop = FALSE])))

  all.large <- h7b_explicit_bw(x, y, c(1e8, 2 / 3, 1))
  for (eval in list(NULL, external)) {
    fit <- h7b_fit(all.large, x, y, eval = eval, se = TRUE)
    expect_identical(fit$grad[, 2:3, drop = FALSE],
                     matrix(0, nrow(fit$grad), 2L))
    expect_identical(fit$gerr[, 2:3, drop = FALSE],
                     matrix(0, nrow(fit$gerr), 2L))
  }
})

test_that("H7B retains two bounded rows and removes estimator helper re-entry", {
  c.source <- paste(
    readLines(test_path("..", "..", "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  r.source <- paste(
    readLines(test_path("..", "..", "R", "np.regression.R"), warn = FALSE),
    collapse = "\n"
  )
  start <- regexpr(
    "if (gradients && !glp.gradient.partial)", r.source, fixed = TRUE
  )[[1L]]
  finish <- regexpr(
    "} else if (gradients) {", r.source, fixed = TRUE
  )[[1L]]
  active.block <- substr(r.source, start, finish - 1L)

  expect_gt(start, 0L)
  expect_gt(finish, start)
  expect_false(grepl("npreghat", active.block, fixed = TRUE))
  expect_match(c.source, "categorical_base_kernel_row", fixed = TRUE)
  expect_match(c.source, "categorical_alternate_kernel_row", fixed = TRUE)
  expect_match(c.source,
               "np_regression_hc0_lp_categorical_standard_error",
               fixed = TRUE)
  c.start <- regexpr(
    "static int np_regression_hc0_lp_categorical_standard_error",
    c.source, fixed = TRUE
  )[[1L]]
  c.tail <- substr(c.source, c.start, nchar(c.source))
  c.finish <- regexpr("typedef struct {", c.tail, fixed = TRUE)[[1L]]
  expect_gt(c.start, 0L)
  expect_gt(c.finish, 0L)
  expect_false(grepl(
    "alloc_matd", substr(c.tail, 1L, c.finish - 1L), fixed = TRUE
  ))
})
