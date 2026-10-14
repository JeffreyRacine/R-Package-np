h6_extract_lp_owner <- function(source) {
  start <- gregexpr(
    "static NP_ALWAYS_INLINE SEXP np_regression_general_lp_fit_body(",
    source, fixed = TRUE
  )[[1L]]
  finish <- gregexpr(
    "static int np_regression_general_lp_fit(\n",
    source, fixed = TRUE
  )[[1L]]
  if (length(start) != 1L || start[[1L]] < 1L ||
      length(finish) != 1L || finish[[1L]] < 1L)
    stop("H6 owner requires unique start and finish anchors", call. = FALSE)
  if (finish[[1L]] <= start[[1L]])
    stop("H6 owner anchors are out of order", call. = FALSE)
  substr(source, start[[1L]], finish[[1L]] - 1L)
}

test_that("general-LP continuous derivative HC0 matches actual hats across bandwidth modes", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  n <- 31L
  txdat <- data.frame(
    x1 = seq(-1.1, 1.15, length.out = n),
    x2 = 0.55 * cos(seq(-0.6, 2.4, length.out = n)),
    u = factor(rep(c("a", "b", "c"), length.out = n))
  )
  tydat <- sin(1.3 * txdat$x1) + 0.45 * txdat$x2 +
    c(a = -0.2, b = 0.15, c = 0.4)[txdat$u] + seq_len(n) / 130
  exdat <- txdat[c(2L, 8L, 15L, 23L, 30L), , drop = FALSE]

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bandwidth <- if (identical(bwtype, "fixed")) {
      c(0.46, 0.39, 0.2)
    } else {
      c(9, 9, 0.2)
    }
    bw <- h6_explicit_lp_bw(
      txdat, tydat, bandwidth, bwtype = bwtype,
      degree = c(2L, 1L), basis = "glp", bernstein = FALSE
    )
    h6_expect_lp_derivative_contract(bw, txdat, tydat)
    h6_expect_lp_derivative_contract(bw, txdat, tydat, exdat)
  }
})

test_that("near-interpolating beta ANN retains finite requested inference", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  n <- 15L
  x <- data.frame(x1 = seq(0.04, 0.96, length.out = n),
    x2 = 0.5 + 0.43 * sin(seq(0.2, 2.8, length.out = n)))
  y <- sin(4 * x$x1) + 0.35 * x$x2 + seq_len(n) / 140
  ex <- x[c(2L, 5L, 8L, 12L, 15L), , drop = FALSE]
  bw <- h6_explicit_lp_bw(x, y, c(7, 7), bwtype = "adaptive_nn",
    degree = c(2L, 1L), basis = "glp", bernstein = FALSE,
    ckertype = "beta", ckerorder = 2L, ckerbound = "fixed",
    ckerlb = c(0, 0), ckerub = c(1, 1))
  point <- npreg(bws = bw, txdat = x, tydat = y, exdat = ex,
                 gradients = TRUE, se = FALSE)
  inference <- npreg(bws = bw, txdat = x, tydat = y, exdat = ex,
                     gradients = TRUE, se = TRUE)
  expect_identical(inference$mean, point$mean)
  expect_identical(inference$grad, point$grad)
  expect_identical(inference$xtra, point$xtra)
  expect_true(all(is.finite(inference$merr)))
  expect_true(all(is.finite(inference$gerr)))
  expect_true(all(inference$merr >= 0))
  expect_true(all(inference$gerr >= 0))
  # Exact retained adjoint/basis/kernel inputs for this original fixture are
  # separately checked by the campaign's 90-digit donor-covariance oracle.
  # Do not replace that proof with a hard-coded platform-dependent SE or a
  # relaxed tolerance against a response-perturbed public forward map.
})

test_that("general-LP derivative HC0 follows accepted ridge and all-large maps", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  ridge.x <- data.frame(
    x = c(-1, -1, -0.35, -0.35, 0.2, 0.2, 0.75, 0.75, 1.1)
  )
  ridge.y <- c(0.2, 0.35, -0.1, 0.15, 0.8, 0.65, 1.4, 1.1, 1.7)
  ridge.eval <- data.frame(x = c(-0.9, -0.25, 0.3, 0.9))
  ridge.bw <- suppressWarnings(h6_explicit_lp_bw(
    ridge.x, ridge.y, 0.19, degree = 3L,
    ckertype = "epanechnikov"
  ))
  h6_expect_lp_derivative_contract(
    ridge.bw, ridge.x, ridge.y, ridge.eval, tolerance = 3e-8,
    oracle.method = "actual"
  )

  large.x <- data.frame(x = seq(-1, 1, length.out = 21L))
  large.y <- sin(large.x$x) + seq_len(nrow(large.x)) / 90
  large.eval <- data.frame(x = c(-0.8, -0.15, 0.45, 0.9))
  large.bw <- h6_explicit_lp_bw(
    large.x, large.y, 1e16, degree = 2L
  )
  h6_expect_lp_derivative_contract(
    large.bw, large.x, large.y, large.eval, tolerance = 2e-9,
    oracle.method = "actual"
  )
})

test_that("balanced q-adjusted multiplier covariance equals general-LP HC0", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 16L
  txdat <- data.frame(x = seq(-1, 1, length.out = n))
  tydat <- sin(1.7 * txdat$x) + seq_len(n) / 60
  exdat <- data.frame(x = c(-0.75, -0.2, 0.35, 0.8))
  bw <- h6_explicit_lp_bw(txdat, tydat, 0.48, degree = 2L)
  result <- h6_expect_lp_derivative_contract(bw, txdat, tydat, exdat)

  hadamard <- matrix(1, 1L, 1L)
  while (nrow(hadamard) < n)
    hadamard <- rbind(
      cbind(hadamard, hadamard),
      cbind(hadamard, -hadamard)
    )
  multipliers <- rbind(hadamard, -hadamard)
  # This is an algebraic covariance identity using q-adjusted donors, not a
  # claim that the separate public wild-bootstrap residual law was changed.
  derivative.donor <- sweep(
    result$oracle$derivative.hat[[1L]],
    2L,
    result$oracle$residual,
    `*`
  )
  bootstrap.deviation <- derivative.donor %*% t(multipliers)
  bootstrap.stderr <- sqrt(rowMeans(bootstrap.deviation^2))

  expect_equal(
    result$fit$gerr[, bw$icon], bootstrap.stderr,
    tolerance = 2e-12
  )
})

test_that("H6 batches adjoint directions without another covariance owner", {
  source <- paste(
    readLines(test_path("..", "..", "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  owner <- h6_extract_lp_owner(source)

  expect_match(
    owner,
    "if(call->do_grad && call->do_gerr)",
    fixed = TRUE
  )
  expect_match(
    owner,
    "np_lp_solve_workspace_solve_adjoint_factored(",
    fixed = TRUE
  )
  expect_match(
    owner,
    "np_regression_hc0_lp_standard_error_reuse(",
    fixed = TRUE
  )
  # Exact-input reuse is only a computational wrapper. No-cache and miss
  # paths must still use the same donor-square calculation and propagate its
  # failures before publishing a reusable result.
  source.lines <- strsplit(source, "\n", fixed = TRUE)[[1L]]
  reuse <- npRmpi_test_extract_c_function(
    source.lines, "np_regression_hc0_lp_standard_error_reuse")
  donor <- gregexpr("np_regression_hc0_lp_standard_error_with_information(",
                   reuse, fixed = TRUE)[[1L]]
  expect_length(donor[donor > 0L], 2L)
  compact.reuse <- gsub("[[:space:]]+", " ", reuse)
  expect_match(compact.reuse,
    paste0("if(reuse->entry == NULL) return ",
           "np_regression_hc0_lp_standard_error_with_information("),
    fixed = TRUE)
  expect_match(compact.reuse,
    "standard_error, &missing)) return 0; np_inference_reuse_put(",
    fixed = TRUE)
  expect_match(reuse, "if(!np_inference_reuse_get(", fixed = TRUE)
  expect_match(reuse, "if(unavailable != NULL && missing) *unavailable = 1;",
               fixed = TRUE)
  expect_lt(regexpr("np_inference_reuse_get(", reuse, fixed = TRUE)[[1L]],
            donor[[2L]])
  expect_lt(donor[[2L]],
            regexpr("np_inference_reuse_put(", reuse, fixed = TRUE)[[1L]])
  expect_match(gsub("[[:space:]]+", " ", owner),
    paste0("if(ordinary_hc0 && call->do_merr && call->all_large_kernel_rows) ",
           "(void)np_inference_reuse_reserve("), fixed = TRUE)
  cleanup <- npRmpi_test_extract_c_function(
    source.lines, "np_regression_general_lp_fit_owner_cleanup")
  expect_match(cleanup, "np_inference_reuse_clear(&owner->inference_reuse);",
               fixed = TRUE)
  expect_false(grepl(
    "!ordinary_hc0 && call->do_grad && call->do_gerr",
    owner,
    fixed = TRUE
  ))
  expect_false(grepl("npreghat", owner, fixed = TRUE))
  expect_false(grepl("inverse", owner, fixed = TRUE))
})
