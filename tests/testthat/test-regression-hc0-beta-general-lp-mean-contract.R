h5b_beta_lp_bw <- function(xdat, ydat, bws, bwtype = "fixed",
                           order = 2L, degree = 2L,
                           basis = "glp", bernstein = FALSE) {
  npregbw(
    xdat = xdat, ydat = ydat, bws = bws,
    bandwidth.compute = FALSE, bwscaling = FALSE,
    bwtype = bwtype, regtype = "lp", degree = degree,
    degree.select = "manual", basis = basis,
    bernstein.basis = bernstein,
    ckertype = "beta", ckerorder = order,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
}

h5b_beta_hc0_oracle <- function(bws, txdat, tydat, exdat = NULL) {
  training_hat <- unclass(suppressWarnings(npreghat(
    bws = bws, txdat = txdat, output = "matrix"
  )))
  evaluation_hat <- if (is.null(exdat)) {
    training_hat
  } else {
    unclass(suppressWarnings(npreghat(
      bws = bws, txdat = txdat, exdat = exdat, output = "matrix"
    )))
  }
  residual <- as.double(tydat) -
    drop(training_hat %*% as.double(tydat))

  sqrt(drop((evaluation_hat^2) %*% (residual^2)))
}

h5b_expect_beta_mean_contract <- function(bws, txdat, tydat, exdat = NULL,
                                           tolerance = 3e-9) {
  args <- list(
    bws = bws, txdat = txdat, tydat = tydat,
    gradients = TRUE
  )
  if (!is.null(exdat))
    args$exdat <- exdat

  without_se <- suppressWarnings(do.call(npreg, c(args, list(se = FALSE))))
  with_se <- suppressWarnings(do.call(npreg, c(args, list(se = TRUE))))
  oracle <- h5b_beta_hc0_oracle(bws, txdat, tydat, exdat)

  expect_identical(with_se$mean, without_se$mean)
  expect_identical(with_se$grad, without_se$grad)
  expect_identical(with_se$xtra, without_se$xtra)
  expect_equal(with_se$merr, oracle, tolerance = tolerance)
  expect_true(all(is.finite(with_se$merr)))
  expect_true(all(with_se$merr >= 0))
  expect_true(all(is.finite(with_se$gerr[, bws$icon, drop = FALSE])))
  expect_true(all(with_se$gerr[, bws$icon, drop = FALSE] >= 0))
  categorical <- setdiff(seq_len(ncol(with_se$gerr)), which(bws$icon))
  if (length(categorical)) {
    expect_true(all(is.finite(with_se$gerr[, categorical, drop = FALSE])))
    expect_true(all(with_se$gerr[, categorical, drop = FALSE] >= 0))
  }
  invisible(with_se)
}

test_that("beta general-LP HC0 spans orders, NN topologies, and mixed predictors", {
  old <- options(
    np.messages = FALSE, np.tree = TRUE,
    np.categorical.compress = FALSE
  )
  on.exit(options(old), add = TRUE)

  n <- 31L
  txdat <- data.frame(
    x = seq(0.025, 0.975, length.out = n),
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(1:4, length.out = n))
  )
  tydat <- sin(4.1 * txdat$x) +
    c(a = -0.18, b = 0.24, c = 0.07)[txdat$u] +
    0.06 * as.integer(txdat$o) + seq_len(n) / 170
  exdat <- txdat[c(2L, 8L, 15L, 23L, 30L), , drop = FALSE]

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bandwidth <- if (identical(bwtype, "fixed")) {
        c(0.23, 0.2, 0.25)
      } else {
        c(9, 0.2, 0.25)
      }
      bw <- h5b_beta_lp_bw(
        txdat, tydat, bandwidth, bwtype = bwtype,
        order = order, degree = 2L,
        bernstein = order %in% c(4L, 8L)
      )
      options(np.categorical.compress = (order + match(
        bwtype, c("fixed", "generalized_nn", "adaptive_nn")
      )) %% 2L == 0L)
      h5b_expect_beta_mean_contract(bw, txdat, tydat, exdat)
    }
  }

  bw <- h5b_beta_lp_bw(
    txdat, tydat, c(9, 0.2, 0.25), bwtype = "adaptive_nn",
    order = 8L, degree = 2L, bernstein = TRUE
  )
  h5b_expect_beta_mean_contract(bw, txdat, tydat)
})

test_that("beta HC0 mean meat is invariant to basis owner and compression", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 29L
  txdat <- data.frame(
    x1 = seq(0.035, 0.965, length.out = n),
    x2 = 0.5 + 0.43 * sin(seq(0.15, 2.95, length.out = n)),
    u = factor(rep(c("a", "b", "c"), length.out = n))
  )
  tydat <- sin(2 * pi * txdat$x1) + 0.31 * txdat$x2 +
    c(a = -0.1, b = 0.16, c = 0.04)[txdat$u] + seq_len(n) / 190
  exdat <- txdat[c(3L, 9L, 16L, 22L, 28L), , drop = FALSE]
  specifications <- list(
    list(basis = "glp", bernstein = FALSE, order = 2L),
    list(basis = "glp", bernstein = TRUE, order = 6L),
    list(basis = "tensor", bernstein = TRUE, order = 8L)
  )

  for (specification in specifications) {
    bw <- h5b_beta_lp_bw(
      txdat, tydat, c(0.27, 0.3, 0.2),
      order = specification$order, degree = c(2L, 1L),
      basis = specification$basis,
      bernstein = specification$bernstein
    )
    options(np.categorical.compress = FALSE)
    dense <- h5b_expect_beta_mean_contract(
      bw, txdat, tydat, exdat, tolerance = 8e-9
    )
    options(np.categorical.compress = TRUE)
    compressed <- h5b_expect_beta_mean_contract(
      bw, txdat, tydat, exdat, tolerance = 8e-9
    )

    expect_identical(compressed$mean, dense$mean)
    expect_identical(compressed$merr, dense$merr)
    expect_identical(compressed$grad, dense$grad)
    expect_identical(compressed$gerr, dense$gerr)
  }
})

test_that("beta general-LP HC0 follows the accepted map under response transforms", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  txdat <- data.frame(x = seq(0.03, 0.97, length.out = 27L))
  tydat <- sin(5 * txdat$x) + txdat$x^2 + seq_len(nrow(txdat)) / 210
  exdat <- data.frame(x = c(0.04, 0.19, 0.43, 0.72, 0.96))
  bw <- h5b_beta_lp_bw(
    txdat, tydat, 0.16, order = 8L, degree = 3L,
    basis = "glp", bernstein = FALSE
  )

  base <- h5b_expect_beta_mean_contract(
    bw, txdat, tydat, exdat, tolerance = 2e-8
  )
  shifted <- suppressWarnings(npreg(
    bws = bw, txdat = txdat, tydat = tydat + 1e7,
    exdat = exdat, gradients = TRUE, se = TRUE
  ))
  scaled <- suppressWarnings(npreg(
    bws = bw, txdat = txdat, tydat = -3.25 * tydat,
    exdat = exdat, gradients = TRUE, se = TRUE
  ))
  constant <- suppressWarnings(npreg(
    bws = bw, txdat = txdat, tydat = rep(7, nrow(txdat)),
    exdat = exdat, gradients = TRUE, se = TRUE
  ))

  expect_lt(max(abs(shifted$merr - base$merr)), 2e-8)
  expect_equal(scaled$merr, 3.25 * base$merr, tolerance = 3e-8)
  expect_identical(constant$merr, rep(0, nrow(exdat)))
  expect_lt(max(abs(shifted$gerr - base$gerr)), 1e-5)
  expect_equal(scaled$gerr, 3.25 * base$gerr, tolerance = 3e-8)
  expect_identical(constant$gerr, matrix(0, nrow(exdat), 1L))
})

test_that("beta general-LP HC0 survives complete absolute-weight underflow", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(6219)
  n <- 100L
  p <- 32L
  log_complement <- matrix(runif(n * p, -1.4, -0.8), nrow = n)
  log_complement <- log_complement - rowMeans(log_complement) - 1.1
  txdat <- as.data.frame(1 - exp(log_complement))
  names(txdat) <- paste0("x", seq_len(p))
  coefficient <- seq_len(p) / p
  tydat <- 1 + drop(as.matrix(txdat) %*% coefficient)
  exdat <- as.data.frame(matrix(0, nrow = 1L, ncol = p))
  names(exdat) <- names(txdat)
  bw <- h5b_beta_lp_bw(
    txdat, tydat, rep(0.1, p), order = 8L,
    degree = rep(1L, p), basis = "glp", bernstein = FALSE
  )
  raw <- npksum(
    bws = bw, txdat = txdat, exdat = exdat,
    return.kernel.weights = TRUE
  )
  without_se <- suppressWarnings(npreg(
    bws = bw, txdat = txdat, tydat = tydat, exdat = exdat,
    gradients = TRUE, se = FALSE
  ))
  with_se <- suppressWarnings(npreg(
    bws = bw, txdat = txdat, tydat = tydat, exdat = exdat,
    gradients = TRUE, se = TRUE
  ))
  oracle <- h5b_beta_hc0_oracle(bw, txdat, tydat, exdat)

  expect_true(all(raw$kw == 0))
  expect_identical(with_se$mean, without_se$mean)
  expect_identical(with_se$grad, without_se$grad)
  expect_lt(max(abs(with_se$merr - oracle)), 3e-14)
  expect_true(all(is.finite(with_se$gerr)))
  expect_true(all(with_se$gerr >= 0))
})

test_that("all-large general-LP point shortcut cedes only HC0 uncertainty", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  n <- 29L
  txdat <- data.frame(x = seq(-1, 1, length.out = n))
  tydat <- sin(2 * txdat$x) + seq_len(n) / 100
  exdat <- data.frame(x = c(-0.8, -0.25, 0.1, 0.55, 0.9))
  bw <- npregbw(
    xdat = txdat, ydat = tydat, bws = 1e8,
    bandwidth.compute = FALSE, bwscaling = FALSE,
    bwtype = "fixed", regtype = "lp", degree = 2L,
    degree.select = "manual", basis = "glp",
    bernstein.basis = FALSE, ckertype = "gaussian", ckerorder = 2L
  )

  h5b_expect_beta_mean_contract(
    bw, txdat, tydat, exdat, tolerance = 1e-12
  )
})
