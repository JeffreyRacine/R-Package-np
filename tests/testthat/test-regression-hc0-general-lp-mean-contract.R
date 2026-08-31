h5a_explicit_lp_bw <- function(xdat, ydat, bws, bwtype = "fixed",
                               degree = 2L, basis = "glp",
                               bernstein = FALSE, ...) {
  npregbw(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = FALSE,
    bwscaling = FALSE,
    bwtype = bwtype,
    regtype = "lp",
    degree = degree,
    degree.select = "manual",
    basis = basis,
    bernstein.basis = bernstein,
    ...
  )
}

h5a_hc0_mean_oracle <- function(bws, txdat, tydat, exdat = NULL) {
  training_hat <- unclass(suppressWarnings(npreghat(
    bws = bws,
    txdat = txdat,
    output = "matrix"
  )))
  evaluation_hat <- if (is.null(exdat)) {
    training_hat
  } else {
    unclass(suppressWarnings(npreghat(
      bws = bws,
      txdat = txdat,
      exdat = exdat,
      output = "matrix"
    )))
  }
  residual <- as.double(tydat) -
    drop(training_hat %*% as.double(tydat))

  sqrt(drop((evaluation_hat^2) %*% (residual^2)))
}

h5a_expect_mean_contract <- function(bws, txdat, tydat, exdat = NULL,
                                     tolerance = 2e-10) {
  args <- list(
    bws = bws,
    txdat = txdat,
    tydat = tydat,
    gradients = TRUE
  )
  if (!is.null(exdat))
    args$exdat <- exdat

  without_se <- do.call(npreg, c(args, list(se = FALSE)))
  with_se <- do.call(npreg, c(args, list(se = TRUE)))
  oracle <- h5a_hc0_mean_oracle(bws, txdat, tydat, exdat)

  expect_identical(with_se$mean, without_se$mean)
  expect_identical(with_se$grad, without_se$grad)
  expect_identical(with_se$xtra, without_se$xtra)
  expect_equal(with_se$merr, oracle, tolerance = tolerance)
  expect_true(all(is.finite(with_se$merr)))
  expect_true(all(with_se$merr >= 0))
  expect_true(all(is.na(with_se$gerr)))
  invisible(with_se)
}

test_that("ordinary general-LP mean SEs equal actual-map HC0 across bandwidth modes", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  n <- 31L
  txdat <- data.frame(
    x = seq(-1.15, 1.2, length.out = n),
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(1:4, length.out = n))
  )
  tydat <- sin(1.4 * txdat$x) +
    c(a = -0.2, b = 0.35, c = 0.1)[txdat$u] +
    0.08 * as.integer(txdat$o) + seq_len(n) / 110
  exdat <- txdat[c(2L, 7L, 14L, 23L, 29L), , drop = FALSE]

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bws <- if (identical(bwtype, "fixed")) {
      c(0.48, 0.2, 0.25)
    } else {
      c(9, 0.2, 0.25)
    }
    bw <- h5a_explicit_lp_bw(
      txdat, tydat, bws, bwtype = bwtype, degree = 2L,
      basis = "glp", bernstein = FALSE
    )

    h5a_expect_mean_contract(bw, txdat, tydat)
    h5a_expect_mean_contract(bw, txdat, tydat, exdat)

    mean_only <- npreg(
      bws = bw, txdat = txdat, tydat = tydat, exdat = exdat,
      gradients = FALSE, se = TRUE
    )
    expect_equal(
      mean_only$merr,
      h5a_hc0_mean_oracle(bw, txdat, tydat, exdat),
      tolerance = 2e-10
    )
  }
})

test_that("raw, shifted-Legendre, and Bernstein LP bases share one HC0 owner", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 29L
  txdat <- data.frame(
    x1 = seq(0.03, 0.97, length.out = n),
    x2 = 0.5 + 0.42 * sin(seq(0.1, 2.9, length.out = n))
  )
  tydat <- sin(2 * pi * txdat$x1) + 0.35 * txdat$x2 +
    seq_len(n) / 130
  exdat <- txdat[c(3L, 8L, 16L, 22L, 27L), , drop = FALSE]
  specifications <- list(
    list(basis = "glp", bernstein = FALSE),
    list(basis = "glp", bernstein = TRUE),
    list(basis = "tensor", bernstein = TRUE)
  )

  for (specification in specifications) {
    bw <- h5a_explicit_lp_bw(
      txdat, tydat, c(0.31, 0.34), degree = c(2L, 1L),
      basis = specification$basis,
      bernstein = specification$bernstein,
      ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
    )
    h5a_expect_mean_contract(bw, txdat, tydat, exdat, tolerance = 8e-10)
  }
})

test_that("general-LP HC0 follows the accepted ridge map and response equivariance", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  txdat <- data.frame(
    x = c(-1, -1, -0.35, -0.35, 0.2, 0.2, 0.75, 0.75, 1.1)
  )
  tydat <- c(0.2, 0.35, -0.1, 0.15, 0.8, 0.65, 1.4, 1.1, 1.7)
  exdat <- data.frame(x = c(-0.9, -0.25, 0.3, 0.9))
  bw <- suppressWarnings(h5a_explicit_lp_bw(
    txdat, tydat, 0.19, degree = 3L, basis = "glp",
    bernstein = FALSE, ckertype = "epanechnikov"
  ))

  base <- suppressWarnings(h5a_expect_mean_contract(
    bw, txdat, tydat, exdat, tolerance = 2e-8
  ))
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

  expect_equal(shifted$merr, base$merr, tolerance = 2e-7)
  expect_equal(scaled$merr, 3.25 * base$merr, tolerance = 2e-8)
  expect_identical(constant$merr, rep(0, nrow(exdat)))
  expect_true(all(is.na(shifted$gerr)))
  expect_true(all(is.na(scaled$gerr)))
  expect_true(all(is.na(constant$gerr)))
})
