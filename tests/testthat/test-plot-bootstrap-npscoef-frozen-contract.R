library(np)

quiet_eval <- function(expr) {
  value <- NULL
  capture.output(value <- force(expr))
  value
}

test_that("npscoef exact helper matches duplicate-row oracle and frozen stays finite", {
  boot.fun <- getFromNamespace(".np_inid_boot_from_scoef", "np")

  xdat <- data.frame(x = c(0.05, 0.05, 0.25, 0.25, 0.60, 0.60, 0.90, 0.90))
  zdat <- data.frame(z = c(0.10, 0.10, 0.35, 0.35, 0.65, 0.65, 0.95, 0.95))
  y <- with(xdat, 1 + x * (1 + zdat$z))
  exdat <- data.frame(x = c(0.10, 0.30, 0.70))
  ezdat <- data.frame(z = c(0.12, 0.40, 0.85))
  counts <- cbind(
    c(2, 0, 2, 0, 1, 1, 1, 1),
    c(0, 2, 0, 2, 1, 1, 1, 1),
    c(1, 1, 3, 0, 0, 2, 1, 0)
  )

  bw <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    bws = 3L,
    bwtype = "adaptive_nn",
    bandwidth.compute = FALSE,
    regtype = "ll"
  )

  exact.out <- boot.fun(
    txdat = xdat,
    ydat = y,
    tzdat = zdat,
    exdat = exdat,
    ezdat = ezdat,
    bws = bw,
    B = ncol(counts),
    counts = counts,
    mode = "exact"
  )
  frozen.out <- boot.fun(
    txdat = xdat,
    ydat = y,
    tzdat = zdat,
    exdat = exdat,
    ezdat = ezdat,
    bws = bw,
    B = ncol(counts),
    counts = counts,
    mode = "frozen"
  )

  exact.oracle <- vapply(seq_len(ncol(counts)), function(j) {
    idx <- np:::.np_counts_to_indices(counts[, j])
    as.vector(npscoef(
      bws = bw,
      txdat = xdat[idx, , drop = FALSE],
      tzdat = zdat[idx, , drop = FALSE],
      tydat = y[idx],
      exdat = exdat,
      ezdat = ezdat,
      iterate = FALSE,
      errors = FALSE
    )$mean)
  }, numeric(nrow(exdat)))
  exact.oracle <- t(exact.oracle)

  t0.oracle <- as.vector(npscoef(
    bws = bw,
    txdat = xdat,
    tzdat = zdat,
    tydat = y,
    exdat = exdat,
    ezdat = ezdat,
    iterate = FALSE,
    errors = FALSE
  )$mean)

  expect_equal(exact.out$t, exact.oracle, tolerance = 1e-10)
  expect_equal(as.vector(exact.out$t0), t0.oracle, tolerance = 1e-10)
  expect_equal(as.vector(frozen.out$t0), t0.oracle, tolerance = 1e-10)
  expect_true(all(is.finite(frozen.out$t)))
  expect_true(max(abs(frozen.out$t)) < 1e4)
})

test_that("npscoef frozen surface plot mode is forwarded", {
  xdat <- data.frame(x = c(0.05, 0.05, 0.25, 0.25, 0.60, 0.60, 0.90, 0.90))
  zdat <- data.frame(z = c(0.10, 0.10, 0.35, 0.35, 0.65, 0.65, 0.95, 0.95))
  y <- with(xdat, 1 + x * (1 + zdat$z))

  bw <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    bws = 3L,
    bwtype = "adaptive_nn",
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  ns <- asNamespace("np")
  orig <- getFromNamespace(".np_inid_boot_from_scoef", "np")
  modes <- character()
  assignInNamespace(".np_inid_boot_from_scoef", function(..., mode = c("exact", "frozen")) {
    mode <- match.arg(mode)
    modes <<- c(modes, mode)
    orig(..., mode = mode)
  }, ns = "np")
  on.exit(assignInNamespace(".np_inid_boot_from_scoef", orig, ns = "np"), add = TRUE)

  expect_no_error(capture.output(plot(
    bw,
    xdat = xdat,
    ydat = y,
    zdat = zdat,
    neval = 6L,
    coef = FALSE,
    plot.behavior = "data",
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.nonfixed = "frozen",
    plot.errors.boot.num = 41L,
    plot.errors.type = "pointwise"
  )))

  expect_gte(length(modes), 1L)
  expect_true(all(modes == "frozen"))
})

test_that("npscoef fixed helper treats exact and frozen identically", {
  boot.fun <- getFromNamespace(".np_inid_boot_from_scoef", "np")

  set.seed(42)
  n <- 30L
  xdat <- data.frame(x = runif(n, -1, 1))
  zdat <- data.frame(z = rnorm(n))
  y <- with(xdat, x^2 + rnorm(n, sd = 0.1))
  exdat <- data.frame(x = seq(-0.9, 0.9, length.out = 9L))
  ezdat <- data.frame(z = seq(-1.0, 1.0, length.out = 9L))
  counts <- rmultinom(n = 5L, size = n, prob = rep.int(1 / n, n))

  bw <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    bws = c(0.6),
    bwtype = "fixed",
    bandwidth.compute = FALSE,
    regtype = "ll"
  )

  exact.out <- boot.fun(
    txdat = xdat,
    ydat = y,
    tzdat = zdat,
    exdat = exdat,
    ezdat = ezdat,
    bws = bw,
    B = ncol(counts),
    counts = counts,
    mode = "exact"
  )
  frozen.out <- boot.fun(
    txdat = xdat,
    ydat = y,
    tzdat = zdat,
    exdat = exdat,
    ezdat = ezdat,
    bws = bw,
    B = ncol(counts),
    counts = counts,
    mode = "frozen"
  )

  expect_equal(exact.out$t0, frozen.out$t0, tolerance = 1e-12)
  expect_equal(exact.out$t, frozen.out$t, tolerance = 1e-12)
})

test_that("npscoef fixed lp helper matches duplicate-row oracle", {
  boot.fun <- getFromNamespace(".np_inid_boot_from_scoef", "np")

  set.seed(20260318)
  n <- 28L
  xdat <- data.frame(x = runif(n, -1, 1))
  zdat <- data.frame(z = rnorm(n))
  y <- with(xdat, x^2 + zdat$z + rnorm(n, sd = 0.1))

  ex.seq <- seq(
    np:::trim.quantiles(xdat[, 1L], 0.05)[1L],
    np:::trim.quantiles(xdat[, 1L], 0.05)[2L],
    length.out = 4L
  )
  ez.seq <- seq(
    np:::trim.quantiles(zdat[, 1L], 0.05)[1L],
    np:::trim.quantiles(zdat[, 1L], 0.05)[2L],
    length.out = 4L
  )
  eval.grid <- expand.grid(x = ex.seq, z = ez.seq)

  bw <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    bws = c(0.7),
    bwtype = "fixed",
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 1,
    bernstein.basis = FALSE
  )

  counts <- cbind(
    c(rep(2L, 6L), rep(1L, n - 12L), rep(0L, 6L)),
    c(rep(0L, 4L), rep(1L, n - 8L), rep(2L, 4L)),
    c(3L, rep(1L, n - 3L), 2L, 0L)
  )

  helper.out <- boot.fun(
    txdat = xdat,
    ydat = y,
    tzdat = zdat,
    exdat = eval.grid["x"],
    ezdat = eval.grid["z"],
    bws = bw,
    B = ncol(counts),
    counts = counts,
    mode = "exact"
  )

  oracle <- vapply(seq_len(ncol(counts)), function(j) {
    idx <- np:::.np_counts_to_indices(counts[, j])
    as.vector(npscoef(
      bws = bw,
      txdat = xdat[idx, , drop = FALSE],
      tzdat = zdat[idx, , drop = FALSE],
      tydat = y[idx],
      exdat = eval.grid["x"],
      ezdat = eval.grid["z"],
      iterate = FALSE,
      errors = FALSE
    )$mean)
  }, numeric(nrow(eval.grid)))
  oracle <- t(oracle)

  expect_equal(helper.out$t, oracle, tolerance = 1e-10)
})

test_that("npscoef frozen nonfixed local-polynomial helpers stay finite", {
  boot.fun <- getFromNamespace(".np_inid_boot_from_scoef", "np")

  set.seed(20260318)
  n <- 40L
  xdat <- data.frame(x = runif(n, -1, 1))
  zdat <- data.frame(z = rnorm(n))
  y <- with(xdat, x^2 + zdat$z + rnorm(n, sd = 0.1))
  exdat <- data.frame(x = seq(
    np:::trim.quantiles(xdat[, 1L], 0.05)[1L],
    np:::trim.quantiles(xdat[, 1L], 0.05)[2L],
    length.out = 5L
  ))
  ezdat <- data.frame(z = seq(
    np:::trim.quantiles(zdat[, 1L], 0.05)[1L],
    np:::trim.quantiles(zdat[, 1L], 0.05)[2L],
    length.out = 5L
  ))
  counts <- rmultinom(n = 5L, size = n, prob = rep.int(1 / n, n))

  configs <- list(
    list(
      label = "adaptive_lp_bernstein",
      bw = npscoefbw(
        xdat = xdat,
        zdat = zdat,
        ydat = y,
        nmulti = 1L,
        regtype = "lp",
        degree = 1L,
        bernstein.basis = TRUE,
        bwtype = "adaptive_nn"
      )
    ),
    list(
      label = "generalized_ll",
      bw = npscoefbw(
        xdat = xdat,
        zdat = zdat,
        ydat = y,
        nmulti = 1L,
        regtype = "ll",
        bwtype = "generalized_nn"
      )
    )
  )

  for (cfg in configs) {
    exact.out <- boot.fun(
      txdat = xdat,
      ydat = y,
      tzdat = zdat,
      exdat = exdat,
      ezdat = ezdat,
      bws = cfg$bw,
      B = ncol(counts),
      counts = counts,
      mode = "exact"
    )
    frozen.out <- boot.fun(
      txdat = xdat,
      ydat = y,
      tzdat = zdat,
      exdat = exdat,
      ezdat = ezdat,
      bws = cfg$bw,
      B = ncol(counts),
      counts = counts,
      mode = "frozen"
    )

    expect_equal(frozen.out$t0, exact.out$t0, tolerance = 1e-10, info = cfg$label)
    expect_true(all(is.finite(frozen.out$t)), info = cfg$label)
    expect_true(max(abs(frozen.out$t)) < 1e4, info = cfg$label)
  }
})
