library(np)

quiet_eval <- function(expr) {
  value <- NULL
  capture.output(value <- force(expr))
  value
}

test_that("npplreg exact helper matches duplicate-row oracle and adaptive frozen stays finite", {
  boot.fun <- getFromNamespace(".np_inid_boot_from_plreg", "np")

  txdat <- data.frame(x = c(0.05, 0.05, 0.25, 0.25, 0.60, 0.60, 0.90, 0.90))
  tzdat <- data.frame(z = c(0.10, 0.10, 0.35, 0.35, 0.65, 0.65, 0.95, 0.95))
  y <- with(txdat, 0.5 + 0.8 * x + sin(2 * pi * tzdat$z))
  exdat <- data.frame(x = c(0.10, 0.30, 0.70))
  ezdat <- data.frame(z = c(0.12, 0.40, 0.85))
  counts <- cbind(
    c(2, 0, 2, 0, 1, 1, 1, 1),
    c(0, 2, 0, 2, 1, 1, 1, 1),
    c(1, 1, 3, 0, 0, 2, 1, 0)
  )

  bw <- npplregbw(
    xdat = txdat,
    zdat = tzdat,
    ydat = y,
    bws = matrix(c(3L, 3L), nrow = 2L, ncol = 1L),
    bwtype = "adaptive_nn",
    bandwidth.compute = FALSE,
    regtype = "ll"
  )

  exact.out <- boot.fun(
    txdat = txdat,
    ydat = y,
    tzdat = tzdat,
    exdat = exdat,
    ezdat = ezdat,
    bws = bw,
    B = ncol(counts),
    counts = counts,
    mode = "exact"
  )
  frozen.out <- boot.fun(
    txdat = txdat,
    ydat = y,
    tzdat = tzdat,
    exdat = exdat,
    ezdat = ezdat,
    bws = bw,
    B = ncol(counts),
    counts = counts,
    mode = "frozen"
  )

  exact.oracle <- vapply(seq_len(ncol(counts)), function(j) {
    idx <- np:::.np_counts_to_indices(counts[, j])
    as.vector(npplreg(
      bws = bw,
      txdat = txdat[idx, , drop = FALSE],
      tydat = y[idx],
      tzdat = tzdat[idx, , drop = FALSE],
      exdat = exdat,
      ezdat = ezdat
    )$mean)
  }, numeric(nrow(exdat)))
  exact.oracle <- t(exact.oracle)

  t0.oracle <- as.vector(npplreg(
    bws = bw,
    txdat = txdat,
    tydat = y,
    tzdat = tzdat,
    exdat = exdat,
    ezdat = ezdat
  )$mean)

  expect_equal(exact.out$t, exact.oracle, tolerance = 1e-10)
  expect_equal(as.vector(exact.out$t0), t0.oracle, tolerance = 1e-10)
  expect_equal(as.vector(frozen.out$t0), t0.oracle, tolerance = 1e-10)
  expect_true(all(is.finite(frozen.out$t)))
  expect_lt(max(abs(frozen.out$t)), 10 * max(abs(exact.out$t)) + 1e-8)
  expect_false(isTRUE(all.equal(
    exact.out$t,
    frozen.out$t,
    tolerance = 1e-10,
    check.attributes = FALSE
  )))
})

test_that("npplreg bootstrap plot payload stays finite across frozen nonfixed bw types", {
  set.seed(42)
  n <- 100
  x <- runif(n, -1, 1)
  z <- rnorm(n)
  y <- x + rnorm(n, sd = 0.25 * sd(x))

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    fit <- npplreg(y ~ x | z, nmulti = 1, bwtype = bwtype)

    for (mode in c("exact", "frozen")) {
      out <- quiet_eval(plot(
        fit,
        plot.errors.method = "bootstrap",
        view = "fixed",
        plot.errors.boot.method = "geom",
        neval = 12,
        plot.errors.boot.num = 61,
        plot.errors.type = "pointwise",
        plot.errors.boot.nonfixed = mode,
        plot.behavior = "data"
      ))

      expect_true(all(is.finite(out$r1$merr)), info = paste(bwtype, mode))
      expect_true(max(abs(out$r1$merr)) < 1e4, info = paste(bwtype, mode))
    }
  }
})

test_that("npplreg frozen bootstrap mode is forwarded into the helper", {
  txdat <- data.frame(x = c(0.05, 0.05, 0.25, 0.25, 0.60, 0.60, 0.90, 0.90))
  tzdat <- data.frame(z = c(0.10, 0.10, 0.35, 0.35, 0.65, 0.65, 0.95, 0.95))
  y <- with(txdat, 0.5 + 0.8 * x + sin(2 * pi * tzdat$z))

  bw <- npplregbw(
    xdat = txdat,
    zdat = tzdat,
    ydat = y,
    bws = matrix(c(3L, 3L), nrow = 2L, ncol = 1L),
    bwtype = "adaptive_nn",
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  ns <- asNamespace("np")
  orig <- getFromNamespace(".np_inid_boot_from_plreg", "np")
  modes <- character()
  assignInNamespace(".np_inid_boot_from_plreg", function(..., mode = c("exact", "frozen")) {
    mode <- match.arg(mode)
    modes <<- c(modes, mode)
    orig(..., mode = mode)
  }, ns = "np")
  on.exit(assignInNamespace(".np_inid_boot_from_plreg", orig, ns = "np"), add = TRUE)

  compute.boot <- getFromNamespace("compute.bootstrap.errors.plbandwidth", "np")

  expect_no_error(compute.boot(
    xdat = txdat,
    ydat = y,
    zdat = tzdat,
    exdat = data.frame(x = seq(0.10, 0.80, length.out = 6L)),
    ezdat = data.frame(z = seq(0.12, 0.92, length.out = 6L)),
    gradients = FALSE,
    slice.index = 0L,
    progress.target = NULL,
    plot.errors.boot.method = "inid",
    plot.errors.boot.nonfixed = "frozen",
    plot.errors.boot.num = 41L,
    plot.errors.center = "estimate",
    plot.errors.type = "pointwise",
    plot.errors.alpha = 0.05,
    bws = bw
  ))

  expect_gte(length(modes), 1L)
  expect_true(all(modes == "frozen"))
})
