test_that("fixed ll gradient helper matches explicit refits without direct regression calls", {
  skip_if_not_installed("np")

  set.seed(20260310)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + 0.4 * x2 + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(15), , drop = FALSE]
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  bw <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "ll",
    bws = c(0.3, 0.3),
    bandwidth.compute = FALSE
  )

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "np")
  np.ns <- asNamespace("np")
  old.calls <- getOption("np.regression.direct.calls")
  options(np.regression.direct.calls = 0L)
  on.exit(options(np.regression.direct.calls = old.calls), add = TRUE)
  trace(
    what = ".np_regression_direct",
    where = np.ns,
    tracer = quote(options(np.regression.direct.calls = getOption("np.regression.direct.calls", 0L) + 1L)),
    print = FALSE
  )
  on.exit(untrace(".np_regression_direct", where = np.ns), add = TRUE)

  options(np.regression.direct.calls = 0L)
  fast.out <- fast.fun(
    xdat = tx,
    exdat = ex,
    bws = bw,
    ydat = y,
    B = B,
    counts = counts,
    gradients = TRUE,
    gradient.order = 1L,
    slice.index = 1L
  )
  fast.calls <- getOption("np.regression.direct.calls")

  explicit.t <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    explicit.t[b, ] <- npreg(
      txdat = tx[idx, , drop = FALSE],
      tydat = y[idx],
      exdat = ex,
      bws = bw,
      gradients = TRUE,
      gradient.order = 1L,
      warn.glp.gradient = FALSE
    )$grad[, 1L]
  }

  fit0 <- npreg(
    txdat = tx,
    tydat = y,
    exdat = ex,
    bws = bw,
    gradients = TRUE,
    gradient.order = 1L,
    warn.glp.gradient = FALSE
  )$grad[, 1L]

  expect_equal(fast.out$t, explicit.t, tolerance = 1e-10)
  expect_equal(as.vector(fast.out$t0), as.vector(fit0), tolerance = 1e-10)
  expect_identical(fast.calls, 0L)
})

test_that("fixed lp gradient helper preserves derivative order and counts drawer semantics", {
  skip_if_not_installed("np")

  set.seed(20260310 + 1L)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  y <- cos(2 * pi * x1) + x2^2 + rnorm(n, sd = 0.06)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(12), , drop = FALSE]
  B <- 7L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  bw <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    degree = c(2L, 2L),
    basis = "glp",
    bernstein.basis = FALSE,
    bws = c(0.3, 0.3),
    bandwidth.compute = FALSE
  )

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "np")
  fast.counts <- fast.fun(
    xdat = tx,
    exdat = ex,
    bws = bw,
    ydat = y,
    B = B,
    counts = counts,
    gradients = TRUE,
    gradient.order = c(2L, 1L),
    slice.index = 1L
  )

  drawer <- function(start, stopi) counts[, start:stopi, drop = FALSE]
  fast.drawer <- fast.fun(
    xdat = tx,
    exdat = ex,
    bws = bw,
    ydat = y,
    B = B,
    counts.drawer = drawer,
    gradients = TRUE,
    gradient.order = c(2L, 1L),
    slice.index = 1L
  )

  explicit.t <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    explicit.t[b, ] <- npreg(
      txdat = tx[idx, , drop = FALSE],
      tydat = y[idx],
      exdat = ex,
      bws = bw,
      gradients = TRUE,
      gradient.order = c(2L, 1L),
      warn.glp.gradient = FALSE
    )$grad[, 1L]
  }

  fit0 <- npreg(
    txdat = tx,
    tydat = y,
    exdat = ex,
    bws = bw,
    gradients = TRUE,
    gradient.order = c(2L, 1L),
    warn.glp.gradient = FALSE
  )$grad[, 1L]

  expect_equal(fast.counts$t, explicit.t, tolerance = 1e-8)
  expect_equal(as.vector(fast.counts$t0), as.vector(fit0), tolerance = 1e-8)
  expect_equal(fast.drawer$t, fast.counts$t, tolerance = 1e-10)
  expect_equal(as.vector(fast.drawer$t0), as.vector(fast.counts$t0), tolerance = 0)
})

test_that("fixed lp plot output records requested gradient order", {
  skip_if_not_installed("np")

  n <- 45
  x <- seq(-1, 1, length.out = n)
  y <- 1 + 2 * x + 3 * x^2
  tx <- data.frame(x = x)

  bw <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    degree = 2L,
    basis = "glp",
    bws = 100,
    bwscaling = FALSE,
    bandwidth.compute = FALSE
  )

  fit <- npreg(
    bws = bw,
    gradients = TRUE,
    gradient.order = 1L,
    warn.glp.gradient = FALSE
  )

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  out <- plot(
    fit,
    gradients = TRUE,
    gradient_order = 2L,
    errors = "none",
    output = "data",
    neval = 20
  )

  expect_equal(as.vector(out[[1L]]$grad), rep(6, 20), tolerance = 1e-8)
  expect_identical(out[[1L]]$gradient.order, 2L)
})
