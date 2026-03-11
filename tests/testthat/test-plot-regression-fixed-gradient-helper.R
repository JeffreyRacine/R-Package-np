test_that("fixed ll gradient helper matches explicit refits without exact fallback", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  old.calls <- getOption("npRmpi.regression.exact.calls")
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(options(npRmpi.regression.exact.calls = old.calls), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

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

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "npRmpi")
  np.ns <- asNamespace("npRmpi")
  options(npRmpi.regression.exact.calls = 0L)
  trace(
    what = ".np_inid_boot_from_regression_exact",
    where = np.ns,
    tracer = quote(options(
      npRmpi.regression.exact.calls = getOption("npRmpi.regression.exact.calls", 0L) + 1L
    )),
    print = FALSE
  )
  on.exit(untrace(".np_inid_boot_from_regression_exact", where = np.ns), add = TRUE)

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
  fast.calls <- getOption("npRmpi.regression.exact.calls")

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
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

  set.seed(20260311)
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

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "npRmpi")
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
  expect_equal(fast.drawer$t, fast.counts$t, tolerance = 0)
  expect_equal(as.vector(fast.drawer$t0), as.vector(fast.counts$t0), tolerance = 0)
})
