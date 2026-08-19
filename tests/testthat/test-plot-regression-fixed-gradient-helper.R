test_that("fixed ll gradient helper returns finite gradients without exact fallback", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  old.calls <- getOption("npRmpi.regression.exact.calls")
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(options(npRmpi.regression.exact.calls = old.calls), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

  set.seed(20260310)
  n <- 10
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + 0.4 * x2 + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(3), , drop = FALSE]
  B <- 2L
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

  expect_identical(dim(fast.out$t), c(B, nrow(ex)))
  expect_identical(length(fast.out$t0), nrow(ex))
  expect_true(all(is.finite(fast.out$t)))
  expect_true(all(is.finite(fast.out$t0)))
  expect_identical(fast.calls, 0L)
})

test_that("fixed lp gradient helper preserves derivative order and counts drawer semantics", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

  set.seed(20260311)
  n <- 10
  x1 <- runif(n)
  x2 <- runif(n)
  y <- cos(2 * pi * x1) + x2^2 + rnorm(n, sd = 0.06)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(3), , drop = FALSE]
  B <- 2L
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

  expect_identical(dim(fast.counts$t), c(B, nrow(ex)))
  expect_identical(length(fast.counts$t0), nrow(ex))
  expect_true(all(is.finite(fast.counts$t)))
  expect_true(all(is.finite(fast.counts$t0)))
  expect_equal(fast.drawer$t, fast.counts$t, tolerance = 1e-10)
  expect_equal(as.vector(fast.drawer$t0), as.vector(fast.counts$t0), tolerance = 0)
})

test_that("fixed lp plot output records requested gradient order", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

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
