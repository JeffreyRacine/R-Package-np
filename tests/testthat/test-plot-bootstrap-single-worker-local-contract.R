test_that("fixed regression wild bootstrap matches remote and single-worker local paths", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  n <- 40L
  B <- 19L
  x <- runif(n, -1, 1)
  y <- x + rnorm(n, sd = 0.25 * sd(x))
  xdat <- data.frame(x = x)
  exdat <- data.frame(x = seq(min(x), max(x), length.out = 17L))
  bw <- npregbw(y ~ x, regtype = "lp", degree = 1L, bwtype = "fixed", nmulti = 1)

  fit.mean <- as.vector(npreghat(
    bws = bw,
    txdat = xdat,
    exdat = xdat,
    y = y,
    output = "apply"
  ))
  H <- npreghat(
    bws = bw,
    txdat = xdat,
    exdat = exdat,
    output = "matrix"
  )

  boot.fun <- get(".np_plot_boot_from_hat_wild", asNamespace("npRmpi"))

  set.seed(20260319)
  remote <- boot.fun(
    H = H,
    ydat = y,
    fit.mean = fit.mean,
    B = B,
    wild = "rademacher",
    progress.label = NULL,
    prefer.local.single_worker = FALSE
  )

  set.seed(20260319)
  local <- boot.fun(
    H = H,
    ydat = y,
    fit.mean = fit.mean,
    B = B,
    wild = "rademacher",
    progress.label = NULL,
    prefer.local.single_worker = TRUE
  )

  expect_equal(local$t0, remote$t0)
  expect_equal(local$t, remote$t)
})

test_that("fixed regression inid bootstrap matches remote and single-worker local paths", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(84)
  n <- 40L
  B <- 19L
  x <- runif(n, -1, 1)
  y <- x + rnorm(n, sd = 0.25 * sd(x))
  xdat <- data.frame(x = x)
  exdat <- data.frame(x = seq(min(x), max(x), length.out = 17L))
  bw <- npregbw(y ~ x, regtype = "lp", degree = 1L, bwtype = "fixed", nmulti = 1)

  boot.fun <- get(".np_inid_boot_from_regression", asNamespace("npRmpi"))

  set.seed(20260319)
  remote <- boot.fun(
    xdat = xdat,
    exdat = exdat,
    bws = bw,
    ydat = y,
    B = B,
    gradients = FALSE,
    gradient.order = 1L,
    slice.index = 0L,
    prefer.local.single_worker = FALSE,
    prep.label = NULL,
    progress.label = NULL
  )

  set.seed(20260319)
  local <- boot.fun(
    xdat = xdat,
    exdat = exdat,
    bws = bw,
    ydat = y,
    B = B,
    gradients = FALSE,
    gradient.order = 1L,
    slice.index = 0L,
    prefer.local.single_worker = TRUE,
    prep.label = NULL,
    progress.label = NULL
  )

  expect_equal(local$t0, remote$t0)
  expect_equal(local$t, remote$t)
})
