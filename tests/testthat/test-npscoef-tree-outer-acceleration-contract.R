test_that("npscoef tree outer acceleration preserves the canonical moments", {
  old <- options(np.messages = FALSE, np.tree = TRUE,
                 np.macMseries.accelerate = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260810L)
  n <- 160L
  zdat <- data.frame(z1 = runif(n, -1, 1), z2 = runif(n, -1, 1))
  xdat <- data.frame(x = rnorm(n))
  ydat <- xdat$x * (1 + sin(pi * zdat$z1) + 0.5 * zdat$z2) +
    rnorm(n, sd = 0.2)
  bw <- npscoefbw(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = c(0.24, 0.24),
    regtype = "lp",
    degree = c(2L, 2L),
    ckertype = "epanechnikov"
  )
  ctx <- np:::.npscoefbw_nomad_context_prepare(xdat, ydat, zdat)
  state <- np:::.npscoefbw_nomad_moment_state(
    ctx = ctx,
    bws = bw,
    leave.one.out = TRUE
  )
  args <- list(
    txdat = state$z.train,
    tydat = state$ytensor,
    weights = state$ytensor,
    bws = state$kernel.bws,
    leave.one.out = TRUE,
    bandwidth.divide = state$bandwidth.divide
  )

  incumbent <- do.call(npksum, args)
  accelerated <- do.call(np:::.npscoef_npksum, args)
  expect_equal(accelerated$ksum, incumbent$ksum, tolerance = 1e-11)

  ## The internal activation is one-call state and cannot affect a later
  ## ordinary npksum invocation, including after the accelerated call.
  ordinary.after <- do.call(npksum, args)
  expect_identical(ordinary.after$ksum, incumbent$ksum)
})

test_that("npscoef tree acceleration leaves dense and nearest-neighbor routes coherent", {
  old <- options(np.messages = FALSE, np.macMseries.accelerate = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260811L)
  n <- 90L
  zdat <- data.frame(z1 = runif(n), z2 = runif(n))
  xdat <- data.frame(x = rnorm(n))
  ydat <- xdat$x * (1 + zdat$z1) + rnorm(n, sd = 0.2)

  run_fixed <- function(tree) {
    options(np.tree = tree)
    npscoefbw(
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      bws = c(0.3, 0.3),
      regtype = "ll",
      ckertype = "epanechnikov"
    )
  }
  dense <- run_fixed(FALSE)
  tree <- run_fixed(TRUE)
  auto <- run_fixed("auto")
  expect_equal(tree$fval, dense$fval, tolerance = 1e-12)
  expect_equal(auto$fval, tree$fval, tolerance = 1e-12)

  options(np.tree = "auto")
  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    nn <- npscoefbw(
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      regtype = "ll",
      bwtype = bwtype,
      ckertype = "epanechnikov",
      nmulti = 1L,
      optim.maxit = 10L,
      optim.maxattempts = 1L
    )
    expect_true(all(is.finite(nn$bw)))
    expect_true(is.finite(nn$fval))
  }
})

test_that("npscoef full-support tree workspace is deterministic at acceptance size", {
  old <- options(np.messages = FALSE, np.tree = TRUE,
                 np.macMseries.accelerate = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260813L)
  n <- 1024L
  zdat <- data.frame(z1 = runif(n, -1, 1), z2 = runif(n, -1, 1))
  xdat <- data.frame(x = rnorm(n))
  ydat <- xdat$x * (1 + sin(pi * zdat$z1) + 0.5 * zdat$z2) +
    rnorm(n, sd = 0.2)
  bw <- npscoefbw(
    xdat = xdat, ydat = ydat, zdat = zdat,
    bws = c(1000, 1000), regtype = "lp", degree = c(3L, 3L),
    ckertype = "epanechnikov"
  )
  ctx <- np:::.npscoefbw_nomad_context_prepare(xdat, ydat, zdat)
  state <- np:::.npscoefbw_nomad_moment_state(ctx, bw, leave.one.out = TRUE)
  args <- list(
    txdat = state$z.train, tydat = state$ytensor,
    weights = state$ytensor, bws = state$kernel.bws,
    leave.one.out = TRUE, bandwidth.divide = state$bandwidth.divide
  )

  ordinary.before <- do.call(npksum, args)
  first <- do.call(np:::.npscoef_npksum, args)
  second <- do.call(np:::.npscoef_npksum, args)
  ordinary.after <- do.call(npksum, args)

  expect_identical(second$ksum, first$ksum)
  expect_equal(first$ksum, ordinary.before$ksum, tolerance = 1e-11)
  expect_identical(ordinary.after$ksum, ordinary.before$ksum)
})
