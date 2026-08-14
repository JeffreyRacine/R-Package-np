test_that("MPI npscoef tree outer acceleration preserves canonical moments", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
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
    xdat = xdat, ydat = ydat, zdat = zdat,
    bws = c(0.24, 0.24), regtype = "lp", degree = c(2L, 2L),
    ckertype = "epanechnikov"
  )
  ctx <- npRmpi:::.npscoefbw_nomad_context_prepare(xdat, ydat, zdat)
  state <- npRmpi:::.npscoefbw_nomad_moment_state(ctx, bw, leave.one.out = TRUE)
  args <- list(
    txdat = state$z.train, tydat = state$ytensor,
    weights = state$ytensor, bws = state$kernel.bws,
    leave.one.out = TRUE, bandwidth.divide = state$bandwidth.divide
  )

  incumbent <- do.call(npksum, args)
  accelerated <- do.call(npRmpi:::.npscoef_npksum, args)
  expect_equal(accelerated$ksum, incumbent$ksum, tolerance = 1e-11)
  expect_identical(do.call(npksum, args)$ksum, incumbent$ksum)
})

test_that("MPI npscoef full-support tree workspace is deterministic at acceptance size", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
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
  ctx <- npRmpi:::.npscoefbw_nomad_context_prepare(xdat, ydat, zdat)
  state <- npRmpi:::.npscoefbw_nomad_moment_state(ctx, bw, leave.one.out = TRUE)
  args <- list(
    txdat = state$z.train, tydat = state$ytensor,
    weights = state$ytensor, bws = state$kernel.bws,
    leave.one.out = TRUE, bandwidth.divide = state$bandwidth.divide
  )

  ordinary.before <- do.call(npksum, args)
  first <- do.call(npRmpi:::.npscoef_npksum, args)
  second <- do.call(npRmpi:::.npscoef_npksum, args)
  ordinary.after <- do.call(npksum, args)

  expect_identical(second$ksum, first$ksum)
  expect_equal(first$ksum, ordinary.before$ksum, tolerance = 1e-11)
  expect_identical(ordinary.after$ksum, ordinary.before$ksum)
})
