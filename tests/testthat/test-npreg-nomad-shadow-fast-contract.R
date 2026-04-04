test_that("npreg shadow NOMAD evaluation reports fast hits for large-h solutions", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 1000
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  xdat <- data.frame(x = x)

  bw <- npregbw(y ~ x, nomad = TRUE)

  shadow <- npRmpi:::.npregbw_nomad_shadow_begin(
    xdat = xdat,
    ydat = y,
    bws = bw,
    start.bw = bw$bw
  )

  out <- npRmpi:::.npregbw_nomad_shadow_eval(
    shadow = shadow,
    bw = bw$bw,
    degree = bw$degree
  )
  npRmpi:::.npregbw_nomad_shadow_end(shadow)

  expect_equal(out[1L], bw$fval, tolerance = 1e-10)
  expect_gt(out[2L], 0)
})
