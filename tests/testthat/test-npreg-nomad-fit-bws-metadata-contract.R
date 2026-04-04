test_that("npreg nomad fit preserves bandwidth telemetry in returned bws", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 1000
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)

  set.seed(42)
  bw <- npregbw(y ~ x, nomad = TRUE)

  set.seed(42)
  fit <- npreg(y ~ x, nomad = TRUE)

  expect_identical(fit$bws$method, bw$method)
  expect_identical(fit$bws$pmethod, bw$pmethod)
  expect_equal(fit$bws$fval, bw$fval, tolerance = 1e-12)
  expect_equal(fit$bws$ifval, bw$ifval, tolerance = 1e-12)
  expect_equal(fit$bws$num.feval, bw$num.feval)
  expect_equal(fit$bws$num.feval.fast, bw$num.feval.fast)
})
