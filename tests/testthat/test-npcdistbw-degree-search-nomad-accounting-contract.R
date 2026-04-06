test_that("npcdist NOMAD accounting is owner-level on fast-path fits", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 300L
  x <- runif(n)
  y <- runif(n)

  fit <- npcdist(
    y ~ x,
    nomad = TRUE,
    search.engine = "nomad",
    cykerbound = "range"
  )

  expect_gt(as.numeric(fit$bws$num.feval.fast[1L]), 0)
  expect_lte(
    as.numeric(fit$bws$num.feval.fast[1L]),
    as.numeric(fit$bws$num.feval[1L])
  )

  ev <- npRmpi:::.npcdistbw_eval_only(
    data.frame(x = x),
    data.frame(y = y),
    NULL,
    fit$bws,
    invalid.penalty = "baseline",
    penalty.multiplier = 10
  )

  expect_equal(as.numeric(ev$num.feval), 1)
  expect_equal(as.numeric(ev$num.feval.fast), 1)
})
