test_that("npcdist NOMAD accounting is owner-level on fast-path fits", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(41)
  n <- 220L
  x <- runif(n)
  y <- pnorm(2 * sin(2 * pi * x) + rnorm(n))

  fit <- npcdist(
    y ~ x,
    nomad = TRUE,
    search.engine = "nomad",
    degree.max = 2L
  )

  expect_true(any(as.integer(fit$bws$degree.engine) > 0L))
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
  expect_equal(
    as.numeric(fit$bws$fval[1L]),
    as.numeric(ev$objective[1L]),
    tolerance = 1e-12
  )
})
