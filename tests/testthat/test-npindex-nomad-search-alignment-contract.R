test_that("npindex NOMAD shortcut stays near the Powell basin on the canonical Ichimura repro", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 500
  dat <- data.frame(
    x1 = runif(n, min = -1, max = 1),
    x2 = runif(n, min = -1, max = 1)
  )
  dat$y <- dat$x1 - dat$x2 + rnorm(n)

  fit_powell <- npRmpi::npindex(
    txdat = dat[c("x1", "x2")],
    tydat = dat$y,
    method = "ichimura",
    regtype = "lp",
    bernstein = TRUE,
    degree = 1L
  )
  fit_nomad <- npRmpi::npindex(
    txdat = dat[c("x1", "x2")],
    tydat = dat$y,
    method = "ichimura",
    nomad = TRUE
  )

  expect_identical(as.integer(fit_nomad$bws$degree), 1L)
  expect_equal(as.numeric(fit_nomad$bws$fval), as.numeric(fit_powell$bws$fval), tolerance = 0.02)
  expect_equal(as.numeric(fit_nomad$bws$beta[2L]), as.numeric(fit_powell$bws$beta[2L]), tolerance = 0.1)
  expect_true(is.finite(fit_nomad$bws$bw))
  expect_lt(as.numeric(fit_nomad$bws$bw), 5)
  expect_lt(as.numeric(fit_nomad$bws$beta[2L]), 0)
})
