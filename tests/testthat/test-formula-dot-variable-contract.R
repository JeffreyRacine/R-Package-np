test_that("formula variables may contain dots without triggering wildcard expansion", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260624L)
  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = factor(rbinom(80L, 1L, 0.5)))
  y.irr <- ordered(rbinom(80L, 1L, 0.5))

  bw.dens <- npcdensbw(y.irr ~ x, data = xdat, bws = c(0.2, 0.2),
                       bandwidth.compute = FALSE)
  bw.dist <- npcdistbw(y.irr ~ x, data = xdat, bws = c(0.2, 0.2),
                       bandwidth.compute = FALSE)

  expect_equal(bw.dens$variableNames$response, "y.irr")
  expect_equal(bw.dist$variableNames$response, "y.irr")
  expect_equal(bw.dens$variableNames$terms, "x")
  expect_equal(bw.dist$variableNames$terms, "x")
})
