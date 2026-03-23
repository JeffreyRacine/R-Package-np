test_that("predict returns log-likelihood metadata for npudens and npcdens", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260323)

  x <- data.frame(x = c(-1.2, -0.5, 0.0, 0.4, 0.9, 1.4))
  y <- c(-0.9, -0.2, 0.1, 0.5, 1.0, 1.6)

  bw.u <- npudensbw(dat = x, bws = 0.6, bandwidth.compute = FALSE)
  fit.u <- npudens(bws = bw.u, tdat = x)
  out.u <- predict(fit.u, se.fit = TRUE)

  expect_false(is.null(out.u$log.likelihood))
  expect_equal(out.u$log.likelihood, fit.u$log_likelihood)

  bw.c <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.6, 0.7),
    bandwidth.compute = FALSE
  )
  fit.c <- npcdens(bws = bw.c, txdat = x, tydat = y)
  out.c <- predict(fit.c, se.fit = TRUE)

  expect_false(is.null(out.c$log.likelihood))
  expect_equal(out.c$log.likelihood, fit.c$log_likelihood)
})
