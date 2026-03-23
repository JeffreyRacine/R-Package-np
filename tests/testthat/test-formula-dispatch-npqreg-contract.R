test_that("named bws formula dispatch matches explicit npqreg bandwidth route", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260323)
  d <- data.frame(
    x = seq(0.1, 0.9, length.out = 9),
    y = seq(0.2, 1.0, length.out = 9)
  )
  nd <- data.frame(x = c(0.2, 0.5, 0.8))

  bw <- npcdistbw(y ~ x, data = d, newdata = nd)
  fit.pos <- npqreg(bws = bw, data = d, newdata = nd, tau = 0.4)
  fit.named <- npqreg(bws = y ~ x, data = d, newdata = nd, tau = 0.4)

  expect_s3_class(fit.named, "qregression")
  expect_equal(length(fit.named$quantile), nrow(nd), tolerance = 0)
  expect_lt(max(abs(as.numeric(fit.named$quantile) - as.numeric(fit.pos$quantile))), 1e-2)
  expect_lt(max(abs(as.numeric(fit.named$quanterr) - as.numeric(fit.pos$quanterr))), 1e-2)
})
