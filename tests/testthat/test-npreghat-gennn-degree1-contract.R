test_that("generalized-nn ll and canonical lp degree-1 share the same public/apply regression route", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260308)
  n <- 90
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + 0.25 * x + rnorm(n, sd = 0.04)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.05, 0.95, length.out = 25))

  bw.ll <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bws = 9,
    bandwidth.compute = FALSE
  )
  bw.lp <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bernstein.basis = FALSE,
    bwtype = "generalized_nn",
    bws = 9,
    bandwidth.compute = FALSE
  )

  fit.ll <- npreg(bws = bw.ll, txdat = tx, tydat = y, exdat = ex, warn.glp.gradient = FALSE)
  fit.lp <- npreg(bws = bw.lp, txdat = tx, tydat = y, exdat = ex, warn.glp.gradient = FALSE)

  hat.apply.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex, y = y, output = "apply")
  hat.apply.lp <- npreghat(bws = bw.lp, txdat = tx, exdat = ex, y = y, output = "apply")
  hat.matrix.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex)
  hat.matrix.lp <- npreghat(bws = bw.lp, txdat = tx, exdat = ex)

  expect_equal(as.vector(fit.ll$mean), as.vector(fit.lp$mean), tolerance = 1e-10)
  expect_equal(as.vector(hat.apply.ll), as.vector(hat.apply.lp), tolerance = 1e-10)
  expect_equal(as.vector(hat.apply.ll), as.vector(hat.matrix.ll %*% y), tolerance = 1e-10)
  expect_equal(as.vector(hat.apply.lp), as.vector(hat.matrix.lp %*% y), tolerance = 1e-10)
})
