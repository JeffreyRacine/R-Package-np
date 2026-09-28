test_that("adaptive external-grid CDF CV gathers every donor weight", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  dat <- data.frame(x = c(0.03, 0.12, 0.24, 0.43, 0.55, 0.69, 0.81, 0.96))
  grid <- data.frame(x = c(0.01, 0.08, 0.31, 0.49, 0.72, 0.91, 0.99))
  n <- nrow(dat)
  # Each fold uses its literal deleted sample. A donor's ANN radius excludes
  # both the held observation and the donor. Keep Gaussian CDF precision.
  leave_one_out <- t(vapply(seq_len(n), function(held) {
    donors <- setdiff(seq_len(n), held)
    radius <- vapply(donors, function(j) {
      other <- setdiff(donors, j)
      sort(abs(dat$x[other] - dat$x[j]))[[3L]]
    }, numeric(1L))
    stopifnot(all(is.finite(radius)), all(radius > 0))
    vapply(grid$x, function(z)
      mean(stats::pnorm((z - dat$x[donors]) / radius *
                        (sqrt(2) * 0.7071067810))), numeric(1L))
  }, numeric(nrow(grid))))
  expected <- mean((outer(dat$x, grid$x, "<=") - leave_one_out)^2)
  bw <- npudistbw(dat = dat, bws = 3, bwtype = "adaptive_nn",
                 bwscaling = FALSE, bandwidth.compute = FALSE)
  first <- npRmpi:::npudistbw.dbandwidth(dat = dat, gdat = grid, bws = bw,
                                        eval.only = TRUE, invalid.penalty = "dbmax")
  second <- npRmpi:::npudistbw.dbandwidth(dat = dat, gdat = grid, bws = bw,
                                         eval.only = TRUE, invalid.penalty = "dbmax")
  expect_equal(as.double(first$fval), expected, tolerance = 2e-12)
  expect_identical(first$fval, second$fval)
})
