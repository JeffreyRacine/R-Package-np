test_that("nonfixed regression NN bandwidth constructor enforces k >= 2 across lc ll and lp", {
  set.seed(42)
  n <- 80
  x <- sort(rnorm(n))
  y <- x^2 + rnorm(n, sd = 0.15)
  xdati <- npRmpi:::untangle(data.frame(x = x))
  ydati <- npRmpi:::untangle(data.frame(y = y))

  expect_error(
    npRmpi:::rbandwidth(bw = 1, bandwidth = 1, regtype = "lc", bwtype = "generalized_nn",
                        bandwidth.compute = FALSE, nobs = n,
                        xdati = xdati, ydati = ydati, xnames = "x", ynames = "y"),
    "nearest-neighbor bandwidth must be at least 2"
  )

  expect_error(
    npRmpi:::rbandwidth(bw = 1, bandwidth = 1, regtype = "ll", bwtype = "generalized_nn",
                        bandwidth.compute = FALSE, nobs = n,
                        xdati = xdati, ydati = ydati, xnames = "x", ynames = "y"),
    "nearest-neighbor bandwidth must be at least 2"
  )

  expect_error(
    npRmpi:::rbandwidth(bw = 1, bandwidth = 1, regtype = "lp", degree = 2L, bwtype = "generalized_nn",
                        bandwidth.compute = FALSE, nobs = n,
                        xdati = xdati, ydati = ydati, xnames = "x", ynames = "y"),
    "nearest-neighbor bandwidth must be at least 2"
  )

  bw.ll <- npRmpi:::rbandwidth(bw = 2, bandwidth = 2, regtype = "ll", bwtype = "generalized_nn",
                               bandwidth.compute = FALSE, nobs = n,
                               xdati = xdati, ydati = ydati, xnames = "x", ynames = "y")
  expect_equal(as.integer(bw.ll$bw[1]), 2L)
})

test_that("semiparametric NN helper floors are aligned at k >= 2", {
  expect_error(
    npRmpi:::.npindex_finalize_bandwidth(1, "generalized_nn", 80L, where = "npindexbw"),
    "nearest-neighbor bandwidth candidate must map to an integer in \\[2,"
  )
  expect_equal(
    npRmpi:::.npindex_finalize_bandwidth(2, "generalized_nn", 80L, where = "npindexbw"),
    2
  )

  expect_equal(
    npRmpi:::.npscoef_finalize_bandwidth(c(1), "adaptive_nn", 80L, where = "npscoefbw"),
    2
  )
  expect_equal(
    npRmpi:::.npscoef_finalize_bandwidth(c(2), "adaptive_nn", 80L, where = "npscoefbw"),
    2
  )
})
