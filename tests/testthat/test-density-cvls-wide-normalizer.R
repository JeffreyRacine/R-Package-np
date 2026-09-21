test_that("density CVLS uses a wide sample-count normalizer", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE, np.largelambda = TRUE)
  on.exit(options(old), add = TRUE)
  # At lambda=1, every normalized two-category kernel is exactly 1/2.
  # The squared-density integral and LOO average are both 1/2, at any n.
  # 66000 crosses both signed-int and unsigned-int product wrap boundaries.
  for (n in c(20L, 66000L)) {
    x <- data.frame(x = factor(rep(1:2, length.out = n)))
    bw <- npudensbw(dat = x, bws = 1, ukertype = "liracine",
                    bwmethod = "cv.ls", bandwidth.compute = FALSE)
    for (compress in c(FALSE, TRUE)) {
      options(np.categorical.compress = compress)
      evaluated <- npudensbw(dat = x, bws = bw, eval.only = TRUE)
      expect_equal(-as.double(evaluated$fval), -0.5, tolerance = 1e-13,
                   info = paste("n =", n, "compression =", compress))
      expect_equal(as.double(evaluated$bw), 1, tolerance = 0)
    }
  }
})
