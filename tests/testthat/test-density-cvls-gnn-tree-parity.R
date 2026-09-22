test_that("unbounded generalized-NN density CVLS is tree-invariant", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(1901)
  dat <- data.frame(x = rnorm(8L))

  for (kernel in c("gaussian", "epanechnikov", "uniform")) {
    values <- vapply(list(FALSE, TRUE, "auto"), function(tree) {
      options(np.tree = tree)
      bws <- npudensbw(dat = dat, bws = 2, bandwidth.compute = FALSE,
                       bwmethod = "cv.ls", bwtype = "generalized_nn",
                       ckertype = kernel)
      out <- npudensbw(dat = dat, bws = bws, bandwidth.compute = TRUE,
                       nmulti = 1L, powell.remin = FALSE,
                       bwsolver = "powell", eval.only = TRUE)
      expect_true(is.finite(out$fval))
      out$fval
    }, numeric(1))
    expect_equal(values, rep(values[[1L]], 3L), tolerance = 2e-12,
                 info = kernel)
  }
})
