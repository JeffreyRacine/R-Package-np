test_that("Ichimura GNN covariance preserves all free predictor rows", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(626L)
  x <- as.data.frame(matrix(runif(120L, -1, 1), 30L, 4L))
  beta <- c(1, .6, -.25, .2)
  for (p in 2:4) {
    xp <- x[seq_len(p)]
    y <- sin(as.vector(as.matrix(xp) %*% beta[seq_len(p)])) + .1*cos(1:30)
    bw <- npindexbw(xdat = xp, ydat = y, method = "ichimura",
      bwtype = "generalized_nn", bws = c(beta[seq_len(p)], 20),
      bandwidth.compute = FALSE)
    fit <- npindex(bws = bw, txdat = xp, tydat = y, gradients = TRUE)
    expect_identical(dim(vcov(fit)), c(p, p))
    expect_true(all(is.finite(vcov(fit))))
    expect_identical(unname(vcov(fit)[1L, ]), rep(0, p))
    expect_identical(unname(vcov(fit)[, 1L]), rep(0, p))
  }
})
