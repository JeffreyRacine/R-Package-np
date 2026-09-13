npsig_streamed_information_contract <- function(package) {
  ns <- asNamespace(package)
  bwfun <- get("npregbw", ns)
  fitfun <- get("npreg", ns)
  tile <- get(".np_npsig_streamed_iid_tile", ns)
  statistic <- get(".np_npsig_statistic", ns)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = c(0, 3, 6, 9, 20, 20.1, 20.2, 20.3),
                  f = factor(rep(c("a", "b"), 4L)))
  y <- c(1, 2, 4, 3, 5, 7, 6, 9)
  responses <- cbind(y, 0, y + c(0, 0, 0, 0, .2, -.1, .1, -.2))
  bw <- bwfun(xdat = x, ydat = y, bws = c(.2, .25),
    bandwidth.compute = FALSE, regtype = "lc", ckertype = "epanechnikov")
  direct <- vapply(seq_len(ncol(responses)), function(column) {
    fit <- suppressWarnings(fitfun(bws = bw, txdat = x,
      tydat = responses[, column], gradients = TRUE, se = TRUE))
    certificate <- attr(fit, ".np.gradient.structural.zero", exact = TRUE)
    expect_true(all(certificate[1:4, 1L]))
    expect_equal(fit$grad[1:4, 1L], numeric(4L))
    expect_equal(fit$gerr[1:4, 1L], numeric(4L))
    statistic(fit, 1L, TRUE)
  }, numeric(1L))
  streamed <- tile(bw, x, 1L, response.matrix = responses,
    null.mean = y, residual.pool = y, pivotal = TRUE)
  expect_equal(streamed, direct, tolerance = 2e-12)
  expect_identical(streamed[[2L]], 0)

  # Certified rows contribute zero; they are not removed from the denominator.
  cluster <- 5:8
  cluster.bw <- bwfun(xdat = x[cluster, ], ydat = y[cluster],
    bws = c(.2, .25), bandwidth.compute = FALSE,
    regtype = "lc", ckertype = "epanechnikov")
  cluster.fit <- fitfun(bws = cluster.bw, txdat = x[cluster, ],
    tydat = y[cluster], gradients = TRUE, se = TRUE)
  expect_equal(streamed[[1L]], statistic(cluster.fit, 1L, TRUE) / 2,
               tolerance = 2e-12)

  # A rounded zero LP derivative is not the LC structural certificate.
  for (degree in 1:2) {
    lp.bw <- bwfun(xdat = x, ydat = y, bws = c(.2, .25),
      bandwidth.compute = FALSE, regtype = "lp", degree = degree,
      ckertype = "epanechnikov")
    fit <- suppressWarnings(fitfun(bws = lp.bw, txdat = x, tydat = y,
                                   gradients = TRUE, se = TRUE))
    expect_error(statistic(fit, 1L, TRUE), "standard error")
    expect_error(tile(lp.bw, x, 1L, response.matrix = matrix(y, ncol = 1L),
      null.mean = y, residual.pool = y, pivotal = TRUE),
      "npsigtest IID tile reduction failed")
  }
}

test_that("streamed pivotal inference shares public residual-information proofs", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  npsig_streamed_information_contract("npRmpi")
})
