test_that("kernel-sum subsets retain optional data and native-row identity", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  withr::local_options(np.messages = FALSE)
  x <- seq(-1, 1, length.out = 18)
  y <- sin(x)
  keep <- seq_along(x) %% 2L == 0L
  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    h <- if (type == "fixed") .4 else 4
    native <- npksum(txdat = data.frame(x = x[keep]), tydat = y[keep],
                     bws = h, bwtype = type)
    implicit <- npksum(y ~ x, subset = keep, bws = h, bwtype = type)
    explicit.null <- npksum(y ~ x, data = NULL, subset = keep, bws = h,
                            bwtype = type)
    expect_identical(implicit$ksum, native$ksum)
    expect_identical(explicit.null$ksum, native$ksum)
    by.expr <- npksum(y ~ x, subset = seq_along(x) %% 2L == 0L,
                      bws = h, bwtype = type)
    expect_identical(by.expr$ksum, native$ksum)
  }
  expect_error(npksum(y ~ x, subset = stop("optional subset sentinel"), bws = .4),
               "optional subset sentinel")
})
