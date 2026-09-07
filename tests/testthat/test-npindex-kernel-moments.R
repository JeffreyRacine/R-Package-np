test_that("single-index raw moments share one kernel traversal", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(621)
  x <- data.frame(index = runif(32, -1, 1))
  y <- cbind(sin(x$index), cos(x$index))
  real.sum <- .np_index_kernel_sum
  calls <- 0L
  counted <- function(...) {
    args <- list(...)
    expect_false(isTRUE(args$return.kernel.weights))
    calls <<- calls + 1L
    real.sum(...)
  }
  testthat::local_mocked_bindings(
    .np_index_kernel_sum = counted, .package = "npRmpi"
  )
  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    for (m in c(1L, 4L)) {
      for (q in c(1L, 2L)) {
        args <- list(txdat = x, exdat = x[seq_len(m), , drop = FALSE],
                     bws = if (type == "fixed") .8 else 24, bwtype = type,
                     ckertype = "gaussian", ckerorder = 2)
        before <- calls
        got <- do.call(.np_index_kernel_moments,
                       c(list(y = y[, seq_len(q), drop = FALSE]), args))
        expect_identical(calls - before, 1L)
        numerator <- do.call(real.sum,
                            c(args, list(tydat = rep(1, nrow(x)),
                                         weights = y[, seq_len(q), drop = FALSE])))$ksum
        denominator <- do.call(real.sum, args)$ksum
        expect_equal(unname(got$numerator),
                     matrix(as.numeric(numerator), nrow = q), tolerance = 1e-14)
        expect_equal(got$denominator, as.numeric(denominator), tolerance = 1e-14)
      }
    }
  }
})
