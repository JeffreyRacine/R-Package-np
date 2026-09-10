test_that("LC single-index apply preserves its matrix oracle without materializing H", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = getOption("np.tree"))
  on.exit(options(old), add = TRUE)
  set.seed(681)
  withCallingHandlers({
    x <- data.frame(x1 = runif(32, -.8, .8), x2 = runif(32, -.6, .6))
    z <- data.frame(index = x$x1 + .3*x$x2)
    y <- cbind(a = sin(z$index), b = cos(z$index))
    for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
      for (kernel in c("gaussian", "epanechnikov", "uniform")) {
        b <- npindexbw(xdat = x, ydat = y[, 1L], bws = c(1, .3, if (type == "fixed") .8 else 24),
                       bandwidth.compute = FALSE, bwtype = type, ckertype = kernel)
        for (m in c(1L, 4L)) {
          ze <- z[seq_len(m), , drop = FALSE]
          for (owner in list(.np_indexhat_lc_mean, .np_indexhat_core)) {
            H <- owner(b, z, ze, output = "matrix")
            actual <- owner(b, z, ze, y = y, output = "apply")
            expect_equal(actual, H %*% y, tolerance = 1e-10)
            expect_equal(owner(b, z, ze, y = y[, 1L], output = "apply"),
                         as.vector(H %*% y[, 1L]), tolerance = 1e-10)
          }
        }
        expect_equal(npindexhat(b, txdat = x, exdat = x[1:4, ], y = y, output = "apply"),
                     npindexhat(b, txdat = x, exdat = x[1:4, ], output = "matrix") %*% y,
                     tolerance = 1e-10)
      }
    }
    # Signed kernels and range bounds retain the supported-row matrix oracle.
    for (kernel in c("gaussian", "epanechnikov")) {
      b <- npindexbw(xdat = x, ydat = y[, 1L], bws = c(1, .3, .8),
                     bandwidth.compute = FALSE, ckertype = kernel,
                     ckerorder = 4, ckerbound = "range")
      for (owner in list(.np_indexhat_lc_mean, .np_indexhat_core)) {
        H <- owner(b, z, z[1:4, , drop = FALSE], output = "matrix")
        expect_equal(owner(b, z, z[1:4, , drop = FALSE], y, "apply"),
                     H %*% y, tolerance = 1e-10)
      }
    }
    b <- npindexbw(xdat = x, ydat = y[, 1L], bws = c(1, .3, .8),
                   bandwidth.compute = FALSE, ckertype = "epanechnikov")
    far <- data.frame(index = 10)
    for (owner in list(.np_indexhat_lc_mean, .np_indexhat_core)) {
      expect_error(owner(b, z, far, y, "apply"),
                   "single-index hat: zero normalizing weight sum", fixed = TRUE)
      expect_error(owner(b, z, far, output = "matrix"),
                   "single-index hat: zero normalizing weight sum", fixed = TRUE)
      for (value in c(NA_real_, Inf)) {
        yy <- y; yy[1, 1] <- value
        H <- owner(b, z, z[1:4, , drop = FALSE], output = "matrix")
        expect_identical(owner(b, z, z[1:4, , drop = FALSE], yy, "apply"), H %*% yy)
      }
    }
    # Actual computation must not ask either matrix owner for returned weights.
    testthat::local_mocked_bindings(
      .np_indexhat_lc_kernel_weights = function(...) stop("unexpected kw"),
      .np_kernel_weights_direct = function(...) stop("unexpected direct kw"),
      .package = "npRmpi"
    )
    for (tree in list(FALSE, TRUE, "auto")) {
      options(np.tree = tree)
      for (owner in list(.np_indexhat_lc_mean, .np_indexhat_core))
        expect_true(all(is.finite(owner(b, z, z[1:4, , drop = FALSE], y, "apply"))))
    }
  }, warning = function(w) {
    if (grepl("ignoring kernel order specified with uniform kernel type",
              conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
})
