test_that("npindexbw never replaces an invalid explicit starting bandwidth", {
  skip_if_not_installed("crs")
  expect_true(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  n <- 31L
  x <- data.frame(x1 = seq_len(n) / n,
                  x2 = cos(seq_len(n)) / 2 + 0.4)
  y <- (x$x1 + 0.5 * x$x2)^2 + sin(seq_len(n)) / 20
  original.default <- getFromNamespace(".npindex_default_start_bandwidth", "npRmpi")
  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    for (degree.search in c(FALSE, TRUE)) {
      default.calls <- 0L
      args <- list(xdat = x, ydat = y,
                    bws = c(1, 0.5, if (bwtype == "fixed") 5e-13 else 1),
                    method = "ichimura", bwtype = bwtype,
                    bandwidth.compute = TRUE, nmulti = 1L,
                    nomad = degree.search, random.seed = 20260902L)
      if (degree.search) {
        args$search.engine <- "nomad"
        args$degree.max <- 1L
        args$nomad.opts <- list(MAX_BB_EVAL = 12L)
      } else {
        args$regtype <- "lc"
      }
      message <- if (bwtype == "fixed") {
        "npindexbw: bandwidth is below the continuous scale-factor lower bound"
      } else {
        "npindexbw: nearest-neighbor bandwidth candidate must map to an integer in [2, 30]"
      }
      expect_error(testthat::with_mocked_bindings(
        do.call(npindexbw, args),
        .npindex_default_start_bandwidth = function(...) {
          default.calls <<- default.calls + 1L
          original.default(...)
        },
        .package = "npRmpi"
      ), message, fixed = TRUE)
      expect_identical(default.calls, 0L)
    }
  }
})
