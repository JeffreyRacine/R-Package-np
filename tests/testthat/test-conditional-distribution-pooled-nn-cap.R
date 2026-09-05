test_that("pooled CDF degree searches preserve the solver's extended NN cap", {
  skip_on_cran()
  withr::local_preserve_seed()
  expect_true(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = TRUE)
  on.exit(options(old), add = TRUE)
  observed <- new.env(parent = emptyenv())
  original <- .npRmpi_bcast_cmd_expr
  testthat::local_mocked_bindings(
    .npRmpi_bcast_cmd_expr = function(expr, ...) {
      if (is.call(expr) && is.call(expr[[1L]]) &&
          identical(expr[[1L]][[1L]], as.name("get")) &&
          identical(expr[[1L]][[2L]], "npRmpiPreparedSearchConditionalDistribution")) {
        observed$template <- expr[[4L]]
        observed$setup <- expr[[5L]]
        observed$upper <- expr[[12L]]
      }
      original(expr, ...)
    }, .package = "npRmpi")
  i <- seq_len(24L)
  x <- data.frame(x = sin(i * sqrt(2)) + i / 24,
                  u = factor(rep(c("a", "b"), 12L)))
  y <- x$x + .3 * (x$u == "b") + sin(i * sqrt(3)) / 2
  for (type in c("generalized_nn", "adaptive_nn")) {
    for (extended in c(TRUE, FALSE)) {
      options(np.extendednn = extended)
      observed$setup <- NULL
      bw <- npcdistbw(xdat = x, ydat = y, bwtype = type, regtype = "lp",
        nomad = TRUE, search.engine = "nomad", nmulti = 2L,
        nomad.nmulti = 1L, degree.min = 0L, degree.max = 1L,
        degree.start = 0L, degree.verify = FALSE, random.seed = 42L,
        ngrid = 7L, nomad.opts = list(MAX_BB_EVAL = 4L))
      expect_type(observed$setup, "list")
      cap <- observed$upper[1:2]
      expect_identical(as.numeric(observed$setup$cont_extendednn_upper),
                       as.numeric(cap))
      expect_identical(as.integer(observed$setup$cont_flat), c(1L, 2L))
      expect_identical(as.integer(observed$setup$cat_flat), 3L)
      expect_true(is.finite(bw$fval))
      # Independent external storage: Y count, X count, X lambda. The expected
      # points come from solver bounds, never from the decoder under test.
      counts <- if (extended) c(22, 23, 24, cap[2L]) else c(22, 23)
      for (count in counts) {
        actual <- .npcdistbw_nomad_point_to_bw(c(8, count, 2000),
          observed$template, observed$setup)
        expect_identical(as.numeric(actual), c(8, count, .2))
      }
      ordinary <- observed$setup
      ordinary$cont_extendednn_upper <- NULL
      expect_identical(as.numeric(.npcdistbw_nomad_point_to_bw(c(8, 80, 2000),
        observed$template, ordinary)), c(8, 23, .2))
    }
  }
})
