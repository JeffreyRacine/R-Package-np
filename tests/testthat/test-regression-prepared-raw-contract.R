test_that("prepared regression raw probes preserve the active owner", {
  expect_true(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  .npRmpi_with_local_regression(local({
    x <- data.frame(x = c(rep(0, 20L), 1:4))
    y <- sin((1:24) / 3)
    for (type in c("generalized_nn", "adaptive_nn")) {
      bw <- npregbw(xdat = x, ydat = y, bws = 20, bwtype = type,
        regtype = "lc", bandwidth.compute = FALSE)
      prep <- .npregbw_prepared_args(x, y, bw, invalid.penalty = "baseline")
      owner <- .npregbw_prepared_begin(x, y, bw, broadcast = FALSE)
      on.exit(.npregbw_prepared_end(owner, broadcast = FALSE), add = TRUE)
      degree <- prep$degree
      valid <- npRmpiPreparedObjectiveEvalRegression(20, degree)
      invalid <- npRmpiPreparedObjectiveEvalRegression(5, degree)
      expect_true(.np_nn_raw_objective_valid(valid[[1L]]))
      expect_true(.np_nn_raw_objective_valid(invalid[[1L]])) # finite search penalty
      expect_identical(npRmpiPreparedObjectiveEvalRegressionRaw(5, degree)[[1L]],
        .Machine$double.xmax)
      expect_equal(npRmpiPreparedObjectiveEvalRegressionRaw(20, degree)[[1L]],
        valid[[1L]], tolerance = 2e-12)
      expect_error(npRmpiPreparedObjectiveEvalRegressionRaw(numeric(), degree),
        "bandwidth vector of unexpected length")
      expect_error(npRmpiPreparedObjectiveEvalRegressionRaw(20, integer()),
        "degree vector of unexpected length")
      expect_identical(npRmpiPreparedObjectiveEvalRegression(20, degree), valid)
      expect_identical(npRmpiPreparedObjectiveEvalRegression(5, degree), invalid)
      .npregbw_prepared_end(owner, broadcast = FALSE)
      owner$active <- FALSE
      expect_error(npRmpiPreparedObjectiveEvalRegressionRaw(20, degree),
        "state is not active")
      replay <- .npregbw_eval_only(x, y, bw, invalid.penalty = "dbmax")
      expect_equal(replay$objective, valid[[1L]], tolerance = 2e-12)
    }
  }))
})
