test_that("automatic degree search resumes from a raw-valid ordinary NN point", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  x <- as.data.frame(lapply(1:2, function(j) c(rep(0, 32L),
    sin((1:16) * sqrt(j + 1)) + (1:16) / (16 * (j + 1)))))
  y <- sin((1:48) / 3) + (1:48) / 90
  for (type in c("generalized_nn", "adaptive_nn")) {
    bw <- do.call(npregbw, list(xdat = x, ydat = y, bwtype = type,
      regtype = "lp", nomad = TRUE, search.engine = "nomad",
      degree.select = "coordinate", degree.min = 0L, degree.max = 2L,
      degree.start = c(1L, 1L), degree.verify = FALSE, nmulti = 1L,
      itmax = 20L, nomad.opts = list(MAX_BB_EVAL = 30L)))
    rs <- bw$nomad.restart.results
    expect_length(rs, 2L)
    expect_identical(rs[[1L]]$start, c(7, 7, 1, 1))
    expect_true(rs[[2L]]$recovery)
    witness <- rs[[2L]]$recovery_witness
    cap <- if (type == "adaptive_nn") 46 else 47
    expect_identical(witness$point, c(cap, cap, 1, 1))
    expect_identical(witness$evaluations, 3L)
    expect_identical(rs[[2L]]$start, witness$point)
    expect_gt(rs[[2L]]$native$compiled_callback_calls, 0)
    raw <- .npregbw_eval_only(x, y, bw, invalid.penalty = "dbmax")$objective
    expect_true(.np_nn_raw_objective_valid(raw))
    expect_equal(as.numeric(bw$fval), as.numeric(raw), tolerance = 2e-12)
    expect_identical(as.numeric(bw$fval),
      as.numeric(rs[[bw$nomad.best.restart]]$native$objective))
    expect_identical(as.numeric(bw$num.feval), sum(vapply(rs,
      function(z) as.numeric(z$native$total_num.feval), numeric(1))) + 3)
  }
})
