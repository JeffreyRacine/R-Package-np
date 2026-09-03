test_that("fixed-degree regression MADS resumes from an ordinary raw-valid NN start", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = c(rep(0, 20L), 1:4))
  y <- sin((1:24) / 3) + (1:24) / 90
  for (type in c("generalized_nn", "adaptive_nn")) {
    for (regtype in c("lc", "ll", "lp")) for (solver in c("mads", "mads+powell")) {
      args <- list(xdat = x, ydat = y, bwtype = type, regtype = regtype,
        bwsolver = solver, nmulti = 1L, itmax = 20L, powell.remin = FALSE,
        nomad.opts = list(MAX_BB_EVAL = 30L))
      if (regtype == "lp") args$degree <- 2L
      bw <- do.call(npregbw, args)
      raw <- npRmpi:::.npregbw_eval_only(x, y, bw, invalid.penalty = "dbmax")$objective
      restarts <- bw$nomad.restart.results
      expect_true(is.finite(raw) && abs(raw) < .Machine$double.xmax)
      expect_equal(as.numeric(bw$fval), as.numeric(raw), tolerance = 2e-12)
      expect_length(restarts, 2L)
      expect_true(isTRUE(restarts[[2L]]$recovery))
      witness <- restarts[[2L]]$recovery_witness
      expect_true(isTRUE(witness$found))
      expect_identical(restarts[[2L]]$start, witness$point)
      expect_true(all(witness$point <= if (type == "adaptive_nn") 22 else 23))
      expect_gt(restarts[[2L]]$native$compiled_callback_calls, 0)
      if (solver == "mads") expect_identical(as.numeric(bw$fval),
        as.numeric(restarts[[bw$nomad.best.restart]]$native$objective))
    }
    expect_error(do.call(npregbw, list(xdat = x, ydat = y, bws = 5,
      bwtype = type, regtype = "lc", bwsolver = "mads", nmulti = 1L,
      itmax = 20L, nomad.opts = list(MAX_BB_EVAL = 30L))),
      "did not return a raw-valid solution")
  }
})
