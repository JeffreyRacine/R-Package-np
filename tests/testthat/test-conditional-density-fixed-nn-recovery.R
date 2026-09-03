test_that("conditional-density MADS recovers only automatic invalid NN starts", {
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x1 = c(rep(0, 32), sin((1:16)*sqrt(2))+(1:16)/32),
    x2 = c(rep(0, 32), sin((1:16)*sqrt(3))+(1:16)/48))
  y <- data.frame(y = sin((1:48)/3)+(1:48)/90)
  for (method in c("cv.ml", "cv.ls")) for (type in c("generalized_nn", "adaptive_nn")) {
    controls <- list(xdat = x, ydat = y, bwtype = type, bwmethod = method,
      regtype = "lp", degree = c(1L,1L), nomad = FALSE, bwsolver = "mads",
      nmulti = 1L, itmax = 20L, powell.remin = FALSE,
      nomad.opts = list(MAX_BB_EVAL = 30L))
    bw <- do.call(npcdensbw, controls)
    restarts <- bw$nomad.restart.results
    expect_length(restarts, 2L)
    expect_true(isTRUE(restarts[[2L]]$recovery))
    expect_identical(restarts[[2L]]$recovery_witness$evaluations, 3L)
    expect_identical(bw$nomad.best.restart, 2L)
    expect_true(restarts[[1L]]$native$compiled_callback_calls > 0L)
    expect_equal(bw$num.feval, 3 + sum(vapply(restarts,
      function(z) z$native$total_num.feval, numeric(1L))), tolerance = 0)
    raw <- .npcdensbw_eval_only(xdat = x, ydat = y, bws = bw, invalid.penalty = "dbmax")
    expect_true(.np_nn_raw_objective_valid(raw$objective))
    expect_equal(bw$fval, raw$objective, tolerance = 2e-12)
    controls$bws <- rep(7, 3)
    controls$nomad.opts <- list(MAX_BB_EVAL = 1L)
    expect_error(do.call(npcdensbw, controls), "did not return a raw-valid solution")
    controls$bws <- rep(if (type == "adaptive_nn") 46 else 47, 3)
    explicit <- do.call(npcdensbw, controls)
    expect_false(any(vapply(explicit$nomad.restart.results,
      function(z) isTRUE(z$recovery), logical(1L))))
    expect_true(.np_nn_raw_objective_valid(explicit$fval))
  }
})
