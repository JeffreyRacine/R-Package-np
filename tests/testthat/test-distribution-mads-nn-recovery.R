test_that("distribution MADS recovers an automatic raw-invalid NN plateau", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  for (type in c("generalized_nn", "adaptive_nn")) {
    for (p in 1:4) {
      x <- as.data.frame(lapply(seq_len(p), function(j) c(rep(0, 20L), (1:4)^j)))
      names(x) <- paste0("x", seq_len(p))
      for (solver in if (p == 1L) c("mads", "mads+powell") else "mads") {
        bw <- do.call(npudistbw, list(dat = x, bwtype = type, bwsolver = solver,
          nmulti = 1L, itmax = 20L, powell.remin = FALSE, ngrid = 7L,
          nomad.opts = list(MAX_BB_EVAL = 30L)))
        raw <- np:::npudistbw.dbandwidth(dat = x, bws = bw, eval.only = TRUE,
          invalid.penalty = "dbmax", nmulti = 1L, ngrid = 7L)$fval
        expect_true(is.finite(raw) && abs(raw) < .Machine$double.xmax)
        expect_equal(as.numeric(raw), as.numeric(bw$fval), tolerance = 2e-12)
        restarts <- bw$nomad.restart.results
        expect_length(restarts, 2L)
        expect_true(isTRUE(restarts[[2L]]$recovery))
        witness <- restarts[[2L]]$recovery_witness
        expect_true(isTRUE(witness$found))
        expect_identical(as.numeric(restarts[[2L]]$start), as.numeric(witness$point))
        expect_true(all(witness$point <= if (type == "adaptive_nn") 22 else 23))
        expect_gt(restarts[[2L]]$native$compiled_callback_calls, 0)
      }
    }
  }
})

test_that("distribution MADS does not replace an explicit invalid NN start", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = c(rep(0, 20L), 1:4))
  for (type in c("generalized_nn", "adaptive_nn")) {
    expect_error(do.call(npudistbw, list(dat = x, bws = 5, bwtype = type,
      bwsolver = "mads", nmulti = 1L, powell.remin = FALSE, ngrid = 7L,
      nomad.opts = list(MAX_BB_EVAL = 30L))), "did not return a raw-valid solution")
  }
})
