test_that("raw certification publishes its score and retains native diagnostics", {
  skip_on_cran()
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = c(rep(0, 20L), 1:4))
  for (family in c("npudensbw", "npudistbw")) for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    args <- list(dat = x, bwtype = type, bwsolver = "mads", nmulti = 1L,
      itmax = 20L, powell.remin = FALSE, nomad.opts = list(MAX_BB_EVAL = 30L))
    if (family == "npudensbw") args$bwmethod <- "cv.ls" else args$ngrid <- 7L
    bw <- do.call(get(family, asNamespace("npRmpi")), args)
    selected <- bw$nomad.restart.results[[bw$nomad.best.restart]]
    certificate <- as.numeric(selected$objective)
    if (family == "npudensbw") certificate <- -certificate
    expect_identical(as.numeric(bw$fval), certificate)
    expect_identical(as.numeric(selected$best_objective),
                     as.numeric(selected$native$objective))
    raw.args <- list(dat = x, bws = bw, eval.only = TRUE, invalid.penalty = "dbmax", nmulti = 1L)
    if (family == "npudistbw") raw.args$ngrid <- 7L
    raw.name <- if (family == "npudensbw") "npudensbw.bandwidth" else "npudistbw.dbandwidth"
    raw <- do.call(get(raw.name, asNamespace("npRmpi")), raw.args)$fval
    expect_true(is.finite(raw) && abs(raw) < .Machine$double.xmax)
    expect_equal(as.numeric(bw$fval), as.numeric(raw), tolerance = 2e-12)
  }
})
