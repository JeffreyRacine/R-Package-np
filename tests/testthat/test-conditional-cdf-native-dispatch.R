test_that("native conditional-CDF calls join the existing rank-symmetric owner", {
  skip_on_cran()
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  i <- seq_len(16L)
  x <- data.frame(x = sin(i * sqrt(2)) + i / 32)
  y <- data.frame(y = sin(i / 3) + i / 90)
  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    for (regtype in c("lc", "lp")) {
      k <- if (type == "fixed") 0.8 else 10
      bw.args <- list(xdat = x, ydat = y, bws = c(k, k),
        bandwidth.compute = FALSE, bwscaling = TRUE, bwtype = type,
        regtype = regtype, nomad = FALSE)
      if (regtype == "lp") bw.args$degree <- 1L
      bw <- do.call(npcdistbw, bw.args)
      prep <- .npcdistbw_nomad_native_prepare_args(
        xdat = x, ydat = y, bws = bw, ngrid = 7L)
      controls <- list(prep = prep, x0 = c(k, k),
        bbin = rep.int(if (type == "fixed") 0L else 1L, 2L),
        lb = rep.int(if (type == "fixed") 0.1 else 1, 2L),
        ub = rep.int(if (type == "fixed") 2 else 14, 2L),
        max.eval = 1L, random.seed = 42L, inner.start.count = 0L,
        option.names = character(), option.values = character())
      command <- substitute(
        do.call(get("npNomadNativeSearchConditionalDistribution",
          envir = asNamespace("npRmpi"), inherits = FALSE), ARGS),
        list(ARGS = controls))
      .npRmpi_sync_prepared_search_options()
      reference <- .npRmpi_bcast_cmd_expr(command, caller.execute = TRUE)
      # A master-owned call must refresh worker options, not inherit stale
      # collective eligibility from an unrelated previous command.
      mpi.bcast.cmd(options(np.tree = TRUE, np.extendednn = TRUE),
        caller.execute = FALSE)
      observed <- do.call(npNomadNativeSearchConditionalDistribution, controls)
      expect_identical(observed, reference, info = paste(type, regtype))
      expect_identical(observed$compiled_callback_failures, 0L)
    }
  }
})
