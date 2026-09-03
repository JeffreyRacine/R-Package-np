test_that("prepared option synchronization preserves unset cache defaults", {
  old <- options(np.objective.cache = NULL)
  on.exit(options(old), add = TRUE)
  captured <- NULL
  local_mocked_bindings(.npRmpi_bcast_cmd_expr = function(expr, comm, caller.execute) {
    captured <<- list(expr = expr, comm = comm, caller.execute = caller.execute)
    invisible(NULL)
  }, .package = "npRmpi")
  expect_null(.npRmpi_sync_prepared_search_options(comm = 2L))
  expect_identical(captured$comm, 2L)
  expect_false(captured$caller.execute)
  expect_null(eval(captured$expr)[["np.objective.cache"]])
  captured <- NULL
  options(np.objective.cache = "invalid")
  expect_error(.npRmpi_sync_prepared_search_options(),
    "option 'np.objective.cache' must be TRUE or FALSE", fixed = TRUE)
  expect_null(captured)
})

test_that("manual prepared degree searches receive current master options", {
  skip_on_cran()
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE,
    np.extendednn = FALSE, np.largeh = TRUE, np.objective.cache = TRUE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = sin((1:18) * sqrt(2)) + (1:18) / 30)
  y <- sin((1:18) / 3) + (1:18) / 90
  for (family in c("npregbw", "npcdensbw", "npcdistbw")) {
    for (extended in c(FALSE, TRUE)) {
      options(np.extendednn = extended, np.largeh = !extended,
        np.tree = extended, np.objective.cache = !extended)
      expected <- .npRmpi_autodispatch_option_snapshot()
      controls <- list(xdat = x,
        ydat = if (family == "npregbw") y else data.frame(y = y),
        bws = rep(14, if (family == "npregbw") 1L else 2L),
        bwtype = if (extended) "adaptive_nn" else "generalized_nn",
        regtype = "lp", nomad = TRUE, search.engine = "nomad",
        degree.select = "coordinate", degree.min = 0L, degree.max = 1L,
        degree.start = 1L, degree.verify = FALSE, nmulti = 1L, itmax = 10L,
        nomad.opts = list(MAX_BB_EVAL = 8L))
      if (family == "npcdistbw") controls$ngrid <- 10L
      bw <- do.call(get(family, envir = asNamespace("npRmpi")), controls)
      expect_true(.np_nn_raw_objective_valid(bw$fval))
      actual <- mpi.bcast.cmd(npRmpi:::mpi.allgather.Robj(
        npRmpi:::.npRmpi_autodispatch_option_snapshot()), caller.execute = TRUE)
      expect_true(all(vapply(seq_len(ncol(actual)), function(rank)
        identical(as.list(actual[, rank]), expected), logical(1L))))
    }
  }
})
