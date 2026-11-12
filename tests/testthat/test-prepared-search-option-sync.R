test_that("prepared option synchronization preserves unset cache defaults", {
  saved <- .npRmpi_option_agreement$snapshot
  on.exit(.npRmpi_option_agreement$snapshot <- saved, add = TRUE)
  old <- options(np.objective.cache = NULL)
  on.exit(options(old), add = TRUE)
  captured <- NULL
  local_mocked_bindings(.npRmpi_bcast_cmd_expr = function(expr, comm, caller.execute) {
    expect_null(.npRmpi_option_agreement$snapshot)
    captured <<- list(expr = expr, comm = comm, caller.execute = caller.execute)
    invisible(NULL)
  }, .package = "npRmpi")
  expect_null(.npRmpi_sync_prepared_search_options(comm = 2L))
  expect_identical(.npRmpi_option_agreement$snapshot,
                   .npRmpi_autodispatch_option_snapshot())
  expect_identical(captured$comm, 2L)
  expect_false(captured$caller.execute)
  expect_null(eval(captured$expr)[["np.objective.cache"]])
  captured <- NULL
  options(np.objective.cache = "invalid")
  expect_error(.npRmpi_sync_prepared_search_options(),
    "option 'np.objective.cache' must be TRUE or FALSE", fixed = TRUE)
  expect_null(captured)
})

test_that("failed prepared option delivery leaves agreement unknown", {
  saved <- .npRmpi_option_agreement$snapshot
  on.exit(.npRmpi_option_agreement$snapshot <- saved, add = TRUE)
  .npRmpi_option_agreement$snapshot <- .npRmpi_autodispatch_option_snapshot()
  local_mocked_bindings(.npRmpi_bcast_cmd_expr = function(...) {
    expect_null(.npRmpi_option_agreement$snapshot)
    stop("prepared delivery sentinel")
  }, .package = "npRmpi")
  expect_error(.npRmpi_sync_prepared_search_options(),
               "prepared delivery sentinel", fixed = TRUE)
  expect_null(.npRmpi_option_agreement$snapshot)
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
      expect_identical(.npRmpi_option_agreement$snapshot, expected)
      # A later task must synchronize back, even if a previous certificate
      # described this restored value. Messages must not be needed to heal it.
      options(np.tree = !extended)
      h <- 0.5
      kernel <- outer(x$x, x$x, function(a, b) dnorm((a - b) / h))
      sums <- npksum(txdat = x, tydat = y, bws = h)$ksum
      expect_equal(as.numeric(sums), as.numeric(kernel %*% y), tolerance = 1e-12)
      workers <- mpi.remote.exec(getOption("np.tree"), simplify = FALSE, ret = TRUE)
      expect_true(all(vapply(workers, identical, logical(1L), !extended)))
    }
  }
})


test_that("SPMD sender failure does not publish an option certificate", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = TRUE,
                 npRmpi.autodispatch.option.sync = "onchange")
  on.exit(options(old), add = TRUE)
  x <- seq(-1, 1, length.out = 25L)
  baseline <- npksum(txdat = x, tydat = sin(x), bws = .4)$ksum
  ns <- asNamespace("npRmpi")
  probe <- new.env(parent = ns)
  probe$dispatch <- .npRmpi_distributed_call_impl
  environment(probe$dispatch) <- probe
  # Isolate only the existing sender seam; real preparation and cleanup still
  # use the installed owners. Never replace transport in the namespace.
  reached <- FALSE
  probe$.npRmpi_bcast_cmd_expr <- function(expr, comm, caller.execute) {
    expect_null(.npRmpi_option_agreement$snapshot)
    reached <<- TRUE
    stop("pre-delivery sentinel")
  }
  options(np.tree = FALSE)
  expect_error(probe$dispatch(quote(npksum(txdat = x, tydat = sin(x), bws = .4)),
                  caller_env = environment(), automatic = TRUE),
               "pre-delivery sentinel", fixed = TRUE)
  expect_true(reached)
  expect_null(.npRmpi_option_agreement$snapshot)
  expect_equal(npksum(txdat = x, tydat = sin(x), bws = .4)$ksum,
               baseline, tolerance = 1e-12)
  expect_identical(.npRmpi_option_agreement$snapshot,
                   .npRmpi_autodispatch_option_snapshot())
})
