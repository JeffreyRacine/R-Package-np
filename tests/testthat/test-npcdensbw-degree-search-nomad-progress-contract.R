test_that("npcdensbw nomad+powell payload does not inject phantom multistart totals", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260320)
  dat <- data.frame(x = runif(20))
  dat$y <- rbeta(nrow(dat), 1, 1)

  trace_env <- new.env(parent = emptyenv())
  trace_env$totals <- integer()
  trace_env$powell.count <- 0L

  trace(
    npRmpi:::.np_progress_bandwidth_set_total,
    tracer = eval(substitute(
      quote({
        assign(
          "totals",
          c(get("totals", envir = .trace_env, inherits = FALSE), as.integer(total)),
          envir = .trace_env
        )
      }),
      list(.trace_env = trace_env)
    )),
    print = FALSE
  )
  trace(
    npRmpi:::.np_nomad_with_powell_progress,
    tracer = eval(substitute(
      quote({
        assign(
          "powell.count",
          get("powell.count", envir = .trace_env, inherits = FALSE) + 1L,
          envir = .trace_env
        )
      }),
      list(.trace_env = trace_env)
    )),
    print = FALSE
  )
  on.exit(untrace(npRmpi:::.np_progress_bandwidth_set_total), add = TRUE)
  on.exit(untrace(npRmpi:::.np_nomad_with_powell_progress), add = TRUE)

  bw <- npRmpi::npcdensbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad+powell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ml",
    nmulti = 2L,
    max.bb.eval = 20L,
    cxkerbound = "range",
    cykerbound = "range"
  )

  expect_s3_class(bw, "conbandwidth")
  expect_identical(get("powell.count", envir = trace_env, inherits = FALSE), 1L)
  expect_false(any(get("totals", envir = trace_env, inherits = FALSE) == 2L))
})
