test_that("npcdensbw nomad+powell payload does not inject phantom multistart totals", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260320)
  dat <- data.frame(x = runif(20))
  dat$y <- rbeta(nrow(dat), 1, 1)

  trace_env <- new.env(parent = emptyenv())
  trace_env$totals <- integer()
  trace_env$powell.count <- 0L

  trace(
    np:::.np_progress_bandwidth_set_total,
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
    np:::.np_nomad_with_powell_progress,
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
  on.exit(untrace(np:::.np_progress_bandwidth_set_total), add = TRUE)
  on.exit(untrace(np:::.np_nomad_with_powell_progress), add = TRUE)

  bw <- np::npcdensbw(
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

test_that("npcdensbw NOMAD plus Powell progress mirrors shared restart detail", {
  skip_if_not_installed("crs")

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.known.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260323)
  dat <- data.frame(x = runif(20))
  dat$y <- rbeta(nrow(dat), 1, 1)

  msgs <- with_np_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "legacy",
      .np_progress_now = degree_progress_time_values(seq(0, 60, by = 0.25))
    ),
    capture_degree_messages_only(
      np::npcdensbw(
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
        max.bb.eval = 12L,
        cxkerbound = "range",
        cykerbound = "range"
      )
    )
  )

  expect_false(any(grepl("nomad\\+powell", msgs, ignore.case = TRUE)))
  expect_false(any(grepl("eval [0-9]+", msgs)))
  expect_false(any(grepl("fval=", msgs, fixed = TRUE)))
  expect_false(any(grepl("%|eta ", msgs)))
  expect_true(any(grepl("^\\[np\\] Selecting degree and bandwidth \\(", msgs)))
  expect_true(any(grepl("^\\[np\\] Refining bandwidth \\(", msgs)))
  expect_true(any(grepl("multistart [12]/2", msgs)))
  expect_true(any(grepl("iteration [0-9]+", msgs)))
  expect_true(any(grepl("iteration [0-9]+ \\([0-9]+\\)", msgs)))
  expect_true(any(grepl("deg \\(", msgs)))
  expect_true(any(grepl("best \\(", msgs)))
})
