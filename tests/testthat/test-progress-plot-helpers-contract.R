with_np_bindings <- function(bindings, code) {
  code <- substitute(code)
  ns <- asNamespace("np")
  old <- lapply(names(bindings), function(name) get(name, envir = ns, inherits = FALSE))
  names(old) <- names(bindings)

  for (name in names(bindings)) {
    was_locked <- bindingIsLocked(name, ns)
    if (was_locked) {
      unlockBinding(name, ns)
    }
    assign(name, bindings[[name]], envir = ns)
    if (was_locked) {
      lockBinding(name, ns)
    }
  }

  on.exit({
    for (name in names(old)) {
      was_locked <- bindingIsLocked(name, ns)
      if (was_locked) {
        unlockBinding(name, ns)
      }
      assign(name, old[[name]], envir = ns)
      if (was_locked) {
        lockBinding(name, ns)
      }
    }
  }, add = TRUE)

  eval(code, envir = parent.frame())
}

capture_messages_only <- function(expr) {
  messages <- character()
  withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  messages
}

normalize_messages <- function(x) {
  sub("\n$", "", x)
}

progress_time_counter <- function(start = 0, by = 0.6) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

test_that("plot helper progress emits append-only bounded messages", {
  begin <- getFromNamespace(".np_plot_progress_begin", "np")
  tick <- getFromNamespace(".np_plot_progress_tick", "np")
  finish <- getFromNamespace(".np_plot_progress_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 3
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      state <- begin(total = 12, label = "Plot bootstrap wild")
      for (i in seq_len(12L)) {
        state <- tick(state, done = i)
      }
      finish(state)
    },
    now = progress_time_counter()
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")

  expect_identical(
    vapply(actual$trace, `[[`, character(1L), "event"),
    c(rep("render", 13L), "finish")
  )
  expect_identical(lines[[1L]], "[np] Plot bootstrap wild...")
  expect_true(any(grepl("^\\[np\\] Plot bootstrap wild 1/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Plot bootstrap wild 3/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Plot bootstrap wild 6/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Plot bootstrap wild 9/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Plot bootstrap wild 12/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_equal(length(lines), 14L)
})

test_that("plot helper stays silent for instant runs below start grace", {
  begin <- getFromNamespace(".np_plot_progress_begin", "np")
  finish <- getFromNamespace(".np_plot_progress_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 1,
    np.plot.progress.max.intermediate = 3
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      state <- begin(total = 5, label = "Plot bootstrap wild")
      finish(state)
    },
    now = progress_time_counter(start = 0, by = 0.2)
  )

  expect_length(actual$trace, 0)
})

test_that("bootstrap execution stage surfaces immediately on begin", {
  begin <- getFromNamespace(".np_plot_bootstrap_progress_begin", "np")
  finish <- getFromNamespace(".np_plot_progress_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 2,
    np.plot.progress.start.grace.sec = 0.75
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      state <- begin(total = 12L, label = "Plot bootstrap (surf 1/1)")
      finish(state)
    },
    now = progress_time_values(c(0, 0.2))
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")
  expect_identical(
    lines[[1L]],
    "[np] Plot bootstrap (surf 1/1) 0/12 (0.0%, elapsed 0.0s, eta 0.0s)"
  )
  expect_identical(vapply(actual$trace, `[[`, character(1L), "event")[[1L]], "render")
})

test_that("single-surface plot targets omit 1/1 suffixes", {
  reg_label <- getFromNamespace(".np_plot_regression_bootstrap_target_label", "np")
  idx_label <- getFromNamespace(".np_plot_singleindex_bootstrap_target_label", "np")
  stage_label <- getFromNamespace(".np_plot_bootstrap_stage_label", "np")

  expect_null(reg_label(list(ndim = 1L, xnames = list("x")), slice.index = 1L))
  expect_null(idx_label(FALSE))
  expect_identical(
    stage_label("Plot bootstrap", method_label = "inid", target_label = idx_label(FALSE)),
    "Plot bootstrap inid"
  )
})

test_that("multi-surface plot targets retain i/n suffixes", {
  reg_label <- getFromNamespace(".np_plot_regression_bootstrap_target_label", "np")
  scoef_label <- getFromNamespace(".np_plot_scoef_bootstrap_target_label", "np")

  expect_identical(
    reg_label(list(ndim = 3L, xnames = list("x", "z", "w")), slice.index = 2L),
    "z 2/3"
  )
  expect_identical(
    scoef_label(list(xndim = 1L, zndim = 1L, xnames = list("x"), znames = list("z")), slice.index = 2L),
    "z 2/2"
  )
})

test_that("plot helper activity renders immediately for long blocking work", {
  begin <- getFromNamespace(".np_plot_activity_begin", "np")
  finish <- getFromNamespace(".np_plot_activity_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.start.grace.sec = 0.75
  )
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    {
      activity <- begin("Constructing bootstrap bands")
      finish(activity)
    },
    force_renderer = "legacy",
    now = progress_time_values(c(0, 1.0))
  )

  actual <- capture_progress_shadow_trace(
    {
      activity <- begin("Constructing bootstrap bands")
      finish(activity)
    },
    now = progress_time_values(c(0, 1.0))
  )

  actual_render_lines <- vapply(
    actual$trace[vapply(actual$trace, `[[`, character(1L), "event") == "render"],
    `[[`,
    character(1L),
    "line"
  )
  legacy_render_lines <- vapply(
    legacy$trace[vapply(legacy$trace, `[[`, character(1L), "event") == "render"],
    `[[`,
    character(1L),
    "line"
  )

  expect_identical(
    actual_render_lines,
    legacy_render_lines
  )
  expect_identical(
    actual_render_lines,
    "[np] Constructing bootstrap bands... elapsed 0.0s"
  )
  expect_identical(vapply(actual$trace, `[[`, character(1L), "event"), c("render", "finish"))
})

test_that("plot helper activity no longer waits for grace before first render", {
  begin <- getFromNamespace(".np_plot_activity_begin", "np")
  finish <- getFromNamespace(".np_plot_activity_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.start.grace.sec = 0.75
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      activity <- begin("Constructing bootstrap bands")
      finish(activity)
    },
    now = progress_time_values(c(0, 0.2))
  )

  expect_identical(vapply(actual$trace, `[[`, character(1L), "event"), c("render", "finish"))
  expect_identical(
    vapply(actual$trace, `[[`, character(1L), "line"),
    rep("[np] Constructing bootstrap bands... elapsed 0.0s", 2L)
  )
})

test_that("first-render helper emits exactly one render activity line", {
  init <- getFromNamespace(".np_plot_first_render_state", "np")
  begin <- getFromNamespace(".np_plot_first_render_begin", "np")
  finish <- getFromNamespace(".np_plot_first_render_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.start.grace.sec = 0.75
  )
  on.exit(options(old_opts), add = TRUE)

  state <- init()
  actual <- capture_progress_shadow_trace(
    {
      begin(state)
      finish(state)
      begin(state)
      finish(state)
    },
    now = progress_time_values(c(0, 0.2))
  )

  expect_false(isTRUE(state$pending))
  expect_identical(vapply(actual$trace, `[[`, character(1L), "event"), c("render", "finish"))
  expect_identical(
    vapply(actual$trace, `[[`, character(1L), "line"),
    rep("[np] Rendering plot surface... elapsed 0.0s", 2L)
  )
})

test_that("plot engine setup does not open a graphics device on the null device", {
  capture.par <- getFromNamespace(".np_plot_capture_par", "np")
  engine.begin <- getFromNamespace(".np_plot_engine_begin", "np")

  skip_if_not(isTRUE(unname(as.integer(dev.cur())) == 1L))

  before <- unname(as.integer(dev.cur()))
  on.exit({
    while (!isTRUE(unname(as.integer(dev.cur())) == 1L)) {
      dev.off()
    }
  }, add = TRUE)

  expect_identical(capture.par(c("mfrow", "cex")), list())
  expect_identical(unname(as.integer(dev.cur())), before)

  state <- engine.begin(plot.par.mfrow = TRUE)
  expect_true(is.list(state))
  expect_identical(unname(as.integer(dev.cur())), before)
})

test_that("plot helper activity yields the line to nested bounded progress", {
  activity.begin <- getFromNamespace(".np_plot_activity_begin", "np")
  activity.end <- getFromNamespace(".np_plot_activity_end", "np")
  progress.begin <- getFromNamespace(".np_plot_progress_begin", "np")
  progress.tick <- getFromNamespace(".np_plot_progress_tick", "np")
  progress.end <- getFromNamespace(".np_plot_progress_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 3
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      activity <- activity.begin("Preparing plot bootstrap inid")
      state <- progress.begin(total = 4L, label = "Plot bootstrap inid")
      for (i in seq_len(4L)) {
        state <- progress.tick(state, done = i)
      }
      progress.end(state)
      activity.end(activity)
    },
    now = progress_time_counter(start = 0, by = 0.6)
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")

  expect_true(any(grepl("^\\[np\\] Preparing plot bootstrap inid\\.\\.\\. elapsed 0\\.0s$", lines)))
  expect_true(any(grepl("^\\[np\\] Plot bootstrap inid 1/4 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Plot bootstrap inid 4/4 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
})

test_that("block-style plot bootstrap chunking is capped only when plot progress is active", {
  chunk_size <- getFromNamespace(".np_inid_chunk_size", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.max.intermediate = 3L,
    np.plot.inid.chunk.size = NULL
  )
  on.exit(options(old_opts), add = TRUE)

  expect_identical(chunk_size(n = 1000L, B = 9999L, progress_cap = FALSE), 8388L)
  expect_identical(chunk_size(n = 1000L, B = 9999L, progress_cap = TRUE), 8388L)

  capped <- with_np_progress_bindings(
    list(.np_progress_is_interactive = function() TRUE),
    chunk_size(n = 1000L, B = 9999L, progress_cap = TRUE)
  )
  expect_true(capped < 2500L)
  expect_true(capped > 0L)
})

test_that("ordinary inid plot bootstrap chunking warms up early when plot progress is active", {
  chunk_size <- getFromNamespace(".np_inid_chunk_size", "np")
  warmup_max <- getFromNamespace(".np_plot_progress_warmup_max_reps", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.max.intermediate = 3L,
    np.plot.progress.warmup.max.reps = 16L,
    np.plot.inid.chunk.size = NULL
  )
  on.exit(options(old_opts), add = TRUE)

  uncapped <- chunk_size(n = 500L, B = 999L, progress_cap = FALSE, progress_enabled = FALSE)
  expect_identical(uncapped, 999L)

  capped <- chunk_size(n = 500L, B = 999L, progress_cap = FALSE, progress_enabled = TRUE)
  expect_true(capped > 0L)
  expect_true(capped <= warmup_max())
})

test_that("ordinary inid plot bootstrap emits intermediate progress updates", {
  helper <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 3L,
    np.plot.progress.warmup.max.reps = 16L,
    np.plot.inid.chunk.size = NULL
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    helper(
      H = diag(4L),
      ydat = c(1, 2, 3, 4),
      B = 9L
    ),
    now = progress_time_counter()
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")
  expect_true(any(grepl("^\\[np\\] Plot bootstrap inid 3/9 \\(", lines)))
  expect_true(any(grepl("^\\[np\\] Plot bootstrap inid 6/9 \\(", lines)))
  expect_true(any(grepl("^\\[np\\] Plot bootstrap inid 9/9 \\(", lines)))
})

test_that("block-style unconditional bootstrap emits intermediate progress updates", {
  helper <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")
  counts.drawer <- function(start, stopi) {
    matrix(1L, nrow = 4L, ncol = stopi - start + 1L)
  }

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 3L,
    np.plot.inid.chunk.size = NULL
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    helper(
      H = diag(4L),
      ydat = c(1, 2, 3, 4),
      B = 9L,
      counts.drawer = counts.drawer
    ),
    now = progress_time_counter()
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")
  expect_true(any(grepl("^\\[np\\] Plot bootstrap block 3/9 \\(", lines)))
  expect_true(any(grepl("^\\[np\\] Plot bootstrap block 6/9 \\(", lines)))
  expect_true(any(grepl("^\\[np\\] Plot bootstrap block 9/9 \\(", lines)))
})

test_that("block bootstrap drawer uses iid fast path when block length is one", {
  drawer_factory <- getFromNamespace(".np_block_counts_drawer", "np")
  boot.ns <- asNamespace("boot")
  calls <- 0L

  trace(
    what = "ts.array",
    where = boot.ns,
    tracer = quote(calls <<- calls + 1L),
    print = FALSE
  )
  on.exit(untrace("ts.array", where = boot.ns), add = TRUE)

  set.seed(20260313)
  drawer <- drawer_factory(n = 8L, B = 11L, blocklen = 1L, sim = "geom", n.sim = 8L)
  out <- drawer(1L, 3L)

  expect_identical(calls, 0L)
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(8L, 3L))
  expect_true(all(colSums(out) == 8L))
  expect_true(all(out >= 0))
})

test_that("block bootstrap drawer defers ts.array setup until chunk demand", {
  drawer_factory <- getFromNamespace(".np_block_counts_drawer", "np")
  boot.ns <- asNamespace("boot")
  assign(".np_test_ts_array_calls", 0L, envir = .GlobalEnv)

  trace(
    what = "ts.array",
    where = boot.ns,
    tracer = quote(assign(
      ".np_test_ts_array_calls",
      get(".np_test_ts_array_calls", envir = .GlobalEnv) + 1L,
      envir = .GlobalEnv
    )),
    print = FALSE
  )
  on.exit({
    untrace("ts.array", where = boot.ns)
    rm(".np_test_ts_array_calls", envir = .GlobalEnv)
  }, add = TRUE)

  set.seed(20260313)
  drawer <- drawer_factory(n = 12L, B = 9L, blocklen = 3L, sim = "geom", n.sim = 12L)

  expect_identical(get(".np_test_ts_array_calls", envir = .GlobalEnv), 0L)

  out1 <- drawer(1L, 2L)
  expect_identical(get(".np_test_ts_array_calls", envir = .GlobalEnv), 1L)
  expect_identical(dim(out1), c(12L, 2L))
  expect_true(all(colSums(out1) == 12L))

  out2 <- drawer(3L, 4L)
  expect_identical(get(".np_test_ts_array_calls", envir = .GlobalEnv), 2L)
  expect_identical(dim(out2), c(12L, 2L))
  expect_true(all(colSums(out2) == 12L))
})

test_that("plot progress chunk controller adapts chunk size toward throttle interval", {
  make_controller <- getFromNamespace(".np_plot_progress_chunk_controller", "np")
  observe_controller <- getFromNamespace(".np_plot_progress_chunk_observe", "np")

  old_opts <- options(np.messages = TRUE, np.plot.progress = TRUE)
  on.exit(options(old_opts), add = TRUE)

  controller <- make_controller(chunk.size = 1000L, progress = list(throttle_sec = 0.5))
  controller$adaptive <- TRUE
  slower <- observe_controller(controller = controller, bsz = 1000L, elapsed.sec = 10)
  expect_true(slower$chunk.size < 1000L)
  expect_true(slower$chunk.size >= 250L)

  slower$adaptive <- TRUE
  faster <- observe_controller(controller = slower, bsz = slower$chunk.size, elapsed.sec = 0.05)
  expect_true(faster$chunk.size > slower$chunk.size)
})

test_that("heavy plot helpers invoke delayed activity notifications", {
  regression.eval <- getFromNamespace(".np_plot_regression_eval", "np")
  unconditional.eval <- getFromNamespace(".np_plot_unconditional_eval", "np")
  singleindex.local <- getFromNamespace(".np_plot_singleindex_local_eval", "np")
  singleindex.asym <- getFromNamespace(".np_plot_singleindex_asymptotic_eval", "np")

  labels <- character()
  record_activity <- function(label) {
    labels <<- c(labels, label)
    TRUE
  }

  dummy_bws <- list(
    type = "fixed",
    xdati = list(icon = TRUE),
    beta = 1,
    bw = 1,
    ckertype = "gaussian",
    ckerorder = 2L,
    xbw = 1,
    ybw = 1,
    ixcon = TRUE,
    ixord = FALSE,
    ixuno = FALSE,
    iycon = TRUE,
    iyord = FALSE,
    iyuno = FALSE,
    xncon = 1L,
    xnord = 0L,
    xnuno = 0L,
    yncon = 1L,
    ynord = 0L,
    ynuno = 0L,
    xndim = 1L,
    ydati = list(all.dlev = list()),
    total.time = 0,
    timing = FALSE,
    xmcv = structure(1, num.row = 1L),
    nconfac = 1,
    ncatfac = 1,
    sdev = 1,
    cxkertype = "gaussian",
    cxkerorder = 2L,
    cykertype = "gaussian",
    cykerorder = 2L,
    cxkerlb = -Inf,
    cxkerub = Inf,
    cykerlb = -Inf,
    cykerub = Inf,
    uxkertype = "aitchisonaitken",
    uykertype = "aitchisonaitken",
    oxkertype = "wangvanryzin",
    oykertype = "wangvanryzin"
  )
  x <- data.frame(x = 1:2)

  with_np_bindings(
    list(
      .np_plot_activity_begin = record_activity,
      .np_plot_activity_end = function(state) invisible(NULL),
      .np_regression_direct = function(...) list(mean = c(1, 2), grad = matrix(c(1, 2), ncol = 1L)),
      .np_ksum_unconditional_eval_exact = function(...) c(1, 2),
      npreg = function(...) list(mean = c(1, 2), merr = c(0.1, 0.2), grad = matrix(c(1, 2), ncol = 1L), gerr = matrix(c(0.1, 0.2), ncol = 1L)),
      npudens = function(...) list(dens = c(1, 2), derr = c(NA_real_, NA_real_)),
      npudist = function(...) list(dist = c(1, 2), derr = c(NA_real_, NA_real_)),
      .np_indexhat_rbw = function(...) dummy_bws,
      .np_plot_singleindex_hat_apply_index = function(...) c(1, 2),
      .npindex_resolve_spec = function(...) list(regtype.engine = "lc"),
      npValidateScalarLogical = function(x, ...) isTRUE(x),
      adjustLevels = function(dat, ...) dat,
      toFrame = function(x) if (is.data.frame(x)) x else data.frame(x = x),
      toMatrix = function(x) as.matrix(x)
    ),
    {
      regression.eval(bws = dummy_bws, xdat = x, ydat = c(1, 2), exdat = x)
      unconditional.eval(xdat = x, exdat = x, bws = dummy_bws, cdf = FALSE)
      singleindex.local(bws = dummy_bws, idx.train = x, idx.eval = x, ydat = c(1, 2))
      singleindex.asym(bws = dummy_bws, txdat = x, tydat = c(1, 2), index.eval = c(1, 2))
    }
  )

  expect_true("Computing regression plot fit" %in% labels)
  expect_true("Computing unconditional density plot fit" %in% labels)
  expect_true("Computing single-index plot fit" %in% labels)
  expect_true("Computing single-index plot asymptotic fit" %in% labels)
})

test_that("plot helper progress is silent by default in noninteractive mode", {
  begin <- getFromNamespace(".np_plot_progress_begin", "np")
  tick <- getFromNamespace(".np_plot_progress_tick", "np")
  finish <- getFromNamespace(".np_plot_progress_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 3
  )
  on.exit(options(old_opts), add = TRUE)

  messages <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() FALSE,
      .np_progress_now = progress_time_counter()
    ),
    capture_messages_only({
      state <- begin(total = 12, label = "Plot bootstrap wild")
      state <- tick(state, done = 12)
      finish(state)
    })
  )

  expect_length(messages, 0)
})

test_that("plot helper progress respects suppressMessages", {
  begin <- getFromNamespace(".np_plot_progress_begin", "np")
  tick <- getFromNamespace(".np_plot_progress_tick", "np")
  finish <- getFromNamespace(".np_plot_progress_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 3
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    suppressMessages({
        state <- begin(total = 12, label = "Plot bootstrap wild")
        state <- tick(state, done = 12)
        finish(state)
      }),
    now = progress_time_counter()
  )

  expect_length(actual$trace, 0)
})

test_that("plot helper progress caps intermediate heartbeats", {
  begin <- getFromNamespace(".np_plot_progress_begin", "np")
  tick <- getFromNamespace(".np_plot_progress_tick", "np")
  finish <- getFromNamespace(".np_plot_progress_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 10,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 2
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      state <- begin(total = 12, label = "Plot bootstrap wild")
      for (i in seq_len(12L)) {
        state <- tick(state, done = i)
      }
      finish(state)
    },
    now = progress_time_counter()
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")

  expect_equal(length(lines), 2L)
  expect_identical(lines[[1L]], "[np] Plot bootstrap wild...")
  expect_true(any(grepl("^\\[np\\] Plot bootstrap wild 12/12 ", lines)))
})

test_that("rotation helper emits bounded progress with eta", {
  begin <- getFromNamespace(".np_plot_rotation_progress_begin", "np")
  tick <- getFromNamespace(".np_plot_rotation_progress_tick", "np")
  finish <- getFromNamespace(".np_plot_rotation_progress_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 3
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      state <- begin(total_frames = 8L)
      for (i in seq_len(8L)) {
        state <- tick(state, done = i)
      }
      finish(state)
    },
    now = progress_time_counter()
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")
  expect_identical(lines[[1L]], "[np] Rotating plot...")
  expect_true(any(grepl("^\\[np\\] Rotating plot 1/8 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Rotating plot 8/8 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
})
