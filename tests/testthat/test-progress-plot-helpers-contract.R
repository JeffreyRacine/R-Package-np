with_nprmpi_bindings <- function(bindings, code) {
  code <- substitute(code)
  ns <- asNamespace("npRmpi")
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

test_that("plot helper progress emits append-only bounded messages on master", {
  begin <- getFromNamespace(".np_plot_progress_begin", "npRmpi")
  tick <- getFromNamespace(".np_plot_progress_tick", "npRmpi")
  finish <- getFromNamespace(".np_plot_progress_end", "npRmpi")

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

  events <- vapply(actual$trace, `[[`, character(1L), "event")
  expect_true(events[[length(events)]] %in% c("render", "finish"))
  expect_true(all(events[-length(events)] == "render"))
  expect_identical(lines[[1L]], "[npRmpi] Plot bootstrap wild...")
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 3/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 6/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 9/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 12/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(length(lines) >= 5L)
})

test_that("plot helper stays silent for instant runs below start grace on master", {
  begin <- getFromNamespace(".np_plot_progress_begin", "npRmpi")
  finish <- getFromNamespace(".np_plot_progress_end", "npRmpi")

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

test_that("bootstrap execution stage surfaces immediately on begin in npRmpi", {
  begin <- getFromNamespace(".np_plot_bootstrap_progress_begin", "npRmpi")
  finish <- getFromNamespace(".np_plot_progress_end", "npRmpi")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 2,
    np.plot.progress.start.grace.sec = 0.75
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      state <- begin(total = 12L, label = "Plot bootstrap (index 1/1)")
      finish(state)
    },
    now = progress_time_values(c(0, 0.2))
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")
  expect_identical(
    lines[[1L]],
    "[npRmpi] Plot bootstrap (index 1/1) 0/12 (0.0%, elapsed 0.0s, eta 0.0s)"
  )
  expect_true(vapply(actual$trace, `[[`, character(1L), "event")[[1L]] == "render")
})

test_that("single-surface plot targets omit 1/1 suffixes in npRmpi", {
  reg_label <- getFromNamespace(".np_plot_regression_bootstrap_target_label", "npRmpi")
  idx_label <- getFromNamespace(".np_plot_singleindex_bootstrap_target_label", "npRmpi")
  stage_label <- getFromNamespace(".np_plot_bootstrap_stage_label", "npRmpi")

  expect_null(reg_label(list(ndim = 1L, xnames = list("x")), slice.index = 1L))
  expect_null(idx_label(FALSE))
  expect_identical(
    stage_label("Plot bootstrap", method_label = "inid", target_label = idx_label(FALSE)),
    "Plot bootstrap inid"
  )
})

test_that("multi-surface plot targets retain i/n suffixes in npRmpi", {
  reg_label <- getFromNamespace(".np_plot_regression_bootstrap_target_label", "npRmpi")
  scoef_label <- getFromNamespace(".np_plot_scoef_bootstrap_target_label", "npRmpi")

  expect_identical(
    reg_label(list(ndim = 3L, xnames = list("x", "z", "w")), slice.index = 2L),
    "z 2/3"
  )
  expect_identical(
    scoef_label(list(xndim = 1L, zndim = 1L, xnames = list("x"), znames = list("z")), slice.index = 2L),
    "z 2/2"
  )
})

test_that("plot helper activity delays its note until grace elapses on master", {
  begin <- getFromNamespace(".np_plot_activity_begin", "npRmpi")
  finish <- getFromNamespace(".np_plot_activity_end", "npRmpi")

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
    now = progress_time_values(c(0, 1.0))
  )

  expect_identical(
    vapply(actual$trace, `[[`, character(1L), "line"),
    rep("[npRmpi] Constructing bootstrap bands... elapsed 0.0s", 2L)
  )
  expect_identical(vapply(actual$trace, `[[`, character(1L), "event"),
                   c("render", "finish"))
})

test_that("plot helper activity stays silent below grace on master", {
  begin <- getFromNamespace(".np_plot_activity_begin", "npRmpi")
  finish <- getFromNamespace(".np_plot_activity_end", "npRmpi")

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

  expect_identical(
    vapply(actual$trace, `[[`, character(1L), "line"),
    rep("[npRmpi] Constructing bootstrap bands... elapsed 0.0s", 2L)
  )
  expect_identical(vapply(actual$trace, `[[`, character(1L), "event"),
                   c("render", "finish"))
})

test_that("first-render helper emits exactly one render activity line on master", {
  init <- getFromNamespace(".np_plot_first_render_state", "npRmpi")
  begin <- getFromNamespace(".np_plot_first_render_begin", "npRmpi")
  finish <- getFromNamespace(".np_plot_first_render_end", "npRmpi")

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
  expect_identical(vapply(actual$trace, `[[`, character(1L), "event"),
                   c("render", "finish"))
  expect_identical(
    vapply(actual$trace, `[[`, character(1L), "line"),
    rep("[npRmpi] Rendering plot surface... elapsed 0.0s", 2L)
  )
})

test_that("plot engine setup does not open a graphics device on the null device", {
  capture.par <- getFromNamespace(".np_plot_capture_par", "npRmpi")

  skip_if_not(isTRUE(unname(as.integer(dev.cur())) == 1L))

  before <- unname(as.integer(dev.cur()))
  on.exit({
    while (!isTRUE(unname(as.integer(dev.cur())) == 1L)) {
      dev.off()
    }
  }, add = TRUE)

  expect_identical(capture.par(c("mfrow", "cex")), list())
  expect_identical(unname(as.integer(dev.cur())), before)
})

test_that("plot activity relinquishes single-line ownership before bounded stage begins", {
  begin.activity <- getFromNamespace(".np_plot_activity_begin", "npRmpi")
  end.activity <- getFromNamespace(".np_plot_activity_end", "npRmpi")
  begin.bounded <- getFromNamespace(".np_plot_progress_begin", "npRmpi")
  end.bounded <- getFromNamespace(".np_plot_progress_end", "npRmpi")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      activity <- begin.activity("Computing regression plot fit")
      end.activity(activity)
      bounded <- begin.bounded(total = 12, label = "Plot bootstrap wild")
      bounded <- getFromNamespace(".np_plot_progress_tick", "npRmpi")(bounded, done = 12L)
      end.bounded(bounded)
    },
    now = progress_time_counter()
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")
  expect_true(any(grepl("^\\[npRmpi\\] Computing regression plot fit\\.\\.\\. elapsed 0\\.0s$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild\\.\\.\\.$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 12/12 ", lines)))
})

test_that("block-style plot bootstrap chunking is capped only when master plot progress is active", {
  chunk_size <- getFromNamespace(".np_inid_chunk_size", "npRmpi")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.max.intermediate = 3L,
    np.plot.inid.chunk.size = NULL
  )
  on.exit(options(old_opts), add = TRUE)

  expect_identical(chunk_size(n = 1000L, B = 9999L, progress_cap = FALSE), 8388L)
  expect_identical(chunk_size(n = 1000L, B = 9999L, progress_cap = TRUE), 8388L)

  capped <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE
    ),
    chunk_size(n = 1000L, B = 9999L, progress_cap = TRUE)
  )
  expect_true(capped < 2500L)
  expect_true(capped > 0L)
})

test_that("ordinary inid plot bootstrap chunking warms up early when master plot progress is active", {
  chunk_size <- getFromNamespace(".np_inid_chunk_size", "npRmpi")
  warmup_max <- getFromNamespace(".np_plot_progress_warmup_max_reps", "npRmpi")

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

test_that("ordinary inid plot bootstrap emits intermediate progress updates in npRmpi", {
  helper <- getFromNamespace(".np_inid_lc_boot_from_hat", "npRmpi")
  progress.begin <- getFromNamespace(".np_plot_bootstrap_progress_begin", "npRmpi")
  progress.tick <- getFromNamespace(".np_plot_progress_tick", "npRmpi")
  progress.end <- getFromNamespace(".np_plot_progress_end", "npRmpi")

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
    with_nprmpi_bindings(
      list(
        .np_progress_is_interactive = function() TRUE,
        .np_progress_is_master = function() TRUE,
        .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
        .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, ...) chunk.size,
        .npRmpi_bootstrap_run_fanout = function(tasks, worker, ncol.out, progress.label, ...) {
          out <- matrix(
            NA_real_,
            nrow = sum(vapply(tasks, function(task) as.integer(task$bsz), integer(1L))),
            ncol = ncol.out
          )
          progress <- progress.begin(total = nrow(out), label = progress.label)
          on.exit(progress.end(progress), add = TRUE)
          row <- 1L
          for (task in tasks) {
            chunk <- worker(task)
            bsz <- nrow(chunk)
            out[row:(row + bsz - 1L), ] <- chunk
            row <- row + bsz
            progress <- progress.tick(progress, done = row - 1L)
          }
          out
        }
      ),
      helper(
        H = diag(4L),
        ydat = c(1, 2, 3, 4),
        B = 9L,
        progress.label = "Plot bootstrap inid"
      )
    ),
    now = progress_time_counter()
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap inid 3/9 \\(", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap inid 6/9 \\(", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap inid 9/9 \\(", lines)))
})

test_that("block bootstrap drawer uses iid fast path when block length is one", {
  drawer_factory <- getFromNamespace(".np_block_counts_drawer", "npRmpi")
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

test_that("block bootstrap drawer defers ts.array setup until chunk demand in npRmpi", {
  drawer_factory <- getFromNamespace(".np_block_counts_drawer", "npRmpi")
  boot.ns <- asNamespace("boot")
  assign(".npRmpi_test_ts_array_calls", 0L, envir = .GlobalEnv)

  trace(
    what = "ts.array",
    where = boot.ns,
    tracer = quote(assign(
      ".npRmpi_test_ts_array_calls",
      get(".npRmpi_test_ts_array_calls", envir = .GlobalEnv) + 1L,
      envir = .GlobalEnv
    )),
    print = FALSE
  )
  on.exit({
    untrace("ts.array", where = boot.ns)
    rm(".npRmpi_test_ts_array_calls", envir = .GlobalEnv)
  }, add = TRUE)

  set.seed(20260313)
  drawer <- drawer_factory(n = 12L, B = 9L, blocklen = 3L, sim = "geom", n.sim = 12L)

  expect_identical(get(".npRmpi_test_ts_array_calls", envir = .GlobalEnv), 0L)

  out1 <- drawer(1L, 2L)
  expect_identical(get(".npRmpi_test_ts_array_calls", envir = .GlobalEnv), 1L)
  expect_identical(dim(out1), c(12L, 2L))
  expect_true(all(colSums(out1) == 12L))

  out2 <- drawer(3L, 4L)
  expect_identical(get(".npRmpi_test_ts_array_calls", envir = .GlobalEnv), 2L)
  expect_identical(dim(out2), c(12L, 2L))
  expect_true(all(colSums(out2) == 12L))
})

test_that("plot progress chunk controller adapts chunk size toward throttle interval in npRmpi", {
  make_controller <- getFromNamespace(".np_plot_progress_chunk_controller", "npRmpi")
  observe_controller <- getFromNamespace(".np_plot_progress_chunk_observe", "npRmpi")

  controller <- make_controller(chunk.size = 1000L, progress = list(throttle_sec = 0.5))
  controller$adaptive <- TRUE
  slower <- observe_controller(controller = controller, bsz = 1000L, elapsed.sec = 10)
  expect_true(slower$chunk.size < 1000L)
  expect_true(slower$chunk.size >= 250L)

  slower$adaptive <- TRUE
  faster <- observe_controller(controller = slower, bsz = slower$chunk.size, elapsed.sec = 0.05)
  expect_true(faster$chunk.size > slower$chunk.size)
})

test_that("npRmpi plot progress can emit on elapsed time before checkpoint", {
  begin <- getFromNamespace(".np_plot_progress_begin", "npRmpi")
  tick <- getFromNamespace(".np_plot_progress_tick", "npRmpi")
  finish <- getFromNamespace(".np_plot_progress_end", "npRmpi")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0.5,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 3L
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    with_nprmpi_bindings(
      list(
        .np_progress_is_interactive = function() TRUE,
        .np_progress_is_master = function() TRUE
      ),
      {
        state <- begin(total = 9L, label = "Plot bootstrap block")
        state <- tick(state = state, done = 1L)
        state <- tick(state = state, done = 2L)
        finish(state)
      }
    ),
    now = progress_time_values(c(0, 0.6, 1.2))
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap block 2/9 \\(", lines)))
})

test_that("block-style plot bootstrap uses progress-capped MPI task partitioning on master", {
  helper <- getFromNamespace(".np_inid_lc_boot_from_hat", "npRmpi")
  counts.drawer <- function(start, stopi) {
    matrix(1L, nrow = 4L, ncol = stopi - start + 1L)
  }
  seen.tasks <- list()

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 3L,
    np.plot.inid.chunk.size = NULL
  )
  on.exit(options(old_opts), add = TRUE)

  out <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
      .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, ...) chunk.size,
      .npRmpi_bootstrap_run_fanout = function(tasks, worker, ncol.out, ...) {
        seen.tasks <<- lapply(tasks, function(task) {
          c(start = as.integer(task$start), bsz = as.integer(task$bsz))
        })
        out <- matrix(NA_real_, nrow = sum(vapply(tasks, function(task) as.integer(task$bsz), integer(1L))), ncol = ncol.out)
        row <- 1L
        for (task in tasks) {
          chunk <- worker(task)
          bsz <- nrow(chunk)
          out[row:(row + bsz - 1L), ] <- chunk
          row <- row + bsz
        }
        out
      }
    ),
    {
      helper(
        H = diag(4L),
        ydat = c(1, 2, 3, 4),
        B = 9L,
        counts.drawer = counts.drawer
      )
    }
  )

  expect_identical(seen.tasks, list(
    c(start = 1L, bsz = 3L),
    c(start = 4L, bsz = 3L),
    c(start = 7L, bsz = 3L)
  ))
  expect_true(is.list(out))
  expect_true(is.matrix(out$t))
  expect_identical(dim(out$t), c(9L, 4L))
})

test_that("heavy plot helpers invoke delayed activity notifications", {
  regression.eval <- getFromNamespace(".np_plot_regression_eval", "npRmpi")
  unconditional.eval <- getFromNamespace(".np_plot_unconditional_eval", "npRmpi")
  singleindex.local <- getFromNamespace(".np_plot_singleindex_local_eval", "npRmpi")
  singleindex.asym <- getFromNamespace(".np_plot_singleindex_asymptotic_eval", "npRmpi")

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

  with_nprmpi_bindings(
    list(
      .np_plot_activity_begin = record_activity,
      .np_plot_activity_end = function(state) invisible(NULL),
      .np_plot_with_local_compiled_eval = function(expr) force(expr),
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

test_that("plot helper progress is silent off master", {
  begin <- getFromNamespace(".np_plot_progress_begin", "npRmpi")
  tick <- getFromNamespace(".np_plot_progress_tick", "npRmpi")
  finish <- getFromNamespace(".np_plot_progress_end", "npRmpi")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 3,
    np.plot.progress.noninteractive = TRUE
  )
  on.exit(options(old_opts), add = TRUE)

  messages <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() FALSE,
      .np_progress_is_master = function() FALSE,
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

test_that("plot helper progress supports the explicit noninteractive override on master", {
  begin <- getFromNamespace(".np_plot_progress_begin", "npRmpi")
  tick <- getFromNamespace(".np_plot_progress_tick", "npRmpi")
  finish <- getFromNamespace(".np_plot_progress_end", "npRmpi")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 3,
    np.plot.progress.noninteractive = TRUE
  )
  on.exit(options(old_opts), add = TRUE)

  messages <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() FALSE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_messages_only({
      state <- begin(total = 12, label = "Plot bootstrap wild")
      for (i in seq_len(12L)) {
        state <- tick(state, done = i)
      }
      finish(state)
    })
  )

  messages <- normalize_messages(messages)
  expect_true(length(messages) >= 5L)
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 12/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
})

test_that("plot helper progress respects suppressMessages", {
  begin <- getFromNamespace(".np_plot_progress_begin", "npRmpi")
  tick <- getFromNamespace(".np_plot_progress_tick", "npRmpi")
  finish <- getFromNamespace(".np_plot_progress_end", "npRmpi")

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

test_that("plot helper progress caps intermediate heartbeats on master", {
  begin <- getFromNamespace(".np_plot_progress_begin", "npRmpi")
  tick <- getFromNamespace(".np_plot_progress_tick", "npRmpi")
  finish <- getFromNamespace(".np_plot_progress_end", "npRmpi")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.max.intermediate = 2
  )
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    {
      state <- begin(total = 12, label = "Plot bootstrap wild")
      for (i in seq_len(12L)) {
        state <- tick(state, done = i)
      }
      finish(state)
    },
    force_renderer = "legacy",
    now = progress_time_counter()
  )

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

  events <- vapply(actual$trace, `[[`, character(1L), "event")
  expect_true(events[[length(events)]] %in% c("render", "finish"))
  expect_true(all(events[-length(events)] == "render"))
  expect_true(length(lines) >= 4L)
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 4/12 ", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 8/12 ", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 12/12 ", lines)))
})

test_that("rotation helper emits bounded progress with eta on master", {
  begin <- getFromNamespace(".np_plot_rotation_progress_begin", "npRmpi")
  tick <- getFromNamespace(".np_plot_rotation_progress_tick", "npRmpi")
  finish <- getFromNamespace(".np_plot_rotation_progress_end", "npRmpi")

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
  expect_identical(lines[[1L]], "[npRmpi] Rotating plot...")
  expect_true(any(grepl("^\\[npRmpi\\] Rotating plot 1/8 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Rotating plot 8/8 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
})
