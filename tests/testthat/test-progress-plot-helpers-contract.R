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

  messages <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
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

  expect_identical(messages[[1L]], "[npRmpi] Plot bootstrap wild...")
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 3/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 6/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 9/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 12/12 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_equal(length(messages), 5L)
  expect_false(any(grepl(intToUtf8(8L), messages, fixed = TRUE)))
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

  messages <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_counter(start = 0, by = 0.2)
    ),
    capture_messages_only({
      state <- begin(total = 5, label = "Plot bootstrap wild")
      finish(state)
    })
  )

  expect_length(messages, 0)
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

  messages <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_values(c(0, 1.0))
    ),
    capture_messages_only({
      activity <- begin("Constructing bootstrap bands")
      finish(activity)
    })
  )

  messages <- normalize_messages(messages)
  expect_identical(messages, "[npRmpi] Constructing bootstrap bands...")
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

  messages <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_values(c(0, 0.2))
    ),
    capture_messages_only({
      activity <- begin("Constructing bootstrap bands")
      finish(activity)
    })
  )

  expect_length(messages, 0)
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
  expect_equal(length(messages), 5L)
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

  messages <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_messages_only(
      suppressMessages({
        state <- begin(total = 12, label = "Plot bootstrap wild")
        state <- tick(state, done = 12)
        finish(state)
      })
    )
  )

  expect_length(messages, 0)
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

  messages <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
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
  expect_equal(length(messages), 4L)
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 4/12 ", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 8/12 ", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Plot bootstrap wild 12/12 ", messages)))
})
