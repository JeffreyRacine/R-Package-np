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

capture_warnings_only <- function(expr) {
  warnings <- character()
  withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  warnings
}

normalize_messages <- function(x) {
  sub("\n$", "", x)
}

progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

capture_single_line_output <- function(pkg, bindings, code) {
  path <- tempfile(fileext = ".log")
  con <- file(path, open = "wb")
  closed <- FALSE

  on.exit({
    if (!closed) {
      close(con)
    }
    unlink(path)
  }, add = TRUE)

  bindings$.np_progress_single_line_connection <- function() con

  if (identical(pkg, "np")) {
    with_np_bindings(bindings, code)
  } else {
    stop("unknown package: ", pkg)
  }

  close(con)
  closed <- TRUE
  readChar(path, nchars = file.info(path)$size, useBytes = TRUE)
}

test_that("progress begin returns disabled state when messages are off", {
  begin <- getFromNamespace(".np_progress_begin", "np")

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  state <- with_np_bindings(
    list(.np_progress_is_interactive = function() TRUE,
         .np_progress_now = function() 10),
    begin("Iterating Landweber-Fridman solve")
  )

  expect_false(state$enabled)
  expect_identical(state$label, "Iterating Landweber-Fridman solve")
  expect_null(state$total)
  expect_false(state$known_total)
})

test_that("known-total formatting is deterministic", {
  fmt <- getFromNamespace(".np_progress_format_known_total", "np")
  state <- list(
    enabled = TRUE,
    pkg_prefix = "[np]",
    label = "Iterating Tikhonov solve",
    total = 5,
    known_total = TRUE,
    started = 10,
    last_emit = 9.5,
    throttle_sec = 0.5,
    last_done = 2,
    domain = "general"
  )

  line <- fmt(state, done = 2, detail = "recomputing alpha", now = 14)

  expect_identical(
    line,
    "[np] Iterating Tikhonov solve 2/5 (40.0%, elapsed 4.0s, eta 6.0s): recomputing alpha"
  )
})

test_that("unknown-total formatting is deterministic", {
  fmt <- getFromNamespace(".np_progress_format_unknown_total", "np")
  state <- list(
    enabled = TRUE,
    pkg_prefix = "[np]",
    label = "Iterating Landweber-Fridman solve",
    total = NULL,
    known_total = FALSE,
    started = 10,
    last_emit = 8,
    throttle_sec = 2,
    last_done = 3,
    domain = "general"
  )

  line <- fmt(state, done = 3, detail = "evaluating stopping rule", now = 14)

  expect_identical(
    line,
    "[np] Iterating Landweber-Fridman solve... iteration 3, elapsed 4.0s: evaluating stopping rule"
  )
})

test_that("known-total progress delays start note and avoids duplicate completion", {
  begin <- getFromNamespace(".np_progress_begin", "np")
  step <- getFromNamespace(".np_progress_step", "np")
  finish <- getFromNamespace(".np_progress_end", "np")

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0.75)
  on.exit(options(old_opts), add = TRUE)

  messages <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_values(c(0, 0.2, 0.8, 0.9))
    ),
    capture_messages_only({
      state <- begin("Bootstrap replications", total = 2)
      state <- step(state, done = 1)
      state <- step(state, done = 2)
      finish(state)
    })
  )

  messages <- normalize_messages(messages)

  expect_identical(messages[[1L]], "[np] Bootstrap replications...")
  expect_false(any(grepl("^\\[np\\] Bootstrap replications 1/2 ", messages)))
  expect_equal(sum(grepl("^\\[np\\] Bootstrap replications 2/2 ", messages)), 1L)
})

test_that("unknown-total progress delays start note until grace elapses", {
  begin <- getFromNamespace(".np_progress_begin", "np")
  step <- getFromNamespace(".np_progress_step", "np")

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.unknown.sec = 1.0)
  on.exit(options(old_opts), add = TRUE)

  messages <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_values(c(0, 0.4, 1.2))
    ),
    capture_messages_only({
      state <- begin("Iterating Landweber-Fridman solve")
      state <- step(state, done = 3, detail = "updating residual smoothing")
      step(state, done = 4, detail = "evaluating stopping rule")
    })
  )

  messages <- normalize_messages(messages)

  expect_identical(messages[[1L]], "[np] Iterating Landweber-Fridman solve...")
  expect_identical(
    messages[[2L]],
    "[np] Iterating Landweber-Fridman solve... iteration 4, elapsed 1.2s: evaluating stopping rule"
  )
  expect_length(messages, 2L)
})

test_that("fast progress stays silent below start grace", {
  begin <- getFromNamespace(".np_progress_begin", "np")
  finish <- getFromNamespace(".np_progress_end", "np")

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0.75)
  on.exit(options(old_opts), add = TRUE)

  messages <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_values(c(0, 0.2))
    ),
    capture_messages_only({
      state <- begin("Bootstrap replications", total = 5)
      finish(state)
    })
  )

  expect_length(messages, 0)
})

test_that("suppressMessages suppresses progress notes", {
  note <- getFromNamespace(".np_progress_note", "np")

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  messages <- with_np_bindings(
    list(.np_progress_is_interactive = function() TRUE),
    capture_messages_only(suppressMessages(note("Preparing IV derivative regression")))
  )

  expect_length(messages, 0)
})

test_that("legacy suppression restores np.messages on success and error", {
  suppress_legacy <- getFromNamespace(".np_progress_with_legacy_suppressed", "np")

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  value <- suppress_legacy({
    expect_false(isTRUE(getOption("np.messages")))
    42
  })

  expect_identical(value, 42)
  expect_true(isTRUE(getOption("np.messages")))

  expect_error(
    suppress_legacy({
      stop("boom")
    }),
    "boom"
  )
  expect_true(isTRUE(getOption("np.messages")))
})

test_that("warning helper prefixes package warnings", {
  warn <- getFromNamespace(".np_warning", "np")

  warnings <- capture_warnings_only(warn("kernel order ignored"))

  expect_identical(warnings, "[np] kernel order ignored")
})

test_that("bandwidth selection helper suppresses legacy output and drives bounded updates", {
  select_bw <- getFromNamespace(".np_progress_select_bandwidth", "np")
  set_total <- getFromNamespace(".np_progress_bandwidth_set_total", "np")
  step_bw <- getFromNamespace(".np_progress_bandwidth_multistart_step", "np")

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old_opts), add = TRUE)

  seen <- NULL
  legacy <- capture_progress_shadow_trace(
    {
      value <- select_bw("Selecting regression bandwidth", {
        seen <<- getOption("np.messages")
        set_total(3)
        step_bw(1, 3)
        step_bw(3, 3)
        7
      })
      expect_identical(value, 7)
    },
    force_renderer = "legacy",
    now = progress_time_values(c(0, 1.1, 2.2, 3.3, 4.4))
  )

  seen <- NULL
  actual <- capture_progress_shadow_trace(
    {
      value <- select_bw("Selecting regression bandwidth", {
        seen <<- getOption("np.messages")
        set_total(3)
        step_bw(1, 3)
        step_bw(3, 3)
        7
      })
      expect_identical(value, 7)
    },
    now = progress_time_values(c(0, 1.1, 2.2, 3.3, 4.4))
  )
  lines <- vapply(actual$trace, `[[`, character(1L), "line")
  legacy_lines <- vapply(
    legacy$trace[vapply(legacy$trace, `[[`, character(1L), "event") == "render"],
    `[[`,
    character(1L),
    "line"
  )
  render_lines <- vapply(
    actual$trace[vapply(actual$trace, `[[`, character(1L), "event") == "render"],
    `[[`,
    character(1L),
    "line"
  )

  expect_false(isTRUE(seen))
  expect_equal(render_lines, legacy_lines)
  expect_true(any(grepl("^\\[np\\] Selecting regression bandwidth\\.\\.\\.$", lines)))
  expect_true(any(grepl("^\\[np\\] Selecting regression bandwidth multistart 1/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Selecting regression bandwidth multistart 3/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(isTRUE(getOption("np.messages")))
})

test_that("compiled progress bridge feeds the bounded bandwidth sink", {
  select_bw <- getFromNamespace(".np_progress_select_bandwidth", "np")

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      value <- select_bw("Selecting regression bandwidth", {
        .Call(
          "C_np_progress_signal",
          "bandwidth_multistart_step",
          "bandwidth",
          1L,
          3L,
          PACKAGE = "np"
        )
        .Call(
          "C_np_progress_signal",
          "bandwidth_multistart_step",
          "bandwidth",
          3L,
          3L,
          PACKAGE = "np"
        )
        11
      })
      expect_identical(value, 11)
    },
    now = progress_time_values(c(0, 1.1, 2.2, 3.3))
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")

  expect_true(any(grepl("^\\[np\\] Selecting regression bandwidth\\.\\.\\.$", lines)))
  expect_true(any(grepl("^\\[np\\] Selecting regression bandwidth multistart 1/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Selecting regression bandwidth multistart 3/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
})

test_that("compiled progress bridge feeds unknown-total bandwidth activity updates", {
  select_bw <- getFromNamespace(".np_progress_select_bandwidth", "np")

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.unknown.sec = 0)
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      value <- select_bw("Selecting conditional density bandwidth", {
        .Call(
          "C_np_progress_signal",
          "bandwidth_activity_step",
          "bandwidth",
          64L,
          0L,
          PACKAGE = "np"
        )
        .Call(
          "C_np_progress_signal",
          "bandwidth_activity_step",
          "bandwidth",
          128L,
          0L,
          PACKAGE = "np"
        )
        13
      })
      expect_identical(value, 13)
    },
    now = progress_time_values(c(0, 1.2, 3.4, 5.6))
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")

  expect_true(any(grepl("^\\[np\\] Selecting conditional density bandwidth\\.\\.\\.$", lines)))
  expect_true(any(grepl("^\\[np\\] Selecting conditional density bandwidth\\.\\.\\. iteration 64, elapsed [0-9]+\\.[0-9]s$", lines)))
  expect_true(any(grepl("^\\[np\\] Selecting conditional density bandwidth\\.\\.\\. iteration 128, elapsed [0-9]+\\.[0-9]s$", lines)))
})

test_that("compiled progress bridge feeds fit updates on the bandwidth surface", {
  fit_begin <- getFromNamespace(".np_fit_progress_begin", "np")
  fit_finish <- getFromNamespace(".np_fit_progress_finish", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      fit_begin("Fitting density", total = 3L, handoff = TRUE, detail = "starting")
      .Call(
        "C_np_progress_signal",
        "fit_step",
        "bandwidth",
        1L,
        3L,
        PACKAGE = "np"
      )
      .Call(
        "C_np_progress_signal",
        "fit_step",
        "bandwidth",
        3L,
        3L,
        PACKAGE = "np"
      )
      fit_finish()
    },
    now = progress_time_values(c(0, 1.1, 2.2, 3.3, 4.4))
  )

  lines <- vapply(actual$trace, `[[`, character(1L), "line")

  expect_true(any(grepl(
    "^\\[np\\] Fitting density 0/3 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[np\\] Fitting density 1/3 \\([0-9]+\\.[0-9]%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[np\\] Fitting density 3/3 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )))
})

test_that("single-line fit drops detail before truncating", {
  fit <- getFromNamespace(".np_progress_fit_single_line", "np")

  line <- "[np] Constructing metric entropy by lag 1/2 (50.0%, elapsed 7.2s, eta 7.2s): lag 1"
  fitted <- fit(line, max_width = 78)

  expect_identical(
    fitted,
    "[np] Constructing metric entropy by lag 1/2 (50.0%, elapsed 7.2s, eta 7.2s)"
  )
})

test_that("single-line fit compacts bandwidth labels before ellipsizing", {
  fit <- getFromNamespace(".np_progress_fit_single_line", "np")

  line <- "[np] Selecting regression bandwidth multistart 2/3 (66.7%, elapsed 6.3s, eta 3.2s)"
  fitted <- fit(line, max_width = 80)

  expect_lte(nchar(fitted, type = "width"), 80L)
  expect_identical(
    fitted,
    "[np] Selecting reg bandwidth multistart 2/3 (66.7%, elapsed 6.3s, eta 3.2s)"
  )
})

test_that("single-line fit keeps smooth coefficient bandwidth lines readable at 80 columns", {
  fit <- getFromNamespace(".np_progress_fit_single_line", "np")

  line <- "[np] Selecting smooth coefficient bandwidth multistart 2/2 (100.0%, elapsed 6.3s, eta 0.0s)"
  fitted <- fit(line, max_width = 80)

  expect_lte(nchar(fitted, type = "width"), 80L)
  expect_identical(
    fitted,
    "[np] Selecting smooth coef bw multistart 2/2 (100.0%, elapsed 6.3s, eta 0.0s)"
  )
})

test_that("single-line fit preserves readable bandwidth fields at 80 columns", {
  fit <- getFromNamespace(".np_progress_fit_single_line", "np")

  line <- "[np] Bandwidth selection (multistart 2/2, iteration 31, elapsed 10.6s, 62.4%, eta 6.4s)"
  fitted <- fit(line, max_width = 80)

  expect_lte(nchar(fitted, type = "width"), 80L)
  expect_identical(
    fitted,
    "[np] Bandwidth selection (2/2, iter 31, elapsed 10.6s, 62.4%, eta 6.4s)"
  )
})

test_that("single-line fit preserves readable coordinator bandwidth fields at 80 columns", {
  fit <- getFromNamespace(".np_progress_fit_single_line", "np")

  line <- "[np] Bandwidth selection (y~z, multistart 2/2, iteration 31, elapsed 10.6s, 62.4%, eta 6.4s)"
  fitted <- fit(line, max_width = 80)

  expect_lte(nchar(fitted, type = "width"), 80L)
  expect_identical(
    fitted,
    "[np] Bandwidth selection (y~z, 2/2, iter 31, elapsed 10.6s, 62.4%, eta 6.4s)"
  )
})

test_that("single-line fit preserves restart and cumulative iteration in NOMAD lines", {
  fit <- getFromNamespace(".np_progress_fit_single_line", "np")

  line <- "[np] Bandwidth selection (multistart 2/2, iteration 31 (128), elapsed 10.6s, deg (1), best (0))"
  fitted <- fit(line, max_width = 80)

  expect_lte(nchar(fitted, type = "width"), 80L)
  expect_identical(
    fitted,
    "[np] Bandwidth selection (2/2, iter 31 (128), elapsed 10.6s, deg (1), best (0))"
  )
})

test_that("single-line fit compacts plot bootstrap labels so prep detail survives", {
  fit <- getFromNamespace(".np_progress_fit_single_line", "np")

  line <- "[np] Preparing plot bootstrap inid 0/50 (0.0%, elapsed 0.0s, eta 0.0s): eval 1/50"
  fitted <- fit(line, max_width = 80)

  expect_lte(nchar(fitted, type = "width"), 80L)
  expect_identical(
    fitted,
    "[np] Prep plot inid 0/50 (0.0%, elapsed 0.0s, eta 0.0s): eval 1/50"
  )
})

test_that("single-line fit preserves readable plot bootstrap counters at 80 columns", {
  fit <- getFromNamespace(".np_progress_fit_single_line", "np")

  line <- "[np] Plot bootstrap (grad index 1/1) 12345/80000 (15.4%, elapsed 8.7s, eta 44.9s)"
  fitted <- fit(line, max_width = 80)

  expect_lte(nchar(fitted, type = "width"), 80L)
  expect_identical(
    fitted,
    "[np] Plot bootstrap (grad idx 1/1) 12345/80000 (15.4%, elapsed 8.7s, eta 44.9s)"
  )
})

test_that("single-line fit preserves both ends when truncation is still required", {
  fit <- getFromNamespace(".np_progress_fit_single_line", "np")

  line <- "[np] Bootstrap replications 123/999 (12.3%, elapsed 12.3s, eta 87.7s)"
  fitted <- fit(line, max_width = 40)

  expect_lte(nchar(fitted, type = "width"), 40L)
  expect_true(startsWith(fitted, "[np]"))
  expect_true(endsWith(fitted, "eta 87.7s)"))
  expect_true(grepl("\\.\\.\\.", fitted))
})

test_that("RStudio capability keeps single-line viable at the boundary", {
  capability <- getFromNamespace(".np_progress_capability", "np")

  actual <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_rstudio_console = function() TRUE
    ),
    capability()
  )

  expect_true(isTRUE(actual$interactive))
  expect_true(isTRUE(actual$rstudio))
  expect_true(isTRUE(actual$single_line_viable))
})

test_that("RStudio width budget reserves redraw margin centrally", {
  output_width <- getFromNamespace(".np_progress_output_width", "np")

  actual <- with_np_bindings(list(.np_progress_is_rstudio_console = function() TRUE), {
    old_options <- options(width = 82)
    on.exit(options(old_options), add = TRUE)
    output_width()
  })

  expect_identical(actual, 78L)
})

test_that("terminal width probe takes precedence over stale option width", {
  output_width <- getFromNamespace(".np_progress_output_width", "np")

  actual <- with_np_bindings(
    list(
      .np_progress_is_rstudio_console = function() FALSE,
      .np_progress_terminal_width_probe = function() 132L
    ),
    {
      old_options <- options(width = 80)
      on.exit(options(old_options), add = TRUE)
      output_width()
    }
  )

  expect_identical(actual, 132L)
})

test_that("single-line finish clears the rendered line without leaving a newline", {
  render <- getFromNamespace(".np_progress_render_single_line", "np")
  line <- "[np] Bootstrap replications 2/2 (100.0%, elapsed 1.0s, eta 0.0s)"

  output <- capture_single_line_output(
    "np",
    list(.np_progress_output_width = function() 80L),
    {
      render(list(render_line = line, last_width = 0L), event = "render")
      render(list(render_line = line, last_width = nchar(line, type = "width")), event = "finish")
    }
  )

  expect_identical(
    output,
    paste0("\r", line, "\r", strrep(" ", nchar(line, type = "width")), "\r")
  )
})

test_that("single-line render uses ANSI erase when the console supports it", {
  render <- getFromNamespace(".np_progress_render_single_line", "np")
  line <- "[np] Rotating plot 18/72 (25.0%, elapsed 10.9s, eta 32.6s)"

  output <- capture_single_line_output(
    "np",
    list(.np_progress_single_line_supports_ansi = function(con) TRUE),
    {
      render(list(render_line = line, last_width = 0L), event = "render")
      render(list(render_line = line, last_width = nchar(line, type = "width")), event = "finish")
    }
  )

  expect_identical(output, paste0("\r\033[2K", line, "\r\033[2K\r"))
})

test_that("single-line render clears stale suffix across state handoff in non-ANSI mode", {
  render <- getFromNamespace(".np_progress_render_single_line", "np")
  output_width <- 120L
  long <- "[np] Calling NOMAD (Nonsmooth Optimization by Mesh Adaptive Direct Search)... elapsed 0.0s"
  short <- "[np] NOMAD search..."

  output <- capture_single_line_output(
    "np",
    list(
      .np_progress_output_width = function() output_width,
      .np_progress_single_line_supports_ansi = function(con) FALSE
    ),
    {
      render(list(render_line = long, last_width = 0L), event = "render")
      render(list(render_line = short, last_width = 0L), event = "render")
    }
  )

  expect_identical(
    output,
    paste0(
      "\r", long,
      "\r", short, strrep(" ", nchar(long, type = "width") - nchar(short, type = "width"))
    )
  )
})

test_that("progress end clears single-line output even when the final state already rendered", {
  begin <- getFromNamespace(".np_progress_begin", "np")
  step <- getFromNamespace(".np_progress_step", "np")
  finish <- getFromNamespace(".np_progress_end", "np")
  reset <- getFromNamespace(".np_progress_reset_registry", "np")

  output <- capture_single_line_output(
    "np",
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_rstudio_console = function() FALSE,
      .np_progress_now = progress_time_values(c(0, 0, 1, 1)),
      .np_progress_output_width = function() 120L
    ),
    {
      reset()
      on.exit(reset(), add = TRUE)
      state <- begin("Bootstrap replications", total = 2, surface = "bootstrap")
      state <- step(state, done = 2)
      finish(state)
    }
  )

  line <- "[np] Bootstrap replications 2/2 (100.0%, elapsed 1.0s, eta 0.0s)"
  expect_identical(
    output,
    paste0(
      "\r",
      strrep(" ", nchar(line, type = "width")),
      "\r"
    )
  )
})

test_that("plot progress can repaint on elapsed time before the first checkpoint", {
  begin <- getFromNamespace(".np_plot_progress_begin", "np")
  tick <- getFromNamespace(".np_plot_progress_tick", "np")
  finish <- getFromNamespace(".np_plot_progress_end", "np")
  reset <- getFromNamespace(".np_progress_reset_registry", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress.start.grace.sec = 0,
    np.plot.progress.interval.sec = 0.5
  )
  on.exit(options(old_opts), add = TRUE)

  output <- capture_single_line_output(
    "np",
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_rstudio_console = function() FALSE,
      .np_progress_now = progress_time_values(c(0, 0, 0.6, 1.2)),
      .np_progress_output_width = function() 120L
    ),
    {
      reset()
      on.exit(reset(), add = TRUE)
      state <- begin(total = 50L, label = "Preparing plot bootstrap inid")
      state <- tick(state, done = 1L)
      state <- tick(state, done = 2L)
      finish(state)
    }
  )

  expect_true(grepl("\\[np\\] Preparing plot bootstrap inid\\.\\.\\.", output))
  expect_true(grepl("\\[np\\] Preparing plot bootstrap inid 2/50", output))
})

test_that("plot progress default cadence matches bandwidth heartbeat", {
  plot_interval <- getFromNamespace(".np_plot_progress_interval_sec", "np")
  begin <- getFromNamespace(".np_progress_begin", "np")

  old_opts <- options(
    np.plot.progress.interval.sec = NULL,
    np.progress.interval.unknown.sec = NULL,
    np.progress.interval.known.sec = NULL
  )
  on.exit(options(old_opts), add = TRUE)

  expect_identical(plot_interval(), 2.0)
  expect_identical(begin("Bandwidth selection")$throttle_sec, 2.0)
  expect_identical(begin("Bootstrap replications", total = 10L)$throttle_sec, 0.5)
})

test_that("interactive bootstrap chunk sizes leave room for intermediate progress", {
  progress_chunk_cap <- getFromNamespace(".np_plot_progress_chunk_cap", "np")
  wild_chunk <- getFromNamespace(".np_wild_chunk_size", "np")
  enabled_chunk <- with_np_bindings(
    list(.np_progress_is_interactive = function() TRUE),
    wild_chunk(n = 500L, B = 9999L)
  )
  disabled_chunk <- with_np_bindings(
    list(.np_progress_is_interactive = function() FALSE),
    wild_chunk(n = 500L, B = 9999L)
  )

  expect_identical(progress_chunk_cap(9999L), 2500L)
  expect_identical(enabled_chunk, 16L)
  expect_identical(disabled_chunk, 9999L)
})

test_that("single-line abort preserves the final line and terminates it", {
  render <- getFromNamespace(".np_progress_render_single_line", "np")
  line <- "[np] Bootstrap replications aborted"

  output <- capture_single_line_output(
    "np",
    list(),
    render(list(render_line = line, last_width = 0L), event = "abort")
  )

  expect_identical(output, paste0("\r", line, "\n"))
})

test_that("progress ownership suppresses nested visibility for legacy renderer too", {
  begin <- getFromNamespace(".np_progress_begin", "np")
  finish <- getFromNamespace(".np_progress_end", "np")
  reset <- getFromNamespace(".np_progress_reset_registry", "np")

  reset()
  on.exit(reset(), add = TRUE)

  states <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_rstudio_console = function() TRUE,
      .np_progress_now = function() 1
    ),
    {
      outer <- begin("Outer task", total = 2)
      inner <- begin("Inner task", total = 2)
      finish(outer)
      list(outer = outer, inner = inner)
    }
  )

  expect_true(isTRUE(states$outer$visible))
  expect_false(isTRUE(states$inner$visible))
})
