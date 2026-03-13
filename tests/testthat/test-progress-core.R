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
  messages <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_values(c(0, 1.1, 2.2, 3.3, 4.4))
    ),
    capture_messages_only({
      value <- select_bw("Selecting regression bandwidth", {
        seen <<- getOption("np.messages")
        set_total(3)
        step_bw(1, 3)
        step_bw(3, 3)
        7
      })
      expect_identical(value, 7)
    })
  )
  messages <- sub("\n$", "", messages, useBytes = TRUE)
  messages <- messages[nzchar(messages)]

  expect_false(isTRUE(seen))
  expect_true(any(grepl("^\\[np\\] Selecting regression bandwidth\\.\\.\\.$", messages)))
  expect_true(any(grepl("^\\[np\\] Selecting regression bandwidth multistart 1/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[np\\] Selecting regression bandwidth multistart 3/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(isTRUE(getOption("np.messages")))
})
