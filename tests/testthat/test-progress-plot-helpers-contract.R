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

test_that("plot helper progress emits append-only bounded messages", {
  begin <- getFromNamespace(".np_plot_progress_begin", "np")
  tick <- getFromNamespace(".np_plot_progress_tick", "np")
  finish <- getFromNamespace(".np_plot_progress_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  messages <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_messages_only({
      state <- begin(total = 9, label = "Plot bootstrap wild")
      state <- tick(state, done = 1)
      state <- tick(state, done = 9)
      finish(state)
    })
  )

  messages <- normalize_messages(messages)

  expect_true(any(grepl("^\\[np\\] Plot bootstrap wild 1/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[np\\] Plot bootstrap wild 9/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_equal(sum(grepl("^\\[np\\] Plot bootstrap wild 9/9 ", messages)), 1L)
  expect_false(any(grepl(intToUtf8(8L), messages, fixed = TRUE)))
})

test_that("plot helper activity emits a single append-only note", {
  begin <- getFromNamespace(".np_plot_activity_begin", "np")
  finish <- getFromNamespace(".np_plot_activity_end", "np")

  old_opts <- options(np.messages = TRUE, np.plot.progress = TRUE)
  on.exit(options(old_opts), add = TRUE)

  messages <- with_np_bindings(
    list(.np_progress_is_interactive = function() TRUE),
    capture_messages_only({
      activity <- begin("Constructing bootstrap bands...")
      finish(activity)
    })
  )

  messages <- normalize_messages(messages)
  expect_identical(messages, "[np] Constructing bootstrap bands...")
})

test_that("plot helper progress is silent by default in noninteractive mode", {
  begin <- getFromNamespace(".np_plot_progress_begin", "np")
  tick <- getFromNamespace(".np_plot_progress_tick", "np")
  finish <- getFromNamespace(".np_plot_progress_end", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.plot.progress = TRUE,
    np.plot.progress.interval.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  messages <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() FALSE,
      .np_progress_now = progress_time_counter()
    ),
    capture_messages_only({
      state <- begin(total = 5, label = "Plot bootstrap wild")
      state <- tick(state, done = 5)
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
    np.plot.progress.interval.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  messages <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_messages_only(
      suppressMessages({
        state <- begin(total = 5, label = "Plot bootstrap wild")
        state <- tick(state, done = 5)
        finish(state)
      })
    )
  )

  expect_length(messages, 0)
})
