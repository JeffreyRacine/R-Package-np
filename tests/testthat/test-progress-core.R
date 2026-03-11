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

test_that("progress begin returns disabled state when messages are off", {
  begin <- getFromNamespace(".np_progress_begin", "npRmpi")

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  state <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = function() 10
    ),
    begin("Iterating Landweber-Fridman solve")
  )

  expect_false(state$enabled)
  expect_identical(state$label, "Iterating Landweber-Fridman solve")
  expect_null(state$total)
  expect_false(state$known_total)
})

test_that("progress begin returns disabled state when not master", {
  begin <- getFromNamespace(".np_progress_begin", "npRmpi")

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  state <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() FALSE,
      .np_progress_now = function() 10
    ),
    begin("Iterating Landweber-Fridman solve")
  )

  expect_false(state$enabled)
})

test_that("known-total formatting is deterministic", {
  fmt <- getFromNamespace(".np_progress_format_known_total", "npRmpi")
  state <- list(
    enabled = TRUE,
    pkg_prefix = "[npRmpi]",
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
    "[npRmpi] Iterating Tikhonov solve 2/5 (40.0%, elapsed 4.0s, eta 6.0s): recomputing alpha"
  )
})

test_that("unknown-total formatting is deterministic", {
  fmt <- getFromNamespace(".np_progress_format_unknown_total", "npRmpi")
  state <- list(
    enabled = TRUE,
    pkg_prefix = "[npRmpi]",
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
    "[npRmpi] Iterating Landweber-Fridman solve... iteration 3, elapsed 4.0s: evaluating stopping rule"
  )
})

test_that("suppressMessages suppresses progress notes", {
  note <- getFromNamespace(".np_progress_note", "npRmpi")

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  messages <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE
    ),
    capture_messages_only(suppressMessages(note("Preparing IV derivative regression")))
  )

  expect_length(messages, 0)
})

test_that("legacy suppression restores np.messages on success and error", {
  suppress_legacy <- getFromNamespace(".np_progress_with_legacy_suppressed", "npRmpi")

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
