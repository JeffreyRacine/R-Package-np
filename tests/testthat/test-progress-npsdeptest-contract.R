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

capture_progress_conditions <- function(expr) {
  messages <- character()
  warnings <- character()

  value <- withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, messages = messages, warnings = warnings)
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

npsdeptest_fun <- function(...) {
  getFromNamespace("npsdeptest", "np")(...)
}

test_that("npsdeptest emits append-only bounded lag and bootstrap progress", {
  set.seed(42)
  y <- arima.sim(n = 50, list(ar = 0.5))

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npsdeptest_fun(y, lag.num = 2, method = "summation", boot.num = 9)
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "sdeptest")
  expect_true(any(grepl("^\\[np\\] Computing bandwidths$", messages)))
  expect_true(any(grepl("^\\[np\\] Constructing metric entropy by lag 1/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): lag 1$", messages)))
  expect_true(any(grepl("^\\[np\\] Constructing metric entropy by lag 2/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): lag 2$", messages)))
  expect_true(any(grepl("^\\[np\\] Bootstrap replications 1/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[np\\] Bootstrap replications 9/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_false(any(grepl("\b", messages, fixed = TRUE)))
})

test_that("npsdeptest progress respects np.messages FALSE", {
  set.seed(42)
  y <- arima.sim(n = 50, list(ar = 0.5))

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npsdeptest_fun(y, lag.num = 2, method = "summation", boot.num = 9)
    )
  )

  expect_length(res$messages, 0)
})

test_that("npsdeptest progress respects suppressMessages", {
  set.seed(42)
  y <- arima.sim(n = 50, list(ar = 0.5))

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      suppressMessages(npsdeptest_fun(y, lag.num = 2, method = "summation", boot.num = 9))
    )
  )

  expect_length(res$messages, 0)
})
