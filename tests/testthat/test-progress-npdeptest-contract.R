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

npdeptest_fun <- function(...) {
  getFromNamespace("npdeptest", "np")(...)
}

test_that("npdeptest emits append-only bounded bootstrap progress", {
  set.seed(42)
  n <- 30
  x <- rnorm(n)
  y <- x + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npdeptest_fun(x, y, method = "summation", boot.num = 9)
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "deptest")
  expect_true(any(grepl("^\\[np\\] Computing bandwidths$", messages)))
  expect_true(any(grepl("^\\[np\\] Constructing metric entropy$", messages)))
  expect_true(any(grepl("^\\[np\\] Bootstrap replications [0-9]+/[0-9]+ \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_false(any(grepl("\b", messages, fixed = TRUE)))
})

test_that("npdeptest progress respects np.messages FALSE", {
  set.seed(42)
  n <- 30
  x <- rnorm(n)
  y <- x + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npdeptest_fun(x, y, method = "summation", boot.num = 9)
    )
  )

  expect_length(res$messages, 0)
})

test_that("npdeptest progress respects suppressMessages", {
  set.seed(42)
  n <- 30
  x <- rnorm(n)
  y <- x + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      suppressMessages(npdeptest_fun(x, y, method = "summation", boot.num = 9))
    )
  )

  expect_length(res$messages, 0)
})
