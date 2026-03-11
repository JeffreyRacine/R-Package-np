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

test_that("npindexbw emits append-only bounded multistart progress", {
  set.seed(42)
  x1 <- runif(30)
  x2 <- runif(30)
  y <- x1 + 0.5 * x2 + rnorm(30, sd = 0.1)

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npindexbw(
        xdat = data.frame(x1 = x1, x2 = x2),
        ydat = y,
        method = "ichimura",
        nmulti = 3,
        optim.maxit = 100
      )
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "sibandwidth")
  expect_true(any(grepl("^\\[np\\] Multistart optimization 1/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): multistart 1$", messages)))
  expect_true(any(grepl("^\\[np\\] Multistart optimization 3/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): multistart 3$", messages)))
  expect_false(any(grepl("\b", messages, fixed = TRUE)))
})

test_that("npindexbw progress respects np.messages FALSE", {
  set.seed(7)
  x1 <- runif(25)
  x2 <- runif(25)
  y <- x1 + 0.5 * x2 + rnorm(25, sd = 0.1)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npindexbw(
        xdat = data.frame(x1 = x1, x2 = x2),
        ydat = y,
        method = "ichimura",
        nmulti = 3,
        optim.maxit = 100
      )
    )
  )

  expect_length(res$messages, 0)
})

test_that("npindexbw progress respects suppressMessages", {
  set.seed(9)
  x1 <- runif(25)
  x2 <- runif(25)
  y <- x1 + 0.5 * x2 + rnorm(25, sd = 0.1)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      suppressMessages(
        npindexbw(
          xdat = data.frame(x1 = x1, x2 = x2),
          ydat = y,
          method = "ichimura",
          nmulti = 3,
          optim.maxit = 100
        )
      )
    )
  )

  expect_length(res$messages, 0)
})
