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

progress_time_counter <- function(start = 0, by = 2.1) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

test_that("npscoefbw emits append-only multistart and objective progress", {
  set.seed(3240)
  n <- 65
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + x * (1 + z) + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npscoefbw(
        xdat = data.frame(x = x),
        zdat = data.frame(z = z),
        ydat = y,
        regtype = "lc",
        nmulti = 2,
        optim.maxit = 10,
        cv.iterate = FALSE
      )
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "scbandwidth")
  expect_true(any(grepl("^\\[np\\] Selecting smooth coefficient bandwidth multistart 1/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[np\\] Selecting smooth coefficient bandwidth multistart 2/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[np\\] Optimizing smooth coefficient bandwidth\\.\\.\\. iteration [0-9]+, elapsed [0-9]+\\.[0-9]s: multistart 1$", messages)))
  expect_false(any(grepl(intToUtf8(8L), messages, fixed = TRUE)))
})

test_that("npscoefbw progress respects np.messages FALSE", {
  set.seed(3240)
  n <- 50
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + x * (1 + z) + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npscoefbw(
        xdat = data.frame(x = x),
        zdat = data.frame(z = z),
        ydat = y,
        regtype = "lc",
        nmulti = 1,
        optim.maxit = 5,
        cv.iterate = FALSE
      )
    )
  )

  expect_length(res$messages, 0)
})

test_that("npscoefbw cv.iterate path emits append-only backfitting progress", {
  set.seed(3240)
  n <- 65
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + x * (1 + z) + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npscoefbw(
        xdat = data.frame(x = x),
        zdat = data.frame(z = z),
        ydat = y,
        regtype = "lc",
        nmulti = 1,
        optim.maxit = 5,
        cv.iterate = TRUE,
        cv.num.iterations = 2,
        backfit.iterate = FALSE
      )
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "scbandwidth")
  expect_true(any(grepl("^\\[np\\] Backfitting smooth coefficient bandwidth 1/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): iteration 1 of 2$", messages)))
  expect_true(any(grepl("^\\[np\\] Optimizing partial residual bandwidth\\.\\.\\. iteration [0-9]+, elapsed [0-9]+\\.[0-9]s: backfitting iteration [0-9]+ of 2, partial residual [0-9]+ of 2, fval ", messages)))
  expect_false(any(grepl(intToUtf8(8L), messages, fixed = TRUE)))
})
