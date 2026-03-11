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

npsigtest_fun <- function(...) {
  getFromNamespace("npsigtest", "np")(...)
}

make_sigtest_fixture <- function(seed = 42, n = 30) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 + rnorm(n, sd = 0.1)
  bw <- npregbw(y ~ x1 + x2, bws = c(0.2, 0.4), bandwidth.compute = FALSE)
  list(x1 = x1, x2 = x2, y = y, bw = bw)
}

test_that("npsigtest joint path emits append-only bounded bootstrap progress", {
  fixture <- make_sigtest_fixture()

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npsigtest_fun(bws = fixture$bw, boot.num = 9, joint = TRUE, index = 1)
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "sigtest")
  expect_true(any(grepl("^\\[np\\] Testing joint significance$", messages)))
  expect_true(any(grepl("^\\[np\\] Bootstrap replications 1/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[np\\] Bootstrap replications 9/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_false(any(grepl("\b", messages, fixed = TRUE)))
})

test_that("npsigtest individual path emits per-variable append-only progress", {
  fixture <- make_sigtest_fixture(seed = 99)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npsigtest_fun(bws = fixture$bw, boot.num = 9, joint = FALSE, index = c(1, 2))
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "sigtest")
  expect_true(any(grepl("^\\[np\\] Testing variable 1 of \\(1,2\\)$", messages)))
  expect_true(any(grepl("^\\[np\\] Testing variable 2 of \\(1,2\\)$", messages)))
  expect_true(sum(grepl("^\\[np\\] Bootstrap replications 1/9 ", messages)) >= 2L)
  expect_true(sum(grepl("^\\[np\\] Bootstrap replications 9/9 ", messages)) >= 2L)
  expect_false(any(grepl("\b", messages, fixed = TRUE)))
})

test_that("npsigtest progress respects np.messages FALSE", {
  fixture <- make_sigtest_fixture(seed = 17)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      npsigtest_fun(bws = fixture$bw, boot.num = 9, joint = TRUE, index = 1)
    )
  )

  expect_length(res$messages, 0)
})

test_that("npsigtest progress respects suppressMessages", {
  fixture <- make_sigtest_fixture(seed = 23)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_conditions(
      suppressMessages(npsigtest_fun(bws = fixture$bw, boot.num = 9, joint = TRUE, index = 1))
    )
  )

  expect_length(res$messages, 0)
})
