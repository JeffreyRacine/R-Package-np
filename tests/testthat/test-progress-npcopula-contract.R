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
  value <- NULL

  message.output <- utils::capture.output(
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
    ),
    type = "message"
  )

  list(value = value, messages = c(messages, message.output), warnings = warnings)
}

normalize_messages <- function(x) {
  sub("\n$", "", x)
}

make_npcopula_fixture <- function(seed = 42, n = 30) {
  set.seed(seed)
  data.frame(
    x = runif(n),
    y = 0.7 * runif(n) + 0.3 * rnorm(n, sd = 0.1)
  )
}

test_that("npcopula sample-realization path emits compact staged progress", {
  skip_if_not(exists("npcopula", mode = "function"), "np package not attached")
  mydat <- make_npcopula_fixture()
  bw <- npudistbw(dat = mydat, bws = c(0.2, 0.2), bandwidth.compute = FALSE)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(.np_progress_is_interactive = function() TRUE),
    capture_progress_conditions(
      npcopula(bws = bw, data = mydat)
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "npcopula")
  expect_true(is.data.frame(as.data.frame(res$value)))
  expect_true(any(grepl("\\[np\\] Copula (distribution|dist) sample", messages)))
  expect_false(any(grepl("^\\[np\\] Computing the marginal of", messages)))
})

test_that("npcopula u-grid density path emits compact staged progress", {
  skip_if_not(exists("npcopula", mode = "function"), "np package not attached")
  mydat <- make_npcopula_fixture(seed = 99)
  bw <- npudensbw(dat = mydat, bws = c(0.2, 0.2), bandwidth.compute = FALSE)
  u <- as.matrix(data.frame(x = c(0.25, 0.75), y = c(0.25, 0.75)))

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(.np_progress_is_interactive = function() TRUE),
    capture_progress_conditions(
      npcopula(
        bws = bw,
        data = mydat,
        u = u,
        n.quasi.inv = 40
      )
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "npcopula")
  expect_true(is.data.frame(as.data.frame(res$value)))
  expect_true(any(grepl("\\[np\\] Copula (density|dens) grid", messages)))
  expect_false(any(grepl("^\\[np\\] Computing the quasi-inverse", messages)))
  expect_false(any(grepl("^\\[np\\] Computing the marginal of", messages)))
})

test_that("npcopula progress respects np.messages FALSE", {
  skip_if_not(exists("npcopula", mode = "function"), "np package not attached")
  mydat <- make_npcopula_fixture(seed = 17)
  bw <- npudistbw(dat = mydat, bws = c(0.2, 0.2), bandwidth.compute = FALSE)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(.np_progress_is_interactive = function() TRUE),
    capture_progress_conditions(
      npcopula(bws = bw, data = mydat)
    )
  )

  expect_length(res$messages, 0)
})

test_that("npcopula progress respects suppressMessages", {
  skip_if_not(exists("npcopula", mode = "function"), "np package not attached")
  mydat <- make_npcopula_fixture(seed = 23)
  bw <- npudistbw(dat = mydat, bws = c(0.2, 0.2), bandwidth.compute = FALSE)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_np_bindings(
    list(.np_progress_is_interactive = function() TRUE),
    capture_progress_conditions(
      suppressMessages(npcopula(bws = bw, data = mydat))
    )
  )

  expect_length(res$messages, 0)
})
