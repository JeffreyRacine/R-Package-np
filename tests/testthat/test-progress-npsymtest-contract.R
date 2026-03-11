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

npsymtest_fun <- function(...) {
  getFromNamespace("npsymtest", "npRmpi")(...)
}

test_that("npsymtest emits append-only bounded bootstrap progress", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  x <- rgamma(30, shape = 2)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_counter(),
      .npRmpi_autodispatch_active = function() FALSE
    ),
    capture_progress_conditions(
      npsymtest_fun(x, method = "summation", boot.num = 9)
    )
  )

  messages <- normalize_messages(res$messages)
  bootstrap_lines <- messages[grepl("^\\[npRmpi\\] Bootstrap replications ", messages)]

  expect_s3_class(res$value, "symtest")
  expect_true(any(grepl("^\\[npRmpi\\] Computing bandwidths$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Bootstrap replications 1/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", bootstrap_lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bootstrap replications 9/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", bootstrap_lines)))
  expect_true(length(unique(bootstrap_lines)) >= 3L)
  expect_false(any(grepl("\b", messages, fixed = TRUE)))
})

test_that("npsymtest progress respects np.messages FALSE", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  x <- rgamma(30, shape = 2)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_counter(),
      .npRmpi_autodispatch_active = function() FALSE
    ),
    capture_progress_conditions(
      npsymtest_fun(x, method = "summation", boot.num = 9)
    )
  )

  expect_length(res$messages, 0)
})

test_that("npsymtest progress respects suppressMessages", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  x <- rgamma(30, shape = 2)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_counter(),
      .npRmpi_autodispatch_active = function() FALSE
    ),
    capture_progress_conditions(
      suppressMessages(npsymtest_fun(x, method = "summation", boot.num = 9))
    )
  )

  expect_length(res$messages, 0)
})
