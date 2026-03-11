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

progress_time_counter <- function(start = 0, by = 2.1) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

make_iv_data <- function(n = 40) {
  set.seed(42)
  z <- runif(n)
  w <- z + rnorm(n, sd = 0.1)
  y <- z^2 + rnorm(n, sd = 0.1)
  list(y = y, z = z, w = w)
}

npregiv_fun <- function(...) {
  getFromNamespace("npregiv", "npRmpi")(...)
}

npregivderiv_fun <- function(...) {
  getFromNamespace("npregivderiv", "npRmpi")(...)
}

test_that("Landweber npregiv emits append-only progress messages", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
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
      npregiv_fun(y = dat$y, z = dat$z, w = dat$w, method = "Landweber-Fridman", iterate.max = 4)
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "npregiv")
  expect_true(any(grepl("^\\[npRmpi\\] Preparing Landweber-Fridman IV regression$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Iterating Landweber-Fridman solve\\.\\.\\. iteration [0-9]+, elapsed [0-9]+\\.[0-9]s", messages)))
  expect_false(any(grepl("\b", messages, fixed = TRUE)))
})

test_that("Tikhonov npregiv emits append-only progress messages", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
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
      npregiv_fun(
        y = dat$y,
        z = dat$z,
        w = dat$w,
        method = "Tikhonov",
        iterate.Tikhonov = TRUE,
        iterate.Tikhonov.num = 3
      )
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "npregiv")
  expect_true(any(grepl("^\\[npRmpi\\] Preparing Tikhonov IV regression$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Iterating Tikhonov solve [0-9]+/[0-9]+ \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)", messages)))
  expect_false(any(grepl("\b", messages, fixed = TRUE)))
})

test_that("npregivderiv emits append-only progress messages", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
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
      npregivderiv_fun(y = dat$y, z = dat$z, w = dat$w, iterate.max = 3)
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "npregivderiv")
  expect_true(any(grepl("^\\[npRmpi\\] Preparing IV derivative regression$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Iterating Landweber-Fridman derivative solve\\.\\.\\. iteration [0-9]+, elapsed [0-9]+\\.[0-9]s", messages)))
  expect_false(any(grepl("\b", messages, fixed = TRUE)))
})

test_that("npregiv proof-slice progress respects np.messages FALSE", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
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
      npregiv_fun(y = dat$y, z = dat$z, w = dat$w, method = "Landweber-Fridman", iterate.max = 4)
    )
  )

  expect_length(res$messages, 0)
})

test_that("npregivderiv proof-slice progress respects suppressMessages", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
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
      suppressMessages(npregivderiv_fun(y = dat$y, z = dat$z, w = dat$w, iterate.max = 3))
    )
  )

  expect_length(res$messages, 0)
})
