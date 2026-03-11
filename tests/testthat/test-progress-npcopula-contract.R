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

make_npcopula_fixture <- function(seed = 42, n = 30) {
  set.seed(seed)
  data.frame(
    x = runif(n),
    y = 0.7 * runif(n) + 0.3 * rnorm(n, sd = 0.1)
  )
}

test_that("npcopula sample-realization path emits append-only progress notes", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  mydat <- make_npcopula_fixture()
  bw <- npudistbw(dat = mydat, bws = c(0.2, 0.2), bandwidth.compute = FALSE)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .npRmpi_autodispatch_active = function() FALSE
    ),
    capture_progress_conditions(
      npcopula(bws = bw, data = mydat)
    )
  )

  messages <- normalize_messages(res$messages)

  expect_s3_class(res$value, "data.frame")
  expect_true(any(grepl("^\\[npRmpi\\] Computing the copula for the sample realizations$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Computing the marginal of x for the sample realizations$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Computing the marginal of y for the sample realizations$", messages)))
  expect_false(any(grepl("\b", messages, fixed = TRUE)))
})

test_that("npcopula u-grid density path emits append-only progress notes", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  mydat <- make_npcopula_fixture(seed = 99)
  bw <- npudensbw(dat = mydat, bws = c(0.2, 0.2), bandwidth.compute = FALSE)
  u <- as.matrix(data.frame(x = c(0.25, 0.75), y = c(0.25, 0.75)))

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .npRmpi_autodispatch_active = function() FALSE
    ),
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

  expect_s3_class(res$value, "data.frame")
  expect_true(any(grepl("^\\[npRmpi\\] Computing the quasi-inverse for the marginal of x$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Computing the quasi-inverse for the marginal of y$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Expanding the u matrix$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Computing the copula density for the expanded grid$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Computing the marginal of x for the expanded grid$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Computing the marginal of y for the expanded grid$", messages)))
  expect_false(any(grepl("\b", messages, fixed = TRUE)))
})

test_that("npcopula progress respects np.messages FALSE", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  mydat <- make_npcopula_fixture(seed = 17)
  bw <- npudistbw(dat = mydat, bws = c(0.2, 0.2), bandwidth.compute = FALSE)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .npRmpi_autodispatch_active = function() FALSE
    ),
    capture_progress_conditions(
      npcopula(bws = bw, data = mydat)
    )
  )

  expect_length(res$messages, 0)
})

test_that("npcopula progress respects suppressMessages", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  mydat <- make_npcopula_fixture(seed = 23)
  bw <- npudistbw(dat = mydat, bws = c(0.2, 0.2), bandwidth.compute = FALSE)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .npRmpi_autodispatch_active = function() FALSE
    ),
    capture_progress_conditions(
      suppressMessages(npcopula(bws = bw, data = mydat))
    )
  )

  expect_length(res$messages, 0)
})
