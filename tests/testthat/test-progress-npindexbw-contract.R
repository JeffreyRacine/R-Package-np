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
  pieces <- unlist(strsplit(paste(x, collapse = "\n"), "[\r\n]+"))
  pieces <- trimws(pieces, which = "right")
  pieces[nzchar(pieces)]
}

progress_time_counter <- function(start = 0, by = 0.6) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

skip_live_route_slice <- function() {
  skip_if_not(
    identical(Sys.getenv("NP_RMPI_PROGRESS_LIVE_ROUTE_TESTS", ""), "true"),
    "live npRmpi route slice is gated to manual session/attach/profile proof artifacts"
  )
}

test_that("npindexbw adopts the generic bandwidth selection line", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  x1 <- runif(30)
  x2 <- runif(30)
  y <- x1 + 0.5 * x2 + rnorm(30, sd = 0.1)

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  res <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_counter(),
      .npRmpi_autodispatch_active = function() FALSE
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
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 1/3\\)$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 1/3, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 2/3, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 2/3, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 3/3, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", messages)))
})

test_that("npindexbw progress respects np.messages FALSE", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(7)
  x1 <- runif(25)
  x2 <- runif(25)
  y <- x1 + 0.5 * x2 + rnorm(25, sd = 0.1)

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
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(9)
  x1 <- runif(25)
  x2 <- runif(25)
  y <- x1 + 0.5 * x2 + rnorm(25, sd = 0.1)

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
