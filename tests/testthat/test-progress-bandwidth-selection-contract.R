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

normalize_messages <- function(x) {
  sub("\n$", "", x)
}

progress_time_counter <- function(start = 0, by = 1.1) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

test_that("npudensbw emits bounded multistart bandwidth progress on master", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  x <- rnorm(35)

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old_opts), add = TRUE)

  res <- NULL
  messages <- capture.output(
    res <- with_nprmpi_bindings(
      list(
        .np_progress_is_interactive = function() TRUE,
        .np_progress_is_master = function() TRUE,
        .np_progress_now = progress_time_counter(),
        .npRmpi_autodispatch_active = function() FALSE
      ),
      npudensbw(
        dat = data.frame(x = x),
        bwmethod = "cv.ml",
        nmulti = 3
      )
    ),
    type = "message"
  )
  messages <- normalize_messages(messages)

  expect_s3_class(res, "bandwidth")
  expect_true(any(grepl("^\\[npRmpi\\] Selecting density bandwidth 1/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): multistart 1 of 3$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Selecting density bandwidth 3/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): multistart 3 of 3$", messages)))
})

test_that("npregbw emits bounded multistart bandwidth progress on master", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(7)
  x <- runif(30)
  y <- sin(2 * pi * x) + rnorm(30, sd = 0.1)

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old_opts), add = TRUE)

  res <- NULL
  messages <- capture.output(
    res <- with_nprmpi_bindings(
      list(
        .np_progress_is_interactive = function() TRUE,
        .np_progress_is_master = function() TRUE,
        .np_progress_now = progress_time_counter(),
        .npRmpi_autodispatch_active = function() FALSE
      ),
      npregbw(
        xdat = data.frame(x = x),
        ydat = y,
        regtype = "lc",
        bwmethod = "cv.aic",
        nmulti = 3
      )
    ),
    type = "message"
  )
  messages <- normalize_messages(messages)

  expect_s3_class(res, "rbandwidth")
  expect_true(any(grepl("^\\[npRmpi\\] Selecting regression bandwidth 1/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): multistart 1 of 3$", messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Selecting regression bandwidth 3/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): multistart 3 of 3$", messages)))
})
