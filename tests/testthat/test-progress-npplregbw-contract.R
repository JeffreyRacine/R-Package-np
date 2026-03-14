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

progress_time_counter <- function(start = 0, by = 1.7) {
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

shadow_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

test_that("npplregbw uses coordinated generic bandwidth selection progress", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  n <- 35
  x <- runif(n, -1, 1)
  z <- rnorm(n)
  y <- x^2 + z + rnorm(n, sd = 0.25 * sd(x))

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_counter(),
      .np_progress_output_width = function() 500L,
      .npRmpi_autodispatch_active = function() FALSE
    ),
    capture_progress_shadow_trace(
      npplregbw(
        xdat = data.frame(x = x),
        zdat = data.frame(z = z),
        ydat = y,
        nmulti = 2
      ),
      force_renderer = "single_line"
    )
  )

  lines <- shadow_lines(actual)

  expect_s3_class(actual$value, "plbandwidth")
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(y~z, multistart 1/2\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(y~z, multistart 2/2, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(x1~z, multistart 1/2, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(x1~z, multistart 2/2, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})

test_that("npplreg formula entry inherits coordinated generic bandwidth progress", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  n <- 35
  dat <- data.frame(x = runif(n, -1, 1), z = rnorm(n))
  dat$y <- dat$x^2 + dat$z + rnorm(n, sd = 0.25 * sd(dat$x))

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_counter(),
      .np_progress_output_width = function() 500L,
      .npRmpi_autodispatch_active = function() FALSE
    ),
    capture_progress_shadow_trace(
      npplreg(
        y ~ x | z,
        data = dat,
        nmulti = 2
      ),
      force_renderer = "single_line"
    )
  )

  lines <- shadow_lines(actual)

  expect_s3_class(actual$value, "plregression")
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(y~z, multistart 1/2\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(x1~z, multistart 1/2, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(x1~z, multistart 2/2, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})
