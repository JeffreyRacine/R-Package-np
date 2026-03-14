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

progress_time_counter <- function(start = 0, by = 2.1) {
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

test_that("npscoefbw adopts the generic bandwidth selection line", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(3240)
  n <- 28
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + x * (1 + z) + rnorm(n, sd = 0.1)

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
      .np_progress_output_width = function() 500L
    ),
    capture_progress_shadow_trace(
      npscoefbw(
        xdat = data.frame(x = x),
        zdat = data.frame(z = z),
        ydat = y,
        regtype = "lc",
        nmulti = 2,
        optim.maxit = 3,
        cv.iterate = FALSE
      ),
      force_renderer = "single_line"
    )
  )

  lines <- shadow_lines(actual)

  expect_s3_class(actual$value, "scbandwidth")
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 1/2\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 1/2, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 2/2, elapsed [0-9]+\\.[0-9]s, 50\\.0%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 2/2, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bandwidth selection \\(multistart 2/2, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})

test_that("npscoefbw progress respects np.messages FALSE", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(3240)
  n <- 24
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + x * (1 + z) + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  silent <- with_nprmpi_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_master = function() TRUE,
      .np_progress_now = progress_time_counter()
    ),
    capture_progress_shadow_trace(
      npscoefbw(
        xdat = data.frame(x = x),
        zdat = data.frame(z = z),
        ydat = y,
        regtype = "lc",
        nmulti = 1,
        optim.maxit = 2,
        cv.iterate = FALSE
      )
    )
  )

  expect_length(silent$trace, 0L)
})

test_that("npscoefbw cv.iterate path retains backfitting progress hooks", {
  src <- installed_function_text("npscoefbw.scbandwidth")

  expect_true(grepl("Backfitting smooth coefficient bandwidth", src, fixed = TRUE))
  expect_true(grepl("Optimizing partial residual bandwidth", src, fixed = TRUE))
  expect_true(grepl("\\.np_progress_begin\\(\"Backfitting smooth coefficient bandwidth\"", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Optimizing partial residual bandwidth\"", src))
})
