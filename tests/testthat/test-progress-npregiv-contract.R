capture_progress_shadow_with_conditions <- function(expr, force_renderer = NULL, now = function() 0) {
  messages <- character()
  warnings <- character()

  shadow <- withCallingHandlers(
    capture_progress_shadow_trace(expr, force_renderer = force_renderer, now = now),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(
    value = shadow$value,
    trace = shadow$trace,
    final_line = shadow$final_line,
    messages = sub("\n$", "", messages),
    warnings = warnings
  )
}

skip_live_route_slice <- function() {
  skip_if_not(
    identical(Sys.getenv("NP_RMPI_PROGRESS_LIVE_ROUTE_TESTS", ""), "true"),
    "live npRmpi route slice is gated to manual session/attach/profile proof artifacts"
  )
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

shadow_lines_matching <- function(shadow, pattern) {
  lines <- vapply(shadow$trace, `[[`, character(1L), "line")
  lines[grepl(pattern, lines)]
}

shadow_signature <- function(shadow, pattern) {
  lines <- vapply(shadow$trace, `[[`, character(1L), "line")
  events <- vapply(shadow$trace, `[[`, character(1L), "event")
  keep <- grepl(pattern, lines) & events == "render"

  data.frame(
    event = events[keep],
    line = lines[keep],
    stringsAsFactors = FALSE
  )
}

test_that("Landweber npregiv single-line progress reports object labels with outer iterations", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  single_line <- capture_progress_shadow_with_conditions(
    npregiv(y = dat$y, z = dat$z, w = dat$w, method = "Landweber-Fridman", iterate.max = 4),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- vapply(single_line$trace, `[[`, character(1L), "line")

  expect_s3_class(single_line$value, "npregiv")
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[y\\|w\\], elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[[^)]*\\|z\\], elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[y-phi\\(z\\)\\|w\\], iteration 1, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[E\\[y-phi\\(z\\)\\|w\\]\\|z\\], iteration 1, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_false(any(grepl("Iterating Landweber-Fridman solve", lines, fixed = TRUE)))
  expect_false(any(grepl("%|eta ", lines)))
})

test_that("Tikhonov npregiv single-line progress restores historical object labels", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  single_line <- capture_progress_shadow_with_conditions(
    npregiv(
      y = dat$y,
      z = dat$z,
      w = dat$w,
      method = "Tikhonov",
      iterate.Tikhonov = TRUE,
      iterate.Tikhonov.num = 3
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- vapply(single_line$trace, `[[`, character(1L), "line")

  expect_s3_class(single_line$value, "npregiv")
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[y\\|w\\], elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[E\\[y\\|w\\]\\|z\\], elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(alpha, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[phi\\(z\\)\\|w\\], iteration [0-9]+, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[E\\[phi\\(z\\)\\|w\\]\\|z\\], iteration [0-9]+, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(phi\\(z\\), iteration [0-9]+, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_false(any(grepl("Iterating Tikhonov solve", lines, fixed = TRUE)))
  expect_false(any(grepl("%|eta ", lines)))
})

test_that("Landweber npregiv without residual smoothing reports the alternate historical objects", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  single_line <- capture_progress_shadow_with_conditions(
    npregiv(
      y = dat$y,
      z = dat$z,
      w = dat$w,
      method = "Landweber-Fridman",
      smooth.residuals = FALSE,
      iterate.max = 4
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- vapply(single_line$trace, `[[`, character(1L), "line")

  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[phi\\(z\\)\\|w\\], iteration 1, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] IV regression \\(E\\[E\\[y\\|w\\]-E\\[phi\\(z\\)\\|w\\]\\|z\\], iteration 1, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
})

test_that("npregivderiv single-line progress matches legacy semantics", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_with_conditions(
    npregivderiv(y = dat$y, z = dat$z, w = dat$w, iterate.max = 3),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  dat <- make_iv_data()
  single_line <- capture_progress_shadow_with_conditions(
    npregivderiv(y = dat$y, z = dat$z, w = dat$w, iterate.max = 3),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  pattern <- "^\\[npRmpi\\] Iterating Landweber-Fridman derivative solve"
  lines <- shadow_lines_matching(single_line, pattern)

  expect_s3_class(single_line$value, "npregivderiv")
  expect_equal(shadow_signature(single_line, pattern), shadow_signature(legacy, pattern))
  expect_true(any(grepl("^\\[npRmpi\\] Preparing IV derivative regression$", single_line$messages)))
  expect_true(any(grepl("^\\[npRmpi\\] Iterating Landweber-Fridman derivative solve\\.\\.\\. iteration [0-9]+, elapsed [0-9]+\\.[0-9]s", lines)))
})

test_that("npregiv proof-slice progress respects np.messages FALSE", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_with_conditions(
    npregiv(y = dat$y, z = dat$z, w = dat$w, method = "Landweber-Fridman", iterate.max = 4),
    now = progress_time_counter()
  )

  expect_length(res$messages, 0)
  expect_length(res$trace, 0)
})

test_that("npregivderiv proof-slice progress respects suppressMessages", {
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  dat <- make_iv_data()
  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_with_conditions(
    suppressMessages(npregivderiv(y = dat$y, z = dat$z, w = dat$w, iterate.max = 3)),
    now = progress_time_counter()
  )

  expect_length(res$messages, 0)
  expect_length(res$trace, 0)
})

test_that("IV source routes use the canonical single-line surface tag", {
  src <- installed_function_text("npregiv")
  expect_true(grepl(".np_progress_select_iv(", src, fixed = TRUE))

  src_deriv <- installed_function_text("npregivderiv")
  expect_true(grepl('surface = "iv_solve"', src_deriv, fixed = TRUE))
})
