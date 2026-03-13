progress_time_counter <- function(start = 0, by = 0.6) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

shadow_surface_signature <- function(shadow) {
  lines <- vapply(shadow$trace, `[[`, character(1L), "line")
  events <- vapply(shadow$trace, `[[`, character(1L), "event")
  keep <- grepl("^\\[npRmpi\\] Constructing metric entropy by lag|^\\[npRmpi\\] Bootstrap replications", lines)

  data.frame(
    event = events[keep],
    line = lines[keep],
    stringsAsFactors = FALSE
  )
}

shadow_lines <- function(shadow) {
  shadow_surface_signature(shadow)$line
}

skip_live_route_slice <- function() {
  skip_if_not(
    identical(Sys.getenv("NP_RMPI_PROGRESS_LIVE_ROUTE_TESTS", ""), "true"),
    "live npRmpi route slice is gated to manual session/attach/profile proof artifacts"
  )
}

npsdeptest_fun <- function(...) {
  getFromNamespace("npsdeptest", "npRmpi")(...)
}

test_that("npsdeptest single-line lag and bootstrap progress match legacy semantics", {
  skip_on_cran()
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  y <- arima.sim(n = 50, list(ar = 0.5))

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    npsdeptest_fun(y, lag.num = 2, method = "summation", boot.num = 9),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  set.seed(42)
  single_line <- capture_progress_shadow_trace(
    npsdeptest_fun(y, lag.num = 2, method = "summation", boot.num = 9),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(single_line)

  expect_s3_class(single_line$value, "sdeptest")
  expect_equal(shadow_surface_signature(single_line), shadow_surface_signature(legacy))
  expect_true(any(grepl("^\\[npRmpi\\] Constructing metric entropy by lag 1/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): lag 1$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Constructing metric entropy by lag 2/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): lag 2$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bootstrap replications [0-9]+/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[npRmpi\\] Bootstrap replications 9/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
})

test_that("npsdeptest progress respects np.messages FALSE", {
  skip_on_cran()
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  y <- arima.sim(n = 50, list(ar = 0.5))

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_trace(
    npsdeptest_fun(y, lag.num = 2, method = "summation", boot.num = 9),
    now = progress_time_counter()
  )

  expect_length(res$trace, 0)
})

test_that("npsdeptest progress respects suppressMessages", {
  skip_on_cran()
  skip_live_route_slice()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  y <- arima.sim(n = 50, list(ar = 0.5))

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_trace(
    suppressMessages(npsdeptest_fun(y, lag.num = 2, method = "summation", boot.num = 9)),
    now = progress_time_counter()
  )

  expect_length(res$trace, 0)
})

test_that("npsdeptest source routes use canonical lag and bootstrap surface tags", {
  src <- installed_function_text("npsdeptest")
  expect_true(grepl('surface = "lag"', src, fixed = TRUE))
  expect_true(grepl('surface = "bootstrap"', src, fixed = TRUE))
})
