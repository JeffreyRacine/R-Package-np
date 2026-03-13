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
  keep <- grepl("^\\[np\\] Constructing metric entropy by lag", lines) |
    grepl("^\\[np\\] Bootstrap replications", lines)

  data.frame(
    event = events[keep],
    line = lines[keep],
    stringsAsFactors = FALSE
  )
}

shadow_lines <- function(shadow) {
  shadow_surface_signature(shadow)$line
}

npsdeptest_fun <- function(...) {
  getFromNamespace("npsdeptest", "np")(...)
}

test_that("npsdeptest single-line lag and bootstrap progress match legacy semantics", {
  set.seed(42)
  y <- arima.sim(n = 50, list(ar = 0.5))

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0
  )
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
  expect_true(any(grepl("^\\[np\\] Constructing metric entropy by lag 1/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): lag 1$", lines)))
  expect_true(any(grepl("^\\[np\\] Constructing metric entropy by lag 2/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): lag 2$", lines)))
  expect_true(any(grepl("^\\[np\\] Bootstrap replications [0-9]+/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bootstrap replications 9/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
})

test_that("npsdeptest progress respects np.messages FALSE", {
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
