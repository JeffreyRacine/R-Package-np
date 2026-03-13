progress_time_counter <- function(start = 0, by = 0.6) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

shadow_bootstrap_signature <- function(shadow) {
  lines <- vapply(shadow$trace, `[[`, character(1L), "line")
  events <- vapply(shadow$trace, `[[`, character(1L), "event")
  keep <- grepl("^\\[np\\] Bootstrap replications", lines)

  data.frame(
    event = events[keep],
    line = lines[keep],
    stringsAsFactors = FALSE
  )
}

shadow_lines <- function(shadow) {
  shadow_bootstrap_signature(shadow)$line
}

npdeptest_fun <- function(...) {
  getFromNamespace("npdeptest", "np")(...)
}

test_that("npdeptest single-line bootstrap progress matches legacy semantics", {
  set.seed(42)
  n <- 30
  x <- rnorm(n)
  y <- x + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    npdeptest_fun(x, y, method = "summation", boot.num = 9),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  set.seed(42)
  single_line <- capture_progress_shadow_trace(
    npdeptest_fun(x, y, method = "summation", boot.num = 9),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(single_line)

  expect_s3_class(single_line$value, "deptest")
  expect_equal(shadow_bootstrap_signature(single_line), shadow_bootstrap_signature(legacy))
  expect_true(any(grepl("^\\[np\\] Bootstrap replications [0-9]+/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bootstrap replications 9/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
})

test_that("npdeptest progress respects np.messages FALSE", {
  set.seed(42)
  n <- 30
  x <- rnorm(n)
  y <- x + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_trace(
    npdeptest_fun(x, y, method = "summation", boot.num = 9),
    now = progress_time_counter()
  )

  expect_length(res$trace, 0)
})

test_that("npdeptest progress respects suppressMessages", {
  set.seed(42)
  n <- 30
  x <- rnorm(n)
  y <- x + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_trace(
    suppressMessages(npdeptest_fun(x, y, method = "summation", boot.num = 9)),
    now = progress_time_counter()
  )

  expect_length(res$trace, 0)
})
