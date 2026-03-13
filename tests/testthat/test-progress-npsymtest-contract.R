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

npsymtest_fun <- function(...) {
  getFromNamespace("npsymtest", "np")(...)
}

test_that("npsymtest single-line bootstrap progress matches legacy semantics", {
  set.seed(42)
  x <- rgamma(30, shape = 2)

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    npsymtest_fun(x, method = "summation", boot.num = 9),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  set.seed(42)
  single_line <- capture_progress_shadow_trace(
    npsymtest_fun(x, method = "summation", boot.num = 9),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(single_line)
  bootstrap_lines <- lines[grepl("^\\[np\\] Bootstrap replications ", lines)]

  expect_s3_class(single_line$value, "symtest")
  expect_equal(shadow_bootstrap_signature(single_line), shadow_bootstrap_signature(legacy))
  expect_true(any(grepl("^\\[np\\] Bootstrap replications [0-9]+/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", bootstrap_lines)))
  expect_true(any(grepl("^\\[np\\] Bootstrap replications 9/9 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", bootstrap_lines)))
  expect_true(length(unique(bootstrap_lines)) >= 3L)
})

test_that("npsymtest progress respects np.messages FALSE", {
  set.seed(42)
  x <- rgamma(30, shape = 2)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_trace(
    npsymtest_fun(x, method = "summation", boot.num = 9),
    now = progress_time_counter()
  )

  expect_length(res$trace, 0)
})

test_that("npsymtest progress respects suppressMessages", {
  set.seed(42)
  x <- rgamma(30, shape = 2)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  res <- capture_progress_shadow_trace(
    suppressMessages(npsymtest_fun(x, method = "summation", boot.num = 9)),
    now = progress_time_counter()
  )

  expect_length(res$trace, 0)
})
