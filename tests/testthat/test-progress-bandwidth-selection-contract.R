progress_time_counter <- function(start = 0, by = 1.1) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

shadow_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

shadow_signature <- function(shadow) {
  data.frame(
    event = vapply(shadow$trace, `[[`, character(1L), "event"),
    line = shadow_lines(shadow),
    stringsAsFactors = FALSE
  )
}

test_that("npudensbw single-line bandwidth progress matches legacy semantics", {
  set.seed(42)
  x <- rnorm(35)

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    npudensbw(
      dat = data.frame(x = x),
      bwmethod = "cv.ml",
      nmulti = 3
    ),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  set.seed(42)
  single_line <- capture_progress_shadow_trace(
    npudensbw(
      dat = data.frame(x = x),
      bwmethod = "cv.ml",
      nmulti = 3
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(single_line)

  expect_s3_class(single_line$value, "bandwidth")
  expect_equal(shadow_signature(single_line), shadow_signature(legacy))
  expect_true(any(grepl("^\\[np\\] Selecting density bandwidth multistart 1/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Selecting density bandwidth multistart 3/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
})

test_that("npregbw single-line bandwidth progress matches legacy semantics", {
  set.seed(7)
  x <- runif(30)
  y <- sin(2 * pi * x) + rnorm(30, sd = 0.1)

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    npregbw(
      xdat = data.frame(x = x),
      ydat = y,
      regtype = "lc",
      bwmethod = "cv.aic",
      nmulti = 3
    ),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  set.seed(7)
  single_line <- capture_progress_shadow_trace(
    npregbw(
      xdat = data.frame(x = x),
      ydat = y,
      regtype = "lc",
      bwmethod = "cv.aic",
      nmulti = 3
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(single_line)

  expect_s3_class(single_line$value, "rbandwidth")
  expect_equal(shadow_signature(single_line), shadow_signature(legacy))
  expect_true(any(grepl("^\\[np\\] Selecting regression bandwidth multistart 1/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Selecting regression bandwidth multistart 3/3 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", lines)))
})
