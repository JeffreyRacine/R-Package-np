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
  trace <- shadow$trace[vapply(shadow$trace, `[[`, character(1L), "event") == "render"]
  data.frame(
    event = vapply(trace, `[[`, character(1L), "event"),
    line = vapply(trace, `[[`, character(1L), "line"),
    stringsAsFactors = FALSE
  )
}

shadow_stage_signature <- function(shadow) {
  sig <- shadow_signature(shadow)
  sig$line <- sub(", elapsed .*", "", sig$line)
  sig
}

test_that("npudensbw bandwidth progress uses the enhanced multistart handoff", {
  set.seed(42)
  x <- rnorm(35)

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0
  )
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
  expect_equal(shadow_stage_signature(single_line), shadow_stage_signature(legacy))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 1/3\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/3, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 3/3, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 3/3, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})

test_that("npudens indirect entry inherits the enhanced density bandwidth progress", {
  set.seed(101)
  x <- rnorm(35)

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npudens(
      tdat = data.frame(x = x),
      bwmethod = "cv.ml",
      nmulti = 3
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(actual)

  expect_s3_class(actual$value, "npdensity")
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 1/3\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/3, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 3/3, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})

test_that("npregbw adopts the generic bandwidth selection line", {
  set.seed(7)
  x <- runif(30)
  y <- sin(2 * pi * x) + rnorm(30, sd = 0.1)

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0
  )
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
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 1/3\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/3, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 3/3, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})
