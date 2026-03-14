progress_time_counter <- function(start = 0, by = 0.6) {
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

test_that("npindexbw adopts the generic bandwidth selection line", {
  set.seed(42)
  x1 <- runif(30)
  x2 <- runif(30)
  y <- x1 + 0.5 * x2 + rnorm(30, sd = 0.1)

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    npindexbw(
      xdat = data.frame(x1 = x1, x2 = x2),
      ydat = y,
      method = "ichimura",
      nmulti = 3,
      optim.maxit = 100
    ),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  set.seed(42)
  single_line <- capture_progress_shadow_trace(
    npindexbw(
      xdat = data.frame(x1 = x1, x2 = x2),
      ydat = y,
      method = "ichimura",
      nmulti = 3,
      optim.maxit = 100
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(single_line)

  expect_s3_class(single_line$value, "sibandwidth")
  expect_equal(shadow_signature(single_line), shadow_signature(legacy))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 1/3\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 1/3, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/3, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/3, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 3/3, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})

test_that("npindexbw single-line progress respects np.messages FALSE and suppressMessages", {
  set.seed(7)
  x1 <- runif(25)
  x2 <- runif(25)
  y <- x1 + 0.5 * x2 + rnorm(25, sd = 0.1)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  silent <- capture_progress_shadow_trace(
    npindexbw(
      xdat = data.frame(x1 = x1, x2 = x2),
      ydat = y,
      method = "ichimura",
      nmulti = 3,
      optim.maxit = 100
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  expect_length(silent$trace, 0L)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(9)
  x1 <- runif(25)
  x2 <- runif(25)
  y <- x1 + 0.5 * x2 + rnorm(25, sd = 0.1)

  suppressed <- capture_progress_shadow_trace(
    suppressMessages(
      npindexbw(
        xdat = data.frame(x1 = x1, x2 = x2),
        ydat = y,
        method = "ichimura",
        nmulti = 3,
        optim.maxit = 100
      )
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  expect_length(suppressed$trace, 0L)
})
