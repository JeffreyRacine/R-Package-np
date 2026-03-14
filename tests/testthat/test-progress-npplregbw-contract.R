progress_time_counter <- function(start = 0, by = 1.7) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

shadow_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

test_that("npplregbw uses coordinated generic bandwidth selection progress", {
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

  actual <- capture_progress_shadow_trace(
    npplregbw(
      xdat = data.frame(x = x),
      zdat = data.frame(z = z),
      ydat = y,
      nmulti = 2
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(actual)

  expect_s3_class(actual$value, "plbandwidth")
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(y~z, multistart 1/2\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(y~z, multistart 2/2, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(x1~z, multistart 1/2, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(x1~z, multistart 2/2, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})

test_that("npplreg formula entry inherits coordinated generic bandwidth progress", {
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

  actual <- capture_progress_shadow_trace(
    npplreg(
      y ~ x | z,
      data = dat,
      nmulti = 2
    ),
    force_renderer = "single_line",
    now = progress_time_counter()
  )

  lines <- shadow_lines(actual)

  expect_s3_class(actual$value, "plregression")
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(y~z, multistart 1/2\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(x1~z, multistart 1/2, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(x1~z, multistart 2/2, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})
