progress_time_counter <- function(start = 0, by = 2.1) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

shadow_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

test_that("npscoefbw uses single-line multistart with legacy objective progress", {
  set.seed(3240)
  n <- 65
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + x * (1 + z) + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  legacy <- capture_progress_shadow_trace(
    npscoefbw(
      xdat = data.frame(x = x),
      zdat = data.frame(z = z),
      ydat = y,
      regtype = "lc",
      nmulti = 2,
      optim.maxit = 10,
      cv.iterate = FALSE
    ),
    force_renderer = "legacy",
    now = progress_time_counter()
  )

  set.seed(3240)
  actual <- capture_progress_shadow_trace(
    npscoefbw(
      xdat = data.frame(x = x),
      zdat = data.frame(z = z),
      ydat = y,
      regtype = "lc",
      nmulti = 2,
      optim.maxit = 10,
      cv.iterate = FALSE
    ),
    now = progress_time_counter()
  )

  bandwidth_lines <- shadow_lines(actual)[grepl("^\\[np\\] Selecting smooth coefficient bandwidth", shadow_lines(actual))]
  legacy_bandwidth_lines <- shadow_lines(legacy)[grepl("^\\[np\\] Selecting smooth coefficient bandwidth", shadow_lines(legacy))]

  expect_s3_class(actual$value, "scbandwidth")
  expect_equal(bandwidth_lines, legacy_bandwidth_lines)
  expect_true(any(grepl("^\\[np\\] Selecting smooth coefficient bandwidth multistart 1/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", bandwidth_lines)))
  expect_true(any(grepl("^\\[np\\] Selecting smooth coefficient bandwidth multistart 2/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", bandwidth_lines)))
  expect_true(any(grepl("^\\[np\\] Optimizing smooth coefficient bandwidth\\.\\.\\. iteration [0-9]+, elapsed [0-9]+\\.[0-9]s: multistart 1$", shadow_lines(actual))))
})

test_that("npscoefbw progress respects np.messages FALSE", {
  set.seed(3240)
  n <- 50
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + x * (1 + z) + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  silent <- capture_progress_shadow_trace(
    npscoefbw(
      xdat = data.frame(x = x),
      zdat = data.frame(z = z),
      ydat = y,
      regtype = "lc",
      nmulti = 1,
      optim.maxit = 5,
      cv.iterate = FALSE
    ),
    now = progress_time_counter()
  )

  expect_length(silent$trace, 0L)
})

test_that("npscoefbw cv.iterate path retains legacy backfitting progress semantics", {
  set.seed(3240)
  n <- 65
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + x * (1 + z) + rnorm(n, sd = 0.1)

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npscoefbw(
      xdat = data.frame(x = x),
      zdat = data.frame(z = z),
      ydat = y,
      regtype = "lc",
      nmulti = 1,
      optim.maxit = 5,
      cv.iterate = TRUE,
      cv.num.iterations = 2,
      backfit.iterate = FALSE
    ),
    now = progress_time_counter()
  )

  lines <- shadow_lines(actual)

  expect_s3_class(actual$value, "scbandwidth")
  expect_true(any(grepl("^\\[np\\] Backfitting smooth coefficient bandwidth 1/2 \\([0-9]+\\.[0-9]%.*, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): iteration 1 of 2$", lines)))
  expect_true(any(grepl("^\\[np\\] Optimizing partial residual bandwidth\\.\\.\\. iteration [0-9]+, elapsed [0-9]+\\.[0-9]s: backfitting iteration [0-9]+ of 2, partial residual [0-9]+ of 2, fval ", lines)))
})
