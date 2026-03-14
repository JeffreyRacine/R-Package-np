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

shadow_render_lines <- function(shadow) {
  trace <- shadow$trace[vapply(shadow$trace, `[[`, character(1L), "event") == "render"]
  vapply(trace, `[[`, character(1L), "line")
}

installed_function_text <- function(name, package = "np") {
  paste(deparse(getFromNamespace(name, package)), collapse = "\n")
}

test_that("npscoefbw adopts the generic bandwidth selection line", {
  set.seed(3240)
  n <- 28
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
      optim.maxit = 3,
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
      optim.maxit = 3,
      cv.iterate = FALSE
    ),
    now = progress_time_counter()
  )

  bandwidth_lines <- shadow_lines(actual)[grepl("^\\[np\\] Bandwidth selection", shadow_lines(actual))]
  render_bandwidth_lines <- shadow_render_lines(actual)[grepl("^\\[np\\] Bandwidth selection", shadow_render_lines(actual))]
  legacy_bandwidth_lines <- shadow_render_lines(legacy)[grepl("^\\[np\\] Bandwidth selection", shadow_render_lines(legacy))]

  expect_s3_class(actual$value, "scbandwidth")
  expect_equal(render_bandwidth_lines, legacy_bandwidth_lines)
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 1/2\\)$", bandwidth_lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 1/2, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s\\)$", bandwidth_lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, elapsed [0-9]+\\.[0-9]s, 50\\.0%, eta [0-9]+\\.[0-9]s\\)$", bandwidth_lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, iteration [0-9]+, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", bandwidth_lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", bandwidth_lines)))
})

test_that("npscoefbw progress respects np.messages FALSE", {
  set.seed(3240)
  n <- 24
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
      optim.maxit = 2,
      cv.iterate = FALSE
    ),
    now = progress_time_counter()
  )

  expect_length(silent$trace, 0L)
})

test_that("npscoefbw cv.iterate path retains backfitting progress hooks", {
  src <- installed_function_text("npscoefbw.scbandwidth")

  expect_true(grepl("Backfitting smooth coefficient bandwidth", src, fixed = TRUE))
  expect_true(grepl("Optimizing partial residual bandwidth", src, fixed = TRUE))
  expect_true(grepl("\\.np_progress_begin\\(\"Backfitting smooth coefficient bandwidth\"", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Optimizing partial residual bandwidth\"", src))
})
