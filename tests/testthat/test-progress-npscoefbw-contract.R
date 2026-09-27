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
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, 50\\.0%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", bandwidth_lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, iteration [0-9]+, [0-9]+\\.[0-9]%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", bandwidth_lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, 100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", bandwidth_lines)))
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

test_that("smooth-coefficient search forwards kernel activity without row counters", {
  x <- data.frame(z = seq(-1, 1, length.out = 96))
  bw <- npscoefbw(xdat = x, zdat = x, ydat = sin(x$z), bws = 0.4,
                  bandwidth.compute = FALSE)
  old <- options(np.messages = TRUE)
  on.exit(options(old), add = TRUE)
  previous <- .np_progress_runtime$fit_forward
  previous.active <- .np_progress_runtime$scoef_search_activity
  result <- capture_progress_shadow_trace(
    .np_progress_select_bandwidth_enhanced("Bandwidth selection", {
      .np_progress_bandwidth_activity_step(done = 7L, force = TRUE)
      value <- .npscoefbw_search_activity(nrow(x), {
        .npscoefbw_search_activity(nrow(x),
          .np_estimator_loo_ksum(txdat = x, tydat = x$z,
                                bws = bw, leave.one.out = TRUE)$ksum)
      })
      expect_identical(.np_progress_runtime$bandwidth_state$last_done, 7L)
      expect_identical(.np_progress_runtime$fit_forward, previous)
      expect_identical(.np_progress_runtime$scoef_search_activity, previous.active)
      expect_error(.npscoefbw_search_activity(nrow(x), stop("scope sentinel")),
                   "scope sentinel", fixed = TRUE)
      expect_identical(.np_progress_runtime$fit_forward, previous)
      expect_identical(.np_progress_runtime$scoef_search_activity, previous.active)
      value
    }),
    force_renderer = "single_line", now = progress_time_counter()
  )
  lines <- shadow_render_lines(result)
  expect_gte(sum(grepl("iteration 7", lines, fixed = TRUE)), 2L)
  expect_false(any(grepl("Fitting", lines, fixed = TRUE)))
  options(np.messages = FALSE)
  expected <- .np_estimator_loo_ksum(txdat = x, tydat = x$z,
                                     bws = bw, leave.one.out = TRUE)$ksum
  expect_identical(result$value, expected)
})

test_that("smooth-coefficient activity has a quiet no-owner path", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  sentinel <- function() invisible(NULL)
  previous <- .np_progress_runtime$fit_forward
  on.exit(.np_progress_runtime$fit_forward <- previous, add = TRUE)
  .np_progress_runtime$fit_forward <- sentinel
  expect_identical(.npscoefbw_search_activity(96L, 17L), 17L)
  expect_identical(.npscoefbw_search_activity(0L, 18L), 18L)
  expect_identical(.np_progress_runtime$fit_forward, sentinel)
  expect_false(withVisible(.npscoefbw_search_activity(96L, invisible(17L)))$visible)
})

test_that("ordinary search renews native activity for successive objectives", {
  set.seed(42)
  x <- data.frame(x = runif(96))
  y <- x$x + rnorm(96)
  calls <- 0L
  original <- .np_with_compiled_fit_progress
  old <- options(np.messages = TRUE)
  on.exit(options(old), add = TRUE)
  with_np_progress_bindings(list(.np_with_compiled_fit_progress = function(...) {
    calls <<- calls + 1L
    original(...)
  }), capture_progress_shadow_trace(
    npscoefbw(xdat = x, zdat = x, ydat = y, bws = .4,
                nmulti = 1L, optim.maxit = 1L),
    force_renderer = "single_line", now = progress_time_counter()
  ))
  expect_gt(calls, 1L)
})

test_that("fits nested in search renew activity at each kernel block", {
  x <- data.frame(x = seq(-1, 1, length.out = 96))
  y <- sin(x$x)
  bw <- npscoefbw(xdat = x, zdat = x, ydat = y, bws = .4,
                  bandwidth.compute = FALSE)
  calls <- 0L
  original <- .np_with_compiled_fit_progress
  old <- options(np.messages = TRUE)
  on.exit(options(old), add = TRUE)
  result <- with_np_progress_bindings(list(.np_with_compiled_fit_progress = function(...) {
    calls <<- calls + 1L
    original(...)
  }), capture_progress_shadow_trace(
    .np_progress_select_bandwidth_enhanced("Bandwidth selection", {
      .np_progress_bandwidth_activity_step(done = 7L, force = TRUE)
      value <- npscoef(bws = bw, txdat = x, tzdat = x, tydat = y,
                        se = TRUE, .np_fit_progress_allow = FALSE)
      expect_identical(.np_progress_runtime$bandwidth_state$last_done, 7L)
      value
    }),
    force_renderer = "single_line", now = progress_time_counter()
  ))
  expect_gte(calls, 2L)
  expect_false(any(grepl("Fitting", shadow_lines(result), fixed = TRUE)))
  options(np.messages = FALSE)
  control <- npscoef(bws = bw, txdat = x, tzdat = x, tydat = y, se = TRUE)
  expect_identical(fitted(result$value), fitted(control))
  expect_identical(se(result$value), se(control))
})
