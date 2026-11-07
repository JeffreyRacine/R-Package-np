test_that("standalone corrected work keeps a live inference display", {
  old <- options(np.messages = TRUE, np.plot.progress = TRUE)
  on.exit(options(old), add = TRUE)
  now <- 0
  z <- capture_progress_shadow_trace({
    value <- withVisible(np:::.np_plot_corrected_activity_run(
      "Constructing bias-corrected intervals", {
        expect_false(is.null(np:::.np_plot_progress_runtime$context))
        now <- 5
        np:::.np_plot_progress_notify(stage = "Constructing intervals",
          done = 1024L, total = 9999L, unit = "draws", force = TRUE)
        now <- 10
        np:::.np_progress_runtime$fit_forward()
        invisible(17L)
      }))
    expect_identical(value, list(value = 17L, visible = FALSE))
  }, force_renderer = "single_line", now = function() now)
  lines <- vapply(z$trace, `[[`, "", "line")
  expect_true(any(grepl("Plot inference.*elapsed 5.0s", lines)))
  expect_true(any(grepl("elapsed 10.0s", lines)))
  expect_false(any(grepl("Plot bootstrap", lines)))
  expect_length(unique(vapply(z$trace, `[[`, "", "id")), 1L)
  expect_equal(sum(vapply(z$trace, function(x) x$event == "finish", TRUE)), 1L)
  expect_null(np:::.np_plot_progress_runtime$context)
  expect_null(np:::.np_progress_runtime$fit_forward)
  expect_null(np:::.np_progress_registry$active_id)
})

test_that("standalone corrected failures and quiet calls leave no display owner", {
  old <- options(np.messages = TRUE, np.plot.progress = TRUE)
  on.exit(options(old), add = TRUE)
  for (condition in list(simpleError("sentinel"),
      structure(list(message = "interrupt sentinel", call = NULL),
                class = c("interrupt", "condition")))) {
    z <- capture_progress_shadow_trace({
      caught <- tryCatch(np:::.np_plot_corrected_activity_run("Corrected", {
        stop(condition)
      }), error = identity, interrupt = identity)
      expect_identical(caught, condition)
    }, force_renderer = "single_line")
    expect_false(any(vapply(z$trace, function(x) x$event == "finish", TRUE)))
    expect_null(np:::.np_plot_progress_runtime$context)
    expect_null(np:::.np_progress_runtime$fit_forward)
    expect_null(np:::.np_progress_registry$active_id)
  }
  options(np.plot.progress = FALSE)
  z <- capture_progress_shadow_trace(
    np:::.np_plot_corrected_activity_run("Corrected", 9L),
    force_renderer = "single_line")
  expect_identical(z$value, 9L)
  expect_length(z$trace, 0L)
})
