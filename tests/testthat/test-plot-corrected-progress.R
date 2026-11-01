test_that("corrected adapters join each planned target exactly once", {
  old <- options(np.messages=TRUE, np.plot.progress=TRUE)
  on.exit(options(old), add=TRUE)
  now <- 0
  z <- capture_progress_shadow_trace({
    npRmpi:::.np_plot_progress_run({
      npRmpi:::.np_plot_progress_plan(2L)
      for (i in 1:2) {
        now <- now + 3
        value <- npRmpi:::.np_plot_corrected_activity_run("Corrected intervals", {
          expect_identical(npRmpi:::.np_plot_progress_runtime$context$target, i)
          now <- now + 3
          npRmpi:::.np_plot_progress_notify(force=TRUE)
          19L
        })
        expect_identical(value, 19L)
      }
    }, enabled=TRUE)
  }, force_renderer="single_line", now=function() now)
  lines <- vapply(z$trace, `[[`, "", "line")
  expect_true(any(grepl("target 1/2", lines)))
  expect_true(any(grepl("target 2/2", lines)))
  expect_false(any(grepl("target 3/2", lines)))
  expect_length(unique(vapply(z$trace, `[[`, "", "id")), 1L)
  expect_null(npRmpi:::.np_plot_progress_runtime$context)
})

test_that("corrected target errors abort the parent rather than announce success", {
  old <- options(np.messages=TRUE,np.plot.progress=TRUE)
  on.exit(options(old), add=TRUE)
  z <- capture_progress_shadow_trace({
    expect_error(npRmpi:::.np_plot_progress_run({
      npRmpi:::.np_plot_progress_plan(1L)
      npRmpi:::.np_plot_corrected_activity_run("Corrected intervals", stop("sentinel"))
    }, enabled=TRUE), "sentinel")
  }, force_renderer="single_line")
  expect_false(any(vapply(z$trace,function(x)x$event=="finish",TRUE)))
  expect_null(npRmpi:::.np_plot_progress_runtime$context)
})
