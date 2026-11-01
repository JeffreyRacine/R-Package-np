test_that("corrected adapters join each planned target exactly once", {
  old <- options(np.messages=TRUE, np.plot.progress=TRUE)
  on.exit(options(old), add=TRUE)
  now <- 0
  z <- capture_progress_shadow_trace({
    np:::.np_plot_progress_run({
      np:::.np_plot_progress_plan(2L)
      for (i in 1:2) {
        now <- now + 3
        value <- np:::.np_plot_corrected_activity_run("Corrected intervals", {
          expect_identical(np:::.np_plot_progress_runtime$context$target, i)
          now <- now + 3
          np:::.np_plot_progress_notify(force=TRUE)
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
  expect_null(np:::.np_plot_progress_runtime$context)
})

test_that("corrected target errors abort the parent rather than announce success", {
  old <- options(np.messages=TRUE,np.plot.progress=TRUE)
  on.exit(options(old), add=TRUE)
  z <- capture_progress_shadow_trace({
    expect_error(np:::.np_plot_progress_run({
      np:::.np_plot_progress_plan(1L)
      np:::.np_plot_corrected_activity_run("Corrected intervals", stop("sentinel"))
    }, enabled=TRUE), "sentinel")
  }, force_renderer="single_line")
  expect_false(any(vapply(z$trace,function(x)x$event=="finish",TRUE)))
  expect_null(np:::.np_plot_progress_runtime$context)
})


test_that("outer ETA uses completed samples and ignores inner counter resets", {
  old <- options(np.messages=TRUE,np.plot.progress=TRUE)
  on.exit(options(old), add=TRUE)
  now <- 0
  z <- capture_progress_shadow_trace({
    np:::.np_plot_progress_run({
      np:::.np_plot_progress_plan(1L)
      np:::.np_plot_corrected_activity_run("Corrected intervals", {
        now <- 10
        np:::.np_plot_progress_outer_run({
          now <- 12
          np:::.np_plot_progress_notify(done=1L,total=39L,unit="samples",force=TRUE)
          now <- 14
          np:::.np_plot_progress_notify(done=2L,total=39L,unit="samples",force=TRUE)
          np:::.np_plot_progress_notify(done=99L,total=99L,unit="rep",force=TRUE)
          np:::.np_plot_progress_notify(done=1L,total=39L,unit="samples",force=TRUE)
          expect_equal(np:::.np_plot_progress_runtime$context$phase$done,2)
          now <- 88
          np:::.np_plot_progress_notify(done=39L,total=39L,unit="samples",force=TRUE)
        })
        expect_null(np:::.np_plot_progress_runtime$context$phase)
        np:::.np_plot_progress_notify(force=TRUE)
      })
    }, enabled=TRUE)
  }, force_renderer="single_line", now=function() now)
  lines <- vapply(z$trace,`[[`,"","line")
  expect_true(any(grepl("samples 1/39.*elapsed 12.0s, outer eta 76.0s", lines)))
  expect_true(any(grepl("samples 2/39.*elapsed 14.0s, outer eta 74.0s", lines)))
  expect_false(any(grepl("rep 99/99",lines)))
  expect_true(any(grepl("samples 39/39.*outer eta finishing",lines)))
  expect_true(any(grepl("detail constructing intervals",lines)))
  expect_null(np:::.np_plot_progress_runtime$context)
})
