progress_completion_clock <- function() { t <- 0; function() { t <<- t + 1; t } }
progress_completion_fn <- function(name) getFromNamespace(name, "np")

test_that("completed starts retain activity until the bandwidth operation returns", {
  old <- options(np.messages=TRUE, np.progress.start.grace.unknown.sec=0)
  on.exit(options(old), add=TRUE)
  select <- progress_completion_fn(".np_progress_select_bandwidth_enhanced")
  multi <- progress_completion_fn(".np_progress_bandwidth_multistart_step")
  activity <- progress_completion_fn(".np_progress_bandwidth_activity_step")
  shadow <- capture_progress_shadow_trace(select("probe", {
    multi(1L, 2L); activity(64L, force=TRUE)
    multi(2L, 2L); activity(128L, force=TRUE)
    17L
  }), force_renderer="single_line", now=progress_completion_clock())
  expect_identical(shadow$value, 17L)
  active <- Filter(function(x) x$event == "render", shadow$trace)
  lines <- vapply(active, `[[`, "", "line")
  expect_false(any(grepl("100[.]0%|eta 0[.]0s", lines)))
  expect_true(any(grepl("multistart 2/2, starts complete; finishing, iteration 128", lines)))
  expect_match(shadow$final_line, "100[.]0%.*eta 0[.]0s")
})

test_that("bandwidth errors abort rather than report successful completion", {
  old <- options(np.messages=TRUE, np.progress.start.grace.unknown.sec=0)
  on.exit(options(old), add=TRUE)
  select <- progress_completion_fn(".np_progress_select_bandwidth_enhanced")
  multi <- progress_completion_fn(".np_progress_bandwidth_multistart_step")
  shadow <- capture_progress_shadow_trace(tryCatch(select("probe", {
    multi(1L, 2L); stop("expected reporting probe")
  }), error=function(e) conditionMessage(e)), force_renderer="single_line",
  now=progress_completion_clock())
  expect_identical(shadow$value, "expected reporting probe")
  expect_false(any(vapply(shadow$trace,function(x)x$event=="finish",logical(1))))
  expect_identical(tail(shadow$trace, 1L)[[1L]]$event, "abort")
  expect_true(is.null(progress_completion_fn(".np_progress_runtime")$bandwidth_state))
  expect_true(getOption("np.messages"))
  expect_false(withVisible(select("probe", invisible(17L)))$visible)
})

test_that("elapsed estimates and completed coordinator groups are not false ETAs", {
  state <- list(started=0, pkg_prefix="[np]", bandwidth_nmulti_total=10L,
    bandwidth_multistart_completed=1L, bandwidth_multistart_current=2L,
    bandwidth_multistart_durations=1, last_done=64L, bandwidth_coordinator_active=FALSE)
  fmt <- progress_completion_fn(".np_progress_bandwidth_format_estimate")
  expect_match(fmt(state, now=11), "iteration 64.*99[.]9%.*eta estimating")
  state$bandwidth_coordinator_active <- TRUE
  state$bandwidth_coordinator_offset <- 0L
  state$bandwidth_coordinator_local_total <- 2L
  state$bandwidth_coordinator_local_current <- 2L
  state$bandwidth_coordinator_group_label <- "first component"
  state$bandwidth_multistart_completed <- 2L
  expect_match(fmt(state, now=11), "first component.*starts complete; finishing.*iteration 64")
  expect_false(grepl("%|eta", fmt(state, now=11)))
})
