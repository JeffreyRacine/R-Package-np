bootstrap_progress_get <- function(name) getFromNamespace(name, "npRmpi")
bootstrap_progress_clock <- function() {
  elapsed <- 0
  function() {
    elapsed <<- elapsed + .25
    elapsed
  }
}

test_that("profile workers cannot own bootstrap progress", {
  old <- options(np.messages = TRUE, np.progress.interval.sec = 0)
  on.exit(options(old), add = TRUE)
  local_mocked_bindings(.np_progress_bandwidth_worker_silent = function() TRUE,
    .package = "npRmpi")
  runtime <- bootstrap_progress_get(".np_progress_runtime")
  prior <- runtime$fit_forward
  actual <- capture_progress_shadow_trace(local({
    context <- bootstrap_progress_get(".np_bootstrap_progress_begin")(3L, "Bootstrap")
    on.exit(bootstrap_progress_get(".np_progress_activity_end")(context), add = TRUE)
    expect_false(context$state$visible)
    bootstrap_progress_get(".np_fit_progress_begin")("Child fit", 2L)
    bootstrap_progress_get(".np_fit_progress_step")(1L)
    bootstrap_progress_get(".np_fit_progress_finish")()
    child <- bootstrap_progress_get(".np_progress_begin")("Child bandwidth", 1L)
    expect_false(child$visible)
    bootstrap_progress_get(".np_progress_show_now")(child, done = 0L)
    bootstrap_progress_get(".np_progress_end")(child)
  }), force_renderer = "single_line", now = bootstrap_progress_clock())
  expect_length(actual$trace, 0L)
  expect_identical(runtime$fit_forward, prior)
})

test_that("disabled bootstrap progress preserves the statistic and parent", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  runtime <- bootstrap_progress_get(".np_progress_runtime")
  parent <- runtime$fit_forward
  statistic <- function(data, indices) sum(data[indices])
  context <- bootstrap_progress_get(".np_bootstrap_progress_begin")(3L, "Bootstrap")
  expect_null(context)
  expect_identical(bootstrap_progress_get(".np_bootstrap_progress_statistic")(
    context, statistic), statistic)
  bootstrap_progress_get(".np_progress_activity_end")(context)
  expect_identical(runtime$fit_forward, parent)
})

test_that("bootstrap progress counts t0 separately and forwards fit activity", {
  old <- options(np.messages = TRUE, np.progress.interval.sec = .1)
  on.exit(options(old), add = TRUE)
  begin <- bootstrap_progress_get(".np_bootstrap_progress_begin")
  end <- bootstrap_progress_get(".np_progress_activity_end")
  wrap <- bootstrap_progress_get(".np_bootstrap_progress_statistic")
  fit.begin <- bootstrap_progress_get(".np_fit_progress_begin")
  fit.step <- bootstrap_progress_get(".np_fit_progress_step")
  fit.end <- bootstrap_progress_get(".np_fit_progress_finish")
  runtime <- bootstrap_progress_get(".np_progress_runtime")
  prior <- runtime$fit_forward
  actual <- capture_progress_shadow_trace(local({
    context <- begin(3L, "Bootstrapping single-index fit")
    on.exit(end(context), add = TRUE)
    calls <- 0L
    statistic <- wrap(context, function(data, indices) {
      calls <<- calls + 1L
      fit.begin("Fitting regression", 4L)
      for (i in 1:4) fit.step(i)
      fit.end()
      sum(data[indices])
    })
    for (i in 0:3) expect_identical(statistic(1:4, 4:1), 10L)
    expect_identical(calls, 4L)
    expect_identical(context$done, 3L)
    end(context, completed = TRUE)
  }), force_renderer = "single_line", now = bootstrap_progress_clock())
  trace <- actual$trace
  lines <- vapply(trace, `[[`, character(1L), "line")
  details <- vapply(trace, function(x) if (is.null(x$detail)) "" else x$detail, "")
  events <- vapply(trace, `[[`, character(1L), "event")
  expect_false(any(grepl("Fitting regression", lines, fixed = TRUE)))
  expect_true(any(details == "initial statistic"))
  for (i in 1:3)
    expect_gte(sum(details == sprintf("replication %d of 3", i)), 4L)
  expect_equal(length(unique(vapply(trace, `[[`, "", "id"))), 1L)
  expect_equal(sum(events == "finish"), 1L)
  expect_identical(runtime$fit_forward, prior)
})

test_that("bootstrap scopes restore the parent and abort exactly once", {
  old <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0,
                 np.progress.interval.sec = 0)
  on.exit(options(old), add = TRUE)
  begin <- bootstrap_progress_get(".np_bootstrap_progress_begin")
  end <- bootstrap_progress_get(".np_progress_activity_end")
  runtime <- bootstrap_progress_get(".np_progress_runtime")
  prior <- runtime$fit_forward
  err <- simpleError("bootstrap witness")
  actual <- capture_progress_shadow_trace(local({
    outer <- begin(3L, "Outer bootstrap")
    on.exit(end(outer), add = TRUE)
    parent <- runtime$fit_forward
    inner <- begin(2L, "Inner bootstrap")
    end(inner, completed = TRUE)
    end(inner)
    expect_identical(runtime$fit_forward, parent)
    seen <- tryCatch(local({
      on.exit(end(outer), add = TRUE)
      bootstrap_progress_get(".np_fit_progress_begin")("Fitting regression", 2L)
      stop(err)
    }), error = identity)
    expect_identical(seen, err)
    expect_identical(runtime$fit_forward, prior)
    end(outer)
    # Fresh ordinary fits must use their own progress again after failure.
    bootstrap_progress_get(".np_fit_progress_begin")("Fresh fit", 2L)
    bootstrap_progress_get(".np_fit_progress_step")(1L)
    bootstrap_progress_get(".np_fit_progress_finish")()
  }), force_renderer = "single_line", now = bootstrap_progress_clock())
  events <- vapply(actual$trace, `[[`, "", "event")
  lines <- vapply(actual$trace, `[[`, "", "line")
  expect_equal(sum(events == "abort"), 1L)
  expect_true(any(grepl("Fresh fit", lines, fixed = TRUE)))
  expect_identical(runtime$fit_forward, prior)
  expect_null(runtime$fit_state)
})
