test_that("plot bootstrap display owns one clock without changing local clocks", {
  old.options <- options(np.messages = TRUE, np.plot.progress = TRUE)
  on.exit(options(old.options), add = TRUE)
  now <- 0
  trace <- capture_progress_shadow_trace({
    result <- npRmpi:::.np_plot_progress_run({
      npRmpi:::.np_plot_progress_plan(2L)
      target <- npRmpi:::.np_plot_progress_target_begin()
      now <- 2
      local <- npRmpi:::.np_plot_bootstrap_progress_begin(10L, "Plot bootstrap")
      expect_equal(local$started, 2)
      throttle <- local$throttle_sec
      now <- 5
      local <- npRmpi:::.np_plot_progress_tick(local, 5L, force = TRUE)
      expect_equal(local$started, 2)
      expect_identical(local$throttle_sec, throttle)
      npRmpi:::.np_plot_progress_end(local)
      npRmpi:::.np_plot_progress_target_end(target)
      target <- npRmpi:::.np_plot_progress_target_begin()
      now <- 8
      local <- npRmpi:::.np_plot_bootstrap_progress_begin(10L, "Plot bootstrap")
      expect_equal(local$started, 8)
      npRmpi:::.np_fit_progress_step(1L)
      now <- 12
      local <- npRmpi:::.np_plot_progress_tick(local, 10L, force = TRUE)
      npRmpi:::.np_plot_progress_end(local)
      npRmpi:::.np_plot_progress_target_end(target)
      17L
    }, enabled = TRUE)
    expect_identical(result, 17L)
  }, force_renderer = "single_line", now = function() now)
  events <- trace$trace
  expect_length(unique(vapply(events, `[[`, "", "id")), 1L)
  expect_true(all(vapply(events, `[[`, 0, "started_at") == 0))
  expect_equal(sum(vapply(events, function(x) x$event == "finish", TRUE)), 1L)
  expect_match(events[[1]]$line, "eta estimating")
  expect_true(any(vapply(events, function(x) grepl("2/2.*elapsed 12.0s", x$line), TRUE)))
  expect_null(npRmpi:::.np_plot_progress_runtime$context)
  expect_null(npRmpi:::.np_progress_runtime$fit_forward)
})

test_that("nested and failed plots restore progress scope without false completion", {
  old.options <- options(np.messages = TRUE, np.plot.progress = TRUE)
  on.exit(options(old.options), add = TRUE)
  now <- 0
  result <- capture_progress_shadow_trace({
    expect_error(npRmpi:::.np_plot_progress_run({
      npRmpi:::.np_plot_progress_plan(3L)
      outer <- npRmpi:::.np_plot_progress_runtime$context
      npRmpi:::.np_plot_progress_run({
        npRmpi:::.np_plot_progress_plan(99L)
        expect_null(npRmpi:::.np_plot_progress_runtime$context)
      }, enabled = TRUE)
      expect_identical(npRmpi:::.np_plot_progress_runtime$context, outer)
      expect_identical(outer$total, 3L)
      now <- 2
      stop("retained failure")
    }, enabled = TRUE), "retained failure")
    expect_null(npRmpi:::.np_plot_progress_runtime$context)
    expect_null(npRmpi:::.np_progress_runtime$fit_forward)
    expect_null(npRmpi:::.np_progress_registry$active_id)
  }, force_renderer = "single_line", now = function() now)
  expect_false(any(vapply(result$trace, function(x) x$event == "finish", TRUE)))
  expect_true(any(vapply(result$trace, function(x) x$event == "abort", TRUE)))
})

test_that("display aggregation leaves chunk-controller decisions unchanged", {
  old.options <- options(np.messages = TRUE, np.plot.progress = TRUE)
  on.exit(options(old.options), add = TRUE)
  sizes <- list()
  capture_progress_shadow_trace({
    for (aggregate in c(FALSE, TRUE)) {
      sizes[[length(sizes) + 1L]] <- npRmpi:::.np_plot_progress_run({
        progress <- npRmpi:::.np_plot_bootstrap_progress_begin(99L, "chunks")
        controller <- npRmpi:::.np_plot_progress_chunk_controller(16L, progress)
        inputs <- controller
        chunks <- integer()
        for (elapsed in c(.1, .4, 2, 4, .8)) {
          chunks <- c(chunks, controller$chunk.size)
          controller <- npRmpi:::.np_plot_progress_chunk_observe(
            controller, controller$chunk.size, elapsed)
        }
        npRmpi:::.np_plot_progress_end(progress)
        list(inputs = inputs, chunks = chunks, final = controller)
      }, enabled = aggregate)
    }
  }, force_renderer = "single_line")
  expect_identical(sizes[[1L]], sizes[[2L]])
})

test_that("quiet and standalone bootstrap helpers retain their local ownership", {
  old.options <- options(np.messages = TRUE, np.plot.progress = TRUE)
  on.exit(options(old.options), add = TRUE)
  trace <- capture_progress_shadow_trace({
    npRmpi:::.np_plot_progress_run({
      expect_null(npRmpi:::.np_plot_progress_runtime$context)
      local <- npRmpi:::.np_plot_bootstrap_progress_begin(3L, "standalone")
      expect_null(local$plot_context)
      npRmpi:::.np_plot_progress_end(local)
    }, enabled = FALSE)
  }, force_renderer = "single_line")
  expect_true(any(vapply(trace$trace, function(x) grepl("standalone", x$line), TRUE)))
  quiet <- capture_progress_shadow_trace(npRmpi:::.np_plot_progress_run({
    expect_null(npRmpi:::.np_plot_progress_runtime$context)
    1L
  }, enabled = TRUE), force_renderer = "single_line", interactive = FALSE)
  expect_length(quiet$trace, 0L)
})

test_that("plot display scope preserves visibility and narrow-width work fields", {
  old.options <- options(np.messages = TRUE, np.plot.progress = TRUE)
  on.exit(options(old.options), add = TRUE)
  now <- 0
  trace <- capture_progress_shadow_trace({
    value <- withVisible(npRmpi:::.np_plot_progress_run({
      npRmpi:::.np_plot_progress_plan(42L)
      target <- npRmpi:::.np_plot_progress_target_begin()
      now <- 12
      npRmpi:::.np_plot_progress_notify(done = 50L, total = 99L, force = TRUE)
      npRmpi:::.np_plot_progress_target_end(target)
      invisible(7L)
    }, enabled = TRUE))
    expect_false(value$visible)
    expect_identical(value$value, 7L)
  }, force_renderer = "single_line", now = function() now)
  lines <- vapply(trace$trace, `[[`, "", "line")
  expect_true(any(grepl("rep 50/99.*elapsed 12.0s.*eta", lines)))
  expect_null(npRmpi:::.np_plot_progress_runtime$context)
})

test_that("verbose target labels cannot displace progress counters at narrow widths", {
  old.options <- options(np.messages = TRUE, np.plot.progress = TRUE)
  on.exit(options(old.options), add = TRUE)
  for (width in c(60L, 80L, 120L)) {
    now <- 0
    result <- with_nprmpi_progress_bindings(
      list(.np_progress_output_width = function() width),
      capture_progress_shadow_trace(npRmpi:::.np_plot_progress_run({
        npRmpi:::.np_plot_progress_plan(42L)
        npRmpi:::.np_plot_progress_target_begin()
        now <- 12
        npRmpi:::.np_plot_progress_notify(
          stage = paste(rep("long target label", 10L), collapse = " "),
          done = 50L, total = 99L, force = TRUE)
      }, enabled = TRUE), force_renderer = "single_line", now = function() now))
    lines <- vapply(result$trace, `[[`, "", "render_line")
    expect_true(all(nchar(lines, type = "width") <= width))
    expect_true(any(grepl("1/42.*rep 50/99.*elap(sed)? 12.0s.*eta", lines)))
  }
})
