test_that("NOMAD Powell handoff emits a refining line immediately", {
  reset <- getFromNamespace(".np_progress_reset_registry", "np")
  reset()
  on.exit(reset(), add = TRUE)

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.unknown.sec = 10,
    np.progress.interval.unknown.sec = 10
  )
  on.exit(options(old_opts), add = TRUE)

  progress_now <- local({
    values <- c(0, 0.5)
    i <- 0L
    function() {
      i <<- min(i + 1L, length(values))
      values[[i]]
    }
  })

  trace <- capture_progress_shadow_trace(
    {
      begin <- getFromNamespace(".np_nomad_progress_begin", "np")
      enter <- getFromNamespace(".np_nomad_progress_enter_powell", "np")
      fmt <- getFromNamespace(".np_progress_format_line", "np")

      state <- begin(
        nmulti = 2L,
        baseline_degree = 1L,
        best_record = list(degree = 1L, objective = 1)
      )
      state$last_done <- 10L
      state$nomad_eval_id <- 10L
      state$last_line <- fmt(state = state, done = state$last_done, now = 0)
      state$rendered <- TRUE
      state$last_render_width <- nchar(state$last_line, type = "width")
      enter(
        state = state,
        degree = 1L,
        best_record = list(degree = 1L, objective = 1)
      )
    },
    force_renderer = "single_line",
    now = progress_now
  )

  events <- vapply(trace$trace, `[[`, character(1), "event")
  lines <- vapply(trace$trace, `[[`, character(1), "line")
  expect_true(any(events == "finish"))
  expect_true(any(grepl("^\\[np\\] Refining bandwidth \\(", lines)))
  expect_true(any(grepl("^\\[np\\] Selecting degree and bandwidth \\(", lines)))
})

test_that("RStudio single-line Powell handoff stays on the single-line surface", {
  reset <- getFromNamespace(".np_progress_reset_registry", "np")
  reset()
  on.exit(reset(), add = TRUE)

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  trace <- with_np_progress_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_rstudio_console = function() TRUE,
      .np_progress_now = function() 0
    ),
    capture_progress_shadow_trace(
      {
        begin <- getFromNamespace(".np_nomad_progress_begin", "np")
        enter <- getFromNamespace(".np_nomad_progress_enter_powell", "np")
        fmt <- getFromNamespace(".np_progress_format_line", "np")
        state <- begin(
          nmulti = 2L,
          baseline_degree = 1L,
          best_record = list(degree = 1L, objective = 1)
        )
        state$last_done <- 10L
        state$nomad_eval_id <- 10L
        state$last_line <- fmt(state = state, done = state$last_done, now = 0)
        state$rendered <- TRUE
        state$last_render_width <- nchar(state$last_line, type = "width")
        enter(
          state = state,
          degree = 1L,
          best_record = list(degree = 1L, objective = 1)
        )
      },
      force_renderer = "single_line",
      now = function() 0
    )
  )

  events <- vapply(trace$trace, `[[`, character(1), "event")
  lines <- vapply(trace$trace, `[[`, character(1), "line")
  expect_true(any(events == "finish"))
  expect_true(any(grepl("^\\[np\\] Refining bandwidth \\(", lines)))
  expect_false(any(grepl("^\\[np\\] Selecting degree and bandwidth \\(", lines[events != "finish"])))
})
