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

      state <- begin(
        nmulti = 2L,
        baseline_degree = 1L,
        best_record = list(degree = 1L, objective = 1)
      )
      enter(
        state = state,
        degree = 1L,
        best_record = list(degree = 1L, objective = 1)
      )
    },
    force_renderer = "legacy",
    now = progress_now
  )

  lines <- vapply(trace$trace, `[[`, character(1), "line")
  expect_true(any(grepl("^\\[np\\] Refining bandwidth \\(", lines)))
  expect_false(any(grepl("^\\[np\\] Selecting degree and bandwidth \\(", lines)))
})

test_that("RStudio single-line Powell handoff switches to a persistent renderer", {
  reset <- getFromNamespace(".np_progress_reset_registry", "np")
  reset()
  on.exit(reset(), add = TRUE)

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  state <- with_np_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_is_rstudio_console = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "single_line",
      .np_progress_now = function() 0
    ),
    {
      begin <- getFromNamespace(".np_nomad_progress_begin", "np")
      enter <- getFromNamespace(".np_nomad_progress_enter_powell", "np")
      enter(
        state = begin(
          nmulti = 2L,
          baseline_degree = 1L,
          best_record = list(degree = 1L, objective = 1)
        ),
        degree = 1L,
        best_record = list(degree = 1L, objective = 1)
      )
    }
  )

  expect_identical(state$renderer, "legacy")
})
