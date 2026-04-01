test_that("NOMAD Powell handoff emits a refining line immediately", {
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
      begin <- getFromNamespace(".np_nomad_progress_begin", "npRmpi")
      enter <- getFromNamespace(".np_nomad_progress_enter_powell", "npRmpi")

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
  expect_true(any(grepl("^\\[npRmpi\\] Refining bandwidth \\(", lines)))
  expect_false(any(grepl("^\\[npRmpi\\] Selecting degree and bandwidth \\(", lines)))
})
