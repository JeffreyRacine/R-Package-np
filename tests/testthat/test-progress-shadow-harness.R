progress_shadow_times <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

test_that("shadow harness captures known-total render events and final line", {
  begin <- getFromNamespace(".np_progress_begin", "np")
  step <- getFromNamespace(".np_progress_step", "np")
  finish <- getFromNamespace(".np_progress_end", "np")

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.known.sec = 0)
  on.exit(options(old_opts), add = TRUE)

  shadow <- capture_progress_shadow_trace(
    {
      state <- begin("Bootstrap replications", total = 2)
      state <- step(state, done = 1)
      finish(state)
    },
    now = progress_shadow_times(c(0, 0.6, 1.2))
  )

  expect_equal(vapply(shadow$trace, `[[`, character(1L), "event"), c("render", "render", "finish"))
  expect_equal(vapply(shadow$trace, `[[`, character(1L), "kind"), rep("known", 3L))
  expect_identical(shadow$trace[[1L]]$line, "[np] Bootstrap replications...")
  expect_match(shadow$trace[[2L]]$line, "^\\[np\\] Bootstrap replications 1/2 ")
  expect_match(shadow$final_line, "^\\[np\\] Bootstrap replications 2/2 ")
})

test_that("shadow harness captures unknown-total delayed activity semantics", {
  begin <- getFromNamespace(".np_progress_begin", "np")
  step <- getFromNamespace(".np_progress_step", "np")

  old_opts <- options(np.messages = TRUE, np.progress.start.grace.unknown.sec = 1.0)
  on.exit(options(old_opts), add = TRUE)

  shadow <- capture_progress_shadow_trace(
    {
      state <- begin("Iterating Landweber-Fridman solve")
      state <- step(state, done = 3, detail = "updating residual smoothing")
      step(state, done = 4, detail = "evaluating stopping rule")
    },
    now = progress_shadow_times(c(0, 0.4, 1.2))
  )

  expect_equal(vapply(shadow$trace, `[[`, character(1L), "event"), c("render", "render"))
  expect_equal(vapply(shadow$trace, `[[`, character(1L), "kind"), rep("unknown", 2L))
  expect_identical(shadow$trace[[1L]]$line, "[np] Iterating Landweber-Fridman solve...")
  expect_identical(
    shadow$final_line,
    "[np] Iterating Landweber-Fridman solve... iteration 4, elapsed 1.2s: evaluating stopping rule"
  )
})
