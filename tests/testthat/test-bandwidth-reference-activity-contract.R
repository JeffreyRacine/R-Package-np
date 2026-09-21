test_that("reference heartbeats do not invent or reset objective counts", {
  ns <- asNamespace("np")
  select <- get(".np_progress_select_bandwidth", ns)
  activity <- get(".np_progress_bandwidth_activity_step", ns)
  runtime <- get(".np_progress_runtime", ns)
  old <- options(np.messages = TRUE, np.progress.bandwidth.enhanced = TRUE,
                 np.progress.start.grace.unknown.sec = 0)
  on.exit(options(old))
  values <- list()
  trace <- capture_progress_shadow_trace({
    select("Bandwidth selection", {
      activity(0L, force = TRUE)
      values$initial <- runtime$bandwidth_state$last_done
      activity(7L, force = TRUE)
      activity(0L, force = TRUE)
      values$retained <- runtime$bandwidth_state$last_done
    })
  }, now = local({time <- 0; function() {time <<- time + 2; time}}))
  expect_null(values$initial)
  expect_identical(values$retained, 7L)
  lines <- vapply(trace$trace, `[[`, "", "line")
  expect_true(any(grepl("iteration 7", lines, fixed = TRUE)))
  expect_false(any(grepl("iteration 0", lines, fixed = TRUE)))
})

test_that("all raw objective calls have an exception-safe activity owner", {
  path <- testthat::test_path("..", "..", "src", "np.c")
  if (!file.exists(path))
    skip("Native source is not present in this installed check")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  # These are the two canonical invocation sites, not eight independently
  # guarded callers. Reference, final verification and counted calls share it.
  raw <- gregexpr("bwmfunc_raw\\(", text)[[1L]]
  expect_length(raw[raw > 0], 2L)
  expect_match(text, "call->value=bwmfunc_raw(call->point)", fixed = TRUE)
  expect_match(text, "if(bwm_progress_eval_active)return bwmfunc_raw(point)",
               fixed = TRUE)
  expect_match(text,
    "R_ExecWithCleanup(bwm_raw_activity_execute,&call,bwm_raw_activity_cleanup,&call)",
    fixed = TRUE)
  expect_match(text, "bwm_progress_eval_active=call->prior_active", fixed = TRUE)
  expect_match(text, "current_eval < 1 && !bwm_progress_eval_active", fixed = TRUE)
})
