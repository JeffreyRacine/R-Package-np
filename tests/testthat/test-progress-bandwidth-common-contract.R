progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

shadow_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

test_that("dark-launched bandwidth engine preserves nmulti=1 iteration heartbeats", {
  select_bw <- getFromNamespace(".np_progress_select_bandwidth", "np")
  activity_bw <- getFromNamespace(".np_progress_bandwidth_activity_step", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.progress.bandwidth.enhanced = TRUE,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      value <- select_bw("Selecting density bandwidth", {
        activity_bw(28L)
        activity_bw(64L)
        7
      })
      expect_identical(value, 7)
    },
    now = progress_time_values(c(0, 1.2, 3.4, 5.6))
  )

  lines <- shadow_lines(actual)

  expect_true(any(grepl("^\\[np\\] Bandwidth selection\\.\\.\\.$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(iteration 28, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(iteration 64, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_false(any(grepl("multistart", lines, fixed = TRUE)))
})

test_that("dark-launched bandwidth engine switches from iteration to estimate mode after first completion", {
  select_bw <- getFromNamespace(".np_progress_select_bandwidth", "np")
  set_total <- getFromNamespace(".np_progress_bandwidth_set_total", "np")
  activity_bw <- getFromNamespace(".np_progress_bandwidth_activity_step", "np")
  step_bw <- getFromNamespace(".np_progress_bandwidth_multistart_step", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.progress.bandwidth.enhanced = TRUE,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      value <- select_bw("Selecting density bandwidth", {
        set_total(2L)
        activity_bw(28L)
        step_bw(1L, 2L)
        activity_bw(56L)
        activity_bw(84L)
        11
      })
      expect_identical(value, 11)
    },
    now = progress_time_values(c(0, 1.0, 4.0, 7.0, 10.0, 12.0))
  )

  lines <- shadow_lines(actual)

  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 1/2, iteration 28, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, elapsed 4\\.0s, 50\\.0%, eta 4\\.0s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, iteration 56, elapsed 7\\.0s, 87\\.5%, eta 1\\.0s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, iteration 84, elapsed 10\\.0s, 99\\.9%, eta 0\\.0s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, elapsed 12\\.0s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})

test_that("dark-launched completion-estimate heartbeats keep the unknown-total throttle cadence", {
  select_bw <- getFromNamespace(".np_progress_select_bandwidth", "np")
  set_total <- getFromNamespace(".np_progress_bandwidth_set_total", "np")
  activity_bw <- getFromNamespace(".np_progress_bandwidth_activity_step", "np")
  step_bw <- getFromNamespace(".np_progress_bandwidth_multistart_step", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.progress.bandwidth.enhanced = TRUE,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      value <- select_bw("Selecting regression bandwidth", {
        set_total(2L)
        activity_bw(20L)
        step_bw(1L, 2L)
        activity_bw(21L)
        activity_bw(22L)
        13
      })
      expect_identical(value, 13)
    },
    now = progress_time_values(c(0, 1.0, 4.0, 5.0, 6.2, 8.4))
  )

  complete_trace <- actual$trace[grepl("eta ", shadow_lines(actual), fixed = TRUE)]
  complete_lines <- vapply(complete_trace, `[[`, character(1L), "line")
  complete_times <- vapply(complete_trace, `[[`, numeric(1L), "now")

  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, elapsed 4\\.0s, 50\\.0%, eta 4\\.0s\\)$", complete_lines)))
  expect_false(any(grepl("elapsed 5\\.0s", complete_lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(multistart 2/2, iteration 22, elapsed 6\\.2s, 77\\.5%, eta 1\\.8s\\)$", complete_lines)))
  expect_true(all(diff(complete_times[seq_len(min(2L, length(complete_times)))]) >= 2.0))
})
