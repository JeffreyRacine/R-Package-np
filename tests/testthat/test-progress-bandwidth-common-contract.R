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

test_that("bandwidth compaction never abbreviates the Bandwidth selection prefix", {
  compact_bandwidth <- getFromNamespace(".np_progress_compact_bandwidth_line", "np")

  compacted <- compact_bandwidth(
    "[np] Bandwidth selection (multistart 2/2, iteration 84, elapsed 10.0s, 99.9%, eta 0.0s)",
    max_width = 74L
  )

  expect_true(grepl("Bandwidth selection \\(", compacted))
  expect_false(grepl("Bandwidth sel \\(", compacted))
})

test_that("unknown-bound NOMAD restart detail drops synthetic percent and eta", {
  nomad_detail <- getFromNamespace(".np_nomad_progress_detail", "np")

  line <- nomad_detail(
    current_degree = 5L,
    best_record = list(degree = 1L),
    iteration = 77L,
    cumulative_iteration = 1234L,
    restart_index = 2L,
    nmulti = 2L,
    restart_durations = c(18.8),
    elapsed = 18.8
  )

  expect_identical(
    line,
    "multistart 2/2, iteration 77 (1234), elapsed 18.8s, deg (5), best (1)"
  )
  expect_false(grepl("%|eta ", line))
})

test_that("NOMAD progress fuses component context with dynamic fields", {
  nomad_detail <- getFromNamespace(".np_nomad_progress_detail", "np")
  powell_detail <- getFromNamespace(".np_nomad_powell_progress_detail", "np")
  npreg_powell_fields <- getFromNamespace(".npregbw_powell_progress_fields", "np")
  fit_line <- getFromNamespace(".np_progress_fit_single_line", "np")
  set_context <- getFromNamespace(".np_progress_bandwidth_set_context", "np")

  set_context("E[y|z] (1/3)")
  on.exit(set_context(NULL), add = TRUE)

  line <- nomad_detail(
    current_degree = 5L,
    best_record = list(degree = 1L),
    iteration = 77L,
    cumulative_iteration = 1234L,
    restart_index = 2L,
    nmulti = 2L,
    restart_durations = c(18.8),
    elapsed = 18.8
  )

  expect_identical(
    line,
    "E[y|z] (1/3), multistart 2/2, iteration 77 (1234), elapsed 18.8s, deg (5), best (1)"
  )

  compact <- fit_line(
    sprintf("[np] Selecting degree/bandwidth (%s)", line),
    max_width = 80L
  )

  expect_true(startsWith(compact, "[np] Degree/bw search (E[y|z] (1/3),"))
  expect_true(grepl("iter 77", compact, fixed = TRUE))
  expect_true(grepl("elap 18.8s", compact, fixed = TRUE))
  expect_true(grepl("deg (5)", compact, fixed = TRUE))
  expect_false(identical(compact, "[np] E[y|z] (1/3)"))

  powell <- powell_detail(
    current_degree = 5L,
    best_record = list(degree = 1L),
    iteration = 9L,
    elapsed = 2.5
  )

  expect_identical(
    powell,
    "E[y|z] (1/3), elapsed 2.5s, iter 9, deg (5), best (1)"
  )

  expect_identical(
    npreg_powell_fields(
      state = list(started = 0, nomad_current_degree = 5L),
      done = 9L,
      now = 2.5
    ),
    c("E[y|z] (1/3)", "elapsed 2.5s", "degree (5)", "iter 9")
  )
})

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

test_that("bandwidth progress can carry a coordinator context label such as degree", {
  select_bw <- getFromNamespace(".np_progress_select_bandwidth", "np")
  set_total <- getFromNamespace(".np_progress_bandwidth_set_total", "np")
  activity_bw <- getFromNamespace(".np_progress_bandwidth_activity_step", "np")
  step_bw <- getFromNamespace(".np_progress_bandwidth_multistart_step", "np")
  set_context <- getFromNamespace(".np_progress_bandwidth_set_context", "np")

  old_opts <- options(
    np.messages = TRUE,
    np.progress.bandwidth.enhanced = TRUE,
    np.progress.start.grace.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    {
      set_context("deg (1,0)")
      on.exit(set_context(NULL), add = TRUE)
      value <- select_bw("Selecting regression bandwidth", {
        set_total(2L)
        activity_bw(20L)
        step_bw(1L, 2L)
        activity_bw(21L)
        13
      })
      expect_identical(value, 13)
    },
    now = progress_time_values(c(0, 1.0, 4.0, 6.2, 8.4))
  )

  lines <- shadow_lines(actual)

  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(deg \\(1,0\\)\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(deg \\(1,0\\), multistart 1/2, iteration 20, elapsed [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(deg \\(1,0\\), multistart 2/2, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(deg \\(1,0\\), multistart 2/2, iteration 21, elapsed [0-9]+\\.[0-9]s, [0-9]+\\.[0-9]%, eta [0-9]+\\.[0-9]s\\)$", lines)))
  expect_true(any(grepl("^\\[np\\] Bandwidth selection \\(deg \\(1,0\\), multistart 2/2, elapsed [0-9]+\\.[0-9]s, 100\\.0%, eta 0\\.0s\\)$", lines)))
})
