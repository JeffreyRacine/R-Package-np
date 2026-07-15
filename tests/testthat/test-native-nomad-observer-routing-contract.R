test_that("all native NOMAD consumers use the observed solve wrapper", {
  source_file <- file.path(npRmpi_namespace_hygiene_root(), "src", "np.c")
  source_text <- paste(readLines(source_file, warn = FALSE), collapse = "\n")

  expect_false(grepl(
    "R_GetCCallable\\(\\\"crs\\\", \\\"crs_nomad_solve\\\"\\)",
    source_text,
    perl = TRUE
  ))
  calls <- gregexpr(
    "np_nomad_solve_with_progress\\s*\\(\\s*&problem",
    source_text,
    perl = TRUE
  )[[1L]]
  expect_identical(length(calls), 7L)

  expect_match(source_text, "crs_nomad_solve_observed", fixed = TRUE)
  expect_match(source_text, "crs_nomad_observer_poll", fixed = TRUE)
  expect_match(source_text, "if (my_rank != 0)", fixed = TRUE)
})

test_that("native NOMAD observer config follows the shared renderer interval", {
  config <- getFromNamespace(
    ".np_progress_nomad_native_observer_config",
    "npRmpi"
  )
  old_opts <- options(np.progress.interval.unknown.sec = 0.375)
  on.exit(options(old_opts), add = TRUE)

  disabled <- config()
  expect_identical(unname(disabled[["enabled"]]), 0)
  expect_equal(unname(disabled[["interval"]]), 0.375)
})

test_that("native NOMAD observer dispatch classifies ordinary errors and interrupts", {
  dispatch <- getFromNamespace(
    ".np_progress_nomad_native_observer_dispatch",
    "npRmpi"
  )

  ordinary <- with_npRmpi_degree_bindings(
    list(.np_progress_bandwidth_activity_step = function(...) stop("observer boom")),
    dispatch(1L, integer(), integer(), NA_real_)
  )
  expect_identical(ordinary[[1L]], 1L)
  expect_match(ordinary[[2L]], "observer boom")

  interrupted <- with_npRmpi_degree_bindings(
    list(.np_progress_bandwidth_activity_step = function(...) {
      stop(structure(
        list(message = "observer interrupt"),
        class = c("interrupt", "condition")
      ))
    }),
    dispatch(1L, integer(), integer(), NA_real_)
  )
  expect_identical(interrupted[[1L]], 2L)
  expect_match(interrupted[[2L]], "observer interrupt")
})

test_that("observer failure reporting cannot turn a completed solve into an error", {
  report <- getFromNamespace(
    ".np_progress_nomad_native_observer_report",
    "npRmpi"
  )
  old_opts <- options(warn = 2)
  on.exit(options(old_opts), add = TRUE)

  expect_message(
    expect_invisible(report("observer boom")),
    "NOMAD progress observer disabled: observer boom"
  )
})
