test_that("density MADS object progress starts only the missing lifecycle", {
  calls <- 0L
  active <- FALSE
  local_mocked_bindings(
    .np_progress_bandwidth_active = function() active,
    .np_progress_select_bandwidth_enhanced = function(label, expr) {
      calls <<- calls + 1L
      force(expr)
    }, .package = "np")
  run <- function(compute = TRUE, solver = "mads", evaluate = FALSE) {
    .np_density_bw_progress(compute, solver, evaluate, 17L)
  }
  expect_identical(run(), 17L)
  expect_identical(run(solver = "mads+powell"), 17L)
  expect_identical(calls, 2L)
  for (solver in list("powell", c("powell", "mads", "mads+powell"), NA_character_, 1))
    expect_identical(run(solver = solver), 17L)
  expect_identical(run(compute = FALSE), 17L)
  expect_identical(run(evaluate = TRUE), 17L)
  active <- TRUE
  expect_identical(run(), 17L)
  expect_identical(calls, 2L)
  active <- FALSE
  frame <- function(x) .np_density_bw_progress(TRUE, "mads", FALSE,
    list(call = match.call(), missing = missing(x)))
  expect_identical(frame()$call, quote(frame()))
  expect_true(frame()$missing)
  expect_error(.np_density_bw_progress(TRUE, "mads", FALSE,
    stop("deliberate progress error")), "deliberate progress error", fixed = TRUE)
})

test_that("density progress restores the existing lifecycle on error", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  was.active <- .np_progress_bandwidth_active()
  expect_error(.np_density_bw_progress(TRUE, "mads", FALSE,
    stop("density entry failure")), "density entry failure", fixed = TRUE)
  expect_identical(.np_progress_bandwidth_active(), was.active)
  expect_false(getOption("np.messages"))
})

