test_that("native progress distinguishes rendering errors from interrupts", {
  signal <- function() .Call("C_np_progress_signal",
    "bandwidth_activity_step", "bandwidth", 0L, 0L, PACKAGE = "np")
  testthat::local_mocked_bindings(
    .np_progress_signal_from_c = function(...) stop("renderer failure"),
    .package = "np")
  expect_null(signal())

  interrupted <- structure(list(message = ""), class = c("interrupt", "condition"))
  testthat::local_mocked_bindings(
    .np_progress_signal_from_c = function(...) signalCondition(interrupted),
    .package = "np")
  caught <- tryCatch(signal(), interrupt = identity)
  expect_s3_class(caught, "interrupt")
  expect_false(inherits(caught, "error"))
})

test_that("native progress uses one typed condition boundary", {
  path <- testthat::test_path("..", "..", "src", "np.c")
  if (!file.exists(path))
    skip("Native source is not present in this installed check")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_false(grepl("R_tryEval", source, fixed = TRUE))
  expect_match(source, "R_tryCatchError(np_progress_evaluate_body", fixed = TRUE)
  expect_match(source, "status = CRS_NOMAD_OBSERVER_OUTCOME_INTERRUPT",
               fixed = TRUE)
  expect_match(source, "goto observer_done;", fixed = TRUE)
})
