# Elapsed time belongs to the specific invocation. Validate it before
# comparing the reproducible fields of two public specification-test results.
without_test_elapsed_time <- function(result) {
  elapsed <- result[["total.time"]]
  expect_type(elapsed, "double")
  expect_length(elapsed, 1L)
  expect_true(is.finite(elapsed) && elapsed >= 0)
  result[["total.time"]] <- NULL
  result
}
