npRmpi_expect_warning_messages <- function(code, expected, fixed = TRUE) {
  expected <- as.character(expected)
  observed <- character()
  value <- withCallingHandlers(
    force(code),
    warning = function(condition) {
      observed <<- c(observed, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(length(observed), length(expected))
  for (index in seq_len(min(length(observed), length(expected)))) {
    expect_match(observed[[index]], expected[[index]], fixed = fixed)
  }

  value
}
