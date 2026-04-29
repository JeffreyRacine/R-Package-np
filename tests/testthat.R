library(testthat)
library(npRmpi)

np_check_filter <- Sys.getenv("NP_CHECK_FILTER", unset = "check-minimal")
np_check_full <- identical(Sys.getenv("NP_CHECK_FULL", unset = ""), "1")

local({
  test_error <- NULL

  tryCatch({
    if (isTRUE(np_check_full) || identical(np_check_filter, "")) {
      test_check("npRmpi")
    } else {
      test_check("npRmpi", filter = np_check_filter)
    }
  }, error = function(e) {
    test_error <<- e
  })

  try(npRmpi.quit(force = TRUE), silent = TRUE)

  if (!is.null(test_error)) {
    message(conditionMessage(test_error))
    quit(save = "no", status = 1, runLast = FALSE)
  }
})
