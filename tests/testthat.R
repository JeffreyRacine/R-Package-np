library(testthat)
library(npRmpi)

np_check_filter <- Sys.getenv("NP_CHECK_FILTER", unset = "check-minimal")
np_check_full <- identical(Sys.getenv("NP_CHECK_FULL", unset = ""), "1")

local({
  on.exit({
    try(npRmpi.quit(force = TRUE), silent = TRUE)
    try(mpi.finalize(), silent = TRUE)
  }, add = TRUE)

  if (isTRUE(np_check_full) || identical(np_check_filter, "")) {
    test_check("npRmpi")
  } else {
    test_check("npRmpi", filter = np_check_filter)
  }
})
