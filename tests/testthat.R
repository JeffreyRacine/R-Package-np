library(testthat)
library(np)

np_check_filter <- Sys.getenv("NP_CHECK_FILTER", unset = "check-core-smoke")
np_check_full <- identical(Sys.getenv("NP_CHECK_FULL", unset = ""), "1")

if (isTRUE(np_check_full) || identical(np_check_filter, "")) {
  test_check("np")
} else {
  test_check("np", filter = np_check_filter)
}
