# This subprocess protects broadcast-replica observer scope, which the local
# observer fixture cannot exercise. Keep this MPI lifecycle gate off CRAN.
test_that("pooled native observers retain broadcast handling", {
  skip_on_cran()
  skip_if_not_installed("crs")
  fixture <- normalizePath(test_path("fixtures", "native-nomad-pooled-observer.R"))
  helpers <- normalizePath(test_path("helper-mpi.R"))
  lines <- c(
    "library(npRmpi)",
    paste0("source(", deparse(helpers), ")"),
    paste0("source(", deparse(fixture), ")")
  )
  child <- npRmpi_run_isolated_contract(lines,
    marker = "NATIVE_POOLED_OBSERVER_PASS", timeout = 45L)
  skip_if(is.null(child), "MPI subprocess environment unavailable")
  expect_identical(child$status, 0L, info = paste(child$output, collapse = "\n"))
  expect_true(child$witnessed, info = paste(child$output, collapse = "\n"))
})
