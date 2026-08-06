test_that("MPI objective collectives retain progress while rank zero waits", {
  source_file <- npRmpi_test_source_path("src", "jksum.c")
  skip_if(is.null(source_file), "source C files unavailable in installed test context")
  src <- paste(readLines(source_file, warn = FALSE), collapse = "\n")

  expect_match(src, "#if MPI_VERSION >= 3", fixed = TRUE)
  expect_match(src, "np_mpi_bandwidth_progress_wait", fixed = TRUE)
  expect_match(src, "np_progress_bandwidth_is_active()", fixed = TRUE)
  expect_match(src, "MPI_Iallgather", fixed = TRUE)
  expect_match(src, "MPI_Iallgatherv", fixed = TRUE)
  expect_match(src, "MPI_Iallreduce", fixed = TRUE)
  expect_match(src, "np_mpi_allgather_in_place_double", fixed = TRUE)
  expect_match(src, "np_mpi_allgatherv_in_place_double", fixed = TRUE)
  expect_match(src, "np_mpi_allreduce_in_place_double", fixed = TRUE)

  ## Fit/evaluation calls and MPI-2 builds retain the incumbent blocking
  ## collectives; only an active bandwidth objective uses the polling wait.
  expect_match(src, "MPI_Allgather(MPI_IN_PLACE", fixed = TRUE)
  expect_match(src, "MPI_Allgatherv(MPI_IN_PLACE", fixed = TRUE)
  expect_match(src, "MPI_Allreduce(MPI_IN_PLACE", fixed = TRUE)
})
