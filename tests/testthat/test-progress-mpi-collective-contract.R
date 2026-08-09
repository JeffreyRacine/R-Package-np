npRmpi_collective_source_region <- function(source, start_marker, stop_marker) {
  start <- regexpr(start_marker, source, fixed = TRUE)[[1L]]
  if (start < 1L)
    stop("collective wrapper start marker is missing: ", start_marker)
  remainder <- substring(source, start)
  stop <- regexpr(stop_marker, remainder, fixed = TRUE)[[1L]]
  if (stop < 2L)
    stop("collective wrapper stop marker is missing: ", stop_marker)
  substring(remainder, 1L, stop - 1L)
}

test_that("MPI typed collective wrappers retain progress and fallbacks", {
  source_file <- npRmpi_test_source_path("src", "jksum.c")
  skip_if(is.null(source_file), "source C files unavailable in installed test context")
  src <- paste(readLines(source_file, warn = FALSE), collapse = "\n")

  progress <- npRmpi_collective_source_region(
    src,
    "#if MPI_VERSION >= 3\nstatic void np_mpi_bandwidth_progress_wait",
    "static void np_mpi_allgather_in_place_double"
  )
  allgather <- npRmpi_collective_source_region(
    src,
    "static void np_mpi_allgather_in_place_double",
    "static void np_mpi_allgatherv_in_place_double"
  )
  allgatherv <- npRmpi_collective_source_region(
    src,
    "static void np_mpi_allgatherv_in_place_double",
    "static void np_mpi_allreduce_in_place_double"
  )
  allreduce <- npRmpi_collective_source_region(
    src,
    "static void np_mpi_allreduce_in_place_double",
    "typedef struct {"
  )

  expect_match(progress, "#if MPI_VERSION >= 3", fixed = TRUE)
  expect_match(progress, "MPI_Wait(request, MPI_STATUS_IGNORE)", fixed = TRUE)
  expect_match(progress, "MPI_Test(request, &complete, MPI_STATUS_IGNORE)", fixed = TRUE)
  expect_match(progress, "np_progress_bandwidth_loop_step()", fixed = TRUE)

  wrappers <- list(
    allgather = list(body = allgather, async = "MPI_Iallgather(",
                     blocking = "MPI_Allgather(MPI_IN_PLACE"),
    allgatherv = list(body = allgatherv, async = "MPI_Iallgatherv(",
                      blocking = "MPI_Allgatherv(MPI_IN_PLACE"),
    allreduce = list(body = allreduce, async = "MPI_Iallreduce(",
                     blocking = "MPI_Allreduce(MPI_IN_PLACE")
  )

  for (name in names(wrappers)) {
    wrapper <- wrappers[[name]]
    expect_match(wrapper$body, "#if MPI_VERSION >= 3", fixed = TRUE,
                 info = name)
    expect_match(wrapper$body, "np_progress_bandwidth_is_active()",
                 fixed = TRUE, info = name)
    expect_match(wrapper$body, wrapper$async, fixed = TRUE, info = name)
    expect_match(wrapper$body, "np_mpi_bandwidth_progress_wait(&request, label)",
                 fixed = TRUE, info = name)
    expect_match(wrapper$body, wrapper$blocking, fixed = TRUE, info = name)
  }

  ## Fit/evaluation calls and MPI-2 builds retain the incumbent blocking
  ## collectives; only an active bandwidth objective uses the polling wait.
  expect_false(grepl("MPI_Iallgatherv\\(", allgather, perl = TRUE))
  expect_false(grepl("MPI_Iallreduce\\(", allgather, perl = TRUE))
  expect_false(grepl("MPI_Iallgather\\(", allgatherv, perl = TRUE))
  expect_false(grepl("MPI_Iallreduce\\(", allgatherv, perl = TRUE))
  expect_false(grepl("MPI_Iallgather\\(", allreduce, perl = TRUE))
  expect_false(grepl("MPI_Iallgatherv\\(", allreduce, perl = TRUE))
})
