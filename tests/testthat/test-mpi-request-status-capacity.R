test_that("growing requests reserves enough statuses for every bulk completion", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  original.request <- mpi.request.maxsize()
  original.status <- mpi.status.maxsize()
  wanted <- max(original.request, original.status) + 1L
  mpi.realloc.request(wanted)
  expect_gte(mpi.request.maxsize(), wanted)
  expect_gte(mpi.status.maxsize(), mpi.request.maxsize())
  ## New slots are MPI_REQUEST_NULL, so no peer or blocking send is needed.
  expect_null(mpi.waitall(wanted))
  expect_true(mpi.testall(wanted))
  for (complete in list(mpi.waitsome, mpi.testsome)) {
    result <- complete(wanted)
    expect_equal(result$count, .Call("mpi_undefined", PACKAGE = "npRmpi"))
    expect_null(result$indices)
  }
  mpi.realloc.status(wanted + 3L)
  mpi.realloc.request(wanted + 2L)
  expect_identical(mpi.status.maxsize(), wanted + 3L)
  expect_identical(mpi.request.maxsize(), wanted + 2L)
  mpi.realloc.request(wanted)
  expect_identical(mpi.request.maxsize(), wanted + 2L)
  expect_identical(mpi.status.maxsize(), wanted + 3L)
})
