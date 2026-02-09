test_that("npRmpi.start() is idempotent", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  # Do not require silence here: some MPI stacks print spawn/status messages.
  expect_true(isTRUE(try({ npRmpi.start(nslaves=1); TRUE }, silent=TRUE)))
  expect_true(mpi.comm.size(1) > 1)
})

test_that("npRmpi.stop() respects reuse/force", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old <- getOption("npRmpi.reuse.slaves", FALSE)
  on.exit(options(npRmpi.reuse.slaves=old), add=TRUE)

  options(npRmpi.reuse.slaves=TRUE)
  expect_true(isTRUE(try({ npRmpi.stop(force=FALSE); TRUE }, silent=TRUE)))
  expect_true(mpi.comm.size(1) > 1) # soft-close keeps pool alive

  expect_true(isTRUE(try({ npRmpi.stop(force=TRUE); TRUE }, silent=TRUE)))
  sz <- try(mpi.comm.size(1), silent=TRUE)
  if (!inherits(sz, "try-error"))
    expect_true(sz < 2)
})

test_that("npRmpi.session.info() returns useful fields", {
  # Can run without slaves; should not error.
  info <- npRmpi.session.info()
  expect_true(is.list(info))
  expect_true(all(c("npRmpi", "Rmpi", "reuse_slaves", "comm", "comm_size", "nslaves") %in% names(info)))
})
