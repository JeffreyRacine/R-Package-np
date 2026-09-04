test_that("invalid communicator handles fail fast with readable errors", {
  expect_error(npRmpi::mpi.comm.size(99999), "communicator=99999 is out of range", fixed = TRUE)
  expect_error(npRmpi::mpi.comm.rank(-1), "communicator=-1 is out of range", fixed = TRUE)
  expect_error(npRmpi::mpi.barrier(99999), "communicator=99999 is out of range", fixed = TRUE)
})

test_that("null communicator queries do not publish undefined scalars", {
  null_comm <- 9L

  expect_equal(npRmpi::mpi.comm.size(null_comm), 0L)
  expect_error(npRmpi::mpi.comm.rank(null_comm), "mpi_comm_rank: communicator is NULL", fixed = TRUE)
  expect_error(npRmpi:::mpi.comm.test.inter(null_comm), "NULL communicator", fixed = TRUE)

  expect_error(
    .Call("mpi_comm_size", null_comm, PACKAGE = "npRmpi"),
    "mpi_comm_size: communicator is NULL",
    fixed = TRUE
  )
  expect_error(
    .Call("mpi_comm_test_inter", null_comm, PACKAGE = "npRmpi"),
    "mpi_comm_test_inter: communicator is NULL",
    fixed = TRUE
  )
})

test_that("invalid status and request handles fail fast with readable errors", {
  expect_error(
    npRmpi:::mpi.get.count(type = 4, status = 99999),
    "status=99999 is out of range",
    fixed = TRUE
  )
  expect_error(
    npRmpi:::mpi.wait(99999, status = 0),
    "request=99999 is out of range",
    fixed = TRUE
  )
  expect_error(
    npRmpi:::mpi.test(99999, status = 0),
    "request=99999 is out of range",
    fixed = TRUE
  )
  expect_error(
    npRmpi:::mpi.waitall(999999),
    "request count=999999 is out of range",
    fixed = TRUE
  )
  expect_error(
    npRmpi:::mpi.testall(999999),
    "request count=999999 is out of range",
    fixed = TRUE
  )
  expect_error(
    npRmpi:::mpi.waitsome(999999),
    "request count=999999 is out of range",
    fixed = TRUE
  )
  expect_error(
    npRmpi:::mpi.testsome(999999),
    "request count=999999 is out of range",
    fixed = TRUE
  )
})

test_that("sendrecv_replace default branch still works on valid self traffic", {
  rank <- npRmpi::mpi.comm.rank(0L)
  expect_equal(
    npRmpi:::mpi.sendrecv.replace(
      3.25,
      type = 5,
      dest = rank,
      sendtag = 11L,
      source = rank,
      recvtag = 11L,
      comm = 0L,
      status = 0L
    ),
    3.25
  )
})
