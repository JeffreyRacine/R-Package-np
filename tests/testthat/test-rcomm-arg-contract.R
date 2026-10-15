test_that("MPI info strings retain values, missing keys and handle reuse", {
  skip_if_not(isTRUE(getOption("npRmpi.mpi.initialized", FALSE)),
              "MPI initialization is required for info handles")
  check_info <- function() {
    .Call("mpi_info_create", 0L, PACKAGE = "npRmpi")
    on.exit(.Call("mpi_info_free", 0L, PACKAGE = "npRmpi"))
    .Call("mpi_info_set", 0L, "np_test_key", "abcdef", PACKAGE = "npRmpi")
    for (capacity in c(6L, 7L, 32L)) {
      expect_identical(.Call("mpi_info_get", 0L, "np_test_key", capacity,
                            PACKAGE = "npRmpi"), "abcdef")
      expect_identical(.Call("mpi_info_get", 0L, "np_test_missing", capacity,
                            PACKAGE = "npRmpi"), "")
    }
    .Call("mpi_info_set", 0L, "np_test_key", "xy", PACKAGE = "npRmpi")
    expect_identical(.Call("mpi_info_get", 0L, "np_test_key", 2L,
                          PACKAGE = "npRmpi"), "xy")
    for (capacity in c(-1L, NA_integer_))
      expect_error(.Call("mpi_info_get", 0L, "np_test_key", capacity,
                         PACKAGE = "npRmpi"), "valuelen must be non-negative")
  }
  check_info()
  check_info()
})

test_that("mpi.comm.spawn validates scalar positive integer nslaves", {
  skip_if_not(is.loaded("mpi_comm_spawn"), "MPI_Comm_spawn is not supported in this build")
  spawn_fun <- get("mpi.comm.spawn", envir = asNamespace("npRmpi"), inherits = FALSE)
  slave_path <- system.file("slavedaemon.R", package = "npRmpi")

  expect_error(
    spawn_fun(slave = slave_path, nslaves = 0),
    "Choose a positive number of slaves."
  )
  expect_error(
    spawn_fun(slave = slave_path, nslaves = NA),
    "Choose a positive number of slaves."
  )
  expect_error(
    spawn_fun(slave = slave_path, nslaves = c(1, 2)),
    "Choose a positive number of slaves."
  )
  expect_error(
    spawn_fun(slave = slave_path, nslaves = 1.5),
    "Choose a positive number of slaves."
  )
})

test_that("low-level Rmpi wrappers reject out-of-range native handle indices", {
  expect_error(.Call("mpi_info_create", -1L, PACKAGE = "npRmpi"), "info=-1 is out of range")
  expect_error(.Call("mpi_info_set", -1L, "k", "v", PACKAGE = "npRmpi"), "info=-1 is out of range")
  expect_error(.Call("mpi_info_get", -1L, "k", 1L, PACKAGE = "npRmpi"), "info=-1 is out of range")
  expect_error(.Call("mpi_info_free", -1L, PACKAGE = "npRmpi"), "info=-1 is out of range")

  expect_error(.Call("mpi_gather", 1L, 1L, integer(1), 0L, -1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_gatherv", 1L, 1L, integer(1), 1L, 0L, -1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_scatter", 1L, 1L, integer(1), 0L, -1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_scatterv", 1L, 1L, 1L, integer(1), 0L, -1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_allgather", 1L, 1L, integer(1), -1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_allgatherv", 1L, 1L, integer(1), 1L, -1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_bcast", integer(1), 1L, 0L, -1L, 1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_reduce", 1L, 1L, 1L, 0L, -1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_allreduce", 1L, 1L, 1L, -1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")

  expect_error(.Call("mpi_cart_create", -1L, 1L, 0L, 0L, 3L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_cart_create", 0L, 1L, 0L, 0L, -1L, PACKAGE = "npRmpi"), "cartesian communicator=-1 is out of range")
  expect_error(.Call("mpi_cartdim_get", -1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_cart_get", -1L, 1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_cart_rank", -1L, 0L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_cart_coords", -1L, 0L, 1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")
  expect_error(.Call("mpi_cart_shift", -1L, 0L, 1L, PACKAGE = "npRmpi"), "communicator=-1 is out of range")

  skip_if_not(is.loaded("mpi_comm_spawn"), "MPI_Comm_spawn is not supported in this build")
  expect_error(.Call("mpi_comm_spawn", "", character(0), 1L, -1L, 0L, 0L, TRUE, PACKAGE = "npRmpi"), "info=-1 is out of range")
  expect_error(.Call("mpi_comm_spawn", "", character(0), 1L, 0L, 0L, -1L, TRUE, PACKAGE = "npRmpi"), "intercommunicator=-1 is out of range")
})
