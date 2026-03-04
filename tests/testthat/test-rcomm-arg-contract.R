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
