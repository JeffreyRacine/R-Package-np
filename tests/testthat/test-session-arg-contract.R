test_that("npRmpi.init validates scalar control arguments before MPI calls", {
  old_allow <- getOption("npRmpi.allow.attached.Rmpi")
  options(npRmpi.allow.attached.Rmpi = TRUE)
  on.exit(options(npRmpi.allow.attached.Rmpi = old_allow), add = TRUE)

  expect_error(npRmpi.init(nslaves = 0, mode = "spawn"), "'nslaves' must be a positive integer")
  expect_error(npRmpi.init(comm = 0, mode = "spawn"), "'comm' must be a positive integer")
  expect_error(npRmpi.init(autodispatch = NA, mode = "spawn"), "'autodispatch' must be TRUE or FALSE")
  expect_error(npRmpi.init(autodispatch.verify.options = c(TRUE, FALSE), mode = "spawn"),
               "'autodispatch.verify.options' must be TRUE or FALSE")
  expect_error(npRmpi.init(np.messages = NA, mode = "spawn"), "'np.messages' must be TRUE or FALSE")
  expect_error(npRmpi.init(nonblock = NA, mode = "spawn"), "'nonblock' must be TRUE or FALSE")
  expect_error(npRmpi.init(quiet = NA, mode = "spawn"), "'quiet' must be TRUE or FALSE")
  expect_error(npRmpi.init(sleep = 0, mode = "spawn"), "'sleep' must be a positive finite numeric scalar")
})

test_that("npRmpi.quit validates scalar control arguments before MPI calls", {
  expect_error(npRmpi.quit(force = NA), "'force' must be TRUE or FALSE")
  expect_error(npRmpi.quit(dellog = NA), "'dellog' must be TRUE or FALSE")
  expect_error(npRmpi.quit(comm = 0), "'comm' must be a positive integer")
})
