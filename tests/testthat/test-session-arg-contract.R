test_that("npRmpi.init validates scalar control arguments before MPI calls", {
  expect_error(npRmpi.init(nslaves = -1, mode = "spawn"), "'nslaves' must be a non-negative integer")
  expect_error(npRmpi.init(nslaves = 0, mode = "spawn"), "'nslaves' must be >= 1 for npRmpi; use package 'np' for serial workflows")
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

test_that("session receive-timeout resolver honors option/env precedence", {
  recv_timeout <- getFromNamespace(".npRmpi_session_recv_timeout", "npRmpi")

  withr::local_options(npRmpi.session.recv.timeout = 12)
  withr::local_envvar(NP_RMPI_SESSION_RECV_TIMEOUT_SEC = "999")
  expect_equal(recv_timeout(), 12)

  withr::local_options(npRmpi.session.recv.timeout = NULL)
  withr::local_envvar(NP_RMPI_SESSION_RECV_TIMEOUT_SEC = "7.5")
  expect_equal(recv_timeout(), 7.5)

  withr::local_options(npRmpi.session.recv.timeout = NULL)
  withr::local_envvar(NP_RMPI_SESSION_RECV_TIMEOUT_SEC = "bad")
  expect_equal(recv_timeout(), 0)
})
