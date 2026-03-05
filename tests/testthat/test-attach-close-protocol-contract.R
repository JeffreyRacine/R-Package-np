test_that("attach close helper primitives are stable", {
  tag_fn <- getFromNamespace(".npRmpi_attach_close_ack_tag", "npRmpi")
  timeout_fn <- getFromNamespace(".npRmpi_attach_close_ack_timeout", "npRmpi")
  next_sid <- getFromNamespace(".npRmpi_attach_next_session_id", "npRmpi")

  expect_identical(tag_fn(1L), 62001L)
  expect_identical(tag_fn(3L), 62003L)

  withr::local_options(list(npRmpi.attach.close.ack.timeout = 2.5))
  expect_equal(timeout_fn(), 2.5)

  withr::local_options(list(npRmpi.attach.close.ack.timeout = NULL))
  withr::local_envvar(c(NP_RMPI_ATTACH_CLOSE_ACK_TIMEOUT_SEC = "7"))
  expect_equal(timeout_fn(), 7)

  withr::local_options(list(npRmpi.attach.session.counter = 41L))
  expect_identical(next_sid(), 42L)
  expect_identical(next_sid(), 43L)
})

test_that("attach worker loop performs ACK send before worker quit", {
  fn <- getFromNamespace(".npRmpi_session_attach_worker_loop", "npRmpi")
  txt <- paste(deparse(body(fn), width.cutoff = 500L), collapse = " ")
  expect_match(txt, "\\.npRmpi_attach_send_close_ack\\(")
  expect_match(txt, "mpi\\.quit\\(")
})

test_that("worker loop routes shutdown detection through helper", {
  fn <- getFromNamespace(".npRmpi_worker_loop", "npRmpi")
  txt <- paste(deparse(body(fn), width.cutoff = 500L), collapse = " ")
  expect_match(txt, "\\.npRmpi_worker_is_break_message\\(")
})

test_that("attach quit path contains idempotent state guard and ACK collection", {
  fn <- getFromNamespace("npRmpi.quit", "npRmpi")
  txt <- paste(deparse(body(fn), width.cutoff = 500L), collapse = " ")
  expect_match(txt, "npRmpi\\.attach\\.close\\.state")
  expect_match(txt, "\\.npRmpi_attach_collect_close_acks\\(")
  expect_match(txt, "\\.npRmpi_attach_state_reset\\(")
})

test_that("attach ACK collector marks complete ACK sets as OK", {
  collect <- getFromNamespace(".npRmpi_attach_collect_close_acks", "npRmpi")

  testthat::local_mocked_bindings(
    mpi.iprobe = function(source, tag, comm = 1L, status = 0L) TRUE,
    mpi.recv.Robj = function(source, tag, comm = 1L, status = 0L) {
      list(kind = "attach_close_ack", rank = as.integer(source), session_id = 11L)
    },
    .env = asNamespace("npRmpi")
  )

  out <- collect(comm = 1L, session_id = 11L, expected = c(1L, 2L), timeout_sec = 0.05, poll_sleep = 0.001)
  expect_true(isTRUE(out$ok))
  expect_identical(out$acked, c(1L, 2L))
  expect_identical(out$missing, integer(0))
})

test_that("attach ACK collector reports missing ranks on partial ACK sets", {
  collect <- getFromNamespace(".npRmpi_attach_collect_close_acks", "npRmpi")

  testthat::local_mocked_bindings(
    mpi.iprobe = function(source, tag, comm = 1L, status = 0L) identical(as.integer(source), 1L),
    mpi.recv.Robj = function(source, tag, comm = 1L, status = 0L) {
      list(kind = "attach_close_ack", rank = as.integer(source), session_id = 12L)
    },
    .env = asNamespace("npRmpi")
  )

  out <- collect(comm = 1L, session_id = 12L, expected = c(1L, 2L), timeout_sec = 0.05, poll_sleep = 0.001)
  expect_false(isTRUE(out$ok))
  expect_identical(out$acked, 1L)
  expect_identical(out$missing, 2L)
})
