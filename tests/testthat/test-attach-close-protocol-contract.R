test_that("attach close helper primitives are stable", {
  tag_fn <- getFromNamespace(".npRmpi_attach_close_ack_tag", "npRmpi")
  release_tag_fn <- getFromNamespace(".npRmpi_attach_close_release_tag", "npRmpi")
  timeout_fn <- getFromNamespace(".npRmpi_attach_close_ack_timeout", "npRmpi")
  next_sid <- getFromNamespace(".npRmpi_attach_next_session_id", "npRmpi")

  expect_identical(tag_fn(1L), 62001L)
  expect_identical(tag_fn(3L), 62003L)
  expect_identical(release_tag_fn(1L), 62101L)
  expect_identical(release_tag_fn(3L), 62103L)

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
  expect_match(txt, "\\.npRmpi_attach_wait_close_release\\(")
  expect_match(txt, "mpi\\.barrier\\(0\\)")
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
  expect_match(txt, "\\.npRmpi_attach_send_close_release\\(")
  expect_match(txt, "isTRUE\\(ack\\.info\\$ok\\) && isTRUE\\(release\\.info\\$ok\\)")
  expect_match(txt, "mpi\\.barrier\\(0\\)")
  expect_match(txt, "\\.npRmpi_attach_state_reset\\(")
})

test_that("attach ACK collector marks complete ACK sets as OK", {
  collect <- getFromNamespace(".npRmpi_attach_collect_close_acks", "npRmpi")

  testthat::with_mocked_bindings({
    out <- collect(comm = 1L, session_id = 11L, expected = c(1L, 2L), timeout_sec = 0.05, poll_sleep = 0.001)
    expect_true(isTRUE(out$ok))
    expect_identical(out$acked, c(1L, 2L))
    expect_identical(out$missing, integer(0))
  },
    mpi.iprobe = function(source, tag, comm = 1L, status = 0L) TRUE,
    mpi.recv.Robj = function(source, tag, comm = 1L, status = 0L) {
      list(kind = "attach_close_ack", rank = as.integer(source), session_id = 11L)
    },
    .package = "npRmpi"
  )
})

test_that("attach ACK collector reports missing ranks on partial ACK sets", {
  collect <- getFromNamespace(".npRmpi_attach_collect_close_acks", "npRmpi")

  testthat::with_mocked_bindings({
    out <- collect(comm = 1L, session_id = 12L, expected = c(1L, 2L), timeout_sec = 0.05, poll_sleep = 0.001)
    expect_false(isTRUE(out$ok))
    expect_identical(out$acked, 1L)
    expect_identical(out$missing, 2L)
  },
    mpi.iprobe = function(source, tag, comm = 1L, status = 0L) identical(as.integer(source), 1L),
    mpi.recv.Robj = function(source, tag, comm = 1L, status = 0L) {
      list(kind = "attach_close_ack", rank = as.integer(source), session_id = 12L)
    },
    .package = "npRmpi"
  )
})

test_that("attach release sender reports per-rank send failures", {
  send_release <- getFromNamespace(".npRmpi_attach_send_close_release", "npRmpi")

  testthat::with_mocked_bindings({
    out <- send_release(comm = 1L, ranks = c(1L, 2L), session_id = 7L, barrier = TRUE)
    expect_false(isTRUE(out$ok))
    expect_identical(out$sent, 1L)
    expect_identical(out$failed, 2L)
  },
    mpi.send.Robj = function(data, dest, tag, comm = 1L) {
      if (identical(as.integer(dest), 2L))
        stop("send failed")
      invisible(TRUE)
    },
    .package = "npRmpi"
  )
})

test_that("attach release waiter accepts matching session payload", {
  wait_release <- getFromNamespace(".npRmpi_attach_wait_close_release", "npRmpi")

  testthat::with_mocked_bindings({
    out <- wait_release(comm = 1L, session_id = 9L, timeout_sec = 0.05, poll_sleep = 0.001)
    expect_true(isTRUE(out$ok))
    expect_true(isTRUE(out$barrier))
    expect_identical(out$payload$session_id, 9L)
  },
    mpi.comm.rank = function(comm = 1L) 1L,
    mpi.iprobe = function(source, tag, comm = 1L, status = 0L) TRUE,
    mpi.recv.Robj = function(source, tag, comm = 1L, status = 0L) {
      list(kind = "attach_close_release", session_id = 9L, barrier = TRUE)
    },
    .package = "npRmpi"
  )
})

test_that("attach release waiter rejects mismatched session payload", {
  wait_release <- getFromNamespace(".npRmpi_attach_wait_close_release", "npRmpi")

  testthat::with_mocked_bindings({
    out <- wait_release(comm = 1L, session_id = 9L, timeout_sec = 0.05, poll_sleep = 0.001)
    expect_false(isTRUE(out$ok))
    expect_false(isTRUE(out$barrier))
  },
    mpi.comm.rank = function(comm = 1L) 1L,
    mpi.iprobe = function(source, tag, comm = 1L, status = 0L) TRUE,
    mpi.recv.Robj = function(source, tag, comm = 1L, status = 0L) {
      list(kind = "attach_close_release", session_id = 10L, barrier = TRUE)
    },
    .package = "npRmpi"
  )
})
