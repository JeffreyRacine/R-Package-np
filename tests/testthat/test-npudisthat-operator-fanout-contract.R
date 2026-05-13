test_that("heavy npudisthat apply fanout uses master plus every active worker", {
  fanout <- getFromNamespace(".npRmpi_hat_operator_row_fanout", "npRmpi")

  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    mpi.comm.rank = function(comm = 1L) 0L,
    .npRmpi_bootstrap_worker_count = function(comm = 1L) 3L,
    .npRmpi_autodispatch_called_from_bcast = function() FALSE,
    .npRmpi_bootstrap_run_fanout = function(tasks,
                                            worker,
                                            ncol.out,
                                            what = "bootstrap",
                                            master_local_chunk = FALSE,
                                            required.bindings = NULL,
                                            ...) {
      captured$tasks <- tasks
      captured$master_local_chunk <- master_local_chunk
      captured$required.bindings <- required.bindings
      matrix(0, nrow = sum(vapply(tasks, function(tt) as.integer(tt$bsz), integer(1L))),
             ncol = as.integer(ncol.out))
    },
    .package = "npRmpi"
  )

  withr::local_options(
    npRmpi.hat.operator.fanout = TRUE,
    npRmpi.hat.operator.fanout.min.work = 1,
    npRmpi.mpi.initialized = TRUE,
    npRmpi.autodispatch.context = FALSE,
    npRmpi.autodispatch.disable = FALSE
  )

  out <- fanout(
    fun.name = "npudisthat",
    bws = list(type = "fixed"),
    tdat = data.frame(x = seq_len(10)),
    edat = data.frame(x = seq_len(8)),
    y = matrix(seq_len(10), ncol = 1L),
    output = "apply",
    target.dist = seq_len(8)
  )

  expect_true(is.matrix(out))
  expect_true(isTRUE(captured$master_local_chunk))
  expect_length(captured$tasks, 4L)
  expect_true(all(vapply(captured$tasks, function(tt) as.integer(tt$bsz), integer(1L)) > 0L))
  expect_equal(sort(unlist(lapply(captured$tasks, `[[`, "rows"), use.names = FALSE)),
               seq_len(8))
  expect_identical(captured$required.bindings$target.dist, seq_len(8))
})

test_that("heavy npudenshat apply fanout uses master plus every active worker", {
  fanout <- getFromNamespace(".npRmpi_hat_operator_row_fanout", "npRmpi")

  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    mpi.comm.rank = function(comm = 1L) 0L,
    .npRmpi_bootstrap_worker_count = function(comm = 1L) 3L,
    .npRmpi_autodispatch_called_from_bcast = function() FALSE,
    .npRmpi_bootstrap_run_fanout = function(tasks,
                                            worker,
                                            ncol.out,
                                            what = "bootstrap",
                                            master_local_chunk = FALSE,
                                            required.bindings = NULL,
                                            ...) {
      captured$tasks <- tasks
      captured$master_local_chunk <- master_local_chunk
      captured$required.bindings <- required.bindings
      matrix(0, nrow = sum(vapply(tasks, function(tt) as.integer(tt$bsz), integer(1L))),
             ncol = as.integer(ncol.out))
    },
    .package = "npRmpi"
  )

  withr::local_options(
    npRmpi.hat.operator.fanout = TRUE,
    npRmpi.hat.operator.fanout.min.work = 1,
    npRmpi.mpi.initialized = TRUE,
    npRmpi.autodispatch.context = FALSE,
    npRmpi.autodispatch.disable = FALSE
  )

  out <- fanout(
    fun.name = "npudenshat",
    bws = list(type = "fixed"),
    tdat = data.frame(x = seq_len(10)),
    edat = data.frame(x = seq_len(8)),
    y = matrix(seq_len(10), ncol = 1L),
    output = "apply"
  )

  expect_true(is.matrix(out))
  expect_true(isTRUE(captured$master_local_chunk))
  expect_length(captured$tasks, 4L)
  expect_true(all(vapply(captured$tasks, function(tt) as.integer(tt$bsz), integer(1L)) > 0L))
  expect_equal(sort(unlist(lapply(captured$tasks, `[[`, "rows"), use.names = FALSE)),
               seq_len(8))
  expect_true(is.function(captured$required.bindings$.np_udenshat_local_chunk))
})

test_that("unconditional hat apply fanout is gated away from small or unsupported calls", {
  fanout <- getFromNamespace(".npRmpi_hat_operator_row_fanout", "npRmpi")

  local_mocked_bindings(
    mpi.comm.rank = function(comm = 1L) 0L,
    .npRmpi_bootstrap_worker_count = function(comm = 1L) 2L,
    .npRmpi_autodispatch_called_from_bcast = function() FALSE,
    .npRmpi_bootstrap_run_fanout = function(...) stop("fanout should not be called"),
    .package = "npRmpi"
  )

  withr::local_options(
    npRmpi.hat.operator.fanout = TRUE,
    npRmpi.hat.operator.fanout.min.work = 1e6,
    npRmpi.mpi.initialized = TRUE,
    npRmpi.autodispatch.context = FALSE,
    npRmpi.autodispatch.disable = FALSE
  )

  expect_null(fanout(
    fun.name = "npudisthat",
    bws = list(type = "fixed"),
    tdat = data.frame(x = seq_len(10)),
    edat = data.frame(x = seq_len(8)),
    y = matrix(seq_len(10), ncol = 1L),
    output = "apply",
    target.dist = seq_len(8)
  ))

  withr::local_options(npRmpi.hat.operator.fanout.min.work = 1)
  expect_null(fanout(
    fun.name = "npudisthat",
    bws = list(type = "fixed"),
    tdat = data.frame(x = seq_len(10)),
    edat = data.frame(x = seq_len(8)),
    y = matrix(seq_len(10), ncol = 1L),
    output = "matrix",
    target.dist = seq_len(8)
  ))
})
