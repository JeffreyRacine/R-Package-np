test_that("heavy npplreghat apply fanout uses master plus every active worker", {
  fanout <- getFromNamespace(".npRmpi_npplreghat_apply_fanout", "npRmpi")

  captured <- new.env(parent = emptyenv())
  fake_hat <- function(bws, txdat, exdat = txdat, y = NULL, output = "apply", ...) {
    if (identical(output, "matrix"))
      return(matrix(0, nrow = nrow(exdat), ncol = nrow(txdat)))
    if (is.null(y))
      return(rep.int(0, nrow(exdat)))
    yy <- as.matrix(y)
    matrix(0, nrow = nrow(exdat), ncol = ncol(yy))
  }

  local_mocked_bindings(
    mpi.comm.rank = function(comm = 1L) 0L,
    .npRmpi_bootstrap_worker_count = function(comm = 1L) 3L,
    .npRmpi_autodispatch_called_from_bcast = function() FALSE,
    npreghat = fake_hat,
    .npRmpi_bootstrap_run_fanout = function(tasks,
                                            worker,
                                            ncol.out,
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

  bws <- list(bw = list(yzbw = list(), list()))
  out <- fanout(
    bws = bws,
    txdat = data.frame(x = seq_len(10)),
    tzdat = data.frame(z = seq_len(10)),
    exdat = data.frame(x = seq_len(8)),
    ezdat = data.frame(z = seq_len(8)),
    y = matrix(seq_len(10), ncol = 1L),
    x.train.num = matrix(seq_len(10), ncol = 1L),
    x.eval.num = matrix(seq_len(8), ncol = 1L)
  )

  expect_true(is.matrix(out))
  expect_true(isTRUE(captured$master_local_chunk))
  expect_length(captured$tasks, 4L)
  expect_true(all(vapply(captured$tasks, function(tt) as.integer(tt$bsz), integer(1L)) > 0L))
  expect_equal(sort(unlist(lapply(captured$tasks, `[[`, "rows"), use.names = FALSE)),
               seq_len(8))
  expect_identical(captured$required.bindings$p, 1L)
})

test_that("npplreghat apply fanout is gated away from small work", {
  fanout <- getFromNamespace(".npRmpi_npplreghat_apply_fanout", "npRmpi")

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
    bws = list(bw = list(yzbw = list(), list())),
    txdat = data.frame(x = seq_len(10)),
    tzdat = data.frame(z = seq_len(10)),
    exdat = data.frame(x = seq_len(8)),
    ezdat = data.frame(z = seq_len(8)),
    y = matrix(seq_len(10), ncol = 1L),
    x.train.num = matrix(seq_len(10), ncol = 1L),
    x.eval.num = matrix(seq_len(8), ncol = 1L)
  ))
})
