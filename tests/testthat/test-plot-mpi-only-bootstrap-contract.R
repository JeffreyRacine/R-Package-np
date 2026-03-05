test_that("bootstrap fanout is rejected inside mpi.bcast.cmd context", {
  fanout_enabled <- getFromNamespace(".npRmpi_bootstrap_fanout_enabled", "npRmpi")

  mpi.bcast.cmd <- function(expr) force(expr)

  expect_error(
    mpi.bcast.cmd(
      fanout_enabled(
        comm = 1L,
        n = 10L,
        B = 10L,
        chunk.size = 2L,
        what = "unit-test"
      )
    ),
    "cannot run inside mpi\\.bcast\\.cmd context"
  )
})

test_that("wild bootstrap helper requires MPI fanout and never falls back locally", {
  wild_boot_t <- getFromNamespace(".np_wild_boot_t", "npRmpi")

  withr::local_options(npRmpi.mpi.initialized = FALSE)

  H <- diag(2)
  fit.mean <- c(0.1, -0.2)
  residuals <- c(0.3, 0.4)

  expect_error(
    wild_boot_t(H = H, fit.mean = fit.mean, residuals = residuals, B = 3L, wild = "rademacher"),
    "requires an active MPI slave pool"
  )
})

test_that("wild fanout does not switch to master-local when workers are active", {
  run_fanout <- getFromNamespace(".npRmpi_bootstrap_run_fanout", "npRmpi")
  tf <- tempfile("npRmpi-boot-transport-", fileext = ".tsv")
  old.trace <- getOption("npRmpi.bootstrap.transport.trace.file")
  options(npRmpi.bootstrap.transport.trace.file = tf)
  on.exit(options(npRmpi.bootstrap.transport.trace.file = old.trace), add = TRUE)

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_master_only_mode = function(comm = 1L) FALSE,
    .npRmpi_bootstrap_worker_count = function(comm = 1L) 2L,
    mpi.applyLB = function(X, FUN, ..., comm = 1L) {
      lapply(X, function(task) do.call(FUN, c(list(task), list(...))))
    },
    .package = "npRmpi"
  )

  tasks <- list(list(start = 1L, bsz = 2L, seed = 123L))
  worker <- function(task) {
    matrix(seq_len(as.integer(task$bsz)), nrow = as.integer(task$bsz), ncol = 1L)
  }

  out <- run_fanout(
    tasks = tasks,
    worker = worker,
    ncol.out = 1L,
    what = "wild",
    profile.where = "unit.test:wild"
  )
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(2L, 1L))

  lines <- readLines(tf, warn = FALSE)
  start.lines <- lines[grepl("event=fanout.start", lines, fixed = TRUE)]
  expect_true(length(start.lines) >= 1L)
  expect_true(any(grepl("master_local=FALSE", start.lines, fixed = TRUE)))
})

test_that("wild fanout master-assist is remote-only when workers are active", {
  run_fanout <- getFromNamespace(".npRmpi_bootstrap_run_fanout", "npRmpi")
  tf <- tempfile("npRmpi-boot-transport-", fileext = ".tsv")
  old.trace <- getOption("npRmpi.bootstrap.transport.trace.file")
  options(npRmpi.bootstrap.transport.trace.file = tf)
  on.exit(options(npRmpi.bootstrap.transport.trace.file = old.trace), add = TRUE)

  qenv <- new.env(parent = emptyenv())
  qenv$queue <- list()

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_master_only_mode = function(comm = 1L) FALSE,
    .npRmpi_bootstrap_worker_count = function(comm = 1L) 1L,
    mpi.any.source = function() 0L,
    mpi.any.tag = function() 0L,
    mpi.bcast.cmd = function(...) invisible(NULL),
    mpi.bcast.Robj = function(...) invisible(NULL),
    mpi.send.Robj = function(obj, dest, tag, comm = 1L) {
      if (is.list(obj) && !is.null(obj$data.arg) && length(obj$data.arg) == 1L) {
        task <- obj$data.arg[[1L]]
        qenv$queue[[length(qenv$queue) + 1L]] <<- list(
          src = as.integer(dest),
          tag = as.integer(tag),
          payload = matrix(
            rep.int(as.integer(task$start), as.integer(task$bsz)),
            nrow = as.integer(task$bsz),
            ncol = 1L
          )
        )
      }
      invisible(NULL)
    },
    mpi.iprobe = function(source, tag, comm = 1L) length(qenv$queue) > 0L,
    mpi.get.sourcetag = function() {
      stopifnot(length(qenv$queue) > 0L)
      c(as.integer(qenv$queue[[1L]]$src), as.integer(qenv$queue[[1L]]$tag))
    },
    mpi.recv.Robj = function(source, tag, comm = 1L) {
      stopifnot(length(qenv$queue) > 0L)
      payload <- qenv$queue[[1L]]$payload
      qenv$queue[[1L]] <<- NULL
      payload
    },
    .package = "npRmpi"
  )

  tasks <- list(
    list(start = 1L, bsz = 2L, seed = 101L),
    list(start = 3L, bsz = 2L, seed = 202L)
  )

  out <- run_fanout(
    tasks = tasks,
    worker = function(task) stop("master-local execution must not run"),
    ncol.out = 1L,
    what = "wild",
    profile.where = "unit.test:remote-only"
  )

  expect_true(is.matrix(out))
  expect_identical(dim(out), c(4L, 1L))
  expect_equal(out[, 1L], c(1, 1, 3, 3))

  lines <- readLines(tf, warn = FALSE)
  done.lines <- lines[grepl("event=fanout.master_assist.done", lines, fixed = TRUE)]
  expect_true(length(done.lines) >= 1L)
  expect_true(any(grepl("local_done=0", done.lines, fixed = TRUE)))
})

test_that("wild fanout fails fast on dispatch timeout when worker replies stall", {
  run_fanout <- getFromNamespace(".npRmpi_bootstrap_run_fanout", "npRmpi")
  withr::local_options(npRmpi.bootstrap.dispatch.timeout.sec = 0.02)

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_master_only_mode = function(comm = 1L) FALSE,
    .npRmpi_bootstrap_worker_count = function(comm = 1L) 1L,
    mpi.any.source = function() 0L,
    mpi.any.tag = function() 0L,
    mpi.bcast.cmd = function(...) invisible(NULL),
    mpi.bcast.Robj = function(...) invisible(NULL),
    mpi.iprobe = function(source, tag, comm = 1L) FALSE,
    .package = "npRmpi"
  )

  tasks <- list(
    list(start = 1L, bsz = 2L, seed = 100L),
    list(start = 3L, bsz = 2L, seed = 200L)
  )
  worker <- function(task) {
    matrix(seq_len(as.integer(task$bsz)), nrow = as.integer(task$bsz), ncol = 1L)
  }

  expect_error(
    run_fanout(
      tasks = tasks,
      worker = worker,
      ncol.out = 1L,
      what = "wild",
      profile.where = "unit.test:dispatch-timeout"
    ),
    "dispatch timeout waiting on worker results"
  )
})
