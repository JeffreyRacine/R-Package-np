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
