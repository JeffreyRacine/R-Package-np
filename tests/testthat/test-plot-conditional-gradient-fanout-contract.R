test_that("conditional gradient bootstrap uses fanout with an active pool", {
  boot <- getFromNamespace(".npRmpi_inid_boot_from_conditional_gradient", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fanout <- FALSE

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, comm = 1L, include.master = TRUE) 2L,
    .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
    .npRmpi_bootstrap_chunk_tasks = function(B, chunk.size) {
      list(
        list(start = 1L, bsz = 2L, seed = 11L),
        list(start = 3L, bsz = 1L, seed = 12L)
      )
    },
    .npRmpi_bootstrap_run_fanout = function(tasks,
                                            worker,
                                            ncol.out,
                                            what = "bootstrap",
                                            master_local_chunk = FALSE,
                                            required.bindings = list(),
                                            ...) {
      called$fanout <- TRUE
      expect_identical(what, "conditional-gradient")
      expect_identical(ncol.out, 2L)
      expect_true(master_local_chunk)
      expect_false("counts.mat" %in% names(required.bindings))
      do.call(rbind, lapply(tasks, worker))
    },
    .np_inid_boot_from_conditional_gradient_local = function(xdat,
                                                            ydat,
                                                            exdat,
                                                            eydat,
                                                            bws,
                                                            B,
                                                            cdf,
                                                            gradient.index,
                                                            counts = NULL,
                                                            counts.drawer = NULL,
                                                            progress.label = NULL) {
      list(
        t = matrix(rep(seq_len(B), each = 2L), nrow = B, ncol = 2L, byrow = TRUE),
        t0 = c(0.25, 0.75)
      )
    },
    .npRmpi_with_local_bootstrap = function(expr) force(expr),
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = TRUE)

  out <- boot(
    xdat = data.frame(x = 1:3),
    ydat = data.frame(y = 1:3),
    exdat = data.frame(x = 1:2),
    eydat = data.frame(y = 1:2),
    bws = list(),
    B = 3L,
    cdf = FALSE,
    gradient.index = 1L
  )

  expect_true(called$fanout)
  expect_identical(out$t0, c(0.25, 0.75))
  expect_identical(dim(out$t), c(3L, 2L))
})

test_that("quantile gradient bootstrap uses fanout with fixed counts", {
  boot <- getFromNamespace(".npRmpi_inid_boot_from_quantile_gradient", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fanout <- FALSE
  counts <- matrix(c(1, 0, 2, 0, 2, 1), nrow = 3L, ncol = 2L)

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, comm = 1L, include.master = TRUE) 1L,
    .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
    .npRmpi_bootstrap_chunk_tasks = function(B, chunk.size) {
      list(
        list(start = 1L, bsz = 1L, seed = 11L),
        list(start = 2L, bsz = 1L, seed = 12L)
      )
    },
    .npRmpi_bootstrap_run_fanout = function(tasks,
                                            worker,
                                            ncol.out,
                                            what = "bootstrap",
                                            required.bindings = list(),
                                            ...) {
      called$fanout <- TRUE
      expect_identical(what, "quantile-gradient")
      expect_identical(ncol.out, 2L)
      expect_true("counts.mat" %in% names(required.bindings))
      do.call(rbind, lapply(tasks, worker))
    },
    .np_inid_boot_from_quantile_gradient_local = function(xdat,
                                                         ydat,
                                                         exdat,
                                                         bws,
                                                         B,
                                                         tau,
                                                         gradient.index,
                                                         counts = NULL,
                                                         counts.drawer = NULL,
                                                         progress.label = NULL) {
      list(
        t = matrix(rep(colSums(counts), each = 2L), nrow = B, ncol = 2L, byrow = TRUE),
        t0 = c(1.25, 1.75)
      )
    },
    .npRmpi_with_local_bootstrap = function(expr) force(expr),
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = TRUE)

  out <- boot(
    xdat = data.frame(x = 1:3),
    ydat = 1:3,
    exdat = data.frame(x = 1:2),
    bws = list(),
    B = 2L,
    tau = 0.5,
    gradient.index = 1L,
    counts = counts
  )

  expect_true(called$fanout)
  expect_identical(out$t0, c(1.25, 1.75))
  expect_identical(dim(out$t), c(2L, 2L))
})

test_that("quantile level bootstrap uses fanout with fixed counts", {
  boot <- getFromNamespace(".npRmpi_inid_boot_from_quantile_level", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fanout <- FALSE
  counts <- matrix(c(1, 0, 2, 0, 2, 1), nrow = 3L, ncol = 2L)

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, comm = 1L, include.master = TRUE) 1L,
    .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
    .npRmpi_bootstrap_chunk_tasks = function(B, chunk.size) {
      list(
        list(start = 1L, bsz = 1L, seed = 11L),
        list(start = 2L, bsz = 1L, seed = 12L)
      )
    },
    .npRmpi_bootstrap_run_fanout = function(tasks,
                                            worker,
                                            ncol.out,
                                            what = "bootstrap",
                                            required.bindings = list(),
                                            ...) {
      called$fanout <- TRUE
      expect_identical(what, "quantile-level")
      expect_identical(ncol.out, 2L)
      expect_true("counts.mat" %in% names(required.bindings))
      do.call(rbind, lapply(tasks, worker))
    },
    .np_inid_boot_from_quantile_level_local = function(xdat,
                                                       ydat,
                                                       exdat,
                                                       bws,
                                                       B,
                                                       tau,
                                                       counts = NULL,
                                                       counts.drawer = NULL,
                                                       progress.label = NULL) {
      list(
        t = matrix(rep(colSums(counts), each = 2L), nrow = B, ncol = 2L, byrow = TRUE),
        t0 = c(2.25, 2.75)
      )
    },
    .npRmpi_with_local_bootstrap = function(expr) force(expr),
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = TRUE)

  out <- boot(
    xdat = data.frame(x = 1:3),
    ydat = 1:3,
    exdat = data.frame(x = 1:2),
    bws = list(),
    B = 2L,
    tau = 0.5,
    counts = counts
  )

  expect_true(called$fanout)
  expect_identical(out$t0, c(2.25, 2.75))
  expect_identical(dim(out$t), c(2L, 2L))
})
