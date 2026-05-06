test_that("conditional fixed localpoly bootstrap core fans out with an active pool", {
  core <- getFromNamespace(".np_inid_boot_from_conditional_localpoly_fixed_core", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fanout <- FALSE

  state <- list(n = 4L, neval = 2L, ngroups = 1L)
  counts <- matrix(
    c(
      1, 0, 2,
      0, 1, 0,
      1, 1, 1,
      2, 2, 1
    ),
    nrow = 4L,
    ncol = 3L
  )
  storage.mode(counts) <- "double"

  local_mocked_bindings(
    .np_inid_boot_from_conditional_localpoly_fixed_group_features = function(state, i) {
      list(rows = seq_len(state$neval))
    },
    .np_inid_boot_from_conditional_localpoly_fixed_t0 = function(state, feat) {
      c(0.25, 0.75)
    },
    .np_inid_boot_from_conditional_localpoly_fixed_fill_chunk = function(state, feat.list, counts.chunk) {
      bsz <- ncol(counts.chunk)
      matrix(seq_len(bsz * state$neval), nrow = bsz, ncol = state$neval)
    },
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
    .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, comm = 1L, include.master = TRUE) 2L,
    .npRmpi_bootstrap_chunk_tasks = function(B, chunk.size) {
      list(
        list(start = 1L, bsz = 2L, seed = 11L),
        list(start = 3L, bsz = 1L, seed = 12L)
      )
    },
    .npRmpi_bootstrap_run_fanout = function(tasks, worker, ncol.out, what = "bootstrap", ...) {
      called$fanout <- TRUE
      expect_identical(what, "inid-conditional-localpoly-fixed")
      expect_identical(ncol.out, state$neval)
      do.call(rbind, lapply(tasks, worker))
    },
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = TRUE)

  out <- core(state = state, B = ncol(counts), counts = counts, progress.label = "unit")

  expect_true(called$fanout)
  expect_identical(out$t0, c(0.25, 0.75))
  expect_identical(dim(out$t), c(3L, 2L))
})

test_that("conditional fixed localpoly bootstrap core remains local without a pool", {
  core <- getFromNamespace(".np_inid_boot_from_conditional_localpoly_fixed_core", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fill <- 0L

  state <- list(n = 3L, neval = 2L, ngroups = 1L)
  counts <- matrix(
    c(
      1, 0,
      1, 2,
      1, 1
    ),
    nrow = 3L,
    ncol = 2L
  )
  storage.mode(counts) <- "double"

  local_mocked_bindings(
    .np_inid_boot_from_conditional_localpoly_fixed_group_features = function(state, i) {
      list(rows = seq_len(state$neval))
    },
    .np_inid_boot_from_conditional_localpoly_fixed_t0 = function(state, feat) {
      c(0.5, 1.5)
    },
    .np_inid_boot_from_conditional_localpoly_fixed_fill_chunk = function(state, feat.list, counts.chunk) {
      called$fill <- called$fill + 1L
      bsz <- ncol(counts.chunk)
      matrix(rep(called$fill, bsz * state$neval), nrow = bsz, ncol = state$neval)
    },
    .npRmpi_bootstrap_run_fanout = function(...) {
      stop("fanout should not be used without an active pool")
    },
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = FALSE, np.messages = FALSE)

  out <- core(state = state, B = ncol(counts), counts = counts, progress.label = "unit")

  expect_identical(called$fill, 1L)
  expect_identical(out$t0, c(0.5, 1.5))
  expect_identical(dim(out$t), c(2L, 2L))
})
