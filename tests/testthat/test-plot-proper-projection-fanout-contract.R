test_that("conditional proper projection fans out matrix bootstrap values with an active pool", {
  project <- getFromNamespace(".npRmpi_project_conditional_proper_values", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fanout <- FALSE
  called$enabled <- FALSE

  values <- matrix(seq_len(6L), nrow = 3L, ncol = 2L)
  plan <- list(dummy = TRUE)

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, comm = 1L, include.master = TRUE) 2L,
    .npRmpi_bootstrap_fanout_enabled = function(...) {
      called$enabled <- TRUE
      TRUE
    },
    .npRmpi_bootstrap_chunk_tasks = function(B, chunk.size) {
      list(
        list(start = 1L, bsz = 2L),
        list(start = 3L, bsz = 1L)
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
      expect_identical(what, "proper-condens-project")
      expect_identical(ncol.out, ncol(values))
      expect_true(master_local_chunk)
      expect_true("values.mat" %in% names(required.bindings))
      do.call(rbind, lapply(tasks, worker))
    },
    .np_condens_project_values_with_plan = function(values, plan, progress.label = NULL) {
      data.matrix(values) + 10
    },
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = TRUE)

  out <- project(values, plan, cdf = FALSE)

  expect_true(called$enabled)
  expect_true(called$fanout)
  expect_equal(out, values + 10)
})

test_that("conditional proper projection uses distribution projection for cdf routes", {
  project <- getFromNamespace(".npRmpi_project_conditional_proper_values", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fanout <- FALSE

  values <- matrix(seq_len(4L), nrow = 2L, ncol = 2L)
  plan <- list(dummy = TRUE)

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, comm = 1L, include.master = TRUE) 1L,
    .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
    .npRmpi_bootstrap_chunk_tasks = function(B, chunk.size) {
      list(
        list(start = 1L, bsz = 1L),
        list(start = 2L, bsz = 1L)
      )
    },
    .npRmpi_bootstrap_run_fanout = function(tasks, worker, ncol.out, what = "bootstrap", ...) {
      called$fanout <- TRUE
      expect_identical(what, "proper-condist-project")
      expect_identical(ncol.out, ncol(values))
      do.call(rbind, lapply(tasks, worker))
    },
    .np_condist_project_values_with_plan = function(values, plan, progress.label = NULL) {
      data.matrix(values) + 20
    },
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = TRUE)

  out <- project(values, plan, cdf = TRUE)

  expect_true(called$fanout)
  expect_equal(out, values + 20)
})

test_that("conditional proper projection remains local for vector and no-pool routes", {
  project <- getFromNamespace(".npRmpi_project_conditional_proper_values", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fanout <- FALSE
  called$local <- 0L

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) FALSE,
    .npRmpi_bootstrap_run_fanout = function(...) {
      called$fanout <- TRUE
      stop("fanout should not be used without an active pool")
    },
    .np_condens_project_values_with_plan = function(values, plan, progress.label = NULL) {
      called$local <- called$local + 1L
      data.matrix(values) + 30
    },
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = TRUE)

  vector.out <- project(c(1, 2), list(dummy = TRUE), cdf = FALSE)
  matrix.out <- project(matrix(1:4, nrow = 2L), list(dummy = TRUE), cdf = FALSE)

  expect_false(called$fanout)
  expect_identical(called$local, 2L)
  expect_equal(as.numeric(vector.out), c(31, 32))
  expect_equal(matrix.out, matrix(1:4, nrow = 2L) + 30)
})
