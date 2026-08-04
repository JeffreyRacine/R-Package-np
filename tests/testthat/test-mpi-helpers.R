test_that("npRmpi.init() is idempotent", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  # Do not require silence here: some MPI stacks print spawn/status messages.
  expect_true(isTRUE(try({ npRmpi.init(nslaves=1); TRUE }, silent=TRUE)))
  expect_true(mpi.comm.size(1) > 1)
})

test_that("npRmpi.quit() respects reuse/force", {
  if (.mpi_suite_pool_owned()) {
    subprocess_env <- npRmpi_subprocess_env()
    expect_false(is.null(subprocess_env))
    if (is.null(subprocess_env))
      return(invisible(NULL))

    result <- npRmpi_run_rscript_subprocess(c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.reuse.slaves = TRUE)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "npRmpi.quit(force = FALSE)",
      "stopifnot(mpi.comm.size(1) > 1L)",
      "stopifnot(isTRUE(getOption('npRmpi.pool.active', FALSE)))",
      "stopifnot(getFromNamespace('.npRmpi_has_active_slave_pool', 'npRmpi')())",
      "getFromNamespace('.npRmpi_require_active_slave_pool', 'npRmpi')(where = 'isolated soft-close probe')",
      "npRmpi.quit(force = TRUE)",
      "size <- try(mpi.comm.size(1), silent = TRUE)",
      "if (!inherits(size, 'try-error')) stopifnot(size < 2L)",
      "cat('NP_RMPI_ISOLATED_QUIT_OK\\n')"
    ), timeout = 90L, env = subprocess_env)

    expect_equal(result$status, 0L, info = paste(result$output, collapse = "\n"))
    expect_true(any(grepl(
      "NP_RMPI_ISOLATED_QUIT_OK", result$output, fixed = TRUE
    )))
    return(invisible(NULL))
  }

  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old <- getOption("npRmpi.reuse.slaves", FALSE)
  on.exit(options(npRmpi.reuse.slaves=old), add=TRUE)

  options(npRmpi.reuse.slaves=TRUE)
  expect_true(isTRUE(try({ npRmpi.quit(force=FALSE); TRUE }, silent=TRUE)))
  expect_true(mpi.comm.size(1) > 1) # soft-close keeps pool alive
  expect_true(isTRUE(getOption("npRmpi.pool.active", FALSE)))
  expect_true(getFromNamespace(".npRmpi_has_active_slave_pool", "npRmpi")())
  expect_silent(getFromNamespace(".npRmpi_require_active_slave_pool", "npRmpi")(
    where = "soft-close probe"
  ))

  expect_true(isTRUE(try({ npRmpi.quit(force=TRUE); TRUE }, silent=TRUE)))
  sz <- try(mpi.comm.size(1), silent=TRUE)
  if (!inherits(sz, "try-error"))
    expect_true(sz < 2)
})

test_that("npRmpi.session.info() returns useful fields", {
  # Can run without slaves; should not error.
  info <- npRmpi.session.info()
  expect_true(is.list(info))
  expect_true(all(c("npRmpi", "Rmpi", "reuse_slaves", "comm", "comm_size", "nslaves") %in% names(info)))
})
