test_that("runtime acceleration certificate follows delivered worker state", {
  e <- new.env(parent = environment(npRmpi.init))
  e$.npRmpi_option_agreement <- new.env(parent = emptyenv())
  e$mpi.comm.rank <- function(comm) 0L
  e$mpi.comm.size <- function(comm) 2L
  sends <- 0L; worker <- NULL; fail <- FALSE
  e$mpi.bcast.cmd <- function(..., rank, comm) {
    sends <<- sends + 1L
    if (fail) stop("injected delivery failure")
    worker <<- list(...)[["np.macMseries.accelerate"]]
  }
  f <- .npRmpi_sync_runtime_options_to_slaves; environment(f) <- e
  old <- options(); on.exit(options(old), add = TRUE)
  options(npRmpi.runtime.option.sync.active = FALSE)
  for (initial in list(FALSE, TRUE, "auto")) {
    options(np.macMseries.accelerate = initial)
    expect_true(f(force = TRUE)); expect_identical(worker, initial)
    snapshot <- options()
    options(np.macMseries.accelerate = if (identical(initial, FALSE)) TRUE else FALSE)
    expect_true(f()); expect_identical(worker, getOption("np.macMseries.accelerate"))
    options(snapshot)
    expect_true(f()); expect_identical(worker, initial)
    before <- sends; expect_false(f()); expect_identical(sends, before)
  }
  fail <- TRUE
  options(np.macMseries.accelerate = TRUE)
  expect_error(f(), "injected delivery failure")
  expect_null(e$.npRmpi_option_agreement$accelerate)
  expect_false(getOption("npRmpi.runtime.option.sync.active"))
  fail <- FALSE
  expect_true(f()); expect_identical(worker, TRUE)
  options(npRmpi.runtime.option.sync.active = TRUE)
  before <- sends; expect_false(f(force=TRUE)); expect_identical(sends,before)
  options(npRmpi.runtime.option.sync.active = FALSE)
  e$mpi.comm.rank <- function(comm) 1L
  expect_false(f(force=TRUE)); expect_identical(sends,before)
  e$mpi.comm.rank <- function(comm) 0L; e$mpi.comm.size <- function(comm) 1L
  expect_false(f(force=TRUE)); expect_identical(sends,before)
})
