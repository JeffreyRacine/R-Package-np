# Unique default-lane contract: rank-different conditional-hat row tasks must
# enter the existing scoped local R/native owner and restore it on failure.
# Numerical MPI coverage stays in the external focused subprocess matrix.
test_that("conditional hat tasks use and restore the scoped local owner", {
  ns <- asNamespace("npRmpi")
  fanout <- get(".npRmpi_cdhat_apply_fanout", ns)
  local.owner <- get(".npRmpi_with_local_regression", ns)
  env <- new.env(parent = ns)
  state <- new.env(parent = emptyenv())
  state$native <- FALSE
  state$fail <- FALSE
  state$entered <- 0L
  env$.Call <- function(name, active, PACKAGE) {
    stopifnot(identical(name, "C_np_set_local_regression_mode"),
              identical(PACKAGE, "npRmpi"))
    old <- state$native
    state$native <- active
    old
  }
  environment(local.owner) <- env
  env$.npRmpi_with_local_regression <- local.owner
  env$mpi.comm.rank <- function(comm) 0L
  env$.npRmpi_bootstrap_worker_count <- function(comm) 2L
  env$.npRmpi_autodispatch_called_from_bcast <- function() FALSE
  env$.npcdhat_exact_apply <- function(bws, txdat, tydat, exdat, eydat,
                                      rhs, operator, x.s, ...) {
    expect_true(isTRUE(getOption("npRmpi.local.regression.mode")))
    expect_true(isTRUE(getOption("npRmpi.autodispatch.disable")))
    expect_true(isTRUE(getOption("npRmpi.autodispatch.context")))
    expect_true(state$native)
    state$entered <- state$entered + 1L
    if (state$fail) stop("deliberate conditional task failure")
    matrix(exdat$x, ncol = ncol(rhs))
  }
  env$.npRmpi_bootstrap_run_fanout <- function(tasks, worker, ...) {
    do.call(rbind, lapply(tasks, worker))
  }
  environment(fanout) <- env
  withr::local_options(list(npRmpi.cdhat.apply.fanout = TRUE,
    npRmpi.cdhat.apply.fanout.min.work = 0,
    npRmpi.mpi.initialized = TRUE, npRmpi.autodispatch.disable = FALSE,
    npRmpi.autodispatch.context = FALSE, npRmpi.local.regression.mode = NULL))
  args <- list(bws = list(), txdat = data.frame(x = 1:8),
    tydat = data.frame(y = 1:8), exdat = data.frame(x = 1:5),
    eydat = data.frame(y = 1:5), rhs = matrix(1, 8, 2), operator = "normal")
  for (native.before in c(FALSE, TRUE)) {
    state$native <- native.before
    state$fail <- FALSE
    expect_equal(do.call(fanout, args), cbind(1:5, 1:5))
    expect_identical(state$native, native.before)
    expect_null(getOption("npRmpi.local.regression.mode"))
    expect_false(getOption("npRmpi.autodispatch.disable"))
    expect_false(getOption("npRmpi.autodispatch.context"))
    state$fail <- TRUE
    expect_error(do.call(fanout, args), "deliberate conditional task failure", fixed = TRUE)
    expect_identical(state$native, native.before)
    expect_null(getOption("npRmpi.local.regression.mode"))
    expect_false(getOption("npRmpi.autodispatch.disable"))
    expect_false(getOption("npRmpi.autodispatch.context"))
  }
  expect_gt(state$entered, 0L)
})
