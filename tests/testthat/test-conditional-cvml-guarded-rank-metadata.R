test_that("conditional guarded CVML counts one evaluation on every active rank", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(301925L)
  X <- as.data.frame(matrix(runif(200, -.8, .8), 100, 2))
  names(X) <- c("x1", "x2")
  Y <- data.frame(y = sin(X$x1) + rnorm(100, sd = .3))
  bw <- npcdensbw(xdat = X, ydat = Y, bws = c(.4, 2, 2),
    bwscaling = FALSE, bandwidth.compute = FALSE, regtype = "lp",
    degree = c(5L, 5L), basis = "glp", bernstein.basis = FALSE,
    bwmethod = "cv.ml")
  args <- list(xdat = X, ydat = Y, bws = bw,
               invalid.penalty = "dbmax", force.local = FALSE)
  command <- substitute(local({
    saved <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE,
      np.categorical.compress = FALSE, npRmpi.local.regression.mode = FALSE)
    on.exit(options(saved), add = TRUE)
    value <- do.call(get(".npcdensbw_eval_only", asNamespace("npRmpi")), ARGS)
    assign(".np_test_conditional_guard_value", value, envir = .GlobalEnv)
    value
  }), list(ARGS = args))
  on.exit({
    npRmpi:::mpi.remote.exec(
      rm(list = intersect(".np_test_conditional_guard_value",
                          ls(envir = .GlobalEnv, all.names = TRUE)),
         envir = .GlobalEnv), comm = 1L)
    if (exists(".np_test_conditional_guard_value", .GlobalEnv, inherits = FALSE))
      rm(".np_test_conditional_guard_value", envir = .GlobalEnv)
  }, add = TRUE, after = FALSE)
  value <- npRmpi:::.npRmpi_bcast_cmd_expr(command, comm = 1L, caller.execute = TRUE)
  # A remote command is also the completion barrier before reading worker state.
  remote <- npRmpi:::mpi.remote.exec(
    get(".np_test_conditional_guard_value", .GlobalEnv), comm = 1L)
  expect_equal(value$num.feval, 1)
  expect_equal(value$num.feval.guarded, 1)
  for (worker in remote) {
    expect_equal(worker$objective, value$objective, tolerance = 1e-10)
    expect_identical(worker$num.feval.guarded, value$num.feval.guarded)
    expect_identical(worker$num.feval, value$num.feval)
  }
  args$force.local <- TRUE
  local <- do.call(npRmpi:::.npcdensbw_eval_only, args)
  expect_equal(local$objective, value$objective, tolerance = 1e-10)
  expect_equal(local$num.feval.guarded, 1)
})
