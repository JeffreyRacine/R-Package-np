auto_beta_arguments <- function(dat, degree.start, bws = NULL) {
  arguments <- list(
    formula = y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 1L,
    degree.start = degree.start,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    ckertype = "beta",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 10L)
  )
  if (!is.null(bws))
    arguments[["bws"]] <- bws
  arguments
}

auto_beta_replay <- function(dat, bandwidth) {
  getFromNamespace(".npregbw_eval_only", "npRmpi")(
    xdat = dat["x"],
    ydat = dat[["y"]],
    bws = bandwidth,
    invalid.penalty = "baseline"
  )[["objective"]][[1L]]
}

test_that("automatic beta regression owns startup and degree transitions", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    npRmpi.autodispatch = TRUE
  )
  on.exit(options(old), add = TRUE)

  x <- seq(0.04, 0.96, length.out = 32L)
  dat <- data.frame(y = sin(2 * pi * x) + 0.04 * cos(7 * pi * x), x = x)

  missing.start <- do.call(
    npregbw,
    auto_beta_arguments(dat, degree.start = 1L)
  )
  missing.replay <- auto_beta_replay(dat, missing.start)
  expect_true(all(is.finite(missing.start[["bw"]])))
  expect_true(all(missing.start[["bw"]] > 0))
  expect_identical(
    as.double(missing.start[["fval"]][[1L]]),
    as.double(missing.replay)
  )

  automatic <- do.call(
    npregbw,
    auto_beta_arguments(dat, degree.start = 0L, bws = 0.17)
  )
  restart <- automatic[["degree.search"]][["restart.results"]][[1L]]
  fixed <- npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 0L,
    bandwidth.compute = FALSE,
    bws = 0.17,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    ckertype = "beta",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )
  fixed.objective <- auto_beta_replay(dat, fixed)
  automatic.replay <- auto_beta_replay(dat, automatic)

  expect_identical(as.integer(restart[["first_degree"]]), 0L)
  expect_identical(
    as.double(restart[["first_objective"]][[1L]]),
    as.double(fixed.objective)
  )
  expect_identical(
    as.double(automatic[["fval"]][[1L]]),
    as.double(automatic.replay)
  )

  local.regression <- getFromNamespace(
    ".npRmpi_with_local_regression", "npRmpi"
  )
  local.automatic <- local.regression(do.call(
    npregbw,
    auto_beta_arguments(dat, degree.start = 0L, bws = 0.17)
  ))
  local.trace <- local.automatic[["degree.search"]][["trace"]]
  local.replay <- local.regression(auto_beta_replay(dat, local.automatic))

  expect_identical(as.integer(local.trace[["degree"]][[1L]]), 0L)
  expect_identical(
    as.double(local.trace[["fval"]][[1L]]),
    as.double(fixed.objective)
  )
  expect_true(any(as.integer(local.trace[["degree"]]) == 1L))
  expect_identical(
    as.double(local.automatic[["fval"]][[1L]]),
    as.double(local.replay)
  )
})
