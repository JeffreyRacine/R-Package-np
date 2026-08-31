test_that("native npreg raw-objective validity excludes the DBL_MAX sentinel", {
  raw_valid <- getFromNamespace(
    ".npregbw_native_raw_objective_valid", "np"
  )

  expect_true(raw_valid(0))
  expect_true(raw_valid(-1))
  expect_false(raw_valid(.Machine$double.xmax))
  expect_false(raw_valid(Inf))
  expect_false(raw_valid(NA_real_))
  expect_false(raw_valid(numeric()))
})

test_that("native npreg rejects a search with no raw-valid restart", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = seq(-1, 1, length.out = 24L))
  expect_error(
    npregbw(
      xdat = x,
      ydat = rep(0, nrow(x)),
      bwmethod = "cv.aic",
      regtype = "lc",
      bwsolver = "mads",
      nmulti = 2L,
      invalid.penalty = "baseline",
      nomad.opts = list(MAX_BB_EVAL = 4L)
    ),
    "native npreg NOMAD route did not return a raw-valid solution",
    fixed = TRUE
  )
})

test_that("native npreg excludes an invalid restart and certifies the next", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.developer.native.nomad.diagnostics = TRUE
  )
  on.exit(options(old), add = TRUE)

  set.seed(915L)
  x <- data.frame(x = seq(-1, 1, length.out = 24L))
  y <- sin(3 * x$x) + seq_len(nrow(x)) / 1000
  fit <- npregbw(
    xdat = x,
    ydat = y,
    bws = 1e-9,
    regtype = "ll",
    bwmethod = "cv.ls",
    bwtype = "fixed",
    ckertype = "epanechnikov",
    bwsolver = "mads",
    nmulti = 2L,
    invalid.penalty = "baseline",
    nomad.opts = list(MAX_BB_EVAL = 1L)
  )

  raw <- np:::.npregbw_eval_only(
    xdat = x,
    ydat = y,
    bws = fit,
    invalid.penalty = "dbmax",
    penalty.multiplier = 10
  )
  first <- fit$nomad.restart.results[[1L]]
  second <- fit$nomad.restart.results[[2L]]

  expect_identical(fit$nomad.best.restart, 2L)
  expect_identical(as.numeric(first$objective), .Machine$double.xmax)
  expect_identical(as.numeric(second$objective), as.numeric(fit$fval))
  expect_identical(as.numeric(fit$fval), as.numeric(raw$objective))
  for (restart in fit$nomad.restart.results) {
    expect_identical(
      as.numeric(restart$native$total_num.feval),
      as.numeric(restart$native$compiled_callback_calls) + 1
    )
  }
})
