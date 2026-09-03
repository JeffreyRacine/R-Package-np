test_that("shared raw-objective validity excludes the DBL_MAX sentinel", {
  package <- packageName(environment(npregbw))
  raw_valid <- getFromNamespace(
    ".np_nn_raw_objective_valid", package
  )

  expect_true(raw_valid(0))
  expect_true(raw_valid(-1))
  expect_false(raw_valid(.Machine$double.xmax))
  expect_false(raw_valid(Inf))
  expect_false(raw_valid(NA_real_))
  expect_false(raw_valid(numeric()))
})

test_that("native npreg raw certification is MPI-session symmetric", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  env <- npRmpi_subprocess_env(
    c("NP_RMPI_NO_REUSE_SLAVES=1", "_R_CHECK_PACKAGE_NAME_=")
  )
  skip_if(is.null(env), "local npRmpi install unavailable for MPI contract")

  out <- npRmpi_run_rscript_subprocess(
    c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(np.messages = FALSE, np.tree = FALSE)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "x <- data.frame(x = seq(-1, 1, length.out = 24L))",
      "invalid <- tryCatch(npregbw(xdat = x, ydat = rep(0, nrow(x)), bwmethod = 'cv.aic', regtype = 'lc', bwsolver = 'mads', nmulti = 2L, invalid.penalty = 'baseline', nomad.opts = list(MAX_BB_EVAL = 4L)), error = identity)",
      "stopifnot(inherits(invalid, 'error'))",
      "stopifnot(grepl('native npreg NOMAD route did not return a raw-valid solution', conditionMessage(invalid), fixed = TRUE))",
      "set.seed(915L)",
      "y <- sin(3 * x$x) + seq_len(nrow(x)) / 1000",
      "fit <- npregbw(xdat = x, ydat = y, bws = 1e-9, regtype = 'll', bwmethod = 'cv.ls', bwtype = 'fixed', ckertype = 'epanechnikov', bwsolver = 'mads', nmulti = 2L, invalid.penalty = 'baseline', nomad.opts = list(MAX_BB_EVAL = 1L))",
      "first <- fit$nomad.restart.results[[1L]]",
      "second <- fit$nomad.restart.results[[2L]]",
      "stopifnot(identical(fit$nomad.best.restart, 2L))",
      "stopifnot(identical(as.numeric(first$objective), .Machine$double.xmax))",
      "stopifnot(identical(as.numeric(second$objective), as.numeric(fit$fval)))",
      "for (restart in fit$nomad.restart.results) stopifnot(identical(as.numeric(restart$native$total_num.feval), as.numeric(restart$native$compiled_callback_calls) + 1))",
      "cat('NATIVE_NPREG_RAW_CERTIFICATION_MPI1_OK\\n')",
      "npRmpi.quit()",
      "quit(save = 'no', status = 0L, runLast = FALSE)"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(out$status, 0L, info = paste(out$output, collapse = "\n"))
  expect_true(
    any(grepl(
      "NATIVE_NPREG_RAW_CERTIFICATION_MPI1_OK",
      out$output,
      fixed = TRUE
    )),
    info = paste(out$output, collapse = "\n")
  )
})
