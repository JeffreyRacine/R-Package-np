test_that("NN degree evaluator preserves errors and reuses its worker pool", {
  # Own a fresh service pool; never commandeer the suite's rank-local pool.
  if(!identical(Sys.getenv("NP_RMPI_SCOEF_STATUS_CHILD", ""), "1")) {
    helper <- normalizePath(testthat::test_path("helper-mpi.R"))
    test <- normalizePath(testthat::test_path("test-scoef-nn-evaluator-status.R"))
    expected <- paste(capture.output(dput(deparse(body(.npscoefbw_eval_pool)))), collapse = "\n")
    out <- npRmpi_run_rscript_subprocess(c(
      "suppressPackageStartupMessages(library(npRmpi))",
      sprintf("stopifnot(identical(deparse(body(npRmpi:::.npscoefbw_eval_pool)), %s))", expected),
      "env <- new.env(parent = asNamespace('npRmpi'))",
      sprintf("sys.source(%s, envir = env)", dQuote(helper, FALSE)),
      sprintf("checks <- testthat::test_file(%s, env = env, reporter = 'silent', stop_on_failure = TRUE)",
              dQuote(test, FALSE)),
      "stopifnot(sum(as.data.frame(checks)$passed) == 10L)",
      "cat('NN_EVALUATOR_STATUS_ISOLATED_PASS\\n')"
    ), env = c(sprintf("R_LIBS=%s", paste(.libPaths(), collapse = .Platform$path.sep)),
      "_R_CHECK_PACKAGE_NAME_=", "NP_RMPI_TEST_SUITE_POOL=",
      "NP_RMPI_TEST_SUITE_LOCAL_MODE=", "NP_RMPI_NO_REUSE_SLAVES=1",
      "NP_RMPI_SCOEF_STATUS_CHILD=1"))
    out$witnessed <- any(grepl("NN_EVALUATOR_STATUS_ISOLATED_PASS", out$output, fixed = TRUE))
    expect_false(is.null(out))
    expect_identical(out$status, 0L, info = paste(out$output, collapse = "\n"))
    expect_true(out$witnessed)
    return(invisible(NULL))
  }
  expect_true(spawn_mpi_slaves(1L))
  pool <- ctx <- NULL
  traced <- FALSE
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit({
    if(!is.null(pool)) .npscoefbw_nomad_pool_stop(pool)
    if(!is.null(ctx)) .npscoefbw_nomad_context_cleanup(ctx)
    if(traced) mpi.bcast.cmd(untrace(".npscoefbw_nomad_moment_state",
      where = asNamespace("npRmpi")), caller.execute = TRUE)
    close_mpi_slaves(force = TRUE)
    options(old)
  })
  i <- seq_len(24L)
  x <- data.frame(w = cos(i / 5) + i / 24)
  z <- data.frame(z = sin(i * sqrt(2)) + i / 24)
  y <- (1 + z$z / 4) * x$w + sin(i * sqrt(5)) / 20
  bw <- .npscoefbw_build_scbandwidth(x, y, z, 12, FALSE,
    list(regtype = "lp", degree = 1L, bwtype = "generalized_nn"))
  bw <- .npscoefbw_normalize_nomad_scbw(bw, z, 12)
  ctx <- .npscoefbw_nomad_context_prepare(x, y, z)
  good <- .npscoefbw_nomad_eval_direct(ctx, bw)
  expect_true(good$raw.valid)
  tied <- data.frame(z = rep(0:1, each = 12L))
  bad.ctx <- .npscoefbw_nomad_lp_context_local(x, y, tied)
  bad <- .npscoefbw_normalize_nomad_scbw(bw, tied, 2)
  rejected <- .npscoefbw_nomad_eval_direct(bad.ctx, bad)
  expect_false(rejected$raw.valid)
  expect_identical(rejected$objective, 10)

  unknown <- simpleError("unexpected NN evaluator failure")
  expect_identical(.npscoefbw_nomad_unknown_nn_error(unknown, bw), unknown)
  fixed <- bw; fixed$type <- "fixed"
  expect_null(.npscoefbw_nomad_unknown_nn_error(unknown, fixed))
  typed <- .np_nn_candidate_invalid_condition("typed invalid", owner = "test")
  expect_null(.npscoefbw_nomad_unknown_nn_error(typed, bw))

  mpi.bcast.cmd({
    assign(".np_test_scoef_status_once", TRUE, envir = .GlobalEnv)
    trace(".npscoefbw_nomad_moment_state", where = asNamespace("npRmpi"),
      print = FALSE, tracer = quote({
        if(npRmpi:::mpi.comm.rank(1L) == 1L &&
           isTRUE(get0(".np_test_scoef_status_once", envir = .GlobalEnv))) {
          assign(".np_test_scoef_status_once", FALSE, envir = .GlobalEnv)
          stop(simpleError("unexpected NN evaluator failure"))
        }
      }))
  }, caller.execute = TRUE)
  traced <- TRUE
  pool <- .npscoefbw_nomad_pool_start(ctx)
  expect_error(.npscoefbw_eval_pool(ctx, bw, pool),
               "unexpected NN evaluator failure", fixed = TRUE)
  reused <- .npscoefbw_eval_pool(ctx, bw, pool)
  expect_true(reused$raw.valid)
  expect_equal(reused$objective, good$objective, tolerance = 1e-11)
})
