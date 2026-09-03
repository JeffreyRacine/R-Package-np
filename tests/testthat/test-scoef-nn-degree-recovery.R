test_that("NN smooth-coefficient degree search recovers and fails closed", {
  # Own a fresh MPI process; never commandeer the suite's rank-local pool.
  if (!identical(Sys.getenv("NP_RMPI_SCOEF_DEGREE_CHILD", ""), "1")) {
    helper <- normalizePath(testthat::test_path("helper-mpi.R"))
    test <- normalizePath(testthat::test_path("test-scoef-nn-degree-recovery.R"))
    expected <- paste(capture.output(dput(deparse(body(.npscoefbw_nomad_search)))),
                      collapse = "\n")
    out <- npRmpi_run_rscript_subprocess(c(
      "suppressPackageStartupMessages(library(npRmpi))",
      sprintf("stopifnot(identical(deparse(body(npRmpi:::.npscoefbw_nomad_search)), %s))",
              expected),
      "env <- new.env(parent = asNamespace('npRmpi'))",
      sprintf("sys.source(%s, envir = env)", dQuote(helper, FALSE)),
      sprintf("testthat::test_file(%s, env = env, reporter = 'silent', stop_on_failure = TRUE)",
              dQuote(test, FALSE)),
      "cat('NN_SCOEF_DEGREE_ISOLATED_PASS\\n')"
    ), env = c(sprintf("R_LIBS=%s", paste(.libPaths(), collapse = .Platform$path.sep)),
      "_R_CHECK_PACKAGE_NAME_=", "NP_RMPI_TEST_SUITE_POOL=",
      "NP_RMPI_TEST_SUITE_LOCAL_MODE=", "NP_RMPI_NO_REUSE_SLAVES=1",
      "NP_RMPI_SCOEF_DEGREE_CHILD=1"))
    out$witnessed <- any(grepl("NN_SCOEF_DEGREE_ISOLATED_PASS", out$output, fixed = TRUE))
    expect_false(is.null(out))
    expect_identical(out$status, 0L, info = paste(out$output, collapse = "\n"))
    expect_true(out$witnessed)
    return(invisible(NULL))
  }

  expect_true(spawn_mpi_slaves(1L))
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE,
                 np.nomad.degree.start.policy = "low_first_full_random")
  on.exit({
    close_mpi_slaves(force = TRUE)
    options(old)
  })
  n <- 24L
  xdat <- data.frame(x = seq_len(n) / n)
  zdat <- data.frame(
    u = factor(rep(c("a", "b"), n / 2L)),
    z = c(rep(0, 16L), seq_len(8L))
  )
  y <- 1 + xdat$x + 0.2 * zdat$z + sin(seq_len(n)) / 10
  common <- list(
    xdat = xdat, ydat = y, zdat = zdat, regtype = "lp", nomad = TRUE,
    degree.min = 0L, degree.max = 1L, nmulti = 1L,
    random.seed = 2026090304L, nomad.opts = list(MAX_BB_EVAL = 1L)
  )

  set.seed(30906L)
  before <- .Random.seed
  gnn <- do.call(npscoefbw, c(common, list(
    bwtype = "generalized_nn", search.engine = "nomad"
  )))
  expect_identical(.Random.seed, before)
  starts <- gnn$degree.search$restart.starts
  expect_length(starts, 2L)
  expect_gt(starts[[2L]][1L], starts[[1L]][1L])
  expect_identical(starts[[2L]][-1L], starts[[1L]][-1L])
  expect_false(gnn$degree.search$certified)
  ctx <- .npscoefbw_nomad_context_prepare(xdat, y, zdat)
  pool <- .npscoefbw_nomad_pool_start(ctx)
  raw <- .npscoefbw_eval_pool(ctx, gnn, pool, invalid.penalty = "large")
  .npscoefbw_nomad_pool_stop(pool)
  .npscoefbw_nomad_context_cleanup(ctx)
  expect_true(raw$raw.valid)
  expect_identical(as.numeric(raw$objective), as.numeric(gnn$fval))

  set.seed(30906L)
  before <- .Random.seed
  ann <- do.call(npscoefbw, c(common, list(
    bwtype = "adaptive_nn", search.engine = "nomad+powell"
  )))
  expect_identical(.Random.seed, before)
  expect_length(ann$degree.search$restart.starts, 2L)
  expect_true(is.finite(ann$fval))
  expect_true(is.finite(ann$degree.search$powell.time))

  explicit <- c(common, list(
    bws = c(0.25, 5), bwtype = "generalized_nn",
    search.engine = "nomad+powell"
  ))
  set.seed(30906L)
  before <- .Random.seed
  expect_error(do.call(npscoefbw, explicit),
               "npscoefbw NOMAD degree search did not return a raw-valid solution",
               fixed = TRUE)
  expect_identical(.Random.seed, before)

  manual <- as.call(c(list(as.name("npscoefbw")), c(common, list(
    bws = c(0.25, 5), bwtype = "adaptive_nn",
    search.engine = "nomad+powell"
  ))))
  expect_error(.npRmpi_bcast_cmd_expr(manual, caller.execute = TRUE),
               "npscoefbw NOMAD degree search did not return a raw-valid solution",
               fixed = TRUE)
  alive <- mpi.bcast.cmd(npRmpi:::mpi.allgather.Robj(npRmpi:::mpi.comm.rank(1L)),
                         caller.execute = TRUE)
  expect_identical(as.integer(alive), 0:1L)
})
