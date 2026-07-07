native_cache_off_opts <- list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 1)
native_cache_unknown_opts <- list(EVAL_USE_CACHE = NA, MAX_BB_EVAL = 1)

test_that("native NOMAD option builders reject cache-off solves before solver entry", {
  builders <- list(
    npRmpi:::.np_nomad_native_option_vectors,
    npRmpi:::.npudensbw_nomad_native_option_vectors,
    npRmpi:::.npcdensbw_nomad_shadow_native_option_vectors,
    npRmpi:::.npcdistbw_nomad_native_option_vectors,
    npRmpi:::.npregbw_nomad_native_option_vectors
  )

  for (builder in builders) {
    expect_error(
      builder(native_cache_off_opts),
      "requires EVAL_USE_CACHE = TRUE",
      fixed = TRUE
    )
    expect_error(
      builder(native_cache_unknown_opts),
      "requires EVAL_USE_CACHE = TRUE",
      fixed = TRUE
    )
    expect_no_error(builder(list(EVAL_USE_CACHE = TRUE, MAX_BB_EVAL = 1)))
  }
})

test_that("native R callback failure is contained in an isolated process", {
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env))

  out <- npRmpi_run_rscript_subprocess(
    c(
      "library(npRmpi)",
      "bad <- npRmpi:::.np_nomad_native_r_callback_search(",
      "  eval.f = function(x) stop('intentional package R callback failure'),",
      "  x0 = 0,",
      "  bbin = 0L,",
      "  lb = -1,",
      "  ub = 1,",
      "  random.seed = 42L,",
      "  option.names = c('MAX_BB_EVAL', 'NB_THREADS_PARALLEL_EVAL'),",
      "  option.values = c('1', '1')",
      ")",
      "stopifnot(identical(bad$value$status, 'error'))",
      "stopifnot(grepl('callback', bad$value$message, ignore.case = TRUE))",
      "good <- npRmpi:::.np_nomad_native_r_callback_search(",
      "  eval.f = function(x) sum((x - 0.25)^2),",
      "  x0 = 0,",
      "  bbin = 0L,",
      "  lb = -1,",
      "  ub = 1,",
      "  random.seed = 42L,",
      "  option.names = c('MAX_BB_EVAL', 'NB_THREADS_PARALLEL_EVAL'),",
      "  option.values = c('10', '1')",
      ")",
      "stopifnot(identical(good$value$status, 'ok'))",
      "stopifnot(identical(good$value$native_status, 0L))",
      "stopifnot(identical(good$value$result_status, 0L))",
      "stopifnot(is.finite(good$value$objective))",
      "cat('NATIVE_R_CALLBACK_SUBPROCESS_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  expect_equal(out$status, 0L, info = paste(out$output, collapse = "\n"))
  expect_true(any(grepl("NATIVE_R_CALLBACK_SUBPROCESS_OK", out$output, fixed = TRUE)),
              info = paste(out$output, collapse = "\n"))
})

test_that("native npreg solver option failure does not poison the next MPI solve", {
  env <- npRmpi_subprocess_env("NP_RMPI_NO_REUSE_SLAVES=1")
  skip_if(is.null(env))

  out <- npRmpi_run_rscript_subprocess(
    c(
      "library(npRmpi)",
      "options(np.messages = FALSE, np.tree = FALSE, np.developer.native.nomad.diagnostics = TRUE)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "set.seed(111)",
      "n <- 40L",
      "dat <- data.frame(x = runif(n))",
      "dat$y <- sin(dat$x) + rnorm(n, sd = 0.2)",
      "err <- tryCatch(npregbw(y ~ x, data = dat, nomad = TRUE, nmulti = 1L, degree.max = 2L, nomad.opts = list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 5)), error = function(e) e)",
      "stopifnot(inherits(err, 'error'))",
      "stopifnot(grepl('requires EVAL_USE_CACHE = TRUE', conditionMessage(err), fixed = TRUE))",
      "fit <- npregbw(y ~ x, data = dat, nomad = TRUE, nmulti = 1L, degree.max = 2L, nomad.opts = list(MAX_BB_EVAL = 8))",
      "stopifnot(is.finite(fit$fval))",
      "stopifnot(!is.null(attr(fit, 'native.nomad.diagnostics')))",
      "quit_status <- tryCatch({ npRmpi.quit(); 0L }, error = function(e) 1L)",
      "stopifnot(identical(quit_status, 0L))",
      "quit(save = 'no', status = 0L, runLast = FALSE)"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(out$status, 0L, info = paste(out$output, collapse = "\n"))
})

test_that("native bandwidth option failures do not poison next MPI solves", {
  env <- npRmpi_subprocess_env("NP_RMPI_NO_REUSE_SLAVES=1")
  skip_if(is.null(env))

  out <- npRmpi_run_rscript_subprocess(
    c(
      "library(npRmpi)",
      "options(np.messages = FALSE, np.tree = FALSE, np.developer.native.nomad.diagnostics = TRUE)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "set.seed(913)",
      "n <- 30L",
      "x <- data.frame(x = runif(n))",
      "y <- data.frame(y = rnorm(n))",
      "check_error <- function(expr) {",
      "  err <- tryCatch(eval.parent(substitute(expr)), error = function(e) e)",
      "  stopifnot(inherits(err, 'error'))",
      "  stopifnot(grepl('requires EVAL_USE_CACHE = TRUE', conditionMessage(err), fixed = TRUE))",
      "}",
      "check_error(npudensbw(dat = x, bwsolver = 'mads', nmulti = 1L, nomad.opts = list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 5)))",
      "stopifnot(is.finite(npudensbw(dat = x, bwsolver = 'mads', nmulti = 1L, nomad.opts = list(MAX_BB_EVAL = 8))$fval))",
      "check_error(npudistbw(dat = x, bwsolver = 'mads', nmulti = 1L, nomad.opts = list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 5)))",
      "stopifnot(is.finite(npudistbw(dat = x, bwsolver = 'mads', nmulti = 1L, nomad.opts = list(MAX_BB_EVAL = 8))$fval))",
      "check_error(npcdensbw(ydat = y, xdat = x, bwsolver = 'mads', nmulti = 1L, nomad.opts = list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 5)))",
      "stopifnot(is.finite(npcdensbw(ydat = y, xdat = x, bwsolver = 'mads', nmulti = 1L, nomad.opts = list(MAX_BB_EVAL = 8))$fval))",
      "check_error(npcdistbw(ydat = y, xdat = x, bwsolver = 'mads', nmulti = 1L, nomad.opts = list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 5)))",
      "stopifnot(is.finite(npcdistbw(ydat = y, xdat = x, bwsolver = 'mads', nmulti = 1L, nomad.opts = list(MAX_BB_EVAL = 8))$fval))",
      "quit_status <- tryCatch({ npRmpi.quit(); 0L }, error = function(e) 1L)",
      "stopifnot(identical(quit_status, 0L))",
      "quit(save = 'no', status = 0L, runLast = FALSE)"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(out$status, 0L, info = paste(out$output, collapse = "\n"))
})
