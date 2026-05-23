test_that("np.powell.cache is included in MPI option sync", {
  expect_true("np.powell.cache" %in% npRmpi:::.npRmpi_autodispatch_option_keys())
})

test_that("np.powell.cache default preserves legacy environment off switch", {
  old <- Sys.getenv("NP_NN_POWELL_CACHE_INSTRUMENT", unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("NP_NN_POWELL_CACHE_INSTRUMENT")
    } else {
      Sys.setenv(NP_NN_POWELL_CACHE_INSTRUMENT = old)
    }
  }, add = TRUE)

  Sys.setenv(NP_NN_POWELL_CACHE_INSTRUMENT = "off")
  expect_false(npRmpi:::.np_powell_cache_default())

  Sys.setenv(NP_NN_POWELL_CACHE_INSTRUMENT = "on")
  expect_true(npRmpi:::.np_powell_cache_default())
})

test_that("np.powell.cache controls continuous NN Powell caching under MPI", {
  env <- npRmpi_subprocess_env("NP_RMPI_NO_REUSE_SLAVES=1")
  skip_if(is.null(env))

  out <- npRmpi_run_rscript_subprocess(
    c(
      "library(npRmpi)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "run_bw <- function(cache) {",
      "  set.seed(42)",
      "  dat <- data.frame(x1 = runif(80), x2 = runif(80))",
      "  dat$y <- dat$x1 + dat$x2 + rnorm(80)",
      "  old <- options(np.messages = FALSE, np.tree = FALSE, np.powell.cache = cache)",
      "  on.exit(options(old), add = TRUE)",
      "  npregbw(y ~ x1 + x2, data = dat, bwmethod = 'cv.ls', bwtype = 'generalized_nn', regtype = 'lc', nmulti = 1)",
      "}",
      "cached <- run_bw(TRUE)",
      "uncached <- run_bw(FALSE)",
      "stopifnot(isTRUE(all.equal(cached$bw, uncached$bw)))",
      "stopifnot(isTRUE(all.equal(cached$fval, uncached$fval, tolerance = 1e-12)))",
      "stopifnot(identical(unname(cached$nn.cache[['enabled']]), 1))",
      "stopifnot(unname(cached$nn.cache[['hits']]) > 0)",
      "stopifnot(cached$num.feval.fast >= unname(cached$nn.cache[['hits']]))",
      "stopifnot(identical(unname(uncached$nn.cache[['enabled']]), 0))",
      "stopifnot(identical(unname(uncached$nn.cache[['hits']]), 0))"
    ),
    timeout = 45L,
    env = env
  )

  expect_equal(out$status, 0L, info = paste(out$output, collapse = "\n"))
})
