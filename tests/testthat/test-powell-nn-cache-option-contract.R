test_that("np.objective.cache is included in MPI option sync", {
  expect_true("np.objective.cache" %in% npRmpi:::.npRmpi_autodispatch_option_keys())
})

test_that("np.objective.cache default ignores legacy environment off switch", {
  old <- Sys.getenv("NP_NN_POWELL_CACHE_INSTRUMENT", unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("NP_NN_POWELL_CACHE_INSTRUMENT")
    } else {
      Sys.setenv(NP_NN_POWELL_CACHE_INSTRUMENT = old)
    }
  }, add = TRUE)

  Sys.setenv(NP_NN_POWELL_CACHE_INSTRUMENT = "off")
  expect_true(npRmpi:::npObjectiveCacheEnabled())

  Sys.setenv(NP_NN_POWELL_CACHE_INSTRUMENT = "on")
  expect_true(npRmpi:::npObjectiveCacheEnabled())
})

test_that("NN objective cache rejects non-finite values", {
  cache <- npRmpi:::.np_r_nn_cache_new(TRUE, key.length = 1L)

  miss <- npRmpi:::.np_r_nn_cache_get(cache, 1L)
  expect_false(npRmpi:::.np_r_nn_cache_put(cache, miss$token, Inf))
  expect_false(npRmpi:::.np_r_nn_cache_get(cache, 1L)$hit)

  miss <- npRmpi:::.np_r_nn_cache_get(cache, 1L)
  expect_false(npRmpi:::.np_r_nn_cache_put(cache, miss$token, NaN))
  expect_false(npRmpi:::.np_r_nn_cache_get(cache, 1L)$hit)

  miss <- npRmpi:::.np_r_nn_cache_get(cache, 1L)
  expect_false(npRmpi:::.np_r_nn_cache_put(cache, miss$token, numeric()))
  expect_false(npRmpi:::.np_r_nn_cache_get(cache, 1L)$hit)

  miss <- npRmpi:::.np_r_nn_cache_get(cache, 1L)
  expect_true(npRmpi:::.np_r_nn_cache_put(cache, miss$token, 3.25))
  hit <- npRmpi:::.np_r_nn_cache_get(cache, 1L)
  expect_true(hit$hit)
  expect_equal(hit$value, 3.25)
})

test_that("C NN objective cache guards non-finite values before insertion", {
  path <- testthat::test_path("..", "..", "src", "np.c")
  skip_if_not(file.exists(path), "package source file unavailable")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "static void bwm_nn_cache_put")
  expect_match(text, "R_FINITE\\(value\\)")
})

test_that("np.objective.cache controls continuous NN Powell caching under MPI", {
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
      "  old <- options(np.messages = FALSE, np.tree = FALSE, np.objective.cache = cache)",
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

test_that("np.objective.cache controls npscoef continuous NN R optimizer caching under MPI", {
  env <- npRmpi_subprocess_env("NP_RMPI_NO_REUSE_SLAVES=1")
  skip_if(is.null(env))

  out <- npRmpi_run_rscript_subprocess(
    c(
      "library(npRmpi)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "run_bw <- function(bwtype, cache) {",
      "  set.seed(20260523)",
      "  n <- 80L",
      "  z1 <- runif(n)",
      "  z2 <- runif(n)",
      "  x <- rnorm(n)",
      "  y <- 1 + (0.5 + sin(2 * pi * z1)) * x + 0.25 * cos(2 * pi * z2) + rnorm(n, sd = 0.2)",
      "  old <- options(np.messages = FALSE, np.objective.cache = cache)",
      "  on.exit(options(old), add = TRUE)",
      "  npscoefbw(xdat = data.frame(x = x), ydat = y, zdat = data.frame(z1 = z1, z2 = z2), regtype = 'lc', bwtype = bwtype, ckertype = 'epanechnikov', nmulti = 1L, optim.maxit = 35L, optim.maxattempts = 1L)",
      "}",
      "for (bwtype in c('generalized_nn', 'adaptive_nn')) {",
      "  cached <- run_bw(bwtype, TRUE)",
      "  uncached <- run_bw(bwtype, FALSE)",
      "  stopifnot(isTRUE(all.equal(cached$bw, uncached$bw, tolerance = 0)))",
      "  stopifnot(isTRUE(all.equal(cached$fval, uncached$fval, tolerance = 0)))",
      "  stopifnot(identical(as.numeric(cached$num.feval[1L]), as.numeric(uncached$num.feval[1L])))",
      "  stopifnot(identical(unname(cached$nn.cache[['enabled']]), 1))",
      "  stopifnot(unname(cached$nn.cache[['hits']]) > 0)",
      "  stopifnot(as.numeric(cached$num.feval.fast[1L]) >= unname(cached$nn.cache[['hits']]))",
      "  stopifnot(as.numeric(cached$num.feval.fast[1L]) <= as.numeric(cached$num.feval[1L]))",
      "  stopifnot(identical(unname(uncached$nn.cache[['enabled']]), 0))",
      "  stopifnot(identical(unname(uncached$nn.cache[['hits']]), 0))",
      "}"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(out$status, 0L, info = paste(out$output, collapse = "\n"))
})

test_that("np.objective.cache controls npindex continuous NN R optimizer caching under MPI", {
  env <- npRmpi_subprocess_env("NP_RMPI_NO_REUSE_SLAVES=1")
  skip_if(is.null(env))

  out <- npRmpi_run_rscript_subprocess(
    c(
      "library(npRmpi)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "run_bw <- function(method, bwtype, cache) {",
      "  set.seed(123 + match(bwtype, c('generalized_nn', 'adaptive_nn')) + if (identical(method, 'kleinspady')) 10L else 0L)",
      "  n <- 70L",
      "  xdat <- data.frame(x1 = runif(n), x2 = runif(n))",
      "  eta <- (xdat$x1 + xdat$x2) / 2",
      "  ydat <- if (identical(method, 'kleinspady')) as.integer(runif(n) < plogis(2 * eta - 1)) else eta + rnorm(n, sd = 0.25)",
      "  old <- options(np.messages = FALSE, np.objective.cache = cache)",
      "  on.exit(options(old), add = TRUE)",
      "  npindexbw(xdat = xdat, ydat = ydat, method = method, regtype = 'lc', bwtype = bwtype, nmulti = 1L, optim.maxit = 35L, optim.maxattempts = 1L)",
      "}",
      "for (method in c('ichimura', 'kleinspady')) {",
      "  for (bwtype in c('generalized_nn', 'adaptive_nn')) {",
      "    cached <- run_bw(method, bwtype, TRUE)",
      "    uncached <- run_bw(method, bwtype, FALSE)",
      "    stopifnot(isTRUE(all.equal(cached$beta, uncached$beta, tolerance = 0)))",
      "    stopifnot(isTRUE(all.equal(cached$bw, uncached$bw, tolerance = 0)))",
      "    stopifnot(isTRUE(all.equal(cached$fval, uncached$fval, tolerance = 0)))",
      "    stopifnot(identical(as.numeric(cached$num.feval[1L]), as.numeric(uncached$num.feval[1L])))",
      "    stopifnot(identical(unname(cached$nn.cache[['enabled']]), 1))",
      "    hits <- unname(cached$nn.cache[['hits']])",
      "    stopifnot(is.numeric(hits), length(hits) == 1L, is.finite(hits), hits >= 0)",
      "    stopifnot(as.numeric(cached$num.feval.fast[1L]) >= hits)",
      "    stopifnot(as.numeric(cached$num.feval.fast[1L]) <= as.numeric(cached$num.feval[1L]))",
      "    stopifnot(identical(unname(uncached$nn.cache[['enabled']]), 0))",
      "    stopifnot(identical(unname(uncached$nn.cache[['hits']]), 0))",
      "  }",
      "}"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(out$status, 0L, info = paste(out$output, collapse = "\n"))
})
