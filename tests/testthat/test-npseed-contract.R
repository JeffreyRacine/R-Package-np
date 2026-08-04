test_that("npseed accepts valid integer-representable scalars", {
  env <- npRmpi_subprocess_env(c(
    "_R_CHECK_PACKAGE_NAME_=",
    "NP_RMPI_TEST_SUITE_POOL="
  ))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess seed contract")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npseed(42)",
      "npseed(-42)",
      "npseed(42.0)",
      "cat('NPSEED_VALID_SCALARS_OK\\n')"
    ),
    timeout = 30L,
    env = env
  )

  info <- paste(res$output, collapse = "\n")
  expect_identical(res$status, 0L, info = info)
  expect_true(any(grepl("NPSEED_VALID_SCALARS_OK", res$output, fixed = TRUE)),
              info = info)
})

test_that("npseed rejects invalid scalar inputs before reaching C", {
  expect_error(npseed(NA_integer_), "'seed' must be a single finite numeric value")
  expect_error(npseed(Inf), "'seed' must be a single finite numeric value")
  expect_error(npseed(NaN), "'seed' must be a single finite numeric value")
  expect_error(npseed(numeric()), "'seed' must be a single finite numeric value")
  expect_error(npseed(c(1, 2)), "'seed' must be a single finite numeric value")
  expect_error(npseed(1.5), "'seed' must be representable as a non-negative integer after abs\\(\\)")
  expect_error(npseed(.Machine$integer.max + 1), "'seed' must be representable as a non-negative integer after abs\\(\\)")
  expect_error(npseed("1"), "'seed' must be a single finite numeric value")
})

test_that("C_np_set_seed rejects direct misuse", {
  expect_error(.Call("C_np_set_seed", NA_integer_, PACKAGE = "npRmpi"), "seed must be finite")
  expect_error(.Call("C_np_set_seed", Inf, PACKAGE = "npRmpi"), "seed must be finite")
  expect_error(.Call("C_np_set_seed", 1.5, PACKAGE = "npRmpi"), "seed must be representable as a non-negative integer after abs\\(\\)")
  expect_error(.Call("C_np_set_seed", c(1L, 2L), PACKAGE = "npRmpi"), "seed must have length 1")
  expect_error(.Call("C_np_set_seed", TRUE, PACKAGE = "npRmpi"), "seed must be numeric")
})

test_that("npseed fails closed for inconsistent active-pool state", {
  old.active <- getOption("npRmpi.pool.active", FALSE)
  on.exit(options(npRmpi.pool.active = old.active), add = TRUE)

  options(npRmpi.pool.active = TRUE)
  set.seed(271828)
  before <- .Random.seed

  expect_error(
    npseed(19),
    "active MPI pool is inconsistent",
    fixed = TRUE
  )
  expect_identical(.Random.seed, before)
})

test_that("npseed synchronizes the C backend across an active MPI pool", {
  env <- npRmpi_subprocess_env(c(
    "_R_CHECK_PACKAGE_NAME_=",
    "NP_RMPI_TEST_SUITE_POOL="
  ))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess MPI seed contract")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.autodispatch = TRUE, np.messages = FALSE)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(8675309)",
      "r.before <- .Random.seed",
      "seq.get <- getFromNamespace('.npRmpi_spmd_seq_get', 'npRmpi')",
      "seq.before <- seq.get()",
      "bad <- try(npseed(NA_real_), silent = TRUE)",
      "stopifnot(inherits(bad, 'try-error'))",
      "stopifnot(grepl('single finite numeric', as.character(bad), fixed = TRUE))",
      "stopifnot(identical(.Random.seed, r.before))",
      "stopifnot(identical(seq.get(), seq.before))",
      "n <- 36L",
      "x <- seq(0.03, 0.97, length.out = n)",
      "u <- factor(rep(c('a', 'b', 'c'), length.out = n))",
      "y <- sin(2 * pi * x) + c(a = 0, b = 0.25, c = -0.15)[u]",
      "dat <- data.frame(y, x, u)",
      "run <- function() npcdistbw(",
      "  y ~ x + u, data = dat, scale.factor.init = 0.55, nmulti = 2L,",
      "  transform.bounds = TRUE, powell.remin = FALSE, itmax = 1L",
      ")",
      "npseed(0)",
      "stopifnot(identical(seq.get(), seq.before))",
      "first <- run()",
      "seq.after.first <- seq.get()",
      "npseed(0)",
      "stopifnot(identical(seq.get(), seq.after.first))",
      "second <- run()",
      "stopifnot(is.finite(first$fval))",
      "stopifnot(is.finite(second$fval))",
      "stopifnot(length(first$fval.history) == 2L)",
      "stopifnot(length(second$fval.history) == 2L)",
      "stopifnot(all(is.finite(first$fval.history)))",
      "stopifnot(all(is.finite(second$fval.history)))",
      "cat('NPSEED_ACTIVE_POOL_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  info <- paste(res$output, collapse = "\n")
  expect_true(res$status %in% c(0L, 137L), info = info)
  expect_true(any(grepl("NPSEED_ACTIVE_POOL_OK", res$output, fixed = TRUE)),
              info = info)
})
