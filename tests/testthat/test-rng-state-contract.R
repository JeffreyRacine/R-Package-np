test_that("categorical padding helper is deterministic and RNG-pristine", {
  num_not_in <- getFromNamespace("numNotIn", "npRmpi")

  set.seed(8675309)
  before <- .Random.seed

  expect_identical(num_not_in(numeric()), 0)
  expect_identical(num_not_in(c(1, 2, 3)), 0)
  expect_false(num_not_in(c(0, 1, 2)) %in% c(0, 1, 2))
  expect_false(num_not_in(c(-3, -2, -1, 0, 1, 2, 3)) %in% c(-3, -2, -1, 0, 1, 2, 3))
  expect_true(is.finite(num_not_in(c(-.Machine$double.xmax, 0, .Machine$double.xmax))))

  expect_identical(.Random.seed, before)
})

test_that("npseed rejects malformed seeds without touching R RNG", {
  npseed <- getFromNamespace("npseed", "npRmpi")

  set.seed(314159)
  before <- .Random.seed

  expect_error(npseed(NA_real_), "single finite numeric")
  expect_error(npseed(1.5), "non-negative integer")
  expect_error(npseed(.Machine$integer.max + 1), "non-negative integer")

  expect_identical(.Random.seed, before)
})

test_that("valid npseed leaves the R RNG unchanged in an isolated process", {
  env <- npRmpi_subprocess_env(c(
    "_R_CHECK_PACKAGE_NAME_=",
    "NP_RMPI_TEST_SUITE_POOL="
  ))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess RNG contract")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "set.seed(314159)",
      "before <- .Random.seed",
      "npseed(0)",
      "stopifnot(identical(.Random.seed, before))",
      "cat('NPSEED_R_RNG_PRISTINE_OK\\n')"
    ),
    timeout = 30L,
    env = env
  )

  info <- paste(res$output, collapse = "\n")
  expect_identical(res$status, 0L, info = info)
  expect_true(any(grepl("NPSEED_R_RNG_PRISTINE_OK", res$output, fixed = TRUE)),
              info = info)
})
