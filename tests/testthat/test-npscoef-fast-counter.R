test_that("npscoefbw fast counter fires only on the z-irrelevant sentinel", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "old_opts <- options(np.messages = FALSE, np.largeh.rel.tol = 0.05, np.disc.upper.rel.tol = 0.05)",
      "on.exit(options(old_opts), add = TRUE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(), silent = TRUE), add = TRUE)",
      "set.seed(77)",
      "n <- 120L",
      "x <- runif(n)",
      "z <- runif(n)",
      "dat_lc <- data.frame(y = rnorm(n, sd = 0.5 * sd(x)), x = x, z = z)",
      "dat_ll <- data.frame(y = x + rnorm(n, sd = 0.5 * sd(x)), x = x, z = z)",
      "bw_lc <- npscoefbw(y ~ x | z, data = dat_lc, regtype = 'lc', bwmethod = 'cv.ls', nmulti = 1)",
      "bw_ll <- npscoefbw(y ~ x | z, data = dat_ll, regtype = 'll', bwmethod = 'cv.ls', nmulti = 1)",
      "bw_lp <- npscoefbw(y ~ x | z, data = dat_ll, regtype = 'lp', degree = 1, bwmethod = 'cv.ls', nmulti = 1)",
      "stopifnot(is.finite(as.numeric(bw_lc$num.feval.fast)), is.finite(as.numeric(bw_ll$num.feval.fast)), is.finite(as.numeric(bw_lp$num.feval.fast)))",
      "stopifnot(as.integer(bw_lc$num.feval.fast) > 0L, as.integer(bw_ll$num.feval.fast) > 0L, as.integer(bw_lp$num.feval.fast) > 0L)",
      "stopifnot(isTRUE(all.equal(as.numeric(bw_ll$fval), as.numeric(bw_lp$fval), tolerance = 1e-10)))",
      "stopifnot(identical(as.integer(bw_ll$num.feval.fast), as.integer(bw_lp$num.feval.fast)))",
      "cat('NPSCOEF_FAST_COUNTER_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPSCOEF_FAST_COUNTER_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
