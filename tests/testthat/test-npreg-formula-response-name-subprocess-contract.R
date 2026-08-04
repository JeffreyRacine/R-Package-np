local_npRmpi_formula_subprocess_env <- function(extra = character()) {
  npRmpi_subprocess_env(extra)
}

test_that("npreg formula fits preserve the response name for plotting metadata in session mode", {
  skip_on_cran()
  env <- local_npRmpi_formula_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(42)",
      "n <- 60L",
      "dat <- data.frame(x = runif(n, -1, 1))",
      "dat$y <- dat$x^2 + rnorm(n, sd = 0.25 * stats::sd(dat$x))",
      "fit <- npreg(y ~ x, data = dat, nmulti = 1, regtype = 'lp', degree = 1, bwtype = 'adaptive_nn')",
      "stopifnot(identical(fit$bws$ynames, 'y'))",
      "png(tempfile(fileext = '.png'), width = 800, height = 600)",
      "plot(fit, errors = 'asymptotic', view = 'fixed')",
      "dev.off()",
      "cat('NPREG_FORMULA_YNAME_OK\\n')"
    ),
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPREG_FORMULA_YNAME_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
