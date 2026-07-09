test_that("public npindex formula NOMAD shortcut completes in session mode", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "options(np.messages = FALSE)",
      "set.seed(42)",
      "n <- 50L",
      "dat <- data.frame(",
      "  x1 = runif(n, min = -1, max = 1),",
      "  x2 = runif(n, min = -1, max = 1)",
      ")",
      "dat$y <- dat$x1 - dat$x2 + rnorm(n)",
      "fit <- npindex(",
      "  y ~ x1 + x2,",
      "  data = dat,",
      "  nomad = TRUE,",
      "  degree.max = 1L,",
      "  nmulti = 1L",
      ")",
      "stopifnot(inherits(fit, 'singleindex'))",
      "stopifnot(inherits(fit$bws, 'sibandwidth'))",
      "stopifnot(is.finite(as.numeric(fit$bws$fval)))",
      "stopifnot(is.finite(as.numeric(fit$bws$bw)))",
      "cat('NPINDEX_PUBLIC_NOMAD_SMOKE_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPINDEX_PUBLIC_NOMAD_SMOKE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
