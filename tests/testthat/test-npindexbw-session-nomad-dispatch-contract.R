test_that("session npindex automatic NOMAD degree search dispatches bandwidth work", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves = 2L, quiet = TRUE)",
      "options(np.messages = FALSE, np.tree = FALSE)",
      "assign('.seen', character(), envir = .GlobalEnv)",
      "trace(npRmpi:::.npRmpi_autodispatch_call,",
      "      tracer = quote({",
      "        head <- mc[[1L]]",
      "        nm <- if (is.symbol(head)) as.character(head) else paste(deparse(head), collapse = ' ')",
      "        .GlobalEnv$.seen <- c(.GlobalEnv$.seen, nm)",
      "      }),",
      "      print = FALSE)",
      "on.exit(try(untrace(npRmpi:::.npRmpi_autodispatch_call), silent = TRUE), add = TRUE)",
      "set.seed(20260504)",
      "n <- 20L",
      "x1 <- runif(n, -1, 1)",
      "x2 <- runif(n, -1, 1)",
      "y <- sin(x1 + 0.5 * x2) + rnorm(n, sd = 0.1)",
      "fit <- npindex(y ~ x1 + x2, nomad = TRUE, degree.max = 1L, nmulti = 1L)",
      "stopifnot(inherits(fit, 'singleindex'))",
      "stopifnot('npindexbw.NULL' %in% .GlobalEnv$.seen)",
      "stopifnot('npindex.sibandwidth' %in% .GlobalEnv$.seen)",
      "cat('NPINDEX_SESSION_NOMAD_AUTODISPATCH_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPINDEX_SESSION_NOMAD_AUTODISPATCH_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindex Klein-Spady NOMAD dispatches bandwidth work with one slave", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "options(np.messages = FALSE, np.tree = FALSE)",
      "assign('.seen', character(), envir = .GlobalEnv)",
      "trace(npRmpi:::.npRmpi_autodispatch_call,",
      "      tracer = quote({",
      "        head <- mc[[1L]]",
      "        nm <- if (is.symbol(head)) as.character(head) else paste(deparse(head), collapse = ' ')",
      "        .GlobalEnv$.seen <- c(.GlobalEnv$.seen, nm)",
      "      }),",
      "      print = FALSE)",
      "on.exit(try(untrace(npRmpi:::.npRmpi_autodispatch_call), silent = TRUE), add = TRUE)",
      "set.seed(20260509)",
      "n <- 30L",
      "x1 <- runif(n, -1, 1)",
      "x2 <- runif(n, -1, 1)",
      "eta <- x1 - 0.75 * x2",
      "prob <- 1 / (1 + exp(-eta))",
      "y <- rbinom(n, 1L, prob)",
      "fit <- npindex(y ~ x1 + x2, method = 'kleinspady', nomad = TRUE, degree.max = 1L, nmulti = 1L)",
      "stopifnot(inherits(fit, 'singleindex'))",
      "stopifnot('npindexbw.NULL' %in% .GlobalEnv$.seen)",
      "stopifnot('npindex.sibandwidth' %in% .GlobalEnv$.seen)",
      "cat('NPINDEX_SESSION_KLEINSPADY_NOMAD_AUTODISPATCH_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPINDEX_SESSION_KLEINSPADY_NOMAD_AUTODISPATCH_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
