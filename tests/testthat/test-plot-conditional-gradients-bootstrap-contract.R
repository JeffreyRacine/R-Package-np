test_that("session-route conditional density/distribution gradient bootstrap inid works for bw and fit objects", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.autodispatch = FALSE, np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(20260226)",
      "n <- 24L",
      "x1 <- factor(sample(c('a', 'b'), n, replace = TRUE))",
      "x2 <- rnorm(n)",
      "y <- rnorm(n)",
      "xdat <- data.frame(x1 = x1, x2 = x2)",
      "ydat <- data.frame(y = y)",
      "bw.cd <- suppressWarnings(npcdensbw(xdat = xdat, ydat = ydat, nmulti = 1L))",
      "bw.cdist <- suppressWarnings(npcdistbw(xdat = xdat, ydat = ydat, nmulti = 1L))",
      "fit.cd <- npcdens(bws = bw.cd)",
      "fit.cdist <- npcdist(bws = bw.cdist)",
      "for (obj in list(bw.cd, bw.cdist)) {",
      "  out <- suppressWarnings(plot(obj, xdat = xdat, ydat = ydat, plot.behavior = 'data', perspective = FALSE, gradients = TRUE, plot.errors.method = 'bootstrap', plot.errors.boot.method = 'inid', plot.errors.boot.num = 5L))",
      "  stopifnot(is.list(out), length(out[[1L]]$bxp) > 0L)",
      "}",
      "for (obj in list(fit.cd, fit.cdist)) {",
      "  out <- suppressWarnings(plot(obj, xdat = xdat, ydat = ydat, plot.behavior = 'data', perspective = FALSE, gradients = TRUE, plot.errors.method = 'bootstrap', plot.errors.boot.method = 'inid', plot.errors.boot.num = 5L))",
      "  stopifnot(is.list(out), length(out[[1L]]$bxp) > 0L)",
      "}",
      "cat('COND_INID_CONTRACT_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("COND_INID_CONTRACT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
