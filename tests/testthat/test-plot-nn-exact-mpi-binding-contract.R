test_that("adaptive conditional exact bootstrap plot binds only required worker state", {
  skip_on_cran()

  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
      "options(np.messages=FALSE, np.tree=FALSE)",
      "set.seed(2468)",
      "x <- data.frame(x = runif(80))",
      "y <- data.frame(y = rnorm(80))",
      "ex <- data.frame(x = seq(0, 1, length.out = 15))",
      "ey <- data.frame(y = seq(min(y$y), max(y$y), length.out = 15))",
      "bw <- npcdensbw(xdat=x, ydat=y, exdat=ex, eydat=ey, bws=c(7, 7), bwtype='adaptive_nn', bandwidth.compute=FALSE)",
      "out <- plot(",
      "  bw,",
      "  xdat=x,",
      "  ydat=y,",
      "  plot.errors.method='bootstrap',",
      "  plot.errors.boot.method='inid',",
      "  plot.errors.boot.num=5,",
      "  plot.errors.boot.nonfixed='exact',",
      "  plot.behavior='data'",
      ")",
      "stopifnot(is.list(out), length(out) > 0)",
      "cat('ADAPTIVE_COND_EXACT_SESSION_OK\\n')"
    ),
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(
    any(grepl("ADAPTIVE_COND_EXACT_SESSION_OK", res$output, fixed = TRUE)),
    info = paste(res$output, collapse = "\n")
  )
})
