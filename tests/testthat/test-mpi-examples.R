test_that("MPI examples from man pages work in an isolated session", {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess MPI examples")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.autodispatch = TRUE, np.messages = FALSE)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "npseed(712)",
      "n <- 50L",
      "x <- runif(n)",
      "y <- x + rnorm(n, sd = 0.1)",
      "mydat <- data.frame(x = x, y = y)",
      "bw <- npregbw(y ~ x, data = mydat)",
      "stopifnot(inherits(bw, 'rbandwidth'))",
      "data(Engel95)",
      "Engel_sub <- Engel95[1:100, ]",
      "Engel_sub <- Engel_sub[order(Engel_sub$logexp), ]",
      "model.iv <- with(Engel_sub, npregiv(",
      "  y = food,",
      "  z = logexp,",
      "  w = logwages,",
      "  method = 'Landweber-Fridman',",
      "  iterate.max = 2L",
      "))",
      "stopifnot(inherits(model.iv, 'npregiv'))",
      "data('USArrests')",
      "dat <- USArrests[1:20, ]",
      "pair_list <- np.pairs(y_vars = c('Murder', 'UrbanPop'), y_dat = dat)",
      "stopifnot(is.list(pair_list))",
      "stopifnot(identical(length(pair_list$pair_kerns), 4L))",
      "cat('MPI_EXAMPLES_SUBPROCESS_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  info <- paste(res$output, collapse = "\n")
  expect_true(res$status %in% c(0L, 137L), info = info)
  expect_true(any(grepl("MPI_EXAMPLES_SUBPROCESS_OK", res$output, fixed = TRUE)),
              info = info)
})
