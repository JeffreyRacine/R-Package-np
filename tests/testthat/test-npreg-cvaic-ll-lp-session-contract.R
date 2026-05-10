test_that("session npreg cv.aic preserves multivariate ll/lp(degree=1) objective parity", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "options(np.messages = FALSE, np.tree = FALSE)",
      "set.seed(20260509)",
      "n <- 90L",
      "x1 <- runif(n, -1, 1)",
      "x2 <- rnorm(n)",
      "xdat <- data.frame(x1 = x1, x2 = x2)",
      "ydat <- sin(2 * pi * x1) + 0.5 * x2 + rnorm(n, sd = 0.15)",
      "set.seed(90210)",
      "bw.ll <- npregbw(xdat = xdat, ydat = ydat, regtype = 'll', bwmethod = 'cv.aic', nmulti = 1L)",
      "set.seed(90210)",
      "bw.lp <- npregbw(xdat = xdat, ydat = ydat, regtype = 'lp', basis = 'glp', degree = rep.int(1L, ncol(xdat)), bwmethod = 'cv.aic', nmulti = 1L)",
      "stopifnot(isTRUE(all.equal(as.numeric(bw.ll$fval), as.numeric(bw.lp$fval), tolerance = 1e-10)))",
      "stopifnot(isTRUE(all.equal(as.numeric(bw.ll$bw), as.numeric(bw.lp$bw), tolerance = 1e-9)))",
      "cat('NPREG_CVAIC_LL_LP1_SESSION_PARITY_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPREG_CVAIC_LL_LP1_SESSION_PARITY_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
