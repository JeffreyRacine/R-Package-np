test_that("session adaptive conditional plot-data smoke completes in subprocess", {
  skip_on_cran()
  marker <- "SESSION_ADAPTIVE_CONDITIONAL_PLOT_OK"
  res <- npRmpi_run_isolated_contract(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "set.seed(42)",
      "n <- 36L",
      "x <- runif(n)",
      "y <- x^2 + rnorm(n, sd=0.1)",
      "d <- data.frame(x=x, y=y)",
      "bw.cd <- npcdensbw(y~x, data=d, bws=c(7, 7), bandwidth.compute=FALSE, regtype='lp', degree=1L, bwtype='adaptive_nn')",
      "bw.cf <- npcdistbw(y~x, data=d, bws=c(7, 7), bandwidth.compute=FALSE, regtype='lp', degree=1L, bwtype='adaptive_nn')",
      "fit.cd <- npcdens(bws=bw.cd, txdat=d['x'], tydat=d['y'], gradients=FALSE)",
      "fit.cf <- npcdist(bws=bw.cf, txdat=d['x'], tydat=d['y'], gradients=FALSE)",
      "png(tempfile(fileext='.png'))",
      "on.exit(dev.off(), add=TRUE)",
      "out.cd <- suppressWarnings(plot(fit.cd, output = 'data'))",
      "out.cf <- suppressWarnings(plot(fit.cf, output = 'data'))",
      "stopifnot(is.list(out.cd), length(out.cd) > 0L)",
      "stopifnot(is.list(out.cf), length(out.cf) > 0L)",
      "cat('SESSION_ADAPTIVE_CONDITIONAL_PLOT_OK\\n')"
    ),
    marker = marker,
    timeout = 20L
  )
  skip_if(is.null(res), "installed npRmpi unavailable for subprocess smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(res$witnessed, info = paste(res$output, collapse = "\n"))
})
