test_that("session-route qreg plots preserve quantreg mode and fitted tau", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(), silent = TRUE), add = TRUE)",
      "set.seed(108)",
      "n <- 20L",
      "x <- runif(n)",
      "y <- sin(2 * pi * x) + rnorm(n, sd = 0.05)",
      "xdat <- data.frame(x = x)",
      "ydat <- data.frame(y = y)",
      "bw <- npcdistbw(xdat = xdat, ydat = ydat, bws = c(0.25, 0.25), bandwidth.compute = FALSE)",
      "fit <- npqreg(bws = bw, txdat = xdat, tydat = ydat, tau = 0.4)",
      "bw.fixed <- plot(bw, xdat = xdat, ydat = ydat, output = 'data', perspective = TRUE, view = 'fixed', quantreg = TRUE, tau = 0.4, neval = 6)",
      "fit.fixed <- plot(fit, output = 'data', perspective = TRUE, view = 'fixed', neval = 6)",
      "bw.slice <- plot(bw, xdat = xdat, ydat = ydat, output = 'data', perspective = FALSE, quantreg = TRUE, tau = 0.4, neval = 6)",
      "fit.slice <- plot(fit, output = 'data', perspective = FALSE, neval = 6)",
      "stopifnot(inherits(fit.fixed, 'np_plot_data') && nrow(fit.fixed$targets) == 1L)",
      "stopifnot(inherits(fit.slice, 'np_plot_data') && nrow(fit.slice$targets) == 1L)",
      "stopifnot(identical(fit.fixed$targets$tau, 0.4))",
      "stopifnot(identical(fit.slice$targets$tau, 0.4))",
      "stopifnot(max(abs(fit.fixed$panels[[1L]]$values$estimate - bw.fixed$panels[[1L]]$values$estimate)) <= 1e-5)",
      "stopifnot(isTRUE(all.equal(fit.fixed$panels[[1L]]$grid, bw.fixed$panels[[1L]]$grid)))",
      "stopifnot(max(abs(fit.slice$panels[[1L]]$values$estimate - bw.slice$panels[[1L]]$values$estimate)) <= 1e-5)",
      "stopifnot(isTRUE(all.equal(fit.slice$panels[[1L]]$grid, bw.slice$panels[[1L]]$grid)))",
      "cat('QREG_PLOT_TAU_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("QREG_PLOT_TAU_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
