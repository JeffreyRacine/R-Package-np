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
      "bw.fixed <- plot(bw, xdat = xdat, ydat = ydat, plot.behavior = 'data', perspective = TRUE, view = 'fixed', quantreg = TRUE, tau = 0.4, neval = 6)",
      "fit.fixed <- plot(fit, plot.behavior = 'data', perspective = TRUE, view = 'fixed', neval = 6)",
      "bw.slice <- plot(bw, xdat = xdat, ydat = ydat, plot.behavior = 'data', perspective = FALSE, quantreg = TRUE, tau = 0.4, neval = 6)",
      "fit.slice <- plot(fit, plot.behavior = 'data', perspective = FALSE, neval = 6)",
      "stopifnot(all(vapply(fit.fixed, inherits, logical(1), 'qregression')))",
      "stopifnot(all(vapply(fit.slice, inherits, logical(1), 'qregression')))",
      "stopifnot(all(vapply(fit.fixed, function(xi) identical(xi$tau, 0.4), logical(1))))",
      "stopifnot(all(vapply(fit.slice, function(xi) identical(xi$tau, 0.4), logical(1))))",
      "stopifnot(max(abs(fit.fixed[[1L]]$quantile - bw.fixed[[1L]]$quantile)) <= 1e-5)",
      "stopifnot(isTRUE(all.equal(fit.fixed[[1L]]$xeval, bw.fixed[[1L]]$xeval)))",
      "stopifnot(max(abs(fit.slice[[1L]]$quantile - bw.slice[[1L]]$quantile)) <= 1e-5)",
      "stopifnot(isTRUE(all.equal(fit.slice[[1L]]$xeval, bw.slice[[1L]]$xeval)))",
      "cat('QREG_PLOT_TAU_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("QREG_PLOT_TAU_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
