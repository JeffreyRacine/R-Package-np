test_that("session-route wild categorical regression gradient helper matches explicit refits", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.autodispatch = FALSE, np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "wild_helper <- getFromNamespace('.np_plot_boot_from_hat_wild_factor_effects', 'npRmpi')",
      "base_helper <- getFromNamespace('.np_plot_boot_from_hat_wild', 'npRmpi')",
      "set.seed(20260312)",
      "n <- 36L",
      "g <- factor(sample(c('a', 'b'), n, replace = TRUE))",
      "x <- runif(n)",
      "y <- 1 + 0.5 * (g == 'b') + sin(2 * pi * x) + rnorm(n, sd = 0.05)",
      "xdat <- data.frame(g = g, x = x)",
      "bw <- npregbw(xdat = xdat, ydat = y, regtype = 'll', bwtype = 'fixed', bws = c(0.25, 0.3), bandwidth.compute = FALSE)",
      "exdat <- plot(bw, xdat = xdat, ydat = y, gradients = TRUE, plot.behavior = 'data', perspective = FALSE)[[1L]]$eval",
      "H <- npreghat(bws = bw, txdat = xdat, exdat = exdat, output = 'matrix')",
      "fit.mean <- as.vector(npreghat(bws = bw, txdat = xdat, exdat = xdat, y = y, output = 'apply'))",
      "B <- 7L",
      "set.seed(11)",
      "base.out <- base_helper(H = H, ydat = y, fit.mean = fit.mean, B = B, wild = 'rademacher')",
      "set.seed(11)",
      "helper.out <- wild_helper(H = H, ydat = y, fit.mean = fit.mean, B = B, wild = 'rademacher')",
      "fit0 <- npreg(txdat = xdat, tydat = y, exdat = exdat, bws = bw, gradients = FALSE, warn.glp.gradient = FALSE)$mean",
      "stopifnot(isTRUE(all.equal(helper.out$t, sweep(base.out$t, 1L, base.out$t[, 1L], '-'), tolerance = 1e-10)))",
      "stopifnot(isTRUE(all.equal(helper.out$t0, base.out$t0 - base.out$t0[1L], tolerance = 1e-10)))",
      "stopifnot(isTRUE(all.equal(as.vector(helper.out$t0), as.vector(fit0 - fit0[1L]), tolerance = 1e-6)))",
      "cat('RBAND_CATGRAD_HELPER_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("RBAND_CATGRAD_HELPER_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session-route categorical regression gradient bootstrap works for default and wild routes", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.autodispatch = FALSE, np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(20260312)",
      "n <- 36L",
      "g <- factor(sample(c('a', 'b'), n, replace = TRUE))",
      "x <- runif(n)",
      "y <- 1 + 0.4 * (g == 'b') + cos(2 * pi * x) + rnorm(n, sd = 0.05)",
      "xdat <- data.frame(g = g, x = x)",
      "bw <- npregbw(xdat = xdat, ydat = y, regtype = 'll', bwtype = 'fixed', bws = c(0.25, 0.3), bandwidth.compute = FALSE)",
      "fit <- npreg(bws = bw, txdat = xdat, tydat = y, gradients = TRUE, errors = TRUE)",
      "for (boot.method in c('default', 'wild')) {",
      "  args <- list(fit, xdat = xdat, ydat = y, plot.behavior = 'data', perspective = FALSE, gradients = TRUE, plot.errors.method = 'bootstrap', plot.errors.boot.num = 11L)",
      "  if (!identical(boot.method, 'default')) args$plot.errors.boot.method <- boot.method",
      "  out <- suppressWarnings(do.call(plot, args))",
      "  stopifnot(is.list(out), length(out) >= 1L, length(out[[1L]]$bxp) > 0L, length(out[[1L]]$bxp$names) == 2L)",
      "  stopifnot(all(is.finite(out[[1L]]$grad)), all(is.finite(out[[1L]]$glerr)), all(is.finite(out[[1L]]$gherr)))",
      "}",
      "cat('RBAND_CATGRAD_CONSUMER_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("RBAND_CATGRAD_CONSUMER_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
