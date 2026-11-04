test_that("session-route wild categorical regression gradient helper matches explicit refits", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      plot_test_support(),
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
      "exdat <- plot_coordinates(plot(bw, xdat = xdat, ydat = y, gradients = TRUE, output = 'data', perspective = FALSE), role='X')",
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

test_that("session-route ordered regression gradient wild bootstrap targets adjacent differences", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      plot_test_support(),
      "options(npRmpi.autodispatch = FALSE, np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "wild_helper <- getFromNamespace('.np_plot_boot_from_hat_wild', 'npRmpi')",
      "ordered_gradient_matrix <- function(x) { dx <- x[, -1L, drop = FALSE] - x[, -ncol(x), drop = FALSE]; cbind(dx[, 1L], dx) }",
      "ordered_gradient_vector <- function(x) { dx <- diff(x); c(dx[1L], dx) }",
      "set.seed(20260709)",
      "n <- 60L",
      "o <- ordered(rep(seq_len(5L), length.out = n), levels = seq_len(5L))",
      "y <- as.numeric(o) + rnorm(n, sd = 0.08)",
      "xdat <- data.frame(o = o)",
      "bw <- npregbw(xdat = xdat, ydat = y, regtype = 'lc', bwtype = 'fixed', bws = 0.25, bandwidth.compute = FALSE)",
      "fit <- npreg(bws = bw, txdat = xdat, tydat = y, gradients = TRUE, se = TRUE)",
      "base <- plot(fit, xdat = xdat, ydat = y, output = 'data', perspective = FALSE, gradients = TRUE)",
      "exdat <- plot_coordinates(base, role='X')",
      "H <- npreghat(bws = bw, txdat = xdat, exdat = exdat, output = 'matrix')",
      "fit.mean <- as.vector(npreghat(bws = bw, txdat = xdat, exdat = xdat, y = y, output = 'apply'))",
      "B <- 31L",
      "set.seed(33)",
      "expected.level <- wild_helper(H = H, ydat = y, fit.mean = fit.mean, B = B, wild = 'rademacher')",
      "expected.t <- ordered_gradient_matrix(expected.level$t)",
      "expected.t0 <- ordered_gradient_vector(expected.level$t0)",
      "set.seed(33)",
      "out <- plot(fit, xdat = xdat, ydat = y, output = 'data', perspective = FALSE, gradients = TRUE, errors = 'bootstrap', band = 'pmzsd', bootstrap = 'wild', boot_control = np_boot_control(wild='rademacher'), B = B, random_seed = 33)",
      "expected.hw <- qnorm(0.05 / 2, lower.tail = FALSE) * apply(expected.t, 2L, stats::sd)",
      "stopifnot(isTRUE(all.equal(plot_target(out)$estimate, as.vector(expected.t0), tolerance = 1e-10)))",
      "stopifnot(isTRUE(all.equal(plot_target(out)$upper - plot_target(out)$estimate, expected.hw, tolerance = 1e-10)))",
      "cat('RBAND_ORDERED_CATGRAD_WILD_TARGET_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("RBAND_ORDERED_CATGRAD_WILD_TARGET_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session-route categorical regression gradient bootstrap works for default and wild routes", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      plot_test_support(),
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
      "fit <- npreg(bws = bw, txdat = xdat, tydat = y, gradients = TRUE, se = TRUE)",
      "for (boot.method in c('default', 'wild')) {",
      "  args <- list(fit, xdat = xdat, ydat = y, output = 'data', perspective = FALSE, gradients = TRUE, errors = 'bootstrap', B = 11L)",
      "  if (!identical(boot.method, 'default')) args$bootstrap <- boot.method",
      "  out <- suppressWarnings(do.call(plot, args))",
      "  stopifnot(inherits(out, 'np_plot_data'), length(out$panels[[1L]]$display$boxplots) > 0L, length(out$panels[[1L]]$display$boxplots[[1L]]$names) == 2L)",
      "  stopifnot(all(is.finite(out$panels[[1L]]$values$estimate)), all(is.finite(out$panels[[1L]]$intervals$lower)), all(is.finite(out$panels[[1L]]$intervals$upper)))",
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

test_that("session-route categorical gradient asymptotic intervals use native HC0 errors", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      plot_test_support(),
      "options(npRmpi.autodispatch = FALSE, np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(20260313)",
      "n <- 36L",
      "g <- factor(sample(c('a', 'b'), n, replace = TRUE))",
      "x <- runif(n)",
      "y <- 1 + 0.4 * (g == 'b') + cos(2 * pi * x) + rnorm(n, sd = 0.05)",
      "xdat <- data.frame(g = g, x = x)",
      "bw <- npregbw(xdat = xdat, ydat = y, regtype = 'll', bwtype = 'fixed', bws = c(0.25, 0.3), bandwidth.compute = FALSE)",
      "fit <- npreg(bws = bw, txdat = xdat, tydat = y, gradients = TRUE, se = TRUE)",
      "out.cat <- suppressWarnings(plot(fit, xdat = xdat, ydat = y, output = 'data', perspective = FALSE, gradients = TRUE, errors = 'asymptotic', data_overlay = FALSE))",
      "stopifnot(is.list(out.cat), inherits(out.cat, 'np_plot_data'))",
      "stopifnot(all(is.finite(plot_target(out.cat)$estimate)), all(is.finite(sweep(plot_band(out.cat, band='pmzsd'), 1L, plot_target(out.cat)$estimate, '-'))))",
      "oracle.cat <- npreg(bws = bw, txdat = xdat, tydat = y, exdat = plot_coordinates(out.cat, role='X'), gradients = TRUE, se = TRUE)",
      "expected.gerr <- cbind(-qnorm(0.975) * oracle.cat$gerr[, 1L], qnorm(0.975) * oracle.cat$gerr[, 1L])",
      "stopifnot(isTRUE(all.equal(as.vector(plot_target(out.cat)$estimate), as.vector(oracle.cat$grad[, 1L]), tolerance = 2e-12)))",
      "stopifnot(isTRUE(all.equal(as.vector(sweep(plot_band(out.cat, band='pmzsd'), 1L, plot_target(out.cat)$estimate, '-')), as.vector(expected.gerr), tolerance = 2e-12)))",
      "xdat2 <- data.frame(x = runif(n), z = runif(n))",
      "y2 <- 1 + sin(2 * pi * xdat2$x) + 0.5 * xdat2$z + rnorm(n, sd = 0.05)",
      "bw2 <- npregbw(xdat = xdat2, ydat = y2, regtype = 'll', bwtype = 'fixed', bws = c(0.25, 0.3), bandwidth.compute = FALSE)",
      "fit2 <- npreg(bws = bw2, txdat = xdat2, tydat = y2, gradients = TRUE, se = TRUE)",
      "out <- suppressWarnings(plot(fit2, xdat = xdat2, ydat = y2, output = 'data', perspective = FALSE, gradients = TRUE, errors = 'asymptotic', data_overlay = FALSE))",
      "stopifnot(is.list(out), inherits(out, 'np_plot_data'))",
      "cat('RBAND_CATGRAD_ASYMPTOTIC_CONTRACT_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("RBAND_CATGRAD_ASYMPTOTIC_CONTRACT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
