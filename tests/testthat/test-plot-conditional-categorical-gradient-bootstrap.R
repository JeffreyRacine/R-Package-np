test_that("session-route conditional categorical gradient helper matches explicit local refits", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.autodispatch = FALSE, np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "helper <- getFromNamespace('.np_inid_boot_from_conditional_gradient_local', 'npRmpi')",
      "eval_fun <- getFromNamespace('.np_plot_conditional_eval', 'npRmpi')",
      "set.seed(20260312)",
      "n <- 24L",
      "x1 <- factor(sample(c('a', 'b'), n, replace = TRUE))",
      "x2 <- rnorm(n)",
      "y <- rnorm(n)",
      "xdat <- data.frame(x1 = x1, x2 = x2)",
      "ydat <- data.frame(y = y)",
      "fixtures <- list(",
      "  list(label = 'dens', cdf = FALSE, bw = suppressWarnings(npcdensbw(xdat = xdat, ydat = ydat, nmulti = 1L))),",
      "  list(label = 'dist', cdf = TRUE, bw = suppressWarnings(npcdistbw(xdat = xdat, ydat = ydat, nmulti = 1L)))",
      ")",
      "counts <- cbind(rep(1L, n), c(rep(2L, 4), rep(0L, 4), rep(1L, n - 8)), c(rep(0L, 3), rep(3L, 3), rep(1L, n - 6)))",
      "for (fixture in fixtures) {",
      "  slice <- suppressWarnings(plot(fixture$bw, xdat = xdat, ydat = ydat, gradients = TRUE, perspective = FALSE, plot.behavior = 'data'))[[1L]]",
      "  boot <- helper(xdat = xdat, ydat = ydat, exdat = slice$xeval, eydat = slice$yeval, bws = fixture$bw, B = ncol(counts), cdf = fixture$cdf, gradient.index = 1L, counts = counts)",
      "  explicit_fit <- function(idx) {",
      "    as.vector(eval_fun(bws = fixture$bw, xdat = xdat[idx, , drop = FALSE], ydat = ydat[idx, , drop = FALSE], exdat = slice$xeval, eydat = slice$yeval, cdf = fixture$cdf, gradients = TRUE)$congrad[, 1L])",
      "  }",
      "  oracle.t0 <- explicit_fit(seq_len(n))",
      "  oracle.t <- vapply(seq_len(ncol(counts)), function(j) explicit_fit(rep.int(seq_len(n), counts[, j])), numeric(length(oracle.t0)))",
      "  stopifnot(isTRUE(all.equal(as.vector(boot$t0), oracle.t0, tolerance = 1e-8)))",
      "  stopifnot(isTRUE(all.equal(boot$t, t(oracle.t), tolerance = 1e-8)))",
      "}",
      "cat('COND_CATGRAD_HELPER_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("COND_CATGRAD_HELPER_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session-route conditional categorical bootstrap gradients work again", {
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
      "n <- 24L",
      "x1 <- factor(sample(c('a', 'b'), n, replace = TRUE))",
      "x2 <- rnorm(n)",
      "y <- rnorm(n)",
      "xdat <- data.frame(x1 = x1, x2 = x2)",
      "ydat <- data.frame(y = y)",
      "cd.bw <- suppressWarnings(npcdensbw(xdat = xdat, ydat = ydat, nmulti = 1L))",
      "cdist.bw <- suppressWarnings(npcdistbw(xdat = xdat, ydat = ydat, nmulti = 1L))",
      "cases <- list(",
      "  list(label = 'npcdens-bw', obj = cd.bw, args = list(xdat = xdat, ydat = ydat)),",
      "  list(label = 'npcdens-fit', obj = npcdens(bws = cd.bw), args = list(xdat = xdat, ydat = ydat)),",
      "  list(label = 'npcdist-bw', obj = cdist.bw, args = list(xdat = xdat, ydat = ydat)),",
      "  list(label = 'npcdist-fit', obj = npcdist(bws = cdist.bw), args = list(xdat = xdat, ydat = ydat))",
      ")",
      "for (case in cases) {",
      "  for (boot.method in c('inid', 'fixed', 'geom')) {",
      "    plot.args <- c(list(case$obj, perspective = FALSE, plot.behavior = 'data', gradients = TRUE, plot.errors.method = 'bootstrap', plot.errors.boot.method = boot.method, plot.errors.boot.num = 7L), case$args)",
      "    out <- suppressWarnings(do.call(plot, plot.args))",
      "    stopifnot(is.list(out), length(out) >= 1L, length(out[[1L]]$bxp) > 0L, length(out[[1L]]$bxp$names) == 2L)",
      "  }",
      "}",
      "cat('COND_CATGRAD_CONSUMER_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("COND_CATGRAD_CONSUMER_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
