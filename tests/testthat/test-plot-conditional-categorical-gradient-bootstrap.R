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
      "  slice <- suppressWarnings(plot(fixture$bw, xdat = xdat, ydat = ydat, gradients = TRUE, perspective = FALSE, output = 'data'))[[1L]]",
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
      "    plot.args <- c(list(case$obj, perspective = FALSE, output = 'data', gradients = TRUE, errors = 'bootstrap', bootstrap = boot.method, B = 7L), case$args)",
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
test_that("conditional gradient coordinates retain physical columns", {
  resolve <- getFromNamespace(".np_plot_resolve_conditional_gradient_index", "npRmpi")
  ordinary <- getFromNamespace(".np_plot_require_conditional_gradient_bootstrap_supported", "npRmpi")
  bias <- getFromNamespace(".np_plot_validate_conditional_gradient_target", "npRmpi")
  label <- getFromNamespace(".np_plot_conditional_bootstrap_target_label", "npRmpi")
  layouts <- list(c(FALSE, FALSE, TRUE, TRUE), c(TRUE, FALSE, TRUE, FALSE),
                  c(FALSE, TRUE), c(TRUE, TRUE), c(FALSE, FALSE))
  for (icon in layouts) {
    bw <- list(xdati = list(icon = icon), xndim = length(icon), yndim = 1L,
               xnames = paste0("x", seq_along(icon)), ynames = "y")
    for (j in seq_along(icon)) {
      expect_identical(resolve(bw, j, "test"), j)
      expect_true(ordinary(bw, j, "test"))
      expect_match(label(bw, 1L, TRUE, j), paste0("grad x", j, " on x1"), fixed = TRUE)
      if (icon[j])
        expect_identical(bias(bw, j, "test"), j)
      else
        expect_error(bias(bw, j, "test"), "smooth-bootstrap gradient bias correction", fixed = TRUE)
    }
  }
  bw <- list(xdati = list(icon = c(FALSE, TRUE)))
  for (j in list(NULL, numeric(), NA_real_, NaN, Inf, -Inf, 0L, -1L,
                 3L, 1.5, c(1L, 2L), "1", TRUE, 1 + 1i)) {
    expect_error(resolve(bw, j, "test"), "invalid conditional gradient coordinate", fixed = TRUE)
  }
  expect_error(resolve(list(xdati = list(icon = c(TRUE, NA))), 1L, "test"),
               "invalid conditional gradient coordinate", fixed = TRUE)
})
