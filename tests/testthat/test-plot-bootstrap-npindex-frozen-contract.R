npindex_frozen_plot_case <- function(label, boot.method) {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- sprintf("NPINDEX_FROZEN_%s_OK", label)
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "set.seed(42)",
    "n <- 80L",
    "x <- runif(n)",
    "z <- rnorm(n)",
    "y <- x + rnorm(n)",
    "fit <- suppressWarnings(npindex(y ~ x + z, regtype='ll', bwtype='adaptive_nn', nmulti=1))",
    "tf <- tempfile(fileext='.pdf')",
    "grDevices::pdf(tf)",
    "on.exit({ try(grDevices::dev.off(), silent=TRUE); unlink(tf) }, add=TRUE)",
    sprintf("args <- list(x=fit, neval=20L, plot.errors.method='bootstrap', plot.errors.boot.method='%s', plot.errors.boot.nonfixed='frozen', plot.errors.boot.num=41L, plot.errors.type='pointwise', gradients=TRUE)", boot.method),
    "if (identical(args$plot.errors.boot.method, 'geom')) args$plot.errors.boot.blocklen <- 4L",
    "capture.output(do.call(plot, args))",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(
    lines = lines,
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
}

test_that("npRmpi npindex frozen inid stays on the exact public scale", {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(42)",
      "n <- 120L",
      "x <- runif(n, -1, 1)",
      "z <- rnorm(n)",
      "y <- x^2 + rnorm(n, sd = 0.25 * stats::sd(x))",
      "fit <- suppressWarnings(npindex(y ~ x + z, bwtype='adaptive_nn', nmulti=1L))",
      "get_obj <- function(mode) suppressWarnings(plot(fit, plot.behavior='data', neval=40L, plot.errors.method='bootstrap', plot.errors.boot.method='inid', plot.errors.boot.nonfixed=mode, plot.errors.boot.num=39L, plot.errors.type='pointwise'))[[1L]]",
      "frozen <- get_obj('frozen')",
      "exact <- get_obj('exact')",
      "ratio <- stats::median(abs(exact$merr[, 1L]) / pmax(abs(frozen$merr[, 1L]), 1e-12), na.rm = TRUE)",
      "stopifnot(max(abs(frozen$mean - exact$mean)) < 1e-12)",
      "stopifnot(is.finite(ratio), ratio > 0.5, ratio < 2.0)",
      "cat('NPINDEX_FROZEN_INID_SCALE_OK\\n')"
    ),
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPINDEX_FROZEN_INID_SCALE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("npRmpi npindex frozen inid still supports gradient slices", {
  npindex_frozen_plot_case("INID", "inid")
})

test_that("npRmpi npindex frozen geom still supports gradient slices", {
  npindex_frozen_plot_case("GEOM", "geom")
})
