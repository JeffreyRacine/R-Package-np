npreg_frozen_plot_case <- function(label, bwtype, boot.method) {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- sprintf("NPREG_FROZEN_%s_OK", label)
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "set.seed(42)",
    "n <- 80L",
    "x <- runif(n)",
    "y <- x + rnorm(n)",
    sprintf("fit <- suppressWarnings(npreg(y ~ x, regtype='ll', degree=c(1), bernstein=TRUE, bwtype='%s', nmulti=1))", bwtype),
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

test_that("npRmpi npreg generalized frozen inid supports continuous gradients", {
  npreg_frozen_plot_case("GENERALIZED_INID", "generalized_nn", "inid")
})

test_that("npRmpi npreg generalized frozen geom supports continuous gradients", {
  npreg_frozen_plot_case("GENERALIZED_GEOM", "generalized_nn", "geom")
})

test_that("npRmpi npreg adaptive frozen inid supports continuous gradients", {
  npreg_frozen_plot_case("ADAPTIVE_INID", "adaptive_nn", "inid")
})

test_that("npRmpi npreg adaptive frozen geom supports continuous gradients", {
  npreg_frozen_plot_case("ADAPTIVE_GEOM", "adaptive_nn", "geom")
})
