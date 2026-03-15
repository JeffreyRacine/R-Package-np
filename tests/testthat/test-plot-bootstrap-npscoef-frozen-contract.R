npscoef_frozen_plot_case <- function(label, bwtype, boot.method) {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- sprintf("NPSCOEF_FROZEN_%s_OK", label)
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "run_case <- function() {",
    "  npRmpi.init(nslaves=1, quiet=TRUE)",
    "  on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "  options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "  set.seed(42)",
    "  n <- 80L",
    "  x <- runif(n)",
    "  z <- runif(n)",
    "  y <- x * (1 + z) + rnorm(n, sd=0.05)",
    sprintf("  fit <- suppressWarnings(npscoef(y ~ x | z, regtype='ll', bwtype='%s', nmulti=1, betas=TRUE))", bwtype),
    "  tf <- tempfile(fileext='.pdf')",
    "  grDevices::pdf(tf)",
    "  on.exit({ try(grDevices::dev.off(), silent=TRUE); unlink(tf) }, add=TRUE)",
    sprintf("  args <- list(x=fit, neval=20L, coef=FALSE, plot.errors.method='bootstrap', plot.errors.boot.method='%s', plot.errors.boot.nonfixed='frozen', plot.errors.boot.num=41L, plot.errors.type='pointwise')", boot.method),
    "  if (identical(args$plot.errors.boot.method, 'geom')) args$plot.errors.boot.blocklen <- 4L",
    "  capture.output(do.call(plot, args))",
    "}",
    "run_case()",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(
    lines = lines,
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
}

test_that("npRmpi npscoef generalized frozen inid is supported", {
  npscoef_frozen_plot_case("GENERALIZED_INID", "generalized_nn", "inid")
})

test_that("npRmpi npscoef generalized frozen geom is supported", {
  npscoef_frozen_plot_case("GENERALIZED_GEOM", "generalized_nn", "geom")
})

test_that("npRmpi npscoef adaptive frozen inid is supported", {
  npscoef_frozen_plot_case("ADAPTIVE_INID", "adaptive_nn", "inid")
})

test_that("npRmpi npscoef adaptive frozen geom is supported", {
  npscoef_frozen_plot_case("ADAPTIVE_GEOM", "adaptive_nn", "geom")
})
