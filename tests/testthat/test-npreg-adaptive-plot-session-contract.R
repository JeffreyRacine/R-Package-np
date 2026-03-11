.npreg_adaptive_mpi_init_env_failure <- function(output) {
  any(grepl("OFI call ep_enable failed", output, fixed = TRUE)) ||
    any(grepl("Fatal error in internal_Init", output, fixed = TRUE)) ||
    any(grepl("MPI_Init", output, fixed = TRUE) & grepl("failed", output, ignore.case = TRUE))
}

test_that("session adaptive-nn regression wild gradient plot-data completes in subprocess", {
  skip_on_cran()

  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(42)",
      "n <- 120L",
      "x <- sort(runif(n))",
      "y <- sin(2*pi*x) + rnorm(n, sd=0.1)",
      "xdat <- data.frame(x=x)",
      "bw <- npregbw(xdat=xdat, ydat=y, regtype='ll', bwtype='adaptive_nn', bws=35L, bandwidth.compute=FALSE)",
      "fit <- npreg(bws=bw, txdat=xdat, tydat=y, warn.glp.gradient=FALSE)",
      "out <- suppressWarnings(plot(",
      "  fit,",
      "  xdat=xdat,",
      "  ydat=y,",
      "  plot.behavior='data',",
      "  perspective=FALSE,",
      "  neval=40L,",
      "  plot.errors.method='bootstrap',",
      "  plot.errors.boot.method='wild',",
      "  plot.errors.boot.num=39L,",
      "  plot.errors.type='all',",
      "  gradients=TRUE))",
      "stopifnot(is.list(out))",
      "stopifnot(inherits(out[[1L]], 'npregression'))",
      "stopifnot(length(out[[1L]]$mean) == 40L)",
      "cat('SESSION_NPREG_ADAPTIVE_PLOT_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  if (res$status != 0L && .npreg_adaptive_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for session-mode smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPREG_ADAPTIVE_PLOT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session adaptive-nn regression object-fed wild gradient plot-data completes in subprocess", {
  skip_on_cran()

  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(42)",
      "n <- 250L",
      "x <- runif(n)",
      "y <- x + rnorm(n)",
      "dat <- data.frame(x=x, y=y)",
      "bw <- npregbw(y ~ x, data=dat, regtype='ll', bwtype='adaptive_nn', bws=35L, bandwidth.compute=FALSE)",
      "fit <- npreg(bws=bw, data=dat, warn.glp.gradient=FALSE)",
      "out <- suppressWarnings(plot(",
      "  fit,",
      "  plot.behavior='data',",
      "  neval=50L,",
      "  plot.errors.method='bootstrap',",
      "  plot.errors.boot.method='wild',",
      "  plot.errors.boot.num=29L,",
      "  plot.errors.type='all',",
      "  gradients=TRUE))",
      "stopifnot(is.list(out))",
      "stopifnot(inherits(out[[1L]], 'npregression'))",
      "stopifnot(length(out[[1L]]$mean) == 50L)",
      "cat('SESSION_NPREG_ADAPTIVE_OBJECT_OK\\n')"
    ),
    timeout = 180L,
    env = env
  )

  if (res$status != 0L && .npreg_adaptive_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for session-mode smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPREG_ADAPTIVE_OBJECT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
