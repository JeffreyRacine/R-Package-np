local_npRmpi_formula_subprocess_env <- function(extra = character()) {
  pkg.root <- tryCatch(
    normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
    error = function(e) ""
  )
  if (!nzchar(pkg.root))
    return(NULL)

  src.copy <- tempfile("npRmpi-formula-subprocess-src-")
  dir.create(src.copy, recursive = TRUE, showWarnings = FALSE)
  pkg.copy <- file.path(src.copy, "np-npRmpi")
  if (!isTRUE(file.copy(pkg.root, src.copy, recursive = TRUE)))
    return(NULL)

  vignette.index <- file.path(pkg.copy, "build", "vignette.rds")
  if (file.exists(vignette.index))
    unlink(vignette.index, force = TRUE)

  build.dir <- tempfile("npRmpi-formula-subprocess-build-")
  dir.create(build.dir, recursive = TRUE, showWarnings = FALSE)

  old.wd <- getwd()
  on.exit(setwd(old.wd), add = TRUE)
  setwd(build.dir)

  build.out <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "build", "--no-build-vignettes", "--no-manual", pkg.copy),
    stdout = TRUE,
    stderr = TRUE
  ))
  build.status <- attr(build.out, "status")
  if (is.null(build.status))
    build.status <- 0L
  if (build.status != 0L)
    return(NULL)

  tarballs <- Sys.glob(file.path(build.dir, "npRmpi_*.tar.gz"))
  if (!length(tarballs))
    return(NULL)

  lib.path <- tempfile("npRmpi-formula-subprocess-lib-")
  dir.create(lib.path, recursive = TRUE, showWarnings = FALSE)

  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-test-load", "-l", lib.path, tail(tarballs, 1L)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(out, "status")
  if (is.null(status))
    status <- 0L

  if (status != 0L) {
    unlink(lib.path, recursive = TRUE, force = TRUE)
    return(NULL)
  }

  c(
    sprintf("R_LIBS=%s", paste(c(lib.path, .libPaths()), collapse = .Platform$path.sep)),
    "NP_RMPI_NO_REUSE_SLAVES=1",
    extra
  )
}

test_that("npreg formula fits preserve the response name for plotting metadata in session mode", {
  skip_on_cran()
  env <- local_npRmpi_formula_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(42)",
      "n <- 60L",
      "dat <- data.frame(x = runif(n, -1, 1))",
      "dat$y <- dat$x^2 + rnorm(n, sd = 0.25 * stats::sd(dat$x))",
      "fit <- npreg(y ~ x, data = dat, nmulti = 1, regtype = 'lp', degree = 1, bwtype = 'adaptive_nn')",
      "stopifnot(identical(fit$bws$ynames, 'y'))",
      "png(tempfile(fileext = '.png'), width = 800, height = 600)",
      "plot(fit, plot.errors.method = 'asymptotic', view = 'fixed')",
      "dev.off()",
      "cat('NPREG_FORMULA_YNAME_OK\\n')"
    ),
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPREG_FORMULA_YNAME_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
