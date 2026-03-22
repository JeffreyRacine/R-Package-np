run_rscript_subprocess <- function(lines, timeout = 45L, env = character()) {
  script <- tempfile("npRmpi-subprocess-", fileext = ".R")
  writeLines(lines, script, useBytes = TRUE)
  on.exit(unlink(script), add = TRUE)

  cmd <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(cmd,
                                  c("--no-save", script),
                                  stdout = TRUE,
                                  stderr = TRUE,
                                  timeout = timeout,
                                  env = env))
  status <- attr(out, "status")
  if (is.null(status))
    status <- 0L
  list(status = as.integer(status), output = out)
}

ensure_subprocess_npRmpi_lib <- local({
  lib.path.cache <- NULL

  function() {
    if (!is.null(lib.path.cache) && dir.exists(lib.path.cache))
      return(lib.path.cache)

    pkg.root <- tryCatch(
      normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
      error = function(e) ""
    )
    if (!nzchar(pkg.root))
      return(NULL)

    lib.path.cache <<- tempfile("npRmpi-subprocess-lib-")
    dir.create(lib.path.cache, recursive = TRUE, showWarnings = FALSE)

    cmd <- file.path(R.home("bin"), "R")
    out <- suppressWarnings(system2(
      cmd,
      c("CMD", "INSTALL", "--no-test-load", "-l", lib.path.cache, pkg.root),
      stdout = TRUE,
      stderr = TRUE
    ))
    status <- attr(out, "status")
    if (is.null(status))
      status <- 0L

    if (status != 0L) {
      warning(paste(out, collapse = "\n"))
      unlink(lib.path.cache, recursive = TRUE, force = TRUE)
      lib.path.cache <<- NULL
      return(NULL)
    }

    lib.path.cache
  }
})

subprocess_env <- function(extra = character()) {
  lib.path <- ensure_subprocess_npRmpi_lib()
  if (is.null(lib.path))
    return(NULL)

  c(
    sprintf("R_LIBS=%s", paste(c(lib.path, .libPaths()), collapse = .Platform$path.sep)),
    extra
  )
}

test_that("session adaptive conditional plot-data smoke completes in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "set.seed(42)",
      "n <- 250",
      "x <- runif(n)",
      "z <- runif(n)",
      "y <- x^2 + rnorm(n, sd=0.1)",
      "d <- data.frame(x=x, z=z, y=y)",
      "fit.cd <- npcdens(y~x, data=d, regtype='lp', degree=1L, bwtype='adaptive_nn')",
      "fit.cf <- npcdist(y~x, data=d, regtype='lp', degree=1L, bwtype='adaptive_nn')",
      "png(tempfile(fileext='.png'))",
      "on.exit(dev.off(), add=TRUE)",
      "out.cd <- suppressWarnings(plot(fit.cd, plot.behavior='data'))",
      "out.cf <- suppressWarnings(plot(fit.cf, plot.behavior='data'))",
      "stopifnot(is.list(out.cd), length(out.cd) > 0L)",
      "stopifnot(is.list(out.cf), length(out.cf) > 0L)",
      "cat('SESSION_ADAPTIVE_CONDITIONAL_PLOT_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_ADAPTIVE_CONDITIONAL_PLOT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
