cleanup_npreg_adaptive_subprocess_mpi <- function() {
  if (.Platform$OS.type != "unix")
    return(invisible(NULL))

  suppressWarnings(try(system2("pkill", c("-f", "slavedaemon.R")), silent = TRUE))
  suppressWarnings(try(system2("pkill", c("-f", "Rslaves.sh")), silent = TRUE))
  Sys.sleep(0.2)
  invisible(NULL)
}

run_npreg_adaptive_rscript_subprocess <- function(lines, timeout = 45L, env = character()) {
  script <- tempfile("npRmpi-adaptive-npreg-", fileext = ".R")
  log <- tempfile("npRmpi-adaptive-npreg-", fileext = ".log")
  writeLines(lines, script, useBytes = TRUE)
  cleanup_npreg_adaptive_subprocess_mpi()
  on.exit(cleanup_npreg_adaptive_subprocess_mpi(), add = TRUE)
  on.exit(unlink(c(script, log)), add = TRUE)

  cmd <- file.path(R.home("bin"), "Rscript")
  status <- suppressWarnings(system2(
    cmd,
    c("--no-save", script),
    stdout = log,
    stderr = log,
    timeout = timeout,
    env = env
  ))
  exit.status <- attr(status, "status")
  if (is.null(exit.status))
    exit.status <- 0L
  out <- if (file.exists(log)) readLines(log, warn = FALSE) else character()
  list(status = as.integer(exit.status), output = out)
}

.npreg_adaptive_mpi_init_env_failure <- function(output) {
  any(grepl("OFI call ep_enable failed", output, fixed = TRUE)) ||
    any(grepl("Fatal error in internal_Init", output, fixed = TRUE)) ||
    any(grepl("MPI_Init", output, fixed = TRUE) & grepl("failed", output, ignore.case = TRUE))
}

ensure_npreg_adaptive_subprocess_lib <- local({
  lib.path.cache <- NULL

  function() {
    if (!is.null(lib.path.cache) && dir.exists(lib.path.cache))
      return(lib.path.cache)

    loaded.pkg <- tryCatch(find.package("npRmpi"), error = function(e) "")
    if (nzchar(loaded.pkg)) {
      loaded.lib <- tryCatch(
        normalizePath(file.path(loaded.pkg, ".."), mustWork = TRUE),
        error = function(e) ""
      )
      if (nzchar(loaded.lib) && dir.exists(loaded.lib)) {
        lib.path.cache <<- loaded.lib
        return(lib.path.cache)
      }
    }

    pkg.root <- tryCatch(
      normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
      error = function(e) ""
    )
    if (!nzchar(pkg.root))
      return(NULL)

    lib.path.cache <<- tempfile("npRmpi-adaptive-npreg-lib-")
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

npreg_adaptive_subprocess_env <- function(extra = character()) {
  lib.path <- ensure_npreg_adaptive_subprocess_lib()
  if (is.null(lib.path))
    return(NULL)

  c(
    sprintf("R_LIBS=%s", paste(c(lib.path, .libPaths()), collapse = .Platform$path.sep)),
    extra
  )
}

test_that("session adaptive-nn regression wild gradient plot-data completes in subprocess", {
  env <- npreg_adaptive_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_npreg_adaptive_rscript_subprocess(
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
  env <- npreg_adaptive_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_npreg_adaptive_rscript_subprocess(
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
