run_rscript_subprocess <- function(lines, timeout = 45L, env = character()) {
  script <- tempfile("npRmpi-subprocess-", fileext = ".R")
  writeLines(lines, script, useBytes = TRUE)
  on.exit(unlink(script), add = TRUE)

  cmd <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(cmd,
                                  c("--vanilla", script),
                                  stdout = TRUE,
                                  stderr = TRUE,
                                  timeout = timeout,
                                  env = env))
  status <- attr(out, "status")
  if (is.null(status))
    status <- 0L
  list(status = as.integer(status), output = out)
}

run_cmd_subprocess <- function(cmd, args = character(), timeout = 60L, env = character()) {
  out <- suppressWarnings(system2(cmd,
                                  args,
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

test_that("session routing smoke completes in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "set.seed(42)",
      "n <- 250",
      "x <- runif(n)",
      "y <- rnorm(n)",
      "bw <- npregbw(y~x, regtype='lc', bwmethod='cv.ls', nmulti=1)",
      "fit <- npreg(bws=bw, gradients=FALSE)",
      "stopifnot(inherits(fit, 'npregression'))",
      "cat('SESSION_ROUTE_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_ROUTE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npcdens user-style example completes with quiet=FALSE", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=FALSE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "set.seed(42)",
      "n <- 1000",
      "x <- rnorm(n)",
      "y <- rnorm(n)",
      "F <- npcdens(y~x)",
      "summary(F$bws)",
      "png(tempfile(fileext='.png'))",
      "plot(F)",
      "dev.off()",
      "stopifnot(inherits(F, 'condensity'))",
      "cat('SESSION_NPCDENS_EXAMPLE_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPCDENS_EXAMPLE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npreg factor example completes with quiet=FALSE", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=FALSE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "set.seed(42)",
      "n <- 250",
      "x <- runif(n)",
      "z1 <- rbinom(n,1,.5)",
      "z2 <- rbinom(n,1,.5)",
      "y <- cos(2*pi*x) + z1 + rnorm(n,sd=.25)",
      "z1 <- factor(z1)",
      "z2 <- factor(z2)",
      "bw <- npregbw(y~x+z1+z2, regtype='lc', bwmethod='cv.ls', nmulti=1)",
      "summary(bw)",
      "fit <- npreg(bws=bw, gradients=FALSE)",
      "summary(fit)",
      "stopifnot(inherits(bw, 'rbandwidth'), inherits(fit, 'npregression'))",
      "cat('SESSION_NPREG_FACTORS_EXAMPLE_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPREG_FACTORS_EXAMPLE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session core density/distribution family smoke completes", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "set.seed(7)",
      "n <- 140",
      "x <- rnorm(n)",
      "y <- rnorm(n)",
      "xd <- data.frame(x = x)",
      "yd <- data.frame(y = y)",
      "bw_ud <- npudensbw(dat = xd, bws = 0.4, bandwidth.compute = FALSE)",
      "fit_ud <- npudens(tdat = xd, bws = bw_ud)",
      "bw_uf <- npudistbw(dat = xd, bws = 0.4, bandwidth.compute = FALSE)",
      "fit_uf <- npudist(tdat = xd, bws = bw_uf)",
      "bw_cf <- npcdistbw(xdat = xd, ydat = yd, bws = c(0.4, 0.4), bandwidth.compute = FALSE)",
      "fit_cf <- npcdist(txdat = xd, tydat = yd, bws = bw_cf)",
      "stopifnot(inherits(fit_ud, 'npdensity'), inherits(fit_uf, 'npdistribution'), inherits(fit_cf, 'condistribution'))",
      "cat('SESSION_DENS_DIST_CORE_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_DENS_DIST_CORE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("skip-init mode fails fast without MPI pool crash", {
  skip_on_cran()
  env <- subprocess_env(extra = "NP_RMPI_SKIP_INIT=1")
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.autodispatch=FALSE)",
      "d <- data.frame(x=c(0.1,0.2))",
      "bws <- list()",
      "err <- tryCatch({",
      "  npRmpi:::npregbw.rbandwidth(xdat=d, ydat=d$x, bws=bws)",
      "  'NO_ERROR'",
      "}, error=function(e) conditionMessage(e))",
      "cat(err, '\\n')"
    ),
    timeout = 20L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("requires an active MPI slave pool", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("manual-broadcast mode smoke completes in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=FALSE)",
      "set.seed(42)",
      "n <- 120",
      "x <- runif(n)",
      "y <- rnorm(n)",
      "d <- data.frame(x=x, y=y)",
      "mpi.bcast.Robj2slave(d)",
      "mpi.bcast.cmd(bw <- npregbw(y~x, data=d, regtype='lc', bwmethod='cv.ls', nmulti=1), caller.execute=TRUE)",
      "mpi.bcast.cmd(fit <- npreg(bws=bw, gradients=FALSE), caller.execute=TRUE)",
      "stopifnot(inherits(fit, 'npregression'))",
      "cat('MANUAL_BCAST_ROUTE_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("MANUAL_BCAST_ROUTE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("attach mode smoke completes under mpiexec when enabled", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_TEST=1 to run attach-mode smoke")
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  script <- tempfile("npRmpi-attach-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(mode='attach', quiet=TRUE)",
    "if (mpi.comm.rank(1L) == 0L) {",
    "  set.seed(42)",
    "  n <- 120",
    "  x <- runif(n)",
    "  y <- rnorm(n)",
    "  bw <- npregbw(y~x, regtype='lc', bwmethod='cv.ls', nmulti=1)",
    "  fit <- npreg(bws=bw, gradients=FALSE)",
    "  stopifnot(inherits(fit, 'npregression'))",
    "  cat('ATTACH_ROUTE_OK\\n')",
    "  npRmpi.quit(mode='attach')",
    "}"
  ), script, useBytes = TRUE)

  env_common <- subprocess_env()
  skip_if(is.null(env_common), "local npRmpi install unavailable for subprocess smoke")
  res <- run_cmd_subprocess(
    mpiexec,
    args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--vanilla", script),
    timeout = 90L,
    env = c(env_common, "FI_TCP_IFACE=en0", "FI_PROVIDER=tcp", "FI_SOCKETS_IFACE=en0")
  )
  if (res$status != 0L) {
    res <- run_cmd_subprocess(
      mpiexec,
      args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--vanilla", script),
      timeout = 90L,
      env = c(env_common, "FI_TCP_IFACE=lo0", "FI_PROVIDER=tcp", "FI_SOCKETS_IFACE=lo0")
    )
  }

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("ATTACH_ROUTE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
