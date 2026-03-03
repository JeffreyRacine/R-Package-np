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

.is_mpi_init_env_failure <- function(output) {
  any(grepl("OFI call ep_enable failed", output, fixed = TRUE)) ||
    any(grepl("Fatal error in internal_Init", output, fixed = TRUE)) ||
    any(grepl("MPI_Init", output, fixed = TRUE) & grepl("failed", output, ignore.case = TRUE))
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

test_that("session master-only mode (nslaves=0) runs estimator workflow", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=0, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(mode='spawn', force=TRUE), silent=TRUE), add=TRUE)",
      "stopifnot(isTRUE(getOption('npRmpi.master.only')))",
      "stopifnot(identical(getOption('npRmpi.autodispatch'), FALSE))",
      "set.seed(42)",
      "n <- 250",
      "x <- runif(n)",
      "y <- rnorm(n)",
      "bw <- npregbw(y~x, regtype='lc', bwmethod='cv.ls', nmulti=1)",
      "fit <- npreg(bws=bw)",
      "pred <- predict(fit, newdata=data.frame(x=seq(0,1,length.out=60)))",
      "stopifnot(inherits(fit, 'npregression'), length(pred)==60, all(is.finite(pred)))",
      "cat('SESSION_MASTER_ONLY_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_MASTER_ONLY_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session master-only mode disables autodispatch to local fallback", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=0, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(mode='spawn', force=TRUE), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE)",
      "set.seed(7)",
      "x <- runif(80)",
      "y <- rnorm(80)",
      "bw <- npregbw(xdat=data.frame(x=x), ydat=y, bws=0.2, bandwidth.compute=FALSE, regtype='ll', bwtype='fixed', ckertype='gaussian')",
      "fit <- npreg(bws=bw, exdat=data.frame(x=seq(0,1,length.out=20)))",
      "stopifnot(inherits(fit, 'npregression'))",
      "cat('SESSION_MASTER_ONLY_AUTOD_FALLBACK_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_MASTER_ONLY_AUTOD_FALLBACK_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("master-only mode", res$output, ignore.case = TRUE)),
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

test_that("session npcdens user-style example completes with default quiet", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "set.seed(42)",
      "n <- 500",
      "x <- rnorm(n)",
      "y <- rnorm(n)",
      "F <- npcdens(y~x)",
      "summary(F$bws)",
      "png(tempfile(fileext='.png'))",
      "plot(F)",
      "dev.off()",
      "stopifnot(inherits(F, 'condensity'))",
      "cat('SESSION_NPCDENS_DEFAULT_QUIET_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPCDENS_DEFAULT_QUIET_OK", res$output, fixed = TRUE)),
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

test_that("session npreghat smoke completes in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=FALSE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "set.seed(31415)",
      "n <- 120",
      "x <- runif(n)",
      "y <- sin(2*pi*x) + rnorm(n, sd=0.1)",
      "bw <- npregbw(y~x, regtype='ll', bwmethod='cv.ls', nmulti=1)",
      "fit <- npreg(bws=bw, gradients=FALSE)",
      "H <- npreghat(bws=bw, txdat=data.frame(x=x))",
      "hy <- as.vector(H %*% y)",
      "stopifnot(inherits(H, 'npreghat'))",
      "stopifnot(max(abs(hy - fitted(fit))) < 1e-6)",
      "cat('SESSION_NPREGHAT_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPREGHAT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session wild selector plot smoke completes in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=FALSE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "set.seed(27182)",
      "n <- 80",
      "x <- runif(n)",
      "y <- sin(2*pi*x) + rnorm(n, sd=0.1)",
      "bw <- npregbw(y~x, bws=0.2, bandwidth.compute=FALSE)",
      "old.chunk <- getOption('np.plot.wild.chunk.size')",
      "on.exit(options(np.plot.wild.chunk.size = old.chunk), add=TRUE)",
      "options(np.plot.wild.chunk.size = 5L)",
      "png(tempfile(fileext='.png'))",
      "on.exit(dev.off(), add=TRUE)",
      "out <- suppressWarnings(plot(",
      "  bw,",
      "  xdat=data.frame(x=x),",
      "  ydat=y,",
      "  plot.behavior='data',",
      "  plot.errors.method='bootstrap',",
      "  plot.errors.boot.method='wild',",
      "  plot.errors.boot.num=9",
      "))",
      "stopifnot(is.list(out), length(out) > 0)",
      "cat('SESSION_WILD_PLOT_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_WILD_PLOT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session wild then inid plot sequence keeps worker pool responsive", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "set.seed(27182)",
      "n <- 120",
      "x <- runif(n)",
      "y <- sin(2*pi*x) + rnorm(n, sd=0.1)",
      "bw <- npregbw(y~x, bws=0.2, bandwidth.compute=FALSE)",
      "out.w <- suppressWarnings(plot(",
      "  bw,",
      "  xdat=data.frame(x=x),",
      "  ydat=y,",
      "  plot.behavior='data',",
      "  plot.errors.method='bootstrap',",
      "  plot.errors.boot.method='wild',",
      "  plot.errors.boot.num=40",
      "))",
      "out.i <- suppressWarnings(plot(",
      "  bw,",
      "  xdat=data.frame(x=x),",
      "  ydat=y,",
      "  plot.behavior='data',",
      "  plot.errors.method='bootstrap',",
      "  plot.errors.boot.method='inid',",
      "  plot.errors.boot.num=40",
      "))",
      "stopifnot(is.list(out.w), length(out.w) > 0)",
      "stopifnot(is.list(out.i), length(out.i) > 0)",
      "r <- npRmpi:::mpi.remote.exec(1+1)",
      "stopifnot(length(r) >= 1L)",
      "cat('SESSION_WILD_INID_SEQUENCE_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_WILD_INID_SEQUENCE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session inid plot smoke completes in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "set.seed(42)",
      "n <- 1000",
      "x <- rnorm(n)",
      "y <- rnorm(n)",
      "g <- npreg(y ~ x)",
      "png(tempfile(fileext='.png'))",
      "on.exit(dev.off(), add=TRUE)",
      "suppressWarnings(plot(g))",
      "suppressWarnings(plot(",
      "  g,",
      "  plot.errors.method='bootstrap',",
      "  plot.errors.boot.method='inid',",
      "  plot.errors.boot.num=999",
      "))",
      "cat('SESSION_INID_PLOT_OK\\n')"
    ),
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_INID_PLOT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session inid density plot smoke completes in subprocess", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_DENSITY_INID_TEST"), "1"),
              "set NP_RMPI_ENABLE_DENSITY_INID_TEST=1 to run density inid session smoke")
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(np.plot.inid.ksum.fastpath.nprmpi=TRUE)",
      "set.seed(77)",
      "n <- 80",
      "x <- rnorm(n)",
      "y <- rnorm(n)",
      "xd <- data.frame(x=x)",
      "yd <- data.frame(y=y)",
      "bw_ud <- npudensbw(dat=xd, bws=0.35, bandwidth.compute=FALSE)",
      "bw_cd <- npcdensbw(xdat=xd, ydat=yd, bws=c(0.45,0.45), bandwidth.compute=FALSE)",
      "png(tempfile(fileext='.png'))",
      "on.exit(dev.off(), add=TRUE)",
      "out.ud <- suppressWarnings(plot(",
      "  bw_ud,",
      "  xdat=xd,",
      "  plot.behavior='data',",
      "  plot.errors.method='bootstrap',",
      "  plot.errors.boot.method='inid',",
      "  plot.errors.boot.num=5",
      "))",
      "out.cd <- suppressWarnings(plot(",
      "  bw_cd,",
      "  xdat=xd,",
      "  ydat=yd,",
      "  plot.behavior='data',",
      "  plot.errors.method='bootstrap',",
      "  plot.errors.boot.method='inid',",
      "  plot.errors.boot.num=3",
      "))",
      "stopifnot(is.list(out.ud), length(out.ud) > 0)",
      "stopifnot(is.list(out.cd), length(out.cd) > 0)",
      "cat('SESSION_INID_DENSITY_PLOT_OK\\n')",
      "cat('SESSION_INID_DENSITY_PLOT_DONE\\n')"
    ),
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_INID_DENSITY_PLOT_OK", res$output, fixed = TRUE)),
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
    args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
    timeout = 90L,
    env = c(
      env_common,
      "R_PROFILE_USER=",
      "R_PROFILE=",
      "FI_TCP_IFACE=en0",
      "FI_PROVIDER=tcp",
      "FI_SOCKETS_IFACE=en0"
    )
  )
  if (res$status != 0L) {
    res <- run_cmd_subprocess(
      mpiexec,
      args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
      timeout = 90L,
      env = c(
        env_common,
        "R_PROFILE_USER=",
        "R_PROFILE=",
        "FI_TCP_IFACE=lo0",
        "FI_PROVIDER=tcp",
        "FI_SOCKETS_IFACE=lo0"
      )
    )
  }

  if (res$status != 0L && .is_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for attach-mode smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("ATTACH_ROUTE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("profile mode smoke completes under mpiexec when enabled", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_PROFILE_TEST"), "1"),
              "set NP_RMPI_ENABLE_PROFILE_TEST=1 to run profile-mode smoke")
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  env_common <- subprocess_env()
  skip_if(is.null(env_common), "local npRmpi install unavailable for subprocess smoke")

  lib.path <- ensure_subprocess_npRmpi_lib()
  skip_if(is.null(lib.path), "local npRmpi install unavailable for subprocess smoke")
  profile.path <- file.path(lib.path, "npRmpi", "Rprofile")
  skip_if(!file.exists(profile.path), "npRmpi profile template unavailable in subprocess lib")

  script <- tempfile("npRmpi-profile-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "if (mpi.comm.rank(0L) == 0L) {",
    "  mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)",
    "  mpi.bcast.cmd(options(np.messages=FALSE), caller.execute=TRUE)",
    "  set.seed(42)",
    "  n <- 120",
    "  x <- runif(n)",
    "  y <- rnorm(n)",
    "  d <- data.frame(x=x, y=y)",
    "  mpi.bcast.Robj2slave(d)",
    "  mpi.bcast.cmd(bw <- npregbw(y~x, data=d, regtype='lc', bwmethod='cv.ls', nmulti=1), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit <- npreg(bws=bw, data=d, gradients=FALSE), caller.execute=TRUE)",
    "  stopifnot(inherits(fit, 'npregression'))",
    "  cat('PROFILE_ROUTE_OK\\n')",
    "  mpi.bcast.cmd(mpi.quit(), caller.execute=TRUE)",
    "}"
  ), script, useBytes = TRUE)

  env_profile <- c(
    env_common,
    sprintf("R_PROFILE_USER=%s", profile.path),
    "R_PROFILE=",
    "NP_RMPI_PROFILE_RECV_TIMEOUT_SEC=90",
    "FI_TCP_IFACE=en0",
    "FI_PROVIDER=tcp",
    "FI_SOCKETS_IFACE=en0"
  )

  res <- run_cmd_subprocess(
    mpiexec,
    args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
    timeout = 90L,
    env = env_profile
  )
  if (res$status != 0L) {
    env_profile_fallback <- c(
      env_common,
      sprintf("R_PROFILE_USER=%s", profile.path),
      "R_PROFILE=",
      "NP_RMPI_PROFILE_RECV_TIMEOUT_SEC=90",
      "FI_TCP_IFACE=lo0",
      "FI_PROVIDER=tcp",
      "FI_SOCKETS_IFACE=lo0"
    )
    res <- run_cmd_subprocess(
      mpiexec,
      args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
      timeout = 90L,
      env = env_profile_fallback
    )
  }

  if (res$status != 0L && .is_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for profile-mode smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("PROFILE_ROUTE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
