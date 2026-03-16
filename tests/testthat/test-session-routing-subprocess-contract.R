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

test_that("session core trio smoke completes in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "set.seed(77)",
      "n <- 80",
      "x <- runif(n)",
      "z <- runif(n)",
      "y <- sin(2*pi*x) + 0.5*z + rnorm(n, sd=0.1)",
      "d1 <- data.frame(x=x)",
      "d2 <- data.frame(x=x, z=z)",
      "dz <- data.frame(z=z)",
      "bw.sc <- npscoefbw(xdat=d1, ydat=y, zdat=dz, regtype='lc', nmulti=1)",
      "fit.sc <- npscoef(bws=bw.sc, gradients=FALSE)",
      "bw.pl <- npplregbw(xdat=d1, ydat=y, zdat=dz, regtype='lc', nmulti=1)",
      "fit.pl <- npplreg(bws=bw.pl, gradients=FALSE)",
      "bw.si <- npindexbw(xdat=d2, ydat=y, regtype='lc', nmulti=1)",
      "fit.si <- npindex(bws=bw.si, gradients=FALSE)",
      "stopifnot(inherits(fit.sc, 'smoothcoefficient'))",
      "stopifnot(inherits(fit.pl, 'plregression'))",
      "stopifnot(inherits(fit.si, 'singleindex'))",
      "cat('SESSION_CORE_TRIO_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_CORE_TRIO_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session generalized-nn shared degree-1 route stays exact after public hats", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(20260308)",
      "n <- 80L",
      "z <- sort(runif(n))",
      "x <- runif(n)",
      "y <- sin(2*pi*z) + 0.5*x + rnorm(n, sd=0.03)",
      "tx <- data.frame(x=x)",
      "tz <- data.frame(z=z)",
      "ez <- data.frame(z=seq(0.05, 0.95, length.out=24L))",
      "bw.reg.ll <- npregbw(xdat=tz, ydat=y, regtype='ll', bwtype='generalized_nn', bws=9, bandwidth.compute=FALSE)",
      "bw.reg.lp <- npregbw(xdat=tz, ydat=y, regtype='lp', basis='glp', degree=1L, bernstein.basis=FALSE, bwtype='generalized_nn', bws=9, bandwidth.compute=FALSE)",
      "fit.reg.ll <- npreg(bws=bw.reg.ll, txdat=tz, tydat=y, exdat=ez, warn.glp.gradient=FALSE)",
      "fit.reg.lp <- npreg(bws=bw.reg.lp, txdat=tz, tydat=y, exdat=ez, warn.glp.gradient=FALSE)",
      "hat.reg.apply.ll <- npreghat(bws=bw.reg.ll, txdat=tz, exdat=ez, y=y, output='apply')",
      "hat.reg.apply.lp <- npreghat(bws=bw.reg.lp, txdat=tz, exdat=ez, y=y, output='apply')",
      "hat.reg.matrix.ll <- npreghat(bws=bw.reg.ll, txdat=tz, exdat=ez)",
      "hat.reg.matrix.lp <- npreghat(bws=bw.reg.lp, txdat=tz, exdat=ez)",
      "bws.pl <- matrix(c(2, 9), nrow=2L, ncol=1L)",
      "bw.pl.ll <- npplregbw(xdat=tx, ydat=y, zdat=tz, regtype='ll', bwtype='generalized_nn', bws=bws.pl, bandwidth.compute=FALSE)",
      "bw.pl.lp <- npplregbw(xdat=tx, ydat=y, zdat=tz, regtype='lp', basis='glp', degree=1L, bernstein.basis=FALSE, bwtype='generalized_nn', bws=bws.pl, bandwidth.compute=FALSE)",
      "fit.pl.ll <- npplreg(bws=bw.pl.ll, txdat=tx, tzdat=tz, tydat=y)",
      "fit.pl.lp <- npplreg(bws=bw.pl.lp, txdat=tx, tzdat=tz, tydat=y)",
      "hat.pl.apply.ll <- npplreghat(bws=bw.pl.ll, txdat=tx, tzdat=tz, y=y, output='apply')",
      "hat.pl.apply.lp <- npplreghat(bws=bw.pl.lp, txdat=tx, tzdat=tz, y=y, output='apply')",
      "hat.pl.matrix.ll <- npplreghat(bws=bw.pl.ll, txdat=tx, tzdat=tz, output='matrix')",
      "hat.pl.matrix.lp <- npplreghat(bws=bw.pl.lp, txdat=tx, tzdat=tz, output='matrix')",
      "tol.direct <- 1e-9",
      "tol.matrix <- 1e-9",
      "checks <- c(",
      "  max(abs(as.numeric(fit.reg.ll$mean) - as.numeric(fit.reg.lp$mean))),",
      "  max(abs(as.numeric(hat.reg.apply.ll) - as.numeric(hat.reg.apply.lp))),",
      "  max(abs(as.numeric(hat.reg.apply.ll) - as.numeric(hat.reg.matrix.ll %*% y))),",
      "  max(abs(as.numeric(hat.reg.apply.lp) - as.numeric(hat.reg.matrix.lp %*% y))),",
      "  max(abs(as.numeric(fit.pl.ll$mean) - as.numeric(fit.pl.lp$mean))),",
      "  max(abs(as.numeric(fit.pl.ll$xcoef) - as.numeric(fit.pl.lp$xcoef))),",
      "  max(abs(as.numeric(hat.pl.apply.ll) - as.numeric(hat.pl.apply.lp))),",
      "  max(abs(as.numeric(hat.pl.apply.ll) - as.numeric(hat.pl.matrix.ll %*% y))),",
      "  max(abs(as.numeric(hat.pl.apply.lp) - as.numeric(hat.pl.matrix.lp %*% y)))",
      ")",
      "stopifnot(all(checks <= tol.direct))",
      "stopifnot(max(abs(as.numeric(hat.pl.matrix.ll %*% y) - as.numeric(hat.pl.matrix.lp %*% y))) <= tol.matrix)",
      "cat('SESSION_SHARED_DEGREE1_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_SHARED_DEGREE1_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session smooth-coefficient ll coef plot-data route completes in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(105)",
      "n <- 60",
      "x <- runif(n)",
      "z <- runif(n, -2, 2)",
      "y <- x * exp(z) * (1 + rnorm(n, sd=0.15))",
      "fit <- npscoef(y ~ x | z, regtype='ll', betas=TRUE)",
      "pdf(file=tempfile(fileext='.pdf'))",
      "on.exit(dev.off(), add=TRUE)",
      "out <- suppressWarnings(plot(",
      "  fit,",
      "  coef=TRUE,",
      "  coef.index=1,",
      "  perspective=FALSE,",
      "  neval=20,",
      "  plot.behavior='plot-data',",
      "  plot.errors.method='none'))",
      "stopifnot(is.list(out))",
      "stopifnot(length(out) > 0L)",
      "stopifnot(all(vapply(out, inherits, logical(1), 'smoothcoefficient')))",
      "cat('SESSION_SCOEF_LL_PLOTDATA_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_SCOEF_LL_PLOTDATA_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindex ichimura plain plot-data returns gradient-bearing payload", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(20260310)",
      "n <- 70L",
      "x1 <- runif(n)",
      "x2 <- runif(n)",
      "y <- sin(x1 + x2) + rnorm(n, sd=0.1)",
      "tx <- data.frame(x1=x1, x2=x2)",
      "bw <- npindexbw(xdat=tx, ydat=y, method='ichimura', regtype='lc', bwtype='fixed', bws=c(1,1,0.85), bandwidth.compute=FALSE)",
      "out <- suppressWarnings(plot(bw, xdat=tx, ydat=y, plot.behavior='data', plot.errors.method='none'))",
      "stopifnot(length(out) == 1L)",
      "stopifnot(all(c('index', 'mean', 'grad') %in% names(out[[1]])))",
      "cat('SESSION_NPINDEX_ICH_PLAIN_PAYLOAD_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_ICH_PLAIN_PAYLOAD_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindex ichimura gradient plot-data returns serial payload contract", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(20260310)",
      "n <- 70L",
      "x1 <- runif(n)",
      "x2 <- runif(n)",
      "y <- sin(x1 + x2) + rnorm(n, sd=0.1)",
      "tx <- data.frame(x1=x1, x2=x2)",
      "bw <- npindexbw(xdat=tx, ydat=y, method='ichimura', regtype='ll', bwtype='fixed', bws=c(1,1,0.85), bandwidth.compute=FALSE)",
      "out <- suppressWarnings(plot(bw, xdat=tx, ydat=y, plot.behavior='data', gradients=TRUE))",
      "stopifnot(length(out) == 1L)",
      "needed <- c('index', 'mean', 'merr', 'grad', 'gerr', 'mean.grad', 'mean.gerr', 'gradients')",
      "stopifnot(all(needed %in% names(out[[1]])))",
      "cat('SESSION_NPINDEX_ICH_GRADIENT_PAYLOAD_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_ICH_GRADIENT_PAYLOAD_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindex ichimura gradient lc fixed plot-data completes locally", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(20260310)",
      "n <- 5L",
      "x1 <- runif(n)",
      "x2 <- runif(n)",
      "y <- sin(x1 + x2) + rnorm(n, sd=0.1)",
      "tx <- data.frame(x1=x1, x2=x2)",
      "bw <- npindexbw(xdat=tx, ydat=y, method='ichimura', regtype='lc', bwtype='fixed', bws=c(1,1,0.85), bandwidth.compute=FALSE)",
      "out <- suppressWarnings(plot(bw, xdat=tx, ydat=y, plot.behavior='data', gradients=TRUE))",
      "needed <- c('index', 'mean', 'merr', 'grad', 'gerr', 'mean.grad', 'mean.gerr', 'gradients')",
      "stopifnot(length(out) == 1L)",
      "stopifnot(all(needed %in% names(out[[1]])))",
      "cat('SESSION_NPINDEX_ICH_GRADIENT_LC_FIXED_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_ICH_GRADIENT_LC_FIXED_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npplreg plain fixed lc plot-data completes locally", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(20260310)",
      "n <- 60L",
      "x <- runif(n)",
      "z <- runif(n)",
      "y <- sin(2*pi*x) + 0.5*z + rnorm(n, sd=0.05)",
      "tx <- data.frame(z=z)",
      "tz <- data.frame(x=x)",
      "bw.fix <- npplregbw(xdat=tx, zdat=tz, ydat=y, regtype='lc', bwmethod='cv.ls', bwtype='fixed', nmulti=1)",
      "bw.ann <- npplregbw(xdat=tx, zdat=tz, ydat=y, regtype='lc', bwmethod='cv.ls', bwtype='adaptive_nn', nmulti=1)",
      "data.out <- suppressWarnings(plot(bw.fix, xdat=tx, ydat=y, zdat=tz, plot.behavior='data', perspective=FALSE, plot.errors.method='none'))",
      "plot.out <- suppressWarnings(plot(bw.fix, xdat=tx, ydat=y, zdat=tz, plot.behavior='plot-data', perspective=FALSE, plot.errors.method='none'))",
      "ann.out <- suppressWarnings(plot(bw.ann, xdat=tx, ydat=y, zdat=tz, plot.behavior='data', perspective=FALSE, plot.errors.method='none'))",
      "stopifnot(is.list(data.out), length(data.out) == 2L, identical(names(data.out), names(plot.out)))",
      "for (i in seq_along(data.out)) {",
      "  stopifnot(inherits(data.out[[i]], 'plregression'))",
      "  stopifnot(inherits(plot.out[[i]], 'plregression'))",
      "  stopifnot(isTRUE(all.equal(data.out[[i]]$mean, plot.out[[i]]$mean, check.attributes=FALSE)))",
      "  stopifnot(isTRUE(all.equal(data.out[[i]]$merr, plot.out[[i]]$merr, check.attributes=FALSE)))",
      "  stopifnot(all(is.na(data.out[[i]]$merr)))",
      "}",
      "stopifnot(identical(as.character(data.out[[1L]]$ptype), 'Fixed'))",
      "stopifnot(identical(as.character(ann.out[[1L]]$ptype), 'Adaptive Nearest Neighbour'))",
      "cat('SESSION_NPPLREG_PLAIN_FIXED_LC_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPPLREG_PLAIN_FIXED_LC_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npsigtest fast-fail contract completes in installed-build subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "options(npRmpi.spmd.timeout.default=2)",
      "set.seed(19)",
      "n <- 20",
      "x1 <- runif(n)",
      "x2 <- runif(n)",
      "y <- rnorm(n)",
      "d <- data.frame(y=y, x1=x1, x2=x2)",
      "bw <- npregbw(y~x1+x2, data=d, bws=c(0.2, 0.4), bandwidth.compute=FALSE)",
      "err <- try(npsigtest(bws=bw, boot.num=9, random.seed=13, index=0), silent=TRUE)",
      "stopifnot(inherits(err, 'try-error'))",
      "msg <- as.character(err)",
      "stopifnot(any(grepl('invalid index provided', msg, fixed=TRUE)))",
      "stopifnot(!any(grepl(\"object 'd' not found\", msg, fixed=TRUE)))",
      "cat('SESSION_NPSIGTEST_SUBPROCESS_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPSIGTEST_SUBPROCESS_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindex residual and evaluation-error branches complete in installed-build subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(20260309)",
      "n <- 56L",
      "x1 <- runif(n)",
      "x2 <- runif(n)",
      "y <- sin(x1 + x2) + rnorm(n, sd=0.05)",
      "tx <- data.frame(x1=x1, x2=x2)",
      "ex <- tx[seq_len(18L), , drop=FALSE]",
      "ey <- sin(ex$x1 + ex$x2) + rnorm(nrow(ex), sd=0.05)",
      "cfgs <- list(",
      "  list(regtype='lc'),",
      "  list(regtype='ll'),",
      "  list(regtype='lp', basis='glp', degree=1L)",
      ")",
      "for (cfg in cfgs) {",
      "  bw.args <- list(xdat=tx, ydat=y, regtype=cfg$regtype, bwtype='adaptive_nn', nmulti=1L)",
      "  if (!is.null(cfg$basis)) {",
      "    bw.args$basis <- cfg$basis",
      "    bw.args$degree <- cfg$degree",
      "    bw.args$bernstein.basis <- FALSE",
      "  }",
      "  bw <- do.call(npindexbw, bw.args)",
      "  fit.is <- npindex(bws=bw, txdat=tx, tydat=y, gradients=FALSE, residuals=TRUE)",
      "  stopifnot(length(fit.is$resid) == nrow(tx))",
      "  stopifnot(max(abs(as.vector(fit.is$resid) - (y - as.vector(fit.is$mean)))) < 1e-8)",
      "  fit.oos <- npindex(bws=bw, txdat=tx, tydat=y, exdat=ex, eydat=ey, gradients=FALSE, errors=TRUE)",
      "  stopifnot(length(fit.oos$mean) == nrow(ex))",
      "  stopifnot(length(fit.oos$merr) == nrow(ex))",
      "  stopifnot(all(is.finite(fit.oos$mean)))",
      "  stopifnot(all(is.finite(fit.oos$merr)))",
      "  fit.grad <- npindex(bws=bw, txdat=tx, tydat=y, exdat=ex, eydat=ey, gradients=TRUE, errors=TRUE)",
      "  stopifnot(identical(dim(fit.grad$grad), c(nrow(ex), ncol(tx))))",
      "  stopifnot(identical(dim(fit.grad$gerr), c(nrow(ex), ncol(tx))))",
      "  stopifnot(all(is.finite(fit.grad$grad)))",
      "  stopifnot(all(is.finite(fit.grad$gerr)))",
      "  cat(sprintf('SESSION_NPINDEX_BRANCH regtype=%s ok\\n', cfg$regtype))",
      "}",
      "cat('SESSION_NPINDEX_BRANCH_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_BRANCH_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session quit/init cycle resets SPMD sequence for LL CV calls", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "mkdat <- function(seed) {",
      "  set.seed(seed)",
      "  n <- 220L",
      "  x1 <- runif(n); x2 <- runif(n)",
      "  z1 <- factor(rbinom(n, 1L, 0.5))",
      "  z2 <- ordered(rbinom(n, 1L, 0.5))",
      "  y <- cos(2*pi*x1) + 0.5*sin(2*pi*x2) + as.numeric(z1) + rnorm(n, sd=0.2)",
      "  data.frame(y=y, x1=x1, x2=x2, z1=z1, z2=z2)",
      "}",
      "runll <- function(seed) {",
      "  d <- mkdat(seed)",
      "  set.seed(321)",
      "  bw <- npregbw(y~x1+x2+z1+z2, regtype='ll', bwmethod='cv.ls', nmulti=1, data=d)",
      "  stopifnot(is.finite(as.numeric(bw$fval)))",
      "  invisible(bw)",
      "}",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "runll(100)",
      "npRmpi.quit(mode='spawn', force=TRUE)",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "runll(101)",
      "stopifnot(identical(as.integer(getOption('npRmpi.spmd.seq_id')), 1L))",
      "npRmpi.quit(mode='spawn', force=TRUE)",
      "cat('SESSION_SEQ_RESET_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_SEQ_RESET_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session soft close/reopen keeps reused slave pool SPMD-synchronized", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE, npRmpi.reuse.slaves=TRUE)",
      "run_lp <- function(seed) {",
      "  set.seed(seed)",
      "  n <- 200L",
      "  x <- runif(n)",
      "  y <- rnorm(n, sd=0.5*sd(x))",
      "  fit <- npreg(y~x, regtype='lp', degree=2L)",
      "  stopifnot(is.finite(as.numeric(fit$bws$fval)))",
      "  invisible(fit)",
      "}",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "run_lp(42)",
      "npRmpi.quit(mode='spawn', force=FALSE)",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "run_lp(43)",
      "npRmpi.quit(mode='spawn', force=TRUE)",
      "cat('SESSION_SOFT_REOPEN_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_SOFT_REOPEN_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session rejects nslaves=0 with serial-workflow remediation", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "err <- try(npRmpi.init(nslaves=0, quiet=TRUE), silent=TRUE)",
      "stopifnot(inherits(err, 'try-error'))",
      "cat(as.character(err), '\\n')",
      "cat('SESSION_NSLAVES0_REJECT_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NSLAVES0_REJECT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("must be >= 1", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("use package 'np' for serial workflows", res$output, fixed = TRUE)),
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

test_that("session nearest-neighbor npreg formula routes complete with summary and plot", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(11)",
      "n <- 100",
      "x <- runif(n)",
      "y <- x^2 + rnorm(n, sd=0.1)",
      "g.gen <- npreg(y~x, regtype='ll', bwtype='generalized_nn')",
      "g.adp <- npreg(y~x, regtype='ll', bwtype='adaptive_nn')",
      "summary(g.gen)",
      "summary(g.adp)",
      "pdf(tempfile(fileext='.pdf'))",
      "plot(g.gen)",
      "plot(g.adp)",
      "dev.off()",
      "stopifnot(inherits(g.gen, 'npregression'), inherits(g.adp, 'npregression'))",
      "stopifnot(length(g.gen$mean) == n, length(g.adp$mean) == n)",
      "cat('SESSION_NPREG_NN_FORMULA_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPREG_NN_FORMULA_OK", res$output, fixed = TRUE)),
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

test_that("session adaptive-nn npreghat matrix owner stays exact in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "ns <- asNamespace('npRmpi')",
      "reg_direct <- get('.np_regression_direct', envir=ns, inherits=FALSE)",
      "with_local <- get('.npRmpi_with_local_regression', envir=ns, inherits=FALSE)",
      "set.seed(20260309)",
      "n <- 70L",
      "x <- sort(runif(n))",
      "y <- sin(2*pi*x) + 0.25*x + rnorm(n, sd=0.04)",
      "tx <- data.frame(x=x)",
      "ex <- data.frame(x=seq(0.05, 0.95, length.out=20L))",
      "probe_one <- function(regtype, selected) {",
      "  bw_args <- list(xdat=tx, ydat=y, regtype=regtype, bwtype='adaptive_nn')",
      "  if (selected) {",
      "    bw_args$nmulti <- 1L",
      "  } else {",
      "    bw_args$bandwidth.compute <- FALSE",
      "    bw_args$bws <- 9L",
      "  }",
      "  if (identical(regtype, 'lp')) {",
      "    bw_args$degree <- 1L",
      "    bw_args$basis <- 'glp'",
      "    bw_args$bernstein.basis <- FALSE",
      "  }",
      "  bw <- do.call(npregbw, bw_args)",
      "  direct <- with_local(reg_direct(",
      "    bws=bw, txdat=tx, tydat=y, exdat=ex,",
      "    gradients=!identical(regtype, 'lc'), gradient.order=1L",
      "  ))",
      "  hat.apply <- npreghat(bws=bw, txdat=tx, exdat=ex, y=y, output='apply')",
      "  hat.matrix <- npreghat(bws=bw, txdat=tx, exdat=ex, output='matrix')",
      "  stopifnot(max(abs(as.vector(hat.apply) - as.vector(direct$mean))) < 1e-8)",
      "  stopifnot(max(abs(as.vector(hat.matrix %*% y) - as.vector(hat.apply))) < 1e-8)",
      "  if (!identical(regtype, 'lc')) {",
      "    grad.apply <- npreghat(bws=bw, txdat=tx, exdat=ex, y=y, output='apply', s=1L)",
      "    grad.matrix <- npreghat(bws=bw, txdat=tx, exdat=ex, output='matrix', s=1L)",
      "    stopifnot(max(abs(as.vector(grad.apply) - as.vector(direct$grad[,1]))) < 1e-6)",
      "    stopifnot(max(abs(as.vector(grad.matrix %*% y) - as.vector(grad.apply))) < 1e-6)",
      "  }",
      "  cat(sprintf('SESSION_NPREGHAT_ADAPTIVE_OWNER route=%s regtype=%s bw=%s\\n', if (selected) 'selected' else 'manual', regtype, paste(bw$bw, collapse=',')))",
      "}",
      "for (selected in c(FALSE, TRUE)) for (regtype in c('lc', 'll', 'lp')) probe_one(regtype, selected)",
      "cat('SESSION_NPREGHAT_ADAPTIVE_OWNER_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPREGHAT_ADAPTIVE_OWNER_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session adaptive-nn npscoefhat selected owner preserves integer support in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(20260308)",
      "n <- 50L",
      "x <- runif(n)",
      "z <- runif(n)",
      "y <- (0.4 + x) * sin(2*pi*z) + rnorm(n, sd=0.04)",
      "tx <- data.frame(x=x)",
      "tz <- data.frame(z=z)",
      "ex <- data.frame(x=seq(0.1, 0.9, length.out=12L))",
      "ez <- data.frame(z=seq(0.1, 0.9, length.out=12L))",
      "tol <- sqrt(.Machine$double.eps)",
      "upper <- n - 1L",
      "for (regtype in c('lc', 'll', 'lp')) {",
      "  bw_args <- list(xdat=tx, zdat=tz, ydat=y, regtype=regtype, bwtype='adaptive_nn', nmulti=1L)",
      "  if (identical(regtype, 'lp')) {",
      "    bw_args$degree <- 1L",
      "    bw_args$basis <- 'glp'",
      "    bw_args$bernstein.basis <- FALSE",
      "  }",
      "  bw <- do.call(npscoefbw, bw_args)",
      "  stopifnot(all(abs(bw$bw - round(bw$bw)) <= tol))",
      "  stopifnot(all(bw$bw >= 1 & bw$bw <= upper))",
      "  hat.apply <- npscoefhat(bws=bw, txdat=tx, tzdat=tz, exdat=ex, ezdat=ez, y=y, output='apply', iterate=FALSE)",
      "  hat.matrix <- npscoefhat(bws=bw, txdat=tx, tzdat=tz, exdat=ex, ezdat=ez, output='matrix', iterate=FALSE)",
      "  stopifnot(max(abs(as.vector(hat.apply) - as.vector(hat.matrix %*% y))) < 1e-8)",
      "  cat(sprintf('SESSION_NPSCOEF_ADAPTIVE_OWNER regtype=%s bw=%s\\n', regtype, paste(bw$bw, collapse=',')))",
      "}",
      "err <- try(npscoefbw(xdat=tx, zdat=tz, ydat=y, regtype='lc', bwtype='adaptive_nn', bandwidth.compute=FALSE, bws=0.13), silent=TRUE)",
      "stopifnot(inherits(err, 'try-error'))",
      "cat('SESSION_NPSCOEF_ADAPTIVE_OWNER_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPSCOEF_ADAPTIVE_OWNER_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindexhat adaptive-nn exact owner route completes in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(105)",
      "n <- 60L",
      "x <- runif(n)",
      "x2 <- runif(n)",
      "y <- x * exp(x2) * (1 + rnorm(n, sd=0.15))",
      "idxdat <- data.frame(y=y, x=x, x2=x2)",
      "fit <- npindex(y ~ x + x2, data=idxdat, method='ichimura', bwtype='adaptive_nn')",
      "hy <- npindexhat(bws=fit$bws, txdat=idxdat[c('x','x2')], exdat=idxdat[c('x','x2')], y=y, output='apply', s=0L)",
      "stopifnot(length(hy) == n)",
      "stopifnot(max(abs(as.vector(hy) - as.vector(fit$mean))) < 1e-8)",
      "cat('SESSION_NPINDEXHAT_ADAPTIVE_EXACT_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEXHAT_ADAPTIVE_EXACT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindexhat adaptive-nn manual owner control stays exact in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(314161)",
      "n <- 70L",
      "x1 <- runif(n)",
      "x2 <- runif(n)",
      "y <- sin(x1 + x2) + rnorm(n, sd=0.06)",
      "tx <- data.frame(x1=x1, x2=x2)",
      "ex <- tx[seq_len(20), , drop=FALSE]",
      "bw <- npindexbw(xdat=tx, ydat=y, bws=c(1, 1, 9), bandwidth.compute=FALSE, regtype='lc', bwtype='adaptive_nn')",
      "fit <- npindex(bws=bw, txdat=tx, tydat=y, exdat=ex, gradients=FALSE)",
      "hy <- npindexhat(bws=bw, txdat=tx, exdat=ex, y=y, output='apply', s=0L)",
      "stopifnot(length(hy) == nrow(ex))",
      "stopifnot(max(abs(as.vector(hy) - as.vector(fit$mean))) < 1e-8)",
      "cat('SESSION_NPINDEXHAT_ADAPTIVE_MANUAL_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEXHAT_ADAPTIVE_MANUAL_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindexhat adaptive-nn ll owner route stays exact in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(20260309)",
      "n <- 60L",
      "x1 <- runif(n)",
      "x2 <- runif(n)",
      "y <- sin(x1 + x2) + rnorm(n, sd=0.05)",
      "tx <- data.frame(x1=x1, x2=x2)",
      "ex <- tx[seq_len(20L), , drop=FALSE]",
      "bw <- npindexbw(xdat=tx, ydat=y, regtype='ll', bwtype='adaptive_nn', nmulti=1L)",
      "fit.mean <- npindex(bws=bw, txdat=tx, tydat=y, exdat=ex, gradients=FALSE)",
      "fit.grad <- npindex(bws=bw, txdat=tx, tydat=y, exdat=ex, gradients=TRUE)",
      "a.mean <- npindexhat(bws=bw, txdat=tx, exdat=ex, y=y, output='apply', s=0L)",
      "H.mean <- npindexhat(bws=bw, txdat=tx, exdat=ex, output='matrix', s=0L)",
      "a.grad <- npindexhat(bws=bw, txdat=tx, exdat=ex, y=y, output='apply', s=1L)",
      "H.grad <- npindexhat(bws=bw, txdat=tx, exdat=ex, output='matrix', s=1L)",
      "stopifnot(max(abs(as.vector(fit.mean$mean) - as.vector(a.mean))) < 1e-8)",
      "stopifnot(max(abs(as.vector(a.mean) - as.vector(H.mean %*% y))) < 1e-8)",
      "stopifnot(max(abs(as.vector(fit.grad$grad[,1L]) - as.vector(a.grad))) < 1e-8)",
      "stopifnot(max(abs(as.vector(a.grad) - as.vector(H.grad %*% y))) < 1e-8)",
      "cat('SESSION_NPINDEXHAT_LL_OWNER_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEXHAT_LL_OWNER_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindex adaptive-nn public route preserves bwtype semantics in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(314161)",
      "n <- 70L",
      "x1 <- runif(n)",
      "x2 <- runif(n)",
      "y <- sin(x1 + x2) + rnorm(n, sd=0.06)",
      "tx <- data.frame(x1=x1, x2=x2)",
      "ex <- tx[seq_len(20), , drop=FALSE]",
      "bw.fixed <- npindexbw(xdat=tx, ydat=y, bws=c(1, 1, 0.85), bandwidth.compute=FALSE, regtype='lc', bwtype='fixed')",
      "bw.adp <- npindexbw(xdat=tx, ydat=y, bws=c(1, 1, 9), bandwidth.compute=FALSE, regtype='lc', bwtype='adaptive_nn')",
      "fit.fixed <- npindex(bws=bw.fixed, txdat=tx, tydat=y, exdat=ex, gradients=FALSE)",
      "fit.adp <- npindex(bws=bw.adp, txdat=tx, tydat=y, exdat=ex, gradients=FALSE)",
      "stopifnot(max(abs(as.vector(fit.fixed$mean) - as.vector(fit.adp$mean))) > 1e-6)",
      "cat('SESSION_NPINDEX_ADAPTIVE_PUBLIC_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_ADAPTIVE_PUBLIC_OK", res$output, fixed = TRUE)),
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

test_that("session npindex inid consumer plot preserves bwtype variants in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE)",
      "set.seed(20260310)",
      "n <- 70",
      "x1 <- runif(n)",
      "x2 <- runif(n)",
      "eta <- 1.2 * x1 - 0.8 * x2",
      "p <- 1 / (1 + exp(-eta))",
      "y <- rbinom(n, 1, p)",
      "tx <- data.frame(x1=x1, x2=x2)",
      "for (bwtype in c('fixed', 'generalized_nn', 'adaptive_nn')) {",
      "  bw.arg <- if (identical(bwtype, 'fixed')) 0.85 else 9",
      "  bw <- npindexbw(",
      "    xdat=tx,",
      "    ydat=y,",
      "    method='kleinspady',",
      "    regtype='lc',",
      "    bwtype=bwtype,",
      "    bws=c(1, 1, bw.arg),",
      "    bandwidth.compute=FALSE",
      "  )",
      "  out <- suppressWarnings(plot(",
      "    bw,",
      "    xdat=tx,",
      "    ydat=y,",
      "    plot.behavior='data',",
      "    plot.errors.method='bootstrap',",
      "    plot.errors.boot.method='inid',",
      "    plot.errors.boot.num=9",
      "  ))",
      "  stopifnot(is.list(out), length(out) > 0L)",
      "  stopifnot(is.matrix(out[[1]]$merr), ncol(out[[1]]$merr) == 2L)",
      "  cat(sprintf('SESSION_NPINDEX_BWTYPE_INID_PLOT_OK bwtype=%s\\n', bwtype))",
      "}"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_BWTYPE_INID_PLOT_OK bwtype=fixed", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_BWTYPE_INID_PLOT_OK bwtype=generalized_nn", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_BWTYPE_INID_PLOT_OK bwtype=adaptive_nn", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindex fixed-block consumer plot preserves bwtype variants in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE)",
      "set.seed(20260310)",
      "n <- 70",
      "x1 <- runif(n)",
      "x2 <- runif(n)",
      "eta <- 1.2 * x1 - 0.8 * x2",
      "p <- 1 / (1 + exp(-eta))",
      "y <- rbinom(n, 1, p)",
      "tx <- data.frame(x1=x1, x2=x2)",
      "for (bwtype in c('fixed', 'generalized_nn', 'adaptive_nn')) {",
      "  bw.arg <- if (identical(bwtype, 'fixed')) 0.85 else 9",
      "  bw <- npindexbw(",
      "    xdat=tx,",
      "    ydat=y,",
      "    method='kleinspady',",
      "    regtype='lc',",
      "    bwtype=bwtype,",
      "    bws=c(1, 1, bw.arg),",
      "    bandwidth.compute=FALSE",
      "  )",
      "  out <- suppressWarnings(plot(",
      "    bw,",
      "    xdat=tx,",
      "    ydat=y,",
      "    plot.behavior='data',",
      "    plot.errors.method='bootstrap',",
      "    plot.errors.boot.method='fixed',",
      "    plot.errors.boot.blocklen=3,",
      "    plot.errors.boot.num=9",
      "  ))",
      "  stopifnot(is.list(out), length(out) > 0L)",
      "  stopifnot(is.matrix(out[[1]]$merr), ncol(out[[1]]$merr) == 2L)",
      "  cat(sprintf('SESSION_NPINDEX_BWTYPE_FIXED_PLOT_OK bwtype=%s\\n', bwtype))",
      "}"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_BWTYPE_FIXED_PLOT_OK bwtype=fixed", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_BWTYPE_FIXED_PLOT_OK bwtype=generalized_nn", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_BWTYPE_FIXED_PLOT_OK bwtype=adaptive_nn", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindex geom-block consumer plot preserves bwtype variants in subprocess", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE)",
      "set.seed(20260310)",
      "n <- 70",
      "x1 <- runif(n)",
      "x2 <- runif(n)",
      "eta <- 1.2 * x1 - 0.8 * x2",
      "p <- 1 / (1 + exp(-eta))",
      "y <- rbinom(n, 1, p)",
      "tx <- data.frame(x1=x1, x2=x2)",
      "for (bwtype in c('fixed', 'generalized_nn', 'adaptive_nn')) {",
      "  bw.arg <- if (identical(bwtype, 'fixed')) 0.85 else 9",
      "  bw <- npindexbw(",
      "    xdat=tx,",
      "    ydat=y,",
      "    method='kleinspady',",
      "    regtype='lc',",
      "    bwtype=bwtype,",
      "    bws=c(1, 1, bw.arg),",
      "    bandwidth.compute=FALSE",
      "  )",
      "  out <- suppressWarnings(plot(",
      "    bw,",
      "    xdat=tx,",
      "    ydat=y,",
      "    plot.behavior='data',",
      "    plot.errors.method='bootstrap',",
      "    plot.errors.boot.method='geom',",
      "    plot.errors.boot.blocklen=3,",
      "    plot.errors.boot.num=9",
      "  ))",
      "  stopifnot(is.list(out), length(out) > 0L)",
      "  stopifnot(is.matrix(out[[1]]$merr), ncol(out[[1]]$merr) == 2L)",
      "  cat(sprintf('SESSION_NPINDEX_BWTYPE_GEOM_PLOT_OK bwtype=%s\\n', bwtype))",
      "}"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_BWTYPE_GEOM_PLOT_OK bwtype=fixed", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_BWTYPE_GEOM_PLOT_OK bwtype=generalized_nn", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_BWTYPE_GEOM_PLOT_OK bwtype=adaptive_nn", res$output, fixed = TRUE)),
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
      "z <- runif(n)",
      "y <- sin(2*pi*x) + 0.5*z + rnorm(n, sd=0.1)",
      "d <- data.frame(x=x, y=y)",
      "d1 <- data.frame(x=x)",
      "d2 <- data.frame(x=x, z=z)",
      "dz <- data.frame(z=z)",
      "xd <- data.frame(u=x)",
      "yd <- data.frame(v=y)",
      "mpi.bcast.Robj2slave(d)",
      "mpi.bcast.Robj2slave(d1)",
      "mpi.bcast.Robj2slave(d2)",
      "mpi.bcast.Robj2slave(dz)",
      "mpi.bcast.Robj2slave(xd)",
      "mpi.bcast.Robj2slave(yd)",
      "mpi.bcast.Robj2slave(y)",
      "mpi.bcast.cmd(bw <- npregbw(y~x, data=d, regtype='lc', bwmethod='cv.ls', nmulti=1), caller.execute=TRUE)",
      "mpi.bcast.cmd(fit <- npreg(bws=bw, gradients=FALSE), caller.execute=TRUE)",
      "mpi.bcast.cmd(bw.sc <- npscoefbw(xdat=d1, ydat=y, zdat=dz, regtype='lc', nmulti=1), caller.execute=TRUE)",
      "mpi.bcast.cmd(fit.sc <- npscoef(bws=bw.sc, gradients=FALSE), caller.execute=TRUE)",
      "mpi.bcast.cmd(bw.pl <- npplregbw(xdat=d1, ydat=y, zdat=dz, regtype='lc', nmulti=1), caller.execute=TRUE)",
      "mpi.bcast.cmd(fit.pl <- npplreg(bws=bw.pl, gradients=FALSE), caller.execute=TRUE)",
      "mpi.bcast.cmd(bw.si <- npindexbw(xdat=d2, ydat=y, regtype='lc', nmulti=1), caller.execute=TRUE)",
      "mpi.bcast.cmd(fit.si <- npindex(bws=bw.si, gradients=FALSE), caller.execute=TRUE)",
      "mpi.bcast.cmd(bw.ud <- npudensbw(dat=xd, bws=0.35, bandwidth.compute=FALSE), caller.execute=TRUE)",
      "mpi.bcast.cmd(fit.ud <- npudens(tdat=xd, bws=bw.ud), caller.execute=TRUE)",
      "mpi.bcast.cmd(bw.uf <- npudistbw(dat=xd, bws=0.35, bandwidth.compute=FALSE), caller.execute=TRUE)",
      "mpi.bcast.cmd(fit.uf <- npudist(tdat=xd, bws=bw.uf), caller.execute=TRUE)",
      "mpi.bcast.cmd(bw.cd <- npcdensbw(xdat=xd, ydat=yd, bws=c(0.45,0.45), bandwidth.compute=FALSE), caller.execute=TRUE)",
      "mpi.bcast.cmd(fit.cd <- npcdens(txdat=xd, tydat=yd, bws=bw.cd), caller.execute=TRUE)",
      "mpi.bcast.cmd(bw.cf <- npcdistbw(xdat=xd, ydat=yd, bws=c(0.45,0.45), bandwidth.compute=FALSE), caller.execute=TRUE)",
      "mpi.bcast.cmd(fit.cf <- npcdist(txdat=xd, tydat=yd, bws=bw.cf), caller.execute=TRUE)",
      "stopifnot(inherits(fit, 'npregression'))",
      "stopifnot(inherits(fit.sc, 'smoothcoefficient'))",
      "stopifnot(inherits(fit.pl, 'plregression'))",
      "stopifnot(inherits(fit.si, 'singleindex'))",
      "stopifnot(inherits(fit.ud, 'npdensity'))",
      "stopifnot(inherits(fit.uf, 'npdistribution'))",
      "stopifnot(inherits(fit.cd, 'condensity'))",
      "stopifnot(inherits(fit.cf, 'condistribution'))",
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
    "  suppressPackageStartupMessages(library(MASS))",
    "  n <- 80",
    "  x <- runif(n)",
    "  z <- runif(n)",
    "  y <- sin(2*pi*x) + 0.5*z + rnorm(n, sd=0.1)",
    "  d1 <- data.frame(x=x)",
    "  d2 <- data.frame(x=x, z=z)",
    "  dz <- data.frame(z=z)",
    "  d.cop <- data.frame(x=x, y=y)",
    "  u.cop <- data.frame(x=c(0.25,0.5,0.75), y=c(0.25,0.5,0.75))",
    "  xd <- data.frame(u=x)",
    "  yd <- data.frame(v=y)",
    "  bw <- npregbw(y~x, regtype='lc', bwmethod='cv.ls', nmulti=1)",
    "  fit <- npreg(bws=bw, gradients=FALSE)",
    "  bw.sc <- npscoefbw(xdat=d1, ydat=y, zdat=dz, regtype='lc', nmulti=1)",
    "  fit.sc <- npscoef(bws=bw.sc, gradients=FALSE)",
    "  bw.pl <- npplregbw(xdat=d1, ydat=y, zdat=dz, regtype='lc', nmulti=1)",
    "  fit.pl <- npplreg(bws=bw.pl, gradients=FALSE)",
    "  bw.si <- npindexbw(xdat=d2, ydat=y, regtype='lc', nmulti=1)",
    "  fit.si <- npindex(bws=bw.si, gradients=FALSE)",
    "  bw.ud <- npudensbw(dat=xd, bws=0.35, bandwidth.compute=FALSE)",
    "  fit.ud <- npudens(tdat=xd, bws=bw.ud)",
    "  bw.uf <- npudistbw(dat=xd, bws=0.35, bandwidth.compute=FALSE)",
    "  fit.uf <- npudist(tdat=xd, bws=bw.uf)",
    "  bw.cd <- npcdensbw(xdat=xd, ydat=yd, bws=c(0.45,0.45), bandwidth.compute=FALSE)",
    "  fit.cd <- npcdens(txdat=xd, tydat=yd, bws=bw.cd)",
    "  bw.cf <- npcdistbw(xdat=xd, ydat=yd, bws=c(0.45,0.45), bandwidth.compute=FALSE)",
    "  fit.cf <- npcdist(txdat=xd, tydat=yd, bws=bw.cf)",
    "  bw.cop <- npudistbw(~x+y, data=d.cop)",
    "  cop <- npcopula(bws=bw.cop, data=d.cop, u=u.cop, n.quasi.inv=60)",
    "  data(birthwt)",
    "  bdat <- birthwt",
    "  bdat$low <- factor(bdat$low)",
    "  bdat$smoke <- factor(bdat$smoke)",
    "  bdat$race <- factor(bdat$race)",
    "  bdat$ht <- factor(bdat$ht)",
    "  bdat$ui <- factor(bdat$ui)",
    "  bdat$ftv <- ordered(bdat$ftv)",
    "  bw.cm <- npcdensbw(low~smoke+race+ht+ui+ftv+age+lwt, data=bdat, nmulti=1)",
    "  fit.cm <- npconmode(bws=bw.cm)",
    "  stopifnot(inherits(fit, 'npregression'))",
    "  stopifnot(inherits(fit.sc, 'smoothcoefficient'))",
    "  stopifnot(inherits(fit.pl, 'plregression'))",
    "  stopifnot(inherits(fit.si, 'singleindex'))",
    "  stopifnot(inherits(fit.ud, 'npdensity'))",
    "  stopifnot(inherits(fit.uf, 'npdistribution'))",
    "  stopifnot(inherits(fit.cd, 'condensity'))",
    "  stopifnot(inherits(fit.cf, 'condistribution'))",
    "  stopifnot(inherits(cop, 'data.frame'))",
    "  stopifnot(nrow(cop) == 9L)",
    "  stopifnot(all(is.finite(cop$copula)))",
    "  stopifnot(inherits(fit.cm, 'conmode'))",
    "  cat('ATTACH_NPCONMODE_ROUTE_OK\\n')",
    "  cat('ATTACH_NPCOPULA_ROUTE_OK\\n')",
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
  expect_true(any(grepl("ATTACH_NPCONMODE_ROUTE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("ATTACH_NPCOPULA_ROUTE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("attach npindexhat adaptive-nn exact owner route completes under mpiexec when enabled", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_TEST=1 to run attach-mode smoke")
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  script <- tempfile("npRmpi-attach-npindexhat-exact-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "is.master <- isTRUE(npRmpi.init(mode='attach', quiet=TRUE, autodispatch=TRUE))",
    "on.exit({",
    "  try(npRmpi.quit(mode='attach'), silent=TRUE)",
    "  if (isTRUE(is.master)) try(Rmpi::mpi.quit(), silent=TRUE)",
    "}, add=TRUE)",
    "options(np.messages=FALSE)",
    "if (isTRUE(is.master)) {",
    "  set.seed(105)",
    "  n <- 60L",
    "  x <- runif(n)",
    "  x2 <- runif(n)",
    "  y <- x * exp(x2) * (1 + rnorm(n, sd=0.15))",
    "  idxdat <- data.frame(y=y, x=x, x2=x2)",
    "  fit <- npindex(y ~ x + x2, data=idxdat, method='ichimura', bwtype='adaptive_nn')",
    "  hy <- npindexhat(bws=fit$bws, txdat=idxdat[c('x','x2')], exdat=idxdat[c('x','x2')], y=y, output='apply', s=0L)",
    "  stopifnot(length(hy) == n)",
    "  stopifnot(max(abs(as.vector(hy) - as.vector(fit$mean))) < 1e-8)",
    "  cat('ATTACH_NPINDEXHAT_ADAPTIVE_EXACT_OK\\n')",
    "}"
  ), script, useBytes = TRUE)

  env_common <- subprocess_env()
  skip_if(is.null(env_common), "local npRmpi install unavailable for subprocess smoke")
  res <- run_cmd_subprocess(
    mpiexec,
    args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
    timeout = 120L,
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
      timeout = 120L,
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
  expect_true(any(grepl("ATTACH_NPINDEXHAT_ADAPTIVE_EXACT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("attach npindexhat adaptive-nn ll owner route stays exact under mpiexec when enabled", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_TEST=1 to run attach-mode smoke")
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  script <- tempfile("npRmpi-attach-npindexhat-ll-owner-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "is.master <- isTRUE(npRmpi.init(mode='attach', quiet=TRUE, autodispatch=TRUE))",
    "on.exit({",
    "  try(npRmpi.quit(mode='attach'), silent=TRUE)",
    "  if (isTRUE(is.master)) try(Rmpi::mpi.quit(), silent=TRUE)",
    "}, add=TRUE)",
    "options(np.messages=FALSE)",
    "if (isTRUE(is.master)) {",
    "  set.seed(20260309)",
    "  n <- 60L",
    "  x1 <- runif(n)",
    "  x2 <- runif(n)",
    "  y <- sin(x1 + x2) + rnorm(n, sd=0.05)",
    "  tx <- data.frame(x1=x1, x2=x2)",
    "  ex <- tx[seq_len(20L), , drop=FALSE]",
    "  bw <- npindexbw(xdat=tx, ydat=y, regtype='ll', bwtype='adaptive_nn', nmulti=1L)",
    "  fit.mean <- npindex(bws=bw, txdat=tx, tydat=y, exdat=ex, gradients=FALSE)",
    "  fit.grad <- npindex(bws=bw, txdat=tx, tydat=y, exdat=ex, gradients=TRUE)",
    "  a.mean <- npindexhat(bws=bw, txdat=tx, exdat=ex, y=y, output='apply', s=0L)",
    "  H.mean <- npindexhat(bws=bw, txdat=tx, exdat=ex, output='matrix', s=0L)",
    "  a.grad <- npindexhat(bws=bw, txdat=tx, exdat=ex, y=y, output='apply', s=1L)",
    "  H.grad <- npindexhat(bws=bw, txdat=tx, exdat=ex, output='matrix', s=1L)",
    "  stopifnot(max(abs(as.vector(fit.mean$mean) - as.vector(a.mean))) < 1e-8)",
    "  stopifnot(max(abs(as.vector(a.mean) - as.vector(H.mean %*% y))) < 1e-8)",
    "  stopifnot(max(abs(as.vector(fit.grad$grad[,1L]) - as.vector(a.grad))) < 1e-8)",
    "  stopifnot(max(abs(as.vector(a.grad) - as.vector(H.grad %*% y))) < 1e-8)",
    "  cat('ATTACH_NPINDEXHAT_LL_OWNER_OK\\n')",
    "}"
  ), script, useBytes = TRUE)

  env_common <- subprocess_env()
  skip_if(is.null(env_common), "local npRmpi install unavailable for subprocess smoke")
  res <- run_cmd_subprocess(
    mpiexec,
    args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
    timeout = 120L,
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
      timeout = 120L,
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
  expect_true(any(grepl("ATTACH_NPINDEXHAT_LL_OWNER_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("attach npindex adaptive-nn public route preserves bwtype semantics under mpiexec when enabled", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_TEST=1 to run attach-mode smoke")
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  script <- tempfile("npRmpi-attach-npindex-public-adaptive-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "is.master <- isTRUE(npRmpi.init(mode='attach', quiet=TRUE, autodispatch=TRUE))",
    "on.exit({",
    "  try(npRmpi.quit(mode='attach'), silent=TRUE)",
    "  if (isTRUE(is.master)) try(Rmpi::mpi.quit(), silent=TRUE)",
    "}, add=TRUE)",
    "options(np.messages=FALSE)",
    "if (isTRUE(is.master)) {",
    "  set.seed(314161)",
    "  n <- 70L",
    "  x1 <- runif(n)",
    "  x2 <- runif(n)",
    "  y <- sin(x1 + x2) + rnorm(n, sd=0.06)",
    "  tx <- data.frame(x1=x1, x2=x2)",
    "  ex <- tx[seq_len(20), , drop=FALSE]",
    "  bw.fixed <- npindexbw(xdat=tx, ydat=y, bws=c(1, 1, 0.85), bandwidth.compute=FALSE, regtype='lc', bwtype='fixed')",
    "  bw.adp <- npindexbw(xdat=tx, ydat=y, bws=c(1, 1, 9), bandwidth.compute=FALSE, regtype='lc', bwtype='adaptive_nn')",
    "  fit.fixed <- npindex(bws=bw.fixed, txdat=tx, tydat=y, exdat=ex, gradients=FALSE)",
    "  fit.adp <- npindex(bws=bw.adp, txdat=tx, tydat=y, exdat=ex, gradients=FALSE)",
    "  stopifnot(max(abs(as.vector(fit.fixed$mean) - as.vector(fit.adp$mean))) > 1e-6)",
    "  cat('ATTACH_NPINDEX_ADAPTIVE_PUBLIC_OK\\n')",
    "}"
  ), script, useBytes = TRUE)

  env_common <- subprocess_env()
  skip_if(is.null(env_common), "local npRmpi install unavailable for subprocess smoke")
  res <- run_cmd_subprocess(
    mpiexec,
    args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
    timeout = 120L,
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
      timeout = 120L,
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
  expect_true(any(grepl("ATTACH_NPINDEX_ADAPTIVE_PUBLIC_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("attach npindex nearest-neighbor exact route stays green under mpiexec when enabled", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_TEST=1 to run attach-mode smoke")
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  script <- tempfile("npRmpi-attach-npindex-nn-exact-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "is.master <- isTRUE(npRmpi.init(mode='attach', quiet=TRUE, autodispatch=TRUE))",
    "on.exit({",
    "  try(npRmpi.quit(mode='attach'), silent=TRUE)",
    "  if (isTRUE(is.master)) try(Rmpi::mpi.quit(), silent=TRUE)",
    "}, add=TRUE)",
    "options(np.messages=FALSE)",
    "if (isTRUE(is.master)) {",
    "  set.seed(314163)",
    "  n <- 70L",
    "  x1 <- rnorm(n)",
    "  x2 <- rnorm(n)",
    "  y <- x1 - x2 + rnorm(n, sd=0.2)",
    "  tx <- data.frame(x1=x1, x2=x2)",
    "  bw.gen <- npindexbw(xdat=tx, ydat=y, regtype='lc', bwtype='generalized_nn', nmulti=1)",
    "  fit.gen <- npindex(bws=bw.gen, txdat=tx, tydat=y, gradients=FALSE)",
    "  hat.gen <- npindexhat(bws=bw.gen, txdat=tx, exdat=tx, y=y, output='apply', s=0L)",
    "  stopifnot(all(is.finite(fit.gen$mean)))",
    "  stopifnot(max(abs(as.vector(hat.gen) - as.vector(fit.gen$mean))) < 1e-8)",
    "  bw.adp <- npindexbw(xdat=tx, ydat=y, regtype='lc', bwtype='adaptive_nn', nmulti=1)",
    "  fit.adp <- npindex(bws=bw.adp, txdat=tx, tydat=y, gradients=FALSE)",
    "  hat.adp <- npindexhat(bws=bw.adp, txdat=tx, exdat=tx, y=y, output='apply', s=0L)",
    "  stopifnot(all(is.finite(fit.adp$mean)))",
    "  stopifnot(max(abs(as.vector(hat.adp) - as.vector(fit.adp$mean))) < 1e-8)",
    "  cat('ATTACH_NPINDEX_NN_EXACT_OK\\n')",
    "}"
  ), script, useBytes = TRUE)

  env_common <- subprocess_env()
  skip_if(is.null(env_common), "local npRmpi install unavailable for subprocess smoke")
  res <- run_cmd_subprocess(
    mpiexec,
    args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
    timeout = 120L,
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
      timeout = 120L,
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
  expect_true(any(grepl("ATTACH_NPINDEX_NN_EXACT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("attach smooth-coefficient ll coef plot-data route completes under mpiexec when enabled", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_TEST=1 to run attach-mode smoke")
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  script <- tempfile("npRmpi-attach-scoef-ll-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(mode='attach', quiet=TRUE)",
    "if (mpi.comm.rank(1L) == 0L) {",
    "  options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "  set.seed(105)",
    "  n <- 60",
    "  x <- runif(n)",
    "  z <- runif(n, -2, 2)",
    "  y <- x * exp(z) * (1 + rnorm(n, sd=0.15))",
    "  fit <- npscoef(y ~ x | z, regtype='ll', betas=TRUE)",
    "  pdf(file=tempfile(fileext='.pdf'))",
    "  on.exit(dev.off(), add=TRUE)",
    "  out <- suppressWarnings(plot(",
    "    fit,",
    "    coef=TRUE,",
    "    coef.index=1,",
    "    perspective=FALSE,",
    "    neval=20,",
    "    plot.behavior='plot-data',",
    "    plot.errors.method='none'))",
    "  stopifnot(is.list(out))",
    "  stopifnot(length(out) > 0L)",
    "  stopifnot(all(vapply(out, inherits, logical(1), 'smoothcoefficient')))",
    "  cat('ATTACH_SCOEF_LL_PLOTDATA_OK\\n')",
    "  npRmpi.quit(mode='attach')",
    "}"
  ), script, useBytes = TRUE)

  env_common <- subprocess_env()
  skip_if(is.null(env_common), "local npRmpi install unavailable for subprocess smoke")
  res <- run_cmd_subprocess(
    mpiexec,
    args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
    timeout = 120L,
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
      timeout = 120L,
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
  expect_true(any(grepl("ATTACH_SCOEF_LL_PLOTDATA_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("attach mode NP transition close stays stable when stress test is enabled", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_STRESS_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_STRESS_TEST=1 to run attach close stress")
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  script <- tempfile("npRmpi-attach-stress-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(mode='attach', quiet=TRUE)",
    "if (mpi.comm.rank(1L) == 0L) {",
    "  suppressPackageStartupMessages(library(MASS))",
    "  data(birthwt)",
    "  birthwt$low <- factor(birthwt$low)",
    "  birthwt$smoke <- factor(birthwt$smoke)",
    "  birthwt$race <- factor(birthwt$race)",
    "  birthwt$ht <- factor(birthwt$ht)",
    "  birthwt$ui <- factor(birthwt$ui)",
    "  birthwt$ftv <- ordered(birthwt$ftv)",
    "  bw <- npcdensbw(low~smoke+race+ht+ui+ftv+age+lwt, data=birthwt, nmulti=1)",
    "  fit <- npconmode(bws=bw)",
    "  stopifnot(inherits(fit, 'conmode'))",
    "  cat('ATTACH_NPCONMODE_ROUTE_OK\\n')",
    "  npRmpi.quit(mode='attach')",
    "}"
  ), script, useBytes = TRUE)

  env_common <- subprocess_env()
  skip_if(is.null(env_common), "local npRmpi install unavailable for subprocess smoke")

  run_once <- function(np) {
    res <- run_cmd_subprocess(
      mpiexec,
      args = c("-n", as.character(np), file.path(R.home("bin"), "Rscript"), "--no-save", script),
      timeout = 120L,
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
        args = c("-n", as.character(np), file.path(R.home("bin"), "Rscript"), "--no-save", script),
        timeout = 120L,
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
    res
  }

  res2 <- run_once(2L)
  if (res2$status != 0L && .is_mpi_init_env_failure(res2$output))
    skip("MPI runtime interface unavailable in this environment for attach stress")
  expect_equal(res2$status, 0L, info = paste(res2$output, collapse = "\n"))
  expect_true(any(grepl("ATTACH_NPCONMODE_ROUTE_OK", res2$output, fixed = TRUE)),
              info = paste(res2$output, collapse = "\n"))

  res3 <- run_once(3L)
  if (res3$status != 0L && .is_mpi_init_env_failure(res3$output))
    skip("MPI runtime interface unavailable in this environment for attach stress")
  expect_equal(res3$status, 0L, info = paste(res3$output, collapse = "\n"))
  expect_true(any(grepl("ATTACH_NPCONMODE_ROUTE_OK", res3$output, fixed = TRUE)),
              info = paste(res3$output, collapse = "\n"))
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
    "  suppressPackageStartupMessages(library(MASS))",
    "  mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)",
    "  mpi.bcast.cmd(options(np.messages=FALSE, npRmpi.autodispatch=FALSE), caller.execute=TRUE)",
    "  set.seed(42)",
    "  n <- 80",
    "  x <- runif(n)",
    "  z <- runif(n)",
    "  y <- sin(2*pi*x) + 0.5*z + rnorm(n, sd=0.1)",
    "  d <- data.frame(x=x, y=y)",
    "  d1 <- data.frame(x=x)",
    "  d2 <- data.frame(x=x, z=z)",
    "  dz <- data.frame(z=z)",
    "  d.cop <- data.frame(x=x, y=y)",
    "  u.cop <- data.frame(x=c(0.25,0.5,0.75), y=c(0.25,0.5,0.75))",
    "  data(birthwt)",
    "  bdat <- birthwt",
    "  bdat$low <- factor(bdat$low)",
    "  bdat$smoke <- factor(bdat$smoke)",
    "  bdat$race <- factor(bdat$race)",
    "  bdat$ht <- factor(bdat$ht)",
    "  bdat$ui <- factor(bdat$ui)",
    "  bdat$ftv <- ordered(bdat$ftv)",
    "  xd <- data.frame(u=x)",
    "  yd <- data.frame(v=y)",
    "  mpi.bcast.Robj2slave(d)",
    "  mpi.bcast.Robj2slave(d1)",
    "  mpi.bcast.Robj2slave(d2)",
    "  mpi.bcast.Robj2slave(dz)",
    "  mpi.bcast.Robj2slave(d.cop)",
    "  mpi.bcast.Robj2slave(u.cop)",
    "  mpi.bcast.Robj2slave(bdat)",
    "  mpi.bcast.Robj2slave(xd)",
    "  mpi.bcast.Robj2slave(yd)",
    "  mpi.bcast.Robj2slave(y)",
    "  mpi.bcast.cmd(bw <- npregbw(y~x, data=d, regtype='lc', bwmethod='cv.ls', nmulti=1), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit <- npreg(bws=bw, data=d, gradients=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(bw.sc <- npscoefbw(xdat=d1, ydat=y, zdat=dz, regtype='lc', nmulti=1), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit.sc <- npscoef(bws=bw.sc, gradients=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(bw.pl <- npplregbw(xdat=d1, ydat=y, zdat=dz, regtype='lc', nmulti=1), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit.pl <- npplreg(bws=bw.pl, gradients=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(bw.si <- npindexbw(xdat=d2, ydat=y, regtype='lc', nmulti=1), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit.si <- npindex(bws=bw.si, gradients=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(bw.ud <- npudensbw(dat=xd, bws=0.35, bandwidth.compute=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit.ud <- npudens(tdat=xd, bws=bw.ud), caller.execute=TRUE)",
    "  mpi.bcast.cmd(bw.uf <- npudistbw(dat=xd, bws=0.35, bandwidth.compute=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit.uf <- npudist(tdat=xd, bws=bw.uf), caller.execute=TRUE)",
    "  mpi.bcast.cmd(bw.cd <- npcdensbw(xdat=xd, ydat=yd, bws=c(0.45,0.45), bandwidth.compute=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit.cd <- npcdens(txdat=xd, tydat=yd, bws=bw.cd), caller.execute=TRUE)",
    "  mpi.bcast.cmd(bw.cf <- npcdistbw(xdat=xd, ydat=yd, bws=c(0.45,0.45), bandwidth.compute=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit.cf <- npcdist(txdat=xd, tydat=yd, bws=bw.cf), caller.execute=TRUE)",
    "  mpi.bcast.cmd(bw.cop <- npudistbw(~x+y, data=d.cop), caller.execute=TRUE)",
    "  mpi.bcast.cmd(cop <- npcopula(bws=bw.cop, data=d.cop, u=u.cop, n.quasi.inv=60), caller.execute=TRUE)",
    "  mpi.bcast.cmd(bw.cm <- npcdensbw(low~smoke+race+ht+ui+ftv+age+lwt, data=bdat, nmulti=1), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit.cm <- npconmode(bws=bw.cm), caller.execute=TRUE)",
    "  stopifnot(inherits(fit, 'npregression'))",
    "  stopifnot(inherits(fit.sc, 'smoothcoefficient'))",
    "  stopifnot(inherits(fit.pl, 'plregression'))",
    "  stopifnot(inherits(fit.si, 'singleindex'))",
    "  stopifnot(inherits(fit.ud, 'npdensity'))",
    "  stopifnot(inherits(fit.uf, 'npdistribution'))",
    "  stopifnot(inherits(fit.cd, 'condensity'))",
    "  stopifnot(inherits(fit.cf, 'condistribution'))",
    "  stopifnot(inherits(cop, 'data.frame'))",
    "  stopifnot(nrow(cop) == 9L)",
    "  stopifnot(all(is.finite(cop$copula)))",
    "  stopifnot(inherits(fit.cm, 'conmode'))",
    "  cat('PROFILE_NPCONMODE_ROUTE_OK\\n')",
    "  cat('PROFILE_NPCOPULA_ROUTE_OK\\n')",
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
  expect_true(any(grepl("PROFILE_NPCONMODE_ROUTE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("PROFILE_NPCOPULA_ROUTE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("profile npindexhat adaptive-nn exact owner route completes under mpiexec when enabled", {
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

  script <- tempfile("npRmpi-profile-npindexhat-exact-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "if (mpi.comm.rank(0L) == 0L) {",
    "  mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)",
    "  mpi.bcast.cmd(options(np.messages=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(set.seed(105), caller.execute=TRUE)",
    "  n <- 60L",
    "  x <- runif(n)",
    "  x2 <- runif(n)",
    "  y <- x * exp(x2) * (1 + rnorm(n, sd=0.15))",
    "  idxdat <- data.frame(y=y, x=x, x2=x2)",
    "  mpi.bcast.Robj2slave(idxdat)",
    "  mpi.bcast.cmd({",
    "    fit <- npindex(y ~ x + x2, data=idxdat, method='ichimura', bwtype='adaptive_nn')",
    "    hy <- npindexhat(bws=fit$bws, txdat=idxdat[c('x','x2')], exdat=idxdat[c('x','x2')], y=idxdat$y, output='apply', s=0L)",
    "    stopifnot(length(hy) == nrow(idxdat))",
    "    stopifnot(max(abs(as.vector(hy) - as.vector(fit$mean))) < 1e-8)",
    "  }, caller.execute=TRUE)",
    "  cat('PROFILE_NPINDEXHAT_ADAPTIVE_EXACT_OK\\n')",
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
    timeout = 120L,
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
      timeout = 120L,
      env = env_profile_fallback
    )
  }

  if (res$status != 0L && .is_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for profile-mode smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("PROFILE_NPINDEXHAT_ADAPTIVE_EXACT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("profile npindexhat adaptive-nn ll owner route stays exact under mpiexec when enabled", {
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

  script <- tempfile("npRmpi-profile-npindexhat-ll-owner-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "if (mpi.comm.rank(0L) == 0L) {",
    "  mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)",
    "  mpi.bcast.cmd(options(np.messages=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(set.seed(20260309), caller.execute=TRUE)",
    "  n <- 60L",
    "  x1 <- runif(n)",
    "  x2 <- runif(n)",
    "  y <- sin(x1 + x2) + rnorm(n, sd=0.05)",
    "  tx <- data.frame(x1=x1, x2=x2)",
    "  ex <- tx[seq_len(20L), , drop=FALSE]",
    "  mpi.bcast.Robj2slave(tx)",
    "  mpi.bcast.Robj2slave(ex)",
    "  mpi.bcast.Robj2slave(y)",
    "  mpi.bcast.cmd({",
    "    bw <- npindexbw(xdat=tx, ydat=y, regtype='ll', bwtype='adaptive_nn', nmulti=1L)",
    "    fit.mean <- npindex(bws=bw, txdat=tx, tydat=y, exdat=ex, gradients=FALSE)",
    "    fit.grad <- npindex(bws=bw, txdat=tx, tydat=y, exdat=ex, gradients=TRUE)",
    "    a.mean <- npindexhat(bws=bw, txdat=tx, exdat=ex, y=y, output='apply', s=0L)",
    "    H.mean <- npindexhat(bws=bw, txdat=tx, exdat=ex, output='matrix', s=0L)",
    "    a.grad <- npindexhat(bws=bw, txdat=tx, exdat=ex, y=y, output='apply', s=1L)",
    "    H.grad <- npindexhat(bws=bw, txdat=tx, exdat=ex, output='matrix', s=1L)",
    "    stopifnot(max(abs(as.vector(fit.mean$mean) - as.vector(a.mean))) < 1e-8)",
    "    stopifnot(max(abs(as.vector(a.mean) - as.vector(H.mean %*% y))) < 1e-8)",
    "    stopifnot(max(abs(as.vector(fit.grad$grad[,1L]) - as.vector(a.grad))) < 1e-8)",
    "    stopifnot(max(abs(as.vector(a.grad) - as.vector(H.grad %*% y))) < 1e-8)",
    "  }, caller.execute=TRUE)",
    "  cat('PROFILE_NPINDEXHAT_LL_OWNER_OK\\n')",
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
    timeout = 120L,
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
      timeout = 120L,
      env = env_profile_fallback
    )
  }

  if (res$status != 0L && .is_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for profile-mode smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("PROFILE_NPINDEXHAT_LL_OWNER_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("profile npindex adaptive-nn public route preserves bwtype semantics under mpiexec when enabled", {
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

  script <- tempfile("npRmpi-profile-npindex-public-adaptive-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "if (mpi.comm.rank(0L) == 0L) {",
    "  mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)",
    "  mpi.bcast.cmd(options(np.messages=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(set.seed(314161), caller.execute=TRUE)",
    "  n <- 70L",
    "  x1 <- runif(n)",
    "  x2 <- runif(n)",
    "  y <- sin(x1 + x2) + rnorm(n, sd=0.06)",
    "  tx <- data.frame(x1=x1, x2=x2)",
    "  ex <- tx[seq_len(20), , drop=FALSE]",
    "  mpi.bcast.Robj2slave(tx)",
    "  mpi.bcast.Robj2slave(ex)",
    "  mpi.bcast.Robj2slave(y)",
    "  mpi.bcast.cmd({",
    "    bw.fixed <- npindexbw(xdat=tx, ydat=y, bws=c(1, 1, 0.85), bandwidth.compute=FALSE, regtype='lc', bwtype='fixed')",
    "    bw.adp <- npindexbw(xdat=tx, ydat=y, bws=c(1, 1, 9), bandwidth.compute=FALSE, regtype='lc', bwtype='adaptive_nn')",
    "    fit.fixed <- npindex(bws=bw.fixed, txdat=tx, tydat=y, exdat=ex, gradients=FALSE)",
    "    fit.adp <- npindex(bws=bw.adp, txdat=tx, tydat=y, exdat=ex, gradients=FALSE)",
    "    stopifnot(max(abs(as.vector(fit.fixed$mean) - as.vector(fit.adp$mean))) > 1e-6)",
    "  }, caller.execute=TRUE)",
    "  cat('PROFILE_NPINDEX_ADAPTIVE_PUBLIC_OK\\n')",
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
    timeout = 120L,
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
      timeout = 120L,
      env = env_profile_fallback
    )
  }

  if (res$status != 0L && .is_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for profile-mode smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("PROFILE_NPINDEX_ADAPTIVE_PUBLIC_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("profile npindex nearest-neighbor exact route stays green under mpiexec when enabled", {
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

  script <- tempfile("npRmpi-profile-npindex-nn-exact-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "if (mpi.comm.rank(0L) == 0L) {",
    "  mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)",
    "  mpi.bcast.cmd(options(np.messages=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(set.seed(314163), caller.execute=TRUE)",
    "  n <- 70L",
    "  x1 <- rnorm(n)",
    "  x2 <- rnorm(n)",
    "  y <- x1 - x2 + rnorm(n, sd=0.2)",
    "  tx <- data.frame(x1=x1, x2=x2)",
    "  mpi.bcast.Robj2slave(tx)",
    "  mpi.bcast.Robj2slave(y)",
    "  mpi.bcast.cmd({",
    "    bw.gen <- npindexbw(xdat=tx, ydat=y, regtype='lc', bwtype='generalized_nn', nmulti=1)",
    "    fit.gen <- npindex(bws=bw.gen, txdat=tx, tydat=y, gradients=FALSE)",
    "    hat.gen <- npindexhat(bws=bw.gen, txdat=tx, exdat=tx, y=y, output='apply', s=0L)",
    "    stopifnot(all(is.finite(fit.gen$mean)))",
    "    stopifnot(max(abs(as.vector(hat.gen) - as.vector(fit.gen$mean))) < 1e-8)",
    "    bw.adp <- npindexbw(xdat=tx, ydat=y, regtype='lc', bwtype='adaptive_nn', nmulti=1)",
    "    fit.adp <- npindex(bws=bw.adp, txdat=tx, tydat=y, gradients=FALSE)",
    "    hat.adp <- npindexhat(bws=bw.adp, txdat=tx, exdat=tx, y=y, output='apply', s=0L)",
    "    stopifnot(all(is.finite(fit.adp$mean)))",
    "    stopifnot(max(abs(as.vector(hat.adp) - as.vector(fit.adp$mean))) < 1e-8)",
    "  }, caller.execute=TRUE)",
    "  cat('PROFILE_NPINDEX_NN_EXACT_OK\\n')",
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
    timeout = 120L,
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
      timeout = 120L,
      env = env_profile_fallback
    )
  }

  if (res$status != 0L && .is_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for profile-mode smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("PROFILE_NPINDEX_NN_EXACT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session npindex nearest-neighbor exact route selects integer support and stays green", {
  skip_on_cran()
  env <- subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  res <- run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
      "set.seed(314163)",
      "n <- 70L",
      "x1 <- rnorm(n)",
      "x2 <- rnorm(n)",
      "y <- x1 - x2 + rnorm(n, sd=0.2)",
      "tx <- data.frame(x1=x1, x2=x2)",
      "bw.gen <- npindexbw(xdat=tx, ydat=y, regtype='lc', bwtype='generalized_nn', nmulti=1)",
      "stopifnot(bw.gen$bw >= 1)",
      "stopifnot(abs(bw.gen$bw - as.integer(bw.gen$bw)) < 1e-12)",
      "fit.gen <- npindex(bws=bw.gen, txdat=tx, tydat=y, gradients=FALSE)",
      "hat.gen <- npindexhat(bws=bw.gen, txdat=tx, exdat=tx, y=y, output='apply', s=0L)",
      "stopifnot(all(is.finite(fit.gen$mean)))",
      "stopifnot(max(abs(as.vector(hat.gen) - as.vector(fit.gen$mean))) < 1e-8)",
      "bw.adp <- npindexbw(xdat=tx, ydat=y, regtype='lc', bwtype='adaptive_nn', nmulti=1)",
      "stopifnot(bw.adp$bw >= 1)",
      "stopifnot(abs(bw.adp$bw - as.integer(bw.adp$bw)) < 1e-12)",
      "fit.adp <- npindex(bws=bw.adp, txdat=tx, tydat=y, gradients=FALSE)",
      "hat.adp <- npindexhat(bws=bw.adp, txdat=tx, exdat=tx, y=y, output='apply', s=0L)",
      "stopifnot(all(is.finite(fit.adp$mean)))",
      "stopifnot(max(abs(as.vector(hat.adp) - as.vector(fit.adp$mean))) < 1e-8)",
      "cat('SESSION_NPINDEX_NN_EXACT_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SESSION_NPINDEX_NN_EXACT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
