local_npRmpi_subprocess_env <- function(extra = character()) {
  pkg.root <- tryCatch(
    normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
    error = function(e) ""
  )
  if (!nzchar(pkg.root))
    return(NULL)

  src.copy <- tempfile("npRmpi-fitted-subprocess-src-")
  dir.create(src.copy, recursive = TRUE, showWarnings = FALSE)
  pkg.copy <- file.path(src.copy, "np-npRmpi")
  if (!isTRUE(file.copy(pkg.root, src.copy, recursive = TRUE)))
    return(NULL)

  vignette.index <- file.path(pkg.copy, "build", "vignette.rds")
  if (file.exists(vignette.index))
    unlink(vignette.index, force = TRUE)

  build.dir <- tempfile("npRmpi-fitted-subprocess-build-")
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

  lib.path <- tempfile("npRmpi-fitted-subprocess-lib-")
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

test_that("npcdens fitted slice repair works in subprocess session mode", {
  skip_on_cran()
  env <- local_npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
      "set.seed(20260322)",
      "x <- runif(36, -1, 1)",
      "y <- sin(2*pi*x) + rnorm(36, sd=0.18)",
      "bw <- npcdensbw(xdat=data.frame(x=x), ydat=data.frame(y=y), bws=c(0.29, 0.23), bandwidth.compute=FALSE, regtype='lp', degree=3L)",
      "fit.raw <- npcdens(bws=bw, txdat=data.frame(x=x), tydat=data.frame(y=y))",
      "fit.slice <- npcdens(bws=bw, txdat=data.frame(x=x), tydat=data.frame(y=y), proper=TRUE, proper.control=list(mode='slice', apply='fitted', slice.grid.size=31L, slice.extend.factor=0))",
      "stopifnot(isTRUE(fit.slice$proper.applied))",
      "stopifnot(identical(fit.slice$proper.info$route, 'slice'))",
      "stopifnot(identical(fit.slice$proper.info$apply.scope, 'fitted'))",
      "stopifnot(isTRUE(all.equal(fit.slice$condens.raw, fit.raw$condens, tolerance=1e-12)))",
      "cat('NPCDENS_FITTED_SLICE_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPCDENS_FITTED_SLICE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("npcdens predict newdata stays unchanged under apply='fitted' in subprocess session mode", {
  skip_on_cran()
  env <- local_npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
      "set.seed(20260322)",
      "x <- runif(70, -1, 1)",
      "y <- sin(2*pi*x) + rnorm(70, sd=0.2)",
      "nd <- rbind(data.frame(y=c(-0.7, -0.15, 0.45), x=rep(-0.35, 3L)), data.frame(y=c(-0.25, 0.2, 0.75, 1.0), x=rep(0.4, 4L)))",
      "ctrl <- list(mode='slice', apply='fitted', slice.grid.size=21L, slice.extend.factor=0)",
      "bw <- npcdensbw(xdat=data.frame(x=x), ydat=data.frame(y=y), bws=c(0.27, 0.21), bandwidth.compute=FALSE, regtype='lp', degree=3L)",
      "fit <- npcdens(bws=bw, txdat=data.frame(x=x), tydat=data.frame(y=y), proper=TRUE)",
      "base <- predict(fit, newdata=nd, proper=FALSE)",
      "req <- predict(fit, newdata=nd, proper.control=ctrl)",
      "stopifnot(isTRUE(all.equal(as.numeric(base), as.numeric(req), tolerance=1e-12)))",
      "cat('NPCDENS_FITTED_PREDICT_OK\\n')"
    ),
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPCDENS_FITTED_PREDICT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
