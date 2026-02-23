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

test_that("session routing smoke completes in subprocess", {
  skip_on_cran()
  pkg_path <- tryCatch(find.package("npRmpi"), error = function(e) "")
  skip_if(!nzchar(pkg_path), "installed npRmpi unavailable for subprocess smoke")

  env <- sprintf("R_LIBS=%s", paste(.libPaths(), collapse = .Platform$path.sep))
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

test_that("skip-init mode fails fast without MPI pool crash", {
  skip_on_cran()
  pkg_path <- tryCatch(find.package("npRmpi"), error = function(e) "")
  skip_if(!nzchar(pkg_path), "installed npRmpi unavailable for subprocess smoke")

  env <- c(
    sprintf("R_LIBS=%s", paste(.libPaths(), collapse = .Platform$path.sep)),
    "NP_RMPI_SKIP_INIT=1"
  )
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
  pkg_path <- tryCatch(find.package("npRmpi"), error = function(e) "")
  skip_if(!nzchar(pkg_path), "installed npRmpi unavailable for subprocess smoke")

  env <- sprintf("R_LIBS=%s", paste(.libPaths(), collapse = .Platform$path.sep))
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
