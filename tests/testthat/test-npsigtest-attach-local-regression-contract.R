run_attach_npsigtest_cmd <- function(cmd, args = character(), timeout = 60L, env = character()) {
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

ensure_attach_npsigtest_lib <- local({
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

    lib.path.cache <<- tempfile("npRmpi-attach-npsigtest-lib-")
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

attach_npsigtest_env <- function(extra = character()) {
  lib.path <- ensure_attach_npsigtest_lib()
  if (is.null(lib.path))
    return(NULL)

  c(
    sprintf("R_LIBS=%s", paste(c(lib.path, .libPaths()), collapse = .Platform$path.sep)),
    extra
  )
}

.is_attach_npsigtest_mpi_init_env_failure <- function(output) {
  any(grepl("OFI call ep_enable failed", output, fixed = TRUE)) ||
    any(grepl("Fatal error in internal_Init", output, fixed = TRUE)) ||
    any(grepl("MPI_Init", output, fixed = TRUE) & grepl("failed", output, ignore.case = TRUE))
}

test_that("attach npsigtest npregression route remains valid under local regression mode", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_TEST=1 to run attach-mode smoke")

  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  env_common <- attach_npsigtest_env()
  skip_if(is.null(env_common), "local npRmpi install unavailable for attach regression")

  script <- tempfile("npRmpi-attach-npsigtest-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)

  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(mode='attach', quiet=TRUE)",
    "options(np.messages=FALSE)",
    "if (mpi.comm.rank(1L) == 0L) {",
    "  set.seed(42)",
    "  n <- 100L",
    "  z <- factor(rbinom(n, 1, 0.5))",
    "  x1 <- rnorm(n)",
    "  x2 <- runif(n, -2, 2)",
    "  y <- x1 + x2 + rnorm(n)",
    "  mydat <- data.frame(z, x1, x2, y)",
    "  rm(z, x1, x2, y)",
    "  model <- npreg(y ~ z + x1 + x2, regtype='ll', bwmethod='cv.aic', data=mydat)",
    "  out <- npsigtest(model, boot.num=9)",
    "  stopifnot(inherits(out, 'sigtest'))",
    "  stopifnot(length(out$P) >= 1L)",
    "  cat('ATTACH_NPSIGTEST_OK\\n')",
    "  npRmpi.quit(mode='attach')",
    "}"
  ), script, useBytes = TRUE)

  run_once <- function(iface) {
    run_attach_npsigtest_cmd(
      mpiexec,
      args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
      timeout = 120L,
      env = c(
        env_common,
        "R_PROFILE_USER=",
        "R_PROFILE=",
        sprintf("FI_TCP_IFACE=%s", iface),
        "FI_PROVIDER=tcp",
        sprintf("FI_SOCKETS_IFACE=%s", iface)
      )
    )
  }

  res <- run_once("en0")
  if (res$status != 0L)
    res <- run_once("lo0")

  if (res$status != 0L && .is_attach_npsigtest_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for attach regression")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("ATTACH_NPSIGTEST_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
