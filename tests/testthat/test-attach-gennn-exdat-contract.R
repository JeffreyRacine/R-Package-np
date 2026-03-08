run_attach_cmd_subprocess <- function(cmd, args = character(), timeout = 60L, env = character()) {
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

ensure_attach_test_npRmpi_lib <- local({
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

    lib.path.cache <<- tempfile("npRmpi-attach-gennn-lib-")
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

attach_test_env <- function(extra = character()) {
  lib.path <- ensure_attach_test_npRmpi_lib()
  if (is.null(lib.path))
    return(NULL)

  c(
    sprintf("R_LIBS=%s", paste(c(lib.path, .libPaths()), collapse = .Platform$path.sep)),
    extra
  )
}

.is_attach_mpi_init_env_failure <- function(output) {
  any(grepl("OFI call ep_enable failed", output, fixed = TRUE)) ||
    any(grepl("Fatal error in internal_Init", output, fixed = TRUE)) ||
    any(grepl("MPI_Init", output, fixed = TRUE) & grepl("failed", output, ignore.case = TRUE))
}

test_that("attach generalized_nn exdat returns full evaluation grid", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_TEST=1 to run attach-mode smoke")

  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  env_common <- attach_test_env()
  skip_if(is.null(env_common), "local npRmpi install unavailable for attach regression")

  script <- tempfile("npRmpi-attach-gennn-exdat-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)

  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(mode='attach', quiet=TRUE)",
    "if (mpi.comm.rank(1L) == 0L) {",
    "  set.seed(123)",
    "  n <- 100",
    "  x <- runif(n)",
    "  y <- x^2 + rnorm(n, sd=0.1)",
    "  xe <- data.frame(x=seq(0, 1, length.out=50))",
    "  fit <- npreg(y~x, bwtype='generalized_nn', exdat=xe)",
    "  stopifnot(length(fit$mean) == 50L)",
    "  stopifnot(sum(fit$mean[26:50] == 0) == 0L)",
    "  stopifnot(all(is.finite(fit$mean)))",
    "  cat('ATTACH_GENNN_EXDAT_OK\\n')",
    "  npRmpi.quit(mode='attach')",
    "}"
  ), script, useBytes = TRUE)

  run_once <- function(iface) {
    run_attach_cmd_subprocess(
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

  if (res$status != 0L && .is_attach_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for attach regression")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("ATTACH_GENNN_EXDAT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
