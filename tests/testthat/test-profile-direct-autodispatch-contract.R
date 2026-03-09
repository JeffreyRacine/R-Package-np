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

test_that("profile direct autodispatch regression routes initialize embedded MPI state automatically", {
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

  script <- tempfile("npRmpi-profile-direct-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "if (mpi.comm.rank(0L) == 0L) {",
    "  set.seed(123)",
    "  n <- 60L",
    "  x <- data.frame(x = runif(n))",
    "  y <- sin(2*pi*x$x) + rnorm(n, sd = 0.1)",
    "  d <- data.frame(x = x$x, y = y)",
    "  bw.cv <- npregbw(xdat = x, ydat = y, regtype = 'lc', bwmethod = 'cv.ls', nmulti = 1)",
    "  fit.cv <- npreg(bws = bw.cv, txdat = x, tydat = y, gradients = FALSE)",
    "  bw.manual <- npregbw(xdat = x, ydat = y, regtype = 'lc', bws = 0.25, bandwidth.compute = FALSE)",
    "  fit.manual <- npreg(bws = bw.manual, txdat = x, tydat = y, gradients = FALSE)",
    "  fit.formula <- npreg(y ~ x, data = d, regtype = 'lc', bwmethod = 'cv.ls', nmulti = 1)",
    "  stopifnot(inherits(bw.cv, 'rbandwidth'))",
    "  stopifnot(inherits(fit.cv, 'npregression'))",
    "  stopifnot(inherits(fit.manual, 'npregression'))",
    "  stopifnot(inherits(fit.formula, 'npregression'))",
    "  cat('PROFILE_DIRECT_AUTODISPATCH_OK\\n')",
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
  expect_true(any(grepl("PROFILE_DIRECT_AUTODISPATCH_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
