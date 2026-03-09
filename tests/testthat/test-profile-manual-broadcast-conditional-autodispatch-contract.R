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

test_that("profile manual-broadcast user-style conditional fixed-bandwidth routes stay stable with autodispatch enabled", {
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

  script <- tempfile("npRmpi-profile-user-conditional-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "if (mpi.comm.rank(0L) == 0L) {",
    "  mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)",
    "  mpi.bcast.cmd(options(np.messages=FALSE), caller.execute=TRUE)",
    "  set.seed(11)",
    "  xd <- data.frame(x = seq(0.1, 0.9, length.out = 18))",
    "  yd <- data.frame(y = sin(xd$x))",
    "  mpi.bcast.Robj2slave(xd)",
    "  mpi.bcast.Robj2slave(yd)",
    "  mpi.bcast.cmd(bw.cd <- npcdensbw(xdat=xd, ydat=yd, bws=c(0.25,0.25), bandwidth.compute=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit.cd <- npcdens(txdat=xd, tydat=yd, bws=bw.cd), caller.execute=TRUE)",
    "  mpi.bcast.cmd(bw.qr <- npcdistbw(xdat=xd, ydat=yd, bws=c(0.25,0.25), bandwidth.compute=FALSE), caller.execute=TRUE)",
    "  mpi.bcast.cmd(fit.qr <- npqreg(bws=bw.qr, txdat=xd, tydat=yd, tau=0.4), caller.execute=TRUE)",
    "  stopifnot(inherits(fit.cd, 'condensity'))",
    "  stopifnot(inherits(fit.qr, 'qregression'))",
    "  cat('PROFILE_USER_CONDITIONAL_AUTODISPATCH_OK\\n')",
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
  expect_true(any(grepl("PROFILE_USER_CONDITIONAL_AUTODISPATCH_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
  expect_false(any(grepl("skipping nested auto-dispatch", res$output, fixed = TRUE)),
               info = paste(res$output, collapse = "\n"))
})
