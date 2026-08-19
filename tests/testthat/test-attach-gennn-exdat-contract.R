run_attach_cmd_subprocess <- function(cmd, args = character(), timeout = 20L,
                                      env = character()) {
  out <- suppressWarnings(system2(
    cmd, args, stdout = TRUE, stderr = TRUE, timeout = timeout, env = env
  ))
  status <- attr(out, "status")
  if (is.null(status))
    status <- 0L
  list(status = as.integer(status), output = out)
}

.is_attach_mpi_init_env_failure <- function(output) {
  any(grepl("OFI call ep_enable failed", output, fixed = TRUE)) ||
    any(grepl("Fatal error in internal_Init", output, fixed = TRUE)) ||
    any(grepl("MPI_Init", output, fixed = TRUE) &
          grepl("failed", output, ignore.case = TRUE))
}

test_that("attach generalized_nn exdat returns full evaluation grid", {
  skip_on_cran()
  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  env_common <- npRmpi_subprocess_env(c(
    "_R_CHECK_PACKAGE_NAME_=",
    "NP_RMPI_TEST_SUITE_POOL="
  ))
  skip_if(is.null(env_common),
          "installed npRmpi unavailable for attach regression")

  script <- tempfile("npRmpi-attach-gennn-exdat-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(mode='attach', quiet=TRUE)",
    "if (mpi.comm.rank(1L) == 0L) {",
    "  set.seed(123)",
    "  n <- 36L",
    "  x <- runif(n)",
    "  y <- x^2 + rnorm(n, sd=0.1)",
    "  d <- data.frame(x=x, y=y)",
    "  xe <- data.frame(x=seq(0, 1, length.out=12L))",
    "  bw <- npregbw(y~x, data=d, bws=7, bandwidth.compute=FALSE, bwtype='generalized_nn')",
    "  fit <- npreg(bws=bw, txdat=d['x'], tydat=d$y, exdat=xe)",
    "  stopifnot(length(fit$mean) == 12L)",
    "  stopifnot(all(is.finite(fit$mean)))",
    "  cat('ATTACH_GENNN_EXDAT_OK\\n')",
    "  npRmpi.quit(mode='attach')",
    "}"
  ), script, useBytes = TRUE)

  run_once <- function(iface) {
    run_attach_cmd_subprocess(
      mpiexec,
      args = c(
        "-n", "2", file.path(R.home("bin"), "Rscript"),
        "--no-save", script
      ),
      timeout = 15L,
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
  if (res$status != 0L && .is_attach_mpi_init_env_failure(res$output))
    res <- run_once("lo0")
  if (res$status != 0L && .is_attach_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for attach regression")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(
    any(grepl("ATTACH_GENNN_EXDAT_OK", res$output, fixed = TRUE)),
    info = paste(res$output, collapse = "\n")
  )
})
