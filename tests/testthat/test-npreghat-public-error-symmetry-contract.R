run_npreghat_error_symmetry_subprocess <- function(cmd,
                                                   args = character(),
                                                   timeout = 90L,
                                                   env = character()) {
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

ensure_npreghat_error_symmetry_npRmpi_lib <- local({
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

    lib.path.cache <<- tempfile("npRmpi-npreghat-error-symmetry-lib-")
    dir.create(lib.path.cache, recursive = TRUE, showWarnings = FALSE)

    out <- suppressWarnings(system2(
      file.path(R.home("bin"), "R"),
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

npreghat_error_symmetry_attach_env <- function(extra = character()) {
  lib.path <- ensure_npreghat_error_symmetry_npRmpi_lib()
  if (is.null(lib.path))
    return(NULL)

  c(
    sprintf("R_LIBS=%s", paste(c(lib.path, .libPaths()), collapse = .Platform$path.sep)),
    extra
  )
}

is_npreghat_error_symmetry_mpi_init_failure <- function(output) {
  any(grepl("OFI call ep_enable failed", output, fixed = TRUE)) ||
    any(grepl("Fatal error in internal_Init", output, fixed = TRUE)) ||
    any(grepl("MPI_Init", output, fixed = TRUE) & grepl("failed", output, ignore.case = TRUE))
}

test_that("public npreghat invalid derivative errors before active-pool broadcast", {
  skip_on_cran()
  skip_on_cran()

  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  env.common <- npreghat_error_symmetry_attach_env()
  skip_if(is.null(env.common), "local npRmpi install unavailable for attach regression")

  script <- tempfile("npRmpi-npreghat-error-symmetry-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)

  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch = FALSE, np.messages = FALSE)",
    "is.master <- isTRUE(npRmpi.init(mode = 'attach', quiet = TRUE, autodispatch = FALSE, np.messages = FALSE))",
    "if (!is.master) q('no', status = 0L)",
    "on.exit(try(npRmpi.quit(force = TRUE, mode = 'attach'), silent = TRUE), add = TRUE)",
    "stopifnot(max(as.integer(mpi.comm.size(1L)) - 1L, 0L) >= 1L)",
    "set.seed(20260630)",
    "n <- 50L",
    "x1 <- runif(n)",
    "x2 <- runif(n)",
    "y <- sin(2 * pi * x1) + x2^2 + rnorm(n, sd = 0.02)",
    "tx <- data.frame(x1 = x1, x2 = x2)",
    "ex <- data.frame(x1 = seq(0.1, 0.9, length.out = 10L), x2 = seq(0.15, 0.85, length.out = 10L))",
    "bw <- npregbw(xdat = tx, ydat = y, regtype = 'lp', degree = c(0L, 2L), basis = 'glp', bws = c(0.3, 0.3), bandwidth.compute = FALSE)",
    "err <- try(npreghat(bws = bw, txdat = tx, exdat = ex, s = c(1L, 0L)), silent = TRUE)",
    "stopifnot(inherits(err, 'try-error'))",
    "stopifnot(grepl(\"exceeds local polynomial degree\", as.character(err)))",
    "H <- npreghat(bws = bw, txdat = tx, exdat = ex, s = c(0L, 1L))",
    "stopifnot(is.matrix(H), nrow(H) == nrow(ex), all(is.finite(H)))",
    "out <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply', s = c(0L, 1L))",
    "stopifnot(length(out) == nrow(ex), all(is.finite(out)))",
    "cat('NPREGHAT_PUBLIC_ERROR_SYMMETRY_OK\\n')"
  ), script, useBytes = TRUE)

  run_once <- function(iface) {
    run_npreghat_error_symmetry_subprocess(
      mpiexec,
      args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
      timeout = 120L,
      env = c(
        env.common,
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

  if (res$status != 0L && is_npreghat_error_symmetry_mpi_init_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for attach regression")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPREGHAT_PUBLIC_ERROR_SYMMETRY_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
