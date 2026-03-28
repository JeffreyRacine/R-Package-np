# Helper for MPI tests
.mpi_check_context <- function() {
  identical(Sys.getenv("_R_CHECK_PACKAGE_NAME_", ""), "npRmpi")
}

.mpi_pool_active <- function() {
  if (.mpi_check_context())
    return(FALSE)
  if (!isTRUE(getOption("npRmpi.mpi.initialized", FALSE)))
    return(FALSE)
  isTRUE(try(mpi.comm.size(1) > 1, silent = TRUE))
}

spawn_mpi_slaves <- function(n = 1L) {
  # R CMD check environments are not a stable MPI runtime target.
  if (.mpi_check_context()) {
    return(FALSE)
  }

  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)

  # Reuse an active pool instead of re-initializing nested MPI sessions.
  if (.mpi_pool_active()) {
    return(TRUE)
  }

  n <- as.integer(n)[1L]
  if (!is.finite(n) || n < 1L) {
    return(FALSE)
  }

  ok <- try({
    npRmpi.init(nslaves = n, quiet = TRUE)
    TRUE
  }, silent = TRUE)

  isTRUE(ok) && .mpi_pool_active()
}

close_mpi_slaves <- function(force = FALSE) {
  if (.mpi_check_context())
    return(invisible())
  if (!.mpi_pool_active())
    return(invisible())

  try(npRmpi.quit(force = force), silent = TRUE)
}

npRmpi_run_rscript_subprocess <- function(lines, timeout = 45L, env = character()) {
  script <- tempfile("npRmpi-subprocess-", fileext = ".R")
  scoped.lines <- c(
    "local({",
    "  on.exit({",
    "    if (exists('npRmpi.quit', mode = 'function')) try(npRmpi.quit(force = TRUE), silent = TRUE)",
    "    if (exists('mpi.finalize', mode = 'function')) try(mpi.finalize(), silent = TRUE)",
    "  }, add = TRUE)",
    paste0("  ", lines),
    "})"
  )
  writeLines(scoped.lines, script, useBytes = TRUE)
  on.exit(unlink(script), add = TRUE)

  cmd <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(
    cmd,
    c("--no-save", script),
    stdout = TRUE,
    stderr = TRUE,
    timeout = timeout,
    env = env
  ))
  status <- attr(out, "status")
  if (is.null(status))
    status <- 0L

  list(status = as.integer(status), output = out)
}

npRmpi_subprocess_env <- local({
  lib.path.cache <- NULL

  function(extra = character()) {
    if (is.null(lib.path.cache) || !dir.exists(lib.path.cache)) {
      installed.pkg <- tryCatch(find.package("npRmpi"), error = function(e) "")
      if (nzchar(installed.pkg) &&
          file.exists(file.path(installed.pkg, "Meta", "package.rds"))) {
        lib.path.cache <<- dirname(installed.pkg)
      }
    }

    if (is.null(lib.path.cache) || !dir.exists(lib.path.cache)) {
      candidates <- unique(Filter(
        nzchar,
        c(
          tryCatch(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), error = function(e) ""),
          tryCatch(normalizePath(getwd(), mustWork = TRUE), error = function(e) ""),
          tryCatch(normalizePath(file.path(getwd(), ".."), mustWork = TRUE), error = function(e) ""),
          tryCatch(normalizePath(file.path(getwd(), "..", ".."), mustWork = TRUE), error = function(e) "")
        )
      ))
      pkg.root <- ""
      for (candidate in candidates) {
        if (file.exists(file.path(candidate, "DESCRIPTION")) &&
            identical(basename(candidate), "np-npRmpi")) {
          pkg.root <- candidate
          break
        }
      }
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
        unlink(lib.path.cache, recursive = TRUE, force = TRUE)
        lib.path.cache <<- NULL
        return(NULL)
      }
    }

    c(
      sprintf("R_LIBS=%s", paste(c(lib.path.cache, .libPaths()), collapse = .Platform$path.sep)),
      "NP_RMPI_NO_REUSE_SLAVES=1",
      extra
    )
  }
})
