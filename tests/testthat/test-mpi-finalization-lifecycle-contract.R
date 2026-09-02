rg0_mpi_status <- function(output) {
  status <- attr(output, "status")
  if (is.null(status)) 0L else as.integer(status)
}

rg0_run_profile <- function(lines, n = 2L, timeout = 45L, iface = "en0") {
  mpiexec <- Sys.which("mpiexec")
  if (!nzchar(mpiexec))
    return(NULL)

  installed <- tryCatch(find.package("npRmpi"), error = function(e) "")
  if (!nzchar(installed))
    return(NULL)
  profile <- file.path(installed, "Rprofile")
  if (!file.exists(profile))
    return(NULL)

  script <- tempfile("npRmpi-finalize-profile-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(lines, script, useBytes = TRUE)

  env <- npRmpi_subprocess_env(c(
    sprintf("R_PROFILE_USER=%s", profile),
    "R_PROFILE=",
    "NP_RMPI_PROFILE_RECV_TIMEOUT_SEC=30",
    sprintf("FI_TCP_IFACE=%s", iface),
    "FI_PROVIDER=tcp",
    sprintf("FI_SOCKETS_IFACE=%s", iface)
  ))
  if (is.null(env))
    return(NULL)

  output <- suppressWarnings(system2(
    mpiexec,
    c("-n", as.character(as.integer(n)), file.path(R.home("bin"), "Rscript"),
      "--no-save", script),
    stdout = TRUE,
    stderr = TRUE,
    timeout = timeout,
    env = env
  ))
  list(status = rg0_mpi_status(output), output = output)
}

test_that("native MPI lifecycle state and finalization are idempotent", {
  env <- npRmpi_subprocess_env(c("R_PROFILE_USER=", "R_PROFILE="))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess proof")

  result <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "ns <- asNamespace('npRmpi')",
      "state <- get('.npRmpi_mpi_runtime_state', ns, inherits=FALSE)",
      "finish <- get('.npRmpi_finalize_before_unload', ns, inherits=FALSE)",
      "before <- state()",
      "stopifnot(isTRUE(before[['initialized']]), !isTRUE(before[['finalized']]))",
      "stopifnot(identical(as.integer(.Call('mpi_finalize', PACKAGE='npRmpi')), 1L))",
      "stopifnot(identical(as.integer(.Call('mpi_finalize', PACKAGE='npRmpi')), 1L))",
      "after <- state()",
      "stopifnot(isTRUE(after[['initialized']]), isTRUE(after[['finalized']]))",
      "stopifnot(isTRUE(finish('test already-finalized entry')))",
      "stopifnot(!isTRUE(getOption('npRmpi.mpi.initialized', TRUE)))",
      "cat('MPI_FINALIZE_IDEMPOTENT_OK\\n')"
    ),
    timeout = 30L,
    env = env,
    cleanup = FALSE
  )

  expect_equal(result$status, 0L, info = paste(result$output, collapse = "\n"))
  expect_true(any(grepl("MPI_FINALIZE_IDEMPOTENT_OK", result$output, fixed = TRUE)),
              info = paste(result$output, collapse = "\n"))
})

test_that("profile exit owns normal and uncaught-error teardown", {
  skip_on_cran()

  clean <- rg0_run_profile(c(
    "cat('MPI_PROFILE_CLEAN_EXIT_OK\\n')"
  ))
  skip_if(is.null(clean), "profile-mode MPI subprocess unavailable")
  if (clean$status != 0L &&
      !any(grepl("MPI_PROFILE_CLEAN_EXIT_OK", clean$output, fixed = TRUE)))
    clean <- rg0_run_profile("cat('MPI_PROFILE_CLEAN_EXIT_OK\\n')", iface = "lo0")

  expect_equal(clean$status, 0L, info = paste(clean$output, collapse = "\n"))
  expect_true(any(grepl("MPI_PROFILE_CLEAN_EXIT_OK", clean$output, fixed = TRUE)),
              info = paste(clean$output, collapse = "\n"))
  expect_false(any(grepl("segfault|BAD TERMINATION|signal 11", clean$output,
                         ignore.case = TRUE)))

  failed <- rg0_run_profile(c(
    "cat('MPI_PROFILE_ERROR_EXIT_WITNESS\\n')",
    "stop('MPI_PROFILE_INTENTIONAL_ERROR', call.=FALSE)"
  ))
  skip_if(is.null(failed), "profile-mode MPI subprocess unavailable")
  if (!any(grepl("MPI_PROFILE_ERROR_EXIT_WITNESS", failed$output, fixed = TRUE)))
    failed <- rg0_run_profile(c(
      "cat('MPI_PROFILE_ERROR_EXIT_WITNESS\\n')",
      "stop('MPI_PROFILE_INTENTIONAL_ERROR', call.=FALSE)"
    ), iface = "lo0")

  expect_equal(failed$status, 1L, info = paste(failed$output, collapse = "\n"))
  expect_true(any(grepl("MPI_PROFILE_ERROR_EXIT_WITNESS", failed$output, fixed = TRUE)),
              info = paste(failed$output, collapse = "\n"))
  expect_true(any(grepl("MPI_PROFILE_INTENTIONAL_ERROR", failed$output, fixed = TRUE)),
              info = paste(failed$output, collapse = "\n"))
  expect_false(any(grepl("segfault|BAD TERMINATION|signal 11", failed$output,
                         ignore.case = TRUE)))
})

test_that("installed profile delegates to the package-owned exit path", {
  profile <- system.file("Rprofile", package = "npRmpi")
  expect_true(nzchar(profile))
  text <- readLines(profile, warn = FALSE)
  expect_true(any(grepl(".npRmpi_finalize_before_unload", text, fixed = TRUE)))
  expect_false(any(grepl(".Call(\"mpi_finalize\"", text, fixed = TRUE)))
  expect_false(any(grepl("options(error", text, fixed = TRUE)))
})
