test_that("spawn sessions finalize cleanly at process exit", {
  env <- npRmpi_subprocess_env(c(
    "_R_CHECK_PACKAGE_NAME_=",
    "NP_RMPI_TEST_SUITE_POOL="
  ))
  skip_if(is.null(env), "installed npRmpi unavailable for process-exit contract")

  for (close.mode in c("noquit", "soft", "force")) {
    close.lines <- switch(
      close.mode,
      noquit = character(),
      soft = "npRmpi.quit()",
      force = "npRmpi.quit(force = TRUE)"
    )
    marker <- paste0("MPI_PROCESS_EXIT_OK_", toupper(close.mode))
    result <- npRmpi_run_rscript_subprocess(
      lines = c(
        "suppressPackageStartupMessages(library(npRmpi))",
        "options(np.messages = FALSE)",
        "npRmpi.init(nslaves = 1L, quiet = TRUE)",
        "value <- npRmpi:::mpi.remote.exec(1L + 1L, simplify = TRUE)",
        "stopifnot(all(unlist(value, use.names = FALSE) == 2L))",
        close.lines,
        sprintf("cat('%s\\n')", marker)
      ),
      timeout = 45L,
      env = env,
      cleanup = FALSE
    )
    info <- paste(result$output, collapse = "\n")

    expect_equal(result$status, 0L, info = info)
    expect_true(any(grepl(marker, result$output, fixed = TRUE)), info = info)
  }
})

test_that("full-suite runner requires clean direct shard exits", {
  runner <- testthat::test_path("..", "testthat.R")
  skip_if_not(file.exists(runner), "source test runner unavailable")
  text <- paste(readLines(runner, warn = FALSE), collapse = "\n")

  expect_false(grepl("setsid", text, fixed = TRUE))
  expect_false(grepl("137L", text, fixed = TRUE))
})
