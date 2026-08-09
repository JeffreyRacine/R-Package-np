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
  parent_runner <- testthat::test_path("..", "testthat.R")
  shard_runner <- testthat::test_path(
    "..", "validation", "run_full_mpi_test_shard.R"
  )
  skip_if_not(
    file.exists(parent_runner) && file.exists(shard_runner),
    "source test runners unavailable"
  )
  parent_text <- paste(readLines(parent_runner, warn = FALSE), collapse = "\n")
  shard_text <- paste(readLines(shard_runner, warn = FALSE), collapse = "\n")

  expect_false(grepl("setsid", parent_text, fixed = TRUE))
  expect_false(grepl("137L", parent_text, fixed = TRUE))
  expect_match(shard_text, "npRmpi_run_full_test_shard <- function(args)",
               fixed = TRUE)
  expect_match(shard_text, "npRmpi_shard_has_failures(result_summary)",
               fixed = TRUE)
  expect_match(shard_text, "npRmpi.quit(force = TRUE)", fixed = TRUE)
  expect_match(shard_text, "file.rename(witness_tmp, witness)", fixed = TRUE)
  expect_match(
    shard_text,
    "quit(save = \"no\", status = exit_status, runLast = FALSE)",
    fixed = TRUE
  )
  expect_false(grepl("writeLines(marker, witness)", shard_text, fixed = TRUE))
})
