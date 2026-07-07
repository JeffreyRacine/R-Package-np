test_that("npreghat contracts pass in an isolated MPI subprocess", {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess npreghat contracts")

  fixture <- normalizePath(
    testthat::test_path("fixtures", "npreghat-contracts-subprocess.R"),
    winslash = "/",
    mustWork = TRUE
  )
  for (block in seq_len(10L)) {
    # Each block owns its own subprocess MPI session; a success marker makes
    # the known local MPI teardown 137 harmless for this contract.
    res <- npRmpi_run_rscript_subprocess(
      lines = c(
        sprintf("source(%s)", shQuote(fixture))
      ),
      timeout = 90L,
      env = c(env, sprintf("NPREGHAT_FIXTURE_BLOCK=%d", block))
    )

    info <- paste(res$output, collapse = "\n")
    expect_true(res$status %in% c(0L, 137L), info = info)
    expect_true(
      any(grepl(sprintf("NPREGHAT_FIXTURE_BLOCK_OK %d", block),
                res$output,
                fixed = TRUE)),
      info = info
    )
  }
})
