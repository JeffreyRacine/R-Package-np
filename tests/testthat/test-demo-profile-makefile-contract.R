test_that("profile demo launcher gives every MPI rank a distinct batch transcript", {
  makefile <- system.file(
    "demo_tools", "makefile", package = "npRmpi", mustWork = TRUE
  )
  lines <- readLines(makefile, warn = FALSE)
  text <- paste(lines, collapse = "\n")

  expect_match(
    text,
    "batch_base=\\$\\$\\{out_file\\}\\.batch;",
    perl = TRUE
  )
  expect_match(
    text,
    "PMI_RANK.*PMIX_RANK.*OMPI_COMM_WORLD_RANK",
    perl = TRUE
  )
  expect_match(
    text,
    "R CMD BATCH --no-save.*\\$\\$2\\.rank\\$\\$\\{rank\\}",
    perl = TRUE
  )
  expect_match(
    text,
    "master_batch=\\$\\$\\{batch_base\\}\\.rank0;",
    perl = TRUE
  )
})
