test_that("profile demo launcher captures the R CMD BATCH transcript explicitly", {
  makefile <- testthat::test_path("..", "..", "demo", "makefile")
  skip_if_not(file.exists(makefile), "demo makefile unavailable in installed test context")
  lines <- readLines(makefile, warn = FALSE)
  text <- paste(lines, collapse = "\n")

  expect_match(
    text,
    "batch_file=\\$\\$\\{out_file\\}\\.batch;",
    perl = TRUE
  )
  expect_match(
    text,
    "R CMD BATCH --no-save \\.\\./\\$\\$\\{d\\}_npRmpi_profile\\.R \"\\$\\$batch_file\" > \"\\$\\$tmp_file\" 2>&1",
    perl = TRUE
  )
  expect_match(
    text,
    "\\{ \\$\\(PROVENANCE_CMD\\); cat \"\\$\\$batch_file\"; if test -s \"\\$\\$tmp_file\"; then cat \"\\$\\$tmp_file\"; fi; \\} > \"\\$\\$out_file\";",
    perl = TRUE
  )
})
