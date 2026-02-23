test_that("npregiv ridging loops avoid superassignment", {
  src_path <- testthat::test_path("..", "..", "R", "npregiv.R")
  expect_true(file.exists(src_path))
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")
  expect_false(grepl("<<-", src, fixed = TRUE))
})
