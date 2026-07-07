test_that("npregiv ridging loops avoid superassignment", {
  src_path <- testthat::test_path("..", "..", "R", "npregiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")
  lines <- grep("<<-", strsplit(src, "\n", fixed = TRUE)[[1L]], value = TRUE, fixed = TRUE)
  expect_false(any(grepl(
    "ridge|doridge|coef\\.mat|mean\\.loo|ghat|WzkWz",
    lines
  )))
})

test_that("npregiv uses shared seed enter helper", {
  src_path <- testthat::test_path("..", "..", "R", "npregiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")
  expect_true(grepl("\\.np_seed_enter\\(random\\.seed\\)", src))
  expect_false(grepl("exists\\(\"\\.Random\\.seed\"", src))
  expect_true(grepl("\\.np_seed_exit\\(seed\\.state, remove_if_absent = TRUE\\)", src))
})
