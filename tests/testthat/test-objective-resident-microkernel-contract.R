test_that("fixed width-three LP CV accumulates one triangle before solving", {
  src_dir <- file.path(testthat::test_path("..", ".."), "src")
  row_file <- file.path(src_dir, "jksum_lp_row.c")
  sum_file <- file.path(src_dir, "jksum.c")
  skip_if_not(file.exists(row_file) && file.exists(sum_file),
              "package C sources unavailable")

  row_text <- paste(readLines(row_file, warn = FALSE), collapse = "\n")
  sum_text <- paste(readLines(sum_file, warn = FALSE), collapse = "\n")
  expect_match(row_text, "attribute_hidden void\\s+np_lp_mirror_dense_moments_row3")
  expect_match(row_text, "s\\[3\\] = s\\[1\\];")
  expect_match(row_text, "s\\[6\\] = s\\[2\\];")
  expect_match(row_text, "s\\[7\\] = s\\[5\\];")
  expect_false(grepl("sj3 += wb1*b0;", row_text, fixed = TRUE))
  expect_false(grepl("sj6 += wb2*b0;", row_text, fixed = TRUE))
  expect_false(grepl("sj7 += wb2*b1;", row_text, fixed = TRUE))
  expect_match(sum_text, "np_lp_mirror_dense_moments_row3\\(moments, num_obs\\)")
})
