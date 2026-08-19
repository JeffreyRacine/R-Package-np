test_that("fixed Gaussian support tracking is elided only after certification", {
  source_file <- file.path(testthat::test_path("..", ".."), "src", "jksum.c")
  skip_if_not(file.exists(source_file), "package C sources unavailable")

  text <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
  expect_match(text, "np_lp_fixed_gaussian_full_support_certified")
  expect_match(text, "kernel_c\\[l\\] != CK_GAUSS2")
  expect_match(text, "operator\\[l\\] != OP_NORMAL")
  expect_match(text, "num_reg_unordered != 0")
  expect_match(text, "num_reg_ordered != 0")
  expect_match(text, "int_cker_bound_extern")
  expect_match(text, "log\\(DBL_MIN\\) \\+ 8\\.0")
  expect_match(
    text,
    "if\\(track_lowsupport &&\\s+np_lp_fixed_gaussian_full_support_certified"
  )
  expect_match(text, "track_lowsupport = 0;")
  expect_match(text, "if\\(track_lowsupport\\)\\{\\s+support_count =")
})
