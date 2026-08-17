test_that("fixed LP objectives retain no private low-support solve policy", {
  source_root <- file.path(testthat::test_path("..", ".."), "src")
  jksum_file <- file.path(source_root, "jksum.c")
  row_file <- file.path(source_root, "jksum_lp_row.c")
  row_header <- file.path(source_root, "jksum_lp_row.h")
  skip_if_not(
    all(file.exists(c(jksum_file, row_file, row_header))),
    "package C sources unavailable"
  )

  jksum_lines <- readLines(jksum_file, warn = FALSE)
  fixed_body <- np_test_extract_c_function(
    jksum_lines, "np_regression_cv_lp_basis_fixed"
  )
  row_source <- paste(
    c(readLines(row_file, warn = FALSE),
      readLines(row_header, warn = FALSE)),
    collapse = "\n"
  )

  expect_match(
    fixed_body,
    "np_lp_solve_workspace_solve_response(&solve_workspace,",
    fixed = TRUE
  )
  expect_false(grepl("dgelsy", fixed_body, fixed = TRUE))
  expect_false(grepl("lowsupport", fixed_body, fixed = TRUE))
  expect_false(grepl("support_count", fixed_body, fixed = TRUE))
  expect_false(grepl("track_lowsupport", row_source, fixed = TRUE))
  expect_false(grepl("support_weight", row_source, fixed = TRUE))
})
