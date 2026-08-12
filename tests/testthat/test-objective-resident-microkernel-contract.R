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

test_that("LP row NEON is an Apple-arm64 moving-row specialization", {
  row_file <- file.path(testthat::test_path("..", ".."),
                        "src", "jksum_lp_row.c")
  skip_if_not(file.exists(row_file), "package C sources unavailable")

  row_text <- paste(readLines(row_file, warn = FALSE), collapse = "\n")
  expect_match(
    row_text,
    "defined\\(__aarch64__\\).*defined\\(NP_USE_ACCELERATE_GAUSS\\)"
  )
  expect_match(row_text, "#define NP_LP_ROW_NEON 1", fixed = TRUE)
  expect_match(row_text, "vfmaq_f64(vld1q_f64(ti), vw, eval_y01)",
               fixed = TRUE)
  expect_match(row_text, "vfmaq_f64(vld1q_f64(si), vw, eval_outer01)",
               fixed = TRUE)
  expect_match(row_text, "tj0 += wb0*yi;", fixed = TRUE)
  expect_match(row_text, "#else\n      ti[0] += w*eval_ybasis[0];",
               fixed = TRUE)
})

test_that("fixed widths five and six retain only the unique Gram triangle", {
  src_dir <- file.path(testthat::test_path("..", ".."), "src")
  row_file <- file.path(src_dir, "jksum_lp_row.c")
  sum_file <- file.path(src_dir, "jksum.c")
  skip_if_not(file.exists(row_file) && file.exists(sum_file),
              "package C sources unavailable")

  row_text <- paste(readLines(row_file, warn = FALSE), collapse = "\n")
  sum_text <- paste(readLines(sum_file, warn = FALSE), collapse = "\n")
  for (width in 5:6) {
    expect_match(
      row_text,
      paste0("NP_LP_DEFINE_SYMMETRIC_RESIDENT_WIDTH\\(", width, "\\)")
    )
    expect_false(grepl(
      paste0("NP_LP_DEFINE_RESIDENT_WIDTH(", width, ")"),
      row_text,
      fixed = TRUE
    ))
  }
  expect_match(row_text, "for\\(b = a; b < \\(WIDTH\\); b\\+\\+\\)")
  expect_match(row_text, "a \\+ 1 < \\(WIDTH\\)")
  expect_match(row_text, "if\\(a < \\(WIDTH\\)\\)")
  expect_match(
    sum_text,
    "!use_sparse_tree && \\(\\(nterms == 5\\) \\|\\| \\(nterms == 6\\)\\)"
  )
  expect_match(sum_text, "NP_REG_CV_LP_RESIDENT_MAX_TERMS = 6",
               fixed = TRUE)
})
