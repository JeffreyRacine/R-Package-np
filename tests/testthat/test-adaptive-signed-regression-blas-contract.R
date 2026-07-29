library(np)

locate_adaptive_signed_regression_source <- function() {
  candidates <- c(
    test_path("..", "..", "src", "jksum.c"),
    test_path("..", "..", "..", "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum.c"),
    file.path(getwd(), "src", "jksum.c"),
    file.path(getwd(), "..", "src", "jksum.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) NULL else hits[[1L]]
}

test_that("adaptive signed regression uses the isolated weighted-design owner", {
  src_file <- locate_adaptive_signed_regression_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  start <- grep(
    "^static NP_NOINLINE NPRegCvLpResult np_regression_cv_lp_basis_adaptive_blas\\(",
    lines
  )
  stop <- grep(
    "^static NPRegCvLpResult np_regression_cv_lp_objective\\(",
    lines
  )
  expect_length(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  body <- paste(lines[start:(stop - 1L)], collapse = "\n")
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(
    compact,
    "if(((kernel_c[l] < 0) || (kernel_c[l] > 4)) && (kernel_c[l] != 8))",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_accel_gauss_adaptive_row_try(kernel_c,",
    fixed = TRUE
  )
  expect_match(
    body,
    "kernel_weighted_sum_np_ctx_ex(kernel_c,",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_accel_gauss_adaptive_higher_row_try(kernel_c,",
    fixed = TRUE
  )
  expect_match(
    body,
    "F77_CALL(dgemm)(&trans_t,",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_blas_dgemv_t_int(num_obs,",
    fixed = TRUE
  )
})
