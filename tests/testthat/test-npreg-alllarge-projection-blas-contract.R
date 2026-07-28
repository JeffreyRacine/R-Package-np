library(np)

locate_alllarge_projection_source <- function() {
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

test_that("regression all-large projection BLAS gate is narrow and bounded", {
  src_file <- locate_alllarge_projection_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  source <- paste(readLines(src_file, warn = FALSE), collapse = "\n")

  expect_match(source, "#define NP_REG_ALLLARGE_BLAS_MIN_ROWS 2048",
               fixed = TRUE)
  expect_match(source, "#define NP_REG_ALLLARGE_BLAS_MIN_TERMS 6",
               fixed = TRUE)
  expect_match(source, "#define NP_REG_ALLLARGE_PROJECTION_MAX_BYTES",
               fixed = TRUE)
  expect_match(source, "((size_t)8*(size_t)1024*(size_t)1024)",
               fixed = TRUE)
  expect_match(
    source,
    "return np_mseries_accelerate_enabled_cache &&",
    fixed = TRUE
  )
  expect_match(source, "#if NP_ACCEL_GAUSS_COMPILED", fixed = TRUE)
})

test_that("regression all-large projection uses bounded BLAS with scalar fallback", {
  src_file <- locate_alllarge_projection_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  function_start <- grep(
    "^static NPRegCvLpResult np_regression_cv_lp_objective\\(",
    lines
  )
  start <- grep(
    "const int basis_stride = np_glp_cv_cache.basis_stride;",
    lines,
    fixed = TRUE
  )
  start <- start[start > function_start][1L]
  stop <- grep("np_fastcv_alllarge_hits++;", lines, fixed = TRUE)
  stop <- stop[stop > start][1L]
  expect_length(function_start, 1L)
  expect_length(start, 1L)
  expect_true(is.finite(stop))
  body <- paste(lines[start:stop], collapse = "\n")
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(body, "np_blas_gram_int(", fixed = TRUE)
  expect_match(body, "np_blas_dgemv_t_int(", fixed = TRUE)
  expect_match(body, "np_blas_dgemv_n_int(", fixed = TRUE)
  expect_match(body, "np_blas_project_inverse_block_int(", fixed = TRUE)
  expect_match(
    compact,
    "projection_max_elements/(size_t)k",
    fixed = TRUE
  )
  expect_match(
    compact,
    "yhat_all = (double *)malloc((size_t)num_obs*sizeof(double));",
    fixed = TRUE
  )
  expect_match(
    compact,
    "hdiag_all = (double *)malloc((size_t)num_obs*sizeof(double));",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(i = 0; i < k; i++){ inverse_workspace.rhs[i] = 0.0;",
    fixed = TRUE
  )
  expect_match(
    body,
    "for(int b = 0; b < k; b++)",
    fixed = TRUE
  )
  expect_false(grepl("num_obs\\s*\\*\\s*num_obs", body))
  expect_false(grepl("num_obs\\)\\s*\\*\\s*\\(size_t\\)num_obs", body))
})
