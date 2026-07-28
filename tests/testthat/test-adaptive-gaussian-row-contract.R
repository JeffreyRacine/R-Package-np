library(np)

locate_adaptive_gaussian_row_source <- function() {
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

adaptive_gaussian_row_source_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_length(start, 1L)
  expect_gte(length(stop), 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("adaptive Gaussian row fusion is narrow and memory bounded", {
  src_file <- locate_adaptive_gaussian_row_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  body <- adaptive_gaussian_row_source_body(
    lines,
    "^static int NP_NOINLINE np_accel_gauss_adaptive_row_try\\(",
    "^/\\*"
  )
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(
    compact,
    "n < NP_CONDITIONAL_X_WEIGHTED_BLAS_MIN_ROWS",
    fixed = TRUE
  )
  expect_match(
    compact,
    "(kernel_c[d] != 0) || (operator[d] != OP_NORMAL)",
    fixed = TRUE
  )
  expect_match(body, "np_accel_gauss_scratch_ensure(n)", fixed = TRUE)
  expect_match(
    body,
    "(np_accel_vdivD|vDSP_vdivD)\\(bandwidth\\[d\\]",
    perl = TRUE
  )
  expect_match(body, "(np_accel_vsqD|vDSP_vsqD)\\(", perl = TRUE)
  expect_match(body, "(np_accel_vaddD|vDSP_vaddD)\\(", perl = TRUE)
  expect_match(body, "vvexp(out, np_accel_gauss_arg", fixed = TRUE)
  expect_match(body, "if(divide_by_bandwidth_product)", fixed = TRUE)
  expect_match(body, "if(!isfinite(out[i]))", fixed = TRUE)
  expect_match(body, "memset(out, 0, (size_t)n*sizeof(double));", fixed = TRUE)
  expect_false(grepl("n\\s*\\*\\s*n", body))
  expect_false(grepl("malloc|calloc|realloc", body))
})

test_that("adaptive Gaussian row fusion has exactly two objective consumers", {
  src_file <- locate_adaptive_gaussian_row_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  text <- paste(lines, collapse = "\n")
  calls <- gregexpr(
    "np_accel_gauss_adaptive_row_try\\(",
    text,
    perl = TRUE
  )[[1L]]

  expect_length(calls[calls > 0L], 3L)

  regression <- adaptive_gaussian_row_source_body(
    lines,
    "^static NPRegCvLpResult np_regression_cv_lp_basis_adaptive_blas\\(",
    paste0(
      "^static (NPRegCvLpResult np_regression_cv_lp_objective|",
      "int np_distribution_cvls_ordered_profile_stream)\\("
    )
  )
  conditional <- adaptive_gaussian_row_source_body(
    lines,
    "^static int .*np_conditional_xrow_from_ctx_impl\\(",
    "^static int np_conditional_xrow_from_ctx\\("
  )

  expect_match(regression, "num_reg_unordered == 0", fixed = TRUE)
  expect_match(regression, "num_reg_ordered == 0", fixed = TRUE)
  expect_match(regression, "!int_cker_bound_extern", fixed = TRUE)
  expect_match(
    regression,
    "np_accel_gauss_adaptive_row_try(kernel_c,",
    fixed = TRUE
  )
  expect_match(
    regression,
    "kernel_weighted_sum_np_ctx_ex(kernel_c,",
    fixed = TRUE
  )

  expect_match(
    conditional,
    "BANDWIDTH_den_extern == BW_ADAP_NN",
    fixed = TRUE
  )
  expect_match(conditional, "num_reg_unordered_extern == 0", fixed = TRUE)
  expect_match(conditional, "num_reg_ordered_extern == 0", fixed = TRUE)
  expect_match(conditional, "!int_cxker_bound_extern", fixed = TRUE)
  expect_match(conditional, "int_TREE_X != NP_TREE_TRUE", fixed = TRUE)
  expect_match(
    conditional,
    "np_shadow_conditional_kernel_row_raw(ctx->kernel_cx,",
    fixed = TRUE
  )
})
