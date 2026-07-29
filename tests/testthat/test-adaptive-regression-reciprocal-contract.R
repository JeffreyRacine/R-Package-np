library(np)

locate_adaptive_regression_reciprocal_source <- function() {
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

adaptive_regression_reciprocal_body <- function(
    lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_length(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("adaptive regression reciprocal workspace is objective-owned and bounded", {
  src_file <- locate_adaptive_regression_reciprocal_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  text <- paste(lines, collapse = "\n")
  workspace <- adaptive_regression_reciprocal_body(
    lines,
    "^static void np_adaptive_bandwidth_reciprocal_workspace_init\\(",
    "^static NP_NOINLINE NPRegCvLpResult np_regression_cv_lp_basis_adaptive_blas\\("
  )
  owner <- adaptive_regression_reciprocal_body(
    lines,
    "^static NP_NOINLINE NPRegCvLpResult np_regression_cv_lp_basis_adaptive_blas\\(",
    "^static NPRegCvLpResult np_regression_cv_lp_objective\\("
  )

  expect_match(
    text,
    "NPAdaptiveBandwidthReciprocalWorkspace;",
    fixed = TRUE
  )
  expect_match(
    workspace,
    "reciprocal_count = ((size_t)ndim + 1)*(size_t)n;",
    fixed = TRUE
  )
  expect_match(
    workspace,
    "workspace->product_reciprocal =",
    fixed = TRUE
  )
  expect_match(
    workspace,
    "if(!(h > 0.0) || !isfinite(h) || !isfinite(reciprocal))",
    fixed = TRUE
  )
  expect_false(grepl("n\\s*\\*\\s*n", workspace))
  expect_false(grepl("malloc|calloc|realloc", workspace))

  prepare_calls <- gregexpr(
    "np_adaptive_bandwidth_reciprocal_workspace_prepare\\(",
    text,
    perl = TRUE
  )[[1L]]
  expect_length(prepare_calls[prepare_calls > 0L], 2L)
  expect_match(
    owner,
    "np_adaptive_bandwidth_reciprocal_workspace_prepare(",
    fixed = TRUE
  )
  expect_match(
    owner,
    "if((kernel_c[l] < 1) || (kernel_c[l] > 3))",
    fixed = TRUE
  )
  expect_match(
    owner,
    "reciprocal_eligible =",
    fixed = TRUE
  )
  expect_match(
    owner,
    "(num_reg_unordered == 0) &&",
    fixed = TRUE
  )
  expect_match(
    owner,
    "(num_reg_ordered == 0) &&",
    fixed = TRUE
  )
  expect_match(
    owner,
    "(!int_cker_bound_extern))",
    fixed = TRUE
  )
  expect_match(
    owner,
    "allocation_count = weighted_count;",
    fixed = TRUE
  )
  expect_match(
    owner,
    "weighted_design + weighted_count,",
    fixed = TRUE
  )
  expect_match(
    owner,
    "reciprocal_workspace.ready ?",
    fixed = TRUE
  )
  expect_match(
    owner,
    "np_adaptive_bandwidth_reciprocal_workspace_clear(&reciprocal_workspace);",
    fixed = TRUE
  )
})

test_that("reciprocal row use is optional and conditional arithmetic stays unchanged", {
  src_file <- locate_adaptive_regression_reciprocal_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  helper <- adaptive_regression_reciprocal_body(
    lines,
    "^static int NP_NOINLINE np_accel_gauss_adaptive_higher_row_try\\(",
    "^/\\*"
  )
  conditional <- adaptive_regression_reciprocal_body(
    lines,
    "^static int .*np_conditional_xrow_from_ctx_impl\\(",
    "^static int np_conditional_xrow_from_ctx\\("
  )
  helper_compact <- gsub("[[:space:]]+", " ", helper)
  conditional_compact <- gsub("[[:space:]]+", " ", conditional)

  expect_match(
    helper_compact,
    "const int use_reciprocal = (bandwidth_reciprocal != NULL)",
    fixed = TRUE
  )
  expect_match(
    helper,
    "vDSP_vmulD(np_accel_gauss_tmp, 1,",
    fixed = TRUE
  )
  expect_match(
    helper,
    "vDSP_vdivD(bandwidth[d], 1, np_accel_gauss_tmp, 1,",
    fixed = TRUE
  )
  expect_match(
    helper,
    "vDSP_vdivD(np_accel_gauss_poly, 1, out, 1,",
    fixed = TRUE
  )
  expect_false(grepl("malloc|calloc|realloc", helper))
  expect_match(
    conditional_compact,
    "ctx->matrix_bandwidth_x, NULL, NULL, num_reg_continuous_extern",
    fixed = TRUE
  )
})
