library(np)

locate_adaptive_xrow_blas_source <- function() {
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

test_that("adaptive X-row context owns one bounded weighted-design slab", {
  src_file <- locate_adaptive_xrow_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  all_source <- paste(lines, collapse = "\n")
  prepare <- np_test_extract_c_function(
    lines, "np_conditional_xrow_ctx_prepare_impl"
  )
  clear <- np_test_extract_c_function(lines, "np_conditional_xrow_ctx_clear")
  compact <- np_test_compact_source(prepare)

  expect_match(all_source, "double *weighted_design;", fixed = TRUE)
  expect_match(
    prepare,
    "np_glp_cv_cache.nterms >= 4",
    fixed = TRUE
  )
  expect_match(
    prepare,
    "np_apple_conditional_x_weighted_blas_profitable(",
    fixed = TRUE
  )
  expect_match(
    compact,
    "ctx->weighted_design = (double *)malloc(weighted_count*sizeof(double));",
    fixed = TRUE
  )
  expect_match(
    clear,
    "if(ctx->weighted_design != NULL) free(ctx->weighted_design);",
    fixed = TRUE
  )
  expect_false(grepl("num_train\\*num_train", prepare))
})

test_that("adaptive X rows preserve signed solve deletion and scalar fallback", {
  src_file <- locate_adaptive_xrow_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  body <- np_test_extract_c_function(
    lines, "np_conditional_xrow_legacy_influence"
  )
  dispatch <- np_test_extract_c_function(
    lines, "np_conditional_xrow_from_ctx_impl"
  )
  compact <- np_test_compact_source(body)

  expect_match(
    body,
    "weighted_row[j] = basis_row[j]*ctx->kw[j];",
    fixed = TRUE
  )
  expect_match(body, "F77_CALL(dgemm)", fixed = TRUE)
  expect_match(body, "F77_CALL(dgemv)", fixed = TRUE)
  expect_match(
    body,
    "row_out[orig_j] = ctx->kw[j]*ctx->mean_row[j];",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_lp_full_row_workspace_solve(&ctx->full_row_workspace",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_lp_delete_denominator(row_out[eval_idx], &denominator)",
    fixed = TRUE
  )
  expect_match(
    compact,
    "if(ctx->weighted_design != NULL){ const char trans_t",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(l = 0; l < k; l++) for(j = 0; j < k; j++) ctx->full_row_workspace.gram[l + j*k] = 0.0;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(j = 0; j < num_train; j++){ double zju = 0.0;",
    fixed = TRUE
  )
  expect_match(
    dispatch,
    "np_conditional_xrow_legacy_influence(",
    fixed = TRUE
  )
  expect_match(
    dispatch,
    "np_regression_xrow_canonical_influence(",
    fixed = TRUE
  )
})

test_that("adaptive objective families share the context without a proof graph", {
  src_file <- locate_adaptive_xrow_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  all_source <- paste(lines, collapse = "\n")
  cv <- np_test_extract_c_function(
    lines, "np_conditional_density_cvml_lp_stream"
  )

  expect_false(grepl("np_shadow_", all_source, fixed = TRUE))
  expect_match(
    cv,
    "if(use_row_ctx){\n    if(np_conditional_xrow_ctx_prepare_ctx(",
    fixed = TRUE
  )
  expect_equal(
    sum(gregexpr("np_conditional_xrow_from_ctx\\(&xctx", cv)[[1L]] > 0L),
    1L
  )
  expect_equal(
    sum(gregexpr("np_conditional_yrow_from_ctx\\(&yctx", cv)[[1L]] > 0L),
    1L
  )
})
