locate_conditional_yrow_reciprocal_source <- function() {
  roots <- c(
    test_path("..", ".."),
    test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  )
  roots <- unique(roots[nzchar(roots)])
  roots <- roots[file.exists(file.path(roots, "src", "jksum.c"))]
  if (!length(roots)) return(NULL)
  file.path(roots[[1L]], "src", "jksum.c")
}

test_that("adaptive conditional response reciprocals are a lazy optional sidecar", {
  path <- locate_conditional_yrow_reciprocal_source()
  skip_if(is.null(path), "package C sources unavailable in this test context")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(
    source,
    "NPConditionalYRowReciprocalCache *reciprocal_cache;",
    fixed = TRUE
  )
  expect_match(source, "int reciprocal_cache_attempted;", fixed = TRUE)
  expect_match(
    source,
    "if(ctx->reciprocal_cache != NULL) free(ctx->reciprocal_cache);",
    fixed = TRUE
  )
  expect_match(
    source,
    "reciprocal_count = ((size_t)ndim + 1)*(size_t)num_train;",
    fixed = TRUE
  )
  expect_match(
    source,
    "(num_train < NP_CONDITIONAL_X_WEIGHTED_BLAS_MIN_ROWS)",
    fixed = TRUE
  )
  expect_match(source, "(num_var_unordered_extern != 0)", fixed = TRUE)
  expect_match(source, "(num_var_ordered_extern != 0)", fixed = TRUE)
  expect_match(source, "int_cyker_bound_extern", fixed = TRUE)
  expect_match(source, "(int_TREE_Y == NP_TREE_TRUE)", fixed = TRUE)

  prepare_start <- regexpr(
    "static int np_conditional_yrow_eval_ctx_prepare",
    source,
    fixed = TRUE
  )
  row_start <- regexpr(
    "static int np_conditional_yrow_from_ctx",
    source,
    fixed = TRUE
  )
  expect_gt(prepare_start, 0L)
  expect_gt(row_start, prepare_start)
  prepare <- substr(source, prepare_start, row_start - 1L)
  expect_false(grepl(
    "np_conditional_yrow_reciprocal_cache_try(ctx)",
    prepare,
    fixed = TRUE
  ))

  row_end <- regexpr(
    "static int np_conditional_y_eval_from_ctx",
    source,
    fixed = TRUE
  )
  expect_gt(row_end, row_start)
  row <- substr(source, row_start, row_end - 1L)
  expect_match(
    row,
    "(!ctx->reciprocal_cache_attempted)",
    fixed = TRUE
  )
  expect_match(
    row,
    "ctx->reciprocal_cache_attempted = 1;",
    fixed = TRUE
  )
  expect_match(
    row,
    "(void)np_conditional_yrow_reciprocal_cache_try(ctx);",
    fixed = TRUE
  )
  expect_match(
    row,
    "np_accel_gauss_adaptive_yrow_reciprocal_try(",
    fixed = TRUE
  )
  expect_match(
    row,
    "ctx->reciprocal_cache->workspace.product_reciprocal",
    fixed = TRUE
  )
})

test_that("adaptive response reciprocal admission retains both fallbacks", {
  path <- locate_conditional_yrow_reciprocal_source()
  skip_if(is.null(path), "package C sources unavailable in this test context")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")
  row_start <- regexpr(
    "static int np_conditional_yrow_from_ctx",
    source,
    fixed = TRUE
  )
  row_end <- regexpr(
    "static int np_conditional_y_eval_from_ctx",
    source,
    fixed = TRUE
  )
  expect_gt(row_start, 0L)
  expect_gt(row_end, row_start)
  row <- substr(source, row_start, row_end - 1L)

  reciprocal <- regexpr(
    "np_accel_gauss_adaptive_yrow_reciprocal_try(",
    row,
    fixed = TRUE
  )
  incumbent <- regexpr(
    "np_accel_gauss_adaptive_row_try(",
    row,
    fixed = TRUE
  )
  generic <- regexpr(
    "np_shadow_conditional_kernel_row(",
    row,
    fixed = TRUE
  )
  expect_gt(reciprocal, 0L)
  expect_gt(incumbent, reciprocal)
  expect_gt(generic, incumbent)
  expect_match(row, "if(!adaptive_gaussian_row &&", fixed = TRUE)
  expect_match(row, "np_conditional_push_bounds(", fixed = TRUE)
  expect_match(row, "np_conditional_pop_bounds(", fixed = TRUE)
})
