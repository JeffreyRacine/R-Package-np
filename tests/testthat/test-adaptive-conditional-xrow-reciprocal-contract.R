locate_mpi_conditional_xrow_source <- function() {
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

test_that("MPI adaptive conditional reciprocals remain an isolated sidecar", {
  path <- locate_mpi_conditional_xrow_source()
  skip_if(is.null(path), "package C sources unavailable in this test context")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(
    source,
    "NPConditionalXRowReciprocalCache *reciprocal_cache;",
    fixed = TRUE
  )
  expect_match(
    source,
    "if(ctx->reciprocal_cache != NULL) free(ctx->reciprocal_cache);",
    fixed = TRUE
  )
  expect_match(source, "(np_glp_cv_cache.nterms < 4)", fixed = TRUE)
  expect_match(source, "(num_train <= 0) || (ndim < 2)", fixed = TRUE)
  expect_match(
    source,
    "reciprocal_count = ((size_t)ndim + 1)*(size_t)num_train;",
    fixed = TRUE
  )
  expect_match(
    source,
    "cache->product_reciprocal =",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_accel_gauss_adaptive_higher_row_reciprocal_try(",
    fixed = TRUE
  )
  guarded_helper <- paste(
    "#if NP_ACCEL_GAUSS_COMPILED",
    "static int NP_NOINLINE",
    "np_accel_gauss_adaptive_higher_row_reciprocal_try(",
    sep = "\n"
  )
  guarded_occurrences <- gregexpr(
    guarded_helper,
    source,
    fixed = TRUE
  )[[1L]]
  expect_length(guarded_occurrences[guarded_occurrences > 0L], 2L)

  prepare_start <- regexpr(
    "static int np_conditional_xrow_ctx_prepare",
    source,
    fixed = TRUE
  )
  row_start <- regexpr(
    "static int NP_NOINLINE NP_HOT_ALIGN np_conditional_xrow_from_ctx_impl",
    source,
    fixed = TRUE
  )
  expect_gt(prepare_start, 0L)
  expect_gt(row_start, prepare_start)
  prepare <- substr(source, prepare_start, row_start - 1L)
  expect_false(grepl(
    "np_conditional_xrow_reciprocal_cache_try(ctx)",
    prepare,
    fixed = TRUE
  ))

  row_end <- regexpr(
    "static int np_conditional_xrow_from_ctx(",
    source,
    fixed = TRUE
  )
  expect_gt(row_end, row_start)
  row <- substr(source, row_start, row_end - 1L)
  expect_match(row, "(eval_idx == 0)", fixed = TRUE)
  expect_match(
    row,
    "(void)np_conditional_xrow_reciprocal_cache_try(ctx);",
    fixed = TRUE
  )
  expect_match(
    row,
    "ctx->reciprocal_cache->storage",
    fixed = TRUE
  )
  expect_match(
    row,
    "ctx->reciprocal_cache->product_reciprocal",
    fixed = TRUE
  )
})

test_that("MPI regression retains its original division-only Gaussian helper", {
  path <- locate_mpi_conditional_xrow_source()
  skip_if(is.null(path), "package C sources unavailable in this test context")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  helper_start <- regexpr(
    "static int NP_NOINLINE np_accel_gauss_adaptive_higher_row_try(",
    source,
    fixed = TRUE
  )
  helper_end <- regexpr(
    "static int NP_NOINLINE np_accel_gauss_product_kind(",
    source,
    fixed = TRUE
  )
  expect_gt(helper_start, 0L)
  expect_gt(helper_end, helper_start)
  helper <- substr(source, helper_start, helper_end - 1L)

  expect_false(grepl("bandwidth_reciprocal", helper, fixed = TRUE))
  expect_match(
    helper,
    "np_accel_vdivD(bandwidth[d], 1, np_accel_gauss_tmp, 1,",
    fixed = TRUE
  )
})
