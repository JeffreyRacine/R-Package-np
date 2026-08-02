locate_objective_pack_source <- function() {
  candidates <- c(
    test_path("..", "..", "src", "jksum.c"),
    test_path("..", "..", "..", "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum.c"),
    file.path(getwd(), "src", "jksum.c"),
    file.path(getwd(), "..", "src", "jksum.c")
  )
  hits <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(hits)) normalizePath(hits[[1L]]) else NULL
}

test_that("generic LP objectives reuse only a matching immutable outer pack", {
  src_file <- locate_objective_pack_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  source <- paste(readLines(src_file, warn = FALSE), collapse = "\n")

  expect_match(source, "const double \\*data;")
  expect_match(source, "} NP_OuterPackCtx;", fixed = TRUE)
  expect_match(source, "outer_pack_ctx->matrix_W == matrix_W")
  expect_match(source, "outer_pack_ctx->ncol_W == ncol_W")
  expect_match(source, "outer_pack_ctx->ncol_Y == ncol_Y")
  expect_match(source, "outer_pack_ctx->num_weights == num_xt")
  expect_match(source, "outer_pack_ctx->symmetric == symmetric")
  expect_match(source, "blas_Apack_owned")
  expect_match(source, "if\\(blas_Apack_owned != NULL\\) free\\(blas_Apack_owned\\)")
})

test_that("objective pack reuse is scoped away from adaptive, tree, and reduced rows", {
  src_file <- locate_objective_pack_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  source <- paste(readLines(src_file, warn = FALSE), collapse = "\n")

  expect_match(
    source,
    "if\\(\\(BANDWIDTH_reg != BW_ADAP_NN\\) &&\\s*!ks_tree_use &&\\s*np_reg_cv_use_symmetric_dropone_path"
  )
  expect_match(
    source,
    "Pack it once at objective scope;\\s*reduced triangular rows advance the source pointers and cannot reuse it"
  )
  expect_match(source, "&objective_pack_ctx,", fixed = TRUE)
  expect_gte(lengths(regmatches(
    source, gregexpr("&objective_pack_ctx,", source, fixed = TRUE)
  )), 2L)
  expect_match(source, "(outer_pack_ctx->data != NULL)", fixed = TRUE)
  expect_match(source, "if\\(objective_Apack != NULL\\) free\\(objective_Apack\\)")
})
