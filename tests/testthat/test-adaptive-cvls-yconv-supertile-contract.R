library(npRmpi)

locate_adaptive_cvls_supertile_source <- function() {
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

adaptive_cvls_supertile_body <- function(lines) {
  start <- grep(
    "^static int np_conditional_density_cvls_lp_adap_block_stream\\(",
    lines
  )
  stop <- grep(
    "^static int np_conditional_density_cvls_categorical_profile_stream\\(",
    lines
  )
  expect_length(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("adaptive CVLS reuses one Y-convolution tile across two X blocks", {
  src_file <- locate_adaptive_cvls_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- adaptive_cvls_supertile_body(readLines(src_file, warn = FALSE))

  expect_match(body, "double **full_blocks[2] = {NULL, NULL};", fixed = TRUE)
  expect_match(body, "i0 += 2*block_size", fixed = TRUE)
  expect_match(body, "for(g = 0; g < 2; g++)", fixed = TRUE)
  expect_match(
    body,
    "np_conditional_yrow_from_ctx(&yconvctx, j, shared_y[jj])",
    fixed = TRUE
  )
  expect_match(
    body,
    "double * const quad_cross = loo_work[0];",
    fixed = TRUE
  )
  expect_match(body, "full_blocks[g][0]", fixed = TRUE)
  expect_match(body, "*cv += quad[g] - 2.0*lin[g];", fixed = TRUE)
})

test_that("adaptive CVLS supertile lowers workspace without changing owners", {
  src_file <- locate_adaptive_cvls_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- adaptive_cvls_supertile_body(readLines(src_file, warn = FALSE))

  expect_match(body, "NPConditionalXRowCtx xctx = {0};", fixed = TRUE)
  expect_match(body, "NPConditionalYRowCtx yctx = {0}, yconvctx = {0};",
               fixed = TRUE)
  expect_match(body, "np_lp_delete_smoother_row(", fixed = TRUE)
  expect_match(body, "np_blas_dgemm_tn_int(", fixed = TRUE)
  expect_false(grepl("alloc_vecd\\(block_size\\*block_size\\)", body))
  expect_false(grepl("malloc|calloc|realloc", body))
})
