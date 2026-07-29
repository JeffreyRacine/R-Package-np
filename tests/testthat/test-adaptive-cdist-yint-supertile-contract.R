library(npRmpi)

locate_adaptive_cdist_supertile_source <- function() {
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

adaptive_cdist_supertile_body <- function(lines) {
  start <- grep(
    "^static int np_conditional_distribution_cvls_lp_adap_block_stream\\(",
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

test_that("adaptive cdist CVLS reuses one Y-integral tile across three X blocks", {
  src_file <- locate_adaptive_cdist_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- adaptive_cdist_supertile_body(readLines(src_file, warn = FALSE))

  expect_match(body, "double **xblocks[3] = {NULL, NULL, NULL};", fixed = TRUE)
  expect_match(body, "xblocks[1] = np_optional_tmatd(", fixed = TRUE)
  expect_match(body, "xblocks[2] = np_optional_tmatd(", fixed = TRUE)
  expect_match(body, "group_width = 2;", fixed = TRUE)
  expect_match(body, "group_width = 3;", fixed = TRUE)
  expect_match(body, "i0 += group_width*block_size", fixed = TRUE)
  expect_match(body, "for(g = 0; g < group_width; g++)", fixed = TRUE)
  expect_match(
    body,
    "np_conditional_y_eval_from_ctx(&yintctx,",
    fixed = TRUE
  )
  expect_match(body, "xblocks[g][0]", fixed = TRUE)
  expect_match(body, "*cv += block_sum[g];", fixed = TRUE)
})

test_that("adaptive cdist CVLS supertile remains bounded and fallback-safe", {
  src_file <- locate_adaptive_cdist_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- adaptive_cdist_supertile_body(readLines(src_file, warn = FALSE))

  expect_match(body, "int group_width = 1;", fixed = TRUE)
  expect_match(body, "if(num_train > block_size)", fixed = TRUE)
  expect_match(body, "NPConditionalXRowCtx xctx = {0};", fixed = TRUE)
  expect_match(body, "NPConditionalYRowCtx yintctx = {0};", fixed = TRUE)
  expect_match(body, "np_blas_dgemm_tn_int(", fixed = TRUE)
  expect_false(grepl("alloc_tmatd\\(num_train, num_train\\)", body))
  expect_false(grepl("malloc|calloc|realloc", body))
})
