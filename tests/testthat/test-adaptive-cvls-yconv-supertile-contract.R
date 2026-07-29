library(np)

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

adaptive_cvls_width3_body <- function(lines) {
  start <- tail(grep(
    "^np_conditional_density_cvls_lp_adap_block3_stream\\(",
    lines
  ), 1L)
  stop <- grep(
    "^#undef NP_CDENS_ADAP_WIDTH3_NOINLINE$",
    lines
  )
  expect_length(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

adaptive_cvls_width4_body <- function(lines) {
  start <- tail(grep(
    "^np_conditional_density_cvls_lp_adap_block4_stream\\(",
    lines
  ), 1L)
  stop <- grep(
    "^#undef NP_CDENS_ADAP_WIDTH4_NOINLINE$",
    lines
  )
  expect_length(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("adaptive CVLS dispatches width three inside its width-two owner", {
  src_file <- locate_adaptive_cvls_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- adaptive_cvls_supertile_body(readLines(src_file, warn = FALSE))

  expect_match(body, "double **full_blocks[2] = {NULL, NULL};", fixed = TRUE)
  expect_match(
    body,
    "np_conditional_density_cvls_lp_adap_block3_stream(",
    fixed = TRUE
  )
  expect_match(
    body,
    "width3_status != NP_CDENS_ADAP_WIDTH3_UNAVAILABLE",
    fixed = TRUE
  )
  expect_match(body, "return width3_status;", fixed = TRUE)
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

test_that("adaptive CVLS width three is an isolated pass-saving sibling", {
  src_file <- locate_adaptive_cvls_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- adaptive_cvls_width3_body(readLines(src_file, warn = FALSE))

  expect_match(
    body,
    "double **full_blocks[3] = {NULL, NULL, NULL};",
    fixed = TRUE
  )
  expect_match(body, "width3_passes >= width2_passes", fixed = TRUE)
  expect_match(
    body,
    "return NP_CDENS_ADAP_WIDTH3_UNAVAILABLE;",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_conditional_density_cvls_lp_adap_block4_stream(",
    fixed = TRUE
  )
  expect_match(
    body,
    "width4_status != NP_CDENS_ADAP_WIDTH4_UNAVAILABLE",
    fixed = TRUE
  )
  expect_match(
    body,
    "full_blocks[2] = np_optional_tmatd(num_obs, block_size);",
    fixed = TRUE
  )
  incumbent_alloc <- regexpr(
    "loo_work = alloc_tmatd(num_obs, block_size);",
    body,
    fixed = TRUE
  )[[1L]]
  optional_alloc <- regexpr(
    "full_blocks[2] = np_optional_tmatd(num_obs, block_size);",
    body,
    fixed = TRUE
  )[[1L]]
  expect_gt(incumbent_alloc, 0L)
  expect_gt(optional_alloc, incumbent_alloc)
  expect_match(body, "i0 += 3*block_size", fixed = TRUE)
  expect_match(body, "for(g = 0; g < 3; g++)", fixed = TRUE)
  expect_match(
    body,
    "double * const quad_cross = loo_work[0];",
    fixed = TRUE
  )
  expect_match(body, "*cv += quad[g] - 2.0*lin[g];", fixed = TRUE)
})

test_that("adaptive CVLS width four is an isolated pass-saving sibling", {
  src_file <- locate_adaptive_cvls_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- adaptive_cvls_width4_body(readLines(src_file, warn = FALSE))

  expect_match(
    body,
    "double **full_blocks[4] = {NULL, NULL, NULL, NULL};",
    fixed = TRUE
  )
  expect_match(body, "width4_passes >= width3_passes", fixed = TRUE)
  expect_match(
    body,
    "return NP_CDENS_ADAP_WIDTH4_UNAVAILABLE;",
    fixed = TRUE
  )
  incumbent_alloc <- regexpr(
    "loo_work = alloc_tmatd(num_obs, block_size);",
    body,
    fixed = TRUE
  )[[1L]]
  optional_alloc <- regexpr(
    "optional_blocks = np_optional_tmatd(num_obs, 2*block_size);",
    body,
    fixed = TRUE
  )[[1L]]
  expect_gt(incumbent_alloc, 0L)
  expect_gt(optional_alloc, incumbent_alloc)
  expect_match(body, "full_blocks[2] = optional_blocks;", fixed = TRUE)
  expect_match(
    body,
    "full_blocks[3] = optional_blocks + block_size;",
    fixed = TRUE
  )
  expect_match(body, "i0 += 4*block_size", fixed = TRUE)
  expect_match(body, "for(g = 0; g < 4; g++)", fixed = TRUE)
  expect_match(
    body,
    "double * const quad_cross = loo_work[0];",
    fixed = TRUE
  )
  expect_match(body, "np_blas_dgemm_tn_int(", fixed = TRUE)
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
