library(npRmpi)

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

adaptive_xrow_source_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_length(start, 1L)
  expect_gte(length(stop), 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("MPI adaptive X-row context owns one bounded rank-local slab", {
  src_file <- locate_adaptive_xrow_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  all_source <- paste(lines, collapse = "\n")
  prepare <- adaptive_xrow_source_body(
    lines,
    "^static int np_conditional_xrow_ctx_prepare\\(",
    "^static int (NP_NOINLINE )?(NP_HOT_ALIGN )?np_conditional_xrow_from_ctx_impl\\("
  )
  clear <- adaptive_xrow_source_body(
    lines,
    "^static void np_conditional_xrow_ctx_clear\\(",
    "^static int np_conditional_xrow_ctx_prepare\\("
  )
  compact <- gsub("[[:space:]]+", " ", prepare)

  expect_match(all_source, "double *weighted_design;", fixed = TRUE)
  expect_match(prepare, "np_glp_cv_cache.nterms >= 4", fixed = TRUE)
  expect_match(
    prepare,
    "np_conditional_x_weighted_blas_profitable(",
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
  expect_false(grepl("MPI_", prepare, fixed = TRUE))
})

test_that("MPI adaptive X rows retain solve deletion and scalar fallback", {
  src_file <- locate_adaptive_xrow_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  body <- adaptive_xrow_source_body(
    lines,
    "^static int (NP_NOINLINE )?(NP_HOT_ALIGN )?np_conditional_xrow_from_ctx_impl\\(",
    "^static int np_conditional_xrow_from_ctx\\("
  )
  compact <- gsub("[[:space:]]+", " ", body)

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
    "np_lp_delete_denominator(row_out[eval_idx], &den)",
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
  expect_false(grepl("MPI_", body, fixed = TRUE))
})

test_that("MPI adaptive acceleration has no dense shadow oracle", {
  src_file <- locate_adaptive_xrow_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  source <- paste(readLines(src_file, warn = FALSE), collapse = "\n")

  expect_false(grepl(
    "np_shadow_conditional_build_x_weights",
    source,
    fixed = TRUE
  ))
  expect_false(grepl(
    "np_shadow_cv_con_density",
    source,
    fixed = TRUE
  ))
})

test_that("MPI adaptive density CVML reuses rank-local row contexts", {
  src_file <- locate_adaptive_xrow_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  body <- adaptive_xrow_source_body(
    lines,
    "^static int np_conditional_density_cvml_lp_prepared_parallel_stream\\(",
    "^int np_conditional_density_cvml_lp_stream\\("
  )
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(
    compact,
    "const int use_xrow_ctx = (BANDWIDTH_den_extern == BW_GEN_NN) || (BANDWIDTH_den_extern == BW_ADAP_NN);",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_conditional_xrow_ctx_prepare(vector_scale_factor, &xctx)",
    fixed = TRUE
  )
  expect_match(
    compact,
    "np_conditional_yrow_ctx_prepare(vector_scale_factor, OP_NORMAL, &yctx)",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_conditional_xrow_from_ctx(&xctx, i, xrow)",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_conditional_yrow_from_ctx(&yctx, i, yrow)",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_objective_outer_buffer_finish(",
    fixed = TRUE
  )
  expect_false(grepl(
    "!int_conditional_prepared_context_extern",
    body,
    fixed = TRUE
  ))
})
