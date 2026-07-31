library(np)

locate_xrow_weighted_blas_source <- function() {
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

xrow_weighted_blas_source_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_length(start, 1L)
  expect_gte(length(stop), 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("conditional CVLS weighted BLAS gate is narrow and memory bounded", {
  src_file <- locate_xrow_weighted_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  all_source <- paste(lines, collapse = "\n")
  gate <- xrow_weighted_blas_source_body(
    lines,
    "^static int np_conditional_x_weighted_blas_profitable\\(",
    "^static void np_blas_dgemv_t_int\\("
  )

  expect_match(
    all_source,
    "#define NP_CONDITIONAL_X_WEIGHTED_BLAS_MIN_ROWS 2048",
    fixed = TRUE
  )
  expect_match(
    all_source,
    "#define NP_CONDITIONAL_X_WEIGHTED_BLAS_MAX_BYTES",
    fixed = TRUE
  )
  expect_match(
    all_source,
    "((size_t)64*(size_t)1024*(size_t)1024)",
    fixed = TRUE
  )
  expect_match(gate, "#if NP_ACCEL_GAUSS_COMPILED", fixed = TRUE)
  expect_match(
    gate,
    "nrows < NP_CONDITIONAL_X_WEIGHTED_BLAS_MIN_ROWS",
    fixed = TRUE
  )
  expect_match(gate, "nterms < 2", fixed = TRUE)
  expect_match(gate, "basis_stride < nrows", fixed = TRUE)
  expect_match(gate, "SIZE_MAX/(size_t)nterms", fixed = TRUE)
  expect_match(
    gate,
    "NP_CONDITIONAL_X_WEIGHTED_BLAS_MAX_BYTES/sizeof(double)",
    fixed = TRUE
  )
  expect_match(gate, "(void)nrows;", fixed = TRUE)
  expect_match(gate, "(void)nterms;", fixed = TRUE)
  expect_match(gate, "(void)basis_stride;", fixed = TRUE)
  expect_false(grepl("num_train\\*num_train", gate))
  expect_false(grepl("nrows\\*nrows", gate))
})

test_that("fixed CVLS weighted BLAS preserves signed row algebra and fallback", {
  src_file <- locate_xrow_weighted_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  body <- xrow_weighted_blas_source_body(
    lines,
    "^static int (NP_NOINLINE )?np_conditional_x_weight_block_pair_stream_core\\(",
    "^static int (NP_NOINLINE )?(NP_HOT_ALIGN )?np_conditional_x_weight_block_pair_gnn_stream_core\\("
  )
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(body, "BANDWIDTH_den_extern != BW_FIXED", fixed = TRUE)
  expect_match(body, "NP_NOINLINE", fixed = TRUE)
  expect_false(grepl("BW_GEN_NN", body, fixed = TRUE))
  expect_match(
    body,
    "np_conditional_x_weighted_blas_profitable(",
    fixed = TRUE
  )
  expect_match(
    compact,
    "weighted_design = (double *)malloc(weighted_count*sizeof(double));",
    fixed = TRUE
  )
  expect_match(
    body,
    "use_weighted_blas = (weighted_design != NULL);",
    fixed = TRUE
  )
  expect_match(
    body,
    "weighted_row[j] = basis_row[j]*kw[j];",
    fixed = TRUE
  )
  expect_match(body, "F77_CALL(dgemm)", fixed = TRUE)
  expect_match(body, "F77_CALL(dgemv)", fixed = TRUE)
  expect_match(
    body,
    "full_rows_out[i][orig_j] = kw[j]*mean_row[j];",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_lp_delete_denominator(full_rows_out[i][eval_idx], &den)",
    fixed = TRUE
  )
  expect_match(body, "if(weighted_design != NULL) free(weighted_design);", fixed = TRUE)
  expect_false(grepl("alloc_vecd\\(num_train\\)", body))
  expect_false(grepl("num_train\\*num_train", body))

  markers <- c(
    "np_conditional_kernel_row_raw",
    "weighted_row[j] = basis_row[j]*kw[j]",
    "F77_CALL(dgemm)",
    "np_lp_full_row_workspace_solve",
    "F77_CALL(dgemv)",
    "full_rows_out[i][orig_j] = kw[j]*mean_row[j]",
    "np_lp_delete_denominator(full_rows_out[i][eval_idx], &den)"
  )
  positions <- vapply(markers, function(marker) {
    regexpr(marker, body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))

  expect_match(
    compact,
    "if(use_weighted_blas){ const char trans_t",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(l = 0; l < k; l++) for(j = 0; j < k; j++) full_row_workspace.gram[l + j*k] = 0.0;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(j = 0; j < num_train; j++){ double zju = 0.0;",
    fixed = TRUE
  )
})

test_that("generalized-NN CVLS reuses signed weighted BLAS with scalar fallback", {
  src_file <- locate_xrow_weighted_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  body <- xrow_weighted_blas_source_body(
    lines,
    "^static int (NP_NOINLINE )?(NP_HOT_ALIGN )?np_conditional_x_weight_block_pair_gnn_stream_core\\(",
    "^static int np_conditional_x_weight_block_pair_selected_stream_core\\("
  )
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(body, "BANDWIDTH_den_extern != BW_GEN_NN", fixed = TRUE)
  expect_match(body, "NP_NOINLINE", fixed = TRUE)
  expect_false(grepl("BANDWIDTH_den_extern != BW_FIXED", body, fixed = TRUE))
  expect_match(
    body,
    "np_conditional_x_weighted_blas_profitable(",
    fixed = TRUE
  )
  expect_match(
    compact,
    "weighted_design = (double *)malloc(weighted_count*sizeof(double));",
    fixed = TRUE
  )
  expect_match(
    body,
    "use_weighted_blas = (weighted_design != NULL);",
    fixed = TRUE
  )
  expect_match(
    body,
    "weighted_row[j] = basis_row[j]*kw[j];",
    fixed = TRUE
  )
  expect_match(body, "F77_CALL(dgemm)", fixed = TRUE)
  expect_match(body, "F77_CALL(dgemv)", fixed = TRUE)
  expect_match(
    body,
    "full_rows_out[i][j] = kw[j]*mean_row[j];",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_lp_delete_denominator(full_rows_out[i][eval_idx], &den)",
    fixed = TRUE
  )
  expect_match(body, "if(weighted_design != NULL) free(weighted_design);", fixed = TRUE)
  expect_false(grepl("alloc_vecd\\(num_train\\)", body))
  expect_false(grepl("num_train\\*num_train", body))

  markers <- c(
    "np_conditional_kernel_row_raw",
    "weighted_row[j] = basis_row[j]*kw[j]",
    "F77_CALL(dgemm)",
    "np_lp_full_row_workspace_solve",
    "F77_CALL(dgemv)",
    "full_rows_out[i][j] = kw[j]*mean_row[j]",
    "np_lp_delete_denominator(full_rows_out[i][eval_idx], &den)"
  )
  positions <- vapply(markers, function(marker) {
    regexpr(marker, body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))

  expect_match(
    compact,
    "if(use_weighted_blas){ const char trans_t",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(l = 0; l < k; l++) for(j = 0; j < k; j++) full_row_workspace.gram[l + j*k] = 0.0;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(j = 0; j < num_train; j++){ double zju = 0.0;",
    fixed = TRUE
  )
})

test_that("shared conditional full-row blocks reuse bounded weighted BLAS", {
  src_file <- locate_xrow_weighted_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  body <- xrow_weighted_blas_source_body(
    lines,
    "^static int np_conditional_x_weight_block_stream_core_impl\\(",
    "^static int np_conditional_x_weight_block_stream_core\\("
  )
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(
    body,
    "np_conditional_x_weighted_blas_profitable(",
    fixed = TRUE
  )
  expect_match(
    compact,
    "weighted_design = (double *)malloc(weighted_count*sizeof(double));",
    fixed = TRUE
  )
  expect_match(body, "use_weighted_blas = (weighted_design != NULL);", fixed = TRUE)
  expect_match(body, "weighted_row[j] = basis_row[j]*kw[j];", fixed = TRUE)
  expect_match(body, "F77_CALL(dgemm)", fixed = TRUE)
  expect_match(body, "F77_CALL(dgemv)", fixed = TRUE)
  expect_match(body, "rows_out[i][orig_j] = kw[j]*mean_row[j];", fixed = TRUE)
  expect_match(body, "if(weighted_design != NULL) free(weighted_design);", fixed = TRUE)
  expect_false(grepl("num_train\\*num_train", body))

  expect_match(
    compact,
    "if(use_weighted_blas){ const char trans_t",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(l = 0; l < k; l++) for(j = 0; j < k; j++) full_row_workspace.gram[l + j*k] = 0.0;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(j = 0; j < num_train; j++){ double zju = 0.0;",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_lp_delete_denominator(rows_out[i][eval_idx], &den)",
    fixed = TRUE
  )
})

test_that("weighted BLAS remains isolated from adjacent conditional routes", {
  src_file <- locate_xrow_weighted_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  density_body <- xrow_weighted_blas_source_body(
    lines,
    "^int np_conditional_density_cvls_lp_stream\\(",
    "^static int np_conditional_distribution_cvls_lp_row_stream\\("
  )
  distribution_body <- xrow_weighted_blas_source_body(
    lines,
    "^int np_conditional_distribution_cvls_lp_stream\\(",
    "^int np_kernel_estimate_density_categorical_leave_one_out_cv\\("
  )

  expect_false(grepl("weighted_design", density_body, fixed = TRUE))
  expect_false(grepl("weighted_design", distribution_body, fixed = TRUE))
  expect_false(grepl(
    "np_conditional_x_weighted_blas_profitable",
    distribution_body,
    fixed = TRUE
  ))
})
