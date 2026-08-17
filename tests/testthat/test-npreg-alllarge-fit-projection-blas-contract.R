library(npRmpi)

locate_alllarge_fit_projection_source <- function() {
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

test_that("MPI all-large fit BLAS route is narrow and bounded", {
  src_file <- locate_alllarge_fit_projection_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  start <- grep(
    "^int kernel_estimate_regression_categorical_tree_np\\(",
    lines
  )
  stop <- grep("^  if\\(estimation_shortcut_done\\)", lines)
  expect_length(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  body <- paste(lines[start:stop], collapse = "\n")
  compact <- gsub("[[:space:]]+", " ", body)
  source <- paste(lines, collapse = "\n")
  source_compact <- gsub("[[:space:]]+", " ", source)

  expect_match(
    source,
    "#define NP_REG_ALLLARGE_FIT_BLOCK_MAX_ROWS 8192",
    fixed = TRUE
  )
  expect_match(compact, "basis_is_contiguous = (!do_grad) &&", fixed = TRUE)
  expect_match(
    compact,
    "np_reg_alllarge_blas_profitable(num_obs_train, glp_nterms, num_obs_train)",
    fixed = TRUE
  )
  expect_match(
    compact,
    "np_reg_alllarge_blas_profitable(num_obs_eval, glp_nterms, num_obs_eval)",
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "basis = np_regression_fit_matrix_alloc( ",
      "num_obs_train, glp_nterms, basis_is_contiguous);"
    ),
    fixed = TRUE
  )
  expect_match(body, "np_blas_gram_int(", fixed = TRUE)
  expect_match(body, "np_blas_dgemv_t_int(", fixed = TRUE)
  expect_match(body, "np_blas_dgemv_n_int(", fixed = TRUE)
  expect_match(body, "np_blas_project_inverse_block_int(", fixed = TRUE)
  expect_match(
    compact,
    "use_fit_projection_blas = basis_is_contiguous;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(i = 0; i < glp_nterms; i++){ inverse_workspace.rhs[i] = 0.0;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "if(basis_is_contiguous) free_tmat(basis); else free_mat(basis, glp_nterms);",
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "R_UnwindProtect( ",
      "np_glp_fit_basis_prepare_execute, ",
      "(void *)&basis_call, ",
      "np_regression_alllarge_lp_fit_cleanup, ",
      "(void *)&all_large_owner, NULL);"
    ),
    fixed = TRUE
  )
  expect_match(
    source_compact,
    "if(!jump) return; if(owner->basis_context != NULL){",
    fixed = TRUE
  )
  expect_match(
    source_compact,
    paste0(
      "free(owner->eval_basis_block); free(owner->projection_block); ",
      "free(owner->fitted_block);"
    ),
    fixed = TRUE
  )
  expect_equal(
    sum(gregexpr(
      "np_progress_fit_loop_step_owned(", body, fixed = TRUE
    )[[1L]] > 0L),
    3L
  )
  expect_match(
    source_compact,
    paste0(
      "free(owner->eval_basis); free(owner->eval_derivative); ",
      "free(owner->coefficient); ",
      "np_lp_full_row_workspace_clear(owner->inverse_workspace);"
    ),
    fixed = TRUE
  )
  expect_match(
    source_compact,
    paste0(
      "static double **np_regression_fit_matrix_alloc( ",
      "const int nrows, const int ncols, const int contiguous)"
    ),
    fixed = TRUE
  )
  expect_match(
    source_compact,
    paste0(
      "if(matrix[column] == NULL) { while(column > 0) ",
      "free(matrix[--column]); free(matrix); return NULL;"
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      ".basis_context = &all_large_owner.basis_context, ",
      ".use_stable_basis = &use_stable_basis, .status = 0"
    ),
    fixed = TRUE
  )
  refresh_position <- regexpr(
    "np_refresh_runtime_tolerances();", compact, fixed = TRUE
  )[[1L]]
  owner_position <- regexpr(
    "np_regression_fit_owner_init(", compact, fixed = TRUE
  )[[1L]]
  terms_position <- regexpr(
    "fast_ok = np_glp_build_terms(", compact, fixed = TRUE
  )[[1L]]
  expect_gt(refresh_position, 0L)
  expect_gt(owner_position, refresh_position)
  expect_gt(terms_position, owner_position)
  expect_false(grepl(
    "np_refresh_mseries_accelerate_option();", body, fixed = TRUE
  ))
  expect_match(
    compact,
    paste0(
      "if(fit_progress_active) np_progress_fit_loop_step_owned( ",
      "i + row + 1, fit_progress_total, ",
      "np_regression_alllarge_lp_fit_cleanup, ",
      "(void *)&all_large_owner);"
    ),
    fixed = TRUE
  )
  expect_false(grepl("num_obs_eval\\s*\\*\\s*num_obs_eval", body))
  expect_false(grepl("num_obs_train\\s*\\*\\s*num_obs_train", body))
})
