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

  expect_match(
    paste(lines, collapse = "\n"),
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
    "basis = basis_is_contiguous ? alloc_tmatd(num_obs_train, glp_nterms) : alloc_matd(num_obs_train, glp_nterms);",
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
    "if(fit_progress_active) np_progress_fit_loop_step(i + row + 1, fit_progress_total);",
    fixed = TRUE
  )
  expect_false(grepl("num_obs_eval\\s*\\*\\s*num_obs_eval", body))
  expect_false(grepl("num_obs_train\\s*\\*\\s*num_obs_train", body))
})
