library(npRmpi)

locate_npreg_adaptive_blas_source <- function() {
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

npreg_adaptive_source_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_length(start, 1L)
  expect_gte(length(stop), 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("adaptive regression BLAS eligibility is narrow and bounded", {
  src_file <- locate_npreg_adaptive_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  body <- npreg_adaptive_source_body(
    lines,
    "^static .*NPRegCvLpResult np_regression_cv_lp_basis_adaptive_blas\\(",
    "^double np_kernel_estimate_regression_categorical_ls_aic\\("
  )
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(
    compact,
    "num_obs < NP_CONDITIONAL_X_WEIGHTED_BLAS_MIN_ROWS",
    fixed = TRUE
  )
  expect_match(compact, "nterms < 4", fixed = TRUE)
  expect_match(
    body,
    "np_mseries_accelerate_enabled_cache",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_conditional_x_weighted_blas_profitable(",
    fixed = TRUE
  )
  expect_match(
    compact,
    "(kernel_c[l] != 0) && (kernel_c[l] != 4) && (kernel_c[l] != 8)",
    fixed = TRUE
  )
  expect_false(grepl("num_obs\\s*\\*\\s*num_obs", body))
  expect_false(grepl("num_obs\\)\\s*\\*\\s*\\(size_t\\)num_obs", body))
})

test_that("adaptive regression BLAS preserves delete-one and full-row algebra", {
  src_file <- locate_npreg_adaptive_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  body <- npreg_adaptive_source_body(
    lines,
    "^static .*NPRegCvLpResult np_regression_cv_lp_basis_adaptive_blas\\(",
    "^double np_kernel_estimate_regression_categorical_ls_aic\\("
  )
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(body, "self_weight = kw[j];", fixed = TRUE)
  expect_match(body, "kw[j] = 0.0;", fixed = TRUE)
  expect_match(body, "F77_CALL(dgemm)", fixed = TRUE)
  expect_match(body, "np_blas_dgemv_t_int(", fixed = TRUE)
  expect_match(
    compact,
    "solve_workspace.rhs_source[l] += self_weight*bl*yj;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "self_weight*bl*eval_basis[i]",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_regression_cv_loss_value(bwm, fit,",
    fixed = TRUE
  )
  expect_match(
    body,
    "trace_rows[j] = self_weight*hii;",
    fixed = TRUE
  )
})

test_that("adaptive regression BLAS MPI ownership and fallback are symmetric", {
  src_file <- locate_npreg_adaptive_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  all_source <- paste(lines, collapse = "\n")
  allocator <- npreg_adaptive_source_body(
    lines,
    "^static double \\*\\*np_regression_adaptive_one_row_matrix\\(",
    "^static .*NPRegCvLpResult np_regression_cv_lp_basis_adaptive_blas\\("
  )
  body <- npreg_adaptive_source_body(
    lines,
    "^static .*NPRegCvLpResult np_regression_cv_lp_basis_adaptive_blas\\(",
    "^double np_kernel_estimate_regression_categorical_ls_aic\\("
  )
  owner <- npreg_adaptive_source_body(
    lines,
    "^double np_kernel_estimate_regression_categorical_ls_aic\\(",
    "^static int np_distribution_cvls_ordered_profile_stream\\("
  )
  owner_compact <- gsub("[[:space:]]+", " ", owner)

  expect_match(allocator, "if(matrix == NULL)", fixed = TRUE)
  expect_match(allocator, "if(storage == NULL)", fixed = TRUE)
  expect_match(
    body,
    "for(j = my_rank; j < num_obs; j += iNum_Processors)",
    fixed = TRUE
  )
  expect_match(
    body,
    "cv_rows = (double *)calloc((size_t)num_obs, sizeof(double));",
    fixed = TRUE
  )
  expect_match(
    body,
    "NP_RMPI_INJECT_REG_ADAPTIVE_BLAS_FAIL_RANK",
    fixed = TRUE
  )
  expect_match(body, "MPI_Allreduce(&local_ready,", fixed = TRUE)
  expect_match(body, "MPI_Allreduce(&local_fail,", fixed = TRUE)
  expect_match(body, "MPI_Allreduce(MPI_IN_PLACE,", fixed = TRUE)
  expect_match(
    body,
    "for(j = 0; j < num_obs; j++){",
    fixed = TRUE
  )
  expect_match(
    all_source,
    "if(!result.ok){\n    result.cv = DBL_MAX;",
    fixed = TRUE
  )
  expect_match(owner, "if(adaptive_result.ok){", fixed = TRUE)
  expect_match(
    owner_compact,
    "(BANDWIDTH_reg == BW_ADAP_NN) && (glp_nterms >= 4)",
    fixed = TRUE
  )
  expect_match(
    owner,
    "const size_t kwm_len =",
    fixed = TRUE
  )
})
