library(npRmpi)

locate_jksum_c <- function() {
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
  if (length(hits) == 0L) {
    return(NULL)
  }
  hits[[1L]]
}

locate_kernelcv_c <- function() {
  candidates <- c(
    test_path("..", "..", "src", "kernelcv.c"),
    test_path("..", "..", "..", "src", "kernelcv.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "kernelcv.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "kernelcv.c"),
    file.path(getwd(), "src", "kernelcv.c"),
    file.path(getwd(), "..", "src", "kernelcv.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) {
    return(NULL)
  }
  hits[[1L]]
}

test_that("canonical LP CV route predicates remain centralized", {
  src_file <- locate_jksum_c()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")

  lines <- readLines(src_file, warn = FALSE)

  expect_true(any(grepl("np_reg_cv_use_symmetric_dropone_path", lines, fixed = TRUE)))
  expect_true(any(grepl("np_reg_cv_use_canonical_lp_fixed_kernel", lines, fixed = TRUE)))
  expect_equal(sum(grepl("np_reg_cv_use_degree1_rawbasis_kernel", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("np_regression_cv_degree1_rawbasis_fixed", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("np_reg_degree1_center_raw_moments_at_eval", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("np_glp_center_raw_moments_at_eval", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("np_glp_fill_shift_raw_from_center", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("np_glp_binom_coeff", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("np_reg_use_canonical_glp_degree1_estimation", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("use_ll_compatible_gnn_se", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("raw_degree1_glp", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("np_glp_fill_basis_raw_train_centered", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("np_glp_fill_basis_eval_raw_centered", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("np_glp_fill_basis_eval_deriv_raw_centered", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("LL_LL", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("const int int_ll_est = int_ll;", lines, fixed = TRUE)), 1L)
  expect_equal(sum(grepl("if(kpow == 2)", lines, fixed = TRUE)), 1L)
  expect_equal(sum(grepl("wbuf[k] = (weights[k] == 0.0) ? 0.0 : wk*wk;", lines, fixed = TRUE)), 1L)

  helper_start <- grep("^static inline int np_reg_cv_use_symmetric_dropone_path\\(", lines)
  helper_stop <- grep("^static int np_lp_fixed_tree_sparse_supported\\(", lines)
  expect_length(helper_start, 1L)
  expect_length(helper_stop, 1L)
  expect_lt(helper_start, helper_stop)
  helper_body <- paste(lines[helper_start:(helper_stop - 1L)], collapse = "\n")
  expect_true(grepl("bwm == RBWM_CVLS", helper_body, fixed = TRUE))
  expect_true(grepl("bwm == RBWM_CVCHECK", helper_body, fixed = TRUE))
  expect_true(grepl("ks_tree_use", helper_body, fixed = TRUE))
  expect_true(grepl("BANDWIDTH_reg == BW_ADAP_NN", helper_body, fixed = TRUE))

  helper_calls <- sum(grepl("np_reg_cv_use_symmetric_dropone_path\\(", lines))
  expect_gte(helper_calls, 2L)

  canonical_lp_calls <- sum(grepl("np_reg_cv_use_canonical_lp_fixed_kernel\\(", lines))
  expect_gte(canonical_lp_calls, 2L)

  canonical_start <- grep(
    "^static inline int np_reg_cv_use_canonical_lp_fixed_kernel\\(",
    lines
  )
  canonical_stop <- grep(
    "^static int np_lp_fixed_tree_sparse_supported\\(",
    lines
  )
  expect_length(canonical_start, 1L)
  expect_length(canonical_stop, 1L)
  expect_lt(canonical_start, canonical_stop)

  canonical_body <- paste(
    lines[canonical_start:(canonical_stop - 1L)],
    collapse = "\n"
  )
  expect_true(grepl("const int nterms", canonical_body, fixed = TRUE))
  expect_true(any(grepl(
    "NP_REG_CV_LP_RESIDENT_MAX_TERMS = 5",
    lines,
    fixed = TRUE
  )))
  expect_true(grepl(
    "return nterms <= NP_REG_CV_LP_RESIDENT_MAX_TERMS;",
    canonical_body,
    fixed = TRUE
  ))
  expect_true(grepl("return bwm == RBWM_CVAIC;", canonical_body, fixed = TRUE))
  expect_false(grepl("use_bernstein", canonical_body, fixed = TRUE))
  expect_false(any(grepl(
    "np_regression_cv_lp_rawbasis_fixed",
    lines,
    fixed = TRUE
  )))
})

test_that("legacy LL compute engines and restoration switches are absent", {
  src_file <- locate_jksum_c()
  skip_if(is.null(src_file), "package C sources unavailable in this test context")

  src_dir <- dirname(src_file)
  source_files <- file.path(src_dir, c(
    "headers.h", "jksum.c", "kernelcv.c", "kernele.c", "np.c"
  ))
  skip_if_not(all(file.exists(source_files)), "complete package C sources unavailable")

  source <- paste(unlist(lapply(source_files, readLines, warn = FALSE)),
                  collapse = "\n")
  forbidden <- c(
    "LL_LL",
    "REGTYPE_LL",
    "REG_OLDREGI",
    "old_reg",
    "kernel_estimate_regression_categorical_leave_one_out",
    "kernel_estimate_regression_categorical_no_stderr",
    "kernel_estimate_regression_categorical_aic_c",
    "kernel_estimate_categorical_gradient_ocg_fast"
  )
  for (token in forbidden)
    expect_false(grepl(token, source, fixed = TRUE), info = token)

  expect_false(grepl(
    "int_ll_extern == LL_LP) ? LL_LP : LL_LC",
    source,
    fixed = TRUE
  ))
})

test_that("fixed resident-row LP CV uses the reusable uncentered solve workspace", {
  src_file <- locate_jksum_c()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")

  lines <- readLines(src_file, warn = FALSE)
  helper_start <- grep(
    "^static NPRegCvLpResult np_regression_cv_lp_basis_fixed\\(",
    lines
  )
  helper_stop <- grep(
    "^double np_kernel_estimate_regression_categorical_ls_aic\\(",
    lines
  )
  expect_length(helper_start, 1L)
  expect_length(helper_stop, 1L)
  expect_lt(helper_start, helper_stop)

  helper_body <- paste(lines[helper_start:(helper_stop - 1L)], collapse = "\n")
  expect_true(grepl(
    "solve_workspace.gram_source[a + b*nterms] = sj[a*nterms+b];",
    helper_body,
    fixed = TRUE
  ))
  expect_true(grepl(
    "fit += eval_basis[a]*solve_workspace.rhs_work[a];",
    helper_body,
    fixed = TRUE
  ))
  expect_true(grepl(
    "hii += eval_basis[a]*solve_workspace.rhs_work[a];",
    helper_body,
    fixed = TRUE
  ))
  expect_true(grepl(
    "np_lp_solve_workspace_solve(&solve_workspace, nterms, 1)",
    helper_body,
    fixed = TRUE
  ))
  expect_false(grepl("mat_solve(", helper_body, fixed = TRUE))
  expect_false(grepl("fit = DELTA[0][0];", helper_body, fixed = TRUE))
  expect_false(grepl("center_raw", helper_body, fixed = TRUE))
  expect_false(grepl("SHIFT", helper_body, fixed = TRUE))
  expect_false(grepl("mat_inv00", helper_body, fixed = TRUE))
})

test_that("packed and nearest-neighbor LP CV avoid legacy solve marshalling", {
  src_file <- locate_jksum_c()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")

  lines <- readLines(src_file, warn = FALSE)
  helper_start <- grep(
    "^double np_kernel_estimate_regression_categorical_ls_aic\\(",
    lines
  )
  helper_stop <- grep(
    "^static int np_distribution_cvls_ordered_profile_stream\\(",
    lines
  )
  expect_length(helper_start, 1L)
  expect_length(helper_stop, 1L)
  expect_lt(helper_start, helper_stop)

  helper_body <- paste(lines[helper_start:(helper_stop - 1L)], collapse = "\n")
  expect_true(grepl(
    "np_lp_solve_workspace_reserve(&solve_workspace, nrc1, 1)",
    helper_body,
    fixed = TRUE
  ))
  expect_true(grepl(
    "solve_workspace.gram_source[i+k*nrc1]",
    helper_body,
    fixed = TRUE
  ))
  expect_true(grepl(
    "mhat += evalv[i]*solve_workspace.rhs_work[i];",
    helper_body,
    fixed = TRUE
  ))
  expect_true(grepl(
    "hii += evalv[i]*solve_workspace.rhs_work[i];",
    helper_body,
    fixed = TRUE
  ))
  expect_false(grepl("mat_solve(", helper_body, fixed = TRUE))
  expect_false(grepl("MATRIX XTKY", helper_body, fixed = TRUE))
  expect_false(grepl("DELTA", helper_body, fixed = TRUE))
  expect_false(grepl("MATRIX KWM", helper_body, fixed = TRUE))
})

test_that("density CV tree-bypass predicate is centralized in one helper", {
  src_file <- locate_jksum_c()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")

  lines <- readLines(src_file, warn = FALSE)
  expect_true(any(grepl("np_den_cv_use_tree_bypass_path", lines, fixed = TRUE)))

  raw_pred <- "gate_x_all_large_fixed || !int_TREE_XY || (BANDWIDTH_den == BW_ADAP_NN)"
  # One raw predicate instance is allowed (the helper definition itself).
  expect_equal(sum(grepl(raw_pred, lines, fixed = TRUE)), 1L)
  # No inline ad hoc usage in if-statements.
  expect_equal(sum(grepl("if\\(gate_x_all_large_fixed \\|\\| !int_TREE_XY \\|\\| \\(BANDWIDTH_den == BW_ADAP_NN\\)\\)", lines)), 0L)
})

test_that("conditional CV large-kernel gating keeps Y and XY activation routes live", {
  src_file <- locate_jksum_c()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")

  lines <- readLines(src_file, warn = FALSE)

  expect_false(any(grepl("Canonical conditional CV gate policy: large-kernel shortcuts are X-only", lines, fixed = TRUE)))
  expect_gte(sum(grepl("^\\s*gate_y_active = 1;\\s*$", lines)), 2L)
  expect_gte(sum(grepl("^\\s*gate_xy_active = 1;\\s*$", lines)), 1L)
  expect_equal(sum(grepl("^\\s*gate_y_active = 0;\\s*$", lines)), 0L)
  expect_equal(sum(grepl("^\\s*gate_xy_active = 0;\\s*$", lines)), 0L)
})

test_that("conditional public LP CV routes stay thin kernelcv dispatches", {
  src_file <- locate_kernelcv_c()
  skip_if(is.null(src_file), "source file src/kernelcv.c unavailable in this test context")

  lines <- readLines(src_file, warn = FALSE)

  expect_equal(sum(grepl("^double cv_func_con_density_categorical_ml\\(", lines)), 1L)
  expect_equal(sum(grepl("^double cv_func_con_density_categorical_ls\\(", lines)), 1L)
  expect_equal(sum(grepl("^double cv_func_con_distribution_categorical_ls\\(", lines)), 1L)
  expect_equal(sum(grepl("^double np_cv_func_con_density_categorical_ml\\(", lines)), 1L)
  expect_equal(sum(grepl("^double np_cv_func_con_density_categorical_ls\\(", lines)), 1L)
  expect_false(any(grepl("np_shadow_proof_cv_con_density_ml\\(", lines)))
  expect_false(any(grepl("np_shadow_proof_cv_con_density_ls\\(", lines)))
  expect_false(any(grepl("np_shadow_proof_cv_con_distribution_ls\\(", lines)))
  expect_equal(sum(grepl("kernel_estimate_regression_categorical_tree_np", lines, fixed = TRUE)), 0L)
})
