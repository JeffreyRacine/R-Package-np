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

test_that("CVLS LL/LP route predicate is centralized in one helper", {
  src_file <- locate_jksum_c()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")

  lines <- readLines(src_file, warn = FALSE)

  expect_true(any(grepl("np_reg_cv_use_symmetric_dropone_path", lines, fixed = TRUE)))
  expect_true(any(grepl("np_reg_cv_use_canonical_glp_fixed_kernel", lines, fixed = TRUE)))
  expect_equal(sum(grepl("np_reg_cv_use_degree1_rawbasis_kernel", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("np_regression_cv_degree1_rawbasis_fixed", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("np_reg_degree1_center_raw_moments_at_eval", lines, fixed = TRUE)), 0L)

  raw_pred <- "(bwm == RBWM_CVLS) || ks_tree_use || (BANDWIDTH_reg == BW_ADAP_NN)"
  # Exactly one raw predicate instance is allowed (the helper definition itself).
  expect_equal(sum(grepl(raw_pred, lines, fixed = TRUE)), 1L)
  # The two-term suffix should also only appear inside the helper definition.
  expect_equal(sum(grepl("ks_tree_use || (BANDWIDTH_reg == BW_ADAP_NN)", lines, fixed = TRUE)), 1L)
  # No inline ad hoc route predicates should remain in if-statements.
  expect_equal(sum(grepl("if\\(\\(bwm == RBWM_CVLS\\)", lines)), 0L)
  expect_equal(sum(grepl("if\\(ks_tree_use \\|\\| \\(BANDWIDTH_reg == BW_ADAP_NN\\)\\)", lines)), 0L)

  helper_calls <- sum(grepl("np_reg_cv_use_symmetric_dropone_path\\(", lines))
  expect_gte(helper_calls, 3L)

  canonical_glp_calls <- sum(grepl("np_reg_cv_use_canonical_glp_fixed_kernel\\(", lines))
  expect_gte(canonical_glp_calls, 3L)
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
