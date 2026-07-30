locate_lp_width_one_source <- function(name) {
  candidates <- c(
    test_path("..", "..", "src", name),
    test_path("..", "..", "..", "src", name),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", name),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", name),
    file.path(getwd(), "src", name),
    file.path(getwd(), "..", "src", name)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits))
    return(NULL)
  hits[[1L]]
}

lp_width_one_region <- function(source, start_text, stop_text) {
  start <- regexpr(start_text, source, fixed = TRUE)
  stop <- regexpr(stop_text, source, fixed = TRUE)
  expect_gt(start, 0L)
  expect_gt(stop, start)
  substr(source, start, stop - 1L)
}

test_that("canonical solve workspace dispatches width one only to scalar algebra", {
  solve_file <- locate_lp_width_one_source("jksum_lp_solve.c")
  skip_if(is.null(solve_file), "source file src/jksum_lp_solve.c unavailable")
  source <- paste(readLines(solve_file, warn = FALSE), collapse = "\n")

  scalar <- lp_width_one_region(
    source,
    "static inline int np_lp_solve_width_one(",
    "int np_lp_solve_workspace_sources_finite("
  )
  solve <- lp_width_one_region(
    source,
    "int np_lp_solve_workspace_solve(",
    "int np_lp_solve_workspace_solve_factored("
  )
  factored <- lp_width_one_region(
    source,
    "int np_lp_solve_workspace_solve_factored(",
    "NPLPWidthOneStatus np_lp_width_one_influence_row("
  )

  expect_true(grepl("rhs_source[i]/gram", scalar, fixed = TRUE))
  expect_true(grepl("R_FINITE(solution)", scalar, fixed = TRUE))
  expect_false(grepl("F77_", scalar, fixed = TRUE))
  expect_false(grepl("kernel_weighted_sum", scalar, fixed = TRUE))
  expect_false(grepl("LL_LC", scalar, fixed = TRUE))

  expect_true(grepl("if(p == 1)", solve, fixed = TRUE))
  expect_true(grepl(
    "np_lp_solve_width_one(workspace->gram_source[0]",
    solve,
    fixed = TRUE
  ))
  expect_equal(lengths(regmatches(
    solve,
    gregexpr("F77_CALL(dgesv)", solve, fixed = TRUE)
  )), 1L)
  expect_lt(
    regexpr("if(p == 1)", solve, fixed = TRUE),
    regexpr("F77_CALL(dgesv)", solve, fixed = TRUE)
  )

  expect_true(grepl("if(p == 1)", factored, fixed = TRUE))
  expect_true(grepl(
    "np_lp_solve_width_one(workspace->gram_work[0]",
    factored,
    fixed = TRUE
  ))
  expect_equal(lengths(regmatches(
    factored,
    gregexpr("F77_CALL(dgetrs)", factored, fixed = TRUE)
  )), 1L)
  expect_lt(
    regexpr("if(p == 1)", factored, fixed = TRUE),
    regexpr("F77_CALL(dgetrs)", factored, fixed = TRUE)
  )
})

test_that("full-row workspace dispatches width one before every LAPACK owner", {
  solve_file <- locate_lp_width_one_source("jksum_lp_solve.c")
  skip_if(is.null(solve_file), "source file src/jksum_lp_solve.c unavailable")
  source <- paste(readLines(solve_file, warn = FALSE), collapse = "\n")

  rcond <- lp_width_one_region(
    source,
    "static int np_lp_full_row_bad_rcond(",
    "int np_lp_full_row_workspace_solve("
  )
  solve <- lp_width_one_region(
    source,
    "int np_lp_full_row_workspace_solve(",
    "static int np_lp_full_row_workspace_factor_invert("
  )
  invert <- lp_width_one_region(
    source,
    "int np_lp_full_row_workspace_invert(",
    "int np_lp_full_row_workspace_invert_retryable("
  )
  retry <- lp_width_one_region(
    source,
    "int np_lp_full_row_workspace_invert_retryable(",
    "int np_lp_full_row_workspace_pack_inverse_rows("
  )

  expect_true(grepl("if(p == 1)", rcond, fixed = TRUE))
  expect_lt(
    regexpr("if(p == 1)", rcond, fixed = TRUE),
    regexpr("F77_CALL(dsyev)", rcond, fixed = TRUE)
  )

  expect_true(grepl("if(p == 1)", solve, fixed = TRUE))
  expect_true(grepl("np_lp_solve_width_one(", solve, fixed = TRUE))
  expect_lt(
    regexpr("if(p == 1)", solve, fixed = TRUE),
    regexpr("F77_CALL(dgesv)", solve, fixed = TRUE)
  )

  expect_true(grepl("if(p == 1)", invert, fixed = TRUE))
  expect_true(grepl("1.0/workspace->gram[0]", invert, fixed = TRUE))
  expect_false(grepl("F77_", invert, fixed = TRUE))

  expect_true(grepl("if(p == 1)", retry, fixed = TRUE))
  expect_true(grepl("1.0/workspace->gram[0]", retry, fixed = TRUE))
  expect_true(grepl("} else {", retry, fixed = TRUE))
  expect_true(grepl(
    "np_lp_full_row_workspace_factor_invert(workspace, p)",
    retry,
    fixed = TRUE
  ))
})

test_that("width-one influence rows retain signed weights without LAPACK", {
  solve_file <- locate_lp_width_one_source("jksum_lp_solve.c")
  skip_if(is.null(solve_file), "source file src/jksum_lp_solve.c unavailable")
  source <- paste(readLines(solve_file, warn = FALSE), collapse = "\n")

  scalar <- lp_width_one_region(
    source,
    "NPLPWidthOneStatus np_lp_width_one_influence_row(",
    "void np_lp_full_row_workspace_init("
  )

  expect_true(grepl(
    "denominator += kw[i]*zi*zi",
    scalar,
    fixed = TRUE
  ))
  expect_true(grepl(
    "kw[i]*basis_train[i]*projection",
    scalar,
    fixed = TRUE
  ))
  expect_false(grepl("positive_weights_only", scalar, fixed = TRUE))
  expect_false(grepl("F77_", scalar, fixed = TRUE))
  expect_false(grepl("sqrt(", scalar, fixed = TRUE))
  expect_false(grepl("malloc(", scalar, fixed = TRUE))
  expect_true(grepl("NP_LP_SOLVE_MAX_RIDGE_STEPS", scalar, fixed = TRUE))
  expect_true(grepl("NP_LP_WIDTH_ONE_NONFINITE", scalar, fixed = TRUE))
  expect_true(grepl("NP_LP_WIDTH_ONE_RIDGE_FAILED", scalar, fixed = TRUE))
})

test_that("compiled width-one hats reuse the scalar influence primitive", {
  hat_file <- locate_lp_width_one_source("reghat_fast.c")
  skip_if(is.null(hat_file), "source file src/reghat_fast.c unavailable")
  source <- paste(readLines(hat_file, warn = FALSE), collapse = "\n")

  scalar_hat <- lp_width_one_region(
    source,
    "static SEXP np_reghat_width_one_matrix(",
    "static int np_reghat_solve_system("
  )
  hat_start <- regexpr("SEXP C_np_reghat_lp_matrix_fast(", source, fixed = TRUE)
  expect_gt(hat_start, 0L)
  hat <- substr(source, hat_start, nchar(source))

  expect_true(grepl("np_lp_width_one_influence_row(", scalar_hat,
                    fixed = TRUE))
  expect_true(grepl("NP_LP_WIDTH_ONE_OK", scalar_hat, fixed = TRUE))
  expect_false(grepl("R_alloc(", scalar_hat, fixed = TRUE))
  expect_false(grepl("F77_CALL(", scalar_hat, fixed = TRUE))
  expect_false(grepl("positive_weights_only", scalar_hat, fixed = TRUE))
  expect_true(grepl("if(nterms == 1)", hat, fixed = TRUE))
  expect_lt(
    regexpr("if(nterms == 1)", hat, fixed = TRUE),
    regexpr("weighted_design = (double *)R_alloc", hat, fixed = TRUE)
  )
  expect_false(grepl("if((nterms == 1)", hat, fixed = TRUE))
})
