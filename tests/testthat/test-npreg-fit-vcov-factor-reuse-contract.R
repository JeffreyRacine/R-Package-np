locate_lp_source <- function(name) {
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
  if (length(hits) == 0L) NULL else hits[[1L]]
}

test_that("response and adjoint LP owners retain validated factorizations", {
  jksum_file <- locate_lp_source("jksum.c")
  solve_file <- locate_lp_source("jksum_lp_solve.c")
  header_file <- locate_lp_source("jksum_lp_solve.h")
  skip_if(
    any(vapply(
      list(jksum_file, solve_file, header_file),
      is.null,
      logical(1L)
    )),
    "LP source files unavailable in this test context"
  )

  jksum_lines <- readLines(jksum_file, warn = FALSE)
  solve_lines <- readLines(solve_file, warn = FALSE)
  header_lines <- readLines(header_file, warn = FALSE)

  expect_equal(
    sum(grepl("np_lp_solve_workspace_solve_factored\\(", jksum_lines)),
    1L
  )
  # Classify retained back-solves by owner rather than pinning a whole-file
  # count. PREPARING and stable categorical contrasts also need adjoints;
  # neither adds a refactorization to accepted covariance directions.
  adjoint.name <- "np_lp_solve_workspace_solve_adjoint_factored("
  owner.names <- c("np_conditional_lp_project_accepted",
                   "np_regression_general_lp_point_at_frame",
                   "np_regression_general_lp_fit_execute")
  owner.bodies <- lapply(owner.names, function(name)
    np_test_extract_c_function(jksum_lines, name))
  names(owner.bodies) <- owner.names
  count.calls <- function(body, name) {
    hit <- gregexpr(name, body, fixed = TRUE)[[1L]]
    sum(hit > 0L)
  }
  # The serial source retains an MPI-guarded branch, but only its serial
  # owner has the categorical SE-off base-adjoint call.
  expect_identical(vapply(owner.bodies, count.calls, integer(1L),
                          name = adjoint.name),
                   setNames(c(1L, 1L, 5L), owner.names))
  conditional <- owner.bodies[[1L]]
  frame <- owner.bodies[[2L]]
  fit <- owner.bodies[[3L]]
  expect_match(conditional, "if(diagnostics->ridge_total > 0.0)", fixed = TRUE)
  expect_match(conditional, "np_lp_solve_workspace_solve_factored(", fixed = TRUE)
  expect_match(frame, "if(projection != NULL)", fixed = TRUE)
  expect_match(fit, "if(hc0_residual_preparing)", fixed = TRUE)
  expect_identical(count.calls(fit, "np_regression_hc0_lp_prepare_row("), 2L)
  expect_match(fit, "if(!ordinary_hc0 && !hc0_residual_preparing)", fixed = TRUE)
  expect_identical(count.calls(fit, "if(!ordinary_hc0 && !hc0_residual_preparing)"),
                   1L)
  expect_match(gsub("[[:space:]]+", " ", fit),
    paste0(adjoint.name, " &owner->solve_workspace, owner->nterms, variance_rhs,"),
    fixed = TRUE)
  expect_true(any(grepl("int factor_ready;", header_lines, fixed = TRUE)))
  expect_true(any(grepl("int factor_p;", header_lines, fixed = TRUE)))

  response_body <- np_test_extract_c_function(
    solve_lines, "np_lp_solve_workspace_try_dgesv"
  )
  expect_true(grepl("workspace->factor_ready = 0;", response_body, fixed = TRUE))
  expect_true(grepl("workspace->factor_ready = 1;", response_body, fixed = TRUE))
  expect_true(grepl("workspace->factor_p = p;", response_body, fixed = TRUE))
  expect_true(grepl("F77_CALL(dgesv)", response_body, fixed = TRUE))

  response_policy <- np_test_extract_c_function(
    solve_lines, "np_lp_solve_workspace_solve_response_ranked"
  )
  adjoint_policy <- np_test_extract_c_function(
    solve_lines, "np_lp_solve_workspace_solve_adjoint_ranked"
  )
  for (body in list(response_policy, adjoint_policy)) {
    expect_true(grepl(
      "np_lp_solve_workspace_admit_retained_factor(", body, fixed = TRUE
    ))
  }
  expect_true(grepl(
    "np_lp_solve_workspace_solve_adjoint_factored(",
    adjoint_policy,
    fixed = TRUE
  ))

  factored_body <- np_test_extract_c_function(
    solve_lines, "np_lp_solve_workspace_solve_factored_with_trans"
  )
  expect_true(grepl("!workspace->factor_ready", factored_body, fixed = TRUE))
  expect_true(grepl("workspace->factor_p != p", factored_body, fixed = TRUE))
  expect_true(grepl("F77_CALL(dgetrs)", factored_body, fixed = TRUE))
  expect_false(grepl("F77_CALL(dgesv)", factored_body, fixed = TRUE))
  expect_false(grepl("F77_CALL(dgetrf)", factored_body, fixed = TRUE))
  expect_false(grepl("gram_source", factored_body, fixed = TRUE))
})
