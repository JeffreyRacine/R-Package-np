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
    2L
  )
  expect_true(any(grepl("int factor_ready;", header_lines, fixed = TRUE)))
  expect_true(any(grepl("int factor_p;", header_lines, fixed = TRUE)))

  response_body <- npRmpi_test_extract_c_function(
    solve_lines, "np_lp_solve_workspace_try_dgesv"
  )
  expect_true(grepl("workspace->factor_ready = 0;", response_body, fixed = TRUE))
  expect_true(grepl("workspace->factor_ready = 1;", response_body, fixed = TRUE))
  expect_true(grepl("workspace->factor_p = p;", response_body, fixed = TRUE))
  expect_true(grepl("F77_CALL(dgesv)", response_body, fixed = TRUE))

  response_policy <- npRmpi_test_extract_c_function(
    solve_lines, "np_lp_solve_workspace_solve_response_ranked"
  )
  adjoint_policy <- npRmpi_test_extract_c_function(
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

  factored_body <- npRmpi_test_extract_c_function(
    solve_lines, "np_lp_solve_workspace_solve_factored_with_trans"
  )
  expect_true(grepl("!workspace->factor_ready", factored_body, fixed = TRUE))
  expect_true(grepl("workspace->factor_p != p", factored_body, fixed = TRUE))
  expect_true(grepl("F77_CALL(dgetrs)", factored_body, fixed = TRUE))
  expect_false(grepl("F77_CALL(dgesv)", factored_body, fixed = TRUE))
  expect_false(grepl("F77_CALL(dgetrf)", factored_body, fixed = TRUE))
  expect_false(grepl("gram_source", factored_body, fixed = TRUE))
})
