locate_general_lp_fit_source <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  roots <- roots[nzchar(roots)]
  paths <- file.path(roots, "src", "jksum.c")
  paths <- paths[file.exists(paths)]
  if (!length(paths))
    return(NULL)
  paths[[1L]]
}

general_lp_owner_region <- function(source) {
  start <- grep("^  NP_REGRESSION_GENERAL_LP_FIT_OK = 0,$", source)
  marker <- grep("^ \\* The all-large shortcut allocates", source)
  stop <- marker - 1L
  stopifnot(length(start) == 1L, length(stop) == 1L, stop > start)
  paste(source[seq.int(start - 1L, stop - 1L)], collapse = "\n")
}

general_lp_callsite_region <- function(source) {
  start <- grep(
    "^  } else if\\(lp_engine_est == NP_LP_ENGINE_GENERAL\\)", source
  )
  stop <- grep("^finish_regression_estimation:$", source)
  stopifnot(length(start) == 1L, length(stop) == 1L, stop > start)
  paste(source[seq.int(start, stop - 1L)], collapse = "\n")
}

test_that("main general-LP fit branch has one complete unwind owner", {
  path <- locate_general_lp_fit_source()
  skip_if(is.null(path), "package sources unavailable")
  source <- readLines(path, warn = FALSE)
  region <- general_lp_owner_region(source)
  callsite <- general_lp_callsite_region(source)

  expect_equal(
    sum(gregexpr("R_UnwindProtect(", region, fixed = TRUE)[[1L]] > 0L),
    1L
  )
  expect_match(region, "NPLPSolveWorkspace solve_workspace;", fixed = TRUE)
  expect_match(
    region,
    "np_lp_solve_workspace_clear(&owner->solve_workspace);",
    fixed = TRUE
  )
  expect_match(
    region,
    "free_mat(owner->basis, owner->nterms);",
    fixed = TRUE
  )
  expect_match(region, "NPRegMpiOwnerChunk mpi_owner_chunk;", fixed = TRUE)
  expect_match(region, "double *mpi_kernel_row;", fixed = TRUE)
  expect_match(
    region,
    "np_reg_mpi_owner_chunk_free(&owner->mpi_owner_chunk);",
    fixed = TRUE
  )
  for (resource in c(
    "response_columns", "basis_columns", "squared_response", "moments",
    "power2_moments", "retained_kernel_row", "coefficient",
    "power2_projection", "eval_basis", "eval_derivative", "terms"
  )) {
    expect_equal(
      sum(gregexpr(
        paste0("free(owner->", resource, ");"),
        region,
        fixed = TRUE
      )[[1L]] > 0L),
      1L,
      info = resource
    )
  }
  expect_false(grepl("double *projection;", region, fixed = TRUE))
  expect_match(
    region,
    "np_lp_solve_workspace_solve_factored(&owner->solve_workspace,",
    fixed = TRUE
  )
  expect_equal(
    sum(gregexpr(
      "np_regression_general_lp_fit(&general_lp_call)",
      callsite,
      fixed = TRUE
    )[[1L]] > 0L),
    1L
  )
  expect_false(grepl("malloc(", callsite, fixed = TRUE))
  expect_false(grepl("alloc_matd(", callsite, fixed = TRUE))
  expect_false(grepl("cleanup_glp_fit", callsite, fixed = TRUE))
})
