locate_lp_objective_source <- function() {
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

lp_objective_owner_region <- function(source) {
  start <- grep("^} NPRegressionCvLPRouteOwner;$", source)
  stop <- grep(
    " * Isolated route-bearing sibling for canonical continuous-kernel activation",
    source,
    fixed = TRUE
  )
  stopifnot(length(start) == 1L, length(stop) == 1L, stop > start)
  paste(source[seq.int(start, stop - 1L)], collapse = "\n")
}

test_that("wider-LP regression objective has one complete unwind owner", {
  path <- locate_lp_objective_source()
  skip_if(is.null(path), "package sources unavailable")
  region <- lp_objective_owner_region(readLines(path, warn = FALSE))

  expect_equal(
    sum(gregexpr("R_UnwindProtect(", region, fixed = TRUE)[[1L]] > 0L),
    1L
  )
  expect_match(region, "NPReghatLPWorkspace lp_workspace;", fixed = TRUE)
  expect_match(region, "int matrix_bandwidth_columns;", fixed = TRUE)
  expect_match(
    region,
    "np_reghat_lp_workspace_clear(&owner->lp_workspace);",
    fixed = TRUE
  )
  expect_match(
    region,
    "free_mat(owner->matrix_bandwidth,\n             owner->matrix_bandwidth_columns);",
    fixed = TRUE
  )
  expect_match(
    region,
    "if(!jump)\n    np_continuous_kernel_beta_prepared_context_release",
    fixed = TRUE
  )
  expect_false(grepl(
    "np_beta_scaled_row_context_clear(&owner->row_context)",
    region,
    fixed = TRUE
  ))
  for (resource in c(
    "operator", "kernel_u", "kernel_o", "lambda", "row",
    "influence_row", "eval_basis"
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
})
