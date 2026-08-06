locate_scalar_fit_source <- function() {
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
  if (!length(paths)) NULL else paths[[1L]]
}

fixed_occurrences <- function(text, pattern) {
  matches <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  sum(matches > 0L)
}

test_that("the generic scalar fit has one complete unwind owner", {
  path <- locate_scalar_fit_source()
  skip_if(is.null(path), "package sources unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")
  start <- regexpr(
    "enum {\n  NP_REGRESSION_SCALAR_FIT_OK = 0,",
    engine,
    fixed = TRUE
  )[[1L]]
  end <- regexpr(
    "enum {\n  NP_REGRESSION_GENERAL_LP_FIT_OK = 0,",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(start, 0L)
  expect_gt(end, start)
  owner <- substr(engine, start, end - 1L)

  expect_identical(fixed_occurrences(owner, "R_UnwindProtect("), 1L)
  expect_match(owner, "NPRegressionScalarFitOwner owner;", fixed = TRUE)
  expect_false(grepl("R_alloc(", owner, fixed = TRUE))
  for (resource in c(
    "mean_columns", "permutation_columns", "conditional_weights",
    "conditional_permutation_weights", "unit_response", "squared_response"
  )) {
    expect_identical(
      fixed_occurrences(owner, paste0("free(owner->", resource, ");")),
      1L,
      info = resource
    )
  }
  expect_match(
    owner,
    "np_regression_conditional_influence_finish(",
    fixed = TRUE
  )
  expect_match(owner, "kernel_weighted_sum_np_ctx(", fixed = TRUE)
})

test_that("the generic scalar branch delegates without duplicate ownership", {
  path <- locate_scalar_fit_source()
  skip_if(is.null(path), "package sources unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")
  start <- regexpr(
    "if(lp_engine_est == NP_LP_ENGINE_SCALAR) { // canonical scalar LP0",
    engine,
    fixed = TRUE
  )[[1L]]
  end <- regexpr(
    "} else if(lp_engine_est == NP_LP_ENGINE_GENERAL)",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(start, 0L)
  expect_gt(end, start)
  branch <- substr(engine, start, end - 1L)

  expect_identical(
    fixed_occurrences(branch, "np_regression_scalar_fit(&scalar_call);"),
    1L
  )
  expect_false(grepl("malloc(", branch, fixed = TRUE))
  expect_false(grepl("free(", branch, fixed = TRUE))
  expect_match(
    branch,
    "if(scalar_fit_status != NP_REGRESSION_SCALAR_FIT_OK)",
    fixed = TRUE
  )
  expect_match(
    branch, "goto finish_regression_estimation;", fixed = TRUE
  )
  expect_match(
    engine,
    "error(\"conditional regression kernel traversal failed\");",
    fixed = TRUE
  )
  expect_match(
    engine,
    "error(\"conditional influence variance construction failed\");",
    fixed = TRUE
  )
})
