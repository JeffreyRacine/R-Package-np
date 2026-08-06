np_regression_prepared_source <- function() {
  source_file <- file.path(np_namespace_hygiene_root(), "src", "np.c")
  expect_true(file.exists(source_file))
  paste(readLines(source_file, warn = FALSE), collapse = "\n")
}

np_regression_prepared_function <- function(source, name) {
  marker <- regexpr(
    sprintf("static void %s\\s*\\(", name),
    source,
    perl = TRUE
  )
  expect_gt(marker[[1L]], 0L)
  tail <- substring(source, marker[[1L]])
  open <- regexpr("\\{", tail, perl = TRUE)[[1L]]
  expect_gt(open, 0L)
  chars <- strsplit(substring(tail, open), "", fixed = TRUE)[[1L]]
  depth <- 0L
  for (i in seq_along(chars)) {
    if (chars[[i]] == "{")
      depth <- depth + 1L
    else if (chars[[i]] == "}") {
      depth <- depth - 1L
      if (depth == 0L)
        return(substring(tail, 1L, open + i - 1L))
    }
  }
  fail(sprintf("could not extract %s", name))
}

test_that("regression bandwidth state has one typed cleanup owner", {
  source <- np_regression_prepared_source()
  owner <- np_regression_prepared_function(
    source, "np_regression_prepared_context_destroy"
  )
  routine <- np_regression_prepared_function(source, "np_regression_bw_mode")

  expect_match(source, "} NPRegressionPreparedCtx;", fixed = TRUE)
  for (field in c(
    "matrix_x_unordered", "matrix_x_ordered", "matrix_x_continuous",
    "response", "lsq_scale", "lsq_q", "num_categories",
    "matrix_categorical_vals", "continuous_stddev", "extendednn_upper",
    "powell_directions", "scale_factor", "scale_factor_startbest",
    "powell_step", "tree_permutation", "tree_lookup", "tree"
  )) {
    expect_match(owner, paste0("context->", field), fixed = TRUE, info = field)
  }
  expect_match(
    routine,
    "NPRegressionPreparedCtx prepared_context = {0};",
    fixed = TRUE
  )
  expect_match(
    routine,
    "np_regression_prepared_context_destroy(&prepared_context);",
    fixed = TRUE
  )
})

