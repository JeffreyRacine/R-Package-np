np_conditional_density_owner_root <- function() {
  root <- Sys.getenv("NP_SOURCE_ROOT", unset = "")
  if (nzchar(root))
    return(normalizePath(root, mustWork = TRUE))
  normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
}

np_conditional_density_owner_source <- function() {
  paste(readLines(
    file.path(np_conditional_density_owner_root(), "src", "np.c"),
    warn = FALSE
  ), collapse = "\n")
}

np_conditional_density_owner_function <- function(source, name) {
  start <- regexpr(
    paste0("static void\\s+", name, "\\s*\\("),
    source, perl = TRUE
  )
  expect_gt(start[[1L]], 0L)
  tail <- substring(source, start[[1L]])
  open <- regexpr("\\{", tail, perl = TRUE)[[1L]]
  expect_gt(open, 0L)
  chars <- strsplit(tail, "", fixed = TRUE)[[1L]]
  depth <- 0L
  close <- NA_integer_
  for (i in seq.int(open, length(chars))) {
    if (chars[[i]] == "{") depth <- depth + 1L
    if (chars[[i]] == "}") depth <- depth - 1L
    if (depth == 0L) {
      close <- i
      break
    }
  }
  expect_false(is.na(close))
  paste0(chars[seq_len(close)], collapse = "")
}

test_that("ordinary conditional-density bandwidth state has one typed owner", {
  source <- np_conditional_density_owner_source()
  destroy <- np_conditional_density_owner_function(
    source, "np_conditional_density_bw_prepared_context_destroy"
  )

  expect_match(source, "} NPConditionalDensityPreparedCtx;", fixed = TRUE)
  for (field in c(
    "matrix_y_unordered_original", "matrix_y_ordered_original",
    "matrix_y_continuous_original", "matrix_x_unordered_original",
    "matrix_x_ordered_original", "matrix_x_continuous_original",
    "matrix_xy_unordered_original", "matrix_xy_ordered_original",
    "matrix_xy_continuous_original", "num_categories",
    "matrix_categorical_vals", "continuous_stddev", "extendednn_upper",
    "powell_directions", "vector_scale_factor", "scale_factor_startbest",
    "powell_step", "ipt_x", "ipt_y", "ipt_xy", "tree_x", "tree_y",
    "tree_xy"
  )) {
    expect_match(destroy, paste0("context->", field), fixed = TRUE, info = field)
  }
  expect_match(
    source,
    "NPConditionalDensityPreparedCtx prepared_context = {0};",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_conditional_density_bw_prepared_context_destroy(&prepared_context);",
    fixed = TRUE
  )
})
