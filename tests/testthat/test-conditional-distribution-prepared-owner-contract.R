np_conditional_distribution_owner_source <- function() {
  root <- Sys.getenv("NP_SOURCE_ROOT", unset = "")
  if (!nzchar(root))
    root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
                          mustWork = TRUE)
  paste(readLines(file.path(root, "src", "np.c"), warn = FALSE),
        collapse = "\n")
}

np_conditional_distribution_owner_function <- function(source, name) {
  start <- regexpr(
    paste0("static void\\s+", name, "\\s*\\([^;]*\\)\\s*\\{"),
    source, perl = TRUE
  )
  expect_gt(start[[1L]], 0L)
  tail <- substring(source, start[[1L]])
  open <- regexpr("\\{", tail, perl = TRUE)[[1L]]
  chars <- strsplit(tail, "", fixed = TRUE)[[1L]]
  depth <- 0L
  for (i in seq.int(open, length(chars))) {
    if (chars[[i]] == "{") depth <- depth + 1L
    if (chars[[i]] == "}") depth <- depth - 1L
    if (depth == 0L)
      return(paste0(chars[seq_len(i)], collapse = ""))
  }
  fail(sprintf("could not extract %s", name))
}

test_that("conditional-distribution bandwidth state has one typed owner", {
  source <- np_conditional_distribution_owner_source()
  destroy <- np_conditional_distribution_owner_function(
    source, "np_conditional_distribution_prepared_context_destroy"
  )

  expect_match(source, "} NPConditionalDistributionPreparedCtx;", fixed = TRUE)
  for (field in c(
    "matrix_y_unordered_train", "matrix_y_ordered_train",
    "matrix_y_continuous_train", "matrix_y_unordered_eval",
    "matrix_y_ordered_eval", "matrix_y_continuous_eval",
    "matrix_x_unordered_train", "matrix_x_ordered_train",
    "matrix_x_continuous_train", "matrix_xy_unordered_train",
    "matrix_xy_ordered_train", "matrix_xy_continuous_train",
    "matrix_categorical_vals", "powell_directions", "scale_factor",
    "scale_factor_startbest", "powell_step", "extendednn_upper",
    "ipt_x", "ipt_y", "ipt_xy"
  )) {
    expect_match(destroy, paste0("context->", field), fixed = TRUE,
                 info = field)
  }
  expect_match(source, "NPConditionalDistributionPreparedCtx local_context = {0};",
               fixed = TRUE)
  expect_match(
    source, "np_conditional_distribution_prepared_context_destroy(prepared_context);",
    fixed = TRUE
  )
})

test_that("native conditional-distribution search retains one prepared owner", {
  source <- np_conditional_distribution_owner_source()

  expect_match(
    source,
    "NPConditionalDistributionPreparedCtx prepared;",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_conditional_distribution_prepared_context_eval(\n      &context->prepared",
    fixed = TRUE
  )
  expect_match(
    source,
    "1, &context.prepared, 1, ndegree > 0);",
    fixed = TRUE
  )
  for (field in c(
    "eval_bandwidth", "eval_degree", "glp_degree",
    "matrix_y_continuous_tree", "matrix_xy_continuous_tree"
  )) {
    expect_match(source, paste0("context->", field), fixed = TRUE,
                 info = field)
  }
})
