np_regression_prepared_source <- function() {
  source_file <- file.path(npRmpi_namespace_hygiene_root(), "src", "np.c")
  expect_true(file.exists(source_file))
  paste(readLines(source_file, warn = FALSE), collapse = "\n")
}

np_regression_prepared_function <- function(source, name, return_type = "void") {
  markers <- gregexpr(
    sprintf("(?:static\\s+)?%s\\s+%s\\s*\\(", return_type, name),
    source,
    perl = TRUE
  )[[1L]]
  expect_false(identical(markers, -1L))
  marker <- NA_integer_
  open <- NA_integer_
  for (candidate in markers) {
    tail <- substring(source, candidate)
    candidate_open <- regexpr("\\{", tail, perl = TRUE)[[1L]]
    candidate_semi <- regexpr(";", tail, fixed = TRUE)[[1L]]
    if (candidate_open > 0L &&
        (candidate_semi <= 0L || candidate_open < candidate_semi)) {
      marker <- candidate
      open <- candidate_open
      break
    }
  }
  expect_false(is.na(marker))
  tail <- substring(source, marker)
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

test_that("MPI regression bandwidth state has one typed cleanup owner", {
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
    "powell_step", "degree", "ckerlb", "ckerub",
    "tree_permutation", "tree_lookup", "tree"
  )) {
    expect_match(owner, paste0("context->", field), fixed = TRUE, info = field)
  }
  expect_match(
    routine,
    "NPRegressionPreparedCtx local_context = {0};",
    fixed = TRUE
  )
  expect_match(
    routine,
    "NPRegressionPreparedCtx *prepared_context =",
    fixed = TRUE
  )
  expect_match(
    routine,
    "np_regression_prepared_context_destroy(prepared_context);",
    fixed = TRUE
  )
})

test_that("MPI native regression search evaluates one retained prepared state", {
  source <- np_regression_prepared_source()
  callback <- np_regression_prepared_function(
    source, "np_regression_prepared_native_search_callback", "int"
  )
  evaluator <- np_regression_prepared_function(
    source, "np_regression_prepared_eval_native_raw", "int"
  )
  prepare <- np_regression_prepared_function(
    source, "np_regression_prepared_prepare_internal", "int"
  )

  expect_match(
    evaluator,
    "np_regression_prepared_context_eval(",
    fixed = TRUE
  )
  expect_false(grepl("NP_NOMAD_CALLBACK_CALLOC", callback, fixed = TRUE))
  expect_false(grepl("NP_NOMAD_CALLBACK_FREE", callback, fixed = TRUE))
  expect_match(
    prepare,
    "&np_regression_prepared, 1, degree_search[0]);",
    fixed = TRUE
  )
  expect_false(grepl("NPRegressionNomadShadowCtx", source, fixed = TRUE))
  expect_false(grepl("np_regression_nomad_shadow", source, fixed = TRUE))
  expect_false(grepl("np_regression_prepared_shadow", source, fixed = TRUE))
})

test_that("MPI prepared-only regression setup elides Powell-only geometry", {
  source <- np_regression_prepared_source()
  routine <- np_regression_prepared_function(source, "np_regression_bw_mode")

  expect_match(
    routine,
    "matrix_y = prepare_only ? NULL :",
    fixed = TRUE
  )
  expect_match(
    routine,
    "vector_scale_factor_startbest = prepare_only ? NULL :",
    fixed = TRUE
  )
  expect_match(
    routine,
    "vsfh = prepare_only ? NULL : alloc_vecd",
    fixed = TRUE
  )
  expect_match(routine, "if(!prepare_only){", fixed = TRUE)
})

test_that("MPI regression prepared API has no obsolete shadow registrations", {
  root <- npRmpi_namespace_hygiene_root()
  init <- paste(readLines(file.path(root, "src", "np_init.c"), warn = FALSE),
                collapse = "\n")
  regression_r <- paste(
    readLines(file.path(root, "R", "np.regression.bw.R"), warn = FALSE),
    collapse = "\n"
  )

  for (name in c(
    "C_np_regression_prepared_prepare",
    "C_np_regression_prepared_eval",
    "C_np_regression_prepared_native_search",
    "C_np_regression_prepared_clear"
  )) {
    expect_match(init, name, fixed = TRUE)
    expect_match(regression_r, name, fixed = TRUE)
  }
  expect_false(grepl("C_np_regression_nomad_shadow", init, fixed = TRUE))
  expect_false(grepl("C_np_regression_nomad_shadow", regression_r, fixed = TRUE))
  expect_false(grepl("npregbw_nomad_shadow", regression_r, fixed = TRUE))
})
