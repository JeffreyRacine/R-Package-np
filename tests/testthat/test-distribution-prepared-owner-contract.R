np_distribution_prepared_source <- function() {
  source_file <- file.path(npRmpi_namespace_hygiene_root(), "src", "np.c")
  expect_true(file.exists(source_file))
  paste(readLines(source_file, warn = FALSE), collapse = "\n")
}

np_distribution_prepared_function <- function(source, name,
                                               return_type = "void") {
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

test_that("distribution bandwidth state has one typed cleanup owner", {
  source <- np_distribution_prepared_source()
  owner <- np_distribution_prepared_function(
    source, "np_distribution_prepared_context_destroy"
  )
  routine <- np_distribution_prepared_function(
    source, "np_distribution_bw_internal"
  )
  public <- np_distribution_prepared_function(source, "np_distribution_bw")

  expect_match(source, "} NPDistributionPreparedCtx;", fixed = TRUE)
  for (field in c(
    "matrix_x_unordered_train", "matrix_x_ordered_train",
    "matrix_x_continuous_train", "matrix_x_unordered_eval",
    "matrix_x_ordered_eval", "matrix_x_continuous_eval",
    "num_categories", "matrix_categorical_vals", "continuous_stddev",
    "extendednn_upper", "powell_directions", "scale_factor",
    "scale_factor_multistart", "scale_factor_startbest", "powell_step",
    "tree_permutation_train", "tree_permutation_eval", "tree"
  )) {
    expect_match(owner, paste0("context->", field), fixed = TRUE, info = field)
  }
  expect_match(owner, "if (!context->cdfontrain)", fixed = TRUE)
  expect_match(
    routine,
    "NPDistributionPreparedCtx local_context = {0};",
    fixed = TRUE
  )
  expect_match(
    routine,
    "NPDistributionPreparedCtx *prepared_context =",
    fixed = TRUE
  )
  expect_match(routine, "if (prepare_only)", fixed = TRUE)
  expect_match(routine, "return;", fixed = TRUE)
  expect_match(
    routine,
    "np_distribution_prepared_context_destroy(prepared_context);",
    fixed = TRUE
  )
  expect_match(public, "np_distribution_bw_internal(", fixed = TRUE)
  expect_match(public, "eval_only, NULL, 0);", fixed = TRUE)
})

test_that("native distribution search evaluates one retained prepared state", {
  source <- np_distribution_prepared_source()
  callback <- np_distribution_prepared_function(
    source, "np_udist_native_search_callback", "int"
  )
  search <- np_distribution_prepared_function(
    source, "C_np_distribution_nomad_native_search", "SEXP"
  )

  expect_match(source, "NPDistributionPreparedCtx prepared;", fixed = TRUE)
  expect_match(
    callback,
    "np_distribution_prepared_context_eval(",
    fixed = TRUE
  )
  expect_false(grepl(
    "np_distribution_nomad_native_eval_once(", callback, fixed = TRUE
  ))
  expect_false(grepl("NP_NOMAD_CALLBACK_CALLOC", callback, fixed = TRUE))
  expect_match(
    search,
    "np_distribution_prepared_context_prepare(",
    fixed = TRUE
  )
  expect_match(
    search,
    "np_distribution_prepared_context_destroy(&context.prepared);",
    fixed = TRUE
  )
})
