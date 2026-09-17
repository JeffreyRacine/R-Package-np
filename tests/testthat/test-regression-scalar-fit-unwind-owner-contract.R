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
  # The cold ANN unknown-information replay has bounded R-transient scratch.
  # It is not manually owned heap or an unprotected SEXP. All other explicit
  # scalar allocations must stay under the existing unwind owner.
  replay.start <- regexpr("if(information_pass && replay_information) {",
                          owner, fixed = TRUE)[[1L]]
  replay.end <- regexpr("const int weighted_sum_status = kernel_weighted_sum_np_ctx_ex(",
                        owner, fixed = TRUE)[[1L]]
  expect_gt(replay.start, 0L)
  expect_gt(replay.end, replay.start)
  replay <- substr(owner, replay.start, replay.end - 1L)
  outside.replay <- paste(substr(owner, 1L, replay.start - 1L),
                           substring(owner, replay.end, nchar(owner)))
  expect_false(grepl("R_alloc(", outside.replay, fixed = TRUE))
  expect_false(grepl("R_allocLD(", outside.replay, fixed = TRUE))
  expect_identical(fixed_occurrences(replay, "R_alloc("), 2L)
  expect_identical(fixed_occurrences(replay, "R_allocLD("), 1L)
  compact <- gsub("[[:space:]]+", " ", owner)
  expect_match(compact,
    paste0("const int replay_information = ordinary_hc0 && call->do_gerr && p_nvar > 0 && ",
           "call->bandwidth_mode == BW_ADAP_NN && call->hc0_context->unknown_count > 0;"),
    fixed = TRUE)
  expect_match(compact,
    "information_pass <= (replay_information || replay_contrast); ++information_pass)",
    fixed = TRUE)
  expect_match(replay,
    "!np_size_mul_checked((size_t)p_nvar,(size_t)call->num_obs_eval,&cells)",
    fixed = TRUE)
  dimension.guard <- regexpr("cells > (size_t)INT_MAX/5U", replay,
                             fixed = TRUE)[[1L]]
  expect_gt(dimension.guard, 0L)
  expect_lt(dimension.guard, regexpr("R_alloc(", replay, fixed = TRUE)[[1L]])
  expect_match(replay, "information.level_denominator = denominators;", fixed = TRUE)
  expect_match(replay, "information.alternate_denominator = alternate_denominators;",
               fixed = TRUE)
  expect_match(replay, "information.certificate = R_allocLD(5U*cells);", fixed = TRUE)
  expect_match(replay, "memset(information.certificate,0,5U*cells*sizeof(long double));",
               fixed = TRUE)
  expect_match(replay, "information_dual.ann_information = &information;", fixed = TRUE)
  expect_match(outside.replay, "information_pass ? &information_dual : ordinary_hc0 ?",
               fixed = TRUE)
  execute <- npRmpi_test_extract_c_function(
    strsplit(engine, "\n", fixed = TRUE)[[1L]],
    "np_regression_scalar_fit_execute")
  expect_identical(fixed_occurrences(execute, "malloc("), 6L)
  for (resource in c(
    "mean_columns", "permutation_columns", "conditional_weights",
    "conditional_permutation_weights", "unit_response", "squared_response"
  )) {
    expect_identical(
      fixed_occurrences(owner, paste0("owner->", resource, " = NULL;")),
      1L,
      info = resource
    )
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
  expect_match(owner, "kernel_weighted_sum_np_ctx_ex(", fixed = TRUE)
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
    "NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_SCALAR, scalar_fit_status, \"conditional regression kernel traversal failed\");",
    fixed = TRUE
  )
  expect_match(
    engine,
    "NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_SCALAR, scalar_fit_status, \"conditional influence variance construction failed\");",
    fixed = TRUE
  )
})

test_that("conditional influence validates pointers before dereference", {
  path <- locate_scalar_fit_source()
  skip_if(is.null(path), "package sources unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")
  start <- regexpr(
    "static int NP_NOINLINE np_regression_conditional_influence_finish(",
    engine,
    fixed = TRUE
  )[[1L]]
  end <- regexpr(
    "#define NP_ACCUMULATE_SQUARE(value_)",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(start, 0L)
  expect_gt(end, start)
  finisher <- substr(engine, start, end - 1L)

  pointer.guard <- regexpr("weighted_sums == NULL", finisher, fixed = TRUE)[[1L]]
  dereference <- regexpr("denominator = weighted_sums[1];", finisher,
                        fixed = TRUE)[[1L]]
  value.guard <- regexpr("!R_FINITE(denominator)", finisher,
                        fixed = TRUE)[[1L]]

  expect_match(finisher, "double denominator;", fixed = TRUE)
  expect_false(grepl("const double denominator = weighted_sums[1];",
                     finisher, fixed = TRUE))
  expect_gt(pointer.guard, 0L)
  expect_gt(dereference, pointer.guard)
  expect_gt(value.guard, dereference)
  expect_identical(fixed_occurrences(finisher,
                                     "denominator = weighted_sums[1];"), 1L)
})
