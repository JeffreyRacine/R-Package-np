# Load only the named harness function; never execute its numerical tests.

np_source_tail_test_function <- function(file, names) {
  expressions <- parse(testthat::test_path(file), keep.source = FALSE)
  selected <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% names &&
      is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function"))
  }, as.list(expressions))
  expect_length(selected, 1L)
  environment <- new.env(parent = parent.frame())
  # Exercise the historical finite default even on R versions with last=NULL.
  environment$substring <- function(text, first, last = 1000000L)
    base::substring(text, first, last)
  # Evaluating this selected function definition binds its actual source body
  # without evaluating sibling test_that blocks or starting estimator/MPI work.
  eval(selected[[1L]], envir = environment)
}

test_that("H6 source extraction is complete and fails closed beyond one million", {
  extract <- np_source_tail_test_function(
    "test-regression-hc0-general-lp-derivative-contract.R", "h6_extract_lp_owner"
  )
  start <- "static SEXP np_regression_general_lp_fit_execute(void *data)"
  finish <- "static int np_regression_general_lp_fit(\n"
  target <- paste0(start, " {\n  required_covariance_owner();\n}\n")
  padding <- paste0(strrep(" ", 1000001L), "\n")
  source <- paste0(padding, target, finish)
  expect_identical(extract(source), target)
  expect_error(extract(paste0(padding, target)), "unique")
  expect_error(extract(paste0(padding, finish)), "unique")
  expect_error(extract(paste0(source, target)), "unique")
  expect_error(extract(paste0(source, finish)), "unique")
  expect_error(extract(paste0(padding, finish, target)), "out of order")
  broken <- sub("required_covariance_owner", "broken_covariance_owner",
                source, fixed = TRUE)
  expect_failure(expect_match(extract(broken), "required_covariance_owner",
                              fixed = TRUE))
})

test_that("C source extractors distinguish declarations and unique full definitions", {
  cases <- list(
    c("test-regression-prepared-owner-contract.R",
      "np_regression_prepared_function"),
    c("test-distribution-prepared-owner-contract.R",
      "np_distribution_prepared_function"),
    c("test-conditional-density-prepared-owner-contract.R",
      "np_conditional_density_owner_function",
      "np_conditional_density_prepared_function"),
    c("test-conditional-distribution-prepared-owner-contract.R",
      "np_conditional_distribution_owner_function"),
    c("test-native-nomad-callback-source-contract.R",
      "np_extract_c_function_body", "npRmpi_extract_c_function_body")
  )
  declaration <- "static void contract_target(void);\n"
  definition <- "static void contract_target(void) {\n  required_call();\n}\n"
  padding <- paste0(strrep(" ", 1000001L), "\n")
  source <- paste0(declaration, padding, definition)
  for (case in cases) {
    extract <- np_source_tail_test_function(case[[1L]], case[-1L])
    expect_match(extract(source, "contract_target"), "required_call();",
                 fixed = TRUE, info = case[[1L]])
    expect_error(extract(paste0(source, definition), "contract_target"),
                 "definition", info = case[[1L]])
    expect_error(extract(paste0(declaration, padding), "contract_target"),
                 "definition", info = case[[1L]])
    expect_error(extract(source, "missing_target"),
                 "missing|definition", info = case[[1L]])
    broken <- sub("required_call", "broken_call", source, fixed = TRUE)
    broken_target <- extract(broken, "contract_target")
    expect_failure(expect_match(broken_target, "required_call();", fixed = TRUE))
  }
})

test_that("profile source regions require a unique definition and ordered endpoint", {
  extract <- np_source_tail_test_function(
    "test-conditional-profile-prepared-index-contract.R",
    "extract_profile_implementation"
  )
  start <- "static void profile_target("
  finish <- "\nstatic void profile_next("
  definition <- "static void profile_target(void) { required_profile_call(); }"
  source <- paste0("static void profile_target(void);\n",
                   strrep(" ", 1000001L), "\n", definition, finish)
  expect_identical(extract(source, start, finish), definition)
  expect_error(extract(paste0(source, definition), start, finish), "one start")
  expect_error(extract(paste0(source, finish), start, finish), "one ordered")
  expect_error(extract(definition, start, finish), "one ordered")
  expect_error(extract(paste0(finish, definition), start, finish), "one ordered")
  broken <- sub("required_profile_call", "broken_profile_call", source, fixed = TRUE)
  expect_failure(expect_match(extract(broken, start, finish),
                              "required_profile_call", fixed = TRUE))
})

test_that("MPI paired source regions reject ambiguous anchors beyond one million", {
  all_large <- np_source_tail_test_function(
    "test-conditional-alllarge-cvls-mpi-ownership-contract.R",
    "extract_conditional_alllarge_region"
  )
  collective <- np_source_tail_test_function(
    "test-progress-mpi-collective-contract.R", "npRmpi_collective_source_region"
  )
  padding <- strrep(" ", 1000001L)
  source <- paste0(padding, "START required_owner(); FINISH")
  for (extract in list(all_large, collective)) {
    expect_identical(extract(source, "START", "FINISH"),
                     "START required_owner(); ")
    expect_error(extract(paste0(source, "START"), "START", "FINISH"), "start")
    expect_error(extract(paste0(padding, "START"), "START", "FINISH"),
                 "finish|stop")
    expect_error(extract(paste0(padding, "FINISH START"), "START", "FINISH"),
                 "finish|stop")
    broken <- sub("required_owner", "broken_owner", source, fixed = TRUE)
    target <- extract(broken, "START", "FINISH")
    expect_failure(expect_match(target, "required_owner", fixed = TRUE))
  }
  expect_error(all_large(paste0(source, "FINISH"), "START", "FINISH"), "finish")
  # A collective wrapper's generic next-section marker may legitimately recur.
  expect_identical(collective(paste0(source, " FINISH"), "START", "FINISH"),
                   "START required_owner(); ")
})
