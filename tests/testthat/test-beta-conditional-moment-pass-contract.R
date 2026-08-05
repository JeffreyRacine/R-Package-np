locate_conditional_moment_source <- function() {
  roots <- unique(c(
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    getwd(),
    file.path(getwd(), "..")
  ))
  paths <- file.path(
    roots[nzchar(roots)],
    "src",
    "continuous_kernel_row.c"
  )
  paths <- paths[file.exists(paths)]
  if (!length(paths))
    return(NULL)
  paths[[1L]]
}

count_fixed <- function(text, pattern) {
  matches <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  if (identical(matches[[1L]], -1L)) 0L else length(matches)
}

test_that("conditional moments fuse validation into the influence pass", {
  path <- locate_conditional_moment_source()
  skip_if(is.null(path), "package source unavailable")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  regression.start <- regexpr(
    "np_continuous_kernel_beta_regression_moment_rows_validated(",
    source,
    fixed = TRUE
  )[[1L]]
  influence.start <- regexpr(
    "np_continuous_kernel_beta_conditional_influence_stderr(",
    source,
    fixed = TRUE
  )[[1L]]
  conditional.start <- regexpr(
    "np_continuous_kernel_beta_conditional_moment_rows_validated(",
    source,
    fixed = TRUE
  )[[1L]]
  expect_gt(regression.start, 0L)
  expect_lt(regression.start, influence.start)
  expect_lt(influence.start, conditional.start)

  regression <- substr(source, regression.start, influence.start - 1L)
  influence <- substr(source, influence.start, conditional.start - 1L)
  conditional <- substr(source, conditional.start, nchar(source))

  expect_match(
    regression,
    "squared_weight_sum += weight * weight;",
    fixed = TRUE
  )
  expect_match(influence, "double *validation_m2", fixed = TRUE)
  expect_match(
    influence,
    "*validation_m2 += weight * residual * residual;",
    fixed = TRUE
  )
  expect_match(
    influence,
    "if(!R_FINITE(*validation_m2))",
    fixed = TRUE
  )
  expect_false(grepl("validation_m2 += weight * delta *", conditional,
                     fixed = TRUE))
  expect_match(conditional, "&validation_m2", fixed = TRUE)
  expect_false(grepl("squared_weight_sum", conditional, fixed = TRUE))

  signed.start <- regexpr(
    "    } else {\n      double weighted_response_sum = 0.0;",
    conditional,
    fixed = TRUE
  )[[1L]]
  signed.end <- regexpr(
    "\n    }\n\n    if(!R_FINITE(total_weight)",
    conditional,
    fixed = TRUE
  )[[1L]]
  expect_gt(signed.start, 0L)
  expect_gt(signed.end, signed.start)
  signed <- substr(conditional, signed.start, signed.end - 1L)
  expect_equal(
    count_fixed(
      signed,
      "for(observation = 0; observation < plan->num_train; ++observation)"
    ),
    1L
  )
})
