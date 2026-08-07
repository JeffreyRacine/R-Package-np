route_unwind_source_root <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  roots <- roots[nzchar(roots)]
  hits <- roots[
    file.exists(file.path(roots, "src", "jksum.c")) &
      file.exists(file.path(roots, "src", "continuous_kernel_row.c"))
  ]
  if (!length(hits)) NULL else hits[[1L]]
}

fixed_occurrences <- function(text, pattern) {
  matches <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  sum(matches > 0L)
}

test_that("the canonical beta row route has one unwind owner", {
  root <- route_unwind_source_root()
  skip_if(is.null(root), "package sources unavailable")
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )

  owner_start <- regexpr("} NPBetaAbsoluteRouteCall;", engine,
                         fixed = TRUE)[[1L]]
  owner_end <- regexpr(
    "This function takes a vector Y and returns a kernel weighted",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(owner_start, 0L)
  expect_gt(owner_end, owner_start)
  owner <- substr(engine, owner_start, owner_end - 1L)

  expect_identical(fixed_occurrences(owner, "R_UnwindProtect("), 1L)
  expect_identical(
    fixed_occurrences(owner, "np_beta_absolute_route_owner_cleanup"),
    2L
  )
  expect_match(owner, "free(owner->route_row);", fixed = TRUE)
  expect_match(
    owner,
    "free_tmat(owner->owned_bandwidth_tmatrix);",
    fixed = TRUE
  )
  expect_match(
    owner,
    paste0(
      "execution.owner.owned_bandwidth_tmatrix =\n",
      "    call->owned_bandwidth_tmatrix;"
    ),
    fixed = TRUE
  )
  expect_match(
    owner,
    "np_continuous_kernel_level_derivative_workspace_release(",
    fixed = TRUE
  )
  expect_match(
    owner,
    "np_continuous_kernel_derivative_accumulator_release(",
    fixed = TRUE
  )
  expect_match(
    owner,
    "np_beta_categorical_factor_context_release(",
    fixed = TRUE
  )
  expect_false(grepl("free(route_row);", engine, fixed = TRUE))
  # The public kernel-sum consumer enters through the MPI-aware dispatcher;
  # the regression-fit consumer enters the same canonical owner directly.
  expect_identical(
    fixed_occurrences(engine, "np_beta_absolute_route(&route_call);"),
    1L
  )
  expect_identical(
    fixed_occurrences(
      engine,
      "np_beta_absolute_route_dispatch(\n        &route_call, suppress_parallel);"
    ),
    1L
  )
  expect_match(
    engine,
    "return np_beta_absolute_route(call);",
    fixed = TRUE
  )

  fit_start <- regexpr(
    "static NP_NOINLINE int np_beta_scalar_regression_fit_block_canonical(",
    engine,
    fixed = TRUE
  )[[1L]]
  fit_end <- regexpr(
    "NP_NOINLINE NP_COLD int np_beta_continuous_bandwidth_prepare_canonical(",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(fit_start, 0L)
  expect_gt(fit_end, fit_start)
  fit <- substr(engine, fit_start, fit_end - 1L)
  expect_match(
    fit,
    ".owned_bandwidth_tmatrix = matrix_bandwidth,",
    fixed = TRUE
  )
  route_call <- regexpr(
    "route_status = np_beta_absolute_route(&route_call);",
    fit,
    fixed = TRUE
  )[[1L]]
  expect_gt(route_call, 0L)
  expect_false(grepl(
    "free_tmat(matrix_bandwidth);",
    substr(fit, route_call, nchar(fit)),
    fixed = TRUE
  ))
})

test_that("nested beta derivative rows borrow route-owned scratch", {
  root <- route_unwind_source_root()
  skip_if(is.null(root), "package sources unavailable")
  row_engine <- paste(
    readLines(file.path(root, "src", "continuous_kernel_row.c"),
              warn = FALSE),
    collapse = "\n"
  )

  factored_start <- regexpr(
    "np_continuous_kernel_beta_derivative_absolute_rows_with_log_factor_validated(",
    row_engine,
    fixed = TRUE
  )[[1L]]
  powered_start <- regexpr(
    "np_continuous_kernel_beta_derivative_powered_rows_validated(",
    row_engine,
    fixed = TRUE
  )[[1L]]
  powered_end <- regexpr(
    "np_continuous_kernel_beta_dual_power_rows_validated(",
    row_engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(factored_start, 0L)
  expect_gt(powered_start, factored_start)
  expect_gt(powered_end, powered_start)
  borrowed <- substr(row_engine, factored_start, powered_end - 1L)

  expect_match(
    borrowed,
    "NPContinuousKernelRowWorkspace *factor_workspace",
    fixed = TRUE
  )
  expect_match(
    borrowed,
    "NPContinuousKernelRowWorkspace *workspace",
    fixed = TRUE
  )
  expect_match(borrowed, "double *row_storage", fixed = TRUE)
  expect_match(borrowed, "double *factor_log_absolute", fixed = TRUE)
  expect_match(borrowed, "signed char *factor_sign", fixed = TRUE)
  expect_false(grepl("malloc(", borrowed, fixed = TRUE))
  expect_false(grepl("free(", borrowed, fixed = TRUE))
  expect_false(grepl(
    "np_continuous_kernel_row_workspace_init(",
    borrowed,
    fixed = TRUE
  ))
  expect_false(grepl(
    "np_continuous_kernel_row_workspace_release(",
    borrowed,
    fixed = TRUE
  ))
})
