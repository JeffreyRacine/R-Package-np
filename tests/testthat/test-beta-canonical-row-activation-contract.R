locate_beta_activation_sources <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  roots <- roots[nzchar(roots)]
  roots <- roots[file.exists(file.path(roots, "src", "jksum.c"))]
  if (!length(roots))
    return(NULL)
  roots[[1L]]
}

test_that("activated beta absolute rows enter the canonical central engine", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")

  ingress <- paste(
    readLines(file.path(root, "src", "np.c"), warn = FALSE),
    collapse = "\n"
  )
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  row_engine <- paste(
    readLines(file.path(root, "src", "continuous_kernel_row.c"),
              warn = FALSE),
    collapse = "\n"
  )

  activation_start <- regexpr(
    "NPContinuousKernelRoute route;",
    ingress,
    fixed = TRUE
  )[[1L]]
  activation_end <- regexpr(
    "if(derivative_dimension >= 0 || p_operator == OP_DERIVATIVE)",
    ingress,
    fixed = TRUE
  )[[1L]]
  expect_gt(activation_start, 0L)
  expect_gt(activation_end, activation_start)
  activation <- substr(ingress, activation_start, activation_end - 1L)

  expect_match(activation, "NPContinuousKernelDerivativeDiagnostics",
               fixed = TRUE)
  expect_match(activation, "kernel_weighted_sum_np_route(", fixed = TRUE)
  expect_match(activation, "&route", fixed = TRUE)
  expect_false(grepl(
    "beta_bandwidth_code == BW_FIXED || !has_overlap",
    activation,
    fixed = TRUE
  ))
  expect_false(grepl("np_beta_kernelsum(", activation, fixed = TRUE))
  expect_false(grepl("np_beta_kernelsum_derivative(", activation,
                     fixed = TRUE))
  derivative_sidecar_calls <- gregexpr(
    "np_beta_kernelsum_derivative(", ingress, fixed = TRUE
  )[[1L]]
  derivative_sidecar_calls <- derivative_sidecar_calls[
    derivative_sidecar_calls > 0L
  ]
  expect_length(derivative_sidecar_calls, 0L)
  expect_false(grepl("KWS_DOTREEI", activation, fixed = TRUE))

  route_start <- regexpr("if(kernel_execution_context != NULL)", engine,
                         fixed = TRUE)[[1L]]
  legacy_start <- regexpr(
    "This function takes a vector Y and returns a kernel weighted",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(route_start, 0L)
  expect_gt(legacy_start, route_start)
  expect_match(engine, "np_beta_absolute_route(", fixed = TRUE)
  expect_match(engine, "np_continuous_kernel_beta_factor_row(", fixed = TRUE)
  expect_match(
    engine,
    "np_continuous_kernel_beta_derivative_absolute_rows_validated(",
    fixed = TRUE
  )
  expect_match(
    engine,
    "np_continuous_kernel_beta_derivative_powered_rows_validated(",
    fixed = TRUE
  )
  expect_match(
    row_engine,
    "if(response_columns == 0 && weight_columns == 0)",
    fixed = TRUE
  )
  expect_match(engine, "np_continuous_kernel_signed_log_restore(", fixed = TRUE)
  expect_match(
    engine,
    "np_continuous_kernel_signed_log_power_restore(",
    fixed = TRUE
  )
  expect_match(
    row_engine,
    "NPContinuousKernelRowStatus np_continuous_kernel_signed_log_power_restore(",
    fixed = TRUE
  )
  scaled_restore_start <- regexpr(
    "NPContinuousKernelRowStatus np_continuous_kernel_scaled_restore(",
    row_engine,
    fixed = TRUE
  )[[1L]]
  signed_restore_start <- regexpr(
    "NPContinuousKernelRowStatus np_continuous_kernel_signed_log_restore(",
    row_engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(scaled_restore_start, 0L)
  expect_gt(signed_restore_start, scaled_restore_start)
  scaled_restore <- substr(
    row_engine, scaled_restore_start, signed_restore_start - 1L
  )
  expect_match(
    scaled_restore,
    "np_continuous_kernel_signed_log_power_restore(",
    fixed = TRUE
  )
  expect_false(grepl(
    "log(fabs(scaled_value)) + (double)power * log_scale",
    scaled_restore,
    fixed = TRUE
  ))
  expect_match(
    engine,
    "np_continuous_kernel_beta_dual_power_rows_validated(",
    fixed = TRUE
  )
  expect_match(
    ingress,
    "kernel_weighted_sum_np_route_power12(",
    fixed = TRUE
  )
  expect_false(grepl(
    "beta kernels do not support the internal dual-power route",
    paste(ingress, engine),
    fixed = TRUE
  ))
  expect_match(
    engine,
    "(!route_has_convolution || matrix_bw_train != NULL)",
    fixed = TRUE
  )
  expect_match(
    engine,
    "(!route_has_convolution || matrix_bw_eval != NULL)",
    fixed = TRUE
  )
})

test_that("legacy callers cannot acquire beta route metadata implicitly", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(
    engine,
    paste0(
      "matrix_categorical_vals, matrix_ordered_indices, weighted_sum,\n",
      "    weighted_permutation_sum, kw, pkw, NULL, NULL);"
    ),
    fixed = TRUE
  )
  expect_match(
    engine,
    "if(!exact_beta_absolute_route)\n      return KWSNP_ERR_BADINVOC;",
    fixed = TRUE
  )
})
