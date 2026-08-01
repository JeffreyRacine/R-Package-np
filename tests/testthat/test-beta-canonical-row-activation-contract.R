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

  ingress_start <- regexpr("SEXP C_np_kernelsum(", ingress,
                           fixed = TRUE)[[1L]]
  ingress_end <- regexpr("SEXP C_np_kernelsum_power12(", ingress,
                         fixed = TRUE)[[1L]]
  expect_gt(ingress_start, 0L)
  expect_gt(ingress_end, ingress_start)
  kernelsum_ingress <- substr(ingress, ingress_start, ingress_end - 1L)

  activation_start <- regexpr(
    "NPContinuousKernelRoute route;",
    kernelsum_ingress,
    fixed = TRUE
  )[[1L]]
  activation_end <- regexpr(
    "if(derivative_dimension >= 0 || p_operator == OP_DERIVATIVE)",
    kernelsum_ingress,
    fixed = TRUE
  )[[1L]]
  expect_gt(activation_start, 0L)
  expect_gt(activation_end, activation_start)
  activation <- substr(
    kernelsum_ingress, activation_start, activation_end - 1L
  )

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
      "    weighted_permutation_sum, kw, pkw, 0, NULL, NULL);"
    ),
    fixed = TRUE
  )
  expect_match(
    engine,
    "if(!exact_beta_absolute_route)\n      return KWSNP_ERR_BADINVOC;",
    fixed = TRUE
  )
})

test_that("beta density and distribution fits have one canonical ingress", {
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

  ingress_start <- regexpr("SEXP C_np_density(", ingress, fixed = TRUE)[[1L]]
  ingress_end <- regexpr(
    "SEXP C_np_beta_conditional_bootstrap(", ingress, fixed = TRUE
  )[[1L]]
  expect_gt(ingress_start, 0L)
  expect_gt(ingress_end, ingress_start)
  density_ingress <- substr(ingress, ingress_start, ingress_end - 1L)

  expect_match(
    density_ingress,
    "active_route, active_diagnostics, 0);",
    fixed = TRUE
  )
  expect_false(grepl("np_beta_kernelsum(", density_ingress, fixed = TRUE))
  expect_false(grepl(
    "if(dens_or_dist == NP_DO_DENS)", density_ingress, fixed = TRUE
  ))

  engine_start <- regexpr(
    "void kernel_estimate_dens_dist_categorical_np(", engine, fixed = TRUE
  )[[1L]]
  engine_end <- regexpr(
    "int np_kernel_estimate_con_density_categorical_leave_one_out_cv(",
    engine, fixed = TRUE
  )[[1L]]
  expect_gt(engine_start, 0L)
  expect_gt(engine_end, engine_start)
  density_engine <- substr(engine, engine_start, engine_end - 1L)

  expect_match(density_engine, "const int exact_beta_route = kernel_route != NULL;",
               fixed = TRUE)
  expect_match(density_engine, "kernel_weighted_sum_np_route_power12(",
               fixed = TRUE)
  expect_match(density_engine, "kernel_weighted_sum_np_route_centered_m2(",
               fixed = TRUE)
  expect_match(density_engine, "np_progress_fit_loop_step);", fixed = TRUE)
  expect_match(
    density_engine,
    "if(beta_route_status != 0)\n      goto cleanup_density_fit;",
    fixed = TRUE
  )
  expect_match(
    density_engine,
    "error(\"canonical beta density/distribution row failed\");",
    fixed = TRUE
  )
})

test_that("signed-log beta rows compose every route segment", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  engine <- paste(
    readLines(file.path(root, "src", "continuous_kernel_row.c"),
              warn = FALSE),
    collapse = "\n"
  )
  multi_start <- regexpr(
    "np_continuous_kernel_beta_log_factor_row_multi_prevalidated(",
    engine, fixed = TRUE
  )[[1L]]
  resident_start <- regexpr(
    "np_continuous_kernel_beta_factor_row_with_log_factor(",
    engine, fixed = TRUE
  )[[1L]]
  resident_owner_start <- regexpr(
    "np_continuous_kernel_beta_log_factor_row(", engine, fixed = TRUE
  )[[1L]]
  expect_gt(multi_start, 0L)
  expect_gt(resident_owner_start, multi_start)
  expect_gt(resident_start, resident_owner_start)
  multi_owner <- substr(engine, multi_start, resident_owner_start - 1L)
  resident_owner <- substr(
    engine, resident_owner_start, resident_start - 1L
  )

  expect_match(multi_owner, "beta_segment_count", fixed = TRUE)
  expect_match(
    multi_owner,
    "segment_index < plan->route->segment_count",
    fixed = TRUE
  )
  expect_match(
    multi_owner,
    "np_continuous_kernel_log_channels_multiply(",
    fixed = TRUE
  )
  expect_match(
    resident_owner,
    "return np_continuous_kernel_beta_log_factor_row_multi_prevalidated(",
    fixed = TRUE
  )
})

test_that("canonical centered moments use an online training-order update", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  engine <- paste(
    readLines(file.path(root, "src", "continuous_kernel_row.c"),
              warn = FALSE),
    collapse = "\n"
  )
  owner_start <- regexpr(
    "np_continuous_kernel_beta_centered_moment_rows_validated(",
    engine, fixed = TRUE
  )[[1L]]
  owner_end <- regexpr(
    "NPContinuousKernelRowStatus np_continuous_kernel_scaled_restore(",
    engine, fixed = TRUE
  )[[1L]]
  expect_gt(owner_start, 0L)
  expect_gt(owner_end, owner_start)
  owner <- substr(engine, owner_start, owner_end - 1L)

  expect_match(owner, "sum[evaluation] += value;", fixed = TRUE)
  expect_match(owner, "delta = value - running_mean;", fixed = TRUE)
  expect_match(
    owner,
    "running_m2 += delta * (value - running_mean);",
    fixed = TRUE
  )
  expect_match(owner, "centered_m2[evaluation] = running_m2;", fixed = TRUE)
  expect_false(grepl("second_moment -", owner, fixed = TRUE))
})

test_that("centered moments have one activated fail-closed route boundary", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  header <- paste(
    readLines(file.path(root, "src", "headers.h"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(
    header, "int kernel_weighted_sum_np_route_centered_m2(", fixed = TRUE
  )
  expect_match(engine, "} NPCenteredMomentCtx;", fixed = TRUE)
  expect_match(
    engine,
    "beta_centered_moment && dual_power_ctx == NULL",
    fixed = TRUE
  )
  expect_match(
    engine,
    "ncol_Y == 0 && ncol_W == 0",
    fixed = TRUE
  )
  expect_match(
    engine,
    "centered_m2 == NULL ? NULL : &centered_moment_ctx",
    fixed = TRUE
  )
  occurrences <- gregexpr(
    "kernel_weighted_sum_np_route_centered_m2(", engine, fixed = TRUE
  )[[1L]]
  expect_length(occurrences[occurrences > 0L], 2L)
})
