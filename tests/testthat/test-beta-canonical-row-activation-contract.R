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

  expect_false(grepl("np_beta_kernelsum(", density_ingress, fixed = TRUE))
  expect_false(grepl(
    "if(dens_or_dist == NP_DO_DENS)", density_ingress, fixed = TRUE
  ))
  expect_false(grepl(
    "currently supports continuous variables only",
    density_ingress,
    fixed = TRUE
  ))
  expect_match(
    density_ingress,
    "active_route, active_diagnostics, categorical_compress);",
    fixed = TRUE
  )

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

test_that("scalar beta regression fits enter the canonical row engine", {
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

  public_start <- regexpr("SEXP C_np_regression(", ingress,
                          fixed = TRUE)[[1L]]
  public_end <- regexpr("SEXP C_np_density(", ingress,
                        fixed = TRUE)[[1L]]
  expect_gt(public_start, 0L)
  expect_gt(public_end, public_start)
  public_regression <- substr(ingress, public_start, public_end - 1L)

  route_start <- regexpr(
    "if(descriptor.family == NP_CKERNEL_FAMILY_BETA)",
    public_regression,
    fixed = TRUE
  )[[1L]]
  common_call <- regexpr(
    "  np_regression(REAL(tuno_r)",
    public_regression,
    fixed = TRUE
  )[[1L]]
  expect_gt(route_start, 0L)
  expect_gt(common_call, route_start)
  route_ingress <- substr(
    public_regression, route_start, common_call - 1L
  )
  expect_match(route_ingress, "beta_route.segment_count = 1;", fixed = TRUE)
  expect_match(route_ingress, "active_route = &beta_route;", fixed = TRUE)
  expect_match(
    route_ingress, "active_diagnostics = &beta_diagnostics;", fixed = TRUE
  )
  expect_match(
    public_regression,
    "active_route, active_diagnostics, categorical_compress);",
    fixed = TRUE
  )
  expect_false(grepl("np_beta_regression_lc(", public_regression,
                     fixed = TRUE))
  expect_false(grepl("beta_regression.h", ingress, fixed = TRUE))
  expect_false(file.exists(file.path(root, "src", "beta_regression.c")))
  expect_false(file.exists(file.path(root, "src", "beta_regression.h")))

  owner_start <- regexpr("void np_regression(", ingress,
                         fixed = TRUE)[[1L]]
  owner_start <- regexpr(
    "void np_regression(",
    substr(ingress, owner_start + 1L, nchar(ingress)),
    fixed = TRUE
  )[[1L]] + owner_start
  owner_end <- regexpr("static void np_kernelsum_common(", ingress,
                       fixed = TRUE)[[1L]]
  expect_gt(owner_start, 0L)
  expect_gt(owner_end, owner_start)
  owner <- substr(ingress, owner_start, owner_end - 1L)
  expect_match(owner, "const NPContinuousKernelRoute *kernel_route",
               fixed = TRUE)
  expect_match(
    owner,
    paste0(
      "&SIGN,\n                                                   kernel_route,\n",
      "                                                   kernel_route_diagnostics,\n",
      "                                                   categorical_compress,\n",
      "                                                   NP_REGRESSION_STDERR_LOCAL_RESIDUAL,\n",
      "                                                   NULL);"
    ),
    fixed = TRUE
  )

  engine_start <- regexpr(
    "int kernel_estimate_regression_categorical_tree_np(",
    engine,
    fixed = TRUE
  )[[1L]]
  engine_end <- regexpr("static int np_conditional_indicator_row_core(",
                        engine, fixed = TRUE)[[1L]]
  expect_gt(engine_start, 0L)
  expect_gt(engine_end, engine_start)
  regression_engine <- substr(engine, engine_start, engine_end - 1L)
  sibling_start <- regexpr(
    "static NP_NOINLINE void np_beta_scalar_regression_fit_canonical(",
    engine, fixed = TRUE
  )[[1L]]
  expect_gt(sibling_start, 0L)
  expect_lt(sibling_start, engine_start)
  regression_sibling <- substr(engine, sibling_start, engine_start - 1L)
  expect_match(regression_engine, "if(NP_UNLIKELY(kernel_route != NULL)",
               fixed = TRUE)
  expect_match(
    regression_engine, "np_beta_scalar_regression_fit_canonical(",
    fixed = TRUE
  )
  expect_match(regression_engine, "np_regression_fit_statistics(",
               fixed = TRUE)
  expect_match(regression_engine, "return 0;", fixed = TRUE)
  expect_match(regression_sibling, "do_grad != do_gerr", fixed = TRUE)
  expect_match(
    regression_sibling, "np_beta_bandwidth_prepare_matrix(", fixed = TRUE
  )
  expect_match(regression_sibling, "NPBetaRegressionMomentCtx", fixed = TRUE)
  expect_match(
    regression_sibling,
    "static NP_NOINLINE int np_beta_regression_lp_moment_row_canonical(",
    fixed = TRUE
  )
  expect_match(
    regression_engine, "np_beta_regression_lp_moment_row_canonical(",
    fixed = TRUE
  )
  expect_match(
    engine, "np_beta_regression_gradient_rows_validated(",
    fixed = TRUE
  )
  expect_match(
    regression_sibling, "&regression_moment_context, kernel_route_diagnostics, NULL);",
    fixed = TRUE
  )
  expect_match(
    regression_sibling, "error(\"canonical beta regression row failed:",
    fixed = TRUE
  )
  expect_false(grepl("NPBetaRegressionMomentCtx", regression_engine,
                     fixed = TRUE))
  expect_false(grepl("np_beta_bandwidth_prepare_matrix(", regression_engine,
                     fixed = TRUE))
  expect_false(grepl("np_beta_regression_lc(", regression_engine,
                     fixed = TRUE))
})

test_that("canonical beta gradient rows separate operator and estimator algebra", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  row_engine <- paste(
    readLines(file.path(root, "src", "continuous_kernel_row_gradient.c"),
              warn = FALSE),
    collapse = "\n"
  )
  estimator_engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )

  observation_start <- regexpr(
    "np_continuous_kernel_beta_level_derivative_observation_bound(",
    row_engine,
    fixed = TRUE
  )[[1L]]
  row_start <- regexpr(
    "np_continuous_kernel_beta_level_derivative_log_row_validated(",
    row_engine,
    fixed = TRUE
  )[[1L]]
  row_end <- nchar(row_engine) + 1L
  expect_gt(observation_start, 0L)
  expect_gt(row_start, observation_start)
  expect_gt(row_end, row_start)
  operator_owner <- substr(
    row_engine, observation_start, row_start - 1L
  )
  row_owner <- substr(row_engine, row_start, row_end - 1L)
  expect_match(operator_owner, "np_beta_log_abs_pdf_order(", fixed = TRUE)
  expect_match(operator_owner, "np_beta_pdf_derivative_order(", fixed = TRUE)
  expect_match(operator_owner, "*common_log_scale = fmax(", fixed = TRUE)
  expect_match(
    row_owner,
    "np_continuous_kernel_beta_level_derivative_observation_bound(",
    fixed = TRUE
  )
  expect_match(
    row_owner, "provider == NULL && omitted_observation == -1",
    fixed = TRUE
  )
  expect_false(grepl("weighted_response", operator_owner, fixed = TRUE))
  expect_false(grepl("gradient_stderr", operator_owner, fixed = TRUE))
  expect_false(grepl("weighted_response", row_owner, fixed = TRUE))
  expect_false(grepl("gradient_stderr", row_owner, fixed = TRUE))

  consumer_start <- regexpr(
    "np_beta_regression_gradient_rows_validated(",
    estimator_engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(consumer_start, 0L)
  consumer <- substr(estimator_engine, consumer_start, nchar(estimator_engine))
  expect_match(consumer, "weighted_response += w * y;", fixed = TRUE)
  expect_match(
    consumer,
    "(regular_response - side_mean * regular_total) / side_weight;",
    fixed = TRUE
  )
  expect_match(consumer, "derivative_coefficient_square_sum", fixed = TRUE)
})

test_that("legacy conditional scalar owner retains dormant route plumbing", {
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

  owner_start <- regexpr(
    "void np_kernel_estimate_con_dens_dist_categorical(",
    engine,
    fixed = TRUE
  )[[1L]]
  owner_end <- regexpr(
    "attribute_hidden int np_fixed_gaussian_density_cvls_pair_dispatch_try(",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(owner_start, 0L)
  expect_gt(owner_end, owner_start)
  owner <- substr(engine, owner_start, owner_end - 1L)

  expect_match(owner, "const NPContinuousKernelRoute *kernel_route",
               fixed = TRUE)
  expect_match(owner, "(void)kernel_route;", fixed = TRUE)
  expect_match(owner, "(void)kernel_route_diagnostics;", fixed = TRUE)
  expect_match(owner, "(void)categorical_compress;", fixed = TRUE)
  expect_false(grepl("NPContinuousKernelExecutionContext", owner,
                     fixed = TRUE))

  conditional_starts <- gregexpr(
    "void np_density_conditional(", ingress, fixed = TRUE
  )[[1L]]
  conditional_start <- tail(conditional_starts[conditional_starts > 0L], 1L)
  density_starts <- gregexpr("void np_density(double", ingress,
                             fixed = TRUE)[[1L]]
  conditional_end <- tail(density_starts[density_starts > conditional_start],
                          1L)
  expect_gt(conditional_start, 0L)
  expect_gt(conditional_end, conditional_start)
  conditional <- substr(ingress, conditional_start, conditional_end - 1L)
  expect_match(
    conditional,
    paste0(
      "pdf_deriv_stderr,\n",
      "                                                 &log_likelihood,\n",
      "                                                 NULL,\n",
      "                                                 NULL,\n",
      "                                                 0);"
    ),
    fixed = TRUE
  )
})

test_that("beta density CVML enters the canonical scaled-row owner", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  ingress <- paste(
    readLines(file.path(root, "src", "kernelcv.c"), warn = FALSE),
    collapse = "\n"
  )
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )

  owner_start <- regexpr(
    "int np_kernel_estimate_density_categorical_leave_one_out_cv(",
    engine,
    fixed = TRUE
  )[[1L]]
  owner_end <- regexpr(
    "int np_kernel_estimate_density_categorical_convolution_cv(",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(owner_start, 0L)
  expect_gt(owner_end, owner_start)
  owner <- substr(engine, owner_start, owner_end - 1L)

  expect_match(owner, "const NPContinuousKernelRoute * const kernel_route",
               fixed = TRUE)
  expect_match(owner, "const int exact_beta_route = kernel_route != NULL;",
               fixed = TRUE)
  expect_match(owner, "np_density_cvml_beta_route(", fixed = TRUE)
  expect_match(
    owner,
    "if(exact_beta_route)\n      goto cleanup_density_leave_one_out_cv;",
    fixed = TRUE
  )
  expect_match(owner, "goto cleanup_density_leave_one_out_cv;",
               fixed = TRUE)
  expect_match(
    engine,
    "np_beta_scaled_row_context_fill_omitting(",
    fixed = TRUE
  )
  expect_match(engine, "np_guarded_cvml_log_contribution(", fixed = TRUE)

  callback_start <- regexpr(
    "double np_cv_func_density_categorical_ml(", ingress, fixed = TRUE
  )[[1L]]
  callback_end <- regexpr(
    "double np_cv_func_density_categorical_ls(", ingress, fixed = TRUE
  )[[1L]]
  expect_gt(callback_start, 0L)
  expect_gt(callback_end, callback_start)
  callback <- substr(ingress, callback_start, callback_end - 1L)
  expect_match(callback, "if(KERNEL_den_extern == NP_CKERNEL_COORDINATE_CODE)",
               fixed = TRUE)
  expect_false(grepl("np_beta_objective_density_ml(", callback,
                     fixed = TRUE))
  expect_match(callback, "active_route = &beta_route;", fixed = TRUE)
  expect_match(
    callback,
    paste0(
      "num_categories_extern,\n",
      "        active_route,\n",
      "        active_diagnostics,\n",
      "        np_density_bw_categorical_compress_extern,\n",
      "        &cv)==1)"
    ),
    fixed = TRUE
  )
})

test_that("beta density CVLS enters canonical quadrature and LOO owners", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  ingress <- paste(
    readLines(file.path(root, "src", "kernelcv.c"), warn = FALSE),
    collapse = "\n"
  )
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )

  owner_start <- regexpr(
    "int np_kernel_estimate_density_categorical_convolution_cv(",
    engine,
    fixed = TRUE
  )[[1L]]
  owner_end <- regexpr(
    "void kernel_estimate_dens_dist_categorical_np(",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(owner_start, 0L)
  expect_gt(owner_end, owner_start)
  owner <- substr(engine, owner_start, owner_end - 1L)

  expect_match(owner, "const NPContinuousKernelRoute * const kernel_route",
               fixed = TRUE)
  expect_match(owner, "const int exact_beta_route = kernel_route != NULL;",
               fixed = TRUE)
  expect_match(owner, "np_density_cvls_beta_cross_term(", fixed = TRUE)
  expect_match(
    owner,
    "if(exact_beta_route)\n      goto cleanup_density_convolution_cv;",
    fixed = TRUE
  )
  expect_match(owner, "goto cleanup_density_convolution_cv;", fixed = TRUE)
  expect_false(grepl("NPContinuousKernelExecutionContext", owner,
                     fixed = TRUE))

  callback_start <- regexpr(
    "double np_cv_func_density_categorical_ls(", ingress, fixed = TRUE
  )[[1L]]
  callback_end <- regexpr(
    "double cv_func_distribution_categorical_ls(", ingress, fixed = TRUE
  )[[1L]]
  expect_gt(callback_start, 0L)
  expect_gt(callback_end, callback_start)
  callback <- substr(ingress, callback_start, callback_end - 1L)
  expect_match(callback, "if(KERNEL_den_extern == NP_CKERNEL_COORDINATE_CODE)",
               fixed = TRUE)
  expect_false(grepl("np_beta_objective_density_ls(", callback,
                     fixed = TRUE))
  expect_match(callback, "active_route = &beta_route;", fixed = TRUE)
  expect_match(
    callback,
    paste0(
      "matrix_categorical_vals_extern,\n",
      "        active_route,\n",
      "        active_diagnostics,\n",
      "        np_density_bw_categorical_compress_extern,\n",
      "        &cv)==1)"
    ),
    fixed = TRUE
  )
})

test_that("beta distribution CVLS enters the canonical row owner", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  ingress <- paste(
    readLines(file.path(root, "src", "kernelcv.c"), warn = FALSE),
    collapse = "\n"
  )
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )

  owner_start <- regexpr(
    "double np_kernel_estimate_distribution_ls_cv(", engine, fixed = TRUE
  )[[1L]]
  owner_end <- regexpr(
    "typedef enum {\n  NP_CONDITIONAL_PROFILE_CV_SUCCESS", engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(owner_start, 0L)
  expect_gt(owner_end, owner_start)
  owner <- substr(engine, owner_start, owner_end - 1L)

  expect_match(owner, "const NPContinuousKernelRoute * const kernel_route",
               fixed = TRUE)
  expect_match(owner, "const int exact_beta_route = kernel_route != NULL;",
               fixed = TRUE)
  expect_match(owner, "np_distribution_cvls_continuous_route(", fixed = TRUE)
  expect_match(owner, "np_distribution_cvls_accumulate_row(", fixed = TRUE)
  expect_match(
    engine,
    paste0(
      "np_continuous_kernel_scaled_restore(\n",
      "         1.0, common_log_scale, 1, &common_scale)"
    ),
    fixed = TRUE
  )
  expect_match(engine, "row_sum = scaled_sum*common_scale;", fixed = TRUE)
  expect_match(engine, "row[observation] *= common_scale;", fixed = TRUE)
  expect_false(grepl("&row[observation]", engine, fixed = TRUE))
  expect_match(owner, "goto cleanup_distribution_ls_cv;", fixed = TRUE)
  expect_false(grepl("NPContinuousKernelExecutionContext", owner,
                     fixed = TRUE))

  callback_start <- regexpr(
    "double cv_func_distribution_categorical_ls(", ingress, fixed = TRUE
  )[[1L]]
  callback_end <- regexpr(
    "double func_con_density_quantile(", ingress, fixed = TRUE
  )[[1L]]
  expect_gt(callback_start, 0L)
  expect_gt(callback_end, callback_start)
  callback <- substr(ingress, callback_start, callback_end - 1L)
  expect_match(callback, "if(KERNEL_den_extern == NP_CKERNEL_COORDINATE_CODE)",
               fixed = TRUE)
  expect_false(grepl("np_beta_objective_distribution_ls(", callback,
                     fixed = TRUE))
  expect_match(callback, "active_route = &beta_route;", fixed = TRUE)
  expect_match(
    callback,
    paste0(
      "matrix_categorical_vals_extern,\n",
      "                                             active_route,\n",
      "                                             active_diagnostics,\n",
      "                                             np_distribution_bw_categorical_compress_extern,\n",
      "                                             &cv)==1)"
    ),
    fixed = TRUE
  )
})

test_that("distribution bandwidth ingress owns categorical compression state", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  distribution_r <- paste(
    readLines(file.path(root, "R", "np.distribution.bw.R"), warn = FALSE),
    collapse = "\n"
  )
  headers <- paste(
    readLines(file.path(root, "src", "headers.h"), warn = FALSE),
    collapse = "\n"
  )
  ingress <- paste(
    readLines(file.path(root, "src", "np.c"), warn = FALSE),
    collapse = "\n"
  )

  option_literal <- "categorical.compress = npStrictLogicalOption("
  option_hits <- gregexpr(option_literal, distribution_r, fixed = TRUE)[[1L]]
  option_hits <- option_hits[option_hits > 0L]
  expect_length(option_hits, 2L)
  expect_match(headers, "#define DBW_CATCOMPI 24", fixed = TRUE)
  expect_match(
    ingress,
    paste0(
      "np_distribution_bw_categorical_compress_extern = ",
      "myopti[DBW_CATCOMPI];"
    ),
    fixed = TRUE
  )
  expect_match(
    ingress,
    "np_distribution_bw_categorical_compress_extern = 0;",
    fixed = TRUE
  )
  expect_match(
    ingress,
    "C_np_distribution_bw: categorical compression must be TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("density bandwidth ingress owns categorical compression state", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  density_r <- paste(
    readLines(file.path(root, "R", "np.density.bw.R"), warn = FALSE),
    collapse = "\n"
  )
  bandwidth_r <- paste(
    readLines(file.path(root, "R", "bandwidth.R"), warn = FALSE),
    collapse = "\n"
  )
  headers <- paste(
    readLines(file.path(root, "src", "headers.h"), warn = FALSE),
    collapse = "\n"
  )
  ingress <- paste(
    readLines(file.path(root, "src", "np.c"), warn = FALSE),
    collapse = "\n"
  )

  option_literal <- "categorical.compress = npStrictLogicalOption("
  option_hits <- gregexpr(option_literal, density_r, fixed = TRUE)[[1L]]
  option_hits <- option_hits[option_hits > 0L]
  expect_length(option_hits, 2L)
  expect_match(headers, "#define BW_CATCOMPI 23", fixed = TRUE)
  expect_match(
    ingress,
    "np_density_bw_categorical_compress_extern = myopti[BW_CATCOMPI];",
    fixed = TRUE
  )
  expect_match(
    ingress,
    "np_density_bw_categorical_compress_extern = 0;",
    fixed = TRUE
  )
  expect_match(
    ingress,
    "C_np_density_bw: categorical compression must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_match(bandwidth_r, "allow.categorical = TRUE", fixed = TRUE)
  allow_hits <- gregexpr(
    "NP_BETA_BW_ALLOW_CATEGORICAL,", ingress, fixed = TRUE
  )[[1L]]
  allow_hits <- allow_hits[allow_hits > 0L]
  reject_hits <- gregexpr(
    "NP_BETA_BW_CONTINUOUS_ONLY,", ingress, fixed = TRUE
  )[[1L]]
  reject_hits <- reject_hits[reject_hits > 0L]
  expect_length(allow_hits, 1L)
  expect_length(reject_hits, 8L)
  contract_hits <- gregexpr(
    "np_density_bw_integer_contract_or_error(", ingress, fixed = TRUE
  )[[1L]]
  contract_hits <- contract_hits[contract_hits > 0L]
  expect_length(contract_hits, 4L)
})

test_that("every beta side enters the common conditional regression owner", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  ingress <- paste(
    readLines(file.path(root, "src", "np.c"), warn = FALSE),
    collapse = "\n"
  )
  public_start <- regexpr("SEXP C_np_density_conditional(", ingress,
                          fixed = TRUE)[[1L]]
  public_end <- regexpr("SEXP C_np_density_bw(", ingress,
                        fixed = TRUE)[[1L]]
  expect_gt(public_start, 0L)
  expect_gt(public_end, public_start)
  public_owner <- substr(ingress, public_start, public_end - 1L)
  expect_match(
    public_owner,
    "x_descriptor.family == NP_CKERNEL_FAMILY_BETA) {",
    fixed = TRUE
  )
  expect_match(public_owner, "active_x_route = &beta_x_route;", fixed = TRUE)
  expect_match(
    public_owner,
    "const NPContinuousKernelRoute *active_y_route = NULL;",
    fixed = TRUE
  )
  expect_match(
    public_owner,
    "y_descriptor.family == NP_CKERNEL_FAMILY_BETA) {",
    fixed = TRUE
  )
  expect_match(public_owner, "active_y_route = &beta_y_route;", fixed = TRUE)
  expect_match(
    public_owner,
    "active_y_route, active_y_diagnostics,",
    fixed = TRUE
  )
  expect_match(
    public_owner,
    "y_descriptor.family == NP_CKERNEL_FAMILY_BETA",
    fixed = TRUE
  )
  expect_false(grepl(
    "x_descriptor.family == NP_CKERNEL_FAMILY_BETA ||\n      y_descriptor.family == NP_CKERNEL_FAMILY_BETA",
    public_owner,
    fixed = TRUE
  ))
  expect_false(grepl("np_beta_conditional_lc(", public_owner, fixed = TRUE))
  expect_false(grepl(
    "np_beta_conditional_lc_gradient(", public_owner, fixed = TRUE
  ))
  expect_false(grepl(
    "beta conditional estimators currently support", public_owner,
    fixed = TRUE
  ))
  expect_match(
    public_owner,
    "np_density_conditional(REAL(tyuno_r), REAL(tyord_r), REAL(tycon_r)",
    fixed = TRUE
  )

  conditional_starts <- gregexpr(
    "void np_density_conditional(", ingress, fixed = TRUE
  )[[1L]]
  conditional_start <- tail(conditional_starts[conditional_starts > 0L], 1L)
  density_starts <- gregexpr("void np_density(double", ingress,
                             fixed = TRUE)[[1L]]
  conditional_end <- tail(density_starts[density_starts > conditional_start],
                          1L)
  conditional <- substr(ingress, conditional_start, conditional_end - 1L)
  expect_match(conditional, "const int beta_y_active = response_kernel_route != NULL;",
               fixed = TRUE)
  expect_match(
    conditional,
    "error(\"np_density_conditional: invalid canonical response-kernel route\")",
    fixed = TRUE
  )
  expect_false(grepl(
    "kernel_route != NULL || response_kernel_route_diagnostics == NULL",
    conditional, fixed = TRUE
  ))
  expect_false(grepl(
    "kernel_weighted_sum_np_route(kernel_cy", conditional, fixed = TRUE
  ))
  expect_match(
    conditional,
    "if(lp_engine_eff == NP_LP_ENGINE_SCALAR && kernel_route == NULL &&\n     response_kernel_route == NULL)",
    fixed = TRUE
  )
  expect_match(
    conditional,
    "NP_REGRESSION_STDERR_CONDITIONAL_INFLUENCE",
    fixed = TRUE
  )
  expect_match(
    conditional,
    "kernel_route,\n                                                               kernel_route_diagnostics,\n                                                               categorical_compress",
    fixed = TRUE
  )
  expect_match(
    conditional,
    "np_beta_prepared_bandwidth_view_init_or_error(",
    fixed = TRUE
  )
  expect_match(
    conditional,
    "prepared_x_bandwidth.evaluation_offset = j;",
    fixed = TRUE
  )
  expect_match(
    conditional, "prepared_x_bandwidth_ptr);", fixed = TRUE
  )
  expect_match(
    conditional, "np_beta_continuous_bandwidth_prepare_canonical(",
    fixed = TRUE
  )
  expect_match(
    conditional, "np_beta_scaled_row_context_prepare(", fixed = TRUE
  )
  expect_match(
    conditional, "np_beta_scaled_row_context_fill(", fixed = TRUE
  )
  expect_match(
    conditional, "np_continuous_kernel_scaled_restore(", fixed = TRUE
  )
  expect_false(grepl(
    "alloc_matd(num_obs_train_extern, num_obs_train_extern)",
    conditional, fixed = TRUE
  ))

  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(
    engine, "np_regression_prepared_bandwidth_copy(", fixed = TRUE
  )
  expect_match(
    engine, "np_regression_conditional_influence_finish(", fixed = TRUE
  )
  expect_match(engine, "double *conditional_kw = NULL;", fixed = TRUE)
  expect_match(engine, "double *conditional_pkw = NULL;", fixed = TRUE)
  expect_match(
    engine,
    "size_t pkw_plane = 0;",
    fixed = TRUE
  )
  expect_match(
    engine,
    "pkw_plane = np_jksum_size_mul_or_die(",
    fixed = TRUE
  )
  expect_match(engine, "(size_t)ii*pkw_plane +", fixed = TRUE)
  expect_false(grepl(
    "ii*num_obs_eval*num_xt + j*num_xt + i", engine, fixed = TRUE
  ))
  expect_match(
    engine,
    "if(prepared_status < 0) {\n      free_tmat(matrix_bandwidth);\n      error(",
    fixed = TRUE
  )
  expect_match(
    engine,
    "if(prepared_status == 0) {\n      const np_beta_bandwidth_prepare_status",
    fixed = TRUE
  )
})

test_that("canonical beta regression moments preserve the sidecar transcript", {
  root <- locate_beta_activation_sources()
  skip_if(is.null(root), "package sources unavailable")
  row_engine <- paste(
    readLines(file.path(root, "src", "continuous_kernel_row.c"),
              warn = FALSE),
    collapse = "\n"
  )

  owner_start <- regexpr(
    "np_continuous_kernel_beta_regression_moment_rows_validated(",
    row_engine,
    fixed = TRUE
  )[[1L]]
  influence_start <- regexpr(
    "np_continuous_kernel_beta_conditional_influence_stderr(",
    row_engine,
    fixed = TRUE
  )[[1L]]
  conditional_start <- regexpr(
    "np_continuous_kernel_beta_conditional_moment_rows_validated(",
    row_engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(owner_start, 0L)
  expect_gt(influence_start, 0L)
  expect_gt(conditional_start, 0L)
  expect_lt(owner_start, influence_start)
  expect_lt(influence_start, conditional_start)
  owner <- substr(row_engine, owner_start, influence_start - 1L)
  influence_owner <- substr(
    row_engine, influence_start, conditional_start - 1L
  )
  conditional_owner <- substr(
    row_engine, conditional_start, nchar(row_engine)
  )

  expect_match(
    owner, "np_continuous_kernel_beta_log_factor_row(", fixed = TRUE
  )
  expect_match(owner, "if(positive_weights)", fixed = TRUE)
  expect_match(owner, "const double new_total_weight", fixed = TRUE)
  expect_match(
    owner,
    "weighted_m2 += weight * delta *\n            (response[observation] - new_mean);",
    fixed = TRUE
  )
  expect_match(owner, "double weighted_response_sum = 0.0;", fixed = TRUE)
  expect_match(
    owner, "weighted_mean = weighted_response_sum / total_weight;",
    fixed = TRUE
  )
  expect_match(owner, "squared_weight_sum += weight * weight;", fixed = TRUE)
  expect_match(
    owner,
    "mean_stderr[evaluation] = sqrt(\n      (weighted_m2 / total_weight) *\n      (squared_weight_sum / (total_weight * total_weight)));",
    fixed = TRUE
  )
  expect_false(grepl(
    "NP_REGRESSION_STDERR_CONDITIONAL_INFLUENCE", owner, fixed = TRUE
  ))
  expect_match(
    influence_owner,
    "weight * (response[observation] - weighted_mean)",
    fixed = TRUE
  )
  expect_match(
    influence_owner,
    "fabs(total_weight) * sqrt((double)(variance_count - 1))",
    fixed = TRUE
  )
  expect_match(
    conditional_owner,
    "np_continuous_kernel_beta_conditional_influence_stderr(",
    fixed = TRUE
  )
  expect_false(grepl("result->row[observation]", owner, fixed = TRUE))
  expect_false(grepl("result->row[observation]", conditional_owner,
                     fixed = TRUE))
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
