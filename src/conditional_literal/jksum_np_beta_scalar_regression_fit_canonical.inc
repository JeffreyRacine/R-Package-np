static NP_NOINLINE int np_beta_scalar_regression_fit_canonical(
  const int BANDWIDTH_reg,
  const int KERNEL_unordered_reg,
  const int KERNEL_ordered_reg,
  const int num_obs_train,
  const int num_obs_eval,
  const int num_reg_unordered,
  const int num_reg_ordered,
  const int num_reg_continuous,
  double **matrix_X_unordered_train,
  double **matrix_X_ordered_train,
  double **matrix_X_continuous_train,
  double **matrix_X_unordered_eval,
  double **matrix_X_ordered_eval,
  double **matrix_X_continuous_eval,
  const double *vector_Y,
  double *vector_scale_factor,
  int *num_categories,
  double **matrix_categorical_vals,
  double *mean,
  double **gradient,
  double *mean_stderr,
  double **gradient_stderr,
  const NPContinuousKernelRoute *kernel_route,
  NPContinuousKernelDerivativeDiagnostics *kernel_route_diagnostics,
  const int categorical_compress,
  const NPRegressionStandardErrorMode standard_error_mode,
  const int ordinary_response,
  const NPContinuousPreparedBandwidthView *prepared_bandwidth,
  const NPNNGeometryContext *nn_geometry_context,
  const NPRegressionHC0Context *hc0_context,
  NPRegressionFailure *failure,
  NPRegressionLPEmptyRows *empty_rows)
{
  const int num_categorical = num_reg_unordered + num_reg_ordered;
  const int num_predictors = num_reg_continuous + num_categorical;
  const int do_grad = gradient != NULL;
  const int do_gerr = gradient_stderr != NULL;
  const int hc0_residual_preparing = hc0_context != NULL &&
    hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_PREPARING;
  const int ordinary_hc0 = hc0_context != NULL &&
    !hc0_residual_preparing;
  const int bwmdim = (BANDWIDTH_reg == BW_GEN_NN) ? num_obs_eval :
    ((BANDWIDTH_reg == BW_ADAP_NN) ? num_obs_train : 1);
  int *operator;
  int *kernel_u = NULL;
  int *kernel_o = NULL;
  double *lambda = NULL;
  double **matrix_bandwidth;
  double *route_row;
  NPBetaRegressionMomentCtx regression_moment_context;
  NPContinuousKernelRowStatus regression_row_status =
    NP_CONTINUOUS_ROW_OK;
  int route_status;
  int coordinate;

  if(np_continuous_kernel_route_validate(kernel_route,
                                         num_reg_continuous) !=
       NP_CKERNEL_ROUTE_OK ||
     !np_continuous_kernel_route_has_beta(kernel_route) ||
     kernel_route->segment_count != 1 ||
     kernel_route->segment[0].coordinate_offset != 0 ||
     kernel_route->segment[0].coordinate_count != num_reg_continuous ||
     num_reg_continuous <= 0 ||
     (do_gerr && (!do_grad || mean_stderr == NULL)) ||
     kernel_route_diagnostics == NULL ||
     (categorical_compress != 0 && categorical_compress != 1) ||
     (standard_error_mode != NP_REGRESSION_STDERR_LOCAL_RESIDUAL &&
      standard_error_mode != NP_REGRESSION_STDERR_CONDITIONAL_INFLUENCE) ||
     (ordinary_hc0 &&
      (hc0_context->scaled_residual == NULL ||
       hc0_context->num_obs_train != num_obs_train ||
       !R_FINITE(hc0_context->residual_scale) ||
       hc0_context->residual_scale < 0.0)))
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_BETA, 1, "canonical beta regression route has an invalid layout");

  operator = (int *)R_alloc((size_t)num_predictors, sizeof(int));
  for(coordinate = 0; coordinate < num_predictors; ++coordinate)
    operator[coordinate] = OP_NORMAL;
  if(num_reg_unordered > 0) {
    kernel_u = (int *)R_alloc((size_t)num_reg_unordered, sizeof(int));
    for(coordinate = 0; coordinate < num_reg_unordered; ++coordinate)
      kernel_u[coordinate] = KERNEL_unordered_reg;
  }
  if(num_reg_ordered > 0) {
    kernel_o = (int *)R_alloc((size_t)num_reg_ordered, sizeof(int));
    for(coordinate = 0; coordinate < num_reg_ordered; ++coordinate)
      kernel_o[coordinate] = KERNEL_ordered_reg;
  }
  if(num_categorical > 0)
    lambda = (double *)R_alloc((size_t)num_categorical, sizeof(double));
  matrix_bandwidth = alloc_tmatd(bwmdim, num_reg_continuous);

  if(num_categorical > 0 &&
     kernel_bandwidth_mean(
       0, BW_FIXED, num_obs_train, num_obs_eval,
       0, 0, 0,
       num_reg_continuous, num_reg_unordered, num_reg_ordered,
       0, vector_scale_factor,
       NULL, NULL,
       matrix_X_continuous_train, matrix_X_continuous_eval,
       NULL, matrix_bandwidth, lambda) == 1) {
    free_tmat(matrix_bandwidth);
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_BETA, 2, "\n** Error: invalid categorical bandwidth.");
  }

  if(BANDWIDTH_reg != BW_FIXED) {
    const int prepared_status = np_regression_prepared_bandwidth_copy(
      prepared_bandwidth, BANDWIDTH_reg, num_obs_train, num_obs_eval,
      num_reg_continuous, matrix_bandwidth);

    if(prepared_status < 0) {
      free_tmat(matrix_bandwidth);
      NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_BETA, 3, "canonical beta regression route has an invalid prepared bandwidth view");
    }
    if(prepared_status == 0) {
      const np_beta_bandwidth_prepare_status bandwidth_status =
        np_beta_bandwidth_prepare_matrix_ctx(
          BANDWIDTH_reg == BW_GEN_NN ?
            NP_BETA_BANDWIDTH_GENERALIZED_NN :
            NP_BETA_BANDWIDTH_ADAPTIVE_NN,
          matrix_X_continuous_train, matrix_X_continuous_eval,
          vector_scale_factor,
          num_obs_train, num_obs_eval, num_reg_continuous,
          matrix_X_continuous_train == matrix_X_continuous_eval,
          BANDWIDTH_reg == BW_GEN_NN, BANDWIDTH_reg == BW_ADAP_NN,
          0, matrix_bandwidth, matrix_bandwidth, nn_geometry_context);

      if(bandwidth_status != NP_BETA_BANDWIDTH_PREPARE_OK) {
        free_tmat(matrix_bandwidth);
        if(bandwidth_status == NP_BETA_BANDWIDTH_PREPARE_ERR_ZERO_RADIUS) {
          const NPNNZeroRadiusInfo info = np_nn_zero_radius_info(
            BANDWIDTH_reg, num_obs_train, num_obs_eval, num_reg_continuous,
            matrix_X_continuous_train, matrix_X_continuous_eval,
            vector_scale_factor, nn_geometry_context, 0);
          if(failure != NULL) {
            np_regression_failure_zero_radius(failure, &info);
            return NP_REGRESSION_FIT_ERR_ZERO_NN_RADIUS;
          }
          np_nn_zero_radius_error(&info);
        }
        NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_BETA, 4, "\n** Error: %s.",
              np_beta_bandwidth_prepare_status_message(bandwidth_status));
      }
    }
  }

  route_row = (double *)np_jksum_malloc_array_or_die(
    (size_t)num_obs_train, sizeof(double), "beta regression row");
  regression_moment_context.response = vector_Y;
  regression_moment_context.hc0_scaled_residual = ordinary_hc0 ?
    hc0_context->scaled_residual : NULL;
  regression_moment_context.hc0_residual_scale = ordinary_hc0 ?
    hc0_context->residual_scale : 0.0;
  regression_moment_context.preserve_mean = ordinary_hc0 ?
    hc0_context->point_already_computed : 0;
  regression_moment_context.mean = mean;
  regression_moment_context.mean_stderr = mean_stderr;
  regression_moment_context.gradient = gradient;
  regression_moment_context.gradient_stderr = gradient_stderr;
  regression_moment_context.compute_errors = mean_stderr != NULL;
  regression_moment_context.compute_gradient = do_grad;
  regression_moment_context.positive_weights =
    kernel_route->segment[0].descriptor.order == 2;
  regression_moment_context.standard_error_mode = standard_error_mode;
  regression_moment_context.ordinary_response = ordinary_response;
  regression_moment_context.status = &regression_row_status;
  regression_moment_context.empty_rows =
    standard_error_mode == NP_REGRESSION_STDERR_CONDITIONAL_INFLUENCE &&
    num_obs_eval == 1 ? empty_rows : NULL;
  regression_moment_context.progress = np_progress_fit_loop_step;
  regression_moment_context.hc0_context = hc0_context;
  {
    const NPBetaAbsoluteRouteCall route_call = {
      .route = kernel_route,
      .bandwidth_mode = BANDWIDTH_reg,
      .num_obs_train = num_obs_train,
      .num_obs_eval = num_obs_eval,
      .num_reg_unordered = num_reg_unordered,
      .num_reg_ordered = num_reg_ordered,
      .num_reg_continuous = num_reg_continuous,
      .leave_one_out = 0,
      .leave_one_out_offset = 0,
      .operator = operator,
      .matrix_X_continuous_train = matrix_X_continuous_train,
      .matrix_X_continuous_eval = matrix_X_continuous_eval,
      .matrix_X_unordered_train = matrix_X_unordered_train,
      .matrix_X_ordered_train = matrix_X_ordered_train,
      .matrix_X_unordered_eval = matrix_X_unordered_eval,
      .matrix_X_ordered_eval = matrix_X_ordered_eval,
      .kernel_unordered = kernel_u,
      .kernel_ordered = kernel_o,
      .categorical_lambda = num_categorical > 0 ? lambda : NULL,
      .num_categories = num_categories,
      .matrix_categorical_vals = matrix_categorical_vals,
      .categorical_compress = categorical_compress,
      .vector_bandwidth = vector_scale_factor,
      .matrix_bandwidth_eval = BANDWIDTH_reg == BW_FIXED ?
        NULL : matrix_bandwidth,
      .matrix_bandwidth_train = BANDWIDTH_reg == BW_FIXED ?
        NULL : matrix_bandwidth,
      .owned_bandwidth_tmatrix = matrix_bandwidth,
      .matrix_Y = NULL,
      .matrix_W = NULL,
      .ncol_Y = 0,
      .ncol_W = 0,
      .matrix_Y_power2 = NULL,
      .matrix_W_power2 = NULL,
      .ncol_Y_power2 = 0,
      .ncol_W_power2 = 0,
      .kernel_power = 1,
      .row = route_row,
      .weighted_sum = NULL,
      .weighted_sum_power2 = NULL,
      .power2_observation_scale = NULL,
      .retain_common_scale = 0,
      .centered_m2 = NULL,
      .kw = NULL,
      .regression_moment_context = &regression_moment_context,
      .route_diagnostics = kernel_route_diagnostics,
      .progress = NULL
    };

    route_status = np_beta_absolute_route(&route_call);
  }

  if(route_status != 0) {
    if(kernel_route_diagnostics->beta_status != NP_BETA_OK)
      NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_BETA, 5, "canonical beta regression row failed in continuous dimension %d: %s",
            kernel_route_diagnostics->bad_coordinate + 1,
            np_beta_status_message(kernel_route_diagnostics->beta_status));
    if(regression_row_status == NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT)
      NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_BETA, 6, "all beta regression weights are zero at an evaluation point");
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_BETA, 7, "canonical beta regression row failed: %s",
          np_continuous_kernel_row_status_message(regression_row_status));
  }
  return 0;
}
