static int np_beta_absolute_route_body(
  const NPBetaAbsoluteRouteCall *call,
  NPBetaAbsoluteRouteOwner *owner)
{
  const NPContinuousKernelRoute * const route = call->route;
  const int bandwidth_mode = call->bandwidth_mode;
  const int num_obs_train = call->num_obs_train;
  const int num_obs_eval = call->num_obs_eval;
  const int num_reg_unordered = call->num_reg_unordered;
  const int num_reg_ordered = call->num_reg_ordered;
  const int num_reg_continuous = call->num_reg_continuous;
  const int leave_one_out = call->leave_one_out;
  const int leave_one_out_offset = call->leave_one_out_offset;
  const int * const operator = call->operator;
  double ** const matrix_X_continuous_train =
    call->matrix_X_continuous_train;
  double ** const matrix_X_continuous_eval = call->matrix_X_continuous_eval;
  double ** const matrix_X_unordered_train = call->matrix_X_unordered_train;
  double ** const matrix_X_ordered_train = call->matrix_X_ordered_train;
  double ** const matrix_X_unordered_eval = call->matrix_X_unordered_eval;
  double ** const matrix_X_ordered_eval = call->matrix_X_ordered_eval;
  const int * const kernel_unordered = call->kernel_unordered;
  const int * const kernel_ordered = call->kernel_ordered;
  const double * const categorical_lambda = call->categorical_lambda;
  const int * const num_categories = call->num_categories;
  double ** const matrix_categorical_vals = call->matrix_categorical_vals;
  const int categorical_compress = call->categorical_compress;
  double * const vector_bandwidth = call->vector_bandwidth;
  double ** const matrix_bandwidth_eval = call->matrix_bandwidth_eval;
  double ** const matrix_bandwidth_train = call->matrix_bandwidth_train;
  double ** const matrix_Y = call->matrix_Y;
  double ** const matrix_W = call->matrix_W;
  const int ncol_Y = call->ncol_Y;
  const int ncol_W = call->ncol_W;
  double ** const matrix_Y_power2 = call->matrix_Y_power2;
  double ** const matrix_W_power2 = call->matrix_W_power2;
  const int ncol_Y_power2 = call->ncol_Y_power2;
  const int ncol_W_power2 = call->ncol_W_power2;
  const int kernel_power = call->kernel_power;
  double * const row = call->row;
  double * const weighted_sum = call->weighted_sum;
  double * const weighted_sum_power2 = call->weighted_sum_power2;
  const int retain_common_scale = call->retain_common_scale;
  double * const centered_m2 = call->centered_m2;
  double * const kw = call->kw;
  const NPBetaRegressionMomentCtx * const regression_moment_context =
    call->regression_moment_context;
  NPContinuousKernelDerivativeDiagnostics * const route_diagnostics =
    call->route_diagnostics;
  NPContinuousKernelProgressFunction progress = call->progress;
  NPContinuousKernelBetaPreparedContext * const beta_prepared =
    &owner->beta_prepared;
  NPContinuousKernelRowWorkspace * const workspace = &owner->workspace;
  NPContinuousKernelRowPlan plan;
  NPContinuousKernelRowResult row_result;
  NPContinuousKernelDerivativeAccumulator * const derivative_accumulator =
    &owner->derivative_accumulator;
  NPContinuousKernelLevelDerivativeWorkspace * const
    regression_gradient_workspace = &owner->regression_gradient_workspace;
  NPContinuousKernelLogFactorProvider categorical_provider;
  NPBetaScaledRowCategoricalContext * const categorical_context =
    &owner->categorical_context;
  const int has_categories =
    num_reg_unordered > 0 || num_reg_ordered > 0;
  int train_is_eval = num_obs_train == num_obs_eval;
  int status = KWSNP_ERR_BADINVOC;
  int evaluation;
  int coordinate;
  int derivative_coordinate = -1;
  int response_column;
  int weight_column;

  if(route_diagnostics != NULL) {
    route_diagnostics->bad_coordinate = -1;
    route_diagnostics->bad_observation = -1;
    route_diagnostics->undefined_count = 0;
    route_diagnostics->beta_status = NP_BETA_OK;
  }
  categorical_provider.function = np_beta_categorical_log_factor;
  categorical_provider.context = categorical_context;

  if(route == NULL ||
     np_continuous_kernel_route_validate(route, num_reg_continuous) !=
       NP_CKERNEL_ROUTE_OK ||
     route->segment_count != 1 ||
     route->segment[0].descriptor.family != NP_CKERNEL_FAMILY_BETA ||
     route->segment[0].coordinate_offset != 0 ||
     route->segment[0].coordinate_count != num_reg_continuous ||
     num_obs_train <= 0 || num_obs_eval <= 0 ||
     num_reg_unordered < 0 || num_reg_ordered < 0 ||
     num_reg_unordered > INT_MAX - num_reg_ordered ||
     num_reg_continuous <= 0 || operator == NULL ||
     matrix_X_continuous_train == NULL ||
     matrix_X_continuous_eval == NULL ||
     (bandwidth_mode != BW_FIXED && bandwidth_mode != BW_GEN_NN &&
      bandwidth_mode != BW_ADAP_NN) ||
     ncol_Y < 0 || ncol_W < 0 ||
     (ncol_Y > 0 && matrix_Y == NULL) ||
     (ncol_W > 0 && matrix_W == NULL) ||
     (weighted_sum == NULL && regression_moment_context == NULL) ||
     leave_one_out_offset < 0)
    return KWSNP_ERR_BADINVOC;
  if(has_categories &&
     ((num_reg_unordered > 0 &&
       (matrix_X_unordered_train == NULL ||
        matrix_X_unordered_eval == NULL || kernel_unordered == NULL)) ||
      (num_reg_ordered > 0 &&
       (matrix_X_ordered_train == NULL ||
        matrix_X_ordered_eval == NULL || kernel_ordered == NULL)) ||
      categorical_lambda == NULL || num_categories == NULL ||
      matrix_categorical_vals == NULL ||
      (categorical_compress != 0 && categorical_compress != 1)))
    return KWSNP_ERR_BADINVOC;

  for(response_column = 0; response_column < ncol_Y; ++response_column)
    if(matrix_Y[response_column] == NULL)
      return KWSNP_ERR_BADINVOC;
  for(weight_column = 0; weight_column < ncol_W; ++weight_column)
    if(matrix_W[weight_column] == NULL)
      return KWSNP_ERR_BADINVOC;

  for(coordinate = 0; coordinate < num_reg_continuous; ++coordinate) {
    if((operator[coordinate] != OP_NORMAL &&
        operator[coordinate] != OP_INTEGRAL &&
        operator[coordinate] != OP_CONVOLUTION &&
        operator[coordinate] != OP_DERIVATIVE) ||
       matrix_X_continuous_train[coordinate] == NULL ||
       matrix_X_continuous_eval[coordinate] == NULL)
      return KWSNP_ERR_BADINVOC;
    if(operator[coordinate] == OP_DERIVATIVE) {
      if(derivative_coordinate >= 0)
        return KWSNP_ERR_BADINVOC;
      derivative_coordinate = coordinate;
    }
    if(bandwidth_mode == BW_FIXED &&
       (vector_bandwidth == NULL ||
        !R_FINITE(vector_bandwidth[coordinate]) ||
        vector_bandwidth[coordinate] <= 0.0))
      return KWSNP_ERR_BADINVOC;
    if(bandwidth_mode == BW_GEN_NN &&
       (matrix_bandwidth_eval == NULL ||
        matrix_bandwidth_eval[coordinate] == NULL))
      return KWSNP_ERR_BADINVOC;
    if(bandwidth_mode == BW_ADAP_NN &&
       (matrix_bandwidth_train == NULL ||
        matrix_bandwidth_train[coordinate] == NULL))
      return KWSNP_ERR_BADINVOC;
    if(bandwidth_mode != BW_FIXED &&
       operator[coordinate] == OP_CONVOLUTION &&
       (matrix_bandwidth_eval == NULL ||
        matrix_bandwidth_train == NULL ||
        matrix_bandwidth_eval[coordinate] == NULL ||
        matrix_bandwidth_train[coordinate] == NULL))
      return KWSNP_ERR_BADINVOC;
    train_is_eval &=
      matrix_X_continuous_train[coordinate] ==
      matrix_X_continuous_eval[coordinate];
  }
  for(coordinate = 0; coordinate < num_reg_unordered; ++coordinate)
    train_is_eval &= matrix_X_unordered_train[coordinate] ==
      matrix_X_unordered_eval[coordinate];
  for(coordinate = 0; coordinate < num_reg_ordered; ++coordinate)
    train_is_eval &= matrix_X_ordered_train[coordinate] ==
      matrix_X_ordered_eval[coordinate];
  if(leave_one_out &&
     (!train_is_eval ||
      num_obs_eval > num_obs_train ||
      leave_one_out_offset > num_obs_train - num_obs_eval))
    return KWSNP_ERR_BADINVOC;
  if(derivative_coordinate < 0 && row == NULL &&
     weighted_sum_power2 == NULL && centered_m2 == NULL)
    return KWSNP_ERR_BADINVOC;
  if(retain_common_scale != 0 && retain_common_scale != 1)
    return KWSNP_ERR_BADINVOC;
  if(centered_m2 != NULL &&
     (derivative_coordinate >= 0 || kernel_power != 1 ||
      weighted_sum_power2 != NULL || kw != NULL ||
      ncol_Y != 0 || ncol_W != 0 ||
      regression_moment_context != NULL))
    return KWSNP_ERR_BADINVOC;
  if(regression_moment_context != NULL &&
     (regression_moment_context->response == NULL ||
      regression_moment_context->mean == NULL ||
      regression_moment_context->status == NULL ||
      (regression_moment_context->compute_errors != 0 &&
       regression_moment_context->compute_errors != 1) ||
      (regression_moment_context->compute_errors &&
       regression_moment_context->mean_stderr == NULL) ||
      (regression_moment_context->compute_gradient != 0 &&
       regression_moment_context->compute_gradient != 1) ||
      (regression_moment_context->compute_gradient &&
       regression_moment_context->gradient == NULL) ||
      (regression_moment_context->compute_gradient &&
       regression_moment_context->compute_errors &&
       regression_moment_context->hc0_scaled_residual == NULL &&
       regression_moment_context->gradient_stderr == NULL) ||
      (regression_moment_context->positive_weights != 0 &&
       regression_moment_context->positive_weights != 1) ||
      (regression_moment_context->preserve_mean != 0 &&
       regression_moment_context->preserve_mean != 1) ||
      (regression_moment_context->hc0_scaled_residual == NULL &&
       (regression_moment_context->preserve_mean != 0 ||
        regression_moment_context->hc0_residual_scale != 0.0)) ||
      (regression_moment_context->hc0_scaled_residual != NULL &&
       (!regression_moment_context->compute_errors ||
        !R_FINITE(regression_moment_context->hc0_residual_scale) ||
        regression_moment_context->hc0_residual_scale < 0.0)) ||
      derivative_coordinate >= 0 || kernel_power != 1 ||
      weighted_sum != NULL || weighted_sum_power2 != NULL ||
      centered_m2 != NULL || kw != NULL || ncol_Y != 0 || ncol_W != 0))
    return KWSNP_ERR_BADINVOC;

  if(bandwidth_mode == BW_FIXED) {
    if((size_t)num_reg_continuous > SIZE_MAX / sizeof(double *))
      return KWSNP_ERR_BADINVOC;
    owner->bandwidth_columns = (double **)malloc(
      (size_t)num_reg_continuous * sizeof(double *));
    if(owner->bandwidth_columns == NULL)
      return KWSNP_ERR_BADINVOC;
    for(coordinate = 0; coordinate < num_reg_continuous; ++coordinate)
      owner->bandwidth_columns[coordinate] = vector_bandwidth + coordinate;
  }

  plan.route = route;
  plan.bandwidth_mode = bandwidth_mode;
  plan.num_train = num_obs_train;
  plan.num_eval = num_obs_eval;
  plan.num_continuous = num_reg_continuous;
  plan.train_is_eval = train_is_eval;
  plan.train = matrix_X_continuous_train;
  plan.evaluation = matrix_X_continuous_eval;
  plan.bandwidth_eval = (bandwidth_mode == BW_FIXED) ?
    owner->bandwidth_columns : matrix_bandwidth_eval;
  plan.bandwidth_train = (bandwidth_mode == BW_FIXED) ?
    owner->bandwidth_columns : matrix_bandwidth_train;
  plan.operator = operator;
  plan.beta_prepared = beta_prepared;
  if(np_continuous_kernel_row_plan_validate(&plan) !=
     NP_CONTINUOUS_ROW_OK)
    goto cleanup;
  if(np_continuous_kernel_beta_prepared_context_prepare(
       beta_prepared, &plan) != NP_CONTINUOUS_ROW_OK)
    goto cleanup;
  row_result.row = row;

  if(has_categories) {
    NPContinuousKernelRowStatus row_status;

    row_status = np_beta_categorical_factor_context_prepare(
      categorical_context, num_obs_train, num_obs_eval,
      num_reg_unordered, num_reg_ordered,
      matrix_X_unordered_train, matrix_X_ordered_train,
      matrix_X_unordered_eval, matrix_X_ordered_eval,
      kernel_unordered, kernel_ordered,
      operator + num_reg_continuous, categorical_lambda,
      num_categories, matrix_categorical_vals, categorical_compress, row);
    if(row_status != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
  }

  if(regression_moment_context != NULL) {
    NPContinuousKernelRowStatus row_status =
      np_beta_regression_moment_rows_dispatch(
        &plan, leave_one_out, leave_one_out_offset,
        has_categories ? &categorical_provider : NULL,
        regression_moment_context->response,
        regression_moment_context->positive_weights,
        regression_moment_context->hc0_scaled_residual,
        regression_moment_context->hc0_residual_scale,
        regression_moment_context->preserve_mean,
        regression_moment_context->standard_error_mode,
        regression_moment_context->hc0_context,
        workspace, &row_result,
        regression_moment_context->mean,
        regression_moment_context->mean_stderr,
        route_diagnostics, regression_moment_context->progress);

    if(regression_moment_context->empty_rows != NULL && plan.num_eval == 1 &&
       row_status == NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT &&
       row_result.total_log_scale == -INFINITY) {
      /* Complete zero X factors, not a finite-scale signed cancellation. */
      np_regression_empty_row_mark(regression_moment_context->empty_rows, 0, 0);
      regression_moment_context->mean[0] = NA_REAL;
      if(regression_moment_context->mean_stderr != NULL)
        regression_moment_context->mean_stderr[0] = NA_REAL;
      if(regression_moment_context->compute_gradient)
        for(int coordinate = 0; coordinate < num_reg_continuous +
            num_reg_unordered + num_reg_ordered; ++coordinate) {
          regression_moment_context->gradient[coordinate][0] = NA_REAL;
          if(regression_moment_context->gradient_stderr != NULL)
            regression_moment_context->gradient_stderr[coordinate][0] = NA_REAL;
        }
      *regression_moment_context->status = NP_CONTINUOUS_ROW_OK;
      status = 0;
      goto cleanup;
    }
    *regression_moment_context->status = row_status;
    if(row_status != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
    if(regression_moment_context->compute_gradient &&
       (!regression_moment_context->preserve_mean ||
        regression_moment_context->gradient_stderr != NULL)) {
      int infinite_count = 0;
      int undefined_count = 0;
      /* Mean uncertainty can retain HC0 residuals independently of gradient
       * uncertainty. The gradient primitives require this channel only when
       * their uncertainty output is requested. */
      const double * const gradient_residual =
        regression_moment_context->gradient_stderr != NULL ?
        regression_moment_context->hc0_scaled_residual : NULL;

      row_status = np_beta_regression_gradient_rows_validated(
        &plan, has_categories ? &categorical_provider : NULL,
        regression_moment_context->response,
        gradient_residual,
        regression_moment_context->hc0_residual_scale,
        regression_moment_context->preserve_mean,
        regression_moment_context->standard_error_mode,
        regression_moment_context->hc0_context,
        regression_gradient_workspace,
        regression_moment_context->gradient,
        regression_moment_context->gradient_stderr,
        route_diagnostics, &infinite_count, &undefined_count);
      *regression_moment_context->status = row_status;
      if(row_status != NP_CONTINUOUS_ROW_OK)
        goto cleanup;
      if(!regression_moment_context->preserve_mean &&
         (infinite_count > 0 || undefined_count > 0))
        warning("beta regression gradient produced %d infinite endpoint value(s) and %d undefined cancellation(s)",
                infinite_count, undefined_count);
      if(has_categories) {
        row_status =
          np_beta_regression_categorical_gradients_validated(
            &plan, leave_one_out, leave_one_out_offset,
            num_reg_unordered, num_reg_ordered,
            categorical_context, regression_moment_context->response,
            gradient_residual,
            regression_moment_context->hc0_residual_scale,
            regression_moment_context->positive_weights,
            regression_moment_context->preserve_mean,
            regression_moment_context->standard_error_mode,
            regression_moment_context->ordinary_response,
            regression_moment_context->hc0_context,
            workspace, &row_result,
            regression_moment_context->mean,
            regression_moment_context->mean_stderr,
            regression_moment_context->gradient,
            regression_moment_context->gradient_stderr,
            route_diagnostics, regression_moment_context->empty_rows);
        *regression_moment_context->status = row_status;
        if(row_status != NP_CONTINUOUS_ROW_OK)
          goto cleanup;
      }
    }
    status = 0;
    goto cleanup;
  }

  if(centered_m2 != NULL) {
    const NPContinuousKernelRowStatus row_status =
      np_continuous_kernel_beta_centered_moment_rows_validated(
        &plan, leave_one_out, leave_one_out_offset,
        has_categories ? &categorical_provider : NULL,
        workspace, &row_result, weighted_sum, centered_m2,
        route_diagnostics, progress);

    if(row_status != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
    status = 0;
    goto cleanup;
  }

  if(weighted_sum_power2 != NULL || retain_common_scale) {
    NPContinuousKernelRowStatus row_status;

    if(derivative_coordinate >= 0 || kernel_power != 1 ||
       (!retain_common_scale && kw != NULL))
      goto cleanup;
    row_status =
      np_continuous_kernel_beta_dual_power_rows_validated(
        &plan, leave_one_out, leave_one_out_offset,
        has_categories ? &categorical_provider : NULL,
        matrix_Y, ncol_Y, matrix_W, ncol_W,
        matrix_Y_power2, ncol_Y_power2, matrix_W_power2, ncol_W_power2,
        workspace, &row_result, weighted_sum, weighted_sum_power2,
        call->power2_observation_scale,
        retain_common_scale, retain_common_scale ? kw : NULL,
        route_diagnostics, progress);

    if(row_status != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
    if(call->support_rank != NULL) {
      /* Only the one-row common-scaled LP owner supplies this request. */
      *call->support_rank = np_lp_rank_upper_bound_from_design(
        row_result.row, num_obs_train, call->support_cap, call->support_design);
    }
    status = 0;
    goto cleanup;
  }

  if(derivative_coordinate >= 0) {
    NPContinuousKernelRowStatus row_status;

    if(kernel_power != 1) {
      size_t derivative_row_count;

      if((size_t)num_obs_train > SIZE_MAX / 2 ||
         !np_size_array_bytes_checked(
           (size_t)num_obs_train * 2, sizeof(double),
           &derivative_row_count))
        goto cleanup;
      owner->derivative_row_storage =
        (double *)malloc(derivative_row_count);
      if(has_categories) {
        owner->derivative_factor_log_absolute = (double *)malloc(
          (size_t)num_obs_train * sizeof(double));
        owner->derivative_factor_sign = (signed char *)malloc(
          (size_t)num_obs_train * sizeof(signed char));
      }
      if(owner->derivative_row_storage == NULL ||
         (has_categories &&
          (owner->derivative_factor_log_absolute == NULL ||
           owner->derivative_factor_sign == NULL)))
        goto cleanup;
    }

    if(has_categories) {
      row_status = kernel_power == 1 ?
        np_continuous_kernel_beta_derivative_absolute_rows_with_log_factor_validated(
          &plan, leave_one_out, leave_one_out_offset,
          derivative_coordinate, &categorical_provider,
          matrix_Y, ncol_Y, matrix_W, ncol_W, derivative_accumulator,
          workspace,
          weighted_sum, kw, route_diagnostics) :
        np_continuous_kernel_beta_derivative_powered_rows_validated(
          &plan, leave_one_out, leave_one_out_offset,
          derivative_coordinate, kernel_power, &categorical_provider,
          matrix_Y, ncol_Y, matrix_W, ncol_W,
          derivative_accumulator, workspace, owner->derivative_row_storage,
          owner->derivative_factor_log_absolute,
          owner->derivative_factor_sign,
          weighted_sum, kw, route_diagnostics);
    } else {
      row_status = kernel_power == 1 ?
        np_continuous_kernel_beta_derivative_absolute_rows_validated(
          &plan, leave_one_out, leave_one_out_offset,
          derivative_coordinate, matrix_Y, ncol_Y, matrix_W, ncol_W,
          derivative_accumulator, weighted_sum, kw, route_diagnostics) :
        np_continuous_kernel_beta_derivative_powered_rows_validated(
          &plan, leave_one_out, leave_one_out_offset,
          derivative_coordinate, kernel_power, NULL,
          matrix_Y, ncol_Y, matrix_W, ncol_W,
          derivative_accumulator, workspace, owner->derivative_row_storage,
          NULL, NULL, weighted_sum, kw, route_diagnostics);
    }

    if(row_status != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
    status = 0;
    goto cleanup;
  }

  for(evaluation = 0; evaluation < num_obs_eval; ++evaluation) {
    const int omitted_observation = call->fold != NULL ?
      call->fold->eval_to_train[evaluation] :
      (leave_one_out ? evaluation + leave_one_out_offset : -1);
    NPContinuousKernelRowStatus row_status;
    int observation;

    if(call->fold != NULL && bandwidth_mode == BW_ADAP_NN &&
       np_nn_adaptive_fold_select_row(num_obs_train, num_reg_continuous,
         call->fold->train, call->fold->primary, call->fold->successor,
         call->fold->scale, omitted_observation, call->fold->selected) !=
       NP_NN_GEOMETRY_OK)
      goto cleanup;

    if(has_categories)
      row_status =
        np_continuous_kernel_beta_factor_row_with_log_factor(
          &plan, evaluation, omitted_observation, &categorical_provider,
          workspace, &row_result);
    else
      row_status = np_continuous_kernel_beta_factor_row(
        &plan, evaluation, omitted_observation, workspace, &row_result);

    if(row_status != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
    if(ncol_Y == 0 && ncol_W == 0) {
      double sum = 0.0;

      for(observation = 0; observation < num_obs_train; ++observation) {
        double value;

        if((kernel_power == 1 ?
            np_continuous_kernel_signed_log_restore(
              workspace->primary_log_absolute[observation],
              workspace->primary_sign[observation], &value) :
            np_continuous_kernel_signed_log_power_restore(
              workspace->primary_log_absolute[observation],
              workspace->primary_sign[observation], kernel_power,
              &value)) != NP_CONTINUOUS_ROW_OK)
          goto cleanup;
        if(kw != NULL) {
          double raw_value;

          if(kernel_power == 1)
            raw_value = value;
          else if(np_continuous_kernel_signed_log_restore(
                    workspace->primary_log_absolute[observation],
                    workspace->primary_sign[observation],
                    &raw_value) != NP_CONTINUOUS_ROW_OK)
            goto cleanup;
          kw[(size_t)evaluation * (size_t)num_obs_train +
             (size_t)observation] = raw_value;
        }
        sum += value;
        if(!R_FINITE(sum))
          goto cleanup;
      }
      weighted_sum[evaluation] = sum;
    } else {
      const int response_extent = (ncol_Y > 0) ? ncol_Y : 1;
      const int weight_extent = (ncol_W > 0) ? ncol_W : 1;
      size_t sum_extent;
      size_t output_base;
      size_t output_offset;

      if((size_t)response_extent > SIZE_MAX / (size_t)weight_extent)
        goto cleanup;
      sum_extent = (size_t)response_extent * (size_t)weight_extent;
      if((size_t)evaluation > SIZE_MAX / sum_extent)
        goto cleanup;
      output_base = (size_t)evaluation * sum_extent;
      for(output_offset = 0; output_offset < sum_extent; ++output_offset)
        weighted_sum[output_base + output_offset] = 0.0;

      for(observation = 0; observation < num_obs_train; ++observation) {
        double value;

        if((kernel_power == 1 ?
            np_continuous_kernel_signed_log_restore(
              workspace->primary_log_absolute[observation],
              workspace->primary_sign[observation], &value) :
            np_continuous_kernel_signed_log_power_restore(
              workspace->primary_log_absolute[observation],
              workspace->primary_sign[observation], kernel_power,
              &value)) != NP_CONTINUOUS_ROW_OK)
          goto cleanup;
        if(kw != NULL) {
          double raw_value;

          if(kernel_power == 1)
            raw_value = value;
          else if(np_continuous_kernel_signed_log_restore(
                    workspace->primary_log_absolute[observation],
                    workspace->primary_sign[observation],
                    &raw_value) != NP_CONTINUOUS_ROW_OK)
            goto cleanup;
          kw[(size_t)evaluation * (size_t)num_obs_train +
             (size_t)observation] = raw_value;
        }
        if(observation == omitted_observation)
          continue;

        for(response_column = 0; response_column < response_extent;
            ++response_column) {
          const double response_value = (ncol_Y > 0) ?
            matrix_Y[response_column][observation] : 1.0;
          for(weight_column = 0; weight_column < weight_extent;
              ++weight_column) {
            const double weight_value = (ncol_W > 0) ?
              matrix_W[weight_column][observation] : 1.0;
            output_offset = output_base +
              (size_t)response_column * (size_t)weight_extent +
              (size_t)weight_column;
            weighted_sum[output_offset] +=
              value * response_value * weight_value;
            if(!R_FINITE(weighted_sum[output_offset]))
              goto cleanup;
          }
        }
      }
    }
    if(progress != NULL)
      progress(evaluation + 1, num_obs_eval);
    if((evaluation & 31) == 0)
      R_CheckUserInterrupt();
  }

  status = 0;

cleanup:
  return status;
}
