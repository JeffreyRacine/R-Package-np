static NPContinuousKernelRowStatus
np_beta_regression_categorical_gradients_validated(
  const NPContinuousKernelRowPlan *plan,
  int leave_one_out,
  int leave_one_out_offset,
  int nunordered,
  int nordered,
  NPBetaScaledRowCategoricalContext *context,
  const double *response,
  const double *hc0_scaled_residual,
  double hc0_residual_scale,
  int positive_weights,
  int preserve_gradient,
  NPRegressionStandardErrorMode standard_error_mode,
  int ordinary_response,
  const NPRegressionHC0Context *hc0_context,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *row_result,
  const double *mean,
  const double *mean_stderr,
  double **gradient,
  double **gradient_stderr,
  NPContinuousKernelDerivativeDiagnostics *diagnostics,
  NPRegressionLPEmptyRows *empty_rows)
{
  NPContinuousKernelRowStatus status = NP_CONTINUOUS_ROW_OK;
  const int do_gerr = gradient_stderr != NULL;
  const int ordinary_hc0 = hc0_scaled_residual != NULL;
  const int paired_influence = do_gerr && standard_error_mode ==
    NP_REGRESSION_STDERR_CONDITIONAL_INFLUENCE;
  NPContinuousKernelLogFactorProvider provider;
  NPContinuousKernelLogFactorProvider level_provider;
  NPBetaScaledRowCategoricalContext level_context;
  double **alternate_unordered = NULL;
  double **alternate_ordered = NULL;
  double * const *eval_unordered;
  double * const *eval_ordered;
  const int *num_categories;
  double * const *category_values;
  double *numeric_storage;
  double *alternate_value;
  double *alternate_mean;
  double *alternate_stderr;
  double *level_coefficient;
  signed char *direction;
  size_t pointer_bytes;
  size_t numeric_count;
  int coordinate;
  int evaluation;

  if(plan == NULL || plan->num_eval <= 0 || plan->num_continuous <= 0 ||
     nunordered < 0 || nordered < 0 ||
     nunordered > INT_MAX - nordered ||
     (nunordered + nordered) <= 0 ||
     context == NULL || response == NULL || workspace == NULL ||
     row_result == NULL || mean == NULL || gradient == NULL ||
     (ordinary_response != 0 && ordinary_response != 1) ||
     (ordinary_response && do_gerr && !ordinary_hc0) ||
     (paired_influence && plan->num_eval != 1) ||
     (preserve_gradient != 0 && preserve_gradient != 1) ||
     (ordinary_hc0 &&
      (!do_gerr || !R_FINITE(hc0_residual_scale) ||
       hc0_residual_scale < 0.0)) ||
     (do_gerr && mean_stderr == NULL) ||
     context->dense_spec.ntrain != plan->num_train ||
     context->dense_spec.neval != plan->num_eval ||
     context->dense_spec.nunordered != nunordered ||
     context->dense_spec.nordered != nordered ||
     context->dense_spec.num_categories == NULL ||
     context->dense_spec.category_values == NULL ||
     (nunordered > 0 && context->dense_spec.eval_unordered == NULL) ||
     (nordered > 0 && context->dense_spec.eval_ordered == NULL) ||
     !np_size_mul_checked((size_t)plan->num_eval, 3U, &numeric_count) ||
     !np_size_add_checked(
       numeric_count, (size_t)plan->num_train, &numeric_count))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;

  eval_unordered = context->dense_spec.eval_unordered;
  eval_ordered = context->dense_spec.eval_ordered;
  num_categories = context->dense_spec.num_categories;
  category_values = context->dense_spec.category_values;
  if(nunordered > 0) {
    if(!np_size_array_bytes_checked((size_t)nunordered,
                                    sizeof(double *), &pointer_bytes))
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
    alternate_unordered = (double **)R_alloc(
      (size_t)nunordered, (int)sizeof(double *));
    memcpy(alternate_unordered, eval_unordered, pointer_bytes);
  }
  if(nordered > 0) {
    if(!np_size_array_bytes_checked((size_t)nordered,
                                    sizeof(double *), &pointer_bytes))
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
    alternate_ordered = (double **)R_alloc(
      (size_t)nordered, (int)sizeof(double *));
    memcpy(alternate_ordered, eval_ordered, pointer_bytes);
  }
  numeric_storage = (double *)R_alloc(numeric_count, (int)sizeof(double));
  direction = (signed char *)R_alloc(
    (size_t)plan->num_eval, (int)sizeof(signed char));
  alternate_value = numeric_storage;
  alternate_mean = numeric_storage + plan->num_eval;
  alternate_stderr = numeric_storage + 2 * plan->num_eval;
  level_coefficient = numeric_storage + 3 * plan->num_eval;
  provider.function = np_beta_categorical_log_factor;
  provider.context = context;
  level_context = *context;
  level_provider.function = np_beta_categorical_log_factor;
  level_provider.context = &level_context;

  if(paired_influence) {
    status = np_beta_conditional_categorical_influence_row(
      workspace, row_result, response, plan->num_train,
      leave_one_out ? leave_one_out_offset : -1, mean[0],
      level_coefficient, NULL);
    if(status != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
  }

  for(coordinate = 0; coordinate < nunordered; ++coordinate) {
    const int output_coordinate = plan->num_continuous + coordinate;

    if(num_categories[coordinate] <= 0 ||
       category_values[coordinate] == NULL) {
      status = NP_CONTINUOUS_ROW_ERR_LAYOUT;
      goto cleanup;
    }
    for(evaluation = 0; evaluation < plan->num_eval; ++evaluation)
      alternate_value[evaluation] = category_values[coordinate][0];
    alternate_unordered[coordinate] = alternate_value;
    context->dense_spec.eval_unordered = alternate_unordered;
    if(context->use_compressed)
      context->compressed_spec.eval_unordered = alternate_unordered;
    if(ordinary_response) {
      status =
        np_continuous_kernel_beta_regression_paired_hc0_rows_validated(
          plan, leave_one_out, leave_one_out_offset,
          &level_provider, &provider, response, hc0_scaled_residual,
          hc0_residual_scale, hc0_context, output_coordinate,
          eval_unordered[coordinate], alternate_value, workspace, row_result,
          level_coefficient, (size_t)plan->num_train,
          preserve_gradient ? NULL : gradient[output_coordinate],
          do_gerr ? gradient_stderr[output_coordinate] : NULL, diagnostics);
    } else {
      status = np_beta_regression_moment_rows_dispatch(
        plan, leave_one_out, leave_one_out_offset, &provider,
        response, positive_weights, NULL, 0.0, 0, standard_error_mode,
        NULL, workspace, row_result,
        alternate_mean, do_gerr ? alternate_stderr : NULL, diagnostics, NULL);
    }
    if(status == NP_CONTINUOUS_ROW_OK && paired_influence)
      status = np_beta_conditional_categorical_influence_row(
        workspace, row_result, response, plan->num_train,
        leave_one_out ? leave_one_out_offset : -1, alternate_mean[0],
        level_coefficient, gradient_stderr[output_coordinate]);
    context->dense_spec.eval_unordered = eval_unordered;
    if(context->use_compressed)
      context->compressed_spec.eval_unordered = eval_unordered;
    alternate_unordered[coordinate] = eval_unordered[coordinate];
    if(empty_rows != NULL && plan->num_eval == 1 &&
       status == NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT &&
       row_result->total_log_scale == -INFINITY) {
      np_regression_empty_row_mark(empty_rows, 0, 1);
      if(empty_rows->component_flags != NULL)
        empty_rows->component_flags[output_coordinate] = 1;
      gradient[output_coordinate][0] = NA_REAL;
      if(do_gerr) gradient_stderr[output_coordinate][0] = NA_REAL;
      status = NP_CONTINUOUS_ROW_OK;
      continue;
    }
    if(status != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
    for(evaluation = 0; evaluation < plan->num_eval; ++evaluation) {
      if(!preserve_gradient && !ordinary_response)
        gradient[output_coordinate][evaluation] =
          mean[evaluation] - alternate_mean[evaluation];
      if(do_gerr && !ordinary_hc0 && !paired_influence)
        gradient_stderr[output_coordinate][evaluation] = sqrt(
          mean_stderr[evaluation] * mean_stderr[evaluation] +
          alternate_stderr[evaluation] * alternate_stderr[evaluation]);
    }
  }

  for(coordinate = 0; coordinate < nordered; ++coordinate) {
    const int category_coordinate = nunordered + coordinate;
    const int output_coordinate =
      plan->num_continuous + nunordered + coordinate;
    const int category_count = num_categories[category_coordinate];
    const double * const support = category_values[category_coordinate];

    if(category_count <= 0 || support == NULL) {
      status = NP_CONTINUOUS_ROW_ERR_LAYOUT;
      goto cleanup;
    }
    if(category_count == 1) {
      for(evaluation = 0; evaluation < plan->num_eval; ++evaluation) {
        if(!preserve_gradient)
          gradient[output_coordinate][evaluation] = 0.0;
        if(do_gerr)
          gradient_stderr[output_coordinate][evaluation] = 0.0;
      }
      continue;
    }
    for(evaluation = 0; evaluation < plan->num_eval; ++evaluation) {
      const double value = eval_ordered[coordinate][evaluation];
      int support_index = 0;

      while(support_index < category_count &&
            support[support_index] != value)
        ++support_index;
      if(support_index == category_count) {
        status = NP_CONTINUOUS_ROW_ERR_LAYOUT;
        goto cleanup;
      }
      if(support_index == 0) {
        alternate_value[evaluation] = support[1];
        direction[evaluation] = -1;
      } else {
        alternate_value[evaluation] = support[support_index - 1];
        direction[evaluation] = 1;
      }
    }
    alternate_ordered[coordinate] = alternate_value;
    context->dense_spec.eval_ordered = alternate_ordered;
    if(context->use_compressed)
      context->compressed_spec.eval_ordered = alternate_ordered;
    if(ordinary_response) {
      status =
        np_continuous_kernel_beta_regression_paired_hc0_rows_validated(
          plan, leave_one_out, leave_one_out_offset,
          &level_provider, &provider, response, hc0_scaled_residual,
          hc0_residual_scale, hc0_context, output_coordinate,
          eval_ordered[coordinate], alternate_value, workspace, row_result,
          level_coefficient, (size_t)plan->num_train,
          preserve_gradient ? NULL : gradient[output_coordinate],
          do_gerr ? gradient_stderr[output_coordinate] : NULL, diagnostics);
    } else {
      status = np_beta_regression_moment_rows_dispatch(
        plan, leave_one_out, leave_one_out_offset, &provider,
        response, positive_weights, NULL, 0.0, 0, standard_error_mode,
        NULL, workspace, row_result,
        alternate_mean, do_gerr ? alternate_stderr : NULL, diagnostics, NULL);
    }
    if(status == NP_CONTINUOUS_ROW_OK && paired_influence)
      status = np_beta_conditional_categorical_influence_row(
        workspace, row_result, response, plan->num_train,
        leave_one_out ? leave_one_out_offset : -1, alternate_mean[0],
        level_coefficient, gradient_stderr[output_coordinate]);
    context->dense_spec.eval_ordered = eval_ordered;
    if(context->use_compressed)
      context->compressed_spec.eval_ordered = eval_ordered;
    alternate_ordered[coordinate] = eval_ordered[coordinate];
    if(empty_rows != NULL && plan->num_eval == 1 &&
       status == NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT &&
       row_result->total_log_scale == -INFINITY) {
      np_regression_empty_row_mark(empty_rows, 0, 1);
      if(empty_rows->component_flags != NULL)
        empty_rows->component_flags[output_coordinate] = 1;
      gradient[output_coordinate][0] = NA_REAL;
      if(do_gerr) gradient_stderr[output_coordinate][0] = NA_REAL;
      status = NP_CONTINUOUS_ROW_OK;
      continue;
    }
    if(status != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
    for(evaluation = 0; evaluation < plan->num_eval; ++evaluation) {
      if(!preserve_gradient)
        gradient[output_coordinate][evaluation] =
          (double)direction[evaluation] *
          (ordinary_response ? gradient[output_coordinate][evaluation] :
           mean[evaluation] - alternate_mean[evaluation]);
      if(do_gerr && !ordinary_hc0 && !paired_influence)
        gradient_stderr[output_coordinate][evaluation] = sqrt(
          mean_stderr[evaluation] * mean_stderr[evaluation] +
          alternate_stderr[evaluation] * alternate_stderr[evaluation]);
    }
  }

cleanup:
  return status;
}
