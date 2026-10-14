NPContinuousKernelRowStatus np_beta_scaled_row_context_prepare(
  NPBetaScaledRowContext *context,
  const NPContinuousKernelRoute *route,
  NPContinuousKernelDerivativeDiagnostics *diagnostics,
  int bandwidth_mode,
  int ntrain,
  int neval,
  int ncontinuous,
  int nunordered,
  int nordered,
  double **train_continuous,
  double **eval_continuous,
  double **train_unordered,
  double **eval_unordered,
  double **train_ordered,
  double **eval_ordered,
  double **bandwidth_train,
  double **bandwidth_eval,
  const int *operator_code,
  const int *kernel_unordered,
  const int *kernel_ordered,
  const double *lambda,
  const int *num_categories,
  double **category_values,
  int categorical_compress,
  double *row)
{
  NPContinuousKernelRowStatus status;
  int coordinate;
  int train_is_eval = ntrain == neval;

  if(context == NULL || route == NULL || diagnostics == NULL ||
     ntrain <= 0 || neval <= 0 || ncontinuous <= 0 ||
     train_continuous == NULL || eval_continuous == NULL ||
     bandwidth_train == NULL || bandwidth_eval == NULL ||
     operator_code == NULL || row == NULL ||
     (categorical_compress != 0 && categorical_compress != 1))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(np_continuous_kernel_route_validate(route, ncontinuous) !=
       NP_CKERNEL_ROUTE_OK ||
     !np_continuous_kernel_route_has_beta(route) ||
     route->segment_count != 1 ||
     route->segment[0].coordinate_offset != 0 ||
     route->segment[0].coordinate_count != ncontinuous)
    return NP_CONTINUOUS_ROW_ERR_ROUTE;

  for(coordinate = 0; coordinate < ncontinuous; ++coordinate)
    train_is_eval &= train_continuous[coordinate] ==
      eval_continuous[coordinate];
  for(coordinate = 0; coordinate < nunordered; ++coordinate)
    train_is_eval &= train_unordered[coordinate] ==
      eval_unordered[coordinate];
  for(coordinate = 0; coordinate < nordered; ++coordinate)
    train_is_eval &= train_ordered[coordinate] == eval_ordered[coordinate];

  context->plan.route = route;
  context->plan.bandwidth_mode = bandwidth_mode;
  context->plan.num_train = ntrain;
  context->plan.num_eval = neval;
  context->plan.num_continuous = ncontinuous;
  context->plan.train_is_eval = train_is_eval;
  context->plan.train = train_continuous;
  context->plan.evaluation = eval_continuous;
  context->plan.bandwidth_train = bandwidth_train;
  context->plan.bandwidth_eval = bandwidth_eval;
  context->plan.operator = operator_code;
  context->plan.beta_prepared = &context->beta_prepared;
  context->row_result.row = row;
  context->diagnostics = diagnostics;
  context->has_categories = nunordered > 0 || nordered > 0;
  context->categorical_provider.function = np_beta_categorical_log_factor;
  context->categorical_provider.context = &context->categorical_context;

  status = np_continuous_kernel_row_plan_validate(&context->plan);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  status = np_continuous_kernel_beta_prepared_context_prepare(
    &context->beta_prepared, &context->plan);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  if(context->has_categories) {
    status = np_beta_categorical_factor_context_prepare(
      &context->categorical_context, ntrain, neval, nunordered, nordered,
      train_unordered, train_ordered, eval_unordered, eval_ordered,
      kernel_unordered, kernel_ordered, operator_code + ncontinuous,
      lambda, num_categories, category_values, categorical_compress, row);
    if(status != NP_CONTINUOUS_ROW_OK)
      return status;
  }
  context->ready = 1;
  return NP_CONTINUOUS_ROW_OK;
}
