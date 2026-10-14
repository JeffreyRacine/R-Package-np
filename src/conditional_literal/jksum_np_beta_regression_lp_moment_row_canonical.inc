static NP_NOINLINE int np_beta_regression_lp_moment_row_canonical(
  const int bandwidth_mode,
  const int num_obs_train,
  const int num_reg_unordered,
  const int num_reg_ordered,
  const int num_reg_continuous,
  const int nterms,
  const int outer_column_count,
  const int *operator,
  int *kernel_unordered,
  int *kernel_ordered,
  double **matrix_X_unordered_train,
  double **matrix_X_ordered_train,
  double **matrix_X_continuous_train,
  double **matrix_X_unordered_eval,
  double **matrix_X_ordered_eval,
  double **matrix_X_continuous_eval,
  double **outer_columns,
  double **basis,
  double *vector_scale_factor,
  double **matrix_bandwidth_train,
  double **matrix_bandwidth_eval,
  double *lambda,
  int *num_categories,
  double **matrix_categorical_vals,
  const int categorical_compress,
  double *weighted_sum,
  double *weighted_sum_power2,
  const double *power2_observation_scale,
  double *scaled_kernel_weights,
  const NPContinuousKernelRoute *kernel_route,
  NPContinuousKernelDerivativeDiagnostics *kernel_route_diagnostics,
  NPLPSupportRankContext * const support_rank)
{
  NP_DualPowerCtx dual_power_context = {
    weighted_sum_power2, 2, basis, basis, nterms, nterms, NULL, 1, NULL,
    power2_observation_scale, NULL, NULL, NULL, NULL
  };
  const NPContinuousKernelExecutionContext execution_context = {
    kernel_route, kernel_route_diagnostics, categorical_compress
  };

  return kernel_weighted_sum_np_ctx_ex(
    NULL, kernel_unordered, kernel_ordered,
    bandwidth_mode, num_obs_train, 1,
    num_reg_unordered, num_reg_ordered, num_reg_continuous,
    0, 0, 1, 0, 0,
    0, 0, 0, 0, operator, OP_NOOP,
    0, 0, NULL, 1, nterms, outer_column_count,
    NP_TREE_FALSE, 0, NULL, NULL, NULL, NULL,
    matrix_X_unordered_train, matrix_X_ordered_train,
    matrix_X_continuous_train,
    matrix_X_unordered_eval, matrix_X_ordered_eval,
    matrix_X_continuous_eval,
    basis, outer_columns, NULL, vector_scale_factor,
    bandwidth_mode == BW_FIXED ? 0 : 1,
    bandwidth_mode == BW_FIXED ? NULL : matrix_bandwidth_train,
    bandwidth_mode == BW_FIXED ? NULL : matrix_bandwidth_eval,
    lambda, num_categories, matrix_categorical_vals, NULL,
    weighted_sum, NULL, scaled_kernel_weights, NULL,
    &dual_power_context,
    NULL, &execution_context, NULL, NULL,
    support_rank);
}
