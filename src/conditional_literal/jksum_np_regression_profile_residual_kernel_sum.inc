static int np_regression_profile_residual_kernel_sum(
  const NPRegressionCategoricalProfileFitCall *call,
  int ntrain, int neval, double **train_u, double **train_o,
  double **eval_u, double **eval_o, double **response_columns,
  int response_column_count, double *weighted_sum,
  NPRegressionResidualKernelContext *residual_context)
{
  NP_DualPowerCtx preparation = {0};
  preparation.residual_preparation = residual_context;
  return kernel_weighted_sum_np_ctx_ex(
    call->kernel_c, call->kernel_u, call->kernel_o,
    call->bandwidth_mode, ntrain, neval,
    call->num_reg_unordered, call->num_reg_ordered, 0,
    0, 0, 1, 1, 0, 0, 0, 0, 0,
    call->operator, OP_NOOP, 0, 0, NULL, 1,
    response_column_count, 0, NP_TREE_FALSE, 0, NULL, NULL, NULL, NULL,
    train_u, train_o, NULL, eval_u, eval_o, NULL,
    response_columns, NULL, NULL, call->vector_scale_factor, 1,
    NULL, NULL, call->lambda, call->num_categories,
    call->matrix_categorical_vals, NULL, weighted_sum, NULL, NULL,
    NULL, &preparation, NULL, NULL, NULL, 0, NULL,
    NULL);
}
