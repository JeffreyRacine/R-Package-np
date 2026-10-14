static NP_ALWAYS_INLINE int np_regression_general_lp_empty_row(
  const int prepared_ann,
  const NPRegressionGeneralLPFitCall *call,
  NPRegressionGeneralLPFitOwner *owner,
  int row, int categorical, int moment_stride, int divide_weights,
  const double *computed_weights)
{
  double row_sum = 0.0;
  int status;
  double **bandwidth = categorical ? call->categorical_matrix_bandwidth :
    call->matrix_bandwidth;
  if(call->empty_rows == NULL ||
     !np_lp_failed_system_is_finite(&owner->solve_workspace, owner->nterms, 1))
    return 0;
  if(computed_weights == NULL) {
    if(owner->retained_kernel_row == NULL)
      owner->retained_kernel_row = (double *)malloc(
        (size_t)call->num_obs_train*sizeof(double));
    if(owner->retained_kernel_row == NULL)
      return 0;
    if(prepared_ann) {
      status = np_ann_fit_moment_row(call, owner, moment_stride, owner->moments,
        NULL, owner->retained_kernel_row, NULL);
    } else if(call->kernel_route != NULL) {
      /* Preserve the canonical response/basis scale used by the failed row. */
      status = np_beta_regression_lp_moment_row_canonical(
        call->bandwidth_mode, call->num_obs_train, call->num_reg_unordered,
        call->num_reg_ordered, call->num_reg_continuous, owner->nterms,
        moment_stride, call->operator, call->kernel_u, call->kernel_o,
        call->matrix_X_unordered_train, call->matrix_X_ordered_train,
        call->matrix_X_continuous_train, owner->eval_unordered,
        owner->eval_ordered, owner->eval_continuous, owner->response_columns,
        owner->basis, call->vector_scale_factor, bandwidth,
        call->bandwidth_mode == BW_GEN_NN ? owner->matrix_bandwidth_eval : bandwidth,
        call->lambda, call->num_categories, call->matrix_categorical_vals,
        call->categorical_compress, owner->moments, NULL, NULL,
        owner->retained_kernel_row, call->kernel_route,
        call->kernel_route_diagnostics,
        NULL);
    } else {
      status = kernel_weighted_sum_np_ctx_ex(
        call->kernel_c, call->kernel_u, call->kernel_o, call->bandwidth_mode,
        call->num_obs_train, 1, call->num_reg_unordered, call->num_reg_ordered,
        call->num_reg_continuous, 0, 0, 1, 1, divide_weights, 0, 0, 0, 0,
        call->operator, OP_NOOP, 0, 0, NULL, 1, 0, 0,
        call->bandwidth_mode == BW_ADAP_NN ? NP_TREE_FALSE : call->tree_enabled,
        0, call->bandwidth_mode == BW_ADAP_NN ? NULL : kdt_extern_X,
        NULL, NULL, NULL, call->matrix_X_unordered_train,
        call->matrix_X_ordered_train, call->matrix_X_continuous_train,
        owner->eval_unordered, owner->eval_ordered, owner->eval_continuous,
        NULL, NULL, NULL, call->vector_scale_factor, 1, bandwidth,
        call->bandwidth_mode == BW_GEN_NN ? owner->matrix_bandwidth_eval : bandwidth,
        call->lambda, call->num_categories, call->matrix_categorical_vals,
        NULL, &row_sum, NULL, owner->retained_kernel_row,
        call->gate_context, NULL, NULL, NULL, NULL,
#ifdef MPI2
        0,
#endif
        NULL,
        NULL);
    }
    if(status != 0)
      return 0;
    computed_weights = owner->retained_kernel_row;
  }
  if(!np_lp_complete_weights_are_zero(computed_weights, call->num_obs_train))
    return 0;
  np_regression_empty_row_mark(call->empty_rows, row, categorical);
  return 1;
}
