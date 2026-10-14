static NP_ALWAYS_INLINE int np_regression_general_lp_point_at_frame(
  const int prepared_ann,
  const NPRegressionGeneralLPFitCall *call,
  NPRegressionGeneralLPFitOwner *owner,
  const int moment_stride,
  const int response_y_offset,
  const int response_basis_offset,
  const double epsilon,
  const int row,
  const int empty_component,
  double *point,
  double *kernel_row,
  double *projection)
{
  NPLPSolvePolicyDiagnostics solve_diagnostics = {0, 0.0};
  NPLPSupportRankContext support_rank = {0, owner->nterms, owner->active_design};
  double pristine_anchor;
  int i, term;

  if(call == NULL || owner == NULL || point == NULL ||
     owner->nterms <= 0 || owner->moments == NULL ||
     owner->eval_basis == NULL)
    return 0;

  if(prepared_ann) {
    if(np_ann_fit_moment_row(call, owner, moment_stride, owner->moments,
         NULL, kernel_row, &support_rank) != 0) return 0;
  } else if(call->kernel_route != NULL) {
    if(np_beta_regression_lp_moment_row_canonical(
         call->bandwidth_mode,
         call->num_obs_train,
         call->num_reg_unordered,
         call->num_reg_ordered,
         call->num_reg_continuous,
         owner->nterms,
         moment_stride,
         call->operator,
         call->kernel_u,
         call->kernel_o,
         call->matrix_X_unordered_train,
         call->matrix_X_ordered_train,
         call->matrix_X_continuous_train,
         owner->eval_unordered,
         owner->eval_ordered,
         owner->eval_continuous,
         owner->response_columns,
         owner->basis,
         call->vector_scale_factor,
         call->categorical_matrix_bandwidth,
         call->bandwidth_mode == BW_GEN_NN ?
           owner->matrix_bandwidth_eval :
           call->categorical_matrix_bandwidth,
         call->lambda,
         call->num_categories,
         call->matrix_categorical_vals,
         call->categorical_compress,
         owner->moments,
         NULL,
         NULL,
         kernel_row,
         call->kernel_route,
         call->kernel_route_diagnostics,
         &support_rank) != 0)
      return 0;
  } else {
    if(kernel_weighted_sum_np_ctx_ex(
         call->kernel_c,
         call->kernel_u,
         call->kernel_o,
         call->bandwidth_mode,
         call->num_obs_train,
         1,
         call->num_reg_unordered,
         call->num_reg_ordered,
         call->num_reg_continuous,
         0, 0, 1, 1, kernel_row != NULL, 0, 0, 0, 0,
         call->operator,
         OP_NOOP,
         0, 0, NULL, 1,
         moment_stride,
         owner->nterms,
         call->bandwidth_mode == BW_ADAP_NN ?
           NP_TREE_FALSE : call->tree_enabled,
         0,
         call->bandwidth_mode == BW_ADAP_NN ? NULL : kdt_extern_X,
         NULL, NULL, NULL,
         call->matrix_X_unordered_train,
         call->matrix_X_ordered_train,
         call->matrix_X_continuous_train,
         owner->eval_unordered,
         owner->eval_ordered,
         owner->eval_continuous,
         owner->response_columns,
         owner->basis_columns,
         NULL,
         call->vector_scale_factor,
         1,
         call->categorical_matrix_bandwidth,
         call->bandwidth_mode == BW_GEN_NN ?
           owner->matrix_bandwidth_eval :
           call->categorical_matrix_bandwidth,
         call->lambda,
         call->num_categories,
         call->matrix_categorical_vals,
         NULL,
         owner->moments,
         NULL,
         kernel_row,
         call->gate_context,
         NULL,
         NULL, NULL, NULL,
#ifdef MPI2
         0,
#endif
         NULL,
         &support_rank) != 0)
      return 0;
  }

  for(i = 0; i < owner->nterms; ++i) {
    const int base = i*moment_stride;

    owner->solve_workspace.rhs_source[i] =
      owner->moments[base + response_y_offset];
    for(term = 0; term < owner->nterms; ++term)
      owner->solve_workspace.gram_source[i + term*owner->nterms] =
        owner->moments[base + term + response_basis_offset];
  }
  pristine_anchor = owner->solve_workspace.gram_source[0];
  if(np_lp_solve_workspace_solve_response_ranked(
       &owner->solve_workspace,
       owner->nterms,
       1,
       epsilon,
       support_rank.count,
       &solve_diagnostics) != NP_LP_SOLVE_POLICY_OK) {
    if(np_regression_general_lp_empty_row(prepared_ann, call, owner, row, 1,
         moment_stride, kernel_row != NULL, kernel_row)) {
      if(empty_component >= 0 && call->empty_rows->component_flags != NULL)
        call->empty_rows->component_flags[empty_component] = 1;
      *point = NA_REAL;
      return 1;
    }
    return 0;
  }

  *point = 0.0;
  for(i = 0; i < owner->nterms; ++i)
    *point += owner->eval_basis[i]*owner->solve_workspace.rhs_work[i];
  if(!R_FINITE(*point))
    return 0;

  if(projection != NULL) {
    for(i = 0; i < owner->nterms; ++i)
      owner->solve_workspace.rhs_source[i] = owner->eval_basis[i];
    if(np_lp_solve_workspace_solve_adjoint_factored(
         &owner->solve_workspace,
         owner->nterms,
         1,
         pristine_anchor,
         &solve_diagnostics) != NP_LP_SOLVE_POLICY_OK)
      return 0;
    memcpy(projection,
           owner->solve_workspace.rhs_work,
           (size_t)owner->nterms*sizeof(double));
  }
  return 1;
}
