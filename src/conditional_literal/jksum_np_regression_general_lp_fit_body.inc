static NP_ALWAYS_INLINE SEXP np_regression_general_lp_fit_body(
  void *data, const int prepared_ann)
{
  NPRegressionGeneralLPFitExecution * const execution =
    (NPRegressionGeneralLPFitExecution *)data;
  const NPRegressionGeneralLPFitCall * const call = execution->call;
  NPRegressionGeneralLPFitOwner * const owner = &execution->owner;
  const int BANDWIDTH_reg = call->bandwidth_mode;
  const int num_obs_train = call->num_obs_train;
  const int num_obs_eval = call->num_obs_eval;
  const int num_reg_unordered = call->num_reg_unordered;
  const int num_reg_ordered = call->num_reg_ordered;
  const int num_reg_continuous = call->num_reg_continuous;
  const int use_bernstein = (int_glp_bernstein_extern != 0);
  const double epsilon = 1.0/(double)MAX(1, num_obs_train);
  const int hc0_residual_preparing = call->hc0_context != NULL &&
    call->hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_PREPARING;
  const int ordinary_hc0 = call->hc0_context != NULL &&
    !hc0_residual_preparing;
  const int preserve_point = ordinary_hc0 &&
    call->hc0_context->point_already_computed;
  const int categorical_pair = call->ordinary_response && call->do_grad &&
    (num_reg_unordered + num_reg_ordered > 0) &&
    !call->categorical_point_invariant;
  const int categorical_hc0 = ordinary_hc0 && call->do_gerr &&
    (num_reg_unordered + num_reg_ordered > 0) &&
    !call->categorical_point_invariant;
  const int hc0_information_row = hc0_residual_preparing || ordinary_hc0;
  /* Only conditional variance consumers retain the squared-moment matrix. */
  const int moment_errors = call->do_merr && !ordinary_hc0;
  NPRegressionResidualPreparation * const preparation =
    hc0_residual_preparing ? call->hc0_context->preparation : NULL;
  const int reuse_fit_kernel_row =
    (BANDWIDTH_reg == BW_FIXED) &&
    (call->tree_enabled != NP_TREE_TRUE) &&
    (num_reg_unordered == 0) &&
    (num_reg_ordered == 0) &&
    (!use_bernstein) &&
    (!int_cker_bound_extern);
  const int fit_tree_active =
    (BANDWIDTH_reg != BW_ADAP_NN) &&
    (call->tree_enabled == NP_TREE_TRUE);
  /*
   * Signed-log beta rows choose one common scale for every response/basis
   * moment.  Keep the historical Y^2 scale anchor even when its standard-
   * error consumer is disabled so requesting uncertainty can never change
   * the fitted mean or gradient.  Peer-kernel routes retain their lean
   * two-column no-error layout.
   */
  const int include_response_square = call->kernel_route != NULL ||
    (!ordinary_hc0 && call->do_merr);
  const int reuse_fit_dual_power =
    moment_errors && (!reuse_fit_kernel_row) && (!fit_tree_active);
  NP_DualPowerCtx fit_dual_power_ctx = {
    NULL, 2, NULL, NULL, 0, 0, NULL, 0, NULL, NULL, NULL, NULL, NULL, NULL
  };
  int variance_nrhs = 1;
  int first_se_nrhs = 0;
  int response_y_offset;
  int response_basis_offset;
  int moment_stride;
  int i, j, l;

  if(hc0_residual_preparing &&
     (!np_regression_residual_preparation_valid(preparation) ||
      preparation->num_obs_train != num_obs_train ||
      preparation->num_obs_eval != num_obs_eval ||
      call->do_merr || call->do_gerr || call->first_se_request != NULL)) {
    execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_HC0;
    return R_NilValue;
  }
  if((vector_glp_degree_extern == NULL) || (num_reg_continuous <= 0)) {
    execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_DEGREE;
    return R_NilValue;
  }
  if(!np_glp_build_terms(num_reg_continuous,
                         vector_glp_degree_extern,
                         int_glp_basis_extern,
                         &owner->terms,
                         &owner->nterms)) {
    execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_TERMS;
    return R_NilValue;
  }
  if(owner->nterms <= 0) {
    execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_DIMENSION;
    return R_NilValue;
  }
  response_y_offset = include_response_square ? 1 : 0;
  response_basis_offset = include_response_square ? 2 : 1;
  moment_stride = owner->nterms + response_basis_offset;
  if(call->do_grad && call->do_gerr)
    for(l = 0; l < num_reg_continuous; ++l)
      if(np_glp_gradient_direction_active(l))
        ++variance_nrhs;
  if(call->first_se_request != NULL)
    for(l = 0; l < num_reg_continuous; ++l)
      first_se_nrhs += call->first_se_request->se[l];
  if(variance_nrhs > owner->nterms) {
    execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_DIMENSION;
    return R_NilValue;
  }
  if(!np_lp_solve_workspace_reserve(&owner->solve_workspace,
                                    owner->nterms,
                                    MAX(variance_nrhs, first_se_nrhs))) {
    execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_SOLVE_ALLOC;
    return R_NilValue;
  }

  owner->basis = alloc_matd(num_obs_train, owner->nterms);
  if(BANDWIDTH_reg == BW_GEN_NN
#ifdef MPI2
     || BANDWIDTH_reg == BW_FIXED
#endif
     )
    owner->matrix_bandwidth_eval = alloc_tmatd(1, num_reg_continuous);
  owner->eval_continuous = alloc_matd(1, num_reg_continuous);
  owner->eval_unordered = alloc_matd(1, num_reg_unordered);
  owner->eval_ordered = alloc_matd(1, num_reg_ordered);
  owner->response_columns = (double **)malloc(
    (size_t)moment_stride*sizeof(double *));
  owner->basis_columns = (double **)malloc(
    (size_t)owner->nterms*sizeof(double *));
  if(prepared_ann) {
    owner->ann_outer = (NP_TreeOuterBlasWorkspace){0};
    owner->ann_current_row = call->prepared_ann_row;
    owner->ann_batch = call->prepared_ann_row == NULL ?
      np_ann_fit_batch_create(call->kernel_c[0],num_obs_train,num_obs_eval,
        num_reg_continuous,call->matrix_X_continuous_train,
        call->matrix_X_continuous_eval,call->matrix_bandwidth) : NULL;
    owner->ann_kernel_row = (double *)malloc((size_t)num_obs_train*sizeof(double));
    owner->ann_active.istart = (int *)malloc((size_t)num_obs_train*sizeof(int));
    owner->ann_active.nlev = (int *)malloc((size_t)num_obs_train*sizeof(int));
    if(owner->ann_kernel_row == NULL || owner->ann_active.istart == NULL ||
       owner->ann_active.nlev == NULL) {
      execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_ALLOC;
      return R_NilValue;
    }
  }
  if(include_response_square)
    owner->squared_response = (double *)malloc(
      (size_t)num_obs_train*sizeof(double));
  owner->moments = (double *)malloc(
    (size_t)moment_stride*(size_t)owner->nterms*sizeof(double));
  if(moment_errors)
    owner->power2_moments = (double *)malloc(
      (size_t)owner->nterms*(size_t)owner->nterms*sizeof(double));
  if((moment_errors && reuse_fit_kernel_row) ||
     (hc0_information_row && !categorical_pair))
    owner->retained_kernel_row = (double *)malloc(
      (size_t)num_obs_train*sizeof(double));
  owner->coefficient = (double *)malloc(
    (size_t)owner->nterms*sizeof(double));
  if(call->do_merr || categorical_pair)
    owner->power2_projection = (double *)malloc(
      (size_t)owner->nterms*sizeof(double));
  owner->eval_basis = (double *)malloc(
    (size_t)owner->nterms*sizeof(double));
  owner->eval_derivative = (double *)malloc(
    (size_t)owner->nterms*sizeof(double));
  if(call->do_grad && (num_reg_unordered + num_reg_ordered > 0))
    owner->categorical_point = (double *)malloc(
      (size_t)(num_reg_unordered + num_reg_ordered)*sizeof(double));
  if(categorical_hc0)
    owner->categorical_stderr = (double *)malloc(
      (size_t)(num_reg_unordered + num_reg_ordered)*sizeof(double));
  if(categorical_pair) {
    owner->categorical_base_kernel_row = (double *)malloc(
      (size_t)num_obs_train*sizeof(double));
    owner->categorical_alternate_kernel_row = (double *)malloc(
      (size_t)num_obs_train*sizeof(double));
  }
  if(use_bernstein)
    owner->basis_context = (NPGLPBasisCtx *)calloc(
      (size_t)num_reg_continuous, sizeof(NPGLPBasisCtx));

  if(owner->basis == NULL ||
     ((BANDWIDTH_reg == BW_GEN_NN
#ifdef MPI2
       || BANDWIDTH_reg == BW_FIXED
#endif
       ) &&
      owner->matrix_bandwidth_eval == NULL) ||
     owner->eval_continuous == NULL ||
     (num_reg_unordered > 0 && owner->eval_unordered == NULL) ||
     (num_reg_ordered > 0 && owner->eval_ordered == NULL) ||
     owner->response_columns == NULL || owner->basis_columns == NULL ||
     (include_response_square && owner->squared_response == NULL) ||
     owner->moments == NULL ||
     (moment_errors && owner->power2_moments == NULL) ||
     (((moment_errors && reuse_fit_kernel_row) ||
       (hc0_information_row && !categorical_pair)) &&
      owner->retained_kernel_row == NULL) ||
     owner->coefficient == NULL ||
     ((call->do_merr || categorical_pair) && owner->power2_projection == NULL) ||
     owner->eval_basis == NULL ||
     owner->eval_derivative == NULL ||
     (call->do_grad && (num_reg_unordered + num_reg_ordered > 0) &&
      owner->categorical_point == NULL) ||
     (categorical_hc0 && owner->categorical_stderr == NULL) ||
     (categorical_pair &&
      (owner->categorical_base_kernel_row == NULL ||
       owner->categorical_alternate_kernel_row == NULL)) ||
     (use_bernstein && owner->basis_context == NULL)) {
    execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_ALLOC;
    return R_NilValue;
  }

  if(moment_errors) {
    fit_dual_power_ctx.weighted_sum = owner->power2_moments;
    fit_dual_power_ctx.matrix_Y = owner->basis;
    fit_dual_power_ctx.matrix_W = owner->basis;
    fit_dual_power_ctx.ncol_Y = owner->nterms;
    fit_dual_power_ctx.ncol_W = owner->nterms;
    fit_dual_power_ctx.observation_scale = NULL;
  }

  /* Optional exact-result reuse; failure to allocate it must not make a
   * valid incumbent donor calculation fail. Required workspaces above retain
   * their existing failure policy. The basis is filled once before the row
   * loop; the residual inputs remain immutable throughout this invocation. */
  if(ordinary_hc0 && call->do_merr && call->all_large_kernel_rows)
    (void)np_inference_reuse_reserve(
      &owner->inference_reuse, owner->nterms, variance_nrhs, num_obs_train,
      owner->basis, call->hc0_context);

  if(include_response_square)
    for(i = 0; i < num_obs_train; ++i)
      owner->squared_response[i] = call->vector_Y[i]*call->vector_Y[i];

  if(use_bernstein) {
    for(l = 0; l < num_reg_continuous; ++l) {
      double xmin = call->matrix_X_continuous_train[l][0];
      double xmax = xmin;

      for(i = 1; i < num_obs_train; ++i) {
        const double xi = call->matrix_X_continuous_train[l][i];
        if(xi < xmin) xmin = xi;
        if(xi > xmax) xmax = xi;
      }
      if(!np_glp_basis_ctx_init(&owner->basis_context[l],
                                vector_glp_degree_extern[l],
                                xmin,
                                xmax,
                                int_glp_basis_extern == 1)) {
        execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_BASIS;
        return R_NilValue;
      }
    }
    np_glp_fill_basis_train(num_reg_continuous,
                            owner->terms,
                            owner->nterms,
                            call->matrix_X_continuous_train,
                            num_obs_train,
                            owner->basis_context,
                            owner->basis);
  } else {
    np_glp_fill_basis_raw_train(num_reg_continuous,
                                owner->terms,
                                owner->nterms,
                                call->matrix_X_continuous_train,
                                num_obs_train,
                                owner->basis);
  }

  if(call->fit_projection != NULL)
    np_fit_projection_basis(call->fit_projection, owner->basis, owner->nterms);

  if(call->prepared_design != NULL)
    owner->active_design = call->prepared_design;
  if(owner->active_design->n == 0)
    np_lp_design_support_prepare(owner->active_design, owner->basis,
                                  num_obs_train, owner->nterms);
  else if(owner->active_design->n != num_obs_train ||
          owner->active_design->p != owner->nterms)
    error("LP fitting design identity shape changed within an invocation");

#ifdef MPI2
	    const int use_mpi_owner_reduce_lp =
	      (iNum_Processors > 1) &&
	      (!np_mpi_local_regression_active()) &&
	      (BANDWIDTH_reg == BW_FIXED);
	    const int owner_chunk_rows_lp =
	      np_reg_mpi_owner_chunk_rows(num_obs_eval,
	                                  1 + call->do_merr + 2*hc0_residual_preparing +
	                                  (call->do_grad ? num_reg_continuous*(1 + (call->do_gerr ? 1 : 0)) : 0));
#endif

  for(j = 0; j < num_obs_eval; ++j) {
    NPLPSupportRankContext support_rank = {0, owner->nterms, owner->active_design};
    double sk, ey, ey2, sigma2hat;
    double * const hc0_kernel_row = categorical_pair ?
      owner->categorical_base_kernel_row :
      (hc0_information_row ? owner->retained_kernel_row : NULL);
    const int divide_retained_weights = categorical_pair ||
      (hc0_information_row && !reuse_fit_kernel_row);
    const double hc0_kernel_divisor =
      reuse_fit_kernel_row && !categorical_pair ? call->bandwidth_product : 1.0;
    double pristine_anchor = 0.0;
    NPLPSolvePolicyDiagnostics solve_diagnostics = {0, 0.0};
    int have_vcov = 0;

#ifdef MPI2
	      if(use_mpi_owner_reduce_lp && ((j % owner_chunk_rows_lp) == 0)){
	        const int chunk_start = j;
	        const int chunk_end = MIN(num_obs_eval, chunk_start + owner_chunk_rows_lp);
	        const int owner_row_width_lp = 1 + call->do_merr + 2*hc0_residual_preparing +
	          (call->do_grad ? num_reg_continuous*(1 + (call->do_gerr ? 1 : 0)) : 0);
	        int local_pos = 0;
	        int owner_solve_failed = 0;

	        np_reg_mpi_owner_chunk_init(&owner->mpi_owner_chunk,
	                                    chunk_start,
	                                    chunk_end,
	                                    owner_row_width_lp,
	                                    "npreg LP");
	        if(call->do_merr || hc0_residual_preparing)
	          owner->mpi_kernel_row = (double *)malloc(
	            (size_t)num_obs_train*sizeof(double));
	        if((call->do_merr || hc0_residual_preparing) && owner->mpi_kernel_row == NULL)
	          error("\n** Error: memory allocation failed.");

	        local_pos = 0;
	        for(int jj = chunk_start + my_rank; jj < chunk_end; jj += iNum_Processors){
	          double sk_owner, ey_owner, ey2_owner, sigma2_owner;
	          double pristine_anchor_owner = 0.0;
	          NPLPSolvePolicyDiagnostics solve_diagnostics_owner = {0, 0.0};
	          int have_vcov_owner = 0;
	          double * const out_owner = owner->mpi_owner_chunk.sendbuf + (size_t)local_pos*(size_t)owner_row_width_lp;

	          for(l = 0; l < num_reg_continuous; l++){
	            owner->eval_continuous[l][0] =
	              call->matrix_X_continuous_eval[l][jj];
	            owner->matrix_bandwidth_eval[l][0] =
	              (BANDWIDTH_reg == BW_GEN_NN) ? call->matrix_bandwidth[l][jj] :
	              ((BANDWIDTH_reg == BW_FIXED) ? call->matrix_bandwidth[l][0] : 0.0);
	          }
	          for(l = 0; l < num_reg_unordered; l++)
	            owner->eval_unordered[l][0] =
	              call->matrix_X_unordered_eval[l][jj];
	          for(l = 0; l < num_reg_ordered; l++)
	            owner->eval_ordered[l][0] =
	              call->matrix_X_ordered_eval[l][jj];

	          if(include_response_square)
	            owner->response_columns[0] = owner->squared_response;
	          owner->response_columns[response_y_offset] = call->vector_Y;
	          for(l = 0; l < owner->nterms; l++){
	            owner->response_columns[l + response_basis_offset] = owner->basis[l];
	            owner->basis_columns[l] = owner->basis[l];
	          }

	          if(call->kernel_route == NULL)
	            kernel_weighted_sum_np_ctx_ex(call->kernel_c,
	                                   call->kernel_u,
	                                   call->kernel_o,
	                                   BANDWIDTH_reg,
	                                   num_obs_train,
	                                   1,
	                                   num_reg_unordered,
	                                   num_reg_ordered,
	                                   num_reg_continuous,
	                                   0,
	                                   0,
	                                   1,
	                                   1,
	                                   1,
	                                   0,
	                                   0,
	                                   0,
	                                   0,
	                                   call->operator,
	                                   OP_NOOP,
	                                   0,
	                                   0,
	                                   NULL,
	                                   1,
	                                   moment_stride,
	                                   owner->nterms,
	                                   (BANDWIDTH_reg == BW_ADAP_NN) ? NP_TREE_FALSE : call->tree_enabled,
	                                   0,
	                                   (BANDWIDTH_reg == BW_ADAP_NN) ? NULL : kdt_extern_X,
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
	                                   call->matrix_bandwidth,
	                                   owner->matrix_bandwidth_eval,
	                                   call->lambda,
	                                   call->num_categories,
	                                   call->matrix_categorical_vals,
	                                   NULL,
	                                   owner->moments,
	                                   NULL,
	                                   owner->mpi_kernel_row,
	                                   call->gate_context,
	                                   NULL, NULL, NULL, NULL, NULL, &support_rank);
	          else {
	            /*
	             * Rank ownership changes only who computes this row.  It must
	             * not change the selected continuous-kernel family: beta rows
	             * retain the same canonical scaled-moment owner used by the
	             * serial and non-owner LP routes.
	             */
	            if(np_beta_regression_lp_moment_row_canonical(
	                 BANDWIDTH_reg,
	                 num_obs_train,
	                 num_reg_unordered,
	                 num_reg_ordered,
	                 num_reg_continuous,
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
	                 call->matrix_bandwidth,
	                 call->matrix_bandwidth,
	                 call->lambda,
	                 call->num_categories,
	                 call->matrix_categorical_vals,
	                 call->categorical_compress,
	                 owner->moments,
	                 moment_errors ? owner->power2_moments : NULL,
	                 ordinary_hc0 ?
	                   call->hc0_context->scaled_residual : NULL,
	                 owner->mpi_kernel_row,
	                 call->kernel_route,
	                 call->kernel_route_diagnostics,
	                 &support_rank) != 0) {
	              owner_solve_failed = 1;
	              break;
	            }
	          }

	          for(i = 0; i < owner->nterms; i++){
	            const int base = i*moment_stride;
	            owner->solve_workspace.rhs_source[i] =
	              owner->moments[base + response_y_offset];
	            for(l = 0; l < owner->nterms; l++)
	              owner->solve_workspace.gram_source[i+l*owner->nterms] =
	                owner->moments[base + l + response_basis_offset];
	          }
	          pristine_anchor_owner = owner->solve_workspace.gram_source[0];
	          if(np_lp_solve_workspace_solve_response_ranked(
	               &owner->solve_workspace,
	               owner->nterms,
	               1,
	               epsilon,
	               support_rank.count,
	               &solve_diagnostics_owner) != NP_LP_SOLVE_POLICY_OK){
	            owner_solve_failed = 1;
	            break;
	          }
	          for(i = 0; i < owner->nterms; i++)
	            owner->coefficient[i] = owner->solve_workspace.rhs_work[i];

	          if(use_bernstein)
	            np_glp_fill_basis_eval(num_reg_continuous,
	                                   owner->terms,
	                                   owner->nterms,
	                                   call->matrix_X_continuous_eval,
	                                   jj,
	                                   owner->basis_context,
	                                   owner->eval_basis);
	          else
	            np_glp_fill_basis_eval_raw(num_reg_continuous,
	                                       owner->terms,
	                                       owner->nterms,
	                                       call->matrix_X_continuous_eval,
	                                       jj,
	                                       owner->eval_basis);

	          out_owner[0] = preserve_point ? call->mean[jj] : 0.0;
	          if(!preserve_point)
	            for(i = 0; i < owner->nterms; i++)
	              out_owner[0] += owner->eval_basis[i]*owner->coefficient[i];

	          if(hc0_residual_preparing) {
	            for(i = 0; i < owner->nterms; ++i)
	              owner->solve_workspace.rhs_source[i] = owner->eval_basis[i];
	            if(np_lp_solve_workspace_solve_adjoint_factored(
	                 &owner->solve_workspace, owner->nterms, 1,
	                 pristine_anchor_owner, &solve_diagnostics_owner) !=
	                 NP_LP_SOLVE_POLICY_OK ||
	               !np_regression_hc0_lp_prepare_row(
	                 owner->basis, owner->nterms, owner->solve_workspace.rhs_work,
	                 owner->mpi_kernel_row, preparation, jj,
	                 solve_diagnostics_owner.ridge_total == 0.0)) {
	              owner_solve_failed = 1;
	              break;
	            }
	            const int donor = preparation->evaluation_to_donor[jj];
	            out_owner[1] = preparation->normalized_residual[donor];
	            out_owner[2] = preparation->information[donor];
	          }

	          sigma2_owner = 0.0;
	          if(call->do_merr && !ordinary_hc0) {
	            sk_owner = copysign(DBL_MIN,
	              owner->moments[response_basis_offset]) +
	              owner->moments[response_basis_offset];
	            ey_owner = owner->moments[response_y_offset]/sk_owner;
	            ey2_owner = owner->moments[0]/sk_owner;
	            sigma2_owner = ey2_owner - ey_owner*ey_owner;
	            {
	              const double v_owner = sigma2_owner *
	                call->kernel_squared_integral /
	                (sk_owner*call->bandwidth_product);
	              out_owner[1] = (v_owner <= 0.0) ? 0.0 : sqrt(v_owner);
	            }
	            sigma2_owner = (sigma2_owner <= 0.0) ? 0.0 : sigma2_owner;
	          } else if(call->do_merr) {
	            out_owner[1] = 0.0;
	          }

	          if(moment_errors && call->kernel_route == NULL) {
	            for(i = 0; i < owner->nterms*owner->nterms; i++)
	              owner->power2_moments[i] = 0.0;

	              for(i = 0; i < num_obs_train; i++){
	                const double w = owner->mpi_kernel_row[i];
	                const double w2 = w*w;
	                if(w2 == 0.0)
	                  continue;
	                for(int a = 0; a < owner->nterms; a++){
	                  const double za = owner->basis[a][i];
	                  for(int b = 0; b < owner->nterms; b++)
	                    owner->power2_moments[a*owner->nterms+b] +=
	                      za*owner->basis[b][i]*w2;
	                }
	              }
	          }
	          {
	            int variance_rhs = 1;
	            if(call->do_merr)
	              for(i = 0; i < owner->nterms; i++)
	                owner->solve_workspace.rhs_source[i] = owner->eval_basis[i];

	            if(call->do_grad){
	              for(l = 0; l < num_reg_continuous; l++){
	                const int grad_order =
	                  (vector_glp_gradient_order_extern != NULL) ?
	                  MAX(1, vector_glp_gradient_order_extern[l]) : 1;
	                const int active = np_glp_gradient_direction_active(l);
	                double * const direction =
	                  (call->do_gerr && active) ?
	                  owner->solve_workspace.rhs_source +
	                    (size_t)variance_rhs*(size_t)owner->nterms :
	                  owner->eval_derivative;
	                double grad_value = 0.0;
	                const int output_offset =
	                  1 + call->do_merr + 2*hc0_residual_preparing +
	                  l*(1 + (call->do_gerr ? 1 : 0));

	                if(use_bernstein)
	                  np_glp_fill_basis_eval_deriv(l,
	                                               grad_order,
	                                               num_reg_continuous,
	                                               owner->terms,
	                                               owner->nterms,
	                                               call->matrix_X_continuous_eval,
	                                               jj,
	                                               owner->basis_context,
	                                               direction);
	                else
	                  np_glp_fill_basis_eval_deriv_raw(l,
	                                                  grad_order,
	                                                  num_reg_continuous,
	                                                  owner->terms,
	                                                  owner->nterms,
	                                                  call->matrix_X_continuous_eval,
	                                                  jj,
	                                                  direction);
	                if(active)
	                  for(i = 0; i < owner->nterms; i++)
	                    grad_value += direction[i]*owner->coefficient[i];
	                out_owner[output_offset] = active ?
	                  (preserve_point ?
	                   call->gradient[l][jj] : grad_value) : NA_REAL;
	                if(call->do_gerr){
	                  out_owner[output_offset + 1] =
	                    active ? 0.0 : NA_REAL;
	                  if(active)
	                    variance_rhs++;
	                }
	              }
	            }

	            if(call->do_merr) {
	              if(variance_rhs != variance_nrhs){
	                owner_solve_failed = 1;
	                break;
	              }
	              if(ordinary_hc0) {
	                if(np_lp_solve_workspace_solve_adjoint_factored(
	                     &owner->solve_workspace,
	                     owner->nterms,
	                     variance_rhs,
	                     pristine_anchor_owner,
	                     &solve_diagnostics_owner) != NP_LP_SOLVE_POLICY_OK) {
	                  owner_solve_failed = 1;
	                  break;
	                }
	                have_vcov_owner = 1;
	              } else if(np_conditional_lp_project_accepted(
	                          &owner->solve_workspace,
	                          owner->nterms,
	                          variance_rhs,
	                          pristine_anchor_owner,
	                          &solve_diagnostics_owner)) {
	                have_vcov_owner = 1;
	              } else {
	                out_owner[1] = NA_REAL;
	                if(call->do_grad && call->do_gerr)
	                  for(l = 0; l < num_reg_continuous; ++l)
	                    if(np_glp_gradient_direction_active(l))
	                      out_owner[2 + call->do_merr + l*2] = NA_REAL;
	              }
	            }

	            if(have_vcov_owner){
	              if(ordinary_hc0) {
	                if(owner->inference_reuse.entry != NULL)
	                  np_inference_reuse_begin_row(
	                    &owner->inference_reuse, owner->basis, call->hc0_context,
	                    owner->mpi_kernel_row, num_obs_train, 1.0,
	                    np_progress_fit_heartbeat);
	                if(!np_regression_hc0_lp_standard_error_reuse(
	                     &owner->inference_reuse, 0,
	                     owner->solve_workspace.rhs_work,
	                     owner->basis, owner->nterms, owner->mpi_kernel_row, 1.0,
	                     call->hc0_context, &out_owner[1], NULL)) {
	                  owner_solve_failed = 1;
	                  break;
	                }
	              } else {
	                const double q = np_lp_variance_quadratic(
	                  owner->solve_workspace.rhs_work,
	                  owner->power2_moments,
	                  owner->power2_projection,
	                  owner->nterms);
	                out_owner[1] =
	                  np_conditional_lp_standard_error(sigma2_owner, q);
	              }

	              if(call->do_grad && call->do_gerr){
	                int rhs = 1;
	                for(l = 0; l < num_reg_continuous; l++){
	                  if(np_glp_gradient_direction_active(l)){
	                    const double * const projection =
	                      owner->solve_workspace.rhs_work +
	                        (size_t)rhs*(size_t)owner->nterms;
	                    if(ordinary_hc0) {
	                      if(!np_regression_hc0_lp_standard_error_reuse(
	                           &owner->inference_reuse, rhs,
	                           projection,
	                           owner->basis, owner->nterms, owner->mpi_kernel_row, 1.0,
	                           call->hc0_context,
	                           &out_owner[2 + call->do_merr + l*2], NULL)) {
	                        owner_solve_failed = 1;
	                        break;
	                      }
	                    } else {
	                      const double grad_q = np_lp_variance_quadratic(
	                        projection,
	                        owner->power2_moments,
	                        owner->power2_projection,
	                        owner->nterms);
	                      out_owner[2 + call->do_merr + l*2] =
	                        np_conditional_lp_standard_error(sigma2_owner, grad_q);
	                    }
	                    rhs++;
	                  }
	                }
	                if(owner_solve_failed)
	                  break;
	              }
	            }
	          }
	          local_pos++;
	        }

	        if(owner_solve_failed){
	          for(int pos = local_pos; pos < owner->mpi_owner_chunk.local_rows; pos++){
	            double * const failed_row =
	              owner->mpi_owner_chunk.sendbuf + (size_t)pos*(size_t)owner_row_width_lp;
	            for(int col = 0; col < owner_row_width_lp; col++)
	              failed_row[col] = 0.0;
	          }
	          owner->mpi_owner_chunk.sendbuf[0] = NAN;
	        }

	        np_reg_mpi_owner_chunk_allgather(&owner->mpi_owner_chunk);

	        owner_solve_failed = 0;
	        for(i = 0; i < iNum_Processors; i++)
	          if((owner->mpi_owner_chunk.recvcounts[i] > 0) &&
	             !isfinite(np_reg_mpi_owner_chunk_recv_ptr(
	               &owner->mpi_owner_chunk,
	               i,
	               0)[0]))
	            owner_solve_failed = 1;

	        if(!owner_solve_failed){
	          for(i = 0; i < iNum_Processors; i++){
	            int pos_i = 0;
	            for(int jj = chunk_start + i; jj < chunk_end; jj += iNum_Processors){
	              const double * const in =
	                np_reg_mpi_owner_chunk_recv_ptr(&owner->mpi_owner_chunk, i, pos_i);
	              int ipos = 0;
	              call->mean[jj] = in[ipos++];
	              if(hc0_residual_preparing) {
	                const int donor = preparation->evaluation_to_donor[jj];
	                preparation->normalized_residual[donor] = in[ipos++];
	                preparation->information[donor] = (int)in[ipos++];
	              }
	              if(call->do_merr) {
	                call->mean_stderr[jj] = in[ipos++];
	                if(ordinary_hc0 && call->hc0_context->mean_unavailable != NULL &&
	                   ISNA(call->mean_stderr[jj]))
	                  call->hc0_context->mean_unavailable[jj] = 1;
	              }
	              if(call->do_grad){
	                for(l = 0; l < num_reg_continuous; l++){
	                  call->gradient[l][jj] = in[ipos++];
	                  if(call->do_gerr) {
	                    call->gradient_stderr[l][jj] = in[ipos++];
	                    if(ordinary_hc0 && np_glp_gradient_direction_active(l) &&
	                       call->hc0_context->gradient_unavailable != NULL &&
	                       ISNA(call->gradient_stderr[l][jj]))
	                      call->hc0_context->gradient_unavailable[l][jj] = 1;
	                  }
	                }
	                for(l = num_reg_continuous;
	                    l < (num_reg_continuous + num_reg_unordered + num_reg_ordered);
	                    l++){
	                  call->gradient[l][jj] = NA_REAL;
	                  if(call->do_gerr)
	                    call->gradient_stderr[l][jj] = NA_REAL;
	                }
	              }
	              pos_i++;
	            }
	          }
	        }

	        np_reg_mpi_owner_chunk_free(&owner->mpi_owner_chunk);
	        free(owner->mpi_kernel_row);
	        owner->mpi_kernel_row = NULL;
	        if(owner_solve_failed){
	          execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_OWNER_SOLVE;
	          return R_NilValue;
	        }
	      }

	      if(use_mpi_owner_reduce_lp){
	        if (call->fit_progress_active)
	          np_progress_fit_loop_step(j + 1, call->fit_progress_total);
	        continue;
	      }
#endif

    for(l = 0; l < num_reg_continuous; ++l) {
      owner->eval_continuous[l][0] = call->matrix_X_continuous_eval[l][j];
      if(BANDWIDTH_reg == BW_GEN_NN)
        owner->matrix_bandwidth_eval[l][0] = call->matrix_bandwidth[l][j];
    }
    for(l = 0; l < num_reg_unordered; ++l)
      owner->eval_unordered[l][0] = call->matrix_X_unordered_eval[l][j];
    for(l = 0; l < num_reg_ordered; ++l)
      owner->eval_ordered[l][0] = call->matrix_X_ordered_eval[l][j];

    if(include_response_square)
      owner->response_columns[0] = owner->squared_response;
    owner->response_columns[response_y_offset] = call->vector_Y;
    for(l = 0; l < owner->nterms; ++l) {
      owner->response_columns[l + response_basis_offset] = owner->basis[l];
      owner->basis_columns[l] = owner->basis[l];
    }

    if(prepared_ann) {
      if(owner->ann_batch != NULL)
        owner->ann_current_row=np_ann_fit_batch_row(owner->ann_batch,j,num_obs_eval,1);
      if(np_ann_fit_moment_row(call, owner, moment_stride, owner->moments,
           moment_errors ? owner->power2_moments : NULL,
           hc0_kernel_row, &support_rank) != 0) {
        execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_ROUTE;
        return R_NilValue;
      }
    } else if(NP_UNLIKELY(call->kernel_route != NULL)) {
      if(np_beta_regression_lp_moment_row_canonical(
           BANDWIDTH_reg,
           num_obs_train,
           num_reg_unordered,
           num_reg_ordered,
           num_reg_continuous,
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
           call->matrix_bandwidth,
           BANDWIDTH_reg == BW_GEN_NN ?
             owner->matrix_bandwidth_eval : call->matrix_bandwidth,
           call->lambda,
           call->num_categories,
           call->matrix_categorical_vals,
	           call->categorical_compress,
	           owner->moments,
	           moment_errors ? owner->power2_moments : NULL,
           ordinary_hc0 ? call->hc0_context->scaled_residual : NULL,
           hc0_kernel_row,
           call->kernel_route,
           call->kernel_route_diagnostics,
           &support_rank) != 0) {
        execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_ROUTE;
        return R_NilValue;
      }
    } else {
      const int moment_status = kernel_weighted_sum_np_ctx_ex(call->kernel_c,
                         call->kernel_u,
                         call->kernel_o,
                         BANDWIDTH_reg,
                         num_obs_train,
                         1,
                         num_reg_unordered,
                         num_reg_ordered,
                         num_reg_continuous,
                         0, 0, 1, 1, divide_retained_weights, 0, 0, 0, 0,
                         call->operator,
                         OP_NOOP,
                         0, 0, NULL, 1,
                         moment_stride,
                         owner->nterms,
                         BANDWIDTH_reg == BW_ADAP_NN ?
                           NP_TREE_FALSE : call->tree_enabled,
                         0,
                         BANDWIDTH_reg == BW_ADAP_NN ? NULL : kdt_extern_X,
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
                         call->matrix_bandwidth,
                         BANDWIDTH_reg == BW_GEN_NN ?
                           owner->matrix_bandwidth_eval : call->matrix_bandwidth,
                         call->lambda,
                         call->num_categories,
                         call->matrix_categorical_vals,
                         NULL,
                         owner->moments,
                         NULL,
                         hc0_kernel_row != NULL ? hc0_kernel_row :
                           ((moment_errors && reuse_fit_kernel_row) ?
                            owner->retained_kernel_row : NULL),
                         call->gate_context,
                         reuse_fit_dual_power ? &fit_dual_power_ctx : NULL,
                         NULL, NULL, NULL,
#ifdef MPI2
                         0,
#endif
                         NULL,
                         &support_rank);
      if((ordinary_hc0 || hc0_residual_preparing) && moment_status != 0) {
        execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_HC0;
        return R_NilValue;
      }
    }

    for(i = 0; i < owner->nterms; ++i) {
      const int base = i*moment_stride;
      owner->solve_workspace.rhs_source[i] =
        owner->moments[base + response_y_offset];
      for(l = 0; l < owner->nterms; ++l)
        owner->solve_workspace.gram_source[i+l*owner->nterms] =
          owner->moments[base + l + response_basis_offset];
    }
    pristine_anchor = owner->solve_workspace.gram_source[0];

    if(np_lp_solve_workspace_solve_response_ranked(
         &owner->solve_workspace,
         owner->nterms,
         1,
         epsilon,
         support_rank.count,
         &solve_diagnostics) !=
       NP_LP_SOLVE_POLICY_OK) {
      if(np_regression_general_lp_empty_row(prepared_ann, call, owner, j, 0,
           moment_stride, divide_retained_weights,
           hc0_kernel_row != NULL ? hc0_kernel_row :
             ((moment_errors && reuse_fit_kernel_row) ? owner->retained_kernel_row : NULL))) {
        call->mean[j] = NA_REAL;
        if(call->do_merr) call->mean_stderr[j] = NA_REAL;
        if(call->do_grad || call->do_gerr)
          for(l = 0; l < num_reg_continuous + num_reg_unordered + num_reg_ordered; ++l) {
            if(call->do_grad) call->gradient[l][j] = NA_REAL;
            if(call->do_gerr) call->gradient_stderr[l][j] = NA_REAL;
          }
#ifdef MPI2
        if(call->fit_progress_active)
#endif
          np_progress_fit_loop_step(j + 1, call->fit_progress_total);
        continue;
      }
      execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_SOLVE;
      return R_NilValue;
    }
    for(i = 0; i < owner->nterms; ++i)
      owner->coefficient[i] = owner->solve_workspace.rhs_work[i];

    if(use_bernstein)
      np_glp_fill_basis_eval(num_reg_continuous,
                             owner->terms,
                             owner->nterms,
                             call->matrix_X_continuous_eval,
                             j,
                             owner->basis_context,
                             owner->eval_basis);
    else
      np_glp_fill_basis_eval_raw(num_reg_continuous,
                                 owner->terms,
                                 owner->nterms,
                                 call->matrix_X_continuous_eval,
                                 j,
                                 owner->eval_basis);
    if(!preserve_point) {
      call->mean[j] = 0.0;
      for(i = 0; i < owner->nterms; ++i)
        call->mean[j] += owner->eval_basis[i]*owner->coefficient[i];
    }
    if(call->fit_projection != NULL) {
      for(i = 0; i < owner->nterms; ++i)
        owner->solve_workspace.rhs_source[i] = owner->eval_basis[i];
      if(np_lp_solve_workspace_solve_adjoint_factored(
           &owner->solve_workspace, owner->nterms, 1,
           pristine_anchor, &solve_diagnostics) != NP_LP_SOLVE_POLICY_OK)
        error("accepted fit projection solve failed");
      np_fit_projection_row(call->fit_projection,
                            owner->solve_workspace.rhs_work, j);
    }
    if(hc0_residual_preparing) {
      for(i = 0; i < owner->nterms; ++i)
        owner->solve_workspace.rhs_source[i] = owner->eval_basis[i];
      if(np_lp_solve_workspace_solve_adjoint_factored(
           &owner->solve_workspace, owner->nterms, 1,
           pristine_anchor, &solve_diagnostics) != NP_LP_SOLVE_POLICY_OK ||
         !np_regression_hc0_lp_prepare_row(
           owner->basis, owner->nterms, owner->solve_workspace.rhs_work,
           hc0_kernel_row, preparation, j,
           solve_diagnostics.ridge_total == 0.0)) {
        execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_HC0;
        return R_NilValue;
      }
    }
    sigma2hat = 0.0;
    if(call->do_merr && !ordinary_hc0) {
      sk = copysign(DBL_MIN, owner->moments[response_basis_offset]) +
        owner->moments[response_basis_offset];
      ey = owner->moments[response_y_offset]/sk;
      ey2 = owner->moments[0]/sk;
      sigma2hat = ey2 - ey*ey;
      call->mean_stderr[j] = sigma2hat <= 0.0 ? 0.0 :
        sqrt(sigma2hat*call->kernel_squared_integral/
             (sk*call->bandwidth_product));
      sigma2hat = sigma2hat <= 0.0 ? 0.0 : sigma2hat;
      if(call->conditional_variance != NULL)
        call->conditional_variance[j] = sigma2hat;
    } else if(call->do_merr) {
      call->mean_stderr[j] = 0.0;
    }

    for(l = 0; l < owner->nterms; ++l) {
      owner->response_columns[l] = owner->basis[l];
      owner->basis_columns[l] = owner->basis[l];
    }

    if(moment_errors && reuse_fit_kernel_row) {
      if(!np_lp_power2_moments_from_kernel_row(
           owner->basis,
           owner->nterms,
           owner->retained_kernel_row,
           num_obs_train,
           call->bandwidth_product,
           ordinary_hc0 ? call->hc0_context->scaled_residual : NULL,
           owner->power2_moments)) {
        execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_HC0;
        return R_NilValue;
      }
    } else if(moment_errors && !reuse_fit_dual_power) {
      kernel_weighted_sum_np_ctx(call->kernel_c,
                         call->kernel_u,
                         call->kernel_o,
                         BANDWIDTH_reg,
                         num_obs_train,
                         1,
                         num_reg_unordered,
                         num_reg_ordered,
                         num_reg_continuous,
                         0, 0, 2, 1, 0, 0, 0, 0, 0,
                         call->operator,
                         OP_NOOP,
                         0, 0, NULL, 1,
                         owner->nterms,
                         owner->nterms,
                         BANDWIDTH_reg == BW_ADAP_NN ?
                           NP_TREE_FALSE : call->tree_enabled,
                         0,
                         BANDWIDTH_reg == BW_ADAP_NN ? NULL : kdt_extern_X,
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
                         call->matrix_bandwidth,
                         BANDWIDTH_reg == BW_GEN_NN ?
                           owner->matrix_bandwidth_eval : call->matrix_bandwidth,
                         call->lambda,
                         call->num_categories,
                         call->matrix_categorical_vals,
                         NULL,
                         owner->power2_moments,
                         NULL, NULL,
                         call->gate_context,
                         NULL);
    }

    {
      int variance_rhs = 1;
      if(call->do_merr)
        for(i = 0; i < owner->nterms; ++i)
          owner->solve_workspace.rhs_source[i] = owner->eval_basis[i];

      if(call->do_grad) {
        for(l = 0; l < num_reg_continuous; ++l) {
          const int grad_order = vector_glp_gradient_order_extern != NULL ?
            MAX(1, vector_glp_gradient_order_extern[l]) : 1;
          const int active = np_glp_gradient_direction_active(l);
          double * const direction =
            (call->do_gerr && active) ?
            owner->solve_workspace.rhs_source +
              (size_t)variance_rhs*(size_t)owner->nterms :
            owner->eval_derivative;

          if(use_bernstein)
            np_glp_fill_basis_eval_deriv(l,
                                         grad_order,
                                         num_reg_continuous,
                                         owner->terms,
                                         owner->nterms,
                                         call->matrix_X_continuous_eval,
                                         j,
                                         owner->basis_context,
                                         direction);
          else
            np_glp_fill_basis_eval_deriv_raw(l,
                                             grad_order,
                                             num_reg_continuous,
                                             owner->terms,
                                             owner->nterms,
                                             call->matrix_X_continuous_eval,
                                             j,
                                             direction);
          if(!preserve_point && active) {
            call->gradient[l][j] = 0.0;
            for(i = 0; i < owner->nterms; ++i)
              call->gradient[l][j] += direction[i]*owner->coefficient[i];
          } else if(!active)
            call->gradient[l][j] = NA_REAL;
          if(call->do_gerr) {
            call->gradient_stderr[l][j] = active ? 0.0 : NA_REAL;
            if(active)
              ++variance_rhs;
          }
        }
      }

      if(call->do_merr) {
        if(variance_rhs != variance_nrhs) {
          execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_DIMENSION;
          return R_NilValue;
        }
        if(ordinary_hc0) {
          if(np_lp_solve_workspace_solve_adjoint_factored(
               &owner->solve_workspace,
               owner->nterms,
               variance_rhs,
               pristine_anchor,
               &solve_diagnostics) != NP_LP_SOLVE_POLICY_OK) {
            execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_HC0;
            return R_NilValue;
          }
          have_vcov = 1;
        } else if(np_conditional_lp_project_accepted(
                    &owner->solve_workspace,
                    owner->nterms,
                    variance_rhs,
                    pristine_anchor,
                    &solve_diagnostics)) {
          have_vcov = 1;
        } else {
          call->mean_stderr[j] = NA_REAL;
          if(call->do_grad && call->do_gerr)
            for(l = 0; l < num_reg_continuous; ++l)
              if(np_glp_gradient_direction_active(l))
                call->gradient_stderr[l][j] = NA_REAL;
        }
      }

      if(have_vcov) {
        if(ordinary_hc0) {
          if(owner->inference_reuse.entry != NULL)
            np_inference_reuse_begin_row(
              &owner->inference_reuse, owner->basis, call->hc0_context,
              hc0_kernel_row, num_obs_train, hc0_kernel_divisor,
              np_progress_fit_heartbeat);
          if(!np_regression_hc0_lp_standard_error_reuse(
               &owner->inference_reuse, 0,
               owner->solve_workspace.rhs_work,
               owner->basis, owner->nterms, hc0_kernel_row, hc0_kernel_divisor,
               call->hc0_context, &call->mean_stderr[j],
               call->hc0_context->mean_unavailable != NULL ?
                 &call->hc0_context->mean_unavailable[j] : NULL)) {
            execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_HC0;
            return R_NilValue;
          }
        } else {
          const double q = np_lp_variance_quadratic(
            owner->solve_workspace.rhs_work,
            owner->power2_moments,
            owner->power2_projection,
            owner->nterms);
          call->mean_stderr[j] =
            np_conditional_lp_standard_error(sigma2hat, q);
        }

        if(call->do_grad && call->do_gerr) {
          int rhs = 1;
          for(l = call->continuous_range.begin;
              l < call->continuous_range.end; ++l) {
            if(np_glp_gradient_direction_active(l)) {
              const double * const projection =
                owner->solve_workspace.rhs_work +
                  (size_t)rhs*(size_t)owner->nterms;
              if(ordinary_hc0) {
                if(!np_regression_hc0_lp_standard_error_reuse(
                     &owner->inference_reuse, rhs,
                     projection,
                     owner->basis, owner->nterms, hc0_kernel_row, hc0_kernel_divisor,
                     call->hc0_context, &call->gradient_stderr[l][j],
                     call->hc0_context->gradient_unavailable != NULL ?
                       &call->hc0_context->gradient_unavailable[l][j] : NULL)) {
                  execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_HC0;
                  return R_NilValue;
                }
              } else {
                const double grad_q = np_lp_variance_quadratic(
                  projection,
                  owner->power2_moments,
                  owner->power2_projection,
                  owner->nterms);
                call->gradient_stderr[l][j] =
                  np_conditional_lp_standard_error(sigma2hat, grad_q);
              }
              ++rhs;
            }
          }
        }
      }

      if(first_se_nrhs > 0) {
        int rhs = 0;
        /* This independent block reuses the accepted factorization, never
         * the level block's success flag or partially solved RHS. */
        for(l = 0; l < num_reg_continuous; ++l) {
          double *direction;
          if(!call->first_se_request->se[l])
            continue;
          call->gradient_stderr[l][j] = NA_REAL;
          direction = owner->solve_workspace.rhs_source +
            (size_t)rhs*(size_t)owner->nterms;
          if(use_bernstein)
            np_glp_fill_basis_eval_deriv(
              l, call->first_se_request->order[l], num_reg_continuous,
              owner->terms, owner->nterms, call->matrix_X_continuous_eval,
              j, owner->basis_context, direction);
          else
            np_glp_fill_basis_eval_deriv_raw(
              l, call->first_se_request->order[l], num_reg_continuous,
              owner->terms, owner->nterms, call->matrix_X_continuous_eval,
              j, direction);
          ++rhs;
        }
        if(np_conditional_lp_project_accepted(
             &owner->solve_workspace, owner->nterms, first_se_nrhs,
             pristine_anchor, &solve_diagnostics)) {
          rhs = 0;
          for(l = 0; l < num_reg_continuous; ++l) {
            double q;
            if(!call->first_se_request->se[l])
              continue;
            q = np_lp_variance_quadratic(
              owner->solve_workspace.rhs_work +
                (size_t)rhs*(size_t)owner->nterms,
              owner->power2_moments, owner->power2_projection, owner->nterms);
            call->gradient_stderr[l][j] =
              np_conditional_lp_standard_error(sigma2hat, q);
            ++rhs;
          }
        }
      }

      if(categorical_pair && !call->categorical_base_requires_refit) {
        /* Mean/continuous coefficients above remain untouched. Reuse an
         * existing accepted mean adjoint whenever SEs or preparation needed
         * it; a point-only categorical request adds only this back-solve. */
        if(!ordinary_hc0 && !hc0_residual_preparing) {
          for(i = 0; i < owner->nterms; ++i)
            owner->solve_workspace.rhs_source[i] = owner->eval_basis[i];
          if(np_lp_solve_workspace_solve_adjoint_factored(
               &owner->solve_workspace, owner->nterms, 1,
               pristine_anchor, &solve_diagnostics) != NP_LP_SOLVE_POLICY_OK) {
            execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_SOLVE;
            return R_NilValue;
          }
        }
        memcpy(owner->power2_projection,
               owner->solve_workspace.rhs_work,
               (size_t)owner->nterms*sizeof(double));
      }

      if(call->do_grad &&
         (!preserve_point || categorical_pair) &&
         (num_reg_unordered + num_reg_ordered > 0)) {
        if(include_response_square)
          owner->response_columns[0] = owner->squared_response;
        owner->response_columns[response_y_offset] = call->vector_Y;
        for(l = 0; l < owner->nterms; ++l) {
          owner->response_columns[l + response_basis_offset] =
            owner->basis[l];
          owner->basis_columns[l] = owner->basis[l];
        }
        if(BANDWIDTH_reg == BW_GEN_NN)
          for(l = 0; l < num_reg_continuous; ++l)
            owner->matrix_bandwidth_eval[l][0] =
              call->categorical_matrix_bandwidth[l][j];
        if(!np_regression_general_lp_categorical_points(
             prepared_ann,
             call,
             owner,
             moment_stride,
             response_y_offset,
             response_basis_offset,
             epsilon,
             call->mean[j],
             j,
             categorical_hc0)) {
          execution->status = categorical_hc0 ?
            NP_REGRESSION_GENERAL_LP_FIT_ERR_HC0 :
            NP_REGRESSION_GENERAL_LP_FIT_ERR_SOLVE;
          return R_NilValue;
        }
      }

      if(call->do_grad) {
        for(l = call->categorical_range.begin;
            l < call->categorical_range.end;
            ++l) {
          if(!preserve_point || categorical_pair)
            call->gradient[l][j] =
              owner->categorical_point[l - num_reg_continuous];
          if(call->do_gerr)
            call->gradient_stderr[l][j] =
              call->categorical_point_invariant ? 0.0 :
              (categorical_hc0 ?
               owner->categorical_stderr[l - num_reg_continuous] : NA_REAL);
        }
      }
    }

#ifdef MPI2
    if(call->fit_progress_active)
#endif
      np_progress_fit_loop_step(j + 1, call->fit_progress_total);
  }

  execution->status = NP_REGRESSION_GENERAL_LP_FIT_OK;
  return R_NilValue;
}
