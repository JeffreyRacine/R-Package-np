int np_regression_lp_hat_matrix(double *vector_scale_factor,
                                int deriv_var,
                                int deriv_order,
                                double *weights_out,
                                double *ridge_used_out,
                                const NPContinuousKernelRoute *kernel_route,
                                NPContinuousKernelDerivativeDiagnostics *kernel_route_diagnostics,
                                int categorical_compress,
                                const NPNNGeometryContext *nn_geometry_context,
                                NPRegressionLPEmptyRows *empty_rows){
  const int num_train = num_obs_train_extern;
  const int num_eval = num_obs_eval_extern;
  const int num_reg_tot =
    num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
  const int bw_rows =
    (BANDWIDTH_den_extern == BW_FIXED) ? 1 :
    ((BANDWIDTH_den_extern == BW_GEN_NN) ? num_eval : num_train);
  const int use_bernstein = (int_glp_bernstein_extern != 0);
  int *kernel_cx = NULL, *kernel_ux = NULL, *kernel_ox = NULL;
  int *x_operator = NULL, *derivative_mask = NULL;
  double *vsfx = NULL, *lambdax = NULL, *kw = NULL, *pkw = NULL;
  double *fallback_lower = NULL, *fallback_upper = NULL;
  double *mean_row = NULL;
  double *out2 = NULL;
  double **matrix_bandwidth_x = NULL, **matrix_bandwidth_eval_one = NULL;
  double **eval_xuno_one = NULL, **eval_xord_one = NULL, **eval_xcon_one = NULL;
  NPLPSolveWorkspace solve_workspace;
  NPLPDesignSupport design_support = {0};
  NPReghatLPWorkspace reghat_workspace;
  NPBetaScaledRowContext beta_row_context;
  double *eval_basis = NULL;
  NPConditionalBoundState bounds_state;
  int i, j, l;
  NPRegressionLPMatrixStatus status = NP_REGRESSION_LP_MATRIX_ERROR;
  NPNNGeometryStatus nn_geometry_status = NP_NN_GEOMETRY_OK;

  if(np_lp_engine_extern != NP_LP_ENGINE_SCALAR &&
     np_lp_engine_extern != NP_LP_ENGINE_GENERAL)
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;
  if((num_train <= 0) || (num_eval <= 0) || (num_reg_continuous_extern <= 0))
    return 1;
  if((deriv_var < 0) || (deriv_var > num_reg_continuous_extern) || (deriv_order < 0))
    return 1;
  if(kernel_route != NULL &&
     (np_continuous_kernel_route_validate(kernel_route,
                                          num_reg_continuous_extern) !=
        NP_CKERNEL_ROUTE_OK ||
      !np_continuous_kernel_route_has_beta(kernel_route) ||
      kernel_route->segment_count != 1 ||
      kernel_route->segment[0].coordinate_offset != 0 ||
      kernel_route->segment[0].coordinate_count !=
        num_reg_continuous_extern ||
      kernel_route_diagnostics == NULL ||
      categorical_compress < 0 || categorical_compress > 1))
    return 1;
  if(np_lp_engine_extern == NP_LP_ENGINE_SCALAR) {
    const int scalar_legacy_derivative =
      kernel_route == NULL && deriv_order == 1 && deriv_var > 0;
    const int scalar_beta_mean =
      kernel_route != NULL && deriv_order == 0 && deriv_var == 0;

    if(!scalar_legacy_derivative && !scalar_beta_mean)
      return 1;
  }

  np_lp_solve_workspace_init(&solve_workspace);
  np_reghat_lp_workspace_init(&reghat_workspace);
  np_beta_scaled_row_context_init(&beta_row_context);
  memset(weights_out, 0, (size_t)num_eval*(size_t)num_train*sizeof(double));
  if(ridge_used_out != NULL)
    memset(ridge_used_out, 0, (size_t)num_eval*sizeof(double));

  vsfx = alloc_vecd(MAX(1, num_reg_tot));
  lambdax = alloc_vecd(MAX(1, num_reg_unordered_extern + num_reg_ordered_extern));
  kw = alloc_vecd(MAX(1, num_train));
  if(np_lp_engine_extern == NP_LP_ENGINE_SCALAR &&
     kernel_route == NULL && deriv_order == 1)
    pkw = alloc_vecd(MAX(1, num_train));
  matrix_bandwidth_x = alloc_tmatd(bw_rows, num_reg_continuous_extern);
  if(kernel_route == NULL) {
    matrix_bandwidth_eval_one =
      alloc_tmatd(1, num_reg_continuous_extern);
    if(num_reg_unordered_extern > 0)
      eval_xuno_one = alloc_matd(1, num_reg_unordered_extern);
    if(num_reg_ordered_extern > 0)
      eval_xord_one = alloc_matd(1, num_reg_ordered_extern);
    if(num_reg_continuous_extern > 0)
      eval_xcon_one = alloc_matd(1, num_reg_continuous_extern);
  }

  if(kernel_route == NULL)
    kernel_cx = (int *)calloc(
      (size_t)MAX(1, num_reg_continuous_extern), sizeof(int));
  kernel_ux = (int *)calloc((size_t)MAX(1, num_reg_unordered_extern), sizeof(int));
  kernel_ox = (int *)calloc((size_t)MAX(1, num_reg_ordered_extern), sizeof(int));
  x_operator = (int *)calloc((size_t)MAX(1, num_reg_tot), sizeof(int));
  if(np_lp_engine_extern == NP_LP_ENGINE_SCALAR && kernel_route == NULL &&
     deriv_order == 1)
    derivative_mask =
      (int *)calloc((size_t)MAX(1, num_reg_tot), sizeof(int));
  if(np_lp_engine_extern == NP_LP_ENGINE_SCALAR && kernel_route == NULL &&
     deriv_order == 1 &&
     (vector_cxkerlb_extern == NULL || vector_cxkerub_extern == NULL)) {
    fallback_lower = alloc_vecd(MAX(1, num_reg_continuous_extern));
    fallback_upper = alloc_vecd(MAX(1, num_reg_continuous_extern));
    for(i = 0; i < num_reg_continuous_extern; ++i) {
      fallback_lower[i] = R_NegInf;
      fallback_upper[i] = R_PosInf;
    }
  }

  if(kernel_route == NULL)
    mean_row = alloc_vecd(MAX(1, num_train));
  if((vsfx == NULL) || (lambdax == NULL) || (kw == NULL) ||
     (np_lp_engine_extern == NP_LP_ENGINE_SCALAR &&
      kernel_route == NULL && deriv_order == 1 && pkw == NULL) ||
     (kernel_route == NULL && mean_row == NULL) ||
     ((num_reg_continuous_extern > 0) && (matrix_bandwidth_x == NULL)) ||
     ((kernel_route == NULL && num_reg_continuous_extern > 0) &&
      (matrix_bandwidth_eval_one == NULL || eval_xcon_one == NULL)) ||
     ((kernel_route == NULL && num_reg_unordered_extern > 0) &&
      eval_xuno_one == NULL) ||
     ((kernel_route == NULL && num_reg_ordered_extern > 0) &&
      eval_xord_one == NULL) ||
     (kernel_route == NULL && kernel_cx == NULL) ||
     (kernel_ux == NULL) || (kernel_ox == NULL) || (x_operator == NULL) ||
     (np_lp_engine_extern == NP_LP_ENGINE_SCALAR && kernel_route == NULL &&
      deriv_order == 1 &&
      (derivative_mask == NULL ||
       ((vector_cxkerlb_extern == NULL || vector_cxkerub_extern == NULL) &&
        (fallback_lower == NULL || fallback_upper == NULL)))))
    goto cleanup_lp_hat;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        vsfx,
                        NULL,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  if(kernel_cx != NULL)
    for(i = 0; i < num_reg_continuous_extern; i++)
      kernel_cx[i] = KERNEL_reg_extern;
  for(i = 0; i < num_reg_unordered_extern; i++) kernel_ux[i] = KERNEL_reg_unordered_extern;
  for(i = 0; i < num_reg_ordered_extern; i++) kernel_ox[i] = KERNEL_reg_ordered_extern;
  for(i = 0; i < num_reg_tot; i++) x_operator[i] = OP_NORMAL;

  if(kernel_route != NULL) {
    if(np_beta_continuous_bandwidth_prepare_canonical_ctx(
        BANDWIDTH_den_extern,
        num_train,
        num_eval,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_reg_continuous_extern,
        matrix_X_continuous_train_extern,
        matrix_X_continuous_eval_extern,
        vsfx,
        matrix_bandwidth_x,
        NULL,
        lambdax,
        NULL,
        nn_geometry_context,
        &nn_geometry_status) == 1) {
      if(nn_geometry_status == NP_NN_GEOMETRY_ZERO_RADIUS)
        status = NP_REGRESSION_LP_MATRIX_ZERO_RADIUS;
      goto cleanup_lp_hat;
    }
  } else if(kernel_bandwidth_mean_ctx(KERNEL_reg_extern,
                           BANDWIDTH_den_extern,
                           num_train,
                           num_eval,
                           0,
                           0,
                           0,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           vsfx,
                           NULL,
                           NULL,
                           matrix_X_continuous_train_extern,
                           matrix_X_continuous_eval_extern,
                           NULL,
                           matrix_bandwidth_x,
                           lambdax,
                           nn_geometry_context,
                           NULL,
                           &nn_geometry_status) == 1) {
    if(nn_geometry_status == NP_NN_GEOMETRY_ZERO_RADIUS)
      status = NP_REGRESSION_LP_MATRIX_ZERO_RADIUS;
    goto cleanup_lp_hat;
  }

  if(kernel_route != NULL &&
     np_beta_scaled_row_context_prepare(
       &beta_row_context, kernel_route, kernel_route_diagnostics,
       BANDWIDTH_den_extern, num_train, num_eval,
       num_reg_continuous_extern, num_reg_unordered_extern,
       num_reg_ordered_extern, matrix_X_continuous_train_extern,
       matrix_X_continuous_eval_extern, matrix_X_unordered_train_extern,
       matrix_X_unordered_eval_extern, matrix_X_ordered_train_extern,
       matrix_X_ordered_eval_extern, matrix_bandwidth_x,
       matrix_bandwidth_x, x_operator, kernel_ux, kernel_ox,
       (num_reg_unordered_extern + num_reg_ordered_extern > 0) ?
         lambdax : NULL,
       num_categories_extern_X, matrix_categorical_vals_extern_X,
       categorical_compress, kw) != NP_CONTINUOUS_ROW_OK)
    goto cleanup_lp_hat;

  if(np_lp_engine_extern == NP_LP_ENGINE_SCALAR) {
    if(kernel_route == NULL) {
      const int which_var = deriv_var - 1;

      derivative_mask[which_var] = 1;

      for(i = 0; i < num_eval; ++i) {
        double sum = 0.0;
        double derivative_sum = 0.0;
        double derivative_sum_output[1] = {0.0};
        int row_status;
        NPPermutationWeightOutput derivative_output =
          np_pkw_output_make(pkw, 1);

        for(l = 0; l < num_reg_unordered_extern; ++l)
          eval_xuno_one[l][0] = matrix_X_unordered_eval_extern[l][i];
        for(l = 0; l < num_reg_ordered_extern; ++l)
          eval_xord_one[l][0] = matrix_X_ordered_eval_extern[l][i];
        for(l = 0; l < num_reg_continuous_extern; ++l) {
          eval_xcon_one[l][0] = matrix_X_continuous_eval_extern[l][i];
          matrix_bandwidth_eval_one[l][0] =
            (BANDWIDTH_den_extern == BW_GEN_NN) ?
              matrix_bandwidth_x[l][i] : matrix_bandwidth_x[l][0];
        }

        np_conditional_push_bounds(int_cxker_bound_extern,
                                   vector_cxkerlb_extern == NULL ?
                                     fallback_lower : vector_cxkerlb_extern,
                                   vector_cxkerub_extern == NULL ?
                                     fallback_upper : vector_cxkerub_extern,
                                   &bounds_state);
        row_status = kernel_weighted_sum_np_ctx(
             kernel_cx, kernel_ux, kernel_ox,
             BANDWIDTH_den_extern, num_train, 1,
             num_reg_unordered_extern, num_reg_ordered_extern,
             num_reg_continuous_extern,
             0, 0, 1, 1, 1, 0, 0, 0, 0,
             x_operator, OP_DERIVATIVE,
             0, 0, derivative_mask, 1,
             0, 0,
             (BANDWIDTH_den_extern == BW_ADAP_NN) ?
               NP_TREE_FALSE : int_TREE_X,
             0,
             (BANDWIDTH_den_extern == BW_ADAP_NN) ? NULL : kdt_extern_X,
             NULL, NULL, NULL,
             matrix_X_unordered_train_extern,
             matrix_X_ordered_train_extern,
             matrix_X_continuous_train_extern,
             eval_xuno_one, eval_xord_one, eval_xcon_one,
             NULL, NULL, NULL, vsfx, 1,
             matrix_bandwidth_x, matrix_bandwidth_eval_one,
             lambdax, num_categories_extern_X,
             matrix_categorical_vals_extern_X, NULL,
             mean_row, derivative_sum_output, kw, NULL,
             &derivative_output);
        np_conditional_pop_bounds(&bounds_state);
        if(row_status != 0)
          goto cleanup_lp_hat;

        for(j = 0; j < num_train; ++j) {
          sum += kw[j];
          derivative_sum += pkw[j];
        }
        if(!R_FINITE(sum) || !R_FINITE(derivative_sum))
          goto cleanup_lp_hat;
        if(sum == 0.0) {
          if(np_lp_record_empty_row(kw, num_train, i, num_eval, num_train,
                                    weights_out, ridge_used_out, empty_rows))
            continue;
          goto cleanup_lp_hat;
        }

        for(j = 0; j < num_train; ++j) {
          const int orig_j = (int_TREE_X == NP_TREE_TRUE) ?
            ipt_extern_X[j] : j;
          weights_out[(size_t)i +
                      (size_t)num_eval*(size_t)orig_j] =
            pkw[j]/sum - kw[j]*derivative_sum/(sum*sum);
        }
      }
      status = NP_REGRESSION_LP_MATRIX_OK;
      goto cleanup_lp_hat;
    }

    for(i = 0; i < num_eval; ++i) {
      double sum = 0.0;
      if(np_beta_scaled_row_context_fill(
           &beta_row_context, i, &sum, NULL) != NP_CONTINUOUS_ROW_OK ||
         !R_FINITE(sum) || sum == 0.0)
        goto cleanup_lp_hat;
      for(j = 0; j < num_train; ++j) {
        const int orig_j = (int_TREE_X == NP_TREE_TRUE) ?
          ipt_extern_X[j] : j;
        weights_out[(size_t)i +
                    (size_t)num_eval*(size_t)orig_j] = kw[j]/sum;
      }
    }
    status = NP_REGRESSION_LP_MATRIX_OK;
    goto cleanup_lp_hat;
  }

  if((vector_glp_degree_extern == NULL) || (num_reg_continuous_extern <= 0))
    goto cleanup_lp_hat;
  if(!np_glp_cv_prepare_extern(NP_LP_ENGINE_GENERAL,
                               num_train,
                               num_reg_continuous_extern,
                               matrix_X_continuous_train_extern))
    goto cleanup_lp_hat;
  if((np_glp_cv_cache.nterms <= 0) ||
     (np_glp_cv_cache.basis == NULL) ||
     (np_glp_cv_cache.terms == NULL))
    goto cleanup_lp_hat;

  np_lp_design_support_prepare(&design_support, np_glp_cv_cache.basis,
                                num_train, np_glp_cv_cache.nterms);
  if(kernel_route != NULL) {
    if(np_reghat_lp_workspace_prepare_columns(
         &reghat_workspace, np_glp_cv_cache.basis, num_train,
         np_glp_cv_cache.nterms) != NP_REGHAT_LP_ROW_OK)
      goto cleanup_lp_hat;
  } else {
    if(!np_lp_solve_workspace_reserve(&solve_workspace,
                                      np_glp_cv_cache.nterms,
                                      1))
      goto cleanup_lp_hat;
    out2 = alloc_vecd(
      MAX(1, np_glp_cv_cache.nterms*np_glp_cv_cache.nterms));
  }
  eval_basis = alloc_vecd(MAX(1, np_glp_cv_cache.nterms));
  if((kernel_route == NULL && out2 == NULL) || (eval_basis == NULL))
    goto cleanup_lp_hat;

  for(i = 0; i < num_eval; i++){
    const int which_var = deriv_var - 1;
    const double ridge_fraction = 1.0/(double)MAX(1, num_train);
    NPLPSolvePolicyDiagnostics solve_diagnostics;
    if(out2 != NULL)
      memset(out2, 0,
             (size_t)np_glp_cv_cache.nterms*
             (size_t)np_glp_cv_cache.nterms*sizeof(double));

    if(kernel_route == NULL) {
      for(l = 0; l < num_reg_unordered_extern; l++)
        eval_xuno_one[l][0] = matrix_X_unordered_eval_extern[l][i];
      for(l = 0; l < num_reg_ordered_extern; l++)
        eval_xord_one[l][0] = matrix_X_ordered_eval_extern[l][i];
      for(l = 0; l < num_reg_continuous_extern; l++) {
        eval_xcon_one[l][0] = matrix_X_continuous_eval_extern[l][i];
        matrix_bandwidth_eval_one[l][0] =
          (BANDWIDTH_den_extern == BW_GEN_NN) ?
            matrix_bandwidth_x[l][i] : matrix_bandwidth_x[l][0];
      }
    }

    if(use_bernstein){
      if(deriv_order > 0){
        np_glp_fill_basis_eval_deriv(which_var,
                                     deriv_order,
                                     num_reg_continuous_extern,
                                     np_glp_cv_cache.terms,
                                     np_glp_cv_cache.nterms,
                                     matrix_X_continuous_eval_extern,
                                     i,
                                     np_glp_cv_cache.basis_ctx,
                                     eval_basis);
      } else {
        np_glp_fill_basis_eval(num_reg_continuous_extern,
                               np_glp_cv_cache.terms,
                               np_glp_cv_cache.nterms,
                               matrix_X_continuous_eval_extern,
                               i,
                               np_glp_cv_cache.basis_ctx,
                               eval_basis);
      }
    } else {
      if(deriv_order > 0){
        np_glp_fill_basis_eval_deriv_raw(which_var,
                                         deriv_order,
                                         num_reg_continuous_extern,
                                         np_glp_cv_cache.terms,
                                         np_glp_cv_cache.nterms,
                                         matrix_X_continuous_eval_extern,
                                         i,
                                         eval_basis);
      } else {
        np_glp_fill_basis_eval_raw(num_reg_continuous_extern,
                                   np_glp_cv_cache.terms,
                                   np_glp_cv_cache.nterms,
                                   matrix_X_continuous_eval_extern,
                                   i,
                                   eval_basis);
      }
    }

    if(kernel_route != NULL) {
      if(np_beta_scaled_row_context_fill(
           &beta_row_context, i, NULL, NULL) != NP_CONTINUOUS_ROW_OK)
        goto cleanup_lp_hat;
      if(np_reghat_lp_workspace_influence_row_ranked(
           &reghat_workspace, kw, eval_basis, weights_out + i,
           (size_t)num_eval,
           np_lp_rank_upper_bound_from_design(
             kw, num_train, np_glp_cv_cache.nterms, &design_support)) != NP_REGHAT_LP_ROW_OK) {
        if(np_lp_failed_system_is_finite(&reghat_workspace.solve_workspace,
             np_glp_cv_cache.nterms, 1) &&
           np_lp_record_empty_row(kw, num_train, i, num_eval, num_train,
             weights_out, ridge_used_out, empty_rows))
          continue;
        goto cleanup_lp_hat;
      }
      continue;
    } else {
      np_conditional_push_bounds(int_cxker_bound_extern,
                                 vector_cxkerlb_extern,
                                 vector_cxkerub_extern,
                                 &bounds_state);
      if(np_conditional_kernel_row(kernel_cx,
                                            kernel_ux,
                                            kernel_ox,
                                            x_operator,
                                            BANDWIDTH_den_extern,
                                            num_train,
                                            num_reg_unordered_extern,
                                            num_reg_ordered_extern,
                                            num_reg_continuous_extern,
                                            matrix_X_unordered_train_extern,
                                            matrix_X_ordered_train_extern,
                                            matrix_X_continuous_train_extern,
                                            eval_xuno_one,
                                            eval_xord_one,
                                            eval_xcon_one,
                                            vsfx,
                                            1,
                                            matrix_bandwidth_x,
                                            matrix_bandwidth_eval_one,
                                            lambdax,
                                            num_categories_extern_X,
                                            matrix_categorical_vals_extern_X,
                                            (BANDWIDTH_den_extern == BW_ADAP_NN) ? NP_TREE_FALSE : int_TREE_X,
                                            (BANDWIDTH_den_extern == BW_ADAP_NN) ? NULL : kdt_extern_X,
                                            kw,
                                            mean_row) != 0){
        np_conditional_pop_bounds(&bounds_state);
        goto cleanup_lp_hat;
      }

      if(kernel_weighted_sum_np_ctx(kernel_cx,
                                  kernel_ux,
                                  kernel_ox,
                                  BANDWIDTH_den_extern,
                                  num_train,
                                  1,
                                  num_reg_unordered_extern,
                                  num_reg_ordered_extern,
                                  num_reg_continuous_extern,
                                  0,
                                  0,
                                  1,
                                  1,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  x_operator,
                                  OP_NOOP,
                                  0,
                                  0,
                                  NULL,
                                  1,
                                  np_glp_cv_cache.nterms,
                                  np_glp_cv_cache.nterms,
                                  (BANDWIDTH_den_extern == BW_ADAP_NN) ? NP_TREE_FALSE : int_TREE_X,
                                  0,
                                  (BANDWIDTH_den_extern == BW_ADAP_NN) ? NULL : kdt_extern_X,
                                  NULL, NULL, NULL,
                                  matrix_X_unordered_train_extern,
                                  matrix_X_ordered_train_extern,
                                  matrix_X_continuous_train_extern,
                                  eval_xuno_one,
                                  eval_xord_one,
                                  eval_xcon_one,
                                  np_glp_cv_cache.basis,
                                  np_glp_cv_cache.basis,
                                  NULL,
                                  vsfx,
                                  1,
                                  matrix_bandwidth_x,
                                  matrix_bandwidth_eval_one,
                                  lambdax,
                                  num_categories_extern_X,
                                  matrix_categorical_vals_extern_X,
                                  NULL,
                                  out2,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL) != 0){
        np_conditional_pop_bounds(&bounds_state);
        goto cleanup_lp_hat;
      }
      np_conditional_pop_bounds(&bounds_state);
    }

    for(l = 0; l < np_glp_cv_cache.nterms; l++){
      const int base = l*np_glp_cv_cache.nterms;
      solve_workspace.rhs_source[l] = eval_basis[l];
      for(j = 0; j < np_glp_cv_cache.nterms; j++)
        solve_workspace.gram_source[l+j*np_glp_cv_cache.nterms] =
          out2[base + j];
    }

    if(np_lp_solve_workspace_solve_adjoint_ranked(
         &solve_workspace,
         np_glp_cv_cache.nterms,
         1,
         ridge_fraction,
         np_lp_rank_upper_bound_from_design(kw, num_train, np_glp_cv_cache.nterms, &design_support),
         &solve_diagnostics) !=
       NP_LP_SOLVE_POLICY_OK) {
      if(empty_rows != NULL &&
         np_lp_failed_system_is_finite(&solve_workspace,
                                       np_glp_cv_cache.nterms, 1) &&
         np_lp_record_empty_row(kw, num_train, i, num_eval, num_train,
                                weights_out, ridge_used_out, empty_rows))
        continue;
      goto cleanup_lp_hat;
    }
    if(ridge_used_out != NULL)
      ridge_used_out[i] = solve_diagnostics.ridge_total;

    for(j = 0; j < num_train; j++){
      double zju = 0.0;
      const int orig_j = (int_TREE_X == NP_TREE_TRUE) ? ipt_extern_X[j] : j;
      for(l = 0; l < np_glp_cv_cache.nterms; l++)
        zju += np_glp_cv_cache.basis[l][j]*solve_workspace.rhs_work[l];
      weights_out[(size_t)i +
                  (size_t)num_eval*(size_t)orig_j] = kw[j]*zju;
    }
  }

  status = NP_REGRESSION_LP_MATRIX_OK;

cleanup_lp_hat:
  np_lp_solve_workspace_clear(&solve_workspace);
  np_reghat_lp_workspace_clear(&reghat_workspace);
  np_beta_scaled_row_context_clear(&beta_row_context);
  if(vsfx != NULL) free(vsfx);
  if(lambdax != NULL) free(lambdax);
  if(kw != NULL) free(kw);
  if(pkw != NULL) free(pkw);
  if(fallback_lower != NULL) free(fallback_lower);
  if(fallback_upper != NULL) free(fallback_upper);
  if(mean_row != NULL) free(mean_row);
  if(out2 != NULL) free(out2);
  if(matrix_bandwidth_x != NULL) free_tmat(matrix_bandwidth_x);
  if(matrix_bandwidth_eval_one != NULL) free_tmat(matrix_bandwidth_eval_one);
  if(eval_xuno_one != NULL) free_mat(eval_xuno_one, num_reg_unordered_extern);
  if(eval_xord_one != NULL) free_mat(eval_xord_one, num_reg_ordered_extern);
  if(eval_xcon_one != NULL) free_mat(eval_xcon_one, num_reg_continuous_extern);
  if(kernel_cx != NULL) free(kernel_cx);
  if(kernel_ux != NULL) free(kernel_ux);
  if(kernel_ox != NULL) free(kernel_ox);
  if(x_operator != NULL) free(x_operator);
  if(derivative_mask != NULL) free(derivative_mask);
  if(eval_basis != NULL) free(eval_basis);
  return status;
}
