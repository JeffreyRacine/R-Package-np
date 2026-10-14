#line 11262 "np.c"
static SEXP np_regression_lp_apply_conditional_impl(SEXP txuno,
                                                     SEXP txord,
                                                     SEXP txcon,
                                                     SEXP exuno,
                                                     SEXP exord,
                                                     SEXP excon,
                                                     SEXP rhs,
                                                     SEXP rbw,
                                                     SEXP bwtype,
                                                     SEXP kernel_x,
                                                     SEXP kernel_xu,
                                                     SEXP kernel_xo,
                                                     SEXP use_tree,
                                                     SEXP glp_degree,
                                                     SEXP glp_gradient_order,
                                                     SEXP glp_bernstein,
                                                     SEXP glp_basis,
                                                     SEXP continuous_kernel_family,
                                                     SEXP continuous_kernel_order,
                                                     SEXP ckerlb,
                                                     SEXP ckerub,
                                                     SEXP categorical_compress,
                                                     SEXP return_hat,
                                                     SEXP train_is_eval,
                                                     SEXP leave_one_out,
                                                     SEXP sigtest_mode,
                                                     SEXP sigtest_coordinate,
                                                     SEXP sigtest_response_ready,
                                                     SEXP sigtest_pivotal,
                                                     SEXP sigtest_null_mean,
                                                     SEXP sigtest_residual_pool,
                                                     SEXP allow_empty_rows,
                                                     SEXP declared_support)
{
  SEXP txuno_r = R_NilValue, txord_r = R_NilValue, txcon_r = R_NilValue;
  SEXP exuno_r = R_NilValue, exord_r = R_NilValue, excon_r = R_NilValue;
  SEXP rhs_r = R_NilValue, rbw_r = R_NilValue, degree_i = R_NilValue;
  SEXP grad_i = R_NilValue, ckerlb_r = R_NilValue, ckerub_r = R_NilValue;
  SEXP sigtest_null_mean_r = R_NilValue;
  SEXP sigtest_residual_pool_r = R_NilValue;
  SEXP out = R_NilValue, ridge_used = R_NilValue;
  SEXP empty_flags = R_NilValue;
  NPRegressionLPEmptyRows empty_rows = {NULL, 0, NULL, NULL, NULL};
  NPRegressionLPEmptyRows *empty_rows_ptr = NULL;
  const int allow_empty_flag = allow_empty_rows != R_NilValue ?
    asLogical(allow_empty_rows) : FALSE;
  int nrow_txuno = 0, ncol_txuno = 0, nrow_txord = 0, ncol_txord = 0, nrow_txcon = 0, ncol_txcon = 0;
  int nrow_exuno = 0, ncol_exuno = 0, nrow_exord = 0, ncol_exord = 0, nrow_excon = 0, ncol_excon = 0;
  int nrow_rhs = 0, ncol_rhs = 0;
  int num_obs_train = 0, num_obs_eval = 0;
  int tree_flag = asLogical(use_tree);
  int categorical_compress_flag = asInteger(categorical_compress);
  int return_hat_flag = asLogical(return_hat);
  const int has_geometry_context = train_is_eval != R_NilValue;
  int train_is_eval_flag = has_geometry_context ? asLogical(train_is_eval) : FALSE;
  int leave_one_out_flag =
    leave_one_out != R_NilValue ? asLogical(leave_one_out) : FALSE;
  int sigtest_mode_flag =
    sigtest_mode != R_NilValue ? asInteger(sigtest_mode) : 0;
  int sigtest_coordinate_flag =
    sigtest_coordinate != R_NilValue ? asInteger(sigtest_coordinate) : -1;
  int sigtest_response_ready_flag =
    sigtest_response_ready != R_NilValue ?
      asInteger(sigtest_response_ready) : 0;
  int sigtest_pivotal_flag =
    sigtest_pivotal != R_NilValue ? asInteger(sigtest_pivotal) : 0;
  int lp_engine = NP_LP_ENGINE_GENERAL;
  int derivative_variable = 0;
  int derivative_order = 0;
  NPRegressionLPMatrixStatus compute_status = NP_REGRESSION_LP_MATRIX_ERROR;
  NPNNZeroRadiusInfo zero_radius_info;
  const NPNNGeometryContext *zero_radius_geometry = NULL;
  const char *failure_message = NULL;
  int int_large_sf_save = int_LARGE_SF;
  double nconfac_save = nconfac_extern;
  double ncatfac_save = ncatfac_extern;
  double *vector_continuous_stddev_save = vector_continuous_stddev_extern;
  int int_cxker_bound_save = int_cxker_bound_extern;
  double *vector_cxkerlb_save = vector_cxkerlb_extern;
  double *vector_cxkerub_save = vector_cxkerub_extern;
  double *shadow_continuous_stddev = NULL;
  double **rhs_cols = NULL;
  np_continuous_kernel_descriptor descriptor;
  NPContinuousKernelRoute beta_route;
  NPContinuousKernelDerivativeDiagnostics beta_diagnostics;
  const NPContinuousKernelRoute *active_route = NULL;
  NPContinuousKernelDerivativeDiagnostics *active_diagnostics = NULL;
  NPNNGeometryContext nn_geometry_context = {
    .mode = NP_NN_QUERY_EXTERNAL,
    .eval_to_train = NULL
  };
  const NPNNGeometryContext *nn_geometry_context_ptr = NULL;
  int nprotect = 0;
  int i;

  txuno_r = PROTECT(coerceVector(txuno, REALSXP));
  nprotect++;
  txord_r = PROTECT(coerceVector(txord, REALSXP));
  nprotect++;
  txcon_r = PROTECT(coerceVector(txcon, REALSXP));
  nprotect++;
  exuno_r = PROTECT(coerceVector(exuno, REALSXP));
  nprotect++;
  exord_r = PROTECT(coerceVector(exord, REALSXP));
  nprotect++;
  excon_r = PROTECT(coerceVector(excon, REALSXP));
  nprotect++;
  rhs_r = PROTECT(coerceVector(rhs, REALSXP));
  nprotect++;
  rbw_r = PROTECT(coerceVector(rbw, REALSXP));
  nprotect++;
  degree_i = PROTECT(coerceVector(glp_degree, INTSXP));
  nprotect++;
  grad_i = PROTECT(coerceVector(glp_gradient_order, INTSXP));
  nprotect++;
  ckerlb_r = PROTECT(coerceVector(ckerlb, REALSXP));
  nprotect++;
  ckerub_r = PROTECT(coerceVector(ckerub, REALSXP));
  nprotect++;
  if(sigtest_mode_flag != 0) {
    sigtest_null_mean_r = PROTECT(coerceVector(sigtest_null_mean, REALSXP));
    nprotect++;
    sigtest_residual_pool_r = PROTECT(
      coerceVector(sigtest_residual_pool, REALSXP));
    nprotect++;
  }

  np_matrix_dims_or_zero(txuno, &nrow_txuno, &ncol_txuno);
  np_matrix_dims_or_zero(txord, &nrow_txord, &ncol_txord);
  np_matrix_dims_or_zero(txcon, &nrow_txcon, &ncol_txcon);
  np_matrix_dims_or_zero(exuno, &nrow_exuno, &ncol_exuno);
  np_matrix_dims_or_zero(exord, &nrow_exord, &ncol_exord);
  np_matrix_dims_or_zero(excon, &nrow_excon, &ncol_excon);
  np_matrix_dims_or_zero(rhs, &nrow_rhs, &ncol_rhs);

  num_obs_train = MAX(nrow_txuno, MAX(nrow_txord, nrow_txcon));
  num_obs_eval = MAX(nrow_exuno, MAX(nrow_exord, nrow_excon));
  descriptor = np_continuous_kernel_descriptor_or_error(
    asInteger(continuous_kernel_family),
    asInteger(kernel_x),
    asInteger(continuous_kernel_order),
    "C_np_regression_lp_apply_conditional");
  if((num_obs_train <= 0) || (num_obs_eval <= 0))
    error("C_np_regression_lp_apply_conditional: zero-row inputs are not supported");
  if((nrow_txuno > 0 && nrow_txuno != num_obs_train) ||
     (nrow_txord > 0 && nrow_txord != num_obs_train) ||
     (nrow_txcon > 0 && nrow_txcon != num_obs_train))
    error("C_np_regression_lp_apply_conditional: training inputs must share the same row count");
  if((nrow_exuno > 0 && nrow_exuno != num_obs_eval) ||
     (nrow_exord > 0 && nrow_exord != num_obs_eval) ||
     (nrow_excon > 0 && nrow_excon != num_obs_eval))
    error("C_np_regression_lp_apply_conditional: evaluation inputs must share the same row count");
  if(return_hat_flag == NA_LOGICAL)
    error("C_np_regression_lp_apply_conditional: invalid return-hat flag");
  if(has_geometry_context && train_is_eval_flag == NA_LOGICAL)
    error("C_np_regression_lp_apply_conditional: invalid training-identity flag");
  if(leave_one_out_flag == NA_LOGICAL)
    error("C_np_regression_lp_apply_conditional: invalid leave-one-out flag");
  if(leave_one_out_flag && !train_is_eval_flag)
    error("C_np_regression_lp_apply_conditional: leave-one-out requires training identity");
  if(train_is_eval_flag && num_obs_eval != num_obs_train)
    error("C_np_regression_lp_apply_conditional: training identity requires equal row counts");
  if(!return_hat_flag && ((nrow_rhs <= 0) || (ncol_rhs <= 0)))
    error("C_np_regression_lp_apply_conditional: rhs must be a non-empty numeric matrix");
  if(!return_hat_flag && nrow_rhs != num_obs_train)
    error("C_np_regression_lp_apply_conditional: rhs row count must match training row count");
  if(sigtest_mode_flag != 0 &&
     (sigtest_mode_flag < 1 || sigtest_mode_flag > 3 ||
      sigtest_coordinate_flag < 0 || return_hat_flag ||
      leave_one_out_flag || !train_is_eval_flag ||
      (sigtest_response_ready_flag != 0 &&
       sigtest_response_ready_flag != 1) ||
      (sigtest_pivotal_flag != 0 && sigtest_pivotal_flag != 1) ||
      ncol_rhs > 8 || tree_flag ||
      XLENGTH(sigtest_null_mean_r) != num_obs_train ||
      XLENGTH(sigtest_residual_pool_r) != num_obs_train))
    error("C_np_regression_lp_apply_conditional: invalid npsigtest tile contract");
  if((asInteger(bwtype) != BW_FIXED) &&
     (asInteger(bwtype) != BW_GEN_NN) &&
     (asInteger(bwtype) != BW_ADAP_NN))
    error("C_np_regression_lp_apply_conditional: unsupported bandwidth type");
  if(tree_flag == NA_LOGICAL)
    error("C_np_regression_lp_apply_conditional: invalid tree flag");
  if(categorical_compress_flag != 0 && categorical_compress_flag != 1)
    error("C_np_regression_lp_apply_conditional: categorical compression must be TRUE or FALSE");
  if((int)XLENGTH(degree_i) != ncol_txcon)
    error("C_np_regression_lp_apply_conditional: glp_degree length mismatch");
  if((XLENGTH(grad_i) != 0) &&
     ((int)XLENGTH(grad_i) != ncol_txcon))
    error("C_np_regression_lp_apply_conditional: glp_gradient_order length mismatch");
  if(XLENGTH(rbw_r) <
     (R_xlen_t)ncol_txcon + (R_xlen_t)ncol_txuno + (R_xlen_t)ncol_txord)
    error("C_np_regression_lp_apply_conditional: bandwidth buffer is too short");
  lp_engine = NP_LP_ENGINE_SCALAR;
  for(i = 0; i < ncol_txcon; ++i) {
    if(INTEGER(degree_i)[i] < 0)
      error("C_np_regression_lp_apply_conditional: invalid LP degree");
    if(INTEGER(degree_i)[i] > 0)
      lp_engine = NP_LP_ENGINE_GENERAL;
  }
  for(i = 0; i < (int)XLENGTH(grad_i); ++i) {
    const int order = INTEGER(grad_i)[i];

    if(order < 0)
      error("C_np_regression_lp_apply_conditional: invalid LP derivative order");
    if(order > 0) {
      if(derivative_variable != 0)
        error("C_np_regression_lp_apply_conditional: only one LP derivative coordinate is supported");
      derivative_variable = i + 1;
      derivative_order = order;
    }
  }
  if(return_hat_flag && lp_engine == NP_LP_ENGINE_SCALAR &&
     derivative_order > 0 &&
     descriptor.family != NP_CKERNEL_FAMILY_LEGACY)
    error("C_np_regression_lp_apply_conditional: scalar derivative hats require a legacy continuous kernel");

  if(descriptor.family == NP_CKERNEL_FAMILY_BETA) {
    if(sigtest_mode_flag != 0)
      error("C_np_regression_lp_apply_conditional: npsigtest tile does not admit beta kernels");
    if(XLENGTH(ckerlb_r) != ncol_txcon ||
       XLENGTH(ckerub_r) != ncol_txcon)
      error("C_np_regression_lp_apply_conditional: beta bounds length mismatch");
    if(tree_flag)
      error("C_np_regression_lp_apply_conditional: beta route cannot use a legacy tree");
    beta_route.segment_count = 1;
    beta_route.segment[0].descriptor = descriptor;
    beta_route.segment[0].coordinate_offset = 0;
    beta_route.segment[0].coordinate_count = ncol_txcon;
    beta_route.segment[0].lower = REAL(ckerlb_r);
    beta_route.segment[0].upper = REAL(ckerub_r);
    beta_diagnostics.bad_coordinate = -1;
    beta_diagnostics.bad_observation = -1;
    beta_diagnostics.undefined_count = 0;
    beta_diagnostics.beta_status = NP_BETA_OK;
    active_route = &beta_route;
    active_diagnostics = &beta_diagnostics;
  } else if(!((XLENGTH(ckerlb_r) == 0 && XLENGTH(ckerub_r) == 0) ||
              (XLENGTH(ckerlb_r) == ncol_txcon &&
               XLENGTH(ckerub_r) == ncol_txcon))) {
    error("C_np_regression_lp_apply_conditional: legacy bounds length mismatch");
  }

  if(allow_empty_flag == NA_LOGICAL ||
     (allow_empty_flag && (!has_geometry_context || train_is_eval_flag ||
       leave_one_out_flag || sigtest_mode_flag != 0 ||
       !((descriptor.family == NP_CKERNEL_FAMILY_LEGACY &&
          ((lp_engine == NP_LP_ENGINE_GENERAL && derivative_order == 0 &&
            asInteger(bwtype) != BW_ADAP_NN) ||
           (return_hat_flag && lp_engine == NP_LP_ENGINE_SCALAR &&
            derivative_order == 1 && derivative_variable > 0))) ||
         (descriptor.family == NP_CKERNEL_FAMILY_BETA &&
          lp_engine == NP_LP_ENGINE_GENERAL && derivative_order == 1 &&
          derivative_variable > 0)))))
    error("C_np_regression_lp_apply_conditional: invalid partial-row policy");
  if(allow_empty_flag) {
    empty_flags = PROTECT(allocVector(INTSXP, num_obs_eval));
    nprotect++;
    memset(INTEGER(empty_flags), 0, (size_t)num_obs_eval*sizeof(int));
    empty_rows.flags = INTEGER(empty_flags);
    empty_rows_ptr = &empty_rows;
  }

  const int support_dimensions[4] = {0, 0, ncol_txuno, ncol_txord};
  const SEXP support_data[4] = {R_NilValue, R_NilValue, txuno_r, txord_r};
  const NPCategoricalSupport support = np_categorical_support_parse(
    declared_support, num_obs_train, support_dimensions, support_data);

  np_native_estimator_state_active = 1;
  num_obs_train_extern = num_obs_train;
  num_obs_eval_extern = num_obs_eval;
  num_var_unordered_extern = 0;
  num_var_ordered_extern = 0;
  num_var_continuous_extern = 0;
  num_reg_unordered_extern = ncol_txuno;
  num_reg_ordered_extern = ncol_txord;
  num_reg_continuous_extern = ncol_txcon;
  int_LARGE_SF = 1;
  nconfac_extern = 1.0;
  ncatfac_extern = 1.0;

  KERNEL_den_extern = CK_GAUSS2;
  KERNEL_den_unordered_extern = UKERNEL_UAA;
  KERNEL_den_ordered_extern = 0;
  KERNEL_reg_extern = asInteger(kernel_x);
  KERNEL_reg_unordered_extern = asInteger(kernel_xu);
  KERNEL_reg_ordered_extern = asInteger(kernel_xo);
  BANDWIDTH_den_extern = asInteger(bwtype);

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train, num_reg_continuous_extern);
  matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval, num_reg_unordered_extern);
  matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval, num_reg_ordered_extern);
  matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval, num_reg_continuous_extern);
  if((num_reg_unordered_extern > 0 &&
      (matrix_X_unordered_train_extern == NULL ||
       matrix_X_unordered_eval_extern == NULL)) ||
     (num_reg_ordered_extern > 0 &&
      (matrix_X_ordered_train_extern == NULL ||
       matrix_X_ordered_eval_extern == NULL)) ||
     (num_reg_continuous_extern > 0 &&
      (matrix_X_continuous_train_extern == NULL ||
       matrix_X_continuous_eval_extern == NULL))) {
    failure_message = "data-matrix allocation failed";
    goto cleanup_lp_apply_wrapper;
  }
  np_fill_column_major_matrix(matrix_X_unordered_train_extern, REAL(txuno_r), num_obs_train, num_reg_unordered_extern);
  np_fill_column_major_matrix(matrix_X_ordered_train_extern, REAL(txord_r), num_obs_train, num_reg_ordered_extern);
  np_fill_column_major_matrix(matrix_X_continuous_train_extern, REAL(txcon_r), num_obs_train, num_reg_continuous_extern);
  np_fill_column_major_matrix(matrix_X_unordered_eval_extern, REAL(exuno_r), num_obs_eval, num_reg_unordered_extern);
  np_fill_column_major_matrix(matrix_X_ordered_eval_extern, REAL(exord_r), num_obs_eval, num_reg_ordered_extern);
  np_fill_column_major_matrix(matrix_X_continuous_eval_extern, REAL(excon_r), num_obs_eval, num_reg_continuous_extern);

  if(num_reg_continuous_extern > 0){
    shadow_continuous_stddev =
      (double *)malloc((size_t)num_reg_continuous_extern*sizeof(double));
    if(shadow_continuous_stddev == NULL) {
      failure_message = "stddev allocation failed";
      goto cleanup_lp_apply_wrapper;
    }
    for(i = 0; i < num_reg_continuous_extern; i++)
      shadow_continuous_stddev[i] =
        standerrd(num_obs_train, matrix_X_continuous_train_extern[i]);
    vector_continuous_stddev_extern = shadow_continuous_stddev;
  } else {
    vector_continuous_stddev_extern = NULL;
  }

  if(tree_flag){
    int_TREE_X = (num_reg_continuous_extern > 0) ? NP_TREE_TRUE : NP_TREE_FALSE;
    int_TREE_Y = int_TREE_XY = NP_TREE_FALSE;
    if(int_TREE_X == NP_TREE_TRUE){
      ipt_extern_X = (int *)malloc((size_t)num_obs_train*sizeof(int));
      ipt_lookup_extern_X = (int *)malloc((size_t)num_obs_train*sizeof(int));
      if((ipt_extern_X == NULL) || (ipt_lookup_extern_X == NULL)) {
        failure_message = "x-tree allocation failed";
        goto cleanup_lp_apply_wrapper;
      }
      for(i = 0; i < num_obs_train; i++) ipt_extern_X[i] = i;
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train, num_reg_continuous_extern,
                   4*num_reg_continuous_extern, ipt_extern_X, &kdt_extern_X);
      for(i = 0; i < num_obs_train; i++) ipt_lookup_extern_X[ipt_extern_X[i]] = i;
      for(int j = 0; j < num_reg_unordered_extern; j++)
        for(i = 0; i < num_obs_train; i++)
          matrix_X_unordered_train_extern[j][i] = REAL(txuno_r)[j*num_obs_train + ipt_extern_X[i]];
      for(int j = 0; j < num_reg_ordered_extern; j++)
        for(i = 0; i < num_obs_train; i++)
          matrix_X_ordered_train_extern[j][i] = REAL(txord_r)[j*num_obs_train + ipt_extern_X[i]];
      for(int j = 0; j < num_reg_continuous_extern; j++)
        for(i = 0; i < num_obs_train; i++)
          matrix_X_continuous_train_extern[j][i] = REAL(txcon_r)[j*num_obs_train + ipt_extern_X[i]];
      if(BANDWIDTH_den_extern == BW_GEN_NN &&
         lp_engine == NP_LP_ENGINE_GENERAL &&
         !return_hat_flag && !leave_one_out_flag) {
        /* The GNN apply owner accumulates in tree donor order. Hat and LOO
           owners map donor identities themselves and retain the original RHS. */
        SEXP ordered_rhs = PROTECT(allocMatrix(REALSXP, num_obs_train, ncol_rhs));
        nprotect++;
        for(int j = 0; j < ncol_rhs; j++)
          for(i = 0; i < num_obs_train; i++)
            REAL(ordered_rhs)[(size_t)j*num_obs_train + i] =
              REAL(rhs_r)[(size_t)j*num_obs_train + ipt_extern_X[i]];
        rhs_r = ordered_rhs;
      }
    } else {
      ipt_extern_X = NULL;
      ipt_lookup_extern_X = NULL;
    kdt_extern_X = NULL;
  }

  } else {
    int_TREE_X = int_TREE_Y = int_TREE_XY = NP_TREE_FALSE;
    ipt_extern_X = ipt_extern_Y = ipt_extern_XY = NULL;
    ipt_lookup_extern_X = ipt_lookup_extern_Y = ipt_lookup_extern_XY = NULL;
    kdt_extern_X = kdt_extern_Y = kdt_extern_XY = NULL;
  }

  if(has_geometry_context) {
    nn_geometry_context_ptr = &nn_geometry_context;
    if(train_is_eval_flag) {
      nn_geometry_context.mode = int_TREE_X == NP_TREE_TRUE ?
        NP_NN_QUERY_TRAINING_MAP : NP_NN_QUERY_TRAINING_IDENTITY;
      nn_geometry_context.eval_to_train =
        int_TREE_X == NP_TREE_TRUE ? ipt_lookup_extern_X : NULL;
    }
  }

  num_categories_extern_X = alloc_vecu(num_reg_unordered_extern + num_reg_ordered_extern);
  matrix_categorical_vals_extern_X = alloc_matd(support.max_levels, num_reg_unordered_extern + num_reg_ordered_extern);
  if((num_reg_unordered_extern + num_reg_ordered_extern > 0) &&
     (num_categories_extern_X == NULL ||
      matrix_categorical_vals_extern_X == NULL)) {
    failure_message = "categorical-metadata allocation failed";
    goto cleanup_lp_apply_wrapper;
  }
  if(int_VERBOSE == 1 || int_DEBUG == 1)
    determine_categorical_vals(num_obs_train, 0, 0,
      num_reg_unordered_extern, num_reg_ordered_extern,
      NULL, NULL, matrix_X_unordered_train_extern,
      matrix_X_ordered_train_extern, num_categories_extern_X,
      matrix_categorical_vals_extern_X);
  np_categorical_support_copy(&support, num_categories_extern_X,
                              matrix_categorical_vals_extern_X);

  np_lp_engine_extern = lp_engine;
  vector_glp_degree_extern = INTEGER(degree_i);
  vector_glp_gradient_order_extern = (XLENGTH(grad_i) > 0) ? INTEGER(grad_i) : NULL;
  int_glp_bernstein_extern = asInteger(glp_bernstein);
  int_glp_basis_extern = asInteger(glp_basis);
  vector_cxkerlb_extern = (XLENGTH(ckerlb_r) == ncol_txcon) ?
    REAL(ckerlb_r) : NULL;
  vector_cxkerub_extern = (XLENGTH(ckerub_r) == ncol_txcon) ?
    REAL(ckerub_r) : NULL;
  int_cxker_bound_extern =
    (vector_cxkerlb_extern != NULL && vector_cxkerub_extern != NULL) ?
      np_has_finite_cker_bounds(vector_cxkerlb_extern,
                                vector_cxkerub_extern,
                                ncol_txcon) : 0;

  if(sigtest_mode_flag != 0) {
    out = PROTECT(allocVector(REALSXP, ncol_rhs));
    nprotect++;
    compute_status = np_regression_lp_sigtest_iid(
      REAL(rbw_r), REAL(rhs_r), ncol_rhs,
      REAL(sigtest_null_mean_r), REAL(sigtest_residual_pool_r),
      sigtest_mode_flag, sigtest_coordinate_flag,
      sigtest_response_ready_flag, sigtest_pivotal_flag, REAL(out));
    if(compute_status != NP_REGRESSION_LP_MATRIX_OK)
      failure_message = "npsigtest IID tile reduction failed";
  } else if(return_hat_flag) {
    out = PROTECT(allocMatrix(REALSXP, num_obs_eval, num_obs_train));
    nprotect++;
    if(descriptor.family == NP_CKERNEL_FAMILY_LEGACY &&
       (lp_engine == NP_LP_ENGINE_GENERAL || leave_one_out_flag)) {
      ridge_used = PROTECT(allocVector(REALSXP, num_obs_eval));
      nprotect++;
    }
    if(leave_one_out_flag){
      if(descriptor.family != NP_CKERNEL_FAMILY_LEGACY ||
         derivative_order != 0)
        error("C_np_regression_lp_apply_conditional: unsupported leave-one-out route");
      compute_status = np_regression_lp_leave_one_out_influence(
        REAL(rbw_r), NULL, 0, REAL(out), NULL,
        (ridge_used != R_NilValue) ? REAL(ridge_used) : NULL);
    } else {
      compute_status = np_regression_lp_hat_matrix(
        REAL(rbw_r), derivative_variable, derivative_order, REAL(out),
        (ridge_used != R_NilValue) ? REAL(ridge_used) : NULL,
        active_route, active_diagnostics, categorical_compress_flag,
        nn_geometry_context_ptr, empty_rows_ptr);
    }
    if(compute_status == NP_REGRESSION_LP_MATRIX_OK &&
       ridge_used != R_NilValue)
      setAttrib(out, install("ridge.used"), ridge_used);
    if(compute_status != NP_REGRESSION_LP_MATRIX_OK)
      failure_message = compute_status == NP_REGRESSION_LP_MATRIX_ZERO_RADIUS ?
        "nearest-neighbor bandwidth has a zero literal radius after occurrence exclusion" :
        (compute_status == NP_REGRESSION_LP_MATRIX_ZERO_MASS ?
         "leave-one-out kernel row has zero effective mass" :
         "LP hat helper failed");
  } else {
    rhs_cols = (double **)malloc((size_t)ncol_rhs*sizeof(double *));
    if(rhs_cols == NULL) {
      failure_message = "rhs column allocation failed";
      goto cleanup_lp_apply_wrapper;
    }
    for(i = 0; i < ncol_rhs; i++)
      rhs_cols[i] = REAL(rhs_r) + ((size_t)i * (size_t)nrow_rhs);

    out = PROTECT(allocMatrix(REALSXP, num_obs_eval, ncol_rhs));
    nprotect++;
    if(leave_one_out_flag){
      if(descriptor.family != NP_CKERNEL_FAMILY_LEGACY ||
         derivative_order != 0)
        error("C_np_regression_lp_apply_conditional: unsupported leave-one-out route");
      ridge_used = PROTECT(allocVector(REALSXP, num_obs_eval));
      nprotect++;
      compute_status = np_regression_lp_leave_one_out_influence(
        REAL(rbw_r), rhs_cols, ncol_rhs, NULL, REAL(out), REAL(ridge_used));
      if(compute_status == NP_REGRESSION_LP_MATRIX_OK)
        setAttrib(out, install("ridge.used"), ridge_used);
    } else {
      compute_status = np_regression_lp_apply_matrix(
        REAL(rbw_r), rhs_cols, ncol_rhs, REAL(out),
        active_route, active_diagnostics, categorical_compress_flag,
        nn_geometry_context_ptr, empty_rows_ptr);
    }
    if(compute_status != NP_REGRESSION_LP_MATRIX_OK)
      failure_message = compute_status == NP_REGRESSION_LP_MATRIX_ZERO_RADIUS ?
        "nearest-neighbor bandwidth has a zero literal radius after occurrence exclusion" :
        (compute_status == NP_REGRESSION_LP_MATRIX_ZERO_MASS ?
         "leave-one-out kernel row has zero effective mass" :
         "LP apply helper failed");
  }

cleanup_lp_apply_wrapper:

  if(compute_status == NP_REGRESSION_LP_MATRIX_ZERO_RADIUS) {
    zero_radius_geometry = nn_geometry_context_ptr;
    zero_radius_info = np_nn_zero_radius_info(
      BANDWIDTH_den_extern, num_obs_train, num_obs_eval, ncol_txcon,
      matrix_X_continuous_train_extern, matrix_X_continuous_eval_extern,
      REAL(rbw_r), descriptor.family == NP_CKERNEL_FAMILY_BETA ?
        NULL : zero_radius_geometry, 0);
  }

  if(kdt_extern_X != NULL) free_kdtree(&kdt_extern_X);
  safe_free(ipt_extern_X); safe_free(ipt_lookup_extern_X);
  if(matrix_X_unordered_train_extern != NULL) free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  if(matrix_X_ordered_train_extern != NULL) free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  if(matrix_X_continuous_train_extern != NULL) free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  if(matrix_X_unordered_eval_extern != NULL) free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
  if(matrix_X_ordered_eval_extern != NULL) free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
  if(matrix_X_continuous_eval_extern != NULL) free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  if(num_categories_extern_X != NULL) safe_free(num_categories_extern_X);
  if(matrix_categorical_vals_extern_X != NULL) free_mat(matrix_categorical_vals_extern_X, num_reg_unordered_extern + num_reg_ordered_extern);
  safe_free(shadow_continuous_stddev);
  matrix_X_unordered_train_extern = matrix_X_ordered_train_extern = matrix_X_continuous_train_extern = NULL;
  matrix_X_unordered_eval_extern = matrix_X_ordered_eval_extern = matrix_X_continuous_eval_extern = NULL;
  num_categories_extern_X = NULL;
  matrix_categorical_vals_extern_X = NULL;
  ipt_extern_X = NULL;
  ipt_lookup_extern_X = NULL;
  np_lp_engine_extern = NP_LP_ENGINE_SCALAR;
  vector_glp_degree_extern = NULL;
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = 0;
  int_glp_basis_extern = 1;
  int_TREE_X = int_TREE_Y = int_TREE_XY = NP_TREE_FALSE;
  vector_continuous_stddev_extern = vector_continuous_stddev_save;
  int_cxker_bound_extern = int_cxker_bound_save;
  vector_cxkerlb_extern = vector_cxkerlb_save;
  vector_cxkerub_extern = vector_cxkerub_save;
  int_LARGE_SF = int_large_sf_save;
  nconfac_extern = nconfac_save;
  ncatfac_extern = ncatfac_save;
  np_glp_cv_clear_extern();
  safe_free(rhs_cols);
  np_native_estimator_state_active = 0;

  if(compute_status == NP_REGRESSION_LP_MATRIX_OK &&
     failure_message == NULL && empty_rows.count > 0)
    setAttrib(out, install(".np.empty.rows"), empty_flags);
  UNPROTECT(nprotect);
  if(compute_status != NP_REGRESSION_LP_MATRIX_OK || failure_message != NULL) {
    if(compute_status == NP_REGRESSION_LP_MATRIX_ZERO_RADIUS)
      np_nn_zero_radius_error(&zero_radius_info);
    if(descriptor.family == NP_CKERNEL_FAMILY_BETA &&
       beta_diagnostics.beta_status != NP_BETA_OK)
      error("C_np_regression_lp_apply_conditional: beta row failed in continuous dimension %d: %s",
            beta_diagnostics.bad_coordinate + 1,
            np_beta_status_message(beta_diagnostics.beta_status));
    error("C_np_regression_lp_apply_conditional: %s",
          failure_message != NULL ? failure_message : "LP computation failed");
  }
  return out;
}
