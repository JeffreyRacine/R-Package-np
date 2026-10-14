int kernel_estimate_regression_categorical_tree_np(
int lp_engine,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_reg,
int num_obs_train,
int num_obs_eval,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_Y,
double *vector_Y_eval,
double *vector_scale_factor,
int *num_categories,
double ** matrix_categorical_vals,
double * const restrict mean,
double **gradient,
double * const restrict mean_stderr,
double **gradient_stderr,
double *R_squared,
double *MSE,
double *MAE,
double *MAPE,
double *CORR,
double *SIGN,
const NPContinuousKernelRoute *kernel_route,
NPContinuousKernelDerivativeDiagnostics *kernel_route_diagnostics,
int categorical_compress,
NPRegressionStandardErrorMode standard_error_mode,
int ordinary_response,
const NPContinuousPreparedBandwidthView *prepared_bandwidth,
const NPNNGeometryContext *nn_geometry_context,
const NPRegressionHC0Context *hc0_context,
NPRegressionLPEmptyRows *empty_rows,
const NPConditionalLPFirstSERequest *first_se_request,
double *conditional_variance,
NPRegressionFailure *failure,
NPLPDesignSupport *prepared_design,
const double *prepared_ann_row,
NPRegressionFitProjection *fit_projection){

  // note that mean has 2*num_obs allocated for npksum
  int i, j, l;

	double INT_KERNEL_P;					 /* Integral of K(z)^2 */
	double K_INT_KERNEL_P;				 /*  K^p */
	/* Integral of K(z-0.5)*K(z+0.5) */
	double INT_KERNEL_PM_HALF = 0.0;
  double DIFF_KER_PPM = 0.0;		 /* Difference between int K(z)^p and int K(z-.5)K(z+.5) */
  double hprod;
  double gfac = 0.0;

  double * lambda = NULL;
  double ** matrix_bandwidth = NULL;
  double ** matrix_bandwidth_deriv = NULL;
  double ** categorical_matrix_bandwidth = NULL;
  int categorical_base_requires_refit = 0;
  int * operator = NULL;
  NP_GateOverrideCtx gate_ctx_local;
  int gate_override_active = 0;
  int *ov_cont_ok = NULL, *ov_disc_uno_ok = NULL, *ov_disc_ord_ok = NULL;
  double *ov_cont_hmin = NULL, *ov_cont_k0 = NULL;
  double *ov_disc_uno_const = NULL, *ov_disc_ord_const = NULL;
  int estimation_shortcut_done = 0;
  int hc0_point_shortcut_done = 0;
  int categorical_point_invariant = 0;
  int all_large_kernel_rows = 0;
  NPRegressionHC0Context hc0_point_context;
  const NPRegressionHC0Context *effective_hc0_context = hc0_context;
  int regression_fit_status = NP_REGRESSION_FIT_OK;
  int zero_radius_external_query = 0;
  int scalar_fit_status = NP_REGRESSION_SCALAR_FIT_OK;
  int general_lp_fit_status = NP_REGRESSION_GENERAL_LP_FIT_OK;
  NPRegressionFitOwner fit_owner;

#ifdef MPI2
  int stride_e = MAX((int)ceil((double) num_obs_eval / (double) iNum_Processors),1);
  int num_obs_eval_alloc = stride_e*iNum_Processors;

#else
  int num_obs_eval_alloc = num_obs_eval;
#endif

  const int do_merr = (mean_stderr != NULL);
  const int do_grad = (gradient != NULL); 
  const int do_gerr = (gradient_stderr != NULL);
  const int lp_engine_est = lp_engine;
  if(ordinary_response != 0 && ordinary_response != 1)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 1,
      "invalid internal regression response target");
  if((lp_engine_est != NP_LP_ENGINE_SCALAR) && (lp_engine_est != NP_LP_ENGINE_GENERAL))
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 1, "invalid internal regression engine");
  if(standard_error_mode != NP_REGRESSION_STDERR_LOCAL_RESIDUAL &&
     standard_error_mode != NP_REGRESSION_STDERR_CONDITIONAL_INFLUENCE)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 2, "invalid internal regression standard-error mode");
  if(do_gerr && (!do_merr || (!do_grad && first_se_request == NULL)))
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 3, "gradient standard errors require gradients and mean standard errors");
  if(conditional_variance != NULL) {
    if(lp_engine_est != NP_LP_ENGINE_GENERAL || !do_merr ||
       standard_error_mode != NP_REGRESSION_STDERR_LOCAL_RESIDUAL ||
       num_obs_eval != 1 || hc0_context != NULL)
      NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 4, "invalid internal conditional variance owner request");
#ifdef MPI2
    if(iNum_Processors > 1 && !np_mpi_local_regression_active())
      NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 5, "conditional variance export requires a local row owner");
#endif
    conditional_variance[0] = NA_REAL;
  }
  if(first_se_request != NULL) {
    int selected = 0;
    if(lp_engine_est != NP_LP_ENGINE_GENERAL ||
       standard_error_mode != NP_REGRESSION_STDERR_LOCAL_RESIDUAL ||
       num_obs_eval != 1 || !do_merr || !do_gerr || do_grad ||
       hc0_context != NULL || num_reg_continuous <= 0 ||
       first_se_request->ncon != num_reg_continuous ||
       first_se_request->order == NULL || first_se_request->se == NULL ||
       vector_glp_degree_extern == NULL)
      NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 6, "invalid internal conditional first-SE owner request");
#ifdef MPI2
    if(iNum_Processors > 1 && !np_mpi_local_regression_active())
      NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 7, "conditional first-SE supplement requires a local row owner");
#endif
    for(l = 0; l < num_reg_continuous; ++l) {
      const int demand = first_se_request->se[l];
      if(first_se_request->order[l] < 1 || (demand != 0 && demand != 1) ||
         (demand && (first_se_request->order[l] != 1 ||
                     vector_glp_degree_extern[l] < 1)))
        NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 8, "invalid internal conditional first-SE direction");
      selected += demand;
    }
    if(selected == 0)
      NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 9, "empty internal conditional first-SE owner request");
    for(l = 0; l < num_reg_continuous + num_reg_unordered + num_reg_ordered; ++l)
      gradient_stderr[l][0] = NA_REAL;
  }
  if(lp_engine_est != NP_LP_ENGINE_SCALAR &&
     standard_error_mode == NP_REGRESSION_STDERR_CONDITIONAL_INFLUENCE)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 10, "conditional influence standard errors require the scalar regression engine");
  if(hc0_context != NULL) {
    const int residual_preparing =
      hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_PREPARING;

    if((lp_engine_est != NP_LP_ENGINE_SCALAR &&
        lp_engine_est != NP_LP_ENGINE_GENERAL) ||
       standard_error_mode != NP_REGRESSION_STDERR_LOCAL_RESIDUAL ||
       hc0_context->donor_to_canonical == NULL ||
       hc0_context->num_obs_train != num_obs_train ||
       (hc0_context->status != NP_REGRESSION_HC0_RESIDUAL_PREPARING &&
        hc0_context->status != NP_REGRESSION_HC0_RESIDUAL_READY &&
        hc0_context->status != NP_REGRESSION_HC0_RESIDUAL_ALL_ZERO) ||
       (residual_preparing && do_merr) ||
       (!residual_preparing && !do_merr) ||
       (!residual_preparing && hc0_context->scaled_residual == NULL) ||
       (!residual_preparing && !R_FINITE(hc0_context->residual_scale)) ||
       (!residual_preparing && hc0_context->residual_scale < 0.0) ||
       (!residual_preparing &&
        hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_READY &&
        hc0_context->residual_scale == 0.0) ||
       (!residual_preparing &&
        hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_ALL_ZERO &&
        hc0_context->residual_scale != 0.0) ||
       (hc0_context->point_already_computed != 0 &&
        hc0_context->point_already_computed != 1))
      NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 11, "invalid internal ordinary-regression HC0 context");
    for(i = 0; i < num_obs_train; i++) {
      if(hc0_context->donor_to_canonical[i] < 0 ||
         hc0_context->donor_to_canonical[i] >= num_obs_train ||
         (!residual_preparing &&
          !R_FINITE(hc0_context->scaled_residual[i])))
        NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 12, "invalid internal ordinary-regression HC0 donor alignment");
    }
  }
  if(kernel_route == NULL &&
     standard_error_mode == NP_REGRESSION_STDERR_CONDITIONAL_INFLUENCE &&
     num_obs_eval != 1)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 13, "legacy conditional influence standard errors require one evaluation row");
  if(prepared_bandwidth != NULL && kernel_route == NULL &&
     BANDWIDTH_reg != BW_GEN_NN && BANDWIDTH_reg != BW_ADAP_NN)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 14, "ordinary prepared regression bandwidths require NN");
  if(NP_UNLIKELY(kernel_route != NULL) &&
     lp_engine_est == NP_LP_ENGINE_SCALAR) {
    const int beta_status = np_beta_scalar_regression_fit_canonical(
      BANDWIDTH_reg, KERNEL_unordered_reg, KERNEL_ordered_reg,
      num_obs_train, num_obs_eval,
      num_reg_unordered, num_reg_ordered, num_reg_continuous,
      matrix_X_unordered_train, matrix_X_ordered_train,
      matrix_X_continuous_train,
      matrix_X_unordered_eval, matrix_X_ordered_eval,
      matrix_X_continuous_eval,
      vector_Y, vector_scale_factor,
      num_categories, matrix_categorical_vals,
      mean, gradient, mean_stderr, gradient_stderr,
      kernel_route, kernel_route_diagnostics, categorical_compress,
      standard_error_mode, ordinary_response, prepared_bandwidth,
      nn_geometry_context, hc0_context, failure,
      failure != NULL ? empty_rows : NULL);
    if(beta_status != 0)
      return beta_status;
    np_regression_fit_statistics(
      num_obs_eval, vector_Y_eval, mean,
      R_squared, MSE, MAE, MAPE, CORR, SIGN);
    return 0;
  }
  if(NP_UNLIKELY(kernel_route != NULL) &&
     (np_continuous_kernel_route_validate(kernel_route,
                                          num_reg_continuous) !=
        NP_CKERNEL_ROUTE_OK ||
      !np_continuous_kernel_route_has_beta(kernel_route) ||
      kernel_route->segment_count != 1 ||
      kernel_route->segment[0].coordinate_offset != 0 ||
      kernel_route->segment[0].coordinate_count != num_reg_continuous ||
      kernel_route_diagnostics == NULL ||
      categorical_compress < 0 || categorical_compress > 1))
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_VALIDATION, 15, "canonical beta regression route has an invalid layout");

  /* Option reads can longjmp; perform them before malloc-backed state exists. */
  np_refresh_runtime_tolerances();
  np_regression_fit_owner_init(
    &fit_owner, num_reg_continuous, &gate_ctx_local);

  {
    const int num_reg_total =
      num_reg_continuous + num_reg_unordered + num_reg_ordered;
    size_t bytes = 0;
    if(num_reg_total < 0 ||
       !np_size_array_bytes_checked((size_t)num_reg_total,
                                    sizeof(int), &bytes)) {
      regression_fit_status = NP_REGRESSION_FIT_ERR_ALLOC;
      goto finish_regression_estimation;
    }
    operator = num_reg_total > 0 ? (int *)malloc(bytes) : NULL;
    fit_owner.operator = operator;
    if(num_reg_total > 0 && operator == NULL) {
      regression_fit_status = NP_REGRESSION_FIT_ERR_ALLOC;
      goto finish_regression_estimation;
    }
  }
  for(i = 0;
      i < (num_reg_continuous+num_reg_unordered+num_reg_ordered);
      i++)
    operator[i] = OP_NORMAL;
  np_gate_ctx_clear(&gate_ctx_local);

  struct th_table * otabs = NULL;
  struct th_entry * ret = NULL;
  int ** matrix_ordered_indices = NULL;

  const int bwmdim = (BANDWIDTH_reg==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_train:1);
  const int fit_progress_total =
    (BANDWIDTH_reg == BW_ADAP_NN) ? num_obs_train : num_obs_eval;

  int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;

  kernel_c = num_reg_continuous > 0 ?
    (int *)malloc((size_t)num_reg_continuous*sizeof(int)) : NULL;
  fit_owner.kernel_c = kernel_c;
  if(num_reg_continuous > 0 && kernel_c == NULL) {
    regression_fit_status = NP_REGRESSION_FIT_ERR_ALLOC;
    goto finish_regression_estimation;
  }

  for(i = 0; i < num_reg_continuous; i++)
    kernel_c[i] = KERNEL_reg;

  kernel_u = num_reg_unordered > 0 ?
    (int *)malloc((size_t)num_reg_unordered*sizeof(int)) : NULL;
  fit_owner.kernel_u = kernel_u;
  if(num_reg_unordered > 0 && kernel_u == NULL) {
    regression_fit_status = NP_REGRESSION_FIT_ERR_ALLOC;
    goto finish_regression_estimation;
  }

  for(i = 0; i < num_reg_unordered; i++)
    kernel_u[i] = KERNEL_unordered_reg;

  kernel_o = num_reg_ordered > 0 ?
    (int *)malloc((size_t)num_reg_ordered*sizeof(int)) : NULL;
  fit_owner.kernel_o = kernel_o;
  if(num_reg_ordered > 0 && kernel_o == NULL) {
    regression_fit_status = NP_REGRESSION_FIT_ERR_ALLOC;
    goto finish_regression_estimation;
  }

  for(i = 0; i < num_reg_ordered; i++)
    kernel_o[i] = KERNEL_ordered_reg;


  // assert(BANDWIDTH_reg == BW_FIXED);

  // Allocate memory for objects 

  {
    const int num_reg_categorical = num_reg_unordered + num_reg_ordered;
    lambda = num_reg_categorical > 0 ?
      (double *)malloc((size_t)num_reg_categorical*sizeof(double)) : NULL;
    fit_owner.lambda = lambda;
    matrix_bandwidth = np_regression_fit_matrix_alloc(
      bwmdim, num_reg_continuous, 1);
    fit_owner.matrix_bandwidth = matrix_bandwidth;
    matrix_bandwidth_deriv = np_regression_fit_matrix_alloc(
      bwmdim, num_reg_continuous, 0);
    fit_owner.matrix_bandwidth_deriv = matrix_bandwidth_deriv;
    if((num_reg_categorical > 0 && lambda == NULL) ||
       (num_reg_continuous > 0 &&
        (matrix_bandwidth == NULL || matrix_bandwidth_deriv == NULL))) {
      regression_fit_status = NP_REGRESSION_FIT_ERR_ALLOC;
      goto finish_regression_estimation;
    }
  }

  {
    NPRegressionFitBandwidthCall bandwidth_call = {
      .use_canonical_beta = kernel_route != NULL,
      .kernel = KERNEL_reg,
      .bandwidth_mode = BANDWIDTH_reg,
      .num_obs_train = num_obs_train,
      .num_obs_eval = num_obs_eval,
      .num_reg_unordered = num_reg_unordered,
      .num_reg_ordered = num_reg_ordered,
      .num_reg_continuous = num_reg_continuous,
      .vector_scale_factor = vector_scale_factor,
      .matrix_X_continuous_train = matrix_X_continuous_train,
      .matrix_X_continuous_eval = matrix_X_continuous_eval,
      .matrix_bandwidth = matrix_bandwidth,
      .matrix_bandwidth_deriv = matrix_bandwidth_deriv,
      .lambda = lambda,
      .prepared_bandwidth = prepared_bandwidth,
      .kernel_integral = &INT_KERNEL_P,
      .kernel_power_integral = &K_INT_KERNEL_P,
      .kernel_half_integral = &INT_KERNEL_PM_HALF,
      .kernel_difference = &DIFF_KER_PPM,
      .nn_geometry_context = nn_geometry_context,
      .nn_geometry_status = NP_NN_GEOMETRY_OK,
      .status = 1
    };
    if(np_regression_fit_bandwidth_prepare_owned(
         &bandwidth_call, &fit_owner) != 0) {
      regression_fit_status =
        bandwidth_call.nn_geometry_status == NP_NN_GEOMETRY_ZERO_RADIUS ?
        NP_REGRESSION_FIT_ERR_ZERO_NN_RADIUS :
        NP_REGRESSION_FIT_ERR_BANDWIDTH;
      goto finish_regression_estimation;
    }
  }

  categorical_matrix_bandwidth = matrix_bandwidth;
  /* Identity selects the radius, not which responses enter the fit. Do not
   * replace ordinary fitted-row geometry with an external-query radius. */
  if(!ordinary_response && lp_engine_est == NP_LP_ENGINE_GENERAL &&
     kernel_route == NULL && BANDWIDTH_reg == BW_GEN_NN && do_grad &&
     (num_reg_unordered + num_reg_ordered > 0) &&
     nn_geometry_context != NULL &&
     nn_geometry_context->mode != NP_NN_QUERY_EXTERNAL) {
    const NPNNGeometryContext external_geometry = {
      .mode = NP_NN_QUERY_EXTERNAL,
      .eval_to_train = NULL
    };
    NPNNGeometryStatus external_geometry_status = NP_NN_GEOMETRY_OK;

    if(kernel_bandwidth_mean_ctx(
         KERNEL_reg,
         BW_GEN_NN,
         num_obs_train,
         num_obs_eval,
         0, 0, 0,
         num_reg_continuous,
         num_reg_unordered,
         num_reg_ordered,
         0,
         vector_scale_factor,
         NULL,
         NULL,
         matrix_X_continuous_train,
         matrix_X_continuous_eval,
         NULL,
         matrix_bandwidth_deriv,
         lambda,
         &external_geometry,
         NULL,
         &external_geometry_status) != 0) {
      regression_fit_status =
        external_geometry_status == NP_NN_GEOMETRY_ZERO_RADIUS ?
        NP_REGRESSION_FIT_ERR_ZERO_NN_RADIUS :
        NP_REGRESSION_FIT_ERR_BANDWIDTH;
      zero_radius_external_query = 1;
      goto finish_regression_estimation;
    }
    categorical_matrix_bandwidth = matrix_bandwidth_deriv;
    categorical_base_requires_refit = 1;
  }

  hprod = 1.0;
  for(l = 0; l < num_reg_continuous; l++)
    hprod *= matrix_bandwidth[l][0];

  gfac = sqrt(DIFF_KER_PPM/K_INT_KERNEL_P);

  // compute hash stuff here if necessary

  if(do_grad && (num_reg_ordered > 0)){
    size_t index_count = 0;
    size_t index_bytes = 0;
    int *tc = NULL;
    otabs = (struct th_table *)calloc(
      (size_t)num_reg_ordered, sizeof(struct th_table));
    matrix_ordered_indices = (int **)malloc(
      (size_t)num_reg_ordered*sizeof(int *));
    if(!np_size_mul_checked((size_t)num_reg_ordered,
                            (size_t)num_obs_eval,
                            &index_count) ||
       !np_size_array_bytes_checked(index_count, sizeof(int), &index_bytes)) {
      regression_fit_status = NP_REGRESSION_FIT_ERR_ALLOC;
      goto finish_regression_estimation;
    }
    tc = index_count > 0 ? (int *)malloc(index_bytes) : NULL;
    fit_owner.ordered_tables = otabs;
    fit_owner.matrix_ordered_indices = matrix_ordered_indices;
    fit_owner.ordered_indices_storage = tc;
    if(otabs == NULL || matrix_ordered_indices == NULL || tc == NULL) {
      regression_fit_status = NP_REGRESSION_FIT_ERR_ALLOC;
      goto finish_regression_estimation;
    }
    for(l = 0; l < num_reg_ordered; l++)
      matrix_ordered_indices[l] = tc + l*num_obs_eval;

    for(l = 0; l < num_reg_ordered; l++){
      if(thcreate_r((size_t)ceil(1.6*num_categories[l+num_reg_unordered]),
                    otabs + l) == TH_ERROR) {
        regression_fit_status = NP_REGRESSION_FIT_ERR_HASH_CREATE;
        goto finish_regression_estimation;
      }
      fit_owner.ordered_tables_initialized = l + 1;

      for(i = 0; i < num_categories[l+num_reg_unordered]; i++){
        struct th_entry centry;
        centry.key.dkey = matrix_categorical_vals[l+num_reg_unordered][i];
        centry.data = i;

        if(thsearch_r(&centry, TH_ENTER, &ret, otabs+l) == TH_FAILURE) {
          regression_fit_status = NP_REGRESSION_FIT_ERR_HASH_INSERT;
          goto finish_regression_estimation;
        }
      }

      // now do lookups
      struct th_entry te;
      te.key.dkey = 0.0;
      te.data = -1;
      ret = &te;

      for(i = 0; i < num_obs_eval; i++){
        /* The placeholder key is not a cached category until the first lookup. */
        if(i == 0 || ret->key.dkey != matrix_X_ordered_eval[l][i]){
          te.key.dkey = matrix_X_ordered_eval[l][i];
          if(thsearch_r(&te, TH_SEARCH, &ret, otabs+l) == TH_FAILURE) {
            regression_fit_status = NP_REGRESSION_FIT_ERR_HASH_LOOKUP;
            goto finish_regression_estimation;
          }
        } 

        matrix_ordered_indices[l][i] = ret->data;
      }
    }
  }

  if((num_reg_continuous + num_reg_unordered + num_reg_ordered) > 0){
    int ok_all = 1;

    if(num_reg_continuous > 0){
      const double rel_tol = np_cont_largeh_rel_tol();
      if(np_cont_largeh_cache_get_or_build(num_obs_train,
                                           num_obs_eval,
                                           num_reg_continuous,
                                           kernel_c,
                                           matrix_X_continuous_train,
                                           matrix_X_continuous_eval,
                                           rel_tol,
                                           &ov_cont_ok,
                                           &ov_cont_hmin,
                                           &ov_cont_k0)){
        fit_owner.ov_cont_from_cache = 1;
        fit_owner.ov_cont_ok = ov_cont_ok;
        fit_owner.ov_cont_hmin = ov_cont_hmin;
        fit_owner.ov_cont_k0 = ov_cont_k0;
      } else {
        ov_cont_ok = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
        ov_cont_hmin = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
        ov_cont_k0 = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
        fit_owner.ov_cont_ok = ov_cont_ok;
        fit_owner.ov_cont_hmin = ov_cont_hmin;
        fit_owner.ov_cont_k0 = ov_cont_k0;
        ok_all = (ov_cont_ok != NULL) && (ov_cont_hmin != NULL) && (ov_cont_k0 != NULL);
        if(ok_all){
          for(i = 0; i < num_reg_continuous; i++){
            const int kern = kernel_c[i];
            double xmin = DBL_MAX, xmax = -DBL_MAX;
            ov_cont_ok[i] = 0; ov_cont_hmin[i] = DBL_MAX; ov_cont_k0[i] = 0.0;
            if(!np_cont_largeh_kernel_supported(kern)) continue;
            for(j = 0; j < num_obs_train; j++){
              const double v = matrix_X_continuous_train[i][j];
              if(!isfinite(v)) continue;
              xmin = MIN(xmin, v); xmax = MAX(xmax, v);
            }
            if(xmax >= xmin){
              const double utol = np_cont_largeh_utol(kern, rel_tol);
              if(utol > 0.0 && isfinite(utol)){
                ov_cont_ok[i] = 1;
                ov_cont_hmin[i] = (xmax - xmin)/utol;
                ov_cont_k0[i] = np_cont_largeh_k0(kern);
              }
            }
          }
        }
      }
    }

    if(ok_all && num_reg_unordered > 0){
      ov_disc_uno_ok = (int *)calloc((size_t)num_reg_unordered, sizeof(int));
      ov_disc_uno_const = (double *)malloc((size_t)num_reg_unordered*sizeof(double));
      fit_owner.ov_disc_uno_ok = ov_disc_uno_ok;
      fit_owner.ov_disc_uno_const = ov_disc_uno_const;
      ok_all = (ov_disc_uno_ok != NULL) && (ov_disc_uno_const != NULL);
      if(ok_all){
        double (* const ukf[])(int, double, int) = {
          np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
          np_score_uaa, np_score_unli_racine
        };
        const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
        for(i = 0; i < num_reg_unordered; i++){
          const int ku = kernel_u[i];
          const int ncat = (num_categories != NULL) ? num_categories[i] : 0;
          const double lam = lambda[i];
          ov_disc_uno_ok[i] = 0; ov_disc_uno_const[i] = 0.0;
          if(ku < 0 || ku >= nuk) continue;
          if(!NP_CF316_DISC_UPPER(ku, lam, ncat)) continue;
          {
            const double ks = ukf[ku](1, lam, ncat);
            const double kd = ukf[ku](0, lam, ncat);
            if(NP_CF316_DISC_CONSTANT(ks, kd)){
              ov_disc_uno_ok[i] = 1;
              ov_disc_uno_const[i] = 0.5*(ks + kd);
            }
          }
        }
      }
    }

    if(ok_all && num_reg_ordered > 0){
      ov_disc_ord_ok = (int *)calloc((size_t)num_reg_ordered, sizeof(int));
      ov_disc_ord_const = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
      fit_owner.ov_disc_ord_ok = ov_disc_ord_ok;
      fit_owner.ov_disc_ord_const = ov_disc_ord_const;
      ok_all = (ov_disc_ord_ok != NULL) && (ov_disc_ord_const != NULL);
      if(ok_all){

        const int nok = 16;
        for(i = 0; i < num_reg_ordered; i++){
          const int oi = i + num_reg_unordered;
          const int ko = kernel_o[i];
          const int ncat = (num_categories != NULL) ? num_categories[oi] : 0;
          const double lam = lambda[oi];
          ov_disc_ord_ok[i] = 0; ov_disc_ord_const[i] = 0.0;
          if(ko < 0 || ko >= nok) continue;
          if(ncat <= 0 || matrix_categorical_vals == NULL) continue;
          if(!NP_CF316_DISC_ORDERED_UPPER(ko, lam)) continue;
          {
            double k0, k1;
        np_ordered_kernel_range(ko,lam,matrix_categorical_vals[oi],ncat,&k0,&k1);
            if(NP_CF316_DISC_CONSTANT(k0, k1)){
              ov_disc_ord_ok[i] = 1;
              ov_disc_ord_const[i] = 0.5*(k0 + k1);
            }
          }
        }
      }
    }

    if(ok_all){
      np_gate_ctx_set(&gate_ctx_local,
                      num_reg_continuous,
                      num_reg_unordered,
                      num_reg_ordered,
                      kernel_c,
                      kernel_u,
                      kernel_o,
                      operator,
                      ov_cont_ok,
                      ov_cont_hmin,
                      ov_cont_k0,
                      ov_disc_uno_ok,
                      ov_disc_uno_const,
                      ov_disc_ord_ok,
                      ov_disc_ord_const);
      gate_override_active = 1;
      fit_owner.gate_override_active = 1;
    }
  }

  {
    int all_large_gate = (BANDWIDTH_reg == BW_FIXED) && gate_override_active;
    if(all_large_gate){
      for(i = 0; i < num_reg_continuous; i++){
        const double h = matrix_bandwidth[i][0];
        if((ov_cont_ok == NULL) || (!ov_cont_ok[i]) || (!isfinite(h)) ||
           (fabs(h) < ov_cont_hmin[i])){
          all_large_gate = 0;
          break;
        }
      }
    }
    if(all_large_gate){
      for(i = 0; i < num_reg_unordered; i++){
        if((ov_disc_uno_ok == NULL) || (!ov_disc_uno_ok[i])){
          all_large_gate = 0;
          break;
        }
      }
    }
    if(all_large_gate){
      for(i = 0; i < num_reg_ordered; i++){
        if((ov_disc_ord_ok == NULL) || (!ov_disc_ord_ok[i])){
          all_large_gate = 0;
          break;
        }
      }
    }

    all_large_kernel_rows =
      all_large_gate && lp_engine_est == NP_LP_ENGINE_GENERAL;
    categorical_point_invariant = all_large_kernel_rows;

    if(all_large_gate &&
       (hc0_context == NULL ||
        hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_PREPARING ||
        (lp_engine_est == NP_LP_ENGINE_GENERAL &&
         !hc0_context->point_already_computed &&
         (hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_READY ||
          hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_ALL_ZERO))) &&
       standard_error_mode == NP_REGRESSION_STDERR_LOCAL_RESIDUAL &&
       ((lp_engine_est == NP_LP_ENGINE_SCALAR) || (lp_engine_est == NP_LP_ENGINE_GENERAL))){
      const int point_only_hc0_shortcut =
        hc0_context != NULL && lp_engine_est == NP_LP_ENGINE_GENERAL &&
        !hc0_context->point_already_computed &&
        (hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_READY ||
         hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_ALL_ZERO);
      const int shortcut_do_merr = do_merr && !point_only_hc0_shortcut;
      const int shortcut_do_gerr = do_gerr && !point_only_hc0_shortcut;
      double kconst = 1.0;
      int kconst_ok = 1;
      const double ridge_eps = 1.0/(double)MAX(1, num_obs_train);
      double sigma2hat = 0.0;
      const double ymean = meand(num_obs_train, vector_Y);

      if(shortcut_do_merr) {
        for(i = 0; i < num_obs_train; i++){
          const double dy = vector_Y[i] - ymean;
          sigma2hat += dy*dy;
        }
        sigma2hat /= (double)MAX(1, num_obs_train);
        if(conditional_variance != NULL)
          conditional_variance[0] = sigma2hat;
      }

      for(i = 0; i < num_reg_continuous; i++){
        const double h = matrix_bandwidth[i][0];
        if(!isfinite(h) || (h == 0.0)){
          kconst_ok = 0;
          break;
        }
        kconst *= ov_cont_k0[i]/h;
      }
      if(kconst_ok){
        for(i = 0; i < num_reg_unordered; i++)
          kconst *= ov_disc_uno_const[i];
      }
      if(kconst_ok){
        for(i = 0; i < num_reg_ordered; i++)
          kconst *= ov_disc_ord_const[i];
      }
      if(!isfinite(kconst) || (kconst <= 0.0))
        kconst_ok = 0;

      if(kconst_ok && (lp_engine_est == NP_LP_ENGINE_SCALAR) && (!do_grad)){
        const double sk = ((double)num_obs_train)*kconst;
        const double sefac = (sk*hprod > 0.0) ?
          sqrt(MAX(0.0, sigma2hat) * K_INT_KERNEL_P / (sk*hprod)) : 0.0;

        if(hc0_context != NULL &&
           hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_PREPARING) {
          NPRegressionResidualPreparation * const preparation =
            hc0_context->preparation;
          const long double anchor = vector_Y[0];
          long double centered_mean = 0.0L;
          for(i = 0; i < num_obs_train; ++i)
            centered_mean += (long double)vector_Y[i] - anchor;
          centered_mean /= num_obs_train;
          const long double norm = sqrtl(
            (long double)(num_obs_train - 1) / num_obs_train);
          for(i = 0; i < num_obs_eval; ++i) {
            const int self = preparation->evaluation_to_donor[i];
            preparation->information[self] = num_obs_train > 1 ?
              NP_RESIDUAL_IDENTIFIED : NP_RESIDUAL_UNIDENTIFIED;
            preparation->normalized_residual[self] = num_obs_train > 1 ?
              (double)(((long double)vector_Y[self] - anchor - centered_mean) /
                       norm) : 0.0;
          }
        }

        for(i = 0; i < num_obs_eval; i++){
          mean[i] = ymean;
          if(shortcut_do_merr)
            mean_stderr[i] = sefac;
          np_progress_fit_loop_step_owned(
            i + 1,
            fit_progress_total,
            np_regression_fit_owner_cleanup,
            (void *)&fit_owner);
        }

        estimation_shortcut_done = 1;
        hc0_point_shortcut_done = point_only_hc0_shortcut;
      } else if(kconst_ok && lp_engine_est == NP_LP_ENGINE_GENERAL &&
                (vector_glp_degree_extern != NULL) && (num_reg_continuous > 0)){
        const int use_bernstein = (int_glp_bernstein_extern != 0);
        int *glp_terms = NULL;
        int glp_nterms = 0;
        double **basis = NULL;
        NPGLPBasisCtx *basis_ctx = NULL;
        double *eval_basis = NULL, *eval_deriv = NULL;
        NPLPFullRowWorkspace inverse_workspace;
        double *beta = NULL;
        double *eval_basis_block = NULL;
        double *projection_block = NULL;
        double *yhat_block = NULL;
        int eval_block_rows = 0;
        int basis_is_contiguous = 0;
        int use_fit_projection_blas = 0;
        int covariance_ridged = 0;
        double *covariance_workspace = NULL;
        double *covariance_work = NULL;
        const int residual_preparing = hc0_context != NULL &&
          hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_PREPARING;
        double *residual_gram = NULL;
        NPAllLargeResidualWorkspace *residual_workspace = NULL;
        NPContrastNumber residual_intercept_ridge = {0};
        int residual_unregularized = 1;
        int fast_ok;
        NPRegressionAllLargeLPFitOwner all_large_owner = {
          .num_reg_continuous = num_reg_continuous,
          .terms = NULL,
          .nterms = 0,
          .basis = NULL,
          .basis_is_contiguous = 0,
          .basis_context = NULL,
          .eval_basis = NULL,
          .eval_derivative = NULL,
          .coefficient = NULL,
          .eval_basis_block = NULL,
          .projection_block = NULL,
          .fitted_block = NULL,
          .covariance_workspace = NULL,
          .residual_gram = NULL,
          .residual_workspace = NULL,
          .inverse_workspace = &inverse_workspace,
          .enclosing_owner = &fit_owner
        };

        np_lp_full_row_workspace_init(&inverse_workspace);
        fast_ok = np_glp_build_terms(num_reg_continuous,
                                     vector_glp_degree_extern,
                                     int_glp_basis_extern,
                                     &glp_terms,
                                     &glp_nterms);
        if(fast_ok && (glp_nterms > 0)){
          basis_is_contiguous =
            (!do_grad) &&
            np_reg_alllarge_blas_profitable(num_obs_train,
                                            glp_nterms,
                                            num_obs_train) &&
            np_reg_alllarge_blas_profitable(num_obs_eval,
                                            glp_nterms,
                                            num_obs_eval);
          basis = np_regression_fit_matrix_alloc(
            num_obs_train, glp_nterms, basis_is_contiguous);
          fast_ok = np_lp_full_row_workspace_reserve(
            &inverse_workspace, glp_nterms, 1);
          beta = (double *)malloc((size_t)glp_nterms*sizeof(double));
          eval_basis = (double *)malloc((size_t)glp_nterms*sizeof(double));
          if(do_grad || first_se_request != NULL)
            eval_deriv = (double *)malloc((size_t)glp_nterms*sizeof(double));
          if(use_bernstein)
            basis_ctx = (NPGLPBasisCtx *)calloc(
              (size_t)num_reg_continuous, sizeof(NPGLPBasisCtx));
          fast_ok = fast_ok && (basis != NULL) && (beta != NULL) &&
            (eval_basis != NULL) &&
            ((!do_grad && first_se_request == NULL) || (eval_deriv != NULL)) &&
            (!use_bernstein || (basis_ctx != NULL));
          all_large_owner.terms = glp_terms;
          all_large_owner.nterms = glp_nterms;
          all_large_owner.basis = basis;
          all_large_owner.basis_is_contiguous = basis_is_contiguous;
          all_large_owner.basis_context = basis_ctx;
          all_large_owner.eval_basis = eval_basis;
          all_large_owner.eval_derivative = eval_deriv;
          all_large_owner.coefficient = beta;
        } else {
          fast_ok = 0;
        }

        if(fast_ok){
          for(l = 0; l < num_reg_continuous; l++){
            if(use_bernstein){
              double xmin = matrix_X_continuous_train[l][0];
              double xmax = matrix_X_continuous_train[l][0];
              for(i = 1; i < num_obs_train; i++){
                const double xi = matrix_X_continuous_train[l][i];
                if(xi < xmin) xmin = xi;
                if(xi > xmax) xmax = xi;
              }
              if(!np_glp_basis_ctx_init(
                  &basis_ctx[l], vector_glp_degree_extern[l], xmin, xmax,
                  int_glp_basis_extern == 1)){
                fast_ok = 0;
                break;
              }
            }
          }
        }

        if(fast_ok){
          if(use_bernstein){
            const NPRegressionAllLargeBernsteinBasisCall basis_call = {
              .num_reg_continuous = num_reg_continuous,
              .terms = glp_terms,
              .nterms = glp_nterms,
              .matrix_X_continuous_train = matrix_X_continuous_train,
              .num_obs_train = num_obs_train,
              .basis_context = basis_ctx,
              .basis = basis
            };
            R_UnwindProtect(
              np_regression_alllarge_bernstein_basis_execute,
              (void *)&basis_call,
              np_regression_alllarge_lp_fit_cleanup,
              (void *)&all_large_owner,
              NULL);
          } else {
            np_glp_fill_basis_raw_train(num_reg_continuous,
                                        glp_terms,
                                        glp_nterms,
                                        matrix_X_continuous_train,
                                        num_obs_train,
                                        basis);
          }
        }

        if(fast_ok){
          const double sk = ((double)num_obs_train)*kconst;
          const double se_default = (sk*hprod > 0.0) ?
            sqrt(MAX(0.0, sigma2hat) * K_INT_KERNEL_P / (sk*hprod)) : 0.0;

          use_fit_projection_blas = basis_is_contiguous;
          if(use_fit_projection_blas){
            const size_t projection_max_elements =
              NP_REG_ALLLARGE_PROJECTION_MAX_BYTES/sizeof(double);
            const size_t block_rows =
              MIN((size_t)NP_REG_ALLLARGE_FIT_BLOCK_MAX_ROWS,
                  MIN((size_t)num_obs_eval,
                      projection_max_elements/(size_t)glp_nterms));

            if((block_rows == 0) || (block_rows > (size_t)INT_MAX)){
              use_fit_projection_blas = 0;
            } else {
              eval_block_rows = (int)block_rows;
              eval_basis_block = (double *)malloc(
                block_rows*(size_t)glp_nterms*sizeof(double));
              if(shortcut_do_merr)
                projection_block = (double *)malloc(
                  block_rows*(size_t)glp_nterms*sizeof(double));
              yhat_block = (double *)malloc(block_rows*sizeof(double));
              if((eval_basis_block == NULL) ||
                 (shortcut_do_merr && projection_block == NULL) ||
                 (yhat_block == NULL)){
                free(eval_basis_block);
                free(projection_block);
                free(yhat_block);
                eval_basis_block = NULL;
                projection_block = NULL;
                yhat_block = NULL;
                eval_block_rows = 0;
                use_fit_projection_blas = 0;
              } else {
                all_large_owner.eval_basis_block = eval_basis_block;
                all_large_owner.projection_block = projection_block;
                all_large_owner.fitted_block = yhat_block;
              }
            }
          }

          if(use_fit_projection_blas){
            np_blas_gram_int(num_obs_train,
                             glp_nterms,
                             basis[0],
                             num_obs_train,
                             inverse_workspace.matrix_copy);
            np_blas_dgemv_t_int(num_obs_train,
                                glp_nterms,
                                basis[0],
                                num_obs_train,
                                vector_Y,
                                inverse_workspace.rhs);
          } else {
            for(i = 0; i < glp_nterms; i++){
              inverse_workspace.rhs[i] = 0.0;
              for(j = 0; j < glp_nterms; j++)
                inverse_workspace.matrix_copy[i + j*glp_nterms] = 0.0;
            }

            for(i = 0; i < num_obs_train; i++){
              const double yi = vector_Y[i];
              for(int a = 0; a < glp_nterms; a++){
                const double za = basis[a][i];
                inverse_workspace.rhs[a] += za*yi;
                for(int b = a; b < glp_nterms; b++){
                  const double zb = basis[b][i];
                  inverse_workspace.matrix_copy[a + b*glp_nterms] += za*zb;
                  if(b != a)
                    inverse_workspace.matrix_copy[b + a*glp_nterms] += za*zb;
                }
              }
            }
          }

          /* eval_basis is not yet an evaluation vector. The inverse retry
           * changes only source diagonals, before packing destroys A0. */
          if(shortcut_do_merr)
            for(i = 0; i < glp_nterms; ++i)
              eval_basis[i] = inverse_workspace.matrix_copy[i + i*glp_nterms];
          if(residual_preparing) {
            size_t elements, gram_bytes;
            if(np_size_mul_checked((size_t)glp_nterms, (size_t)glp_nterms,
                                   &elements) &&
               np_size_mul_checked(elements, sizeof(double), &gram_bytes)) {
              residual_gram = (double *)malloc(gram_bytes);
              residual_workspace = (NPAllLargeResidualWorkspace *)calloc(
                1, sizeof(NPAllLargeResidualWorkspace));
              all_large_owner.residual_gram = residual_gram;
              all_large_owner.residual_workspace = residual_workspace;
              if(residual_gram != NULL)
                memcpy(residual_gram, inverse_workspace.matrix_copy, gram_bytes);
            }
            /* An unavailable preparation remains INVALID. Allocation failure
             * does not change fast_ok, point ownership or the BLAS decision. */
          }
          fast_ok = np_lp_full_row_workspace_invert_retryable(
            &inverse_workspace, glp_nterms, ridge_eps, 64);
          if(fast_ok && residual_gram != NULL) {
            residual_intercept_ridge = np_contrast_two_sum(
              inverse_workspace.matrix_copy[0], -residual_gram[0]);
            for(i = 0; i < glp_nterms; ++i)
              if(inverse_workspace.matrix_copy[i + i*glp_nterms] !=
                 residual_gram[i + i*glp_nterms])
                residual_unregularized = 0;
          }
          if(fast_ok && shortcut_do_merr) {
            for(i = 0; i < glp_nterms; ++i)
              if(inverse_workspace.matrix_copy[i + i*glp_nterms] != eval_basis[i])
                covariance_ridged = 1;
            if(covariance_ridged) {
              size_t gram_elements, work_elements, work_bytes;
              if(np_size_mul_checked((size_t)glp_nterms, (size_t)glp_nterms,
                                     &gram_elements) &&
                 np_size_add_checked(gram_elements, (size_t)glp_nterms,
                                     &work_elements) &&
                 np_size_mul_checked(work_elements, sizeof(double), &work_bytes)) {
                covariance_workspace = (double *)malloc(work_bytes);
                all_large_owner.covariance_workspace = covariance_workspace;
                if(covariance_workspace != NULL) {
                  memcpy(covariance_workspace, inverse_workspace.matrix_copy,
                         gram_elements*sizeof(double));
                  for(i = 0; i < glp_nterms; ++i)
                    covariance_workspace[i + i*glp_nterms] = eval_basis[i];
                  covariance_work = covariance_workspace + gram_elements;
                }
              }
              /* Unavailable covariance never changes fast_ok or its point owner. */
            }
          }
          if(fast_ok)
            fast_ok = np_lp_full_row_workspace_pack_inverse_rows(
              &inverse_workspace, glp_nterms);

          if(fast_ok){
            if(use_fit_projection_blas){
              np_blas_dgemv_n_int(glp_nterms,
                                  glp_nterms,
                                  inverse_workspace.gram,
                                  glp_nterms,
                                  inverse_workspace.rhs,
                                  beta);
            } else {
              for(i = 0; i < glp_nterms; i++){
                double s = 0.0;
                for(j = 0; j < glp_nterms; j++)
                  s += inverse_workspace.matrix_copy[i*glp_nterms + j] *
                    inverse_workspace.rhs[j];
                beta[i] = s;
              }
            }

            if(fit_projection != NULL) {
              np_fit_projection_basis(fit_projection, basis, glp_nterms);
              fit_projection->constant = 1;
            }

            if(residual_preparing && residual_gram != NULL &&
               residual_workspace != NULL) {
              const NPRegressionAllLargeResidualCall residual_call = {
                .output = hc0_context->preparation,
                .nterms = glp_nterms,
                .basis = basis,
                .inverse = inverse_workspace.gram,
                .gram = residual_gram,
                .intercept_ridge = residual_intercept_ridge,
                .unregularized = residual_unregularized,
                .workspace = residual_workspace
              };
              R_UnwindProtect(
                np_regression_alllarge_residual_execute, (void *)&residual_call,
                np_regression_alllarge_lp_fit_cleanup, (void *)&all_large_owner,
                NULL);
            }

            if(use_fit_projection_blas){
              for(i = 0; i < num_obs_eval; i += eval_block_rows){
                const int nblock = MIN(eval_block_rows, num_obs_eval - i);

                /*
                 * Keep the fixed leading dimension for every basis block.
                 * The projection helper returns its nblock-by-k result with
                 * leading dimension nblock, including the final short block.
                 */
                for(int row = 0; row < nblock; row++){
                  if(use_bernstein){
                    np_glp_fill_basis_eval(num_reg_continuous,
                                           glp_terms,
                                           glp_nterms,
                                           matrix_X_continuous_eval,
                                           i + row,
                                           basis_ctx,
                                           eval_basis);
                  } else {
                    np_glp_fill_basis_eval_raw(num_reg_continuous,
                                               glp_terms,
                                               glp_nterms,
                                               matrix_X_continuous_eval,
                                               i + row,
                                               eval_basis);
                  }
                  for(j = 0; j < glp_nterms; j++)
                    eval_basis_block[j*eval_block_rows + row] =
                      eval_basis[j];
                }

                np_blas_dgemv_n_int(nblock,
                                    glp_nterms,
                                    eval_basis_block,
                                    eval_block_rows,
                                    beta,
                                    yhat_block);
                if(shortcut_do_merr && !covariance_ridged)
                  np_blas_project_inverse_block_int(
                    nblock,
                    glp_nterms,
                    eval_basis_block,
                    eval_block_rows,
                    inverse_workspace.gram,
                    projection_block);

                for(int row = 0; row < nblock; row++){
                  double q = 0.0;
                  if(shortcut_do_merr && !covariance_ridged)
                    for(j = 0; j < glp_nterms; j++)
                      q += eval_basis_block[j*eval_block_rows + row] *
                        projection_block[j*nblock + row];
                  mean[i + row] = yhat_block[row];
                  if(fit_projection != NULL)
                    np_fit_projection_inverse(fit_projection, inverse_workspace.gram,
                      eval_basis_block + row, eval_block_rows, i + row);
                  if(shortcut_do_merr) {
                    const double mv = sigma2hat*q;
                    mean_stderr[i + row] = covariance_ridged ?
                      np_conditional_alllarge_standard_error(
                        inverse_workspace.matrix_copy, eval_basis_block + row,
                        eval_block_rows, covariance_workspace,
                        inverse_workspace.rhs, covariance_work, glp_nterms,
                        sigma2hat) :
                      (mv > 0.0 && isfinite(mv)) ?
                      sqrt(mv) : se_default;
                  }
                  np_progress_fit_loop_step_owned(
                    i + row + 1,
                    fit_progress_total,
                    np_regression_alllarge_lp_fit_cleanup,
                    (void *)&all_large_owner);
                }
              }
            } else {
              for(i = 0; i < num_obs_eval; i++){
                double yhat = 0.0;
                double q = 0.0;

                if(use_bernstein){
                  np_glp_fill_basis_eval(num_reg_continuous,
                                         glp_terms,
                                         glp_nterms,
                                         matrix_X_continuous_eval,
                                         i,
                                         basis_ctx,
                                         eval_basis);
                } else {
                  np_glp_fill_basis_eval_raw(num_reg_continuous,
                                             glp_terms,
                                             glp_nterms,
                                             matrix_X_continuous_eval,
                                             i,
                                             eval_basis);
                }

                for(j = 0; j < glp_nterms; j++)
                  yhat += eval_basis[j]*beta[j];
                if(shortcut_do_merr && !covariance_ridged)
                  for(j = 0; j < glp_nterms; j++){
                    const double zj = eval_basis[j];
                    for(int b = 0; b < glp_nterms; b++)
                      q += zj *
                        inverse_workspace.matrix_copy[j*glp_nterms + b] *
                        eval_basis[b];
                  }

                mean[i] = yhat;
                if(fit_projection != NULL)
                  np_fit_projection_inverse(fit_projection, inverse_workspace.gram,
                                            eval_basis, 1, i);
                if(shortcut_do_merr) {
                  const double mv = sigma2hat*q;
                  mean_stderr[i] = covariance_ridged ?
                    np_conditional_alllarge_standard_error(
                      inverse_workspace.matrix_copy, eval_basis, 1,
                      covariance_workspace, inverse_workspace.rhs,
                      covariance_work, glp_nterms, sigma2hat) :
                    (mv > 0.0 && isfinite(mv)) ?
                    sqrt(mv) : se_default;
                }

                if(do_grad || first_se_request != NULL){
                  const int nvars = num_reg_continuous + num_reg_unordered + num_reg_ordered;
                  for(l = 0; l < num_reg_continuous; l++){
                    const int grad_order =
                      first_se_request != NULL ? first_se_request->order[l] :
                      ((vector_glp_gradient_order_extern != NULL) ?
                       MAX(1, vector_glp_gradient_order_extern[l]) : 1);
                    const int active = first_se_request != NULL ?
                      first_se_request->se[l] : np_glp_gradient_direction_active(l);
                    double qg = 0.0;
                    double dg = 0.0;

                    if(!active) {
                      if(do_grad) gradient[l][i] = NA_REAL;
                      if(shortcut_do_gerr)
                        gradient_stderr[l][i] = NA_REAL;
                      continue;
                    }

                    if(use_bernstein){
                      np_glp_fill_basis_eval_deriv(l,
                                                   grad_order,
                                                   num_reg_continuous,
                                                   glp_terms,
                                                   glp_nterms,
                                                   matrix_X_continuous_eval,
                                                   i,
                                                   basis_ctx,
                                                   eval_deriv);
                    } else {
                      np_glp_fill_basis_eval_deriv_raw(l,
                                                       grad_order,
                                                       num_reg_continuous,
                                                       glp_terms,
                                                       glp_nterms,
                                                       matrix_X_continuous_eval,
                                                       i,
                                                       eval_deriv);
                    }

                    if(do_grad) {
                      for(j = 0; j < glp_nterms; j++)
                        dg += eval_deriv[j]*beta[j];
                      gradient[l][i] = dg;
                    }

                    if(shortcut_do_gerr){
                      for(j = 0; !covariance_ridged && j < glp_nterms; j++){
                        const double dj = eval_deriv[j];
                        for(int b = 0; b < glp_nterms; b++)
                          qg += dj *
                            inverse_workspace.matrix_copy[j*glp_nterms + b] *
                            eval_deriv[b];
                      }
                      {
                        const double gv = sigma2hat*qg;
                        gradient_stderr[l][i] = covariance_ridged ?
                          np_conditional_alllarge_standard_error(
                            inverse_workspace.matrix_copy, eval_deriv, 1,
                            covariance_workspace, inverse_workspace.rhs,
                            covariance_work, glp_nterms, sigma2hat) :
                          (gv > 0.0 && isfinite(gv)) ?
                          sqrt(gv) : 0.0;
                      }
                    }
                  }

                  for(l = num_reg_continuous; do_grad && l < nvars; l++){
                    /* The all-large gate's accepted point map is invariant
                     * to every categorical endpoint. */
                    gradient[l][i] = 0.0;
                    if(shortcut_do_gerr) gradient_stderr[l][i] = NA_REAL;
                  }
                }
                np_progress_fit_loop_step_owned(
                  i + 1,
                  fit_progress_total,
                  np_regression_alllarge_lp_fit_cleanup,
                  (void *)&all_large_owner);
              }
            }
            estimation_shortcut_done = 1;
            hc0_point_shortcut_done = point_only_hc0_shortcut;
          }
        }

        if(use_bernstein && (basis_ctx != NULL)){
          for(l = 0; l < num_reg_continuous; l++) np_glp_basis_ctx_free(&basis_ctx[l]);
          free(basis_ctx);
        }
        free(eval_basis_block);
        free(projection_block);
        free(yhat_block);
        free(covariance_workspace);
        free(residual_gram);
        np_alllarge_residual_workspace_clear(residual_workspace);
        free(residual_workspace);
        if(eval_basis != NULL) free(eval_basis);
        if(eval_deriv != NULL) free(eval_deriv);
        free(beta);
        np_lp_full_row_workspace_clear(&inverse_workspace);
        if(basis != NULL){
          if(basis_is_contiguous)
            free_tmat(basis);
          else
            free_mat(basis, glp_nterms);
        }
        if(glp_terms != NULL) free(glp_terms);
      }
    }
  }

  if(estimation_shortcut_done) {
    if(!hc0_point_shortcut_done)
      goto finish_regression_estimation;
    hc0_point_context = *hc0_context;
    hc0_point_context.point_already_computed = 1;
    effective_hc0_context = &hc0_point_context;
  }

  // Conduct the estimation 

  if(lp_engine_est == NP_LP_ENGINE_SCALAR) { // canonical scalar LP0
    // Nadaraya-Watson
    // Generate bandwidth vector given scale factors, nearest neighbors, or lambda 
    const NPRegressionScalarFitCall scalar_call = {
      .kernel_c = kernel_c,
      .kernel_u = kernel_u,
      .kernel_o = kernel_o,
      .lp_engine = lp_engine_est,
      .bandwidth_mode = BANDWIDTH_reg,
      .num_obs_train = num_obs_train,
      .num_obs_eval = num_obs_eval,
      .num_obs_eval_alloc = num_obs_eval_alloc,
      .num_reg_unordered = num_reg_unordered,
      .num_reg_ordered = num_reg_ordered,
      .num_reg_continuous = num_reg_continuous,
      .matrix_X_unordered_train = matrix_X_unordered_train,
      .matrix_X_ordered_train = matrix_X_ordered_train,
      .matrix_X_continuous_train = matrix_X_continuous_train,
      .matrix_X_unordered_eval = matrix_X_unordered_eval,
      .matrix_X_ordered_eval = matrix_X_ordered_eval,
      .matrix_X_continuous_eval = matrix_X_continuous_eval,
      .vector_Y = vector_Y,
      .vector_scale_factor = vector_scale_factor,
      .lambda = lambda,
      .num_categories = num_categories,
      .operator = operator,
      .matrix_categorical_vals = matrix_categorical_vals,
      .matrix_bandwidth = matrix_bandwidth,
      .matrix_ordered_indices = matrix_ordered_indices,
      .mean = mean,
      .gradient = gradient,
      .mean_stderr = mean_stderr,
      .gradient_stderr = gradient_stderr,
      .tree_enabled = int_TREE_X,
      .do_merr = do_merr,
      .do_grad = do_grad,
      .do_gerr = do_gerr,
      .standard_error_mode = standard_error_mode,
      .ordinary_response = ordinary_response,
      .hc0_context = hc0_context,
      .gate_context = &gate_ctx_local,
      .enclosing_owner = &fit_owner,
      .kernel_squared_integral = K_INT_KERNEL_P,
      .bandwidth_product = hprod,
      .gradient_stderr_factor = gfac,
      .empty_rows = failure != NULL ? empty_rows : NULL
    };
    scalar_fit_status = np_regression_scalar_fit(&scalar_call);

    if(scalar_fit_status == NP_REGRESSION_SCALAR_FIT_PROFILE)
      goto finish_regression_estimation;
    if(scalar_fit_status != NP_REGRESSION_SCALAR_FIT_OK)
      goto finish_regression_estimation;
  } else if(lp_engine_est == NP_LP_ENGINE_GENERAL) { // local polynomial (regtype = "lp")
    /* Reuse the NN search certificate for ordinary and conditional GNN fitting.
       Full support proves traversal cannot discard a donor; dense moments
       still evaluate the varying kernel weights. Fixed owners are untouched. */
    const int fit_tree_enabled =
      BANDWIDTH_reg == BW_GEN_NN &&
      int_TREE_X == NP_TREE_TRUE && kernel_route == NULL &&
      !int_cker_bound_extern &&
      np_reg_gnn_tree_cannot_prune(num_obs_eval, num_reg_continuous,
                                   kdt_extern_X, kernel_c, operator,
                                   matrix_X_continuous_eval, matrix_bandwidth) ?
        NP_TREE_FALSE : int_TREE_X;
    const NPRegressionGeneralLPFitCall general_lp_call = {
      .prepared_ann_row = prepared_ann_row,
      .fit_projection = fit_projection,
      .kernel_c = kernel_c,
      .kernel_u = kernel_u,
      .kernel_o = kernel_o,
      .bandwidth_mode = BANDWIDTH_reg,
      .num_obs_train = num_obs_train,
      .num_obs_eval = num_obs_eval,
      .num_reg_unordered = num_reg_unordered,
      .num_reg_ordered = num_reg_ordered,
      .num_reg_continuous = num_reg_continuous,
      .matrix_X_unordered_train = matrix_X_unordered_train,
      .matrix_X_ordered_train = matrix_X_ordered_train,
      .matrix_X_continuous_train = matrix_X_continuous_train,
      .matrix_X_unordered_eval = matrix_X_unordered_eval,
      .matrix_X_ordered_eval = matrix_X_ordered_eval,
      .matrix_X_continuous_eval = matrix_X_continuous_eval,
      .vector_Y = vector_Y,
      .vector_scale_factor = vector_scale_factor,
      .lambda = lambda,
      .num_categories = num_categories,
      .operator = operator,
      .matrix_categorical_vals = matrix_categorical_vals,
      .matrix_bandwidth = matrix_bandwidth,
      .categorical_matrix_bandwidth = categorical_matrix_bandwidth,
      .categorical_base_requires_refit =
        categorical_base_requires_refit,
      .categorical_point_invariant = categorical_point_invariant,
      .all_large_kernel_rows = all_large_kernel_rows,
      .ordinary_response = ordinary_response,
      .mean = mean,
      .gradient = gradient,
      .mean_stderr = mean_stderr,
      .gradient_stderr = gradient_stderr,
      .tree_enabled = fit_tree_enabled,
      .do_merr = do_merr,
      .do_grad = do_grad,
      .do_gerr = do_gerr,
      .categorical_compress = categorical_compress,
      .fit_progress_total = fit_progress_total,
#ifdef MPI2
      .fit_progress_active = fit_progress_active,
#endif
      .gate_context = &gate_ctx_local,
      .kernel_squared_integral = K_INT_KERNEL_P,
      .bandwidth_product = hprod,
      .kernel_route = kernel_route,
      .kernel_route_diagnostics = kernel_route_diagnostics,
      .hc0_context = effective_hc0_context,
      .enclosing_owner = &fit_owner,
      .empty_rows = empty_rows,
      .first_se_request = first_se_request,
      .conditional_variance = conditional_variance,
      .prepared_design = prepared_design,
      .continuous_range = np_regression_gradient_range(
        hc0_context != NULL ? hc0_context->gradient_request : NULL,
        0, num_reg_continuous),
      .unordered_range = np_regression_gradient_range(
        hc0_context != NULL ? hc0_context->gradient_request : NULL,
        num_reg_continuous, num_reg_continuous + num_reg_unordered),
      .ordered_range = np_regression_gradient_range(
        hc0_context != NULL ? hc0_context->gradient_request : NULL,
        num_reg_continuous + num_reg_unordered,
        num_reg_continuous + num_reg_unordered + num_reg_ordered),
      .categorical_range = np_regression_gradient_range(
        hc0_context != NULL ? hc0_context->gradient_request : NULL,
        num_reg_continuous,
        num_reg_continuous + num_reg_unordered + num_reg_ordered)
    };

    general_lp_fit_status = np_regression_general_lp_fit(&general_lp_call);
  }

finish_regression_estimation:
  np_regression_fit_owner_clear(&fit_owner);

  if(regression_fit_status == NP_REGRESSION_FIT_ERR_ALLOC)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_FIT, regression_fit_status, "\n** Error: memory allocation failed.");
  if(regression_fit_status == NP_REGRESSION_FIT_ERR_BANDWIDTH)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_FIT, regression_fit_status, "\n** Error: invalid bandwidth.");
  if(regression_fit_status == NP_REGRESSION_FIT_ERR_ZERO_NN_RADIUS) {
    const NPNNGeometryContext external_geometry = {.mode = NP_NN_QUERY_EXTERNAL};
    const NPNNZeroRadiusInfo info = np_nn_zero_radius_info(
      BANDWIDTH_reg, num_obs_train, num_obs_eval, num_reg_continuous,
      matrix_X_continuous_train, matrix_X_continuous_eval, vector_scale_factor,
      kernel_route != NULL ? NULL :
        (zero_radius_external_query ? &external_geometry : nn_geometry_context), 0);
    if(failure != NULL) {
      np_regression_failure_zero_radius(failure, &info);
      return NP_REGRESSION_FIT_ERR_ZERO_NN_RADIUS;
    }
    np_nn_zero_radius_error(&info);
  }
  if(regression_fit_status == NP_REGRESSION_FIT_ERR_HASH_CREATE)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_FIT, regression_fit_status, "hash table creation failed");
  if(regression_fit_status == NP_REGRESSION_FIT_ERR_HASH_INSERT)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_FIT, regression_fit_status, "insertion into hash table failed");
  if(regression_fit_status == NP_REGRESSION_FIT_ERR_HASH_LOOKUP)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_FIT, regression_fit_status, "hash table lookup failed (which should be impossible)");
  if(regression_fit_status != NP_REGRESSION_FIT_OK)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_FIT, regression_fit_status, "invalid internal regression fit status");

  if(scalar_fit_status == NP_REGRESSION_SCALAR_FIT_ERR_DIMENSION)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_SCALAR, scalar_fit_status, "conditional influence kernel-weight dimensions overflow");
  if(scalar_fit_status == NP_REGRESSION_SCALAR_FIT_ERR_CONDITIONAL_ALLOC)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_SCALAR, scalar_fit_status, "memory allocation failed in conditional influence workspace");
  if(scalar_fit_status == NP_REGRESSION_SCALAR_FIT_ERR_ALLOC)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_SCALAR, scalar_fit_status, "\n** Error: memory allocation failed.");
  if(scalar_fit_status == NP_REGRESSION_SCALAR_FIT_ERR_TRAVERSAL)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_SCALAR, scalar_fit_status, "conditional regression kernel traversal failed");
  if(scalar_fit_status == NP_REGRESSION_SCALAR_FIT_ERR_INFLUENCE)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_SCALAR, scalar_fit_status, "conditional influence variance construction failed");
  if(scalar_fit_status == NP_REGRESSION_SCALAR_FIT_ERR_WORKSPACE_DIMENSION)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_SCALAR, scalar_fit_status, "scalar regression workspace dimensions overflow");
  if(scalar_fit_status == NP_REGRESSION_SCALAR_FIT_ERR_HC0)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_SCALAR, scalar_fit_status, "ordinary-regression HC0 variance construction failed");
  if(scalar_fit_status != NP_REGRESSION_SCALAR_FIT_OK &&
     scalar_fit_status != NP_REGRESSION_SCALAR_FIT_PROFILE)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_SCALAR, scalar_fit_status, "invalid internal scalar regression fit status");

  if(general_lp_fit_status == NP_REGRESSION_GENERAL_LP_FIT_ERR_ROUTE) {
    if(kernel_route_diagnostics != NULL &&
       kernel_route_diagnostics->beta_status != NP_BETA_OK)
      NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "canonical beta LP row failed in continuous dimension %d: %s",
            kernel_route_diagnostics->bad_coordinate + 1,
            np_beta_status_message(kernel_route_diagnostics->beta_status));
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "canonical beta LP row failed");
  }
  if(general_lp_fit_status == NP_REGRESSION_GENERAL_LP_FIT_ERR_DEGREE)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "glp degree vector unavailable");
  if(general_lp_fit_status == NP_REGRESSION_GENERAL_LP_FIT_ERR_TERMS)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "failed to build glp basis terms");
  if(general_lp_fit_status == NP_REGRESSION_GENERAL_LP_FIT_ERR_DIMENSION)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "invalid glp basis dimension");
  if(general_lp_fit_status == NP_REGRESSION_GENERAL_LP_FIT_ERR_SOLVE_ALLOC)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "memory allocation failed in glp solve workspace");
  if(general_lp_fit_status == NP_REGRESSION_GENERAL_LP_FIT_ERR_ALLOC)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "memory allocation failed in glp path");
  if(general_lp_fit_status == NP_REGRESSION_GENERAL_LP_FIT_ERR_BASIS)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "failed to initialize glp Bernstein basis");
  if(general_lp_fit_status == NP_REGRESSION_GENERAL_LP_FIT_ERR_SOLVE)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "LP solve failed in glp path");
  if(general_lp_fit_status == NP_REGRESSION_GENERAL_LP_FIT_ERR_OWNER_SOLVE)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "LP solve failed in MPI owner-row path");
  if(general_lp_fit_status == NP_REGRESSION_GENERAL_LP_FIT_ERR_HC0)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "ordinary-regression HC0 variance construction failed in glp path");
  if(general_lp_fit_status != NP_REGRESSION_GENERAL_LP_FIT_OK)
    NP_REGRESSION_RETURN_FAILURE(NP_REGRESSION_FAILURE_GENERAL, general_lp_fit_status, "invalid internal general LP fit status");

  np_regression_fit_statistics(
    num_obs_eval, vector_Y_eval, mean,
    R_squared, MSE, MAE, MAPE, CORR, SIGN);

  return(0);
}
