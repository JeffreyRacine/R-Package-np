static int NP_OUTER_PACK_ADJACENT_HOT_ALIGN kernel_weighted_sum_np_ctx_ex(
int * KERNEL_reg,
int * KERNEL_unordered_reg,
int * KERNEL_ordered_reg,
const int BANDWIDTH_reg,
const int num_obs_train,
const int num_obs_eval,
const int num_reg_unordered,
const int num_reg_ordered,
const int num_reg_continuous,
const int leave_one_out,
const int leave_one_out_offset,
const int kernel_pow,
const int bandwidth_divide,
const int bandwidth_divide_weights,
const int symmetric,
const int gather_scatter,
const int drop_one_train,
const int drop_which_train,
const int * const operator,
const int permutation_operator,
int do_score,
int do_ocg,
int * bpso,
const int suppress_parallel,
const int ncol_Y,
const int ncol_W,
const int int_TREE,
const int do_partial_tree,
KDT * const kdt,
NL * const inl,
int * const nld,
int * const idx,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double **matrix_Y,
double **matrix_W,
double * sgn,
double *vector_scale_factor,
int bandwidth_provided,
double ** matrix_bw_train,
double ** matrix_bw_eval,
double * lambda_pre,
int *num_categories,
double **matrix_categorical_vals,
int ** matrix_ordered_indices,
double * const weighted_sum,
double * const weighted_permutation_sum,
double * const kw,
const NP_GateOverrideCtx * const gate_override_ctx,
const NP_DualPowerCtx * const dual_power_ctx,
const NP_OuterPackCtx * const outer_pack_ctx,
const NPContinuousKernelExecutionContext * const kernel_execution_context,
const NPCenteredMomentCtx * const centered_moment_ctx,
NPPermutationWeightOutput * const pkw_output,
  NPLPSupportRankContext * const support_rank){
  if(support_rank != NULL) {
    if(num_obs_eval != 1 || !suppress_parallel || support_rank->cap <= 0 ||
       leave_one_out || drop_one_train || symmetric || gather_scatter ||
       do_score || do_ocg || permutation_operator != OP_NOOP ||
       kernel_pow != 1)
      return KWSNP_ERR_BADINVOC;
    support_rank->count = 0;
    np_lp_design_support_begin(support_rank->design);
  }
  const NP_GateOverrideCtx * const gate_ctx_raw =
    (gate_override_ctx != NULL) ? gate_override_ctx : &np_gate_override_ctx;
  const NP_GateOverrideCtx gate_ctx_empty = {0};
  const NP_GateOverrideCtx * const gate_ctx =
    np_gate_ctx_is_sane(gate_ctx_raw) ? gate_ctx_raw : &gate_ctx_empty;
  const int caller_override_active =
    ((gate_ctx != NULL) && (gate_ctx->active == NP_GATE_CTX_OVERRIDE));
  const int disable_gate_features =
    (((!np_partial_gate_features_enabled) && (!caller_override_active))) ||
    ((gate_ctx != NULL) && (gate_ctx->active == NP_GATE_CTX_DISABLE));
  assert(np_gate_ctx_is_sane(gate_ctx));

  const NPNNKernelFold * const fold = kernel_execution_context == NULL ? NULL :
    kernel_execution_context->fold;
  if(fold != NULL) {
    if(fold->num_train != num_obs_train || fold->num_eval != num_obs_eval ||
       fold->num_continuous != num_reg_continuous || num_reg_continuous <= 0 ||
       fold->eval_to_train == NULL || leave_one_out || drop_one_train ||
       symmetric || gather_scatter || do_score || do_ocg ||
       permutation_operator != OP_NOOP || dual_power_ctx != NULL ||
       centered_moment_ctx != NULL || support_rank != NULL ||
       (kernel_pow != 1 && kernel_pow != 2))
      return KWSNP_ERR_BADINVOC;
    for(int d = 0; d < num_reg_continuous + num_reg_unordered + num_reg_ordered; ++d)
      if(operator[d] != OP_NORMAL) return KWSNP_ERR_BADINVOC;
  }

  if(kernel_execution_context != NULL && kernel_execution_context->route != NULL) {
    const NPContinuousKernelRoute * const kernel_route =
      kernel_execution_context->route;
    NPContinuousKernelDerivativeDiagnostics * const
      kernel_route_diagnostics = kernel_execution_context->diagnostics;
    int route_has_convolution = 0;
    int route_has_derivative = 0;
    int route_coordinate;
    const int beta_scaled_outer =
      dual_power_ctx != NULL && dual_power_ctx->retain_common_scale;
    const int beta_dual_power =
      dual_power_ctx != NULL && dual_power_ctx->weighted_sum != NULL &&
      dual_power_ctx->kernel_pow == 2;
    /* A fit can retain progress while declining the optional second moment. */
    const int beta_point_request =
      dual_power_ctx != NULL && dual_power_ctx->weighted_sum == NULL &&
      dual_power_ctx->kernel_pow == 2 &&
      dual_power_ctx->matrix_Y == NULL && dual_power_ctx->matrix_W == NULL &&
      dual_power_ctx->ncol_Y == 0 && dual_power_ctx->ncol_W == 0 &&
      !dual_power_ctx->retain_common_scale &&
      dual_power_ctx->observation_scale == NULL &&
      dual_power_ctx->regression_derivative == NULL;
    const int beta_retain_common_scale =
      beta_scaled_outer;
    const int beta_centered_moment =
      centered_moment_ctx != NULL &&
      centered_moment_ctx->centered_m2 != NULL;
    double ** const beta_power2_Y = beta_dual_power ?
      ((dual_power_ctx->matrix_Y != NULL) ?
       dual_power_ctx->matrix_Y : matrix_Y) : NULL;
    double ** const beta_power2_W = beta_dual_power ?
      ((dual_power_ctx->matrix_W != NULL) ?
       dual_power_ctx->matrix_W : matrix_W) : NULL;
    const int beta_power2_ncol_Y = beta_dual_power ?
      ((dual_power_ctx->matrix_Y != NULL) ?
       dual_power_ctx->ncol_Y : ncol_Y) : 0;
    const int beta_power2_ncol_W = beta_dual_power ?
      ((dual_power_ctx->matrix_W != NULL) ?
       dual_power_ctx->ncol_W : ncol_W) : 0;
    const int beta_has_categories =
      num_reg_unordered > 0 || num_reg_ordered > 0;
    const int beta_categories_valid =
      num_reg_unordered >= 0 && num_reg_ordered >= 0 &&
      num_reg_unordered <= INT_MAX - num_reg_ordered &&
      (!beta_has_categories ||
       (((num_reg_unordered == 0) ||
         (KERNEL_unordered_reg != NULL &&
          matrix_X_unordered_train != NULL &&
          matrix_X_unordered_eval != NULL)) &&
        ((num_reg_ordered == 0) ||
         (KERNEL_ordered_reg != NULL &&
          matrix_X_ordered_train != NULL &&
          matrix_X_ordered_eval != NULL)) &&
        lambda_pre != NULL && num_categories != NULL &&
        matrix_categorical_vals != NULL &&
        (kernel_execution_context->categorical_compress == 0 ||
         kernel_execution_context->categorical_compress == 1)));
    for(route_coordinate = 0; route_coordinate < num_reg_continuous;
        ++route_coordinate) {
      route_has_convolution |=
        operator != NULL && operator[route_coordinate] == OP_CONVOLUTION;
      route_has_derivative |=
        operator != NULL && operator[route_coordinate] == OP_DERIVATIVE;
    }
    const int exact_beta_absolute_route =
      np_continuous_kernel_route_has_beta(kernel_route) &&
      (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN ||
       BANDWIDTH_reg == BW_ADAP_NN) &&
      bandwidth_divide == 0 && bandwidth_divide_weights == 0 &&
      beta_categories_valid &&
      num_reg_continuous > 0 && permutation_operator == OP_NOOP &&
      !do_score && !do_ocg &&
      int_TREE == NP_TREE_FALSE && !do_partial_tree &&
      !symmetric && !gather_scatter && !drop_one_train &&
      sgn == NULL &&
      ((BANDWIDTH_reg == BW_FIXED && !bandwidth_provided &&
        matrix_bw_train == NULL && matrix_bw_eval == NULL &&
        vector_scale_factor != NULL) ||
       (BANDWIDTH_reg == BW_GEN_NN && bandwidth_provided &&
        matrix_bw_eval != NULL &&
        (!route_has_convolution || matrix_bw_train != NULL)) ||
       (BANDWIDTH_reg == BW_ADAP_NN && bandwidth_provided &&
        matrix_bw_train != NULL &&
        (!route_has_convolution || matrix_bw_eval != NULL))) &&
      ((!beta_has_categories && lambda_pre == NULL) ||
       beta_has_categories) &&
      (dual_power_ctx == NULL ||
       ((beta_dual_power || beta_scaled_outer || beta_point_request) &&
        kernel_pow == 1 && !route_has_derivative &&
        (kw == NULL || beta_retain_common_scale))) &&
      (centered_moment_ctx == NULL ||
       (dual_power_ctx == NULL &&
        kernel_pow == 1 && !route_has_derivative && kw == NULL &&
        ncol_Y == 0 && ncol_W == 0)) &&
      outer_pack_ctx == NULL &&
      weighted_sum != NULL;
    double *route_row = NULL;
    int route_status;

    if(!exact_beta_absolute_route)
      return KWSNP_ERR_BADINVOC;
    if((!route_has_derivative && !beta_dual_power &&
        !beta_centered_moment) || beta_has_categories ||
       beta_retain_common_scale) {
      if((size_t)num_obs_train > SIZE_MAX / sizeof(double))
        return KWSNP_ERR_BADINVOC;
      route_row = (double *)malloc((size_t)num_obs_train * sizeof(double));
      if(route_row == NULL)
        return KWSNP_ERR_BADINVOC;
    }
    {
      const NPBetaAbsoluteRouteCall route_call = {
        .route = kernel_route,
        .bandwidth_mode = BANDWIDTH_reg,
        .num_obs_train = num_obs_train,
        .num_obs_eval = num_obs_eval,
        .num_reg_unordered = num_reg_unordered,
        .num_reg_ordered = num_reg_ordered,
        .num_reg_continuous = num_reg_continuous,
        .leave_one_out = leave_one_out,
        .leave_one_out_offset = leave_one_out_offset,
        .operator = operator,
        .matrix_X_continuous_train = matrix_X_continuous_train,
        .matrix_X_continuous_eval = matrix_X_continuous_eval,
        .matrix_X_unordered_train = matrix_X_unordered_train,
        .matrix_X_ordered_train = matrix_X_ordered_train,
        .matrix_X_unordered_eval = matrix_X_unordered_eval,
        .matrix_X_ordered_eval = matrix_X_ordered_eval,
        .kernel_unordered = KERNEL_unordered_reg,
        .kernel_ordered = KERNEL_ordered_reg,
        .categorical_lambda = lambda_pre,
        .num_categories = num_categories,
        .matrix_categorical_vals = matrix_categorical_vals,
        .categorical_compress =
          kernel_execution_context->categorical_compress,
        .vector_bandwidth = vector_scale_factor,
        .matrix_bandwidth_eval = matrix_bw_eval,
        .matrix_bandwidth_train = matrix_bw_train,
        .matrix_Y = matrix_Y,
        .matrix_W = matrix_W,
        .ncol_Y = ncol_Y,
        .ncol_W = ncol_W,
        .matrix_Y_power2 = beta_power2_Y,
        .matrix_W_power2 = beta_power2_W,
        .ncol_Y_power2 = beta_power2_ncol_Y,
        .ncol_W_power2 = beta_power2_ncol_W,
        .kernel_power = kernel_pow,
        .row = route_row,
        .weighted_sum = weighted_sum,
        .weighted_sum_power2 =
          beta_dual_power ? dual_power_ctx->weighted_sum : NULL,
        .power2_observation_scale = beta_dual_power ?
          dual_power_ctx->observation_scale : NULL,
        .retain_common_scale = beta_retain_common_scale,
        .support_rank = support_rank != NULL ? &support_rank->count : NULL,
        .support_cap = support_rank != NULL ? support_rank->cap : 0,
        .support_design = support_rank != NULL ? support_rank->design : NULL,
        .centered_m2 = beta_centered_moment ?
          centered_moment_ctx->centered_m2 : NULL,
        .kw = kw,
        .regression_moment_context = NULL,
        .fold = kernel_execution_context->fold,
        .route_diagnostics = kernel_route_diagnostics,
        .progress = dual_power_ctx != NULL ? dual_power_ctx->progress :
          (centered_moment_ctx != NULL ? centered_moment_ctx->progress : NULL)
      };

      route_status = np_beta_absolute_route(&route_call);
    }
    return route_status;
  }
  
  /* This function takes a vector Y and returns a kernel weighted
     leave-one-out sum. By default Y should be a vector of ones
     (simply compute the kernel sum). This function will allow users
     to `roll their own' with mixed data leave-one out kernel sums. */

  /*
   * Every categorical kernel consumes its category-count metadata.  Reject
   * an incomplete internal context before entering the hot loop so an owner
   * allocation defect cannot become a null dereference.
   */
  if(((num_reg_unordered > 0) || (num_reg_ordered > 0)) &&
     (num_categories == NULL))
    return 1;

  /* Declarations */

  int i, ii, j, kk, k, l, mstep, js, je, num_obs_eval_alloc, sum_element_length, ip;
  int status = 0;
  int do_psum, swap_xxt;
  NPNNGeometryStatus nn_geometry_status = NP_NN_GEOMETRY_OK;

  int * permutation_kernel = NULL;
  int doscoreocg = do_score || do_ocg;
  int do_perm = permutation_operator != OP_NOOP; 

  
  const int no_bpso = (NULL == bpso);

  int p_nvar;

  if(!np_runtime_tol_cache_ready)
    np_refresh_runtime_tolerances();
  if((outer_pack_ctx == NULL) || !outer_pack_ctx->runtime_options_frozen)
    np_refresh_mseries_accelerate_option();

  if(no_bpso){
    bpso = (int *)malloc((num_reg_unordered + num_reg_ordered + num_reg_continuous)*sizeof(int));

    for(i = 0; i < num_reg_continuous; i++)
      bpso[i] = do_perm;

    for(i = num_reg_continuous; i < num_reg_continuous+num_reg_unordered+num_reg_ordered; i++)
      bpso[i] = doscoreocg;

    p_nvar = (do_perm ? num_reg_continuous : 0) +
      (doscoreocg ? num_reg_unordered + num_reg_ordered : 0);

  } else {
    // it's possible that do_perm and do_ocg could be set, but no continuous or categorical
    // variables are enabled via bpso, e.g. conditional density gradients

    for(i = 0, ii = 0, l = 0; i < num_reg_continuous; i++){
      ii |= bpso[i];
      l += bpso[i];
    }

    do_perm &= ii;

    for(i = num_reg_continuous, ii = 0; i < num_reg_continuous + num_reg_unordered + num_reg_ordered; i++){
      ii |= bpso[i];
      l += bpso[i];
    }

    do_ocg &= ii;

    doscoreocg = do_score || do_ocg;
    p_nvar = l;
  }

  int * ps_ukernel = NULL, * ps_okernel = NULL;

  int ps_ok_nli = (num_reg_ordered != 0) && (KERNEL_ordered_reg[0] != 1);

  /* Trees are currently not compatible with all operations */
  int np_ks_tree_use = (int_TREE == NP_TREE_TRUE);
  int np_ks_tree_requested;
  int fixed_tree_common_geometry = 0;
  NPTreeCapability fixed_tree_capability = NP_TREE_CAPABILITY_UNKNOWN;
  int any_convolution = 0;
  int all_continuous_normal = 1;
  int is_adaptive = (BANDWIDTH_reg == BW_ADAP_NN);
  const NPNNKernelFold * const adaptive_fold = is_adaptive ? fold : NULL;
  XL fold_primary = {0}, fold_successor = {0};
  if(adaptive_fold != NULL) {
    if(!bandwidth_provided || adaptive_fold->primary == NULL ||
       adaptive_fold->successor == NULL || adaptive_fold->scale == NULL ||
       (kw != NULL && bandwidth_divide != bandwidth_divide_weights))
      return KWSNP_ERR_BADINVOC;
    fold_primary.istart = (int *)R_alloc((size_t)num_obs_eval, sizeof(int));
    fold_primary.nlev = (int *)R_alloc((size_t)num_obs_eval, sizeof(int));
    fold_successor.istart = (int *)R_alloc((size_t)num_obs_eval, sizeof(int));
    fold_successor.nlev = (int *)R_alloc((size_t)num_obs_eval, sizeof(int));
  }
  int tree_pkw_sparse = 0;

  int lod = 0;

  const int nws = (weighted_sum == NULL);
  const int do_dual_power =
    (dual_power_ctx != NULL) && (dual_power_ctx->weighted_sum != NULL);
  const NPRegressionHC0DerivativeMomentCtx * const hc0_derivative_ctx =
    dual_power_ctx != NULL ? dual_power_ctx->regression_derivative : NULL;
  const int do_hc0_derivative = hc0_derivative_ctx != NULL;
  NPRegressionResidualKernelContext * const residual_preparation =
    dual_power_ctx != NULL ? dual_power_ctx->residual_preparation : NULL;
  const NPRegressionHC0Context * const uncertainty_context =
    dual_power_ctx != NULL ? dual_power_ctx->uncertainty_context : NULL;
  NPRegressionScalarANNInformation * const ann_information =
    dual_power_ctx != NULL ? dual_power_ctx->ann_information : NULL;
  const NPRegressionScalarContrastContext * const scalar_contrast =
    dual_power_ctx != NULL ? dual_power_ctx->scalar_contrast : NULL;
  int information_invalid = 0;
  if(scalar_contrast != NULL &&
     (((BANDWIDTH_reg == BW_ADAP_NN) != (scalar_contrast->adaptive != NULL)) ||
      permutation_operator != OP_DERIVATIVE ||
      scalar_contrast->num_eval != num_obs_eval ||
      scalar_contrast->num_eval_alloc < num_obs_eval ||
      scalar_contrast->response == NULL || scalar_contrast->gradient == NULL ||
      (scalar_contrast->want_variance &&
       (scalar_contrast->uncertainty == NULL || scalar_contrast->gradient_stderr == NULL))))
    return KWSNP_ERR_BADINVOC;
  if(ann_information != NULL &&
     (BANDWIDTH_reg != BW_ADAP_NN ||
      permutation_operator != OP_DERIVATIVE ||
      ann_information->derivative == NULL ||
      ann_information->derivative->num_eval != num_obs_eval ||
      ann_information->derivative->num_predictors != p_nvar ||
      ann_information->derivative->uncertainty_context == NULL ||
      ann_information->derivative->uncertainty_context->gradient_unavailable == NULL ||
      ann_information->level_denominator == NULL ||
      ann_information->alternate_denominator == NULL ||
      ann_information->certificate == NULL || do_dual_power || do_hc0_derivative))
    return KWSNP_ERR_BADINVOC;
  if(residual_preparation != NULL) {
    const int profile = residual_preparation->profile_eval_to_train != NULL;
    if(residual_preparation->output == NULL ||
       residual_preparation->state == NULL ||
       residual_preparation->num_rows != num_obs_eval ||
       leave_one_out || drop_one_train || gather_scatter || symmetric ||
       do_score || (permutation_operator != OP_NOOP &&
                    permutation_operator != OP_DERIVATIVE) ||
       kernel_pow != 1 || bandwidth_divide != 1 || ncol_W != 0 ||
       ncol_Y != 2 || weighted_sum == NULL || do_dual_power ||
       (!profile &&
        (residual_preparation->output->num_obs_train != num_obs_train ||
         residual_preparation->output->num_obs_eval != num_obs_eval)) ||
       (profile && (num_reg_continuous != 0 ||
                    residual_preparation->profile_count == NULL ||
                    residual_preparation->profile_anchor == NULL ||
                    residual_preparation->profile_centered_sum == NULL ||
                    residual_preparation->profile_self_weight == NULL)))
      return KWSNP_ERR_BADINVOC;
    for(i = 0; i < num_reg_continuous + num_reg_unordered + num_reg_ordered; ++i)
      if(operator[i] != OP_NORMAL) return KWSNP_ERR_BADINVOC;
  }
  NPCategoricalLeadingMomentCtx * const categorical_leading =
    outer_pack_ctx != NULL ? outer_pack_ctx->categorical_leading_moments : NULL;
  if(categorical_leading != NULL &&
     (BANDWIDTH_reg == BW_ADAP_NN || num_reg_continuous <= 0 ||
      num_reg_unordered + num_reg_ordered <= 0 || leave_one_out ||
      drop_one_train || gather_scatter || symmetric || do_score ||
      (!categorical_leading->conditional_base_row &&
       (do_ocg || permutation_operator != OP_NOOP)) ||
      (categorical_leading->conditional_base_row &&
       permutation_operator != OP_NOOP && permutation_operator != OP_DERIVATIVE) ||
      kernel_pow != 1 || ncol_W != 0 ||
      ncol_Y != 0 || weighted_sum == NULL || categorical_leading->sum == NULL))
    return KWSNP_ERR_BADINVOC;
  NPConditionalLeadingRatioCtx * const conditional_ratio =
    outer_pack_ctx != NULL ? outer_pack_ctx->conditional_leading_ratio : NULL;
  NPConditionalANNDirectCtx * const conditional_ann =
    outer_pack_ctx != NULL ? outer_pack_ctx->conditional_ann_direct : NULL;
  if(conditional_ann != NULL) {
    if(BANDWIDTH_reg != BW_ADAP_NN ||
       num_reg_continuous+num_reg_unordered+num_reg_ordered <= 0 ||
       leave_one_out || drop_one_train || gather_scatter || symmetric ||
       do_score || kernel_pow != 1 || bandwidth_divide != 1 ||
       ncol_W != 0 || ncol_Y != 0 || weighted_sum == NULL || kw != NULL ||
       do_dual_power || categorical_leading != NULL || conditional_ratio != NULL ||
       outer_pack_ctx->categorical_density_moments != NULL ||
       outer_pack_ctx->row_tile_sink != NULL ||
       conditional_ann->n != num_obs_train || conditional_ann->m != num_obs_eval ||
       conditional_ann->response == NULL || conditional_ann->level == NULL ||
       conditional_ann->denominator == NULL || conditional_ann->mean == NULL ||
       conditional_ann->response->yncon+conditional_ann->response->ynuno+
         conditional_ann->response->ynord <= 0 ||
       (conditional_ann->response->yncon > 0 &&
        (conditional_ann->response->ykc == NULL ||
         conditional_ann->response->train_c == NULL ||
         conditional_ann->response->eval_c == NULL ||
         conditional_ann->response->bandwidth == NULL)) ||
       conditional_ann->response_value == NULL ||
       conditional_ann->response_donor == NULL ||
       (permutation_operator != OP_NOOP && permutation_operator != OP_DERIVATIVE) ||
       (conditional_ann->directions != 0 &&
        (conditional_ann->directions != num_reg_continuous ||
         conditional_ann->derivative == NULL || !do_perm ||
         conditional_ann->derivative_denominator == NULL ||
         conditional_ann->gradient == NULL ||
         p_nvar < num_reg_continuous || weighted_permutation_sum == NULL)))
      return KWSNP_ERR_BADINVOC;
    for(i=0;i<num_reg_continuous;++i)
      if(operator[i] != OP_NORMAL ||
         (conditional_ann->directions > 0 && !bpso[i]))
        return KWSNP_ERR_BADINVOC;
#ifdef NP_ANN_DIRECT_DIAGNOSTIC
    ++conditional_ann->entries;
#endif
  }
  if(conditional_ratio != NULL &&
     (BANDWIDTH_reg == BW_ADAP_NN || leave_one_out || drop_one_train ||
      gather_scatter || symmetric || do_score || kernel_pow != 1 ||
      ncol_W != 0 || ncol_Y != 0 || weighted_sum == NULL ||
      conditional_ratio->sum == NULL || conditional_ratio->numerator == NULL ||
      (permutation_operator != OP_NOOP && permutation_operator != OP_DERIVATIVE)))
    return KWSNP_ERR_BADINVOC;
  NPCategoricalDensityMomentCtx * const categorical_moments =
    outer_pack_ctx != NULL ? outer_pack_ctx->categorical_density_moments : NULL;
  if(categorical_moments != NULL &&
     (num_reg_continuous != 0 || num_reg_unordered + num_reg_ordered <= 0 ||
      leave_one_out || drop_one_train || gather_scatter || do_score || do_ocg ||
      permutation_operator != OP_NOOP || kernel_pow != 1 || ncol_W != 0 ||
      ncol_Y > 1 || weighted_sum == NULL || categorical_moments->m2 == NULL))
    return KWSNP_ERR_BADINVOC;
  const int dual_ncol_Y =
    (do_dual_power && dual_power_ctx->matrix_Y != NULL) ?
    dual_power_ctx->ncol_Y : ncol_Y;
  const int dual_ncol_W =
    (do_dual_power && dual_power_ctx->matrix_W != NULL) ?
    dual_power_ctx->ncol_W : ncol_W;
  double ** const dual_matrix_Y =
    (do_dual_power && dual_power_ctx->matrix_Y != NULL) ?
    dual_power_ctx->matrix_Y : matrix_Y;
  double ** const dual_matrix_W =
    (do_dual_power && dual_power_ctx->matrix_W != NULL) ?
    dual_power_ctx->matrix_W : matrix_W;
  
  assert(!(do_score && do_ocg));
  assert(!(gather_scatter && is_adaptive));
  assert(!(do_dual_power && nws));
  if(do_hc0_derivative &&
     (!do_dual_power || permutation_operator != OP_DERIVATIVE ||
      p_nvar != num_reg_continuous + num_reg_unordered + num_reg_ordered ||
      p_nvar <= 0 ||
      leave_one_out || drop_one_train || gather_scatter ||
      hc0_derivative_ctx->num_continuous != num_reg_continuous ||
      hc0_derivative_ctx->num_predictors != p_nvar ||
      hc0_derivative_ctx->num_eval != num_obs_eval ||
      hc0_derivative_ctx->response_column_count !=
        MAX(ncol_Y, 1)*MAX(ncol_W, 1)))
    return KWSNP_ERR_BADINVOC;

  for(i = 0; (i < (num_reg_unordered + num_reg_ordered + num_reg_continuous)); i++){
    any_convolution |= (operator[i] == OP_CONVOLUTION);
  }
  for(i = 0; i < num_reg_continuous; ++i)
    all_continuous_normal &= operator[i] == OP_NORMAL;

  if(any_convolution && is_adaptive) np_ks_tree_use = 0;

  np_ks_tree_use &= (num_reg_continuous != 0);
  np_ks_tree_requested = np_ks_tree_use;

  int p_ipow = 0, * bpow = NULL;

  NL nls = {.node = NULL, .n = 0, .nalloc = 0};
  XL xl = {.istart = NULL, .nlev = NULL, .n = 0, .nalloc = 0};

  XL * pxl=  np_ks_tree_use ? &xl : NULL;
  XL * p_pxl =  NULL;

  if((np_ks_tree_use) && (p_nvar > 0)){
    p_pxl = (XL *)np_jksum_malloc_array_or_die(
      (size_t)p_nvar, sizeof(XL), "kernel weighted sum permutation ranges");
    for(k = 0; k < p_nvar; k++){
      p_pxl[k].istart = NULL;
      p_pxl[k].nlev = NULL;
      p_pxl[k].n = 0;
      p_pxl[k].nalloc = 0;
    }
  }

  // root node
  nls.node = (int *)np_jksum_malloc_array_or_die(1, sizeof(int), "kernel weighted sum nls.node");
  nls.node[0] = 0;
  nls.nalloc = nls.n = 1;

#ifdef MPI2
  // switch parallelisation strategies based on biggest stride
 
  int stride = MAX((int)ceil((double) num_obs_eval / (double) iNum_Processors),1);
  num_obs_eval_alloc = suppress_parallel ? num_obs_eval : (stride*iNum_Processors);

  int * igatherv = NULL, * idisplsv = NULL;

  igatherv = (int *)malloc(iNum_Processors*sizeof(int));
  if(igatherv == NULL) error("memory allocation failed");
  idisplsv = (int *)malloc(iNum_Processors*sizeof(int));
  if(idisplsv == NULL) error("memory allocation failed");
#else
  num_obs_eval_alloc = num_obs_eval;
#endif

  // we now allow one to use derivative, integral, convolution or standard kernels on a regressor by regressor basis
  int * KERNEL_reg_np = NULL, * KERNEL_unordered_reg_np = NULL, * KERNEL_ordered_reg_np = NULL;

  KERNEL_reg_np = (int *)malloc(sizeof(int)*num_reg_continuous);
  KERNEL_unordered_reg_np = (int *)malloc(sizeof(int)*num_reg_unordered);
  KERNEL_ordered_reg_np = (int *)malloc(sizeof(int)*num_reg_ordered);

  for(l = 0; l < num_reg_continuous; l++)
    KERNEL_reg_np[l] = KERNEL_reg[l] + OP_CFUN_OFFSETS[operator[l]];

  for(l = num_reg_continuous; l < (num_reg_continuous + num_reg_unordered); l++)
    KERNEL_unordered_reg_np[l - num_reg_continuous] = KERNEL_unordered_reg[l - num_reg_continuous] + OP_UFUN_OFFSETS[operator[l]];

  /*
   * Ordered standard, convolution, derivative, and integral operators share
   * the canonical table-offset dispatch. Focused categorical-profile and
   * ordered-distribution contracts protect the convolution/integral routes.
   */
  for(l = (num_reg_continuous+num_reg_unordered); l < (num_reg_continuous + num_reg_unordered + num_reg_ordered); l++)
    KERNEL_ordered_reg_np[l - (num_reg_continuous + num_reg_unordered)] = KERNEL_ordered_reg[l - (num_reg_continuous + num_reg_unordered)] + OP_OFUN_OFFSETS[operator[l]];

  /*
   * Validate categorical operator dispatch once at the map-construction
   * boundary. The vector microkernels are hot and may then index their
   * canonical tables directly without per-row fallback branches.
   */
  for(l = 0; l < num_reg_unordered; l++)
    if((KERNEL_unordered_reg_np[l] < 0) ||
       (KERNEL_unordered_reg_np[l] > 5))
      error("unsupported unordered kernel/operator combination");
  for(l = 0; l < num_reg_ordered; l++)
    if((KERNEL_ordered_reg_np[l] < 0) ||
       (KERNEL_ordered_reg_np[l] > 15))
      error("unsupported ordered kernel/operator combination");

  const int num_xt = is_adaptive?num_obs_eval:num_obs_train;
  const int progress_total = is_adaptive ? num_obs_train : num_obs_eval;
  size_t pkw_plane = 0;
  const int ws_step = is_adaptive? 0 :
                                 (MAX(ncol_Y, 1) * MAX(ncol_W, 1));
  const int dual_ws_step = is_adaptive ? 0 :
    (MAX(dual_ncol_Y, 1) * MAX(dual_ncol_W, 1));

  double *lambda = NULL, **matrix_bandwidth = NULL, **matrix_alt_bandwidth = NULL, **m = NULL;
  double *tprod = NULL, dband, *ws, *ws2 = NULL, * p_ws, * tprod_mp = NULL, * p_dband = NULL;
  double *perm_kbuf = NULL;
  double **bounded_cdf_lower_fixed = NULL;
  double **bounded_cdf_den_fixed = NULL;
  double *tree_convolution_alt_max = NULL;
  double *blas_Apack_owned = NULL;
  const double *blas_Apack = NULL;
  NP_TreeOuterBlasWorkspace tree_outer_workspace = {0};
  int use_disc_profile_cache = 0, disc_nprof = 0, disc_mark_token = 1;
  int disc_profile_from_override = 0;
  int disc_profile_from_global_cache = 0;
  int *disc_prof_id = NULL, *disc_prof_rep = NULL, *disc_prof_mark = NULL;
  int *disc_prof_list = NULL, *disc_active_idx = NULL;
  uint64_t *disc_prof_hash = NULL;
  double *disc_prof_val = NULL, *disc_ord_cl = NULL, *disc_ord_ch = NULL;
  int *disc_uno_const_ok = NULL, *disc_ord_const_ok = NULL;
  double *disc_uno_const = NULL, *disc_ord_const = NULL;
  int *cont_largeh_ok = NULL;
  double *cont_largeh_hmin = NULL, *cont_largeh_k0 = NULL;
  double cont_largeh_rel_tol = 1e-3;
  int *cont_largeh_active = NULL, *cont_largeh_active_fixed = NULL;
  int *fixed_tree_nonpruning_dimension = NULL;
  int *tree_active_dims = NULL;
  int cont_largeh_all_fixed = 0, cont_largeh_fixed_ready = 0;
  int cont_largeh_any_fixed = 0;
  int fused_gaussian_fixed_higher_convolution_eligible = 0;
#if NP_ACCEL_GAUSS_COMPILED
  /* 0 = none, 1 = ordinary Gaussian, 2 = higher-order Gaussian. */
  int fused_gaussian_product_eligible = 0;
  int fused_gaussian_convolution_kind = 0;
  double fused_gaussian_product_coef = 1.0;
#endif
  int tree_alllarge_bypass = 0;
  int cont_largeh_from_override = 0;
  int cont_largeh_from_global_cache = 0;
  int disc_uno_from_override = 0, disc_ord_from_override = 0;
  const int lean_reg_cont_loop =
    (!doscoreocg) &&
    (p_nvar == 0) &&
    (!int_cker_bound_extern) &&
    (!any_convolution);

  double * const * const xtc = is_adaptive?
    matrix_X_continuous_eval:matrix_X_continuous_train;
  double * const * const xtu = is_adaptive?
    matrix_X_unordered_eval:matrix_X_unordered_train;
  double * const * const xto = is_adaptive?
    matrix_X_ordered_eval:matrix_X_ordered_train;

  double * const * const xc = is_adaptive?
    matrix_X_continuous_train:matrix_X_continuous_eval;
  double * const * const xu = is_adaptive?
    matrix_X_unordered_train:matrix_X_unordered_eval;
  double * const * const xo = is_adaptive?
    matrix_X_ordered_train:matrix_X_ordered_eval;

  if (num_obs_eval == 0) {
    status = KWSNP_ERR_NOEVAL;
    goto cleanup;
  }

  do_psum = BANDWIDTH_reg == BW_ADAP_NN;
  swap_xxt = BANDWIDTH_reg == BW_ADAP_NN;
  /* Allocate memory for objects */

  mstep = (BANDWIDTH_reg==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_train:1);


  if(bandwidth_provided){
    if(BANDWIDTH_reg == BW_GEN_NN){
      matrix_bandwidth = matrix_bw_eval;
      if(any_convolution)
        matrix_alt_bandwidth = matrix_bw_train;
    }
    else if (is_adaptive){
      if (any_convolution){
        matrix_alt_bandwidth = matrix_bw_eval;        
      } 
      matrix_bandwidth = matrix_bw_train;
    } else {
      matrix_bandwidth = matrix_bw_train;
      if(any_convolution)
        matrix_alt_bandwidth = matrix_bw_train;
    }
    lambda = lambda_pre;
  } else {
    matrix_bandwidth = alloc_tmatd(mstep, num_reg_continuous);  
    lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);
    if(any_convolution && (BANDWIDTH_reg == BW_FIXED))
      matrix_alt_bandwidth = matrix_bandwidth;
  } 


  tprod = alloc_vecd((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_eval:num_obs_train);

  sum_element_length = MAX(ncol_Y, 1) * 
    MAX(ncol_W, 1);
  const int dual_sum_element_length =
    MAX(dual_ncol_Y, 1) * MAX(dual_ncol_W, 1);

  /* assert(!(BANDWIDTH_reg == BW_ADAP_NN)); */
  /* Conduct the estimation */

  /* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

  bpow = (int *) malloc(num_reg_continuous*sizeof(int));
  if(bpow == NULL){
    status = KWSNP_ERR_BADINVOC;
    goto cleanup;
  }

  for(i = 0; i < num_reg_continuous; i++){
    bpow[i] = (bandwidth_divide ? 1 : 0) + ((operator[i] == OP_DERIVATIVE) ? 1 : ((operator[i] == OP_INTEGRAL) ? -1 : 0));
  }

  if(p_nvar > 0){
    if(do_perm){
      permutation_kernel = (int *) malloc(num_reg_continuous*sizeof(int));
      for(i = 0; i < num_reg_continuous; i++){
        permutation_kernel[i] = KERNEL_reg[i] + OP_CFUN_OFFSETS[permutation_operator];
      }
    }
    if(do_score){
      ps_ukernel = (int *) malloc(num_reg_unordered*sizeof(int));
      ps_okernel = (int *) malloc(num_reg_ordered*sizeof(int));

      for(i = 0; i < num_reg_unordered; i++){
        ps_ukernel[i] = KERNEL_unordered_reg[i] +
          (operator[num_reg_continuous+i] == OP_CONVOLUTION ? 6 :
           OP_UFUN_OFFSETS[OP_DERIVATIVE]);
      }

      for(i = 0; i < num_reg_ordered; i++){
        const int op = operator[num_reg_continuous+num_reg_unordered+i];
        ps_okernel[i] = KERNEL_ordered_reg[i] +
          (op == OP_CONVOLUTION ? 16 : op == OP_INTEGRAL ? 20 :
           OP_OFUN_OFFSETS[OP_DERIVATIVE]);
      }
    } else if(do_ocg) {
      ps_ukernel = KERNEL_unordered_reg;
      ps_okernel = KERNEL_ordered_reg;
    }

    {
      const size_t tprod_rows = (size_t)((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_eval:num_obs_train);
      const size_t tprod_count = np_jksum_size_mul_or_die(tprod_rows, (size_t)p_nvar, "tprod_mp");
      tprod_mp = (double *)malloc(np_jksum_size_mul_or_die(tprod_count, sizeof(double), "tprod_mp"));
    }
    if(tprod_mp == NULL){
      status = KWSNP_ERR_BADINVOC;
      goto cleanup;
    }

    p_dband = (double *)malloc(p_nvar*sizeof(double));
    if(p_dband == NULL){
      status = KWSNP_ERR_BADINVOC;
      goto cleanup;
    }

    p_ipow = (bandwidth_divide ? 1 : 0) + ((permutation_operator == OP_DERIVATIVE) ? 1 : ((permutation_operator == OP_INTEGRAL) ? -1 : 0));

    perm_kbuf = (double *)malloc((size_t)num_xt*sizeof(double));
    if(perm_kbuf == NULL){
      status = KWSNP_ERR_BADINVOC;
      goto cleanup;
    }

  }

  if(!bandwidth_provided){
    if(kernel_bandwidth_mean_ctx((num_reg_continuous != 0) ? KERNEL_reg[0]: 0,
                             BANDWIDTH_reg,
                             num_obs_train,
                             num_obs_eval,
                             0,
                             0,
                             0,
                             num_reg_continuous,
                             num_reg_unordered,
                             num_reg_ordered,
                             suppress_parallel, // suppress_parallel if requested
                             vector_scale_factor,
                             matrix_X_continuous_train,				 /* Not used */
                             matrix_X_continuous_eval,				 /* Not used */
                             matrix_X_continuous_train,
                             matrix_X_continuous_eval,
                             matrix_bandwidth,						 /* Not used */
                             matrix_bandwidth,
                             lambda,
                             NULL, NULL, &nn_geometry_status)==1){

      status = nn_geometry_status == NP_NN_GEOMETRY_ZERO_RADIUS ?
        KWSNP_ERR_ZERO_NN_RADIUS : KWSNP_ERR_BADBW;
      goto cleanup;
    }
  }

#ifdef NP_ANN_DIRECT_DIAGNOSTIC
  if(np_ann_direct_expected_radius != NULL && BANDWIDTH_reg == BW_ADAP_NN) {
    np_ann_direct_radius_equal=
      num_reg_continuous == np_ann_direct_radius_columns &&
      num_obs_train == np_ann_direct_radius_rows;
    if(np_ann_direct_radius_equal)
      for(i=0;i<num_reg_continuous;++i)
        np_ann_direct_radius_equal &= memcmp(matrix_bandwidth[i],
          np_ann_direct_expected_radius[i],(size_t)num_obs_train*sizeof(double)) == 0;
    np_ann_direct_expected_radius=NULL;
  }
#endif

  if(!bandwidth_provided && any_convolution &&
     (num_reg_continuous > 0) &&
     (BANDWIDTH_reg == BW_GEN_NN)){
    if(np_gnn_convolution_training_bandwidth_prepare(
         (num_reg_continuous != 0) ? KERNEL_reg[0] : 0,
         num_obs_train, num_reg_unordered, num_reg_ordered,
         num_reg_continuous, suppress_parallel, vector_scale_factor,
         matrix_X_continuous_train, lambda,
         &matrix_alt_bandwidth, &nn_geometry_status) != 0){
      status = nn_geometry_status == NP_NN_GEOMETRY_ZERO_RADIUS ?
        KWSNP_ERR_ZERO_NN_RADIUS : KWSNP_ERR_BADBW;
      goto cleanup;
    }
  }

  if(!bandwidth_provided && (is_adaptive && any_convolution)){ // need additional bandwidths
    matrix_alt_bandwidth = alloc_tmatd(num_obs_eval, num_reg_continuous);  

    // Adaptive convolution requires an auxiliary BW_GEN_NN bandwidth matrix.
    nn_geometry_status = NP_NN_GEOMETRY_OK;
    if(kernel_bandwidth_mean_ctx((num_reg_continuous != 0) ? KERNEL_reg[0]: 0,
                             BW_GEN_NN, // this is not an error!
                             num_obs_train,
                             num_obs_eval,
                             0,
                             0,
                             0,
                             num_reg_continuous,
                             num_reg_unordered,
                             num_reg_ordered,
                             suppress_parallel, // suppress_parallel if requested
                             vector_scale_factor,
                             NULL,				 /* Not used */
                             NULL,				 /* Not used */
                             matrix_X_continuous_train,
                             matrix_X_continuous_eval,
                             NULL,						 /* Not used */
                             matrix_alt_bandwidth,
                             lambda,
                             NULL, NULL, &nn_geometry_status)==1){

      status = nn_geometry_status == NP_NN_GEOMETRY_ZERO_RADIUS ?
        KWSNP_ERR_ZERO_NN_RADIUS : KWSNP_ERR_BADBW;
      goto cleanup;
    }

  }

  if(np_ks_tree_use && BANDWIDTH_reg == BW_GEN_NN && any_convolution)
    tree_convolution_alt_max = np_tree_convolution_maxima(
      num_reg_continuous, KERNEL_reg_np, num_xt, matrix_alt_bandwidth);

  /*
   * A tree is an invocation-local arithmetic owner, not a persistent property
   * of the prepared data.  When fixed primary and requested permutation
   * kernels have full support over an exactly shared train/evaluation
   * geometry, the tree cannot remove any pair.  Select the canonical dense
   * sibling without changing the caller's tree readiness or data order.
   *
   * This admission rule is deliberately structural: it contains no sample-
   * size, host, elapsed-time, or optimizer-dependent crossover.  Unknown,
   * separate-evaluation, partial-tree, and genuinely prunable invocations
   * retain the existing owner.
   */
  fixed_tree_common_geometry =
    np_ks_tree_use &&
    BANDWIDTH_reg == BW_FIXED &&
    !do_partial_tree &&
    np_common_continuous_resident_geometry(
      num_obs_train, num_obs_eval, num_reg_continuous,
      matrix_X_continuous_train, matrix_X_continuous_eval);
  if(fixed_tree_common_geometry)
    fixed_tree_capability = np_fixed_tree_capability_from_common_root(
       BANDWIDTH_reg,
       num_reg_continuous,
       kdt,
       KERNEL_reg_np,
       do_perm ? permutation_kernel : NULL,
       do_perm ? bpso : NULL,
       matrix_bandwidth);
  if(fixed_tree_capability == NP_TREE_CANNOT_PRUNE){
    np_ks_tree_use = 0;
    pxl = NULL;
    free(p_pxl);
    p_pxl = NULL;
  }
  if(fixed_tree_common_geometry && np_ks_tree_use &&
     all_continuous_normal && !do_perm && (p_nvar == 0) &&
     !int_cker_bound_extern){
    int inactive_dimensions = 0;

    fixed_tree_nonpruning_dimension =
      (int *)calloc((size_t)num_reg_continuous, sizeof(int));
    if(fixed_tree_nonpruning_dimension != NULL){
      for(i = 0; i < num_reg_continuous; ++i){
        const NPTreeCapability coordinate_capability =
          np_fixed_tree_coordinate_capability_from_common_root(
            i, num_reg_continuous, kdt, KERNEL_reg_np, matrix_bandwidth);
        fixed_tree_nonpruning_dimension[i] =
          coordinate_capability == NP_TREE_CANNOT_PRUNE;
        inactive_dimensions += fixed_tree_nonpruning_dimension[i];
      }
      if(inactive_dimensions == 0){
        free(fixed_tree_nonpruning_dimension);
        fixed_tree_nonpruning_dimension = NULL;
      }
    }
    {
      const char * const dimension_oracle =
        getenv("NP_KWS_TREE_DIMENSION_ORACLE");
      if(dimension_oracle != NULL && dimension_oracle[0] == '1' &&
         dimension_oracle[1] == '\0')
        REprintf("NP_KWS_TREE_DIMENSION_ORACLE dimensions=%d inactive=%d active=%d\n",
                 num_reg_continuous,
                 inactive_dimensions,
                 num_reg_continuous - inactive_dimensions);
    }
  }
  {
    const char * const owner_oracle =
      getenv("NP_KWS_TREE_OWNER_ORACLE");
    if(owner_oracle != NULL && owner_oracle[0] == '1' &&
       owner_oracle[1] == '\0')
      REprintf("NP_KWS_TREE_OWNER_ORACLE requested=%d common=%d capability=%d owner=%s\n",
               np_ks_tree_requested,
               fixed_tree_common_geometry,
               (int)fixed_tree_capability,
               np_ks_tree_use ? "tree" : "dense");
  }

  if((BANDWIDTH_reg == BW_FIXED) && int_cker_bound_extern &&
     (num_reg_continuous > 0)){
    int need_bounded_cdf_cache = 0;
    for(i = 0; i < num_reg_continuous; i++)
      need_bounded_cdf_cache |= (operator[i] == OP_INTEGRAL) ||
        (do_perm && (permutation_operator == OP_INTEGRAL));
    if(need_bounded_cdf_cache){
      bounded_cdf_lower_fixed = alloc_tmatd(num_obs_train, num_reg_continuous);
      bounded_cdf_den_fixed = alloc_tmatd(num_obs_train, num_reg_continuous);
      for(i = 0; i < num_reg_continuous; i++){
        if((operator[i] == OP_INTEGRAL) ||
           (do_perm && (permutation_operator == OP_INTEGRAL))){
          const int k0 = KERNEL_reg_np[i] % 10;
          const double h = matrix_bandwidth[i][0];
          const double lb = vector_ckerlb_extern[i];
          const double ub = vector_ckerub_extern[i];
          for(int q = 0; q < num_obs_train; q++){
            const double train = xtc[i][q];
            const double zl = isfinite(lb) ? ((lb-train)/h) : R_NegInf;
            const double zu = isfinite(ub) ? ((ub-train)/h) : R_PosInf;
            bounded_cdf_lower_fixed[i][q] =
              np_cker_centered_primitive(k0, zl);
            bounded_cdf_den_fixed[i][q] = np_cker_interval_mass(k0, zl, zu);
          }
        }
      }
    }
  }

  if((!disable_gate_features) && (num_reg_continuous > 0)){
    if((gate_ctx->active == NP_GATE_CTX_OVERRIDE) &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->cont_ok != NULL &&
       gate_ctx->cont_hmin != NULL &&
       gate_ctx->cont_k0 != NULL){
      cont_largeh_ok = (int *)gate_ctx->cont_ok;
      cont_largeh_hmin = (double *)gate_ctx->cont_hmin;
      cont_largeh_k0 = (double *)gate_ctx->cont_k0;
      cont_largeh_from_override = 1;
    } else if(gate_ctx->active &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->cont_ok != NULL &&
       gate_ctx->cont_hmin != NULL &&
       gate_ctx->cont_k0 != NULL){
      cont_largeh_ok = (int *)gate_ctx->cont_ok;
      cont_largeh_hmin = (double *)gate_ctx->cont_hmin;
      cont_largeh_k0 = (double *)gate_ctx->cont_k0;
      cont_largeh_from_override = 1;
    }

    if(!cont_largeh_from_override){
      cont_largeh_rel_tol = np_cont_largeh_rel_tol();
      if(np_cont_largeh_cache_get_or_build(num_obs_train,
                                           num_obs_eval,
                                           num_reg_continuous,
                                           KERNEL_reg_np,
                                           matrix_X_continuous_train,
                                           matrix_X_continuous_eval,
                                           cont_largeh_rel_tol,
                                           &cont_largeh_ok,
                                           &cont_largeh_hmin,
                                           &cont_largeh_k0)){
        cont_largeh_from_global_cache = 1;
      } else {
        if(!np_cont_largeh_build_params(num_obs_train,
                                        num_obs_eval,
                                        num_reg_continuous,
                                        KERNEL_reg_np,
                                        matrix_X_continuous_train,
                                        matrix_X_continuous_eval,
                                        cont_largeh_rel_tol,
                                        &cont_largeh_ok,
                                        &cont_largeh_hmin,
                                        &cont_largeh_k0)){
          cont_largeh_ok = NULL;
          cont_largeh_hmin = NULL;
          cont_largeh_k0 = NULL;
        }
      }
    }
  }

  if((!disable_gate_features) && (num_reg_continuous > 0)){
    cont_largeh_active = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
    if(BANDWIDTH_reg == BW_FIXED){
      cont_largeh_active_fixed = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
      if(cont_largeh_active != NULL && cont_largeh_active_fixed != NULL){
        cont_largeh_all_fixed = 1;
        cont_largeh_any_fixed = 0;
        for(i = 0; i < num_reg_continuous; i++){
          const int is_largeh = (cont_largeh_ok != NULL) &&
            cont_largeh_ok[i] &&
            (operator[i] != OP_CONVOLUTION) &&
            isfinite(matrix_bandwidth[i][0]) &&
            (fabs(matrix_bandwidth[i][0]) >= cont_largeh_hmin[i]);

          cont_largeh_active_fixed[i] = is_largeh;
          cont_largeh_active[i] = is_largeh;
          cont_largeh_all_fixed &= is_largeh;
          cont_largeh_any_fixed |= is_largeh;
        }
        cont_largeh_fixed_ready = 1;
      }
    }
  }

  if(np_ks_tree_use &&
     cont_largeh_fixed_ready &&
     cont_largeh_all_fixed &&
     (p_nvar == 0) &&
     (!do_partial_tree)){
    /* All continuous kernels collapse to constants for every eval point.
       Bypass per-j tree list construction and use dense accumulation. */
    tree_alllarge_bypass = 1;
    pxl = NULL;
  }

  if(np_ks_tree_use && (!tree_alllarge_bypass) && (num_reg_continuous > 0)){
    tree_active_dims = (int *)malloc((size_t)num_reg_continuous*sizeof(int));
  }

  if(any_convolution &&
     (!doscoreocg) &&
     (!do_perm) &&
     (p_nvar == 0) &&
     (!int_cker_bound_extern) &&
     (BANDWIDTH_reg == BW_FIXED) &&
     (num_reg_continuous > 0)){
    fused_gaussian_fixed_higher_convolution_eligible = 1;
    for(i = 0; i < num_reg_continuous; i++){
      if((operator[i] != OP_CONVOLUTION) ||
         ((KERNEL_reg[i] != 1) && (KERNEL_reg[i] != 2))){
        fused_gaussian_fixed_higher_convolution_eligible = 0;
        break;
      }
    }
  }

#if NP_ACCEL_GAUSS_COMPILED
  if(lean_reg_cont_loop &&
     ((BANDWIDTH_reg == BW_FIXED) || (BANDWIDTH_reg == BW_GEN_NN) ||
      (BANDWIDTH_reg == BW_ADAP_NN)) &&
     (num_reg_continuous >= 2) &&
     (num_xt >= 256) &&
     (pxl == NULL) &&
     ((BANDWIDTH_reg != BW_FIXED) || (!cont_largeh_any_fixed)) &&
     np_mseries_accelerate_enabled_cache &&
     ((KERNEL_reg_np[0] >= 0) && (KERNEL_reg_np[0] <= 3)) &&
     np_accel_gauss_scratch_ensure(num_xt)){
    fused_gaussian_product_eligible =
      np_accel_gauss_product_kind(
        KERNEL_reg_np, num_reg_continuous, &fused_gaussian_product_coef);
  }

  if((num_reg_continuous >= 2) &&
     (BANDWIDTH_reg == BW_FIXED) &&
     any_convolution &&
     (!doscoreocg) &&
     (!do_perm) &&
     (p_nvar == 0) &&
     (!int_cker_bound_extern) &&
     (num_xt >= 256) &&
     np_mseries_accelerate_enabled_cache){
    fused_gaussian_convolution_kind = 1;
    for(i = 0; i < num_reg_continuous; i++){
      if((operator[i] != OP_CONVOLUTION) || (KERNEL_reg[i] != 0)){
        fused_gaussian_convolution_kind = 0;
        break;
      }
    }
  }

  if((num_reg_continuous >= 2) &&
     ((BANDWIDTH_reg == BW_GEN_NN) ||
      (BANDWIDTH_reg == BW_ADAP_NN)) &&
     any_convolution &&
     (!doscoreocg) &&
     (!do_perm) &&
     (p_nvar == 0) &&
     (!int_cker_bound_extern) &&
     (num_xt >= 256) &&
     np_mseries_accelerate_enabled_cache){
    fused_gaussian_convolution_kind = 2;
    for(i = 0; i < num_reg_continuous; i++){
      if((operator[i] != OP_CONVOLUTION) || (KERNEL_reg[i] != 0)){
        fused_gaussian_convolution_kind = 0;
        break;
      }
    }
  }
#endif

  if((!disable_gate_features) && (num_reg_unordered > 0)){
#if !NP_CF316_LITERAL
    if((gate_ctx->active == NP_GATE_CTX_OVERRIDE) &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->disc_uno_ok != NULL &&
       gate_ctx->disc_uno_const != NULL){
      disc_uno_const_ok = (int *)gate_ctx->disc_uno_ok;
      disc_uno_const = (double *)gate_ctx->disc_uno_const;
      disc_uno_from_override = 1;
    } else if(gate_ctx->active &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->disc_uno_ok != NULL &&
       gate_ctx->disc_uno_const != NULL){
      disc_uno_const_ok = (int *)gate_ctx->disc_uno_ok;
      disc_uno_const = (double *)gate_ctx->disc_uno_const;
      disc_uno_from_override = 1;
    }

#endif
    if(!disc_uno_from_override){
    disc_uno_const_ok = (int *)calloc((size_t)num_reg_unordered, sizeof(int));
    disc_uno_const = (double *)malloc((size_t)num_reg_unordered*sizeof(double));

    if(disc_uno_const_ok != NULL && disc_uno_const != NULL){
      double (* const ukf[])(int, double, int) = {
        np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
        np_score_uaa, np_score_unli_racine
      };
      const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));

      for(i = 0; i < num_reg_unordered; i++){
        const int ku = KERNEL_unordered_reg_np[i];
        const int ncat = (num_categories != NULL) ? num_categories[i] : 0;
        const double lam = lambda[i];

        disc_uno_const_ok[i] = 0;
        disc_uno_const[i] = 0.0;

        if(ku < 0 || ku >= nuk) continue;
        if(!NP_CF316_DISC_UPPER(ku, lam, ncat)) continue;

        const double kn_same = ukf[ku](1, lam, ncat);
        const double kn_diff = ukf[ku](0, lam, ncat);
        if(NP_CF316_DISC_CONSTANT(kn_same, kn_diff)){
          disc_uno_const_ok[i] = 1;
          disc_uno_const[i] = 0.5*(kn_same + kn_diff);
        }
      }
    } else {
      if(disc_uno_const_ok != NULL){ free(disc_uno_const_ok); disc_uno_const_ok = NULL; }
      if(disc_uno_const != NULL){ free(disc_uno_const); disc_uno_const = NULL; }
    }
    }
  }

  if((!disable_gate_features) && (num_reg_ordered > 0)){
#if !NP_CF316_LITERAL
    if((gate_ctx->active == NP_GATE_CTX_OVERRIDE) &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->disc_ord_ok != NULL &&
       gate_ctx->disc_ord_const != NULL){
      disc_ord_const_ok = (int *)gate_ctx->disc_ord_ok;
      disc_ord_const = (double *)gate_ctx->disc_ord_const;
      disc_ord_from_override = 1;
    } else if(gate_ctx->active &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->disc_ord_ok != NULL &&
       gate_ctx->disc_ord_const != NULL){
      disc_ord_const_ok = (int *)gate_ctx->disc_ord_ok;
      disc_ord_const = (double *)gate_ctx->disc_ord_const;
      disc_ord_from_override = 1;
    }

#endif
    if(!disc_ord_from_override){
    disc_ord_const_ok = (int *)calloc((size_t)num_reg_ordered, sizeof(int));
    disc_ord_const = (double *)malloc((size_t)num_reg_ordered*sizeof(double));

    if(disc_ord_const_ok != NULL && disc_ord_const != NULL){

      const int nok = 16;

      for(i = 0; i < num_reg_ordered; i++){
        const int oi = i + num_reg_unordered;
        const int ko = KERNEL_ordered_reg_np[i];
        const int ncat = (num_categories != NULL) ? num_categories[oi] : 0;
        const double lam = lambda[oi];

        disc_ord_const_ok[i] = 0;
        disc_ord_const[i] = 0.0;

        if(ko < 0 || ko >= nok) continue;
        if(ncat <= 0 || matrix_categorical_vals == NULL) continue;
        if(!NP_CF316_DISC_ORDERED_UPPER(ko, lam)) continue;

        double k0, k1;
        np_ordered_kernel_range(ko,lam,matrix_categorical_vals[oi],ncat,&k0,&k1);
        if(NP_CF316_DISC_CONSTANT(k0, k1)){
          disc_ord_const_ok[i] = 1;
          disc_ord_const[i] = 0.5*(k0 + k1);
        }
      }
    } else {
      if(disc_ord_const_ok != NULL){ free(disc_ord_const_ok); disc_ord_const_ok = NULL; }
      if(disc_ord_const != NULL){ free(disc_ord_const); disc_ord_const = NULL; }
    }
    }
  }
  
  if (leave_one_out && drop_one_train) {
    REprintf("\n error, leave one out estimator and drop-one estimator can't be enabled simultaneously");
    status = KWSNP_ERR_BADINVOC;
    goto cleanup;

  }

  if ((num_obs_train < (num_obs_eval + leave_one_out_offset)) && leave_one_out){
    
    REprintf("\nnumber of training points must be >= number of evaluation points to use the 'leave one out' estimator");
    status = KWSNP_ERR_BADINVOC;
    goto cleanup;
  }

  if(!gather_scatter && !nws)
    for(i = 0; i < num_obs_eval_alloc*sum_element_length; i++){
      weighted_sum[i] = 0.0;
    }

  if(!gather_scatter && do_dual_power)
    for(i = 0; i < num_obs_eval_alloc*dual_sum_element_length; i++){
      dual_power_ctx->weighted_sum[i] = 0.0;
    }

  if(categorical_leading != NULL)
    memset(categorical_leading->sum, 0,
           (size_t)num_obs_eval_alloc * sizeof(double));
  if(conditional_ratio != NULL)
    memset(conditional_ratio->sum, 0,
           (size_t)num_obs_eval_alloc * sizeof(double));

  if(categorical_moments != NULL) {
    memset(categorical_moments->m2, 0,
           (size_t)num_obs_eval_alloc * sizeof(double));
    categorical_moments->count = 0.0;
    if(is_adaptive) {
      memset(categorical_moments->anchor, 0, (size_t)num_obs_eval*sizeof(double));
      memset(categorical_moments->mean_delta, 0, (size_t)num_obs_eval*sizeof(double));
    }
  }

  if(!gather_scatter){
    for(i = 0; i < num_obs_eval_alloc*sum_element_length*p_nvar; i++){
      weighted_permutation_sum[i] = 0.0;
    }
  }

  if(do_hc0_derivative) {
    if(hc0_derivative_ctx->num_eval_alloc != num_obs_eval_alloc) {
      status = KWSNP_ERR_BADINVOC;
      goto cleanup;
    }
    for(i = 0; i < p_nvar; ++i) {
      if(hc0_derivative_ctx->gradient[i] == NULL ||
         hc0_derivative_ctx->gradient_stderr[i] == NULL) {
        status = KWSNP_ERR_BADINVOC;
        goto cleanup;
      }
      memset(hc0_derivative_ctx->gradient_stderr[i], 0,
             (size_t)num_obs_eval_alloc*sizeof(double));
      if(is_adaptive)
        memset(hc0_derivative_ctx->gradient[i], 0,
               (size_t)num_obs_eval_alloc*sizeof(double));
    }
  }

  if(pkw_output != NULL){
    if(pkw_output->data == NULL || pkw_output->block_count <= 0 ||
       pkw_output->block_count != p_nvar){
      status = KWSNP_ERR_BADINVOC;
      goto cleanup;
    }
    pkw_plane = np_jksum_size_mul_or_die(
      (size_t)progress_total, (size_t)num_xt,
      "permuted kernel-weight plane");
    pkw_output->plane_elements = pkw_plane;
    pkw_output->total_elements = np_jksum_size_mul_or_die(
      (size_t)pkw_output->block_count, pkw_plane,
      "derivative kernel-weight output");
  }

  if(np_ks_tree_use && pkw_plane > 0){
    for(i = 0; i < num_reg_continuous; i++){
      const int kernel = KERNEL_reg_np[i];
      const int p_kernel = (do_perm && bpso[i]) ? permutation_kernel[i] : kernel;
      tree_pkw_sparse =
        (fabs(cksup[kernel][0]) < 0.5*DBL_MAX) ||
        (fabs(cksup[kernel][1]) < 0.5*DBL_MAX) ||
        (fabs(cksup[p_kernel][0]) < 0.5*DBL_MAX) ||
        (fabs(cksup[p_kernel][1]) < 0.5*DBL_MAX);
      if(tree_pkw_sparse)
        break;
    }

  }

  if(pkw_output != NULL){
    pkw_output->state = tree_pkw_sparse ?
      NP_PKW_OUTPUT_TREE_SPARSE :
      (is_adaptive ?
       NP_PKW_OUTPUT_ADAPTIVE_FULL : NP_PKW_OUTPUT_DENSE_FULL);

    if(tree_pkw_sparse)
      memset(pkw_output->data, 0,
             np_jksum_size_mul_or_die(pkw_output->total_elements,
                                      sizeof(double),
                                      "derivative kernel-weight initialization"));
  }

  if (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN){
#ifdef MPI2
    if((!suppress_parallel) && (!gather_scatter)){
      js = stride * my_rank;
      je = MIN(num_obs_eval - 1, js + stride - 1);

      ws = nws ? NULL : weighted_sum + js*sum_element_length;
      ws2 = do_dual_power ?
        dual_power_ctx->weighted_sum + js*dual_sum_element_length : NULL;

      if(weighted_permutation_sum != NULL){
        p_ws = weighted_permutation_sum +js*sum_element_length;
      } else {
        p_ws = NULL;
      }

      // hopefully this now handles the case where there are more processors than observations
      for(i = 0, ii = 0; ii*stride < num_obs_eval; ii++){
        igatherv[ii] = stride*sum_element_length;
        idisplsv[ii] = i;
        i += stride*sum_element_length;
      }
      
      if(i < num_obs_eval*sum_element_length){
        const int de1 = (num_obs_eval - (ii-1)*stride)*sum_element_length;
        igatherv[ii] = de1;
        idisplsv[ii++] = i;
        i += de1;
      }

      for(; ii < iNum_Processors; ii++){
        igatherv[ii] = 0;
        idisplsv[ii] = i;
      }


    } else {
      js = 0;
      je = num_obs_eval - 1;
      ws = weighted_sum;
      ws2 = do_dual_power ? dual_power_ctx->weighted_sum : NULL;
      p_ws = weighted_permutation_sum;

      for(ii = 0; ii < iNum_Processors; ii++){
        igatherv[ii] = -1;
        idisplsv[ii] = -1;
      }
    }
#else
    js = 0;
    je = num_obs_eval - 1;
    ws = weighted_sum;
    ws2 = do_dual_power ? dual_power_ctx->weighted_sum : NULL;
    p_ws = weighted_permutation_sum;
#endif
    

  } else {
#ifdef MPI2
    if((!suppress_parallel) && (!gather_scatter)){
      js = stride * my_rank;
      je = MIN(num_obs_train - 1, js + stride - 1);
      ws = weighted_sum;
      ws2 = do_dual_power ? dual_power_ctx->weighted_sum : NULL;
      p_ws = weighted_permutation_sum;

      for(ii = 0; ii < iNum_Processors; ii++){
        igatherv[ii] = -1;
        idisplsv[ii] = -1;
      }

    } else {
      js = 0;
      je = num_obs_train - 1;
      ws = weighted_sum;
      ws2 = do_dual_power ? dual_power_ctx->weighted_sum : NULL;
      p_ws = weighted_permutation_sum;

      for(ii = 0; ii < iNum_Processors; ii++){
        igatherv[ii] = -1;
        idisplsv[ii] = -1;
      }

    }
#else
    js = 0;
    je = num_obs_train - 1;
    ws = weighted_sum;
    ws2 = do_dual_power ? dual_power_ctx->weighted_sum : NULL;
    p_ws = weighted_permutation_sum;
#endif
  }

  if((!disable_gate_features) && (!doscoreocg) && (p_nvar == 0) &&
     ((num_reg_unordered + num_reg_ordered) > 1) &&
     ((num_reg_ordered > 0) || ((num_reg_unordered + num_reg_ordered) > 2)) &&
     (num_xt >= 128)){
    if(gate_ctx->disc_prof_id != NULL &&
       gate_ctx->disc_prof_rep != NULL &&
       gate_ctx->disc_nprof > 0 &&
       gate_ctx->disc_prof_num_xt == num_xt){
      disc_prof_id = (int *)gate_ctx->disc_prof_id;
      disc_prof_rep = (int *)gate_ctx->disc_prof_rep;
      disc_nprof = gate_ctx->disc_nprof;
      disc_profile_from_override = 1;
      disc_profile_from_global_cache = 0;
    } else {
      disc_profile_from_override = 0;
      if(np_disc_profile_cache_get_or_build(num_xt,
                                            num_reg_unordered,
                                            num_reg_ordered,
                                            xtu,
                                            xto,
                                            &disc_prof_id,
                                            &disc_prof_rep,
                                            &disc_nprof)){
        disc_profile_from_global_cache = 1;
      } else {
        disc_prof_id = NULL;
        disc_prof_rep = NULL;
        disc_nprof = 0;
        disc_profile_from_global_cache = 0;
      }
    }

    if((disc_nprof > 0) && (4*disc_nprof <= 3*num_xt)){
      disc_prof_mark = (int *)calloc((size_t)disc_nprof, sizeof(int));
      disc_prof_list = (int *)malloc((size_t)disc_nprof*sizeof(int));
      disc_active_idx = (int *)malloc((size_t)num_xt*sizeof(int));
      disc_prof_val = (double *)malloc((size_t)disc_nprof*sizeof(double));

      if(num_reg_ordered > 0){
        disc_ord_cl = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
        disc_ord_ch = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
      }

      if(disc_prof_mark != NULL && disc_prof_list != NULL &&
         disc_active_idx != NULL && disc_prof_val != NULL &&
         ((num_reg_ordered == 0) || (disc_ord_cl != NULL && disc_ord_ch != NULL))){
        for(i = 0; i < num_reg_ordered; i++){
          const int oi = i + num_reg_unordered;
          disc_ord_cl[i] = (matrix_categorical_vals != NULL) ? matrix_categorical_vals[oi][0] : 0.0;
          disc_ord_ch[i] = (matrix_categorical_vals != NULL) ? matrix_categorical_vals[oi][num_categories[oi] - 1] : 0.0;
        }
        use_disc_profile_cache = 1;
      }
    }

    if(!use_disc_profile_cache){
      if(disc_prof_mark != NULL){ free(disc_prof_mark); disc_prof_mark = NULL; }
      if(disc_prof_list != NULL){ free(disc_prof_list); disc_prof_list = NULL; }
      if(disc_active_idx != NULL){ free(disc_active_idx); disc_active_idx = NULL; }
      if(disc_prof_val != NULL){ free(disc_prof_val); disc_prof_val = NULL; }
      if(disc_ord_cl != NULL){ free(disc_ord_cl); disc_ord_cl = NULL; }
      if(disc_ord_ch != NULL){ free(disc_ord_ch); disc_ord_ch = NULL; }
      if((disc_prof_id != NULL) && (!disc_profile_from_override) && (!disc_profile_from_global_cache)){ free(disc_prof_id); disc_prof_id = NULL; }
      if((disc_prof_rep != NULL) && (!disc_profile_from_override) && (!disc_profile_from_global_cache)){ free(disc_prof_rep); disc_prof_rep = NULL; }
      disc_nprof = 0;
    }
    if(use_disc_profile_cache && ((disc_prof_id == NULL) || (disc_prof_rep == NULL)))
      error("discrete profile cache setup failed");
  }

  const int leave_or_drop = leave_one_out || (drop_one_train && (BANDWIDTH_reg != BW_ADAP_NN));
  if(drop_one_train && (BANDWIDTH_reg != BW_ADAP_NN)) lod = drop_which_train;

  m = matrix_bandwidth;

  if((!nws) && (!do_psum) && (!gather_scatter) && (sgn == NULL) &&
     (!np_ks_tree_use) && (ncol_W > 0) &&
     np_outer_weighted_sum_blas_eligible(MAX(ncol_W, 1), MAX(ncol_Y, 1),
                                         num_xt, symmetric, NULL, NULL, NULL) &&
     (outer_pack_ctx != NULL) &&
     (outer_pack_ctx->data != NULL) &&
     (outer_pack_ctx->matrix_W == matrix_W) &&
     (outer_pack_ctx->ncol_W == ncol_W) &&
     (outer_pack_ctx->ncol_Y == ncol_Y) &&
     (outer_pack_ctx->num_weights == num_xt) &&
     (outer_pack_ctx->symmetric == symmetric))
    blas_Apack = outer_pack_ctx->data;

  if((blas_Apack == NULL) &&
     (!nws) && (!do_psum) && (!gather_scatter) && (sgn == NULL) &&
     (!np_ks_tree_use) && (ncol_W > 0) &&
     np_outer_weighted_sum_blas_eligible(MAX(ncol_W, 1), MAX(ncol_Y, 1),
                                         num_xt, symmetric, NULL, NULL, NULL)){
    blas_Apack_owned =
      np_outer_weighted_sum_pack_A(matrix_W, 1, MAX(ncol_W, 1),
                                   MAX(ncol_Y, 1), num_xt, symmetric);
    blas_Apack = blas_Apack_owned;
  }

  const int interrupt_total = (je >= js) ? (je - js + 1) : 0;
  const int interrupt_stride = np_hot_loop_interrupt_stride(interrupt_total);
  NP_KernelRowTileSink * const row_tile_sink =
    (outer_pack_ctx == NULL) ? NULL : outer_pack_ctx->row_tile_sink;

  if((row_tile_sink != NULL) &&
     ((row_tile_sink->consume == NULL) ||
      (row_tile_sink->weights == NULL) ||
      (row_tile_sink->weight_elements == 0U) ||
      (row_tile_sink->capacity_rows <= 0) ||
      (num_xt <= 0) ||
      ((size_t)row_tile_sink->capacity_rows >
       SIZE_MAX/(size_t)num_xt) ||
      ((size_t)row_tile_sink->capacity_rows*(size_t)num_xt >
       row_tile_sink->weight_elements) ||
      (row_tile_sink->count_rows != 0) ||
      (weighted_sum == NULL) ||
      (kw != row_tile_sink->weights) ||
      np_ks_tree_use ||
      (!suppress_parallel))){
    status = KWSNP_ERR_BADINVOC;
    goto cleanup;
  }

    /* do sums */
  for(j=js; j <= je; j++,
      ws = (ws==NULL ? NULL : ws+ws_step),
      ws2 = (ws2==NULL ? NULL : ws2+dual_ws_step),
      p_ws=(p_ws == NULL ? NULL : p_ws+ws_step)){
    if(((j - js) & 31) == 0)
      np_progress_bandwidth_loop_step();

    np_hot_loop_check_interrupt(j - js, interrupt_total, interrupt_stride);

    dband = 1.0;
    const int jbw = (BANDWIDTH_reg != BW_FIXED) ? j:0;
    int tprod_has_vals = 0;
    int deferred_const_active = 0;
    double deferred_const = 1.0;
    int all_cont_largeh = (num_reg_continuous > 0);
    int any_cont_largeh = 0;
    int tree_active_n = num_reg_continuous;
    const int tree_use_active_dims =
      np_ks_tree_use &&
      (!tree_alllarge_bypass) &&
      (p_nvar == 0) &&
      (!do_partial_tree) &&
      ((cont_largeh_active != NULL) ||
       (fixed_tree_nonpruning_dimension != NULL)) &&
      (tree_active_dims != NULL);

    if(cont_largeh_fixed_ready){
      all_cont_largeh = cont_largeh_all_fixed;
      any_cont_largeh = cont_largeh_any_fixed;
    } else if(cont_largeh_active != NULL){
      for(i = 0; i < num_reg_continuous; i++){
        const int is_largeh = (cont_largeh_ok != NULL) &&
          cont_largeh_ok[i] &&
          (operator[i] != OP_CONVOLUTION) &&
          isfinite(m[i][jbw]) &&
          (fabs(m[i][jbw]) >= cont_largeh_hmin[i]);

        cont_largeh_active[i] = is_largeh;
        all_cont_largeh &= is_largeh;
        any_cont_largeh |= is_largeh;
      }
    } else {
      all_cont_largeh = 0;
      any_cont_largeh = 0;
    }

    if(tree_use_active_dims){
      tree_active_n = 0;
      for(i = 0; i < num_reg_continuous; i++){
        const int approximate_constant =
          cont_largeh_active != NULL && cont_largeh_active[i];
        const int exact_nonpruning =
          fixed_tree_nonpruning_dimension != NULL &&
          fixed_tree_nonpruning_dimension[i];
        if(!approximate_constant && !exact_nonpruning){
          tree_active_dims[tree_active_n++] = i;
        }
      }
    }

    for (ii = 0; ii < p_nvar; ii++)
      p_dband[ii] = 1.0;

    // if we are consistently dropping one obs from training data, and we are doing a parallel sum, then that means
    // we need to skip here

    if(leave_one_out) lod = j + leave_one_out_offset;

    // do a hail mary, then generate the interaction list
    // anything but a fixed bandwidth is not yet supported
    // that includes convolutions 

    if(np_ks_tree_use && (!tree_alllarge_bypass)){

      if(kw != NULL)
        for(i = 0; i < num_xt; i++)
          tprod[i] = 0.0;

      double bb[kdt->ndim*2];

      // reset the interaction node list
      xl.n = 0;
      if(all_cont_largeh || (tree_use_active_dims && (tree_active_n == 0))){
        /* Global large-h: all continuous kernels are effectively K(0), so support is full. */
        merge_end_xl(pxl, &kdt->kdn[0]);
      } else {
#ifdef NP_ANN_DIRECT_DIAGNOSTIC
        if(conditional_ann != NULL) ++conditional_ann->traversals;
#endif
        if(!do_partial_tree){
          if(tree_use_active_dims && (tree_active_n < num_reg_continuous)){
            for(kk = 0; kk < tree_active_n; kk++){
              const int id = tree_active_dims[kk];
              const double sf = np_tree_convolution_scale(KERNEL_reg_np[id],
                adaptive_fold == NULL ? m[id][jbw] : adaptive_fold->successor[id][jbw],
                tree_convolution_alt_max, id);
              if(!is_adaptive){
                bb[2*id] = -cksup[KERNEL_reg_np[id]][1];
                bb[2*id+1] = -cksup[KERNEL_reg_np[id]][0];
              } else {
                bb[2*id] = cksup[KERNEL_reg_np[id]][0];
                bb[2*id+1] = cksup[KERNEL_reg_np[id]][1];
              }
              bb[2*id] = (fabs(bb[2*id]) == DBL_MAX) ? bb[2*id] : (xc[id][j] + bb[2*id]*sf);
              bb[2*id+1] = (fabs(bb[2*id+1]) == DBL_MAX) ? bb[2*id+1] : (xc[id][j] + bb[2*id+1]*sf);
              np_nn_rect_box(KERNEL_reg_np[id], BANDWIDTH_reg != BW_FIXED,
                             &bb[2*id], &bb[2*id+1]);
            }

            boxSearchNLPartial(kdt, &nls, bb, NULL, pxl, tree_active_dims, tree_active_n);
          } else {
            for(i = 0; i < num_reg_continuous; i++){
              const double sf = np_tree_convolution_scale(KERNEL_reg_np[i],
                adaptive_fold == NULL ? m[i][jbw] : adaptive_fold->successor[i][jbw],
                tree_convolution_alt_max, i);
              if(!is_adaptive){
                bb[2*i] = -cksup[KERNEL_reg_np[i]][1];
                bb[2*i+1] = -cksup[KERNEL_reg_np[i]][0];
              }else{
                bb[2*i] = cksup[KERNEL_reg_np[i]][0];
                bb[2*i+1] = cksup[KERNEL_reg_np[i]][1];
              }
              bb[2*i] = (fabs(bb[2*i]) == DBL_MAX) ? bb[2*i] : (xc[i][j] + bb[2*i]*sf);
              bb[2*i+1] = (fabs(bb[2*i+1]) == DBL_MAX) ? bb[2*i+1] : (xc[i][j] + bb[2*i+1]*sf);
              np_nn_rect_box(KERNEL_reg_np[i], BANDWIDTH_reg != BW_FIXED,
                             &bb[2*i], &bb[2*i+1]);
            }

            boxSearchNL(kdt, &nls, bb, NULL, pxl);
          }
        } else {
          for(i = 0; i < num_reg_continuous; i++){
            const double sf = np_tree_convolution_scale(KERNEL_reg_np[i],
              adaptive_fold == NULL ? m[i][jbw] : adaptive_fold->successor[i][jbw],
              tree_convolution_alt_max, i);
            if(!is_adaptive){
              bb[2*nld[i]] = -cksup[KERNEL_reg_np[i]][1];
              bb[2*nld[i]+1] = -cksup[KERNEL_reg_np[i]][0];
            }else{
              bb[2*nld[i]] = cksup[KERNEL_reg_np[i]][0];
              bb[2*nld[i]+1] = cksup[KERNEL_reg_np[i]][1];
            }

            bb[2*nld[i]] = (fabs(bb[2*nld[i]]) == DBL_MAX) ? bb[2*nld[i]] : (xc[i][j] + bb[2*nld[i]]*sf);
            bb[2*nld[i]+1] = (fabs(bb[2*nld[i]+1]) == DBL_MAX) ? bb[2*nld[i]+1] : (xc[i][j] + bb[2*nld[i]+1]*sf);
            np_nn_rect_box(KERNEL_reg_np[i], BANDWIDTH_reg != BW_FIXED,
                           &bb[2*nld[i]], &bb[2*nld[i]+1]);
          }
          if(idx == NULL)
            boxSearchNLPartial(kdt, inl, bb, NULL, pxl, nld, num_reg_continuous);
          else {
            boxSearchNLPartialIdx(kdt, inl, bb, NULL, pxl, nld, num_reg_continuous, idx);
          }
        }
      }

      if(do_perm && (permutation_operator == OP_INTEGRAL)){
        for(ii = 0, k = 0; ii < num_reg_continuous; ii++){
          if(bpso[ii]){
            // reset the interaction node list
            p_pxl[k].n = 0;

            if(all_cont_largeh){
              mirror_xl(pxl, p_pxl + k);
            } else if(!do_partial_tree){
              for(i = 0; i < num_reg_continuous; i++){
                const int knp = (i == ii) ? permutation_kernel[i] : KERNEL_reg_np[i];
                if(!is_adaptive){
                  bb[2*i] = -cksup[knp][1];
                  bb[2*i+1] = -cksup[knp][0];
                }else{
                  bb[2*i] = cksup[knp][0];
                  bb[2*i+1] = cksup[knp][1];
                }
                const double sf = np_tree_convolution_scale(knp,m[i][jbw],
                  tree_convolution_alt_max,i);
                bb[2*i] = (fabs(bb[2*i]) == DBL_MAX) ? bb[2*i] : (xc[i][j] + bb[2*i]*sf);
                bb[2*i+1] = (fabs(bb[2*i+1]) == DBL_MAX) ? bb[2*i+1] : (xc[i][j] + bb[2*i+1]*sf);
                np_nn_rect_box(knp, BANDWIDTH_reg != BW_FIXED,
                               &bb[2*i], &bb[2*i+1]);
              }

              boxSearchNL(kdt, &nls, bb, NULL, p_pxl + k);
            } else {
              for(i = 0; i < num_reg_continuous; i++){
                const int knp = (i == ii) ? permutation_kernel[i] : KERNEL_reg_np[i];
                if(!is_adaptive){
                  bb[2*nld[i]] = -cksup[knp][1];
                  bb[2*nld[i]+1] = -cksup[knp][0];
                }else{
                  bb[2*nld[i]] = cksup[knp][0];
                  bb[2*nld[i]+1] = cksup[knp][1];
                }
                const double sf = np_tree_convolution_scale(knp,m[i][jbw],
                  tree_convolution_alt_max,i);
                bb[2*nld[i]] = (fabs(bb[2*nld[i]]) == DBL_MAX) ? bb[2*nld[i]] : (xc[i][j] + bb[2*nld[i]]*sf);
                bb[2*nld[i]+1] = (fabs(bb[2*nld[i]+1]) == DBL_MAX) ? bb[2*nld[i]+1] : (xc[i][j] + bb[2*nld[i]+1]*sf);
                np_nn_rect_box(knp, BANDWIDTH_reg != BW_FIXED,
                               &bb[2*nld[i]], &bb[2*nld[i]+1]);
              }
              if(idx == NULL)
                boxSearchNLPartial(kdt, inl, bb, NULL, p_pxl + k, nld, num_reg_continuous);
              else{
                boxSearchNLPartialIdx(kdt, inl, bb, NULL, p_pxl + k, nld, num_reg_continuous, idx);
              }
            }

            k++;
          }
        }
        
        for(ii = k; ii < p_nvar; ii++){
          p_pxl[ii] = pxl[0];
        }
                  
      } else {
        for(ii = 0; ii < p_nvar; ii++){
          p_pxl[ii] = pxl[0];
        }
      } 

      /* No continuous support from the tree: all kernel products are zero. */
      if(xl.n == 0){
        int any_p_support = 0;
        for(ii = 0; ii < p_nvar; ii++){
          any_p_support |= (p_pxl[ii].n > 0);
        }

        if(!any_p_support){
          if(kw != NULL){
            for(i = 0; i < num_xt; i++)
              kw[j*num_xt + i] = 0.0;
          }
          continue;
        }
      }

    }
    /* continuous first */

    /* for the first iteration, no weights */
    /* for the rest, the accumulated products are the weights */
    if(adaptive_fold != NULL) {
      /* Fold requests have no permutations; initialize their common tail index. */
      for(i = 0, l = 0, k = 0; i < num_reg_continuous; ++i, ++l) {
        if(np_nn_fold_ckernelv(adaptive_fold, i, j, KERNEL_reg_np[i],
             xtc[i], tprod_has_vals, tprod, pxl, bandwidth_divide,
             vector_ckerlb_extern == NULL ? R_NegInf : vector_ckerlb_extern[i],
             vector_ckerub_extern == NULL ? R_PosInf : vector_ckerub_extern[i],
             &fold_primary, &fold_successor) != 0) {
          status = KWSNP_ERR_BADINVOC;
          goto cleanup;
        }
        tprod_has_vals = 1;
      }
    } else if(lean_reg_cont_loop){
#if NP_ACCEL_GAUSS_COMPILED
      int fused_gaussian_product = 0;

      if((fused_gaussian_product_eligible == 1) && (!any_cont_largeh))
        fused_gaussian_product =
          np_accel_gauss_product_try(xtc, xc, m, num_reg_continuous,
                                     num_xt, j, jbw,
                                     fused_gaussian_product_coef, tprod);
      else if((fused_gaussian_product_eligible == 2) && (!any_cont_largeh))
        fused_gaussian_product =
          np_accel_gauss_higher_product_try(
            KERNEL_reg_np, xtc, xc, m, num_reg_continuous,
            num_xt, j, jbw, fused_gaussian_product_coef, tprod);

      if(fused_gaussian_product){
        tprod_has_vals = 1;
        for(i = 0, l = 0, ip = 0, k = 0;
            i < num_reg_continuous;
            i++, l++, ip += do_perm){
          dband *= ipow(m[i][jbw], bpow[i]);
          k += bpso[l];
        }
      } else
#endif
      {
        for(i = 0, l = 0, ip = 0, k = 0;
            i < num_reg_continuous;
            i++, l++, ip += do_perm){
          const int use_largeh = any_cont_largeh &&
            (cont_largeh_active != NULL) ? cont_largeh_active[i] : 0;

          if(use_largeh){
            deferred_const *= cont_largeh_k0[i];
            deferred_const_active = 1;
          } else {
            np_ckernelv(KERNEL_reg_np[i], xtc[i], num_xt, tprod_has_vals,
                        xc[i][j], m[i][jbw], tprod, pxl, swap_xxt,
                        BANDWIDTH_reg != BW_FIXED, 1, 1.0,
                        0, R_NegInf, R_PosInf, NULL, NULL);
            tprod_has_vals = 1;
          }

          dband *= ipow(m[i][jbw], bpow[i]);
          k += bpso[l];
        }
      }
    } else {
      int fused_gaussian_convolution = 0;

#if NP_ACCEL_GAUSS_COMPILED
      if(fused_gaussian_convolution_kind != 0){
        if(fused_gaussian_convolution_kind == 1)
          fused_gaussian_convolution =
            np_accel_gauss_convolution_product_try(
              xtc, xc, m, matrix_alt_bandwidth, bpow,
              num_reg_continuous, num_xt, j, jbw, tprod);
        else
          fused_gaussian_convolution =
            np_accel_gauss_nn_convolution_product_try(
              xtc, xc, m, matrix_alt_bandwidth, bpow,
              num_reg_continuous, num_xt, j, jbw, tprod);
      } else
#endif
      if(fused_gaussian_fixed_higher_convolution_eligible)
        fused_gaussian_convolution =
          np_fixed_gaussian_convolution_product_try(
            KERNEL_reg, xtc, xc, m, matrix_alt_bandwidth, bpow,
            num_reg_continuous, num_xt, j, jbw, tprod);

      if(fused_gaussian_convolution){
        tprod_has_vals = 1;
        for(i = 0, l = 0, ip = 0, k = 0;
            i < num_reg_continuous;
            i++, l++, ip += do_perm){
          dband *= ipow(m[i][jbw], bpow[i]);
          k += bpso[l];
        }
      } else
      {
        for(i = 0, l = 0, ip = 0, k = 0;
            i < num_reg_continuous;
            i++, l++, ip += do_perm){
        const int kbase_i = KERNEL_reg_np[i] % 10;
        const int p_kbase_i = (do_perm ? permutation_kernel[i] : KERNEL_reg_np[i]) % 10;
        const int use_bounds_i = int_cker_bound_extern &&
          vector_ckerlb_extern != NULL &&
          vector_ckerub_extern != NULL &&
          (kbase_i >= 0 && kbase_i <= 9) &&
          ((isfinite(vector_ckerlb_extern[i]) && (fabs(vector_ckerlb_extern[i]) < 0.5*DBL_MAX)) ||
           (isfinite(vector_ckerub_extern[i]) && (fabs(vector_ckerub_extern[i]) < 0.5*DBL_MAX)));
        const int use_p_bounds_i = int_cker_bound_extern &&
          vector_ckerlb_extern != NULL &&
          vector_ckerub_extern != NULL &&
          (p_kbase_i >= 0 && p_kbase_i <= 9) &&
          ((isfinite(vector_ckerlb_extern[i]) && (fabs(vector_ckerlb_extern[i]) < 0.5*DBL_MAX)) ||
           (isfinite(vector_ckerub_extern[i]) && (fabs(vector_ckerub_extern[i]) < 0.5*DBL_MAX)));
        const int bounded_integral = use_bounds_i &&
          (operator[l] == OP_INTEGRAL);
        const int p_bounded_integral = use_p_bounds_i && do_perm &&
          (permutation_operator == OP_INTEGRAL);
        /* Fixed CVLS/AIC may omit the common continuous row mass before
         * symmetric pair reuse. Keep bounded support and owner selection. */
        const double invnorm = (use_bounds_i && !bounded_integral &&
          !(outer_pack_ctx != NULL && outer_pack_ctx->fixed_cv_cancel_row_mass)) ?
          np_cker_bounded_norm(KERNEL_reg_np[i], xc[i][j], m[i][jbw],
                               vector_ckerlb_extern[i],
                               vector_ckerub_extern[i]).inverse_mass : 1.0;
        const double p_invnorm = (use_p_bounds_i && !p_bounded_integral) ?
          np_cker_bounded_norm((do_perm ? permutation_kernel[i] : KERNEL_reg_np[i]),
                               xc[i][j], m[i][jbw],
                               vector_ckerlb_extern[i],
                               vector_ckerub_extern[i]).inverse_mass : 1.0;

        if(operator[l] != OP_CONVOLUTION){
          if(p_nvar == 0){
            const int use_largeh = !bounded_integral && any_cont_largeh &&
              (cont_largeh_active != NULL) ? cont_largeh_active[i] : 0;

            if(use_largeh){
              deferred_const *= cont_largeh_k0[i]*invnorm;
              deferred_const_active = 1;
            } else {
              np_ckernelv(KERNEL_reg_np[i], xtc[i], num_xt, tprod_has_vals,
                          xc[i][j], m[i][jbw], tprod, pxl, swap_xxt,
                          BANDWIDTH_reg != BW_FIXED, 1, invnorm, bounded_integral,
                          vector_ckerlb_extern[i], vector_ckerub_extern[i],
                          bounded_cdf_lower_fixed == NULL ? NULL :
                            bounded_cdf_lower_fixed[i],
                          bounded_cdf_den_fixed == NULL ? NULL :
                            bounded_cdf_den_fixed[i]);
              tprod_has_vals = 1;
            }
          } else {
            const int use_largeh = any_cont_largeh && (cont_largeh_active != NULL) ? cont_largeh_active[i] : 0;

            np_p_ckernelv(KERNEL_reg_np[i],
                          (do_perm ? permutation_kernel[i] : KERNEL_reg_np[i]),
                          k, p_nvar, xtc[i], num_xt, l, xc[i][j], m[i][jbw],
                          tprod, tprod_mp, pxl,
                          (p_pxl==NULL?NULL : p_pxl+k), swap_xxt,
                          BANDWIDTH_reg != BW_FIXED, bpso[l],
                          do_score, perm_kbuf, use_largeh,
                          (use_largeh ? cont_largeh_k0[i] : 0.0),
                          invnorm, p_invnorm, bounded_integral,
                          p_bounded_integral, vector_ckerlb_extern[i],
                          vector_ckerub_extern[i],
                          bounded_cdf_lower_fixed == NULL ? NULL :
                            bounded_cdf_lower_fixed[i],
                          bounded_cdf_den_fixed == NULL ? NULL :
                            bounded_cdf_den_fixed[i]);
          }
        }
        else if(p_nvar == 0){
          np_convol_ckernelv(KERNEL_reg[i], xtc[i], num_xt, tprod_has_vals, xc[i][j],
                             matrix_alt_bandwidth[i], (BANDWIDTH_reg == BW_FIXED),
                             m[i][jbw],
                             tprod, bpow[i]);
          tprod_has_vals = 1;
        } else {
          /* Every permutation plane inherits this coordinate's convolution
           * except the plane replacing it. Reuse the existing scratch row and
           * canonical vector owners; clearing tprod_mp would not compute it. */
          np_convol_ckernelv(KERNEL_reg[i], xtc[i], num_xt, 0, xc[i][j],
                            matrix_alt_bandwidth[i], (BANDWIDTH_reg == BW_FIXED),
                            m[i][jbw], perm_kbuf, bpow[i]);
          for(ii = 0; ii < num_xt; ++ii)
            tprod[ii] = l ? tprod[ii]*perm_kbuf[ii] : perm_kbuf[ii];
          for(kk = 0; kk < p_nvar; ++kk) {
            double * const plane = tprod_mp + (size_t)kk*num_xt;
            if(bpso[l] && kk == k && permutation_operator != OP_CONVOLUTION) {
              np_ckernelv(permutation_kernel[i], xtc[i], num_xt, l,
                          xc[i][j], m[i][jbw], plane,
                          p_pxl == NULL ? NULL : p_pxl + kk, swap_xxt,
                          BANDWIDTH_reg != BW_FIXED, 1, p_invnorm,
                          p_bounded_integral, vector_ckerlb_extern[i],
                          vector_ckerub_extern[i],
                          bounded_cdf_lower_fixed == NULL ? NULL :
                            bounded_cdf_lower_fixed[i],
                          bounded_cdf_den_fixed == NULL ? NULL :
                            bounded_cdf_den_fixed[i]);
            } else {
              for(ii = 0; ii < num_xt; ++ii)
                plane[ii] = l ? plane[ii]*perm_kbuf[ii] : perm_kbuf[ii];
            }
          }
        }
        dband *= ipow(m[i][jbw], bpow[i]);

        if(do_perm){
          for(ii = 0, kk = 0; ii < num_reg_continuous; ii++){
            if(bpso[ii]){
              if (i != ii){
                p_dband[kk] *= ipow(m[i][jbw], bpow[i]);              
              } else {
                p_dband[kk] *= ipow(m[i][jbw], p_ipow);
                if(((BANDWIDTH_reg == BW_FIXED) && (int_LARGE_SF == 0) && do_score)){
                  p_dband[kk] *= vector_scale_factor[ii]/(m[i][jbw]);
                }
              }
              kk++;
            }
          }
        }
        k += bpso[l];
      }
    }
    }


    for(ii = k; ii < p_nvar; ii++){
      p_dband[ii] = dband;
    }

    if(use_disc_profile_cache){
      int nactive = 0, nplist = 0;
      double disc_profile_const = 1.0;
      int disc_profile_const_active = 0;
      const int profile_base_empty = !tprod_has_vals;

      if(disc_uno_const_ok != NULL){
        for(kk = 0; kk < num_reg_unordered; kk++){
          if(!disc_uno_const_ok[kk]) continue;
          disc_profile_const *= disc_uno_const[kk];
          disc_profile_const_active = 1;
        }
      }

      if(disc_ord_const_ok != NULL){
        for(kk = 0; kk < num_reg_ordered; kk++){
          const int opidx = num_reg_continuous + num_reg_unordered + kk;
          if(!disc_ord_const_ok[kk]) continue;
          if(!(ps_ok_nli || (operator[opidx] != OP_CONVOLUTION))) continue;
          disc_profile_const *= disc_ord_const[kk];
          disc_profile_const_active = 1;
        }
      }

      if(profile_base_empty){
        if(pxl == NULL){
          for(i = 0; i < num_xt; i++)
            disc_active_idx[nactive++] = i;
        } else {
          for(ii = 0; ii < pxl->n; ii++){
            const int istart = pxl->istart[ii];
            const int iend = istart + pxl->nlev[ii];
            for(i = istart; i < iend; i++)
              disc_active_idx[nactive++] = i;
          }
        }
      } else if(pxl == NULL){
        for(i = 0; i < num_xt; i++){
          if(tprod[i] == 0.0) continue;
          disc_active_idx[nactive++] = i;
        }
      } else {
        for(ii = 0; ii < pxl->n; ii++){
          const int istart = pxl->istart[ii];
          const int iend = istart + pxl->nlev[ii];
          for(i = istart; i < iend; i++){
            if(tprod[i] == 0.0) continue;
            disc_active_idx[nactive++] = i;
          }
        }
      }

      if(disc_mark_token == INT_MAX){
        memset(disc_prof_mark, 0, (size_t)disc_nprof*sizeof(int));
        disc_mark_token = 1;
      } else {
        disc_mark_token++;
      }

      for(i = 0; i < nactive; i++){
        const int pid = disc_prof_id[disc_active_idx[i]];
        if(disc_prof_mark[pid] == disc_mark_token) continue;
        disc_prof_mark[pid] = disc_mark_token;
        disc_prof_list[nplist++] = pid;
      }

      double (* const ukf[])(int, double, int) = {
        np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
        np_score_uaa, np_score_unli_racine
      };

      for(i = 0; i < nplist; i++){
        const int pid = disc_prof_list[i];
        const int ridx = disc_prof_rep[pid];
        double dprod = 1.0;

        for(kk = 0; kk < num_reg_unordered; kk++){
          if((disc_uno_const_ok != NULL) && disc_uno_const_ok[kk]) continue;
          const int ku = KERNEL_unordered_reg_np[kk];
          dprod *= ukf[ku]((xtu[kk][ridx] == xu[kk][j]), lambda[kk], num_categories[kk]);
          if(dprod == 0.0) break;
        }

        if(dprod != 0.0){
          for(kk = 0; kk < num_reg_ordered; kk++){
            const int opidx = num_reg_continuous + num_reg_unordered + kk;
            if((disc_ord_const_ok != NULL) && disc_ord_const_ok[kk] &&
               (ps_ok_nli || (operator[opidx] != OP_CONVOLUTION)))
              continue;
            const int ko = KERNEL_ordered_reg_np[kk];
            const double c1 = swap_xxt ? xo[kk][j] : xto[kk][ridx];
            const double c2 = swap_xxt ? xto[kk][ridx] : xo[kk][j];
            dprod *= np_ordered_kernel_eval(ko,c1,c2,lambda[num_reg_unordered+kk],
              matrix_categorical_vals[num_reg_unordered+kk],num_categories[num_reg_unordered+kk],disc_ord_cl[kk],disc_ord_ch[kk]);
            if(dprod == 0.0) break;
          }
        }

        disc_prof_val[pid] = dprod;
      }

      if(profile_base_empty){
        for(i = 0; i < nactive; i++){
          const int idx_i = disc_active_idx[i];
          tprod[idx_i] = disc_prof_val[disc_prof_id[idx_i]];
        }
        tprod_has_vals = (nactive > 0);
      } else {
        for(i = 0; i < nactive; i++){
          const int idx_i = disc_active_idx[i];
          tprod[idx_i] *= disc_prof_val[disc_prof_id[idx_i]];
        }
      }

      if(disc_profile_const_active){
        deferred_const *= disc_profile_const;
        deferred_const_active = 1;
      }

      for(ii = 0; ii < (num_reg_unordered + num_reg_ordered); ii++){
        k += bpso[l + ii];
      }
      l += (num_reg_unordered + num_reg_ordered);
      ip += doscoreocg*(num_reg_unordered + num_reg_ordered);
    } else {
      /* unordered second */
      for(i=0; i < num_reg_unordered; i++, l++, ip += doscoreocg){
        if(p_nvar > 0){
          const int p_kernel = doscoreocg ?
            ps_ukernel[i] : KERNEL_unordered_reg_np[i];
          np_p_ukernelv(KERNEL_unordered_reg_np[i], p_kernel, k, p_nvar, xtu[i], num_xt, l, xu[i][j],
                        lambda[i], num_categories[i], matrix_categorical_vals[i][0], tprod, tprod_mp, pxl, (p_pxl == NULL ? NULL : p_pxl + k), swap_xxt, (bpso[l] ? do_ocg : 0), bpso[l], perm_kbuf);
        } else {
          if((p_nvar == 0) && (disc_uno_const_ok != NULL) && disc_uno_const_ok[i]){
            deferred_const *= disc_uno_const[i];
            deferred_const_active = 1;
          } else if(disc_uno_const_ok != NULL && disc_uno_const_ok[i]){
            np_ckernelv_mul_const(disc_uno_const[i], num_xt, l, tprod, pxl);
          } else {
            np_ukernelv(KERNEL_unordered_reg_np[i], xtu[i], num_xt, tprod_has_vals, xu[i][j], 
                        lambda[i], num_categories[i], tprod, pxl,
                        NP_CF316_SKIP_UPPER((disable_gate_features ? 1 : (disc_uno_const_ok != NULL))));
            tprod_has_vals = 1;
          }
        }
        k += bpso[l];
      }

      /* ordered third */
      for(i=0; i < num_reg_ordered; i++, l++, ip += doscoreocg){
        if(p_nvar == 0){
          if((disc_ord_const_ok != NULL) && disc_ord_const_ok[i] &&
             (ps_ok_nli || (operator[l] != OP_CONVOLUTION))){
            deferred_const *= disc_ord_const[i];
            deferred_const_active = 1;
          } else if((disc_ord_const_ok != NULL) && disc_ord_const_ok[i] &&
                    (ps_ok_nli || (operator[l] != OP_CONVOLUTION))){
            np_ckernelv_mul_const(disc_ord_const[i], num_xt, l, tprod, pxl);
          } else if(ps_ok_nli || (operator[l] != OP_CONVOLUTION)){
            np_okernelv(KERNEL_ordered_reg_np[i], xto[i], num_xt, tprod_has_vals,
                        xo[i][j], lambda[num_reg_unordered+i], 
                        (matrix_categorical_vals != NULL) ? matrix_categorical_vals[i+num_reg_unordered] : NULL, 
                        (num_categories != NULL) ? num_categories[i+num_reg_unordered] : 0,
                        tprod, pxl, swap_xxt);      
            tprod_has_vals = 1;
          } else {
            np_convol_okernelv(KERNEL_ordered_reg[i], xto[i], num_xt, tprod_has_vals,
                               xo[i][j], lambda[num_reg_unordered+i], 
                               num_categories[i+num_reg_unordered],
                               matrix_categorical_vals[i+num_reg_unordered],
                               tprod, swap_xxt);
            tprod_has_vals = 1;
          }
        } else {
          const int p_kernel = doscoreocg ?
            ps_okernel[i] : KERNEL_ordered_reg_np[i];
          const int this_do_ocg = bpso[l] ? do_ocg : 0;
          const int swapped_index = (this_do_ocg && !swap_xxt) ?
            matrix_ordered_indices[i][j] : 0;
          np_p_okernelv(KERNEL_ordered_reg_np[i], p_kernel, k, p_nvar, xto[i], num_xt, l,
                        xo[i][j], lambda[num_reg_unordered+i], 
                        (matrix_categorical_vals != NULL) ? matrix_categorical_vals[i+num_reg_unordered] : NULL, 
                        (num_categories != NULL) ? num_categories[i+num_reg_unordered] : 0,
                        tprod, tprod_mp, pxl, (p_pxl == NULL ? NULL : p_pxl + k), swap_xxt, this_do_ocg, bpso[l],
                        matrix_ordered_indices == NULL ? NULL : matrix_ordered_indices[i], swapped_index,
                        perm_kbuf);
        }
        k += bpso[l];
      }
    }
    if((p_nvar == 0) && deferred_const_active){
      np_ckernelv_mul_const(deferred_const, num_xt, tprod_has_vals, tprod, pxl);
      tprod_has_vals = 1;
    }

    /* expand matrix outer product, multiply by kernel weights, etc, do sum */

    if(fold != NULL && !is_adaptive)
      tprod[fold->eval_to_train[j]] = 0.0;

    if (!(drop_one_train && do_psum && (j == drop_which_train))){
      if(support_rank != NULL && support_rank->count < support_rank->cap) {
        /* Fixed/GNN: the product is the complete donor row. ANN: the
         * outer loop visits donors and its one inner entry is this query.
         * Tree scratch outside pxl is not part of the accumulated row. */
        const int ranges = pxl == NULL ? 1 : pxl->n;
        for(int range = 0; range < ranges &&
            support_rank->count < support_rank->cap; ++range) {
          const int begin = pxl == NULL ? 0 : pxl->istart[range];
          const int end = pxl == NULL ? num_xt :
            begin + pxl->nlev[range];
          for(int donor = begin; donor < end &&
              support_rank->count < support_rank->cap; ++donor)
            if(tprod[donor] != 0.0)
              support_rank->count += np_lp_design_support_add(
                support_rank->design, is_adaptive ? j : donor);
        }
      }
      if(!nws){
        np_outer_weighted_sum(matrix_W, sgn, ncol_W,
                              matrix_Y, ncol_Y,
                              tprod, num_xt,
                              leave_or_drop, lod,
                              kernel_pow,
                              do_psum, j,
                              symmetric,
                              gather_scatter,
                              1, dband,
                              ws, pxl,
                              outer_pack_ctx != NULL && outer_pack_ctx->tree_outer_blas,
                              blas_Apack,
                              &tree_outer_workspace);
      }

      if(residual_preparation != NULL)
        np_regression_residual_kernel_row(
          residual_preparation, BANDWIDTH_reg, j, tprod, num_xt, dband, pxl);
      if(uncertainty_context != NULL && uncertainty_context->unknown_count > 0 &&
         !np_regression_scalar_level_information(
           uncertainty_context, is_adaptive, j, num_xt, tprod, pxl))
        information_invalid = 1;

      if(categorical_moments != NULL)
        np_categorical_density_moment_row(
          categorical_moments, tprod, num_xt, j, is_adaptive);

      if(categorical_leading != NULL)
        np_categorical_leading_moment_row(
          categorical_leading, tprod, num_xt, j, dband, pxl,
          num_reg_unordered, num_reg_ordered, KERNEL_unordered_reg_np,
          KERNEL_ordered_reg_np, xtu, xto, xu, xo, lambda, num_categories,
          matrix_categorical_vals, disc_uno_const_ok, disc_uno_const,
          disc_ord_const_ok, disc_ord_const, use_disc_profile_cache,
          disc_prof_id, disc_prof_val);

      if(conditional_ratio != NULL)
        np_conditional_leading_ratio_row(
          conditional_ratio, tprod, num_xt, j, dband, *ws, pxl,
          num_reg_unordered, num_reg_ordered, KERNEL_unordered_reg_np,
          KERNEL_ordered_reg_np, xtu, xto, xu, xo, lambda, num_categories,
          matrix_categorical_vals, disc_uno_const_ok, disc_uno_const,
          disc_ord_const_ok, disc_ord_const, use_disc_profile_cache,
          disc_prof_id, disc_prof_val);

      if(do_dual_power){
        if(dual_power_ctx->observation_scale != NULL) {
          if(dual_power_ctx->kernel_pow != 2 || sgn != NULL ||
             leave_or_drop || symmetric || gather_scatter ||
             !np_regression_hc0_power2_outer_sum(
               dual_matrix_W, dual_ncol_W,
               dual_matrix_Y, dual_ncol_Y,
               tprod, num_xt, do_psum, j, dband,
               dual_power_ctx->observation_scale, ws2, pxl)) {
            status = KWSNP_ERR_BADINVOC;
            goto cleanup;
          }
        } else {
          np_outer_weighted_sum(dual_matrix_W, sgn, dual_ncol_W,
                                dual_matrix_Y, dual_ncol_Y,
                                tprod, num_xt,
                                leave_or_drop, lod,
                                dual_power_ctx->kernel_pow,
                                do_psum, j,
                                symmetric,
                                gather_scatter,
                                1, dband,
                                ws2, pxl,
                                outer_pack_ctx != NULL && outer_pack_ctx->tree_outer_blas,
                                blas_Apack,
                                &tree_outer_workspace);
        }
      }


        for(ii = 0; ii < p_nvar; ii++){
          np_outer_weighted_sum(matrix_W, sgn, ncol_W, 
                                matrix_Y, ncol_Y,
                                tprod_mp+ii*num_xt, num_xt,
                                leave_or_drop, lod,
                                kernel_pow,
                                do_psum, j,
                                symmetric,
                                gather_scatter,
                                1, p_dband[ii],
                                p_ws + ii*num_obs_eval*sum_element_length, (p_pxl == NULL) ? NULL : (p_pxl+ii),
                                outer_pack_ctx != NULL && outer_pack_ctx->tree_outer_blas,
                                blas_Apack,
                                &tree_outer_workspace);
        }
        if(conditional_ann != NULL) {
#ifdef NP_ANN_DIRECT_DIAGNOSTIC
          conditional_ann->tree_selected=np_ks_tree_use;
          conditional_ann->partial_tree=do_partial_tree;
#endif
          np_conditional_ann_direct_row(conditional_ann,j,tprod,dband,
                                        tprod_mp,p_dband,pxl,p_pxl);
        }
        if(do_hc0_derivative && !is_adaptive &&
           hc0_derivative_ctx->uncertainty_context != NULL &&
           hc0_derivative_ctx->uncertainty_context->unknown_count > 0 &&
           !np_regression_scalar_direction_information(
             hc0_derivative_ctx,j,num_xt,tprod,tprod_mp,dband,p_dband,ws,p_ws,pxl))
          information_invalid = 1;
        if(ann_information != NULL &&
           !np_regression_scalar_ann_information_row(
             ann_information,j,num_xt,tprod,tprod_mp,dband,p_dband,pxl))
          information_invalid = 1;
        if(scalar_contrast != NULL) {
          const int contrast_ok=is_adaptive ?
            np_regression_scalar_ann_contrast_row(scalar_contrast,j,num_xt,
              tprod,tprod_mp,dband,p_dband,pxl) :
            np_regression_scalar_contrast_row(scalar_contrast,j,num_xt,tprod,tprod_mp,pxl);
          if(!contrast_ok) information_invalid = 1;
        }
        if(do_hc0_derivative &&
           !np_regression_hc0_derivative_moments_accumulate(
             hc0_derivative_ctx, BANDWIDTH_reg, j, num_xt,
             tprod, dband, tprod_mp, p_dband,
             ws, p_ws, pxl)) {
          status = KWSNP_ERR_BADINVOC;
          goto cleanup;
        }

    }

    /* Convolution embeds its other-side bandwidth divisor in the product.
     * Sums and normalized private consumers above keep that representation;
     * raw exports must not inherit the sum's bandwidth.divide option. */
    if((kw != NULL || pkw_output != NULL) && !bandwidth_divide_weights &&
       bandwidth_divide && any_convolution) {
      for(int coordinate = 0, plane = 0; coordinate < num_reg_continuous;
          plane += bpso[coordinate], ++coordinate) {
        if(operator[coordinate] != OP_CONVOLUTION) continue;
        for(int row = 0; row < num_xt; ++row) {
          const double scale = matrix_alt_bandwidth[coordinate][
            BANDWIDTH_reg == BW_FIXED ? 0 : row];
          if(kw != NULL) tprod[row] *= scale;
          if(pkw_output != NULL)
            for(int block = 0; block < p_nvar; ++block)
              if(!bpso[coordinate] || block != plane)
                tprod_mp[(size_t)block*num_xt + row] *= scale;
        }
      }
    }

    if(kw != NULL){
      if(np_ks_tree_use && any_convolution &&
         !np_kernel_row_zero_outside_support(tprod,num_xt,pxl)) {
        status = KWSNP_ERR_BADINVOC;
        goto cleanup;
      }
      if(NP_UNLIKELY(row_tile_sink != NULL)){
        NP_KernelRowTileSink * const sink = row_tile_sink;
        double * const tile_row =
          sink->weights + (size_t)sink->count_rows*(size_t)num_xt;

        if(sink->count_rows == 0)
          sink->eval_start = j;
        if(bandwidth_divide_weights)
          for(i = 0; i < num_xt; i++)
            tile_row[i] = tprod[i]/dband;
        else
          for(i = 0; i < num_xt; i++)
            tile_row[i] = tprod[i];
        sink->count_rows++;

        if((sink->count_rows == sink->capacity_rows) || (j == je)){
          const double * const tile_weighted_sum =
            weighted_sum +
            (size_t)sink->eval_start*(size_t)sum_element_length;
          if(sink->consume(sink->state,
                           sink->eval_start,
                           sink->count_rows,
                           sink->weights,
                           num_xt,
                           tile_weighted_sum,
                           sum_element_length) != 0){
            status = KWSNP_ERR_BADINVOC;
            goto cleanup;
          }
          sink->count_rows = 0;
        }
      } else {
        // if using adaptive bandwidths, kw is returned transposed
        if(bandwidth_divide_weights)
          for(i = 0; i < num_xt; i++)
            kw[j*num_xt + i] = tprod[i]/dband;
        else
          for(i = 0; i < num_xt; i++)
            kw[j*num_xt + i] = tprod[i];
      }
    }

    if(pkw_output != NULL){
      for(ii = 0; ii < pkw_output->block_count; ii++){
        if(np_pkw_output_is_sparse(pkw_output) && (p_pxl != NULL)){
          const size_t dst_offset =
            (size_t)ii*pkw_plane +
            (size_t)j*(size_t)num_xt;
          const size_t src_offset = (size_t)ii*(size_t)num_xt;

          for(int support = 0; support < p_pxl[ii].n; support++){
            const int istart = p_pxl[ii].istart[support];
            const int nlev = p_pxl[ii].nlev[support];
            const int iend = istart + nlev;

            if(istart < 0 || nlev < 0 || iend < istart || iend > num_xt){
              status = KWSNP_ERR_BADINVOC;
              goto cleanup;
            }

            if(bandwidth_divide_weights){
              for(i = istart; i < iend; i++)
                pkw_output->data[dst_offset + (size_t)i] =
                  tprod_mp[src_offset + (size_t)i]/p_dband[ii];
            } else {
              for(i = istart; i < iend; i++)
                pkw_output->data[dst_offset + (size_t)i] =
                  tprod_mp[src_offset + (size_t)i];
            }
          }
        } else if(bandwidth_divide_weights)
          for(i = 0; i < num_xt; i++)
            pkw_output->data[(size_t)ii*pkw_plane +
                             (size_t)j*(size_t)num_xt +
                             (size_t)i] =
              tprod_mp[ii*num_xt + i]/p_dband[ii];
        else
          for(i = 0; i < num_xt; i++)
            pkw_output->data[(size_t)ii*pkw_plane +
                             (size_t)j*(size_t)num_xt +
                             (size_t)i] =
              tprod_mp[ii*num_xt + i];
      }
    }

    np_progress_fit_loop_step(j + 1, progress_total);
    
  }

  if((!gather_scatter) && (!suppress_parallel)){ 
    // gather_scatter is used only by the row-partitioned regression CV owner;
    // that owner excludes adaptive NN.
#ifdef MPI2
    if(!nws){
      if (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN){
        MPI_Allgather(MPI_IN_PLACE, stride * sum_element_length, MPI_DOUBLE, weighted_sum, stride * sum_element_length, MPI_DOUBLE, comm[1]);
      } else if(BANDWIDTH_reg == BW_ADAP_NN){
        MPI_Allreduce(MPI_IN_PLACE, weighted_sum, num_obs_eval*sum_element_length, MPI_DOUBLE, MPI_SUM, comm[1]);
      }
    }

    if(categorical_leading != NULL)
      MPI_Allgather(MPI_IN_PLACE, stride, MPI_DOUBLE,
                    categorical_leading->sum, stride, MPI_DOUBLE, comm[1]);
    if(conditional_ratio != NULL)
      MPI_Allgather(MPI_IN_PLACE, stride, MPI_DOUBLE,
                    conditional_ratio->sum, stride, MPI_DOUBLE, comm[1]);

    if(categorical_moments != NULL) {
      if(!is_adaptive) {
        MPI_Allgather(MPI_IN_PLACE, stride, MPI_DOUBLE,
                      categorical_moments->m2, stride, MPI_DOUBLE, comm[1]);
      } else {
        /* The original sum is only an anchor.  Keep its correction separate
         * to avoid manufacturing variance from a rounded constant mean. */
        const double local_count = categorical_moments->count;
        for(i = 0; i < num_obs_eval; ++i) {
          const double anchor = weighted_sum[i] / (double)num_obs_train;
          categorical_moments->merge_delta[i] = local_count == 0.0 ? 0.0 :
            local_count * ((categorical_moments->anchor[i] - anchor) +
                           categorical_moments->mean_delta[i]);
        }
        MPI_Allreduce(MPI_IN_PLACE, categorical_moments->merge_delta,
                      num_obs_eval, MPI_DOUBLE, MPI_SUM, comm[1]);
        if(local_count != 0.0)
          for(i = 0; i < num_obs_eval; ++i) {
            const double anchor = weighted_sum[i] / (double)num_obs_train;
            const double centered =
              ((categorical_moments->anchor[i] - anchor) +
               categorical_moments->mean_delta[i]) -
              categorical_moments->merge_delta[i] / (double)num_obs_train;
            categorical_moments->m2[i] += local_count * centered * centered;
          }
        MPI_Allreduce(MPI_IN_PLACE, categorical_moments->m2,
                      num_obs_eval, MPI_DOUBLE, MPI_SUM, comm[1]);
      }
    }

    if(do_dual_power){
      if (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN){
        MPI_Allgather(MPI_IN_PLACE, stride * dual_sum_element_length, MPI_DOUBLE,
                      dual_power_ctx->weighted_sum, stride * dual_sum_element_length,
                      MPI_DOUBLE, comm[1]);
      } else if(BANDWIDTH_reg == BW_ADAP_NN){
        MPI_Allreduce(MPI_IN_PLACE, dual_power_ctx->weighted_sum,
                      num_obs_eval*dual_sum_element_length, MPI_DOUBLE, MPI_SUM,
                      comm[1]);
      }
    }

    if(kw != NULL){
      MPI_Allgather(MPI_IN_PLACE, stride * num_xt, MPI_DOUBLE, kw, stride * num_xt, MPI_DOUBLE, comm[1]);          
    }

    if(p_nvar > 0){
      if (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN){
        for(ii = 0; ii < p_nvar; ii++){
          MPI_Allgatherv(MPI_IN_PLACE, igatherv[my_rank], MPI_DOUBLE, weighted_permutation_sum + ii*num_obs_eval*sum_element_length, igatherv, idisplsv, MPI_DOUBLE, comm[1]);
        }
      } else if(BANDWIDTH_reg == BW_ADAP_NN){
        MPI_Allreduce(MPI_IN_PLACE, weighted_permutation_sum, p_nvar*num_obs_eval*sum_element_length, MPI_DOUBLE, MPI_SUM, comm[1]);
      }
    }

    if(do_hc0_derivative) {
      if(BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN) {
        for(ii = 0; ii < p_nvar; ++ii)
          MPI_Allgather(
            MPI_IN_PLACE, stride, MPI_DOUBLE,
            hc0_derivative_ctx->gradient_stderr[ii], stride,
            MPI_DOUBLE, comm[1]);
      } else {
        const size_t plane_count =
          (size_t)p_nvar*(size_t)num_obs_eval_alloc;
        const int planes_are_adjacent =
          hc0_derivative_ctx->gradient_stderr[0] ==
            hc0_derivative_ctx->gradient[0] + plane_count;

        if(plane_count > (size_t)INT_MAX ||
           (planes_are_adjacent && plane_count > (size_t)INT_MAX/2U)) {
          status = KWSNP_ERR_BADINVOC;
          goto cleanup;
        }
        for(ii = 1; ii < p_nvar; ++ii) {
          if(hc0_derivative_ctx->gradient_stderr[ii] !=
             hc0_derivative_ctx->gradient_stderr[0] +
               (size_t)ii*(size_t)num_obs_eval_alloc ||
             hc0_derivative_ctx->gradient[ii] !=
               hc0_derivative_ctx->gradient[0] +
                 (size_t)ii*(size_t)num_obs_eval_alloc) {
            status = KWSNP_ERR_BADINVOC;
            goto cleanup;
          }
        }
        if(planes_are_adjacent) {
          MPI_Allreduce(MPI_IN_PLACE, hc0_derivative_ctx->gradient[0],
                        (int)(2U*plane_count), MPI_DOUBLE, MPI_SUM, comm[1]);
        } else {
          MPI_Allreduce(MPI_IN_PLACE, hc0_derivative_ctx->gradient[0],
                        (int)plane_count, MPI_DOUBLE, MPI_SUM, comm[1]);
          MPI_Allreduce(MPI_IN_PLACE, hc0_derivative_ctx->gradient_stderr[0],
                        (int)plane_count, MPI_DOUBLE, MPI_SUM, comm[1]);
        }
      }
    }
    if(scalar_contrast != NULL && !is_adaptive) {
      for(ii=scalar_contrast->categorical_range.begin;
          ii<scalar_contrast->categorical_range.end;++ii)
        MPI_Allgather(MPI_IN_PLACE,stride,MPI_DOUBLE,
                      scalar_contrast->gradient[ii],stride,MPI_DOUBLE,comm[1]);
    }
#endif
  }

  if(residual_preparation != NULL &&
     !np_regression_residual_kernel_reduce(
       residual_preparation, BANDWIDTH_reg, suppress_parallel))
    status = KWSNP_ERR_BADINVOC;
  if(ann_information != NULL &&
     !np_regression_scalar_ann_information_finish(ann_information,suppress_parallel))
    information_invalid = 1;
  if(scalar_contrast != NULL && is_adaptive &&
     !np_regression_scalar_ann_contrast_finish(scalar_contrast,suppress_parallel))
    information_invalid = 1;
  if(ann_information != NULL || scalar_contrast != NULL ||
     (uncertainty_context != NULL && uncertainty_context->unknown_count > 0)) {
#ifdef MPI2
    if(!suppress_parallel)
      MPI_Allreduce(MPI_IN_PLACE,&information_invalid,1,MPI_INT,MPI_MAX,comm[1]);
#endif
    if(information_invalid) status = KWSNP_ERR_BADINVOC;
  }

cleanup:
#ifdef MPI2
  free(igatherv);
  free(idisplsv);
#endif

  free(KERNEL_reg_np);
  free(KERNEL_unordered_reg_np);
  free(KERNEL_ordered_reg_np);
  
  if(!bandwidth_provided){
    free(lambda);
    free_tmat(matrix_bandwidth);
  }

  if(!bandwidth_provided && ((BANDWIDTH_reg == BW_ADAP_NN) && any_convolution))
    free_tmat(matrix_alt_bandwidth);
  if(!bandwidth_provided && (BANDWIDTH_reg == BW_GEN_NN) &&
     any_convolution && (num_reg_continuous > 0))
    free_tmat(matrix_alt_bandwidth);

  if(bounded_cdf_lower_fixed != NULL) free_tmat(bounded_cdf_lower_fixed);
  if(bounded_cdf_den_fixed != NULL) free_tmat(bounded_cdf_den_fixed);
  free(tree_convolution_alt_max);

  free(tprod);
  free(bpow);
  if(blas_Apack_owned != NULL) free(blas_Apack_owned);
  if(tree_outer_workspace.scratch != NULL)
    free(tree_outer_workspace.scratch);
  
  clean_xl(pxl);
  clean_nl(&nls);

  if(p_nvar > 0){
    free(tprod_mp);
    free(perm_kbuf);

    if(do_perm)
      free(permutation_kernel);

    if(do_score){
      free(ps_ukernel);
      free(ps_okernel);
    }

    if(p_pxl != NULL){
      if(permutation_operator == OP_INTEGRAL){
        for(l = 0, k = 0; l < num_reg_continuous; l++){
          if(bpso[l]){
            clean_xl(p_pxl+k);
            k++;
          }
        }
      }
      free(p_pxl);
    }

    free(p_dband);
  }

  if((disc_prof_id != NULL) && (!disc_profile_from_override) && (!disc_profile_from_global_cache)) free(disc_prof_id);
  if((disc_prof_rep != NULL) && (!disc_profile_from_override) && (!disc_profile_from_global_cache)) free(disc_prof_rep);
  if(disc_prof_hash != NULL) free(disc_prof_hash);
  if(disc_prof_mark != NULL) free(disc_prof_mark);
  if(disc_prof_list != NULL) free(disc_prof_list);
  if(disc_active_idx != NULL) free(disc_active_idx);
  if(disc_prof_val != NULL) free(disc_prof_val);
  if(disc_ord_cl != NULL) free(disc_ord_cl);
  if(disc_ord_ch != NULL) free(disc_ord_ch);
  if((disc_uno_const_ok != NULL) && (!disc_uno_from_override)) free(disc_uno_const_ok);
  if((disc_uno_const != NULL) && (!disc_uno_from_override)) free(disc_uno_const);
  if((disc_ord_const_ok != NULL) && (!disc_ord_from_override)) free(disc_ord_const_ok);
  if((disc_ord_const != NULL) && (!disc_ord_from_override)) free(disc_ord_const);
  if((cont_largeh_ok != NULL) && (!cont_largeh_from_override) && (!cont_largeh_from_global_cache)) free(cont_largeh_ok);
  if((cont_largeh_hmin != NULL) && (!cont_largeh_from_override) && (!cont_largeh_from_global_cache)) free(cont_largeh_hmin);
  if((cont_largeh_k0 != NULL) && (!cont_largeh_from_override) && (!cont_largeh_from_global_cache)) free(cont_largeh_k0);
  if(cont_largeh_active != NULL) free(cont_largeh_active);
  if(cont_largeh_active_fixed != NULL) free(cont_largeh_active_fixed);
  if(fixed_tree_nonpruning_dimension != NULL) free(fixed_tree_nonpruning_dimension);
  if(tree_active_dims != NULL) free(tree_active_dims);

  if(no_bpso)
    free(bpso);


  return(status);
}
