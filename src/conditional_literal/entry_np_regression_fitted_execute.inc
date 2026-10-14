#line 22600 "np.c"
static SEXP np_regression_fitted_execute(void *data)
{
  NPRegressionFittedCall *call = (NPRegressionFittedCall *)data;
  NPRegressionFittedOwner *owner = &call->owner;
  double * tuno = call->tuno;
  double * tord = call->tord;
  double * tcon = call->tcon;
  double * ty = call->ty;
  double * euno = call->euno;
  double * eord = call->eord;
  double * econ = call->econ;
  double * ey = call->ey;
  double * rbw = call->rbw;
  double * mcv = call->mcv;
  double * padnum = call->padnum;
  double * nconfac = call->nconfac;
  double * ncatfac = call->ncatfac;
  double * mysd = call->mysd;
  int * myopti = call->myopti;
  int * glp_degree = call->glp_degree;
  int * glp_gradient_order = call->glp_gradient_order;
  int * glp_bernstein = call->glp_bernstein;
  int * glp_basis = call->glp_basis;
  double * cm = call->cm;
  double * cmerr = call->cmerr;
  double * g = call->g;
  double *gerr = call->gerr;
  double * xtra = call->xtra;
  double * ckerlb = call->ckerlb;
  double * ckerub = call->ckerub;
  const NPContinuousKernelRoute *kernel_route = call->kernel_route;
  NPContinuousKernelDerivativeDiagnostics *kernel_route_diagnostics = call->kernel_route_diagnostics;
  int categorical_compress = call->categorical_compress;
  owner->active = 1;


  double * vector_scale_factor = NULL, * ecm = NULL, * ecmerr = NULL, ** eg = NULL, **egerr = NULL;
  double **eg_hc0_storage = NULL;
  double RS, MSE, MAE, MAPE, CORR, SIGN, pad_num;

  int i,j, num_var;
  int ey_is_ty, do_merr, do_grad, do_gerr;
  int train_is_eval, num_obs_eval_alloc, max_lev;
  int ordinary_hc0_active = 0;
  double *ordinary_hc0_residual = NULL;
  int *ordinary_hc0_information = NULL;
  NPRegressionHC0Context ordinary_hc0_context;

  int * ipt = NULL, * ipe = NULL;  // point permutation, see tree.c
  /* match integer options with their globals */

  num_reg_continuous_extern = myopti[REG_NCONI];
  num_reg_unordered_extern = myopti[REG_NUNOI];
  num_reg_ordered_extern = myopti[REG_NORDI];
  np_reset_y_side_extern();

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;
  const NPRegressionGradientRange gradient_range =
    np_regression_gradient_range(call->gradient_request, 0, num_var);

  train_is_eval = myopti[REG_TISEI];
  ey_is_ty = myopti[REG_EY];

  num_obs_train_extern = myopti[REG_TNOBSI];
  num_obs_eval_extern = myopti[REG_ENOBSI];

  if(train_is_eval && (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_regression): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    error("\n(np_regression): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
  }

  KERNEL_reg_extern = myopti[REG_CKRNEVI];
  KERNEL_reg_unordered_extern = myopti[REG_UKRNEVI];
  KERNEL_reg_ordered_extern = myopti[REG_OKRNEVI];

  vector_ckerlb_extern = ckerlb;
  vector_ckerub_extern = ckerub;
  int_cker_bound_extern = np_has_finite_cker_bounds(ckerlb, ckerub, num_reg_continuous_extern);

  int_LARGE_SF = myopti[REG_LSFI];
  int_MINIMIZE_IO = myopti[REG_MINIOI];
  BANDWIDTH_reg_extern = myopti[REG_BWI];

  do_merr = cmerr != NULL;
  do_grad = g != NULL;
  do_gerr = gerr != NULL;
  if(do_gerr && (!do_merr || !do_grad))
    error("C_np_regression: gradient standard errors require gradients and mean standard errors");
  if(do_grad != (myopti[REG_GRAD] != 0))
    error("C_np_regression: inconsistent gradient output request");
  np_lp_engine_extern = np_regression_engine_or_error(
    myopti[REG_LL], "C_np_regression");
  ordinary_hc0_active = do_merr &&
    (np_lp_engine_extern == NP_LP_ENGINE_SCALAR ||
     np_lp_engine_extern == NP_LP_ENGINE_GENERAL);
  vector_glp_degree_extern = glp_degree;
  vector_glp_gradient_order_extern = glp_gradient_order;
  int_glp_bernstein_extern = *glp_bernstein;
  int_glp_basis_extern = *glp_basis;

  max_lev = myopti[REG_MLEVI];
  pad_num = *padnum;

  nconfac_extern = *nconfac;
  ncatfac_extern = *ncatfac;

  int_TREE_X = myopti[REG_DOTREEI];
  int_TREE_PROFILE_X = myopti[REG_DOTREEI];

#ifdef MPI2
  num_obs_eval_alloc = MAX((int)ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif


  /* Allocate memory for objects */

  matrix_X_unordered_train_extern = owner->x_train[0] = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = owner->x_train[1] = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = owner->x_train[2] = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  vector_Y_extern = owner->y_train = alloc_vecd(num_obs_train_extern);

  if(!train_is_eval){
    matrix_X_unordered_eval_extern = owner->x_eval[0] = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
    matrix_X_ordered_eval_extern = owner->x_eval[1] = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
    matrix_X_continuous_eval_extern = owner->x_eval[2] = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);

    if(!ey_is_ty)
      vector_Y_eval_extern = owner->y_eval = alloc_vecd(num_obs_eval_extern);
    else
      vector_Y_eval_extern = NULL;

  } else {
    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;

    if(!ey_is_ty)
      vector_Y_eval_extern = owner->y_eval = alloc_vecd(num_obs_eval_extern);
    else
      vector_Y_eval_extern = vector_Y_extern;

  }

  ecm = owner->mean = alloc_vecd(num_obs_eval_alloc);
  if(do_merr)
    ecmerr = owner->mean_error = alloc_vecd(num_obs_eval_alloc);
  if(do_grad && do_gerr && ordinary_hc0_active) {
    eg_hc0_storage = owner->hc0_gradient = alloc_tmatd(num_obs_eval_alloc, 2*num_var);
    eg = eg_hc0_storage;
    egerr = eg_hc0_storage + num_var;
  } else if(do_grad) {
    eg = owner->gradient = alloc_matd(num_obs_eval_alloc, num_var);
  }
  if(do_gerr && !ordinary_hc0_active)
    egerr = owner->gradient_error = alloc_matd(num_obs_eval_alloc, num_var);

  num_categories_extern = owner->category_counts = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  vector_scale_factor = owner->scale = alloc_vecd(num_var + 1);
  matrix_categorical_vals_extern = owner->categories = alloc_matd(max_lev, num_reg_unordered_extern + num_reg_ordered_extern);

  vector_continuous_stddev_extern = mysd;
  /* train */

  for( j=0;j<num_reg_unordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+i];

  for( i=0;i<num_obs_train_extern;i++ )
    vector_Y_extern[i] = ty[i];

  /* eval */
  if(!train_is_eval){
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+i];
  }

  if (!ey_is_ty)
    for(i=0;i<num_obs_eval_extern;i++)
      vector_Y_eval_extern[i] = ey[i];

  /*  bandwidths/scale factors */

  for( i=0; i<num_var; i++ )
    vector_scale_factor[i+1] = rbw[i];

  /* fix up categories */

  for(j=0; j < (num_reg_unordered_extern + num_reg_ordered_extern); j++){
    i = 0;
    do {
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;
  }

  ipt = owner->training_order = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  if(!train_is_eval) {
    ipe = owner->evaluation_order = (int *)malloc(num_obs_eval_extern*sizeof(int));
    if(!(ipe != NULL))
      error("!(ipe != NULL)");

    for(i = 0; i < num_obs_eval_extern; i++){
      ipe[i] = i;
    }
  } else {
    ipe = ipt;
  }

  // attempt tree build, if enabled
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    if((BANDWIDTH_reg_extern != BW_ADAP_NN) || ((BANDWIDTH_reg_extern == BW_ADAP_NN) && train_is_eval)){
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern,
                   4*num_reg_continuous_extern, ipt, &owner->outer_tree);
      kdt_extern_X = owner->outer_tree;

      //put training data into tree-order using the index array

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+ipt[i]];


      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+ipt[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+ipt[i]];

      /* response variable */
      for( i=0;i<num_obs_train_extern;i++ )
        vector_Y_extern[i] = ty[ipt[i]];

      /* A separately supplied response at training locations owns the same
       * evaluation occurrences, hence the same tree permutation. */
      if(train_is_eval && !ey_is_ty)
        for(i = 0; i < num_obs_eval_extern; ++i)
          vector_Y_eval_extern[i] = ey[ipt[i]];

    } else {
      build_kdtree(matrix_X_continuous_eval_extern, num_obs_eval_extern, num_reg_continuous_extern,
                   4*num_reg_continuous_extern, ipe, &owner->outer_tree);
      kdt_extern_X = owner->outer_tree;

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+ipe[i]];

      if(!ey_is_ty)
        for(i=0;i<num_obs_eval_extern;i++)
          vector_Y_eval_extern[i] = ey[ipe[i]];

    }
  }


  /* Conduct estimation */

  /*
     nb - KERNEL_(|un)ordered_den are set to zero upon declaration
     - they have only one kernel type each at the moment
  */

  const NPNNGeometryContext nn_geometry_context = {
    .mode = train_is_eval ?
      NP_NN_QUERY_TRAINING_IDENTITY : NP_NN_QUERY_EXTERNAL,
    .eval_to_train = NULL
  };

  memset(&ordinary_hc0_context, 0, sizeof(ordinary_hc0_context));
  if(ordinary_hc0_active) {
    const NPNNGeometryContext training_geometry_context = {
      .mode = NP_NN_QUERY_TRAINING_IDENTITY,
      .eval_to_train = NULL
    };
    double *training_mean = NULL;
    double ** const training_gradient = train_is_eval ? eg : NULL;
    long double residual_scale = 0.0L;
    int invalid_preparation = 0;
    const int temporary_training_tree =
      int_TREE_X == NP_TREE_TRUE &&
      BANDWIDTH_reg_extern == BW_ADAP_NN &&
      !train_is_eval && kernel_route == NULL;
    KDT *outer_evaluation_kdt = NULL;
    NPRegressionHC0Context residual_preparation_context;
    NPRegressionResidualPreparation residual_preparation;
    int *training_identity = NULL;

    if(temporary_training_tree) {
      outer_evaluation_kdt = kdt_extern_X;
      build_kdtree(
        matrix_X_continuous_train_extern,
        num_obs_train_extern,
        num_reg_continuous_extern,
        4*num_reg_continuous_extern,
        ipt,
        &owner->training_tree);
      for(j = 0; j < num_reg_unordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_X_unordered_train_extern[j][i] =
            tuno[j*num_obs_train_extern + ipt[i]];
      for(j = 0; j < num_reg_ordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_X_ordered_train_extern[j][i] =
            tord[j*num_obs_train_extern + ipt[i]];
      for(j = 0; j < num_reg_continuous_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_X_continuous_train_extern[j][i] =
            tcon[j*num_obs_train_extern + ipt[i]];
      for(i = 0; i < num_obs_train_extern; i++)
        vector_Y_extern[i] = ty[ipt[i]];
      kdt_extern_X = owner->training_tree;
    }

    ordinary_hc0_residual =
      (double *)R_alloc((size_t)num_obs_train_extern, sizeof(double));
    ordinary_hc0_information =
      (int *)R_alloc((size_t)num_obs_train_extern, sizeof(int));
    training_identity =
      (int *)R_alloc((size_t)num_obs_train_extern, sizeof(int));
    training_mean = train_is_eval ? ecm :
      (double *)R_alloc((size_t)num_obs_train_extern, sizeof(double));
    for(i = 0; i < num_obs_train_extern; ++i) {
      training_identity[i] = i;
      ordinary_hc0_residual[i] = 0.0;
      ordinary_hc0_information[i] = NP_RESIDUAL_INVALID;
    }
    residual_preparation.response = vector_Y_extern;
    residual_preparation.evaluation_to_donor = training_identity;
    residual_preparation.num_obs_train = num_obs_train_extern;
    residual_preparation.num_obs_eval = num_obs_train_extern;
    residual_preparation.normalized_residual = ordinary_hc0_residual;
    residual_preparation.information = ordinary_hc0_information;
    memset(&residual_preparation_context, 0,
           sizeof(residual_preparation_context));
    residual_preparation_context.donor_to_canonical = ipt;
    residual_preparation_context.num_obs_train = num_obs_train_extern;
    residual_preparation_context.status =
      NP_REGRESSION_HC0_RESIDUAL_PREPARING;
    residual_preparation_context.gradient_request = call->gradient_request;
    residual_preparation_context.preparation = &residual_preparation;

    kernel_estimate_regression_categorical_tree_np(
      np_lp_engine_extern,
      KERNEL_reg_extern,
      KERNEL_reg_unordered_extern,
      KERNEL_reg_ordered_extern,
      BANDWIDTH_reg_extern,
      num_obs_train_extern,
      num_obs_train_extern,
      num_reg_unordered_extern,
      num_reg_ordered_extern,
      num_reg_continuous_extern,
      matrix_X_unordered_train_extern,
      matrix_X_ordered_train_extern,
      matrix_X_continuous_train_extern,
      matrix_X_unordered_train_extern,
      matrix_X_ordered_train_extern,
      matrix_X_continuous_train_extern,
      vector_Y_extern,
      vector_Y_extern,
      &vector_scale_factor[1],
      num_categories_extern,
      matrix_categorical_vals_extern,
      training_mean,
      training_gradient,
      NULL,
      NULL,
      &RS,
      &MSE,
      &MAE,
      &MAPE,
      &CORR,
      &SIGN,
      kernel_route,
      kernel_route_diagnostics,
      categorical_compress,
      NP_REGRESSION_STDERR_LOCAL_RESIDUAL,
      1, /* ordinary regression response, including SE preparation */
      NULL,
      &training_geometry_context,
      &residual_preparation_context, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if(temporary_training_tree) {
      kdt_extern_X = outer_evaluation_kdt;
      free_kdtree(&owner->training_tree);
    }

    for(i = 0; i < num_obs_train_extern; i++) {
      const long double residual = (long double)ordinary_hc0_residual[i];
      const long double magnitude = fabsl(residual);

      if(!isfinite((double)vector_Y_extern[i]) ||
         !isfinite((double)training_mean[i]) || !isfinite(magnitude) ||
         (ordinary_hc0_information[i] != NP_RESIDUAL_IDENTIFIED &&
          ordinary_hc0_information[i] != NP_RESIDUAL_UNIDENTIFIED))
        invalid_preparation = 1;
      if(magnitude > residual_scale)
        residual_scale = magnitude;
    }

    /* A failed/incomplete producer must not strand a peer in the subsequent
     * covariance traversal. Local fanout work does not enter collectives. */
#ifdef MPI2
    if(!np_mpi_local_regression_active() && comm != NULL &&
       comm[1] != MPI_COMM_NULL && iNum_Processors > 1)
      MPI_Allreduce(MPI_IN_PLACE, &invalid_preparation, 1, MPI_INT,
                    MPI_MAX, comm[1]);
#endif
    if(invalid_preparation)
      error("ordinary-regression normalized residual preparation is invalid or incomplete");

    if(residual_scale > (long double)DBL_MAX)
      error("ordinary-regression HC0 residual scale exceeds double range");

    ordinary_hc0_context.scaled_residual = ordinary_hc0_residual;
    ordinary_hc0_context.residual_information = ordinary_hc0_information;
    for(i = 0; i < num_obs_train_extern; ++i)
      if(ordinary_hc0_information[i] == NP_RESIDUAL_UNIDENTIFIED)
        ++ordinary_hc0_context.unknown_count;
    ordinary_hc0_context.donor_to_canonical = ipt;
    ordinary_hc0_context.num_obs_train = num_obs_train_extern;
    ordinary_hc0_context.residual_scale = (double)residual_scale;
    ordinary_hc0_context.point_already_computed = train_is_eval;
    ordinary_hc0_context.gradient_request = call->gradient_request;

    if(residual_scale == 0.0L) {
      ordinary_hc0_context.status = NP_REGRESSION_HC0_RESIDUAL_ALL_ZERO;
      for(i = 0; i < num_obs_train_extern; i++)
        ordinary_hc0_residual[i] = 0.0;
    } else {
      ordinary_hc0_context.status = NP_REGRESSION_HC0_RESIDUAL_READY;
      for(i = 0; i < num_obs_train_extern; i++) {
        ordinary_hc0_residual[i] =
          (double)((long double)ordinary_hc0_residual[i] / residual_scale);
        if(!R_FINITE(ordinary_hc0_residual[i]))
          error("ordinary-regression HC0 scaled residual is not finite");
      }
    }

    if(temporary_training_tree) {
      for(i = 0; i < num_obs_train_extern; i++) {
        while(ipt[i] != i) {
          const int destination = ipt[i];
          const double residual = ordinary_hc0_residual[i];
          const int information = ordinary_hc0_information[i];
          const int permutation = ipt[i];

          ordinary_hc0_residual[i] =
            ordinary_hc0_residual[destination];
          ordinary_hc0_residual[destination] = residual;
          ordinary_hc0_information[i] =
            ordinary_hc0_information[destination];
          ordinary_hc0_information[destination] = information;
          ipt[i] = ipt[destination];
          ipt[destination] = permutation;
        }
      }
      for(j = 0; j < num_reg_unordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_X_unordered_train_extern[j][i] =
            tuno[j*num_obs_train_extern + i];
      for(j = 0; j < num_reg_ordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_X_ordered_train_extern[j][i] =
            tord[j*num_obs_train_extern + i];
      for(j = 0; j < num_reg_continuous_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_X_continuous_train_extern[j][i] =
            tcon[j*num_obs_train_extern + i];
      for(i = 0; i < num_obs_train_extern; i++)
        vector_Y_extern[i] = ty[i];
    }

    /*
     * A training-point HC0 pass may reuse the internal continuous-gradient
     * planes as adaptive-NN cross moments.  Preserve the already-computed
     * point estimates directly in their final caller-owned storage first;
     * this moves the ordinary output copy rather than allocating another
     * O(n p) workspace.
     */
    if(train_is_eval && do_grad && do_gerr)
      for(j = gradient_range.begin; j < gradient_range.end; ++j)
        for(i = 0; i < num_obs_eval_extern; ++i)
          g[j*num_obs_eval_extern + ipe[i]] = eg[j][i];

    np_progress_fit_set_offset(num_obs_train_extern);
  }

  /* Dependency flags are separate from the finite partial covariance sums.
   * They exist only when a requested covariance can encounter unknown donor
   * variance; never allocate them on the uncertainty-off/common known path. */
  if(ordinary_hc0_active && ordinary_hc0_context.unknown_count > 0) {
    ordinary_hc0_context.mean_unavailable =
      (int *)R_alloc((size_t)num_obs_eval_extern, sizeof(int));
    memset(ordinary_hc0_context.mean_unavailable, 0,
           (size_t)num_obs_eval_extern * sizeof(int));
    if(do_gerr) {
      ordinary_hc0_context.gradient_unavailable =
        (int **)R_alloc((size_t)num_var, sizeof(int *));
      ordinary_hc0_context.gradient_structural_zero =
        (int **)R_alloc((size_t)num_var, sizeof(int *));
      memset(ordinary_hc0_context.gradient_unavailable, 0,
             (size_t)num_var * sizeof(int *));
      memset(ordinary_hc0_context.gradient_structural_zero, 0,
             (size_t)num_var * sizeof(int *));
      *call->gradient_zero_out = (int *)R_alloc(
        (size_t)num_var * num_obs_eval_extern, sizeof(int));
      memset(*call->gradient_zero_out, 0,
             (size_t)num_var * num_obs_eval_extern * sizeof(int));
      for(j = gradient_range.begin; j < gradient_range.end; ++j) {
        ordinary_hc0_context.gradient_unavailable[j] =
          (int *)R_alloc((size_t)num_obs_eval_extern, sizeof(int));
        memset(ordinary_hc0_context.gradient_unavailable[j], 0,
               (size_t)num_obs_eval_extern * sizeof(int));
        ordinary_hc0_context.gradient_structural_zero[j] =
          (int *)R_alloc((size_t)num_obs_eval_extern, sizeof(int));
        memset(ordinary_hc0_context.gradient_structural_zero[j], 0,
               (size_t)num_obs_eval_extern * sizeof(int));
      }
    }
  }

  if(call->fit_projection != NULL) {
    call->fit_projection->donor_map = ipt;
    call->fit_projection->eval_map = ipe;
  }
  if(call->empty_rows != NULL)
    call->empty_rows->row_map = ipe;
  kernel_estimate_regression_categorical_tree_np(np_lp_engine_extern,
                                                   KERNEL_reg_extern,
                                                   KERNEL_reg_unordered_extern,
                                                   KERNEL_reg_ordered_extern,
                                                   BANDWIDTH_reg_extern,
                                                   num_obs_train_extern,
                                                   num_obs_eval_extern,
                                                   num_reg_unordered_extern,
                                                   num_reg_ordered_extern,
                                                   num_reg_continuous_extern,
                                                   /* Train */
                                                   matrix_X_unordered_train_extern,
                                                   matrix_X_ordered_train_extern,
                                                   matrix_X_continuous_train_extern,
                                                   /* Eval */
                                                   matrix_X_unordered_eval_extern,
                                                   matrix_X_ordered_eval_extern,
                                                   matrix_X_continuous_eval_extern,
                                                   vector_Y_extern,
                                                   vector_Y_eval_extern,
                                                   &vector_scale_factor[1],
                                                   num_categories_extern,
                                                   matrix_categorical_vals_extern,
                                                   ecm,
                                                   eg,
                                                   ecmerr,
                                                   egerr,
                                                   &RS,
                                                   &MSE,
                                                   &MAE,
                                                   &MAPE,
                                                   &CORR,
                                                   &SIGN,
                                                   kernel_route,
                                                   kernel_route_diagnostics,
                                                   categorical_compress,
                                                   NP_REGRESSION_STDERR_LOCAL_RESIDUAL,
                                                   1, /* ordinary regression response */
                                                   NULL,
                                                   &nn_geometry_context,
                                                   ordinary_hc0_active ?
                                                     &ordinary_hc0_context : NULL,
                                                   call->empty_rows, NULL, NULL, NULL, NULL, NULL, call->fit_projection);


  if(ordinary_hc0_active && ordinary_hc0_context.unknown_count > 0) {
#ifdef MPI2
    const int distributed = !np_mpi_local_regression_active() &&
      comm != NULL && comm[1] != MPI_COMM_NULL && iNum_Processors > 1;
    if(distributed) {
      MPI_Allreduce(MPI_IN_PLACE, ordinary_hc0_context.mean_unavailable,
                    num_obs_eval_extern, MPI_INT, MPI_MAX, comm[1]);
      if(do_gerr)
        for(j = gradient_range.begin; j < gradient_range.end; ++j) {
          MPI_Allreduce(MPI_IN_PLACE,
                        ordinary_hc0_context.gradient_unavailable[j],
                        num_obs_eval_extern, MPI_INT, MPI_MAX, comm[1]);
          MPI_Allreduce(MPI_IN_PLACE,
                        ordinary_hc0_context.gradient_structural_zero[j],
                        num_obs_eval_extern, MPI_INT, MPI_MAX, comm[1]);
        }
    }
#endif
    for(i = 0; i < num_obs_eval_extern; ++i) {
      if(ordinary_hc0_context.mean_unavailable[i]) {
        ecmerr[i] = NA_REAL;
        *call->variance_unavailable = 1;
      }
      if(do_gerr)
        for(j = gradient_range.begin; j < gradient_range.end; ++j) {
          (*call->gradient_zero_out)[j*num_obs_eval_extern + ipe[i]] =
            ordinary_hc0_context.gradient_structural_zero[j][i];
          if(ordinary_hc0_context.gradient_unavailable[j][i]) {
            egerr[j][i] = NA_REAL;
            *call->variance_unavailable = 1;
          }
        }
    }
  }

  for(i=0;i<num_obs_eval_extern;i++)
    cm[ipe[i]] = ecm[i];

  if(do_merr)
    for(i=0;i<num_obs_eval_extern;i++)
      cmerr[ipe[i]] = ecmerr[i];


  if(do_grad){
    if(!(ordinary_hc0_active && train_is_eval && do_gerr))
      for(j=gradient_range.begin;j<gradient_range.end;j++)
        for(i=0;i<num_obs_eval_extern;i++)
          g[j*num_obs_eval_extern+ipe[i]]=eg[j][i];

    /* The canonical mixed beta route defines and returns discrete-contrast
       standard errors.  Preserve the historical non-beta bridge unchanged;
       its narrower copy contract is tracked as separate legacy debt. */
    if(do_gerr) {
      if(ordinary_hc0_active) {
        for(j=gradient_range.begin;j<gradient_range.end;j++)
          for(i=0;i<num_obs_eval_extern;i++)
            gerr[j*num_obs_eval_extern+ipe[i]] = egerr[j][i];
      } else {
        const int gradient_stderr_count =
          (kernel_route != NULL) ? num_var : num_reg_continuous_extern;

        for(j=0;j<gradient_stderr_count;j++)
          for(i=0;i<num_obs_eval_extern;i++)
            gerr[j*num_obs_eval_extern+ipe[i]]=egerr[j][i];
      }
    }
  }


  /* write the return values */


  xtra[0] = RS;
  xtra[1] = MSE;
  xtra[2] = MAE;
  xtra[3] = MAPE;
  xtra[4] = CORR;
  xtra[5] = SIGN;

  return R_NilValue;
}
