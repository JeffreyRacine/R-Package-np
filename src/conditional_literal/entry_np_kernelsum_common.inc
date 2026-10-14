#line 24333 "np.c"
static void np_kernelsum_common(double * tuno, double * tord, double * tcon,
                                double * ty, double * weights,
                                double * euno, double * eord, double * econ,
                                double * bw,
                                double * mcv, double * padnum,
                                int * operator,
                                int * myopti, double * kpow,
                                double * weighted_sum,
                                double * weighted_sum_power2,
                                double * weighted_p_sum,
                                double * kernel_weights,
                                double * permutation_kernel_weights,
                                double * ckerlb, double * ckerub,
                                const int *eval_to_train, int fold_requested){

  int * ipt = NULL, * ipe = NULL;  // point permutation, see tree.c
      
  /* the ys are the weights */

  double * vector_scale_factor, * ksum, * ksum2 = NULL, * p_ksum = NULL, pad_num, * kw = NULL, * pkw = NULL;
  double **query_bandwidth = NULL;
  double *query_lambda = NULL;
  int i,j,k, num_var, num_obs_eval_alloc;
  int no_y, leave_one_out, train_is_eval, do_divide_bw;
  int do_divide_returned_weights, do_divide_returned_bw;
  int max_lev, no_weights, sum_element_length, return_kernel_weights;
  int p_operator, do_score, do_ocg, p_nvar = 0;
  int tree_outer_blas_requested = 0;
  int query_bandwidth_provided = 0;
  NPNNKernelFold fold = {0};
  int suppress_parallel_ksum = 0;

  int use_tree = 0;
  int allocated_X_train = 1, allocated_X_eval = 1;
  int allocated_Y = 1, allocated_W = 1;
  int ksum_is_output = 0, ksum2_is_output = 0, pksum_is_output = 0, kw_is_output = 0, pkw_is_output = 0;

  struct th_table * otabs = NULL;
  struct th_entry * ret = NULL;
  int ** matrix_ordered_indices = NULL;

  int ncol_Y, ncol_W;

  int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;

  int npks_err = 0;

  /* OCG needs an index in the retained ordered-category table. Ordinary
     numeric ordered evaluation may extrapolate and must remain unrestricted.
     Check before native allocations (and before any MPI work is entered). */
  if(myopti[KWS_POCGI] && myopti[KWS_NORDI] > 0){
    const int rows = myopti[KWS_TISEI] ? myopti[KWS_TNOBSI] : myopti[KWS_ENOBSI];
    const int levels = myopti[KWS_MLEVI];
    const double *eval = myopti[KWS_TISEI] ? tord : eord;
    if(levels <= 0 || rows < 0)
      error("npksum: invalid ordered contrast metadata");
    for(int coordinate = 0; coordinate < myopti[KWS_NORDI]; ++coordinate){
      const double *categories = mcv +
        (size_t)(myopti[KWS_NUNOI] + coordinate)*(size_t)levels;
      int count = 0;
      while(count < levels && categories[count] != *padnum) ++count;
      for(int row = 0; row < rows; ++row){
        const double value = eval[(size_t)coordinate*(size_t)rows + (size_t)row];
        int index = 0;
        while(index < count && categories[index] != value) ++index;
        if(index == count)
          error("npksum: compute.ocg requires ordered evaluation values in the bandwidth category table");
      }
    }
  }

  /* match integer options with their globals */

  num_reg_continuous_extern = myopti[KWS_NCONI];
  num_reg_unordered_extern = myopti[KWS_NUNOI];
  num_reg_ordered_extern = myopti[KWS_NORDI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[KWS_TNOBSI];
  num_obs_eval_extern = myopti[KWS_ENOBSI];

  KERNEL_reg_extern = myopti[KWS_CKRNEVI];
  KERNEL_reg_unordered_extern = myopti[KWS_UKRNEVI];
  KERNEL_reg_ordered_extern = myopti[KWS_OKRNEVI];
  vector_ckerlb_extern = ckerlb;
  vector_ckerub_extern = ckerub;
  int_cker_bound_extern = np_has_finite_cker_bounds(ckerlb, ckerub, num_reg_continuous_extern);

  int_LARGE_SF = myopti[KWS_LSFI];
  int_MINIMIZE_IO = myopti[KWS_MINIOI];
  BANDWIDTH_reg_extern = myopti[KWS_BWI];

  train_is_eval = myopti[KWS_TISEI];
  // no_y = myopti[KWS_NOYI];
  leave_one_out = myopti[KWS_LOOI];
  do_divide_bw = myopti[KWS_BDIVI];
  do_divide_returned_weights = myopti[KWS_BDIVWI];
  do_divide_returned_bw =
    do_divide_returned_weights ? do_divide_bw : 0;
  
  max_lev = myopti[KWS_MLEVI];
  pad_num = *padnum;

  /* the y and weight matrices will be contained in these variables */
  ncol_Y = myopti[KWS_YNCOLI];
  ncol_W = myopti[KWS_WNCOLI];

  int_TREE_X = myopti[KWS_DOTREEI];
  tree_outer_blas_requested =
    (myopti[KWS_TREEBLASI] == 1) ? 1 : 0;
  return_kernel_weights = myopti[KWS_RKWI];
  p_operator = myopti[KWS_POPI];
  do_score = myopti[KWS_PSCOREI];
  do_ocg = myopti[KWS_POCGI];
  suppress_parallel_ksum = myopti[KWS_SPARI];

  nconfac_extern = ncatfac_extern = 0.0;

  no_y = (ncol_Y == 0);
  no_weights = (ncol_W == 0);

  sum_element_length = (no_y ? 1 : ncol_Y)*(no_weights ? 1 : ncol_W);

  use_tree = (int_TREE_X == NP_TREE_TRUE);

#ifdef MPI2
  num_obs_eval_alloc =
    np_mpi_padded_count_or_error(num_obs_eval_extern, "np_kernelsum");
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif

  if(train_is_eval && (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_kernelsum): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    error("\n(np_kernelsum): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
  }

  /* allocate */

  if(use_tree){
    matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
    matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
    matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);
  } else {
    allocated_X_train = 0;
    matrix_X_unordered_train_extern = (num_reg_unordered_extern > 0) ?
      (double **)R_alloc((size_t)num_reg_unordered_extern, sizeof(double*)) : NULL;
    matrix_X_ordered_train_extern = (num_reg_ordered_extern > 0) ?
      (double **)R_alloc((size_t)num_reg_ordered_extern, sizeof(double*)) : NULL;
    matrix_X_continuous_train_extern = (num_reg_continuous_extern > 0) ?
      (double **)R_alloc((size_t)num_reg_continuous_extern, sizeof(double*)) : NULL;

    for(j = 0; j < num_reg_unordered_extern; j++)
      matrix_X_unordered_train_extern[j] = tuno + j*num_obs_train_extern;
    for(j = 0; j < num_reg_ordered_extern; j++)
      matrix_X_ordered_train_extern[j] = tord + j*num_obs_train_extern;
    for(j = 0; j < num_reg_continuous_extern; j++)
      matrix_X_continuous_train_extern[j] = tcon + j*num_obs_train_extern;
  }
  
  /* for the moment we will just allocate a vector of ones */
  /* vector_Y_extern = (no_y)?NULL:alloc_vecd(num_obs_train_extern); */

  if(use_tree){
    matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, ncol_Y);
    matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, ncol_W);
  } else {
    allocated_Y = 0;
    allocated_W = 0;
    matrix_Y_continuous_train_extern = (ncol_Y > 0) ?
      (double **)R_alloc((size_t)ncol_Y, sizeof(double*)) : NULL;
    matrix_Y_ordered_train_extern = (ncol_W > 0) ?
      (double **)R_alloc((size_t)ncol_W, sizeof(double*)) : NULL;

    for(j = 0; j < ncol_Y; j++)
      matrix_Y_continuous_train_extern[j] = ty + j*num_obs_train_extern;
    for(j = 0; j < ncol_W; j++)
      matrix_Y_ordered_train_extern[j] = weights + j*num_obs_train_extern;
  }

  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_categorical_vals_extern = alloc_matd(max_lev, num_reg_unordered_extern + num_reg_ordered_extern);

  vector_scale_factor = alloc_vecd(num_var + 1);
  if(!use_tree && (num_obs_eval_alloc == num_obs_eval_extern)){
    ksum = weighted_sum;
    ksum_is_output = 1;
  } else {
    if(np_int_product_overflows(num_obs_eval_alloc, sum_element_length))
      error("C_np_kernelsum: requested kernel-sum buffer is too large");
    ksum = alloc_vecd(num_obs_eval_alloc*sum_element_length);
  }

  if(weighted_sum_power2 != NULL){
    if(!use_tree && (num_obs_eval_alloc == num_obs_eval_extern)){
      ksum2 = weighted_sum_power2;
      ksum2_is_output = 1;
    } else {
      if(np_int_product_overflows(num_obs_eval_alloc, sum_element_length))
        error("C_np_kernelsum_power12: requested power-two kernel-sum buffer is too large");
      ksum2 = alloc_vecd(num_obs_eval_alloc*sum_element_length);
    }
  }

  if((p_operator != OP_NOOP) || do_score || do_ocg){
    p_nvar = ((p_operator != OP_NOOP) ? num_reg_continuous_extern : 0) + ((do_score || do_ocg) ? num_reg_unordered_extern + num_reg_ordered_extern : 0);
    if(!use_tree && (num_obs_eval_alloc == num_obs_eval_extern)){
      p_ksum = weighted_p_sum;
      pksum_is_output = 1;
    } else {
      if(np_int_product3_overflows(num_obs_eval_alloc, sum_element_length, p_nvar))
        error("C_np_kernelsum: requested derivative kernel-sum buffer is too large");
      p_ksum = alloc_vecd(num_obs_eval_alloc*sum_element_length*p_nvar);
    }
  }

  if(!train_is_eval){
    if(use_tree){
      matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
      matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
      matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);
    } else {
      allocated_X_eval = 0;
      matrix_X_unordered_eval_extern = (num_reg_unordered_extern > 0) ?
        (double **)R_alloc((size_t)num_reg_unordered_extern, sizeof(double*)) : NULL;
      matrix_X_ordered_eval_extern = (num_reg_ordered_extern > 0) ?
        (double **)R_alloc((size_t)num_reg_ordered_extern, sizeof(double*)) : NULL;
      matrix_X_continuous_eval_extern = (num_reg_continuous_extern > 0) ?
        (double **)R_alloc((size_t)num_reg_continuous_extern, sizeof(double*)) : NULL;

      for(j = 0; j < num_reg_unordered_extern; j++)
        matrix_X_unordered_eval_extern[j] = euno + j*num_obs_eval_extern;
      for(j = 0; j < num_reg_ordered_extern; j++)
        matrix_X_ordered_eval_extern[j] = eord + j*num_obs_eval_extern;
      for(j = 0; j < num_reg_continuous_extern; j++)
        matrix_X_continuous_eval_extern[j] = econ + j*num_obs_eval_extern;
    }
  } else {
    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;
  }

  /* train */

  if(use_tree){
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+i];

    for( j = 0; j < ncol_Y; j++ )
      for( i = 0; i < num_obs_train_extern; i++ )
        matrix_Y_continuous_train_extern[j][i] = ty[j*num_obs_train_extern+i];

    for( j = 0; j < ncol_W; j++ )
      for( i = 0; i < num_obs_train_extern; i++ )
        matrix_Y_ordered_train_extern[j][i] = weights[j*num_obs_train_extern+i];
  }

  if(use_tree && !train_is_eval){
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

  if(use_tree){
    ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
    if(!(ipt != NULL))
      error("!(ipt != NULL)");

    for(i = 0; i < num_obs_train_extern; i++){
      ipt[i] = i;
    }

    if(!train_is_eval) {
      ipe = (int *)malloc(num_obs_eval_extern*sizeof(int));
      if(!(ipe != NULL))
        error("!(ipe != NULL)");

      for(i = 0; i < num_obs_eval_extern; i++){
        ipe[i] = i;
      }
    } else {
      ipe = ipt;
    }
  } else {
    ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
    if(!(ipt != NULL))
      error("!(ipt != NULL)");

    for(i = 0; i < num_obs_train_extern; i++){
      ipt[i] = i;
    }

    if(!train_is_eval){
      ipe = (int *)malloc(num_obs_eval_extern*sizeof(int));
      if(!(ipe != NULL))
        error("!(ipe != NULL)");

      for(i = 0; i < num_obs_eval_extern; i++){
        ipe[i] = i;
      }
    } else {
      ipe = ipt;
    }
  }

  // attempt tree build, if enabled 
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    if((BANDWIDTH_reg_extern != BW_ADAP_NN) || ((BANDWIDTH_reg_extern == BW_ADAP_NN) && train_is_eval)){
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipt, &kdt_extern_X);

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+ipt[i]];
    
    
      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+ipt[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+ipt[i]];

      for( j = 0; j < ncol_Y; j++ )
        for( i = 0; i < num_obs_train_extern; i++ )
          matrix_Y_continuous_train_extern[j][i] = ty[j*num_obs_train_extern+ipt[i]];

      for( j = 0; j < ncol_W; j++ )
        for( i = 0; i < num_obs_train_extern; i++ )
          matrix_Y_ordered_train_extern[j][i] = weights[j*num_obs_train_extern+ipt[i]];

    } else {
      build_kdtree(matrix_X_continuous_eval_extern, num_obs_eval_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipe, &kdt_extern_X);


      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+ipe[i]];

    }
  }



  /* bandwidths */
  for( i=0; i<num_var; i++ )
    vector_scale_factor[i+1] = bw[i];
  
  /* fix up categories */

  for(j=0; j<num_reg_unordered_extern; j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;
  }

  if(do_ocg && (num_reg_ordered_extern > 0)){
    if(num_obs_eval_extern < 0 ||
       (size_t)num_obs_eval_extern > SIZE_MAX / (size_t)num_reg_ordered_extern ||
       ((size_t)num_reg_ordered_extern * (size_t)num_obs_eval_extern) > SIZE_MAX / sizeof(int))
      error("ordered index buffer size overflow");

    size_t ordered_index_count = (size_t)num_reg_ordered_extern * (size_t)num_obs_eval_extern;
    otabs = (struct th_table *)malloc(num_reg_ordered_extern*sizeof(struct th_table));
    matrix_ordered_indices = (int **)malloc(num_reg_ordered_extern*sizeof(int *));
    int * tc = (int *)malloc(ordered_index_count*sizeof(int));
    if(tc == NULL)
      error("failed to allocate ordered index buffer");
    for(i = 0; i < num_reg_ordered_extern; i++)
      matrix_ordered_indices[i] = tc + i*num_obs_eval_extern;
  }

  for(j=num_reg_unordered_extern, k = 0; j < (num_reg_unordered_extern+num_reg_ordered_extern); j++, k++){
    i = 0;

    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;

    if(do_ocg){
      if(thcreate_r((size_t)ceil(1.6*num_categories_extern[j]), otabs + k) == TH_ERROR)
        error("hash table creation failed");

      for(i = 0; i < num_categories_extern[j]; i++){
        struct th_entry centry;
        centry.key.dkey = matrix_categorical_vals_extern[j][i];
        centry.data = i;

        if(thsearch_r(&centry, TH_ENTER, &ret, otabs+k) == TH_FAILURE)
          error("insertion into hash table failed");
      }
      
      // now do lookups
      struct th_entry te;
      te.key.dkey = pad_num;
      te.data = -1;

      ret = &te;
      for(i = 0; i < num_obs_eval_extern; i++){
        if(i == 0 || ret->key.dkey != matrix_X_ordered_eval_extern[k][i]){
          te.key.dkey = matrix_X_ordered_eval_extern[k][i];
          if(thsearch_r(&te, TH_SEARCH, &ret, otabs+k) == TH_FAILURE)
            error("hash table lookup failed (which should be impossible)");

        } 

        matrix_ordered_indices[k][i] = ret->data;
      }
    }
  }

  if(return_kernel_weights){
    if(!use_tree && (BANDWIDTH_reg_extern != BW_ADAP_NN) && (num_obs_eval_alloc == num_obs_eval_extern)){
      kw = kernel_weights;
      kw_is_output = 1;
    } else {
      if(np_int_product_overflows(num_obs_train_extern, num_obs_eval_extern))
        error("C_np_kernelsum: requested kernel-weight buffer is too large");
      kw = alloc_vecd(num_obs_train_extern*num_obs_eval_extern);
    }
  }

  if(return_kernel_weights && (p_nvar > 0)){
    if(np_int_product3_overflows(num_obs_train_extern, num_obs_eval_extern, p_nvar))
      error("C_np_kernelsum: requested derivative kernel-weight buffer is too large");
    pkw = alloc_vecd(num_obs_train_extern*num_obs_eval_extern*p_nvar);
  }
  
  kernel_c = (int *)malloc(sizeof(int)*num_reg_continuous_extern);

  for(i = 0; i < num_reg_continuous_extern; i++)
    kernel_c[i] = KERNEL_reg_extern;

  kernel_u = (int *)malloc(sizeof(int)*num_reg_unordered_extern);

  for(i = 0; i < num_reg_unordered_extern; i++)
    kernel_u[i] = KERNEL_reg_unordered_extern;

  kernel_o = (int *)malloc(sizeof(int)*num_reg_ordered_extern);

  for(i = 0; i < num_reg_ordered_extern; i++)
    kernel_o[i] = KERNEL_reg_ordered_extern;

  int *bpso = NULL;
  {
    const int bpso_len = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
    if(bpso_len > 0){
      const int do_perm = (p_operator != OP_NOOP);
      const int doscoreocg = (do_score || do_ocg);

      bpso = (int *)R_alloc((size_t)bpso_len, sizeof(int));
      for(i = 0; i < num_reg_continuous_extern; i++)
        bpso[i] = do_perm;
      for(i = num_reg_continuous_extern; i < bpso_len; i++)
        bpso[i] = doscoreocg;
    }
  }

  /*
   * The kernel-sum ingress owns whether evaluation rows are the training
   * occurrences themselves. Preserve that identity when generalized-NN
   * radii are constructed. The existing O(n p) bandwidth workspace is
   * prepared once and consumed by the unchanged serial/MPI weighted-sum
   * owner; rank-local calls honor the caller's suppress-parallel contract.
   */
  const int *query_map = eval_to_train;
  if(eval_to_train != NULL && use_tree) {
    int *inverse_train = (int *)R_alloc((size_t)num_obs_train_extern, sizeof(int));
    int *permuted_map = (int *)R_alloc((size_t)num_obs_eval_extern, sizeof(int));
    for(i = 0; i < num_obs_train_extern; ++i) inverse_train[ipt[i]] = i;
    for(i = 0; i < num_obs_eval_extern; ++i)
      permuted_map[i] = inverse_train[eval_to_train[ipe[i]]];
    query_map = permuted_map;
  }
  if(fold_requested) {
    query_lambda = (double *)R_alloc(
      (size_t)MAX(1, num_reg_unordered_extern + num_reg_ordered_extern), sizeof(double));
    npks_err = np_kernelsum_fold_prepare(&fold, BANDWIDTH_reg_extern,
      num_obs_train_extern, num_obs_eval_extern, num_reg_continuous_extern,
      num_reg_unordered_extern, num_reg_ordered_extern,
      matrix_X_continuous_train_extern, matrix_X_continuous_eval_extern,
      query_map, &vector_scale_factor[1], query_lambda, 0, suppress_parallel_ksum);
    if(npks_err != 0) goto cleanup_kernelsum_fold;
    query_bandwidth = fold.primary;
    query_bandwidth_provided = 1;
  } else if(BANDWIDTH_reg_extern == BW_GEN_NN && (train_is_eval || eval_to_train != NULL) &&
     num_reg_continuous_extern > 0){
    const NPNNGeometryContext query_geometry = {
      .mode = query_map == NULL ? NP_NN_QUERY_TRAINING_IDENTITY :
                                 NP_NN_QUERY_TRAINING_MAP,
      .eval_to_train = query_map,
      .adaptive_successor = NULL
    };
    NPNNGeometryStatus query_status = NP_NN_GEOMETRY_OK;

    query_bandwidth = alloc_tmatd(num_obs_eval_extern,
                                  num_reg_continuous_extern);
    query_lambda = alloc_vecd(MAX(
      1, num_reg_unordered_extern + num_reg_ordered_extern));
    if(query_bandwidth == NULL || query_lambda == NULL){
      free_tmat(query_bandwidth);
      free(query_lambda);
      error("C_np_kernelsum: generalized nearest-neighbor training geometry allocation failed");
    }
    if(kernel_bandwidth_mean_ctx(
         KERNEL_reg_extern,
         BANDWIDTH_reg_extern,
         num_obs_train_extern,
         num_obs_eval_extern,
         0,
         0,
         0,
         num_reg_continuous_extern,
         num_reg_unordered_extern,
         num_reg_ordered_extern,
         suppress_parallel_ksum,
         &vector_scale_factor[1],
         NULL,
         NULL,
         matrix_X_continuous_train_extern,
         matrix_X_continuous_eval_extern,
         NULL,
         query_bandwidth,
         query_lambda,
         &query_geometry,
         NULL,
         &query_status) != 0){
      free_tmat(query_bandwidth);
      free(query_lambda);
      if(query_status == NP_NN_GEOMETRY_ZERO_RADIUS) {
        const NPNNZeroRadiusInfo info = np_nn_zero_radius_info(
          BANDWIDTH_reg_extern, num_obs_train_extern, num_obs_eval_extern,
          num_reg_continuous_extern, matrix_X_continuous_train_extern,
          matrix_X_continuous_eval_extern, &vector_scale_factor[1], &query_geometry, 0);
        np_nn_zero_radius_error(&info);
      }
      error("C_np_kernelsum: generalized nearest-neighbor training geometry failed with status %d",
            (int)query_status);
    }
    query_bandwidth_provided = 1;
  }
  
  
  if(weighted_sum_power2 != NULL){
    npks_err = kernel_weighted_sum_np_power12(kernel_c,
                                      kernel_u,
                                      kernel_o,
                                      BANDWIDTH_reg_extern,
                                      num_obs_train_extern,
                                      num_obs_eval_extern,
                                      num_reg_unordered_extern,
                                      num_reg_ordered_extern,
                                      num_reg_continuous_extern,
                                      leave_one_out,
                                      0,
                                      do_divide_bw,
                                      do_divide_returned_bw,
                                      0, //not symmetric
                                      0, //disable 'twisting'
                                      0, // do not drop train
                                      0, // do not drop train
                                      operator,
                                      p_operator,
                                      do_score,
                                      do_ocg, // no ocg (for now)
                                      bpso,
                                      suppress_parallel_ksum,
                                      ncol_Y,
                                      ncol_W,
                                      int_TREE_X,
                                      0,
                                      kdt_extern_X,
                                      NULL, NULL, NULL,
                                      matrix_X_unordered_train_extern,
                                      matrix_X_ordered_train_extern,
                                      matrix_X_continuous_train_extern,
                                      matrix_X_unordered_eval_extern,
                                      matrix_X_ordered_eval_extern,
                                      matrix_X_continuous_eval_extern,
                                      /* ys matrix */
                                      matrix_Y_continuous_train_extern,
                                      /* weights matrix */
                                      matrix_Y_ordered_train_extern,
                                      NULL,
                                      &vector_scale_factor[1],
                                      query_bandwidth_provided,
                                      query_bandwidth,
                                      query_bandwidth,
                                      query_lambda,
                                      num_categories_extern,
                                      matrix_categorical_vals_extern,
                                      matrix_ordered_indices,
                                      ksum,
                                      ksum2,
                                      p_ksum,
                                      kw,
                                      pkw);
  } else {
    int_TREE_OUTER_BLAS = tree_outer_blas_requested;
    npks_err = kernel_weighted_sum_np_fold_route(kernel_c,
                                      kernel_u,
                                      kernel_o,
                                      BANDWIDTH_reg_extern,
                                      num_obs_train_extern,
                                      num_obs_eval_extern,
                                      num_reg_unordered_extern,
                                      num_reg_ordered_extern,
                                      num_reg_continuous_extern,
                                      leave_one_out,
                                      0,
                                      (int)(*kpow),
                                      do_divide_bw,
                                      do_divide_returned_bw,
                                      0, //not symmetric
                                      0, //disable 'twisting'
                                      0, // do not drop train
                                      0, // do not drop train
                                      operator,
                                      p_operator,
                                      do_score,
                                      do_ocg,
                                      bpso,
                                      suppress_parallel_ksum,
                                      ncol_Y,
                                      ncol_W,
                                      int_TREE_X,
                                      0,
                                      kdt_extern_X,
                                      NULL, NULL, NULL,
                                      matrix_X_unordered_train_extern,
                                      matrix_X_ordered_train_extern,
                                      matrix_X_continuous_train_extern,
                                      matrix_X_unordered_eval_extern,
                                      matrix_X_ordered_eval_extern,
                                      matrix_X_continuous_eval_extern,
                                      matrix_Y_continuous_train_extern,
                                      matrix_Y_ordered_train_extern,
                                      NULL,
                                      &vector_scale_factor[1],
                                      query_bandwidth_provided,
                                      query_bandwidth,
                                      query_bandwidth,
                                      query_lambda,
                                      num_categories_extern,
                                      matrix_categorical_vals_extern,
                                      matrix_ordered_indices,
                                      ksum,
                                      p_ksum,
                                      kw,
                                      pkw, 0, NULL, NULL,
                                      fold_requested ? &fold : NULL);
  }
  if(query_bandwidth != NULL && !fold_requested)
    free_tmat(query_bandwidth);
  if(!fold_requested) free(query_lambda);

  if(npks_err != 0){
    if(fold_requested) goto cleanup_kernelsum_fold;
    if(npks_err == KWSNP_ERR_ZERO_NN_RADIUS) {
      const NPNNZeroRadiusInfo info = np_nn_zero_radius_info(
        BANDWIDTH_reg_extern, num_obs_train_extern, num_obs_eval_extern,
        num_reg_continuous_extern, matrix_X_continuous_train_extern,
        matrix_X_continuous_eval_extern, &vector_scale_factor[1], NULL, 0);
      np_nn_zero_radius_error(&info);
    }
    error("kernel_weighted_sum_np failed with code %d", npks_err);
  }


  if(!npks_err){
    if(use_tree || !ksum_is_output){
      for(j = 0; j < num_obs_eval_extern; j++)
        for(i = 0; i < sum_element_length; i++)
          weighted_sum[ipe[j]*sum_element_length + i] = ksum[j*sum_element_length+i];
    }

    if((weighted_sum_power2 != NULL) && (use_tree || !ksum2_is_output)){
      for(j = 0; j < num_obs_eval_extern; j++)
        for(i = 0; i < sum_element_length; i++)
          weighted_sum_power2[ipe[j]*sum_element_length + i] = ksum2[j*sum_element_length+i];
    }

    if(return_kernel_weights){
      if(use_tree || (BANDWIDTH_reg_extern == BW_ADAP_NN) || !kw_is_output){
        if(BANDWIDTH_reg_extern != BW_ADAP_NN){ // adaptive weights are currently returned transposed...
          for(j = 0; j < num_obs_eval_extern; j++)
            for(i = 0; i < num_obs_train_extern; i++)
              kernel_weights[ipe[j]*num_obs_train_extern + ipt[i]] = kw[j*num_obs_train_extern + i];
        } else {
          for(j = 0; j < num_obs_train_extern; j++)
            for(i = 0; i < num_obs_eval_extern; i++)
              kernel_weights[ipe[i]*num_obs_train_extern + ipt[j]] = kw[j*num_obs_eval_extern + i];      
        }
      }
    }

    if(return_kernel_weights && (p_nvar > 0) && (pkw != NULL)){
      if(BANDWIDTH_reg_extern != BW_ADAP_NN){
        for(k = 0; k < p_nvar; k++){
          const int koff = k*num_obs_train_extern*num_obs_eval_extern;
          for(j = 0; j < num_obs_eval_extern; j++)
            for(i = 0; i < num_obs_train_extern; i++)
              permutation_kernel_weights[koff + ipe[j]*num_obs_train_extern + ipt[i]] =
                pkw[koff + j*num_obs_train_extern + i];
        }
      } else {
        for(k = 0; k < p_nvar; k++){
          const int koff = k*num_obs_train_extern*num_obs_eval_extern;
          for(j = 0; j < num_obs_train_extern; j++)
            for(i = 0; i < num_obs_eval_extern; i++)
              permutation_kernel_weights[koff + ipe[i]*num_obs_train_extern + ipt[j]] =
                pkw[koff + j*num_obs_eval_extern + i];
        }
      }
    }

    if(p_nvar > 0){
      if(use_tree || !pksum_is_output){
        for(k = 0; k < p_nvar; k++){
          const int kidx = k*num_obs_eval_extern*sum_element_length;
          for(j = 0; j < num_obs_eval_extern; j++)
            for(i = 0; i < sum_element_length; i++)
              weighted_p_sum[kidx + ipe[j]*sum_element_length + i] = p_ksum[kidx + j*sum_element_length + i];
        }
      }
    }

  }
cleanup_kernelsum_fold:
  /* Fold-geometry/accumulator errors release the existing ingress workspace
     before signalling. The new geometry itself is call-local R storage. */

  if(allocated_X_train){
    free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  }
  
  if(!train_is_eval && allocated_X_eval){
    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);
  if(allocated_Y)
    free_mat(matrix_Y_continuous_train_extern, ncol_Y);
  if(allocated_W)
    free_mat(matrix_Y_ordered_train_extern, ncol_W);

  safe_free(num_categories_extern);
  safe_free(vector_scale_factor);
  if(!ksum_is_output)
    safe_free(ksum);
  if((ksum2 != NULL) && !ksum2_is_output)
    safe_free(ksum2);

  if(!kw_is_output)
    safe_free(kw);
  if(!pkw_is_output)
    safe_free(pkw);

  safe_free(ipt);

  free(kernel_c);
  free(kernel_u);
  free(kernel_o);

  if(!train_is_eval)
    safe_free(ipe);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  if(p_nvar > 0){
    if(!pksum_is_output)
      safe_free(p_ksum);
  }

  if(do_ocg && (num_reg_ordered_extern > 0)){

    for(k = 0; k < num_reg_ordered_extern; k++){
      thdestroy_r(otabs+k);
    }
    free(otabs);
    free(matrix_ordered_indices[0]);
    free(matrix_ordered_indices);
  }

  int_cker_bound_extern = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;
  np_reset_y_side_extern();
  np_clear_estimator_extern_aliases();
  num_obs_train_extern = 0;
  num_obs_eval_extern = 0;

  if(fold_requested && npks_err != 0)
    error("raw delete-one kernel geometry or accumulation failed with code %d", npks_err);
  return;
}
