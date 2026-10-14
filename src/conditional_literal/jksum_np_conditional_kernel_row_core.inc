static int np_conditional_kernel_row_core(const int *kernel_c,
                                                 const int *kernel_u,
                                                 const int *kernel_o,
                                                 const int *operator,
                                                 const int BANDWIDTH_den,
                                                 const int num_train,
                                                 const int num_unordered,
                                                 const int num_ordered,
                                                 const int num_continuous,
                                                 double **matrix_train_unordered,
                                                 double **matrix_train_ordered,
                                                 double **matrix_train_continuous,
                                                 double **matrix_eval_unordered_one,
                                                 double **matrix_eval_ordered_one,
                                                 double **matrix_eval_continuous_one,
                                                 double *vector_scale_factor,
                                                 int bandwidth_provided,
                                                 double **matrix_bandwidth_train,
                                                 double **matrix_bandwidth_eval_one,
                                                 double *lambda,
                                                 int *num_categories,
                                                 double **matrix_categorical_vals,
                                                 const int int_tree,
                                                 KDT *kdt,
                                                 const int bandwidth_divide_weights,
                                                 double *kw_out,
                                                 double *mean_out){
  /* ANN swaps training/evaluation roles in the kernel traversal. This adapter
     supplies exactly one evaluation point, so a training-indexed tree cannot
     describe its vector operand. Retain tree traversal with the matching
     one-point geometry. R scratch is bounded per call and unwinds on errors. */
  KDN eval_node = {0};
  KDT eval_tree = {0};
  const void *scratch = NULL;
  if(BANDWIDTH_den == BW_ADAP_NN && int_tree == NP_TREE_TRUE &&
     num_continuous > 0){
    scratch = vmaxget();
    double *bounds = (double *)R_alloc((size_t)num_continuous, 2*sizeof(double));
    for(int coordinate = 0; coordinate < num_continuous; ++coordinate)
      bounds[2*coordinate] = bounds[2*coordinate+1] =
        matrix_eval_continuous_one[coordinate][0];
    eval_node.bb = bounds;
    eval_node.childl = eval_node.childu = KD_NOCHILD;
    eval_node.nlev = 1;
    eval_node.istart = 0;
    eval_tree.kdn = &eval_node;
    eval_tree.bb = bounds;
    eval_tree.ndim = num_continuous;
    eval_tree.nbucket = eval_tree.numnode = 1;
    kdt = &eval_tree;
  }
  const int status = kernel_weighted_sum_np_ctx((int *)kernel_c,
                                    (int *)kernel_u,
                                    (int *)kernel_o,
                                    BANDWIDTH_den,
                                    num_train,
                                    1,
                                    num_unordered,
                                    num_ordered,
                                    num_continuous,
                                    0,
                                    0,
                                    1,
                                    1,
                                    bandwidth_divide_weights,
                                    0,
                                    0,
                                    0,
                                    0,
                                    operator,
                                    OP_NOOP,
                                    0,
                                    0,
                                    NULL,
                                    1,
                                    0,
                                    0,
                                    int_tree,
                                    0,
                                    kdt,
                                    NULL,
                                    NULL,
                                    NULL,
                                    matrix_train_unordered,
                                    matrix_train_ordered,
                                    matrix_train_continuous,
                                    matrix_eval_unordered_one,
                                    matrix_eval_ordered_one,
                                    matrix_eval_continuous_one,
                                    NULL,
                                    NULL,
                                    NULL,
                                    vector_scale_factor,
                                    bandwidth_provided,
                                    matrix_bandwidth_train,
                                    matrix_bandwidth_eval_one,
                                    lambda,
                                    num_categories,
                                    matrix_categorical_vals,
                                    NULL,
                                    mean_out,
                                    NULL,
                                    kw_out,
                                    NULL,
                                    NULL);
  if(scratch != NULL) vmaxset(scratch);
  return status;
}
