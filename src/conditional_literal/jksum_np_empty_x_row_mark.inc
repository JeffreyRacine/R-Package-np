static int np_empty_x_row_mark(const NPEmptyXRowContext *ctx, int row,
                              int contrast, NPRegressionLPEmptyRows *rows)
{
  if(rows == NULL) return 0;
  const int dims[3] = {ctx->nu, ctx->no, ctx->nc};
  const size_t nptr = (size_t)ctx->nu + ctx->no + 2*(size_t)ctx->nc;
  double **storage = nptr > 0 ? (double **)malloc(nptr*sizeof(double *)) : NULL;
  double *weights = (double *)malloc((size_t)ctx->n*sizeof(double));
  double **eval[3] = {NULL, NULL, NULL}, **bandwidth_eval = NULL;
  double alternate = 0.0, sum = 0.0;
  size_t offset = 0;
  int empty = 0;
  if((nptr > 0 && storage == NULL) || weights == NULL) goto cleanup;
  for(int kind = 0; kind < 3; ++kind) {
    if(dims[kind]) eval[kind] = storage + offset;
    for(int j = 0; j < dims[kind]; ++j)
      eval[kind][j] = ctx->eval[kind][j] + row;
    offset += dims[kind];
  }
  if(ctx->nc) {
    bandwidth_eval = storage + offset;
    for(int j = 0; j < ctx->nc; ++j)
      bandwidth_eval[j] = ctx->bandwidth[j] +
        (ctx->bw == BW_GEN_NN ? row : 0);
  }
  if(contrast >= 0) {
    const int cat = contrast - ctx->nc;
    if(cat < 0 || cat >= ctx->nu + ctx->no || ctx->categories[cat] <= 0)
      goto cleanup;
    if(cat < ctx->nu) {
      alternate = ctx->catvals[cat][0];
      eval[0][cat] = &alternate;
    } else {
      const int ordered = cat - ctx->nu;
      const double current = *eval[1][ordered];
      int index = 0;
      while(index < ctx->categories[cat] && ctx->catvals[cat][index] != current)
        ++index;
      if(index == ctx->categories[cat]) goto cleanup;
      alternate = ctx->catvals[cat][index == 0 ?
        (ctx->categories[cat] > 1 ? 1 : 0) : index-1];
      eval[1][ordered] = &alternate;
    }
  }
  const int status = kernel_weighted_sum_np_ctx_ex(
    ctx->kc, ctx->ku, ctx->ko, ctx->bw, ctx->n, 1,
    ctx->nu, ctx->no, ctx->nc, 0, 0, 1, 1, 0, 0, 0, 0, 0,
    ctx->op, OP_NOOP, 0, 0, NULL, 1, 0, 0,
    NP_TREE_FALSE, 0, NULL, NULL, NULL, NULL,
    ctx->train[0], ctx->train[1], ctx->train[2],
    eval[0], eval[1], eval[2], NULL, NULL, NULL,
    ctx->vsf, 1, ctx->bandwidth, bandwidth_eval,
    ctx->lambda, ctx->categories, ctx->catvals, NULL,
    &sum, NULL, weights, NULL, NULL, NULL, NULL, NULL,
#ifdef MPI2
    0,
#endif
    NULL,
    NULL);
  if(status == 0 && np_lp_complete_weights_are_zero(weights, ctx->n)) {
    np_regression_empty_row_mark(rows, row, contrast >= 0);
    if(contrast >= 0 && rows->component_flags != NULL)
      rows->component_flags[contrast] = 1;
    empty = 1;
  }
cleanup:
  free(weights);
  free(storage);
  return empty;
}
