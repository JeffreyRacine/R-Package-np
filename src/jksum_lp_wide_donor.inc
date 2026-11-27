/* Shared donor traversal. The standard row object always instantiates the
 * original update. The qualified Apple wide-row object instantiates tiling. */
static void np_lp_accumulate_dense_row_generic(
    const NPLPDenseRowContext *ctx)
{
  int i, a, b;
  const int nterms = ctx->nterms;

  for(i = 0; i < ctx->nsub; i++){
    const int orig_ii = ctx->row_j + 1 + i;
    const int ii = ctx->use_tree ? ctx->tree_lookup[orig_ii] : orig_ii;
    const int widx = ctx->use_tree ? ii : i;
    const double weight = ctx->weights[widx];
    const double reverse_weight = np_lp_reverse_pair_weight(weight,
      ctx->rly_log_normalizer,ii,ctx->eval_idx);

    if(weight == 0.0)
      continue;

    if(ctx->track_lowsupport){
      np_lp_dense_support_add(ctx->row_j, nterms, ctx->support_count);
      np_lp_dense_support_add(orig_ii, nterms, ctx->support_count);
    }

    {
      const double yi = ctx->response[ii];
      double * const fixed_moments =
        ctx->moments + (size_t)ctx->row_j*(size_t)nterms*(size_t)nterms;
      double * const moving_moments =
        ctx->moments + (size_t)orig_ii*(size_t)nterms*(size_t)nterms;
      double * const fixed_rhs =
        ctx->rhs + (size_t)ctx->row_j*(size_t)nterms;
      double * const moving_rhs =
        ctx->rhs + (size_t)orig_ii*(size_t)nterms;

#if NP_LP_COMPILE_TILED_DONOR
      /* The basis is stored by term. Gather contiguous donor tiles before
         updating full rows. Each entry retains its original donor order. */
      for(a = 0; a < nterms; ++a){
        const double weighted_bia = weight*ctx->basis[a][ii];
        fixed_rhs[a] += weighted_bia*yi;
        moving_rhs[a] += reverse_weight*ctx->eval_ybasis[a];
      }
      for(int first = 0; first < nterms; ){
        double donor[32];
        const int count = nterms - first < 32 ? nterms - first : 32;
        for(b = 0; b < count; ++b)
          donor[b] = ctx->basis[first+b][ii];
        for(a = 0; a < nterms; ++a){
          const double weighted_bia = weight*ctx->basis[a][ii];
          double * const row = fixed_moments + (size_t)a*nterms + first;
          for(b = 0; b < count; ++b)
            row[b] += weighted_bia*donor[b];
        }
        first += count;
      }
      for(a = 0; a < nterms; ++a){
        const size_t offset = (size_t)a*nterms;
        for(b = 0; b < nterms; ++b)
          moving_moments[offset+b] += reverse_weight*ctx->eval_outer[offset+b];
      }
#else
      for(a = 0; a < nterms; a++){
        const double bia = ctx->basis[a][ii];
        const double weighted_bia = weight*bia;
        fixed_rhs[a] += weighted_bia*yi;
        moving_rhs[a] += reverse_weight*ctx->eval_ybasis[a];
        for(b = 0; b < nterms; b++){
          fixed_moments[a*nterms+b] +=
            weighted_bia*ctx->basis[b][ii];
          moving_moments[a*nterms+b] +=
            reverse_weight*ctx->eval_outer[a*nterms+b];
        }
      }
#endif
    }
  }
}
