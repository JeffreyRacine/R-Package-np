/* Fixed-row-resident accumulation for dense raw-basis regression CV. */

#include <stddef.h>

#include "jksum_lp_row.h"

static inline void np_lp_dense_support_add(const int row,
                                            const int orig_idx,
                                            const int data_idx,
                                            const double weight,
                                            const int nterms,
                                            int *support_count,
                                            int *support_orig,
                                            int *support_data,
                                            double *support_weight)
{
  const int count = support_count[row];
  const size_t off = (size_t)row*(size_t)nterms;

  if(count > nterms)
    return;

  if(count == nterms){
    support_count[row] = nterms + 1;
    return;
  }

  if(count < nterms){
    support_orig[off + (size_t)count] = orig_idx;
    support_data[off + (size_t)count] = data_idx;
    support_weight[off + (size_t)count] = weight;
  }

  support_count[row] = count + 1;
}

/*
  A one-column LP basis is the degree-zero constant basis.  Keep its resident
  row scalar and implicit: no basis loads or general outer-product algebra.
  The wider generated siblings below retain their incumbent transcript.
*/
static void np_lp_accumulate_dense_row_1(
    const NPLPDenseRowContext *ctx)
{
  int i;
  double * const stored_moment = ctx->moments + (size_t)ctx->row_j;
  double * const stored_rhs = ctx->rhs + (size_t)ctx->row_j;
  const double eval_response = ctx->response[ctx->eval_idx];
  double fixed_moment = *stored_moment;
  double fixed_rhs = *stored_rhs;

  if(!ctx->use_tree){
    for(i = 0; i < ctx->nsub; i++){
      const int orig_ii = ctx->row_j + 1 + i;
      const double weight = ctx->weights[i];

      if(weight == 0.0)
        continue;

      if(ctx->track_lowsupport){
        np_lp_dense_support_add(ctx->row_j, orig_ii, orig_ii, weight, 1,
                                ctx->support_count, ctx->support_orig,
                                ctx->support_data, ctx->support_weight);
        np_lp_dense_support_add(orig_ii, ctx->row_j, ctx->eval_idx, weight, 1,
                                ctx->support_count, ctx->support_orig,
                                ctx->support_data, ctx->support_weight);
      }

      fixed_rhs += weight*ctx->response[orig_ii];
      ctx->rhs[orig_ii] += weight*eval_response;
      fixed_moment += weight;
      ctx->moments[orig_ii] += weight;
    }
  } else {
    for(i = 0; i < ctx->nsub; i++){
      const int orig_ii = ctx->row_j + 1 + i;
      const int ii = ctx->tree_lookup[orig_ii];
      const double weight = ctx->weights[ii];

      if(weight == 0.0)
        continue;

      if(ctx->track_lowsupport){
        np_lp_dense_support_add(ctx->row_j, orig_ii, ii, weight, 1,
                                ctx->support_count, ctx->support_orig,
                                ctx->support_data, ctx->support_weight);
        np_lp_dense_support_add(orig_ii, ctx->row_j, ctx->eval_idx, weight, 1,
                                ctx->support_count, ctx->support_orig,
                                ctx->support_data, ctx->support_weight);
      }

      fixed_rhs += weight*ctx->response[ii];
      ctx->rhs[orig_ii] += weight*eval_response;
      fixed_moment += weight;
      ctx->moments[orig_ii] += weight;
    }
  }

  *stored_rhs = fixed_rhs;
  *stored_moment = fixed_moment;
}

#define NP_LP_DEFINE_RESIDENT_WIDTH(WIDTH)                                  \
static void np_lp_accumulate_dense_row_##WIDTH(                             \
    const NPLPDenseRowContext *ctx)                                         \
{                                                                            \
  int i, a, b;                                                               \
  double fixed_moments[(WIDTH)*(WIDTH)];                                     \
  double fixed_rhs[(WIDTH)];                                                 \
  double * const stored_moments =                                           \
    ctx->moments + (size_t)ctx->row_j*(size_t)(WIDTH)*(size_t)(WIDTH);       \
  double * const stored_rhs =                                               \
    ctx->rhs + (size_t)ctx->row_j*(size_t)(WIDTH);                           \
                                                                             \
  for(a = 0; a < (WIDTH); a++){                                              \
    fixed_rhs[a] = stored_rhs[a];                                            \
    for(b = 0; b < (WIDTH); b++)                                             \
      fixed_moments[a*(WIDTH)+b] = stored_moments[a*(WIDTH)+b];              \
  }                                                                          \
                                                                             \
  for(i = 0; i < ctx->nsub; i++){                                            \
    const int orig_ii = ctx->row_j + 1 + i;                                 \
    const int ii = ctx->use_tree ? ctx->tree_lookup[orig_ii] : orig_ii;      \
    const int widx = ctx->use_tree ? ii : i;                                \
    const double weight = ctx->weights[widx];                               \
    double *moving_moments;                                                  \
    double *moving_rhs;                                                      \
    double yi;                                                               \
                                                                             \
    if(weight == 0.0)                                                        \
      continue;                                                              \
                                                                             \
    if(ctx->track_lowsupport){                                               \
      np_lp_dense_support_add(ctx->row_j, orig_ii, ii, weight, (WIDTH),      \
                              ctx->support_count, ctx->support_orig,          \
                              ctx->support_data, ctx->support_weight);        \
      np_lp_dense_support_add(orig_ii, ctx->row_j, ctx->eval_idx, weight,    \
                              (WIDTH), ctx->support_count,                    \
                              ctx->support_orig, ctx->support_data,           \
                              ctx->support_weight);                           \
    }                                                                        \
                                                                             \
    yi = ctx->response[ii];                                                   \
    moving_moments =                                                        \
      ctx->moments + (size_t)orig_ii*(size_t)(WIDTH)*(size_t)(WIDTH);        \
    moving_rhs = ctx->rhs + (size_t)orig_ii*(size_t)(WIDTH);                 \
    for(a = 0; a < (WIDTH); a++){                                            \
      const double bia = ctx->basis[a][ii];                                 \
      const double weighted_bia = weight*bia;                               \
      fixed_rhs[a] += weighted_bia*yi;                                      \
      if((WIDTH) >= 3){                                                      \
        moving_rhs[a] += weight*ctx->eval_ybasis[a];                         \
        for(b = 0; b < (WIDTH); b++){                                        \
          fixed_moments[a*(WIDTH)+b] +=                                     \
            weighted_bia*ctx->basis[b][ii];                                 \
          moving_moments[a*(WIDTH)+b] +=                                    \
            weight*ctx->eval_outer[a*(WIDTH)+b];                            \
        }                                                                    \
      } else {                                                               \
        const double bja = ctx->basis[a][ctx->eval_idx];                    \
        const double weighted_bja = weight*bja;                             \
        moving_rhs[a] += weighted_bja*ctx->response[ctx->eval_idx];          \
        for(b = 0; b < (WIDTH); b++){                                        \
          fixed_moments[a*(WIDTH)+b] +=                                     \
            weighted_bia*ctx->basis[b][ii];                                 \
          moving_moments[a*(WIDTH)+b] +=                                    \
            weighted_bja*ctx->basis[b][ctx->eval_idx];                      \
        }                                                                    \
      }                                                                      \
    }                                                                        \
  }                                                                          \
                                                                             \
  for(a = 0; a < (WIDTH); a++){                                              \
    stored_rhs[a] = fixed_rhs[a];                                            \
    for(b = 0; b < (WIDTH); b++)                                             \
      stored_moments[a*(WIDTH)+b] = fixed_moments[a*(WIDTH)+b];              \
  }                                                                          \
}

NP_LP_DEFINE_RESIDENT_WIDTH(2)
NP_LP_DEFINE_RESIDENT_WIDTH(4)
NP_LP_DEFINE_RESIDENT_WIDTH(5)
NP_LP_DEFINE_RESIDENT_WIDTH(6)
NP_LP_DEFINE_RESIDENT_WIDTH(7)
NP_LP_DEFINE_RESIDENT_WIDTH(8)
NP_LP_DEFINE_RESIDENT_WIDTH(9)
NP_LP_DEFINE_RESIDENT_WIDTH(10)
NP_LP_DEFINE_RESIDENT_WIDTH(11)
NP_LP_DEFINE_RESIDENT_WIDTH(12)
NP_LP_DEFINE_RESIDENT_WIDTH(13)
NP_LP_DEFINE_RESIDENT_WIDTH(14)
NP_LP_DEFINE_RESIDENT_WIDTH(15)
NP_LP_DEFINE_RESIDENT_WIDTH(16)

#undef NP_LP_DEFINE_RESIDENT_WIDTH

void np_lp_accumulate_dense_resident_row3(
    const int row_j,
    const int nsub,
    const int use_tree,
    const int eval_idx,
    const int track_lowsupport,
    const int *tree_lookup,
    const double *weights,
    double * const *basis,
    const double *response,
    double *moments,
    double *rhs,
    const double *eval_ybasis,
    const double *eval_outer,
    int *support_count,
    int *support_orig,
    int *support_data,
    double *support_weight)
{
  /*
   * Keep the unique upper triangle resident in the non-MPI pairwise route.
   * The caller mirrors it after every unordered pair has been accumulated.
   */
  int i;
  double * const sj = moments + (size_t)row_j*9;
  double * const tj = rhs + (size_t)row_j*3;
  double sj0 = sj[0], sj1 = sj[1], sj2 = sj[2];
  double sj4 = sj[4], sj5 = sj[5], sj8 = sj[8];
  double tj0 = tj[0], tj1 = tj[1], tj2 = tj[2];

  for(i = 0; i < nsub; i++){
    const int orig_ii = row_j + 1 + i;
    const int ii = use_tree ? tree_lookup[orig_ii] : orig_ii;
    const int widx = use_tree ? ii : i;
    const double w = weights[widx];

    if(w == 0.0)
      continue;

    if(track_lowsupport){
      np_lp_dense_support_add(row_j, orig_ii, ii, w, 3,
                              support_count, support_orig,
                              support_data, support_weight);
      np_lp_dense_support_add(orig_ii, row_j, eval_idx, w, 3,
                              support_count, support_orig,
                              support_data, support_weight);
    }

    {
      const double yi = response[ii];
      const double b0 = basis[0][ii];
      const double b1 = basis[1][ii];
      const double b2 = basis[2][ii];
      const double wb0 = w*b0;
      const double wb1 = w*b1;
      const double wb2 = w*b2;
      double * const si = moments + (size_t)orig_ii*9;
      double * const ti = rhs + (size_t)orig_ii*3;

      tj0 += wb0*yi;
      tj1 += wb1*yi;
      tj2 += wb2*yi;
      ti[0] += w*eval_ybasis[0];
      ti[1] += w*eval_ybasis[1];
      ti[2] += w*eval_ybasis[2];

      sj0 += wb0*b0;
      sj1 += wb0*b1;
      sj2 += wb0*b2;
      sj4 += wb1*b1;
      sj5 += wb1*b2;
      sj8 += wb2*b2;

      si[0] += w*eval_outer[0];
      si[1] += w*eval_outer[1];
      si[2] += w*eval_outer[2];
      si[4] += w*eval_outer[4];
      si[5] += w*eval_outer[5];
      si[8] += w*eval_outer[8];
    }
  }

  tj[0] = tj0;
  tj[1] = tj1;
  tj[2] = tj2;
  sj[0] = sj0;
  sj[1] = sj1;
  sj[2] = sj2;
  sj[4] = sj4;
  sj[5] = sj5;
  sj[8] = sj8;
}

attribute_hidden void
np_lp_mirror_dense_moments_row3(double *moments, const int nrows)
{
  int row;

  for(row = 0; row < nrows; row++){
    double * const s = moments + (size_t)row*9;
    s[3] = s[1];
    s[6] = s[2];
    s[7] = s[5];
  }
}

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

    if(weight == 0.0)
      continue;

    if(ctx->track_lowsupport){
      np_lp_dense_support_add(ctx->row_j, orig_ii, ii, weight, nterms,
                              ctx->support_count, ctx->support_orig,
                              ctx->support_data, ctx->support_weight);
      np_lp_dense_support_add(orig_ii, ctx->row_j, ctx->eval_idx, weight,
                              nterms, ctx->support_count, ctx->support_orig,
                              ctx->support_data, ctx->support_weight);
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

      for(a = 0; a < nterms; a++){
        const double bia = ctx->basis[a][ii];
        const double weighted_bia = weight*bia;
        fixed_rhs[a] += weighted_bia*yi;
        moving_rhs[a] += weight*ctx->eval_ybasis[a];
        for(b = 0; b < nterms; b++){
          fixed_moments[a*nterms+b] +=
            weighted_bia*ctx->basis[b][ii];
          moving_moments[a*nterms+b] +=
            weight*ctx->eval_outer[a*nterms+b];
        }
      }
    }
  }
}

void np_lp_accumulate_dense_resident_row(const NPLPDenseRowContext *ctx)
{
  switch(ctx->nterms){
  case 1: np_lp_accumulate_dense_row_1(ctx); return;
  case 2: np_lp_accumulate_dense_row_2(ctx); return;
  case 3:
    np_lp_accumulate_dense_resident_row3(
      ctx->row_j, ctx->nsub, ctx->use_tree, ctx->eval_idx,
      ctx->track_lowsupport, ctx->tree_lookup, ctx->weights, ctx->basis,
      ctx->response, ctx->moments, ctx->rhs, ctx->eval_ybasis,
      ctx->eval_outer, ctx->support_count, ctx->support_orig,
      ctx->support_data, ctx->support_weight);
    return;
  case 4: np_lp_accumulate_dense_row_4(ctx); return;
  case 5: np_lp_accumulate_dense_row_5(ctx); return;
  case 6: np_lp_accumulate_dense_row_6(ctx); return;
  case 7: np_lp_accumulate_dense_row_7(ctx); return;
  case 8: np_lp_accumulate_dense_row_8(ctx); return;
  case 9: np_lp_accumulate_dense_row_9(ctx); return;
  case 10: np_lp_accumulate_dense_row_10(ctx); return;
  case 11: np_lp_accumulate_dense_row_11(ctx); return;
  case 12: np_lp_accumulate_dense_row_12(ctx); return;
  case 13: np_lp_accumulate_dense_row_13(ctx); return;
  case 14: np_lp_accumulate_dense_row_14(ctx); return;
  case 15: np_lp_accumulate_dense_row_15(ctx); return;
  case 16: np_lp_accumulate_dense_row_16(ctx); return;
  default: np_lp_accumulate_dense_row_generic(ctx); return;
  }
}
