#ifndef NP_JKSUM_LP_ROW_H
#define NP_JKSUM_LP_ROW_H

#include <math.h>
#include <stddef.h>
#include <R_ext/Visibility.h>

#if defined(__aarch64__) && defined(NP_USE_ACCELERATE_GAUSS) && NP_USE_ACCELERATE_GAUSS
#include <arm_neon.h>
#define NP_LP_ROW_NEON 1
#else
#define NP_LP_ROW_NEON 0
#endif

/*
 * Exact denominator contract for the full-row/delete-one LP identity.
 * Every finite nonzero signed denominator is valid and remains unchanged;
 * zero and non-finite values denote a failed deleted system.
 */
static inline int np_lp_delete_denominator(const double leverage,
                                           double *denominator)
{
  if(denominator == NULL)
    return 0;
  *denominator = 1.0 - leverage;
  return isfinite(*denominator) && (*denominator != 0.0);
}

typedef struct {
  int nterms;
  int row_j;
  int nsub;
  int use_tree;
  int eval_idx;
  int track_lowsupport;
  const int *tree_lookup;
  const double *weights;
  double * const *basis;
  const double *response;
  double *moments;
  double *rhs;
  const double *eval_ybasis;
  const double *eval_outer;
  int *support_count;
} NPLPDenseRowContext;

/*
 * Width six is the established resident/packed crossover boundary and is
 * common for a bivariate quadratic generalized basis. Keep its sparse
 * unordered-pair transcript compile-time visible so the resident row stays
 * in registers and the contiguous moving row can use two-lane SIMD on
 * AArch64. Tree support, pair ownership, and sum order are unchanged.
 */
#if NP_LP_ROW_NEON
static inline void np_lp_accumulate_sparse_pair_resident6(
    double * const *basis,
    const double *response,
    double *moments,
    double *rhs,
    double *row_moments,
    double *row_rhs,
    const double *eval_ybasis,
    const double *eval_outer,
    const int orig_ii,
    const int tree_ii,
    const double weight)
{
  enum { nterms = 6 };
  const double yi = response[tree_ii];
  double * const moving_moments = moments +
    (size_t)orig_ii*(size_t)nterms*(size_t)nterms;
  double * const moving_rhs = rhs + (size_t)orig_ii*(size_t)nterms;
  int a, b;

  for(a = 0; a < nterms; a++){
    const double bia = basis[a][tree_ii];
    const double weighted_bia = weight*bia;
    const int aoff = a*nterms;

    row_rhs[a] += weighted_bia*yi;
    for(b = a; b < nterms; b++)
      row_moments[aoff+b] += weighted_bia*basis[b][tree_ii];
  }

  {
    const float64x2_t vw = vdupq_n_f64(weight);

    for(a = 0; a + 1 < nterms; a += 2)
      vst1q_f64(moving_rhs + a,
                vfmaq_f64(vld1q_f64(moving_rhs + a), vw,
                          vld1q_f64(eval_ybasis + a)));
    if(a < nterms)
      moving_rhs[a] += weight*eval_ybasis[a];

    for(a = 0; a < nterms; a++){
      const int end = a*nterms + nterms;
      int pos = a*nterms + a;

      for(; pos + 1 < end; pos += 2)
        vst1q_f64(moving_moments + pos,
                  vfmaq_f64(vld1q_f64(moving_moments + pos), vw,
                            vld1q_f64(eval_outer + pos)));
      if(pos < end)
        moving_moments[pos] += weight*eval_outer[pos];
    }
  }
}

/* Rank-owned sparse rows omit the symmetric partner update. */
static inline void np_lp_accumulate_sparse_row_resident6(
    double * const *basis,
    const double *response,
    double *row_moments,
    double *row_rhs,
    const int tree_ii,
    const double weight)
{
  enum { nterms = 6 };
  const double yi = response[tree_ii];
  int a, b;

  for(a = 0; a < nterms; a++){
    const double bia = basis[a][tree_ii];
    const double weighted_bia = weight*bia;
    const int aoff = a*nterms;

    row_rhs[a] += weighted_bia*yi;
    for(b = a; b < nterms; b++)
      row_moments[aoff+b] += weighted_bia*basis[b][tree_ii];
  }
}
#endif

/*
  A rank-owned MPI row is accumulated independently of every other row.
  Width one has a dedicated scalar caller and is deliberately excluded from
  this wider-row context; width seven and above retain the generic caller loop.
*/
typedef struct {
  int nterms;
  int row_j;
  int nobs;
  int use_tree;
  int track_lowsupport;
  const int *tree_lookup;
  const double *weights;
  double * const *basis;
  const double *response;
  double *moments;
  double *rhs;
  int *support_count;
} NPLPOwnedRowContext;

void np_lp_accumulate_dense_resident_row(const NPLPDenseRowContext *ctx);

attribute_hidden int
np_lp_accumulate_owned_resident_row(const NPLPOwnedRowContext *ctx);

attribute_hidden void np_lp_mirror_dense_moments_row3(double *moments,
                                                       int nrows);

void np_lp_accumulate_dense_resident_row3(
  int row_j,
  int nsub,
  int use_tree,
  int eval_idx,
  int track_lowsupport,
  const int *tree_lookup,
  const double *weights,
  double * const *basis,
  const double *response,
  double *moments,
  double *rhs,
  const double *eval_ybasis,
  const double *eval_outer,
  int *support_count);

#endif
