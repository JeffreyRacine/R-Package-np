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

typedef struct { double value; double log_value; } NPRlyNormalizer;

/* Borrowed view for one resident-row accumulation call. All pointed-to
 * storage remains caller-owned and must outlive the call. Inputs (including
 * basis, despite its historical non-const pointer type) are read-only here.
 * moments/rhs and, when enabled, support_count are updated additively for
 * both endpoints of each visited pair; the caller initializes these buffers.
 * No pointer is retained and no heap storage is allocated by the row helpers.
 * The enclosing jksum.c owner handles triangle completion, any rank reduction,
 * solve policy and cleanup. See src/README.md for the surrounding route map. */
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
  const NPRlyNormalizer *rly_log_normalizer;
} NPLPDenseRowContext;

/* Reverse a directed pair without recomputing its continuous or categorical
 * kernel. NULL preserves the symmetric transcript. The log-domain branch
 * avoids overflowing the ratio when its final weighted value is representable. */
static inline double np_lp_reverse_pair_weight(const double weight,
    const NPRlyNormalizer *log_normalizer, const int donor, const int evaluation)
{
  if(log_normalizer == NULL || weight == 0.0) return weight;
  const double numerator = log_normalizer[donor].value;
  const double denominator = log_normalizer[evaluation].value;
  if(isfinite(numerator) && isfinite(denominator))
    return weight*(numerator/denominator);
  const double difference = log_normalizer[donor].log_value-log_normalizer[evaluation].log_value;
  if(fabs(difference) < 350.0) return weight*exp(difference);
  return copysign(exp(log(fabs(weight))+difference),weight);
}

/*
 * Width six is the established resident/packed crossover boundary and is
 * common for a bivariate quadratic generalized basis.  Keep its sparse
 * unordered-pair transcript compile-time visible so the resident row stays
 * in registers and the contiguous moving row can use two-lane SIMD on
 * AArch64.  Tree support, pair ownership, and sum order are unchanged.
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
    const double weight,
    const double reverse_weight)
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
    const float64x2_t vw = vdupq_n_f64(reverse_weight);

    for(a = 0; a + 1 < nterms; a += 2)
      vst1q_f64(moving_rhs + a,
                vfmaq_f64(vld1q_f64(moving_rhs + a), vw,
                          vld1q_f64(eval_ybasis + a)));
    if(a < nterms)
      moving_rhs[a] += reverse_weight*eval_ybasis[a];

    for(a = 0; a < nterms; a++){
      const int end = a*nterms + nterms;
      int pos = a*nterms + a;

      for(; pos + 1 < end; pos += 2)
        vst1q_f64(moving_moments + pos,
                  vfmaq_f64(vld1q_f64(moving_moments + pos), vw,
                            vld1q_f64(eval_outer + pos)));
      if(pos < end)
        moving_moments[pos] += reverse_weight*eval_outer[pos];
    }
  }
}
#endif

void np_lp_accumulate_dense_resident_row(const NPLPDenseRowContext *ctx);

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
  int *support_count,
  const NPRlyNormalizer *rly_log_normalizer);

#endif
