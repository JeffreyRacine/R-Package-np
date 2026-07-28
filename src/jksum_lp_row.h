#ifndef NP_JKSUM_LP_ROW_H
#define NP_JKSUM_LP_ROW_H

#include <math.h>
#include <R_ext/Visibility.h>

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
  int *support_orig;
  int *support_data;
  double *support_weight;
} NPLPDenseRowContext;

void np_lp_accumulate_dense_resident_row(const NPLPDenseRowContext *ctx);

attribute_hidden void np_lp_mirror_dense_moments_row3(double *moments,
                                                       int nrows);

attribute_hidden int np_lp_delete_smoother_row(const double *full_row,
                                                int n,
                                                int eval_idx,
                                                double *loo_row);

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
  int *support_orig,
  int *support_data,
  double *support_weight);

#endif
