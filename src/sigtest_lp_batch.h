#ifndef NP_SIGTEST_LP_BATCH_H
#define NP_SIGTEST_LP_BATCH_H

#include <Rinternals.h>
#include "regression_contrast.h"

/* Invocation-owned, rank-local state. The fit supplies one accepted row at
 * a time. Storage is O(n * draws), with at most 64 responses per tile. */
typedef struct NPSigtestLPBatch {
  int n, draws, residual_rows, contrast_rows, ready, unsupported;
  const double *response;
  long double *weight, *rescale;
  NPContrastNumber *coefficient;
  double *residual, *scale, *gradient, *stderr;
} NPSigtestLPBatch;

NPSigtestLPBatch *np_sigtest_lp_begin(SEXP responses, int n);
void np_sigtest_lp_residual(NPSigtestLPBatch *tile, int self,
                            const NPResidualOffDiagonal *geometry);
void np_sigtest_lp_standardize(NPSigtestLPBatch *tile);
void np_sigtest_lp_contrast(NPSigtestLPBatch *tile, int row, int anchor, double sign);
void np_sigtest_lp_attach(NPSigtestLPBatch *tile, SEXP result);

#endif
