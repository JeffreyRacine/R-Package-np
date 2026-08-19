#ifndef NP_REGHAT_FAST_H
#define NP_REGHAT_FAST_H

#include <stddef.h>

#include "jksum_lp_solve.h"

typedef enum {
  NP_REGHAT_LP_ROW_OK = 0,
  NP_REGHAT_LP_ROW_INVALID,
  NP_REGHAT_LP_ROW_NONFINITE,
  NP_REGHAT_LP_ROW_RIDGE_FAILED,
  NP_REGHAT_LP_ROW_MEMORY
} NPReghatLPRowStatus;

/*
 * Reusable bounded workspace for a wider local-polynomial influence row.
 * The packed design and all scratch storage are O(n * p + p^2); no
 * evaluation-by-training matrix is retained here.  Width one remains owned by
 * np_lp_width_one_influence_row() and must not enter this LAPACK workspace.
 */
typedef struct {
  int ntrain;
  int nterms;
  size_t design_capacity;
  double *design;
  double *weighted_design;
  double *prediction;
  NPLPSolveWorkspace solve_workspace;
} NPReghatLPWorkspace;

void np_reghat_lp_workspace_init(NPReghatLPWorkspace *workspace);
void np_reghat_lp_workspace_clear(NPReghatLPWorkspace *workspace);

/* Copy column-fragmented package basis storage into one BLAS-ready matrix. */
NPReghatLPRowStatus np_reghat_lp_workspace_prepare_columns(
  NPReghatLPWorkspace *workspace,
  double * const *basis_columns,
  int ntrain,
  int nterms);

/*
 * Compute one exact wider-LP influence row from caller-supplied kernel
 * weights.  output_stride supports both contiguous row blocks and strided R
 * column-major matrices.  The common scaling of a complete weight row cancels
 * from a nonsingular solve, so signed-log kernel owners may pass their
 * representable common-scaled row directly.
 */
NPReghatLPRowStatus np_reghat_lp_workspace_influence_row(
  NPReghatLPWorkspace *workspace,
  const double *weights,
  const double *basis_eval,
  double *row_out,
  size_t output_stride);

#endif
