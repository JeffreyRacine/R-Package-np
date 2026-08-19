#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

#include "jksum_lp_solve.h"
#include "reghat_fast.h"

/*
 * Exact compiled counterpart of the R loop in
 * .npreghat_exact_lp_matrix_from_kernel_weights().  Width one is scalar and
 * independent of R's matprod mode; wider systems retain the BLAS-backed Gram
 * and projection route while delegating rank admission, bounded ridging, and
 * the retained factor transcript to the canonical LP solve workspace.
 */

static int np_matrix_dims(SEXP x, int *nr, int *nc)
{
  SEXP dim = getAttrib(x, R_DimSymbol);
  if((TYPEOF(dim) != INTSXP) || (XLENGTH(dim) != 2))
    return 0;
  *nr = INTEGER(dim)[0];
  *nc = INTEGER(dim)[1];
  return (*nr >= 0) && (*nc >= 0);
}

/*
 * Exact compiled counterpart of
 *
 *   sweep(t(kw), 1L, denominator, "/")
 *
 * Callers retain ownership of the denominator calculation and any zero-floor
 * policy.  This helper changes only the transpose/allocation traversal.  The
 * scalar division for every output element is therefore identical to the R
 * route, while small tiles keep both the strided input and contiguous output
 * writes resident in cache.
 */
SEXP C_np_lc_hat_normalize(SEXP kw, SEXP denominator)
{
  const int tile = 32;
  int ntrain = 0, neval = 0;
  const double *kw_p = NULL;
  const double *denominator_p = NULL;
  double *out_p = NULL;
  SEXP out = R_NilValue;

  if((TYPEOF(kw) != REALSXP) || (TYPEOF(denominator) != REALSXP) ||
     !np_matrix_dims(kw, &ntrain, &neval) ||
     (XLENGTH(denominator) != (R_xlen_t)neval))
    error("invalid LC hat normalization input");

  kw_p = REAL(kw);
  denominator_p = REAL(denominator);
  out = PROTECT(allocMatrix(REALSXP, neval, ntrain));
  out_p = REAL(out);

  for(int jb = 0; jb < neval; jb += tile){
    const int jend = (jb + tile < neval) ? jb + tile : neval;
    for(int ib = 0; ib < ntrain; ib += tile){
      const int iend = (ib + tile < ntrain) ? ib + tile : ntrain;
      for(int i = ib; i < iend; i++)
        for(int j = jb; j < jend; j++)
          out_p[j + (size_t)neval*(size_t)i] =
            kw_p[i + (size_t)ntrain*(size_t)j]/denominator_p[j];
    }
    R_CheckUserInterrupt();
  }

  UNPROTECT(1);
  return out;
}

static SEXP np_reghat_width_one_matrix(SEXP kw,
                                       SEXP wtrain,
                                       SEXP weval,
                                       const int ntrain,
                                       const int neval)
{
  SEXP out = PROTECT(allocMatrix(REALSXP, neval, ntrain));

  for(int j = 0; j < neval; j++){
    const double * const weights = REAL(kw) + (size_t)j*(size_t)ntrain;
    NPLPWidthOneStatus status;

    if((j == 0) || (j + 1 == neval) || ((j % 32) == 0))
      R_CheckUserInterrupt();

    status = np_lp_width_one_influence_row(
      REAL(wtrain),
      ntrain,
      weights,
      REAL(weval)[j],
      REAL(out) + j,
      (size_t)neval
    );

    if(status == NP_LP_WIDTH_ONE_OK)
      continue;
    if(status == NP_LP_WIDTH_ONE_NONFINITE)
      error("LP solve failed in compiled hat-matrix path: non-finite system");
    if(status == NP_LP_WIDTH_ONE_RIDGE_FAILED)
      error("LP solve failed in compiled hat-matrix path after bounded ridging");
    error("invalid width-one compiled hat-matrix input");
  }

  UNPROTECT(1);
  return out;
}

static NPReghatLPRowStatus np_reghat_lp_prediction_raw(
  const int ntrain,
  const int nterms,
  const double * const design,
  const double * const weights,
  const double * const basis_eval,
  double * const weighted_design,
  NPLPSolveWorkspace * const solve_workspace,
  double * const prediction,
  const int check_interrupt)
{
  const double alpha = 1.0;
  const double beta = 0.0;
  const double ridge_fraction = 1.0/(double)ntrain;
  const char trans_t = 'T';
  const char trans_n = 'N';
  const int one = 1;
  NPLPSolvePolicyDiagnostics diagnostics;
  int term;
  int i;

  if((ntrain <= 0) || (nterms <= 1) || (design == NULL) ||
     (weights == NULL) || (basis_eval == NULL) ||
     (weighted_design == NULL) || (solve_workspace == NULL) ||
     (prediction == NULL))
    return NP_REGHAT_LP_ROW_INVALID;

  for(term = 0; term < nterms; term++) {
    const double * const src = design + (size_t)term*(size_t)ntrain;
    double * const dst = weighted_design + (size_t)term*(size_t)ntrain;

    for(i = 0; i < ntrain; i++)
      dst[i] = src[i]*weights[i];
  }

  F77_CALL(dgemm)(&trans_t, &trans_n,
                  &nterms, &nterms, &ntrain,
                  &alpha, design, &ntrain,
                  weighted_design, &ntrain,
                  &beta, solve_workspace->gram_source,
                  &nterms FCONE FCONE);
  memcpy(solve_workspace->rhs_source, basis_eval,
         (size_t)nterms*sizeof(double));

  if(check_interrupt)
    R_CheckUserInterrupt();
  switch(np_lp_solve_workspace_solve_adjoint_ranked(
           solve_workspace, nterms, 1, ridge_fraction,
           NP_LP_RANK_UPPER_BOUND_UNKNOWN, &diagnostics)) {
  case NP_LP_SOLVE_POLICY_OK:
    break;
  case NP_LP_SOLVE_POLICY_NONFINITE:
    return NP_REGHAT_LP_ROW_NONFINITE;
  case NP_LP_SOLVE_POLICY_RIDGE_EXHAUSTED:
  case NP_LP_SOLVE_POLICY_FINAL_FAILED:
    return NP_REGHAT_LP_ROW_RIDGE_FAILED;
  default:
    return NP_REGHAT_LP_ROW_INVALID;
  }

  F77_CALL(dgemv)(&trans_n, &ntrain, &nterms, &alpha,
                  design, &ntrain, solve_workspace->rhs_work, &one,
                  &beta, prediction, &one FCONE);
  return NP_REGHAT_LP_ROW_OK;
}

void np_reghat_lp_workspace_init(NPReghatLPWorkspace *workspace)
{
  if(workspace != NULL) {
    memset(workspace, 0, sizeof(*workspace));
    np_lp_solve_workspace_init(&workspace->solve_workspace);
  }
}

void np_reghat_lp_workspace_clear(NPReghatLPWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  free(workspace->design);
  free(workspace->weighted_design);
  free(workspace->prediction);
  np_lp_solve_workspace_clear(&workspace->solve_workspace);
  np_reghat_lp_workspace_init(workspace);
}

static NPReghatLPRowStatus np_reghat_lp_workspace_reserve(
  NPReghatLPWorkspace *workspace,
  int ntrain,
  int nterms)
{
  size_t design_elements;
  double *design = NULL;
  double *weighted_design = NULL;
  double *prediction = NULL;
  NPLPSolveWorkspace solve_workspace;

  if((workspace == NULL) || (ntrain <= 0) || (nterms <= 1) ||
     ((size_t)ntrain > SIZE_MAX/(size_t)nterms))
    return NP_REGHAT_LP_ROW_INVALID;
  design_elements = (size_t)ntrain*(size_t)nterms;
  if((design_elements > SIZE_MAX/sizeof(double)) ||
     ((size_t)ntrain > SIZE_MAX/sizeof(double)))
    return NP_REGHAT_LP_ROW_INVALID;
  if((workspace->ntrain == ntrain) && (workspace->nterms == nterms) &&
     (workspace->design_capacity >= design_elements) &&
     workspace->design != NULL && workspace->weighted_design != NULL &&
     workspace->prediction != NULL &&
     workspace->solve_workspace.gram_source != NULL &&
     workspace->solve_workspace.rhs_source != NULL)
    return NP_REGHAT_LP_ROW_OK;

  np_lp_solve_workspace_init(&solve_workspace);
  design = (double *)malloc(design_elements*sizeof(double));
  weighted_design = (double *)malloc(design_elements*sizeof(double));
  prediction = (double *)malloc((size_t)ntrain*sizeof(double));
  if(design == NULL || weighted_design == NULL || prediction == NULL ||
     !np_lp_solve_workspace_reserve(&solve_workspace, nterms, 1)) {
    free(design);
    free(weighted_design);
    free(prediction);
    np_lp_solve_workspace_clear(&solve_workspace);
    return NP_REGHAT_LP_ROW_MEMORY;
  }

  np_reghat_lp_workspace_clear(workspace);
  workspace->ntrain = ntrain;
  workspace->nterms = nterms;
  workspace->design_capacity = design_elements;
  workspace->design = design;
  workspace->weighted_design = weighted_design;
  workspace->prediction = prediction;
  workspace->solve_workspace = solve_workspace;
  return NP_REGHAT_LP_ROW_OK;
}

NPReghatLPRowStatus np_reghat_lp_workspace_prepare_columns(
  NPReghatLPWorkspace *workspace,
  double * const *basis_columns,
  int ntrain,
  int nterms)
{
  NPReghatLPRowStatus status;
  int term;

  if(basis_columns == NULL)
    return NP_REGHAT_LP_ROW_INVALID;
  status = np_reghat_lp_workspace_reserve(workspace, ntrain, nterms);
  if(status != NP_REGHAT_LP_ROW_OK)
    return status;
  for(term = 0; term < nterms; term++) {
    if(basis_columns[term] == NULL)
      return NP_REGHAT_LP_ROW_INVALID;
    memcpy(workspace->design + (size_t)term*(size_t)ntrain,
           basis_columns[term], (size_t)ntrain*sizeof(double));
  }
  return NP_REGHAT_LP_ROW_OK;
}

NPReghatLPRowStatus np_reghat_lp_workspace_influence_row(
  NPReghatLPWorkspace *workspace,
  const double *weights,
  const double *basis_eval,
  double *row_out,
  size_t output_stride)
{
  NPReghatLPRowStatus status;
  int observation;

  if((workspace == NULL) || (workspace->ntrain <= 0) ||
     (workspace->nterms <= 1))
    return NP_REGHAT_LP_ROW_INVALID;
  if(row_out == NULL || output_stride == 0U)
    return NP_REGHAT_LP_ROW_INVALID;
  status = np_reghat_lp_prediction_raw(
    workspace->ntrain, workspace->nterms, workspace->design,
    weights, basis_eval, workspace->weighted_design,
    &workspace->solve_workspace, workspace->prediction, 0);
  if(status != NP_REGHAT_LP_ROW_OK)
    return status;
  for(observation = 0; observation < workspace->ntrain; ++observation) {
    const double value =
      weights[observation]*workspace->prediction[observation];

    if(!isfinite(value))
      return NP_REGHAT_LP_ROW_NONFINITE;
    row_out[(size_t)observation*output_stride] = value;
  }
  return NP_REGHAT_LP_ROW_OK;
}

typedef struct {
  SEXP kw;
  SEXP weval;
  int ntrain;
  int neval;
  int nterms;
  NPReghatLPWorkspace workspace;
} NPReghatMatrixExecution;

static void np_reghat_matrix_execution_cleanup(void *data, Rboolean jump)
{
  NPReghatMatrixExecution * const execution =
    (NPReghatMatrixExecution *)data;
  (void)jump;
  np_reghat_lp_workspace_clear(&execution->workspace);
}

static SEXP np_reghat_matrix_execution_run(void *data)
{
  NPReghatMatrixExecution * const execution =
    (NPReghatMatrixExecution *)data;
  const int ntrain = execution->ntrain;
  const int neval = execution->neval;
  const int nterms = execution->nterms;
  NPReghatLPWorkspace * const workspace = &execution->workspace;
  SEXP out = PROTECT(allocMatrix(REALSXP, neval, ntrain));

  for(int j = 0; j < neval; j++){
    const double * const weights =
      REAL(execution->kw) + (size_t)j*(size_t)ntrain;
    NPReghatLPRowStatus row_status;

    if((j == 0) || (j + 1 == neval) || ((j % 32) == 0))
      R_CheckUserInterrupt();
    for(int term = 0; term < nterms; term++)
      workspace->solve_workspace.rhs_work[term] =
        REAL(execution->weval)[j + (size_t)neval*(size_t)term];
    row_status = np_reghat_lp_prediction_raw(
      ntrain, nterms, workspace->design, weights,
      workspace->solve_workspace.rhs_work, workspace->weighted_design,
      &workspace->solve_workspace, workspace->prediction, 0);
    if(row_status == NP_REGHAT_LP_ROW_NONFINITE)
      error("LP solve failed in compiled hat-matrix path: non-finite system");
    if(row_status == NP_REGHAT_LP_ROW_RIDGE_FAILED)
      error("LP solve failed in compiled hat-matrix path after bounded ridging");
    if(row_status != NP_REGHAT_LP_ROW_OK)
      error("invalid wider-LP compiled hat-matrix input");
    for(int i = 0; i < ntrain; i++)
      REAL(out)[j + (size_t)neval*(size_t)i] =
        weights[i]*workspace->prediction[i];
  }

  UNPROTECT(1);
  return out;
}

SEXP C_np_reghat_lp_matrix_fast(SEXP kw, SEXP wtrain, SEXP weval)
{
  int ntrain = 0, neval = 0, kw_neval = 0;
  int wtrain_n = 0, nterms = 0, weval_n = 0, weval_p = 0;
  NPReghatMatrixExecution execution;

  memset(&execution, 0, sizeof(execution));
  np_reghat_lp_workspace_init(&execution.workspace);

  if((TYPEOF(kw) != REALSXP) || (TYPEOF(wtrain) != REALSXP) ||
     (TYPEOF(weval) != REALSXP) ||
     !np_matrix_dims(kw, &ntrain, &kw_neval) ||
     !np_matrix_dims(wtrain, &wtrain_n, &nterms) ||
     !np_matrix_dims(weval, &weval_n, &weval_p) ||
     (ntrain <= 0) || (kw_neval <= 0) || (nterms <= 0) ||
     (wtrain_n != ntrain) || (weval_n != kw_neval) ||
     (weval_p != nterms))
    return R_NilValue;

  neval = kw_neval;
  if(nterms == 1)
    return np_reghat_width_one_matrix(kw, wtrain, weval, ntrain, neval);

  if(((size_t)ntrain > SIZE_MAX/(size_t)nterms) ||
     ((size_t)neval > SIZE_MAX/(size_t)ntrain) ||
     np_reghat_lp_workspace_reserve(&execution.workspace, ntrain, nterms) !=
       NP_REGHAT_LP_ROW_OK)
    return R_NilValue;
  memcpy(execution.workspace.design, REAL(wtrain),
         (size_t)ntrain*(size_t)nterms*sizeof(double));
  execution.kw = kw;
  execution.weval = weval;
  execution.ntrain = ntrain;
  execution.neval = neval;
  execution.nterms = nterms;
  return R_UnwindProtect(
    np_reghat_matrix_execution_run, &execution,
    np_reghat_matrix_execution_cleanup, &execution, NULL);
}
