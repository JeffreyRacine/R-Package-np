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
 * independent of R's matprod mode; wider systems retain the incumbent
 * BLAS-backed route.  Matrix products, LAPACK condition checks, ridge
 * increments, and final row construction deliberately retain the incumbent
 * operation order for those wider systems.
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

static NPReghatLPRowStatus np_reghat_lp_prediction_raw(
  const int ntrain,
  const int nterms,
  const double * const design,
  const double * const weights,
  const double * const basis_eval,
  double * const weighted_design,
  double * const prediction,
  NPLPSolveWorkspace * const solve_workspace,
  NPLPSolvePolicyDiagnostics * const diagnostics)
{
  const double alpha = 1.0;
  const double beta = 0.0;
  const double epsilon = 1.0/(double)ntrain;
  const char trans_t = 'T';
  const char trans_n = 'N';
  const int one = 1;
  int term;
  int i;

  if((ntrain <= 0) || (nterms <= 1) || (design == NULL) ||
     (weights == NULL) || (basis_eval == NULL) ||
     (weighted_design == NULL) || (prediction == NULL) ||
     (solve_workspace == NULL) ||
     (solve_workspace->gram_source == NULL) ||
     (solve_workspace->rhs_source == NULL) ||
     (solve_workspace->rhs_work == NULL))
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

  switch(np_lp_solve_workspace_solve_adjoint(solve_workspace,
                                              nterms,
                                              1,
                                              epsilon,
                                              diagnostics)) {
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
     np_lp_solve_workspace_reserve(&workspace->solve_workspace, nterms, 1))
    return NP_REGHAT_LP_ROW_OK;

  design = (double *)malloc(design_elements*sizeof(double));
  weighted_design = (double *)malloc(design_elements*sizeof(double));
  prediction = (double *)malloc((size_t)ntrain*sizeof(double));
  if(design == NULL || weighted_design == NULL || prediction == NULL ||
     !np_lp_solve_workspace_reserve(&workspace->solve_workspace, nterms, 1)) {
    free(design);
    free(weighted_design);
    free(prediction);
    return NP_REGHAT_LP_ROW_MEMORY;
  }

  free(workspace->design);
  free(workspace->weighted_design);
  free(workspace->prediction);
  workspace->ntrain = ntrain;
  workspace->nterms = nterms;
  workspace->design_capacity = design_elements;
  workspace->design = design;
  workspace->weighted_design = weighted_design;
  workspace->prediction = prediction;
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
    workspace->prediction, &workspace->solve_workspace, NULL);
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
  int ntrain;
  int nterms;
  const double *design;
  double *weighted_design;
  double *prediction;
  double *eval_basis;
  double *influence_row;
  NPLPSolveWorkspace solve_workspace;
} NPReghatCallWorkspace;

static int np_reghat_call_dims(SEXP kw,
                               SEXP wtrain,
                               SEXP weval,
                               int *ntrain,
                               int *neval,
                               int *nterms)
{
  int kw_neval = 0;
  int wtrain_n = 0;
  int weval_n = 0;
  int weval_p = 0;

  if((TYPEOF(kw) != REALSXP) || (TYPEOF(wtrain) != REALSXP) ||
     (TYPEOF(weval) != REALSXP) ||
     !np_matrix_dims(kw, ntrain, &kw_neval) ||
     !np_matrix_dims(wtrain, &wtrain_n, nterms) ||
     !np_matrix_dims(weval, &weval_n, &weval_p) ||
     (*ntrain <= 0) || (kw_neval <= 0) || (*nterms <= 0) ||
     (wtrain_n != *ntrain) || (weval_n != kw_neval) ||
     (weval_p != *nterms))
    return 0;
  if(((size_t)*nterms > SIZE_MAX/sizeof(double)) ||
     ((size_t)*ntrain > SIZE_MAX/((size_t)*nterms*sizeof(double))) ||
     ((size_t)*nterms > SIZE_MAX/((size_t)*nterms*sizeof(double))))
    return 0;
  *neval = kw_neval;
  return 1;
}

static void np_reghat_call_workspace_bind(NPReghatCallWorkspace *workspace,
                                           const double *design,
                                           int ntrain,
                                           int nterms)
{
  memset(workspace, 0, sizeof(*workspace));
  workspace->ntrain = ntrain;
  workspace->nterms = nterms;
  workspace->design = design;
  workspace->eval_basis = (double *)R_alloc((size_t)nterms, sizeof(double));
  workspace->influence_row = (double *)R_alloc((size_t)ntrain, sizeof(double));

  if(nterms <= 1)
    return;
  workspace->weighted_design = (double *)R_alloc(
    (size_t)ntrain*(size_t)nterms, sizeof(double));
  workspace->prediction = (double *)R_alloc((size_t)ntrain, sizeof(double));
  workspace->solve_workspace.p_capacity = nterms;
  workspace->solve_workspace.nrhs_capacity = 1;
  workspace->solve_workspace.gram_capacity = (size_t)nterms*(size_t)nterms;
  workspace->solve_workspace.rhs_capacity = (size_t)nterms;
  workspace->solve_workspace.gram_source = (double *)R_alloc(
    workspace->solve_workspace.gram_capacity, sizeof(double));
  workspace->solve_workspace.rhs_source = (double *)R_alloc(
    workspace->solve_workspace.rhs_capacity, sizeof(double));
  workspace->solve_workspace.gram_work = (double *)R_alloc(
    workspace->solve_workspace.gram_capacity, sizeof(double));
  workspace->solve_workspace.rhs_work = (double *)R_alloc(
    workspace->solve_workspace.rhs_capacity, sizeof(double));
  workspace->solve_workspace.ipiv = (int *)R_alloc(
    (size_t)nterms, sizeof(int));
  workspace->solve_workspace.rcond_work = (double *)R_alloc(
    4U*(size_t)nterms, sizeof(double));
  workspace->solve_workspace.rcond_iwork = (int *)R_alloc(
    (size_t)nterms, sizeof(int));
}

static NPReghatLPRowStatus np_reghat_call_influence_row(
  NPReghatCallWorkspace *workspace,
  const double *weights,
  const double *weval,
  int eval_index,
  int neval,
  NPLPSolvePolicyDiagnostics *diagnostics)
{
  NPReghatLPRowStatus status;

  diagnostics->ridge_steps = 0;
  diagnostics->ridge_total = 0.0;
  for(int term = 0; term < workspace->nterms; term++)
    workspace->eval_basis[term] =
      weval[eval_index + (size_t)neval*(size_t)term];

  if(workspace->nterms == 1) {
    const NPLPWidthOneStatus width_status = np_lp_width_one_influence_row(
      workspace->design, workspace->ntrain, weights,
      workspace->eval_basis[0], workspace->influence_row, 1U, diagnostics);

    if(width_status == NP_LP_WIDTH_ONE_OK)
      return NP_REGHAT_LP_ROW_OK;
    if(width_status == NP_LP_WIDTH_ONE_NONFINITE)
      return NP_REGHAT_LP_ROW_NONFINITE;
    if(width_status == NP_LP_WIDTH_ONE_RIDGE_FAILED)
      return NP_REGHAT_LP_ROW_RIDGE_FAILED;
    return NP_REGHAT_LP_ROW_INVALID;
  }

  status = np_reghat_lp_prediction_raw(
    workspace->ntrain, workspace->nterms, workspace->design, weights,
    workspace->eval_basis, workspace->weighted_design, workspace->prediction,
    &workspace->solve_workspace, diagnostics);
  if(status != NP_REGHAT_LP_ROW_OK)
    return status;
  for(int observation = 0; observation < workspace->ntrain; observation++) {
    const double value = weights[observation]*workspace->prediction[observation];

    if(!isfinite(value))
      return NP_REGHAT_LP_ROW_NONFINITE;
    workspace->influence_row[observation] = value;
  }
  return NP_REGHAT_LP_ROW_OK;
}

static void np_reghat_check_row_status(NPReghatLPRowStatus status,
                                       const char *path)
{
  if(status == NP_REGHAT_LP_ROW_OK)
    return;
  if(status == NP_REGHAT_LP_ROW_NONFINITE)
    error("LP solve failed in compiled %s path: non-finite system", path);
  if(status == NP_REGHAT_LP_ROW_RIDGE_FAILED)
    error("LP solve failed in compiled %s path after bounded ridging", path);
  error("invalid compiled LP %s input", path);
}

SEXP C_np_reghat_lp_matrix_fast(SEXP kw, SEXP wtrain, SEXP weval)
{
  int ntrain = 0, neval = 0, nterms = 0;
  NPReghatCallWorkspace workspace;
  SEXP out = R_NilValue;
  SEXP ridge = R_NilValue;

  if(!np_reghat_call_dims(kw, wtrain, weval,
                          &ntrain, &neval, &nterms) ||
     ((size_t)neval > SIZE_MAX/(size_t)ntrain))
    return R_NilValue;

  np_reghat_call_workspace_bind(&workspace, REAL(wtrain), ntrain, nterms);

  out = PROTECT(allocMatrix(REALSXP, neval, ntrain));
  ridge = PROTECT(allocVector(REALSXP, neval));
  for(int j = 0; j < neval; j++){
    const double * const weights = REAL(kw) + (size_t)j*(size_t)ntrain;
    NPLPSolvePolicyDiagnostics diagnostics;

    if((j == 0) || (j + 1 == neval) || ((j % 32) == 0))
      R_CheckUserInterrupt();
    np_reghat_check_row_status(
      np_reghat_call_influence_row(
        &workspace, weights, REAL(weval), j, neval, &diagnostics),
      "hat-matrix");
    REAL(ridge)[j] = diagnostics.ridge_total;
    for(int i = 0; i < ntrain; i++)
      REAL(out)[j + (size_t)neval*(size_t)i] = workspace.influence_row[i];
  }

  setAttrib(out, install("ridge.used"), ridge);
  UNPROTECT(2);
  return out;
}

SEXP C_np_reghat_lp_apply_fast(SEXP kw, SEXP wtrain, SEXP weval, SEXP y)
{
  int ntrain = 0, neval = 0, nterms = 0;
  int y_n = 0, nrhs = 0;
  NPReghatCallWorkspace workspace;
  SEXP out = R_NilValue;
  SEXP ridge = R_NilValue;

  if(!np_reghat_call_dims(kw, wtrain, weval,
                          &ntrain, &neval, &nterms) ||
     (TYPEOF(y) != REALSXP) ||
     !np_matrix_dims(y, &y_n, &nrhs) ||
     (nrhs <= 0) || (y_n != ntrain) ||
     ((size_t)neval > SIZE_MAX/(size_t)nrhs))
    return R_NilValue;

  np_reghat_call_workspace_bind(&workspace, REAL(wtrain), ntrain, nterms);

  out = PROTECT(allocMatrix(REALSXP, neval, nrhs));
  ridge = PROTECT(allocVector(REALSXP, neval));
  for(int j = 0; j < neval; j++) {
    const double * const weights = REAL(kw) + (size_t)j*(size_t)ntrain;
    NPLPSolvePolicyDiagnostics diagnostics = {0, 0.0};
    NPReghatLPRowStatus row_status = NP_REGHAT_LP_ROW_OK;

    if((j == 0) || (j + 1 == neval) || ((j % 32) == 0))
      R_CheckUserInterrupt();
    row_status = np_reghat_call_influence_row(
      &workspace, weights, REAL(weval), j, neval, &diagnostics);
    np_reghat_check_row_status(row_status, "hat-apply");

    REAL(ridge)[j] = diagnostics.ridge_total;
    for(int rhs = 0; rhs < nrhs; rhs++) {
      double value = 0.0;
      const double * const response = REAL(y) + (size_t)rhs*(size_t)ntrain;

      for(int i = 0; i < ntrain; i++)
        value += workspace.influence_row[i]*response[i];
      if(!isfinite(value))
        error("LP solve failed in compiled hat-apply path: non-finite result");
      REAL(out)[j + (size_t)neval*(size_t)rhs] = value;
    }
  }

  setAttrib(out, install("ridge.used"), ridge);
  UNPROTECT(2);
  return out;
}
