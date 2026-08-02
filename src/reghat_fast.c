#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
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

static int np_reghat_solve_system(const int nterms,
                                  const double * const matrix,
                                  const double * const rhs,
                                  double * const matrix_work,
                                  double * const solution,
                                  int * const pivot,
                                  double * const condition_work,
                                  int * const condition_iwork)
{
  const char norm = '1';
  const int nrhs = 1;
  int info = 0;
  double anorm = 0.0;
  double rcond = 0.0;

  memcpy(matrix_work, matrix,
         (size_t)nterms*(size_t)nterms*sizeof(double));
  memcpy(solution, rhs, (size_t)nterms*sizeof(double));
  anorm = F77_CALL(dlange)(&norm, &nterms, &nterms, matrix_work, &nterms,
                           condition_work FCONE);
  F77_CALL(dgesv)(&nterms, &nrhs, matrix_work, &nterms, pivot,
                  solution, &nterms, &info);
  if(info != 0)
    return 0;

  F77_CALL(dgecon)(&norm, &nterms, matrix_work, &nterms, &anorm, &rcond,
                   condition_work, condition_iwork, &info FCONE);
  if((info != 0) || !isfinite(rcond) || (rcond < DBL_EPSILON))
    return 0;
  for(int term = 0; term < nterms; term++)
    if(!isfinite(solution[term]))
      return 0;
  return 1;
}

static int np_reghat_sources_finite(const int nterms,
                                    const double * const matrix,
                                    const double * const rhs)
{
  const size_t matrix_elements = (size_t)nterms*(size_t)nterms;

  for(size_t i = 0; i < matrix_elements; i++)
    if(!isfinite(matrix[i]))
      return 0;
  for(int i = 0; i < nterms; i++)
    if(!isfinite(rhs[i]))
      return 0;
  return 1;
}

static NPReghatLPRowStatus np_reghat_lp_prediction_raw(
  const int ntrain,
  const int nterms,
  const double * const design,
  const double * const weights,
  const double * const basis_eval,
  double * const weighted_design,
  double * const gram,
  double * const gram_work,
  double * const rhs,
  double * const solution,
  double * const prediction,
  double * const condition_work,
  int * const pivot,
  int * const condition_iwork,
  const int check_interrupt)
{
  const double alpha = 1.0;
  const double beta = 0.0;
  const double epsilon = 1.0/(double)ntrain;
  const char trans_t = 'T';
  const char trans_n = 'N';
  const int one = 1;
  double nepsilon = 0.0;
  int term;
  int i;

  if((ntrain <= 0) || (nterms <= 1) || (design == NULL) ||
     (weights == NULL) || (basis_eval == NULL) ||
     (weighted_design == NULL) || (gram == NULL) ||
     (gram_work == NULL) || (rhs == NULL) || (solution == NULL) ||
     (prediction == NULL) || (condition_work == NULL) ||
     (pivot == NULL) || (condition_iwork == NULL))
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
                  &beta, gram, &nterms FCONE FCONE);
  if(rhs != basis_eval)
    memcpy(rhs, basis_eval, (size_t)nterms*sizeof(double));

  if(!np_reghat_solve_system(nterms, gram, rhs, gram_work, solution,
                             pivot, condition_work, condition_iwork)) {
    int ridge_step;
    int solved = 0;

    if(!np_reghat_sources_finite(nterms, gram, rhs))
      return NP_REGHAT_LP_ROW_NONFINITE;
    for(ridge_step = 0;
        ridge_step < NP_LP_SOLVE_MAX_RIDGE_STEPS;
        ridge_step++) {
      for(term = 0; term < nterms; term++)
        gram[term + (size_t)nterms*(size_t)term] += epsilon;
      nepsilon += epsilon;
      if(check_interrupt)
        R_CheckUserInterrupt();
      if(np_reghat_solve_system(nterms, gram, rhs, gram_work, solution,
                                pivot, condition_work, condition_iwork)) {
        solved = 1;
        break;
      }
    }
    if(!solved)
      return NP_REGHAT_LP_ROW_RIDGE_FAILED;

    {
      double denom = gram[0];

      if(!isfinite(denom) || (fabs(denom) < DBL_MIN))
        denom = DBL_MIN;
      solution[0] *= 1.0 + nepsilon/denom;
    }
  }

  F77_CALL(dgemv)(&trans_n, &ntrain, &nterms, &alpha,
                  design, &ntrain, solution, &one,
                  &beta, prediction, &one FCONE);
  return NP_REGHAT_LP_ROW_OK;
}

void np_reghat_lp_workspace_init(NPReghatLPWorkspace *workspace)
{
  if(workspace != NULL)
    memset(workspace, 0, sizeof(*workspace));
}

void np_reghat_lp_workspace_clear(NPReghatLPWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  free(workspace->design);
  free(workspace->weighted_design);
  free(workspace->gram);
  free(workspace->gram_work);
  free(workspace->rhs);
  free(workspace->solution);
  free(workspace->prediction);
  free(workspace->condition_work);
  free(workspace->pivot);
  free(workspace->condition_iwork);
  np_reghat_lp_workspace_init(workspace);
}

static NPReghatLPRowStatus np_reghat_lp_workspace_reserve(
  NPReghatLPWorkspace *workspace,
  int ntrain,
  int nterms)
{
  size_t design_elements;
  size_t gram_elements;
  double *design = NULL;
  double *weighted_design = NULL;
  double *gram = NULL;
  double *gram_work = NULL;
  double *rhs = NULL;
  double *solution = NULL;
  double *prediction = NULL;
  double *condition_work = NULL;
  int *pivot = NULL;
  int *condition_iwork = NULL;

  if((workspace == NULL) || (ntrain <= 0) || (nterms <= 1) ||
     ((size_t)ntrain > SIZE_MAX/(size_t)nterms) ||
     ((size_t)nterms > SIZE_MAX/(size_t)nterms))
    return NP_REGHAT_LP_ROW_INVALID;
  design_elements = (size_t)ntrain*(size_t)nterms;
  gram_elements = (size_t)nterms*(size_t)nterms;
  if((design_elements > SIZE_MAX/sizeof(double)) ||
     (gram_elements > SIZE_MAX/sizeof(double)) ||
     ((size_t)nterms > SIZE_MAX/(4U*sizeof(double))) ||
     ((size_t)nterms > SIZE_MAX/sizeof(int)) ||
     ((size_t)ntrain > SIZE_MAX/sizeof(double)))
    return NP_REGHAT_LP_ROW_INVALID;
  if((workspace->ntrain == ntrain) && (workspace->nterms == nterms) &&
     (workspace->design_capacity >= design_elements) &&
     (workspace->gram_capacity >= gram_elements) &&
     workspace->design != NULL && workspace->weighted_design != NULL &&
     workspace->gram != NULL && workspace->gram_work != NULL &&
     workspace->rhs != NULL && workspace->solution != NULL &&
     workspace->prediction != NULL && workspace->condition_work != NULL &&
     workspace->pivot != NULL && workspace->condition_iwork != NULL)
    return NP_REGHAT_LP_ROW_OK;

  design = (double *)malloc(design_elements*sizeof(double));
  weighted_design = (double *)malloc(design_elements*sizeof(double));
  gram = (double *)malloc(gram_elements*sizeof(double));
  gram_work = (double *)malloc(gram_elements*sizeof(double));
  rhs = (double *)malloc((size_t)nterms*sizeof(double));
  solution = (double *)malloc((size_t)nterms*sizeof(double));
  prediction = (double *)malloc((size_t)ntrain*sizeof(double));
  condition_work = (double *)malloc(4U*(size_t)nterms*sizeof(double));
  pivot = (int *)malloc((size_t)nterms*sizeof(int));
  condition_iwork = (int *)malloc((size_t)nterms*sizeof(int));
  if(design == NULL || weighted_design == NULL || gram == NULL ||
     gram_work == NULL || rhs == NULL || solution == NULL ||
     prediction == NULL || condition_work == NULL || pivot == NULL ||
     condition_iwork == NULL) {
    free(design);
    free(weighted_design);
    free(gram);
    free(gram_work);
    free(rhs);
    free(solution);
    free(prediction);
    free(condition_work);
    free(pivot);
    free(condition_iwork);
    return NP_REGHAT_LP_ROW_MEMORY;
  }

  np_reghat_lp_workspace_clear(workspace);
  workspace->ntrain = ntrain;
  workspace->nterms = nterms;
  workspace->design_capacity = design_elements;
  workspace->gram_capacity = gram_elements;
  workspace->design = design;
  workspace->weighted_design = weighted_design;
  workspace->gram = gram;
  workspace->gram_work = gram_work;
  workspace->rhs = rhs;
  workspace->solution = solution;
  workspace->prediction = prediction;
  workspace->condition_work = condition_work;
  workspace->pivot = pivot;
  workspace->condition_iwork = condition_iwork;
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
    workspace->gram, workspace->gram_work, workspace->rhs,
    workspace->solution, workspace->prediction, workspace->condition_work,
    workspace->pivot, workspace->condition_iwork, 0);
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

SEXP C_np_reghat_lp_matrix_fast(SEXP kw, SEXP wtrain, SEXP weval)
{
  int ntrain = 0, neval = 0, kw_neval = 0;
  int wtrain_n = 0, nterms = 0, weval_n = 0, weval_p = 0;
  double *weighted_design = NULL;
  double *gram = NULL;
  double *gram_work = NULL;
  double *rhs = NULL;
  double *solution = NULL;
  double *prediction = NULL;
  double *condition_work = NULL;
  int *pivot = NULL;
  int *condition_iwork = NULL;
  SEXP out = R_NilValue;

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

  if(((size_t)nterms > SIZE_MAX/sizeof(double)) ||
     ((size_t)ntrain > SIZE_MAX/((size_t)nterms*sizeof(double))) ||
     ((size_t)nterms > SIZE_MAX/((size_t)nterms*sizeof(double))) ||
     ((size_t)nterms > SIZE_MAX/(4*sizeof(double))) ||
     ((size_t)neval > SIZE_MAX/(size_t)ntrain))
    return R_NilValue;

  weighted_design = (double *)R_alloc((size_t)ntrain*(size_t)nterms,
                                      sizeof(double));
  gram = (double *)R_alloc((size_t)nterms*(size_t)nterms, sizeof(double));
  gram_work = (double *)R_alloc((size_t)nterms*(size_t)nterms, sizeof(double));
  rhs = (double *)R_alloc((size_t)nterms, sizeof(double));
  solution = (double *)R_alloc((size_t)nterms, sizeof(double));
  prediction = (double *)R_alloc((size_t)ntrain, sizeof(double));
  condition_work = (double *)R_alloc((size_t)4*(size_t)nterms, sizeof(double));
  pivot = (int *)R_alloc((size_t)nterms, sizeof(int));
  condition_iwork = (int *)R_alloc((size_t)nterms, sizeof(int));

  out = PROTECT(allocMatrix(REALSXP, neval, ntrain));
  for(int j = 0; j < neval; j++){
    const double * const weights = REAL(kw) + (size_t)j*(size_t)ntrain;
    NPReghatLPRowStatus row_status;

    if((j == 0) || (j + 1 == neval) || ((j % 32) == 0))
      R_CheckUserInterrupt();

    for(int term = 0; term < nterms; term++)
      rhs[term] = REAL(weval)[j + (size_t)neval*(size_t)term];
    row_status = np_reghat_lp_prediction_raw(
      ntrain, nterms, REAL(wtrain), weights, rhs, weighted_design, gram,
      gram_work, rhs, solution, prediction, condition_work, pivot,
      condition_iwork, 1);
    if(row_status == NP_REGHAT_LP_ROW_NONFINITE)
      error("LP solve failed in compiled hat-matrix path: non-finite system");
    if(row_status == NP_REGHAT_LP_ROW_RIDGE_FAILED)
      error("LP solve failed in compiled hat-matrix path after bounded ridging");
    if(row_status != NP_REGHAT_LP_ROW_OK)
      error("invalid wider-LP compiled hat-matrix input");
    for(int i = 0; i < ntrain; i++)
      REAL(out)[j + (size_t)neval*(size_t)i] = weights[i]*prediction[i];
  }

  UNPROTECT(1);
  return out;
}
