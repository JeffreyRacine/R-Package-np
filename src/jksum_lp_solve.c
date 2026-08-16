/* Reusable contiguous LAPACK workspace for canonical local-polynomial solves. */

#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <R_ext/Arith.h>
#include <R_ext/Lapack.h>

#include "jksum_lp_solve.h"

static int np_lp_size_product(size_t a, size_t b, size_t *result)
{
  if((result == NULL) || ((b != 0U) && (a > SIZE_MAX/b)))
    return 0;
  *result = a*b;
  return 1;
}

static int np_lp_double_bytes(size_t elements, size_t *bytes)
{
  return np_lp_size_product(elements, sizeof(double), bytes);
}

static int np_lp_solve_workspace_shape(
  const NPLPSolveWorkspace *workspace,
  int p,
  int nrhs,
  size_t *gram_elements,
  size_t *rhs_elements)
{
  size_t gram_count, rhs_count;

  if((workspace == NULL) || (p <= 0) || (nrhs <= 0) ||
     (workspace->p_capacity < p) || (workspace->nrhs_capacity < nrhs) ||
     (workspace->gram_source == NULL) || (workspace->rhs_source == NULL) ||
     (workspace->gram_work == NULL) || (workspace->rhs_work == NULL) ||
     (workspace->ipiv == NULL) ||
     !np_lp_size_product((size_t)p, (size_t)p, &gram_count) ||
     !np_lp_size_product((size_t)p, (size_t)nrhs, &rhs_count) ||
     (gram_count > workspace->gram_capacity) ||
     (rhs_count > workspace->rhs_capacity))
    return 0;

  if(gram_elements != NULL)
    *gram_elements = gram_count;
  if(rhs_elements != NULL)
    *rhs_elements = rhs_count;
  return 1;
}

/*
 * A one-term local-polynomial system is an exact scalar problem.  Keep this
 * specialization inside the shared solve workspace so all of its width-one
 * consumers use the same algebra and cannot fall through to LAPACK or a
 * legacy LC evaluator.  Widths greater than one retain the existing LAPACK
 * transcript.
 */
static inline int np_lp_solve_width_one(const double gram,
                                        const double *rhs_source,
                                        double *rhs_work,
                                        int nrhs)
{
  int i;

  if((rhs_source == NULL) || (rhs_work == NULL) || (nrhs <= 0))
    return 0;

  for(i = 0; i < nrhs; i++){
    const double solution = rhs_source[i]/gram;
    if(!R_FINITE(solution))
      return 0;
    rhs_work[i] = solution;
  }

  return 1;
}

int np_lp_solve_workspace_sources_finite(
  const NPLPSolveWorkspace *workspace,
  int p,
  int nrhs)
{
  size_t gram_elements, rhs_elements, i;

  if(!np_lp_solve_workspace_shape(workspace,
                                  p,
                                  nrhs,
                                  &gram_elements,
                                  &rhs_elements))
    return 0;
  for(i = 0; i < gram_elements; i++)
    if(!R_FINITE(workspace->gram_source[i]))
      return 0;
  for(i = 0; i < rhs_elements; i++)
    if(!R_FINITE(workspace->rhs_source[i]))
      return 0;
  return 1;
}

void np_lp_solve_workspace_init(NPLPSolveWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  memset(workspace, 0, sizeof(*workspace));
}

void np_lp_solve_workspace_clear(NPLPSolveWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  free(workspace->gram_source);
  free(workspace->rhs_source);
  free(workspace->gram_work);
  free(workspace->rhs_work);
  free(workspace->ipiv);
  np_lp_solve_workspace_init(workspace);
}

int np_lp_solve_workspace_reserve(NPLPSolveWorkspace *workspace,
                                  int p,
                                  int nrhs)
{
  size_t gram_elements, rhs_elements, gram_bytes, rhs_bytes, pivot_bytes;
  double *gram_source = NULL, *rhs_source = NULL;
  double *gram_work = NULL, *rhs_work = NULL;
  int *ipiv = NULL;

  if((workspace == NULL) || (p <= 0) || (nrhs <= 0))
    return 0;
  workspace->factor_ready = 0;
  workspace->factor_p = 0;
  if((workspace->p_capacity >= p) && (workspace->nrhs_capacity >= nrhs))
    return 1;
  if(!np_lp_size_product((size_t)p, (size_t)p, &gram_elements) ||
     !np_lp_size_product((size_t)p, (size_t)nrhs, &rhs_elements) ||
     !np_lp_double_bytes(gram_elements, &gram_bytes) ||
     !np_lp_double_bytes(rhs_elements, &rhs_bytes) ||
     !np_lp_size_product((size_t)p, sizeof(int), &pivot_bytes))
    return 0;

  gram_source = (double *)malloc(gram_bytes);
  rhs_source = (double *)malloc(rhs_bytes);
  gram_work = (double *)malloc(gram_bytes);
  rhs_work = (double *)malloc(rhs_bytes);
  ipiv = (int *)malloc(pivot_bytes);
  if((gram_source == NULL) || (rhs_source == NULL) ||
     (gram_work == NULL) || (rhs_work == NULL) || (ipiv == NULL)){
    free(gram_source);
    free(rhs_source);
    free(gram_work);
    free(rhs_work);
    free(ipiv);
    return 0;
  }

  np_lp_solve_workspace_clear(workspace);
  workspace->p_capacity = p;
  workspace->nrhs_capacity = nrhs;
  workspace->gram_capacity = gram_elements;
  workspace->rhs_capacity = rhs_elements;
  workspace->gram_source = gram_source;
  workspace->rhs_source = rhs_source;
  workspace->gram_work = gram_work;
  workspace->rhs_work = rhs_work;
  workspace->ipiv = ipiv;
  return 1;
}

int np_lp_solve_workspace_solve(NPLPSolveWorkspace *workspace,
                                int p,
                                int nrhs)
{
  size_t gram_elements, rhs_elements;
  int info = 0;
  size_t i;

  if(workspace != NULL){
    workspace->factor_ready = 0;
    workspace->factor_p = 0;
  }
  if((workspace == NULL) || (p <= 0) || (nrhs <= 0) ||
     (workspace->p_capacity < p) || (workspace->nrhs_capacity < nrhs) ||
     (workspace->gram_source == NULL) || (workspace->rhs_source == NULL) ||
     (workspace->gram_work == NULL) || (workspace->rhs_work == NULL) ||
     (workspace->ipiv == NULL) ||
     !np_lp_size_product((size_t)p, (size_t)p, &gram_elements) ||
     !np_lp_size_product((size_t)p, (size_t)nrhs, &rhs_elements) ||
     (gram_elements > workspace->gram_capacity) ||
     (rhs_elements > workspace->rhs_capacity))
    return 0;

  if(p == 1){
    workspace->gram_work[0] = workspace->gram_source[0];
    if(!np_lp_solve_width_one(workspace->gram_source[0],
                              workspace->rhs_source,
                              workspace->rhs_work,
                              nrhs))
      return 0;
    workspace->factor_ready = 1;
    workspace->factor_p = 1;
    return 1;
  }

  memcpy(workspace->gram_work,
         workspace->gram_source,
         gram_elements*sizeof(double));
  memcpy(workspace->rhs_work,
         workspace->rhs_source,
         rhs_elements*sizeof(double));
  F77_CALL(dgesv)(&p, &nrhs,
                  workspace->gram_work, &p,
                  workspace->ipiv,
                  workspace->rhs_work, &p,
                  &info);
  if(info != 0)
    return 0;
  for(i = 0; i < rhs_elements; i++)
    if(!R_FINITE(workspace->rhs_work[i]))
      return 0;
  workspace->factor_ready = 1;
  workspace->factor_p = p;
  return 1;
}

static int np_lp_solve_workspace_solve_factored_with_trans(
  NPLPSolveWorkspace *workspace,
  int p,
  int nrhs,
  const char trans)
{
  size_t rhs_elements;
  int info = 0;
  size_t i;

  if((trans != 'N' && trans != 'T') ||
     (workspace == NULL) || (p <= 0) || (nrhs <= 0) ||
     (!workspace->factor_ready) || (workspace->factor_p != p) ||
     (workspace->p_capacity < p) || (workspace->nrhs_capacity < nrhs) ||
     (workspace->rhs_source == NULL) || (workspace->gram_work == NULL) ||
     (workspace->rhs_work == NULL) || (workspace->ipiv == NULL) ||
     !np_lp_size_product((size_t)p, (size_t)nrhs, &rhs_elements) ||
     (rhs_elements > workspace->rhs_capacity))
    return 0;

  if(p == 1)
    return np_lp_solve_width_one(workspace->gram_work[0],
                                 workspace->rhs_source,
                                 workspace->rhs_work,
                                 nrhs);

  memcpy(workspace->rhs_work,
         workspace->rhs_source,
         rhs_elements*sizeof(double));
  F77_CALL(dgetrs)(&trans, &p, &nrhs,
                   workspace->gram_work, &p,
                   workspace->ipiv,
                   workspace->rhs_work, &p,
                   &info
                   FCONE);
  if(info != 0)
    return 0;
  for(i = 0; i < rhs_elements; i++)
    if(!R_FINITE(workspace->rhs_work[i]))
      return 0;
  return 1;
}

int np_lp_solve_workspace_solve_factored(NPLPSolveWorkspace *workspace,
                                         int p,
                                         int nrhs)
{
  return np_lp_solve_workspace_solve_factored_with_trans(
    workspace, p, nrhs, 'N');
}

/*
 * Select the factor/ridge state from the Gram system alone.  Response shape
 * and magnitude must not choose a different statistical regularization.  The
 * retained LU is consumed by both response and adjoint solves.
 */
static int np_lp_solve_workspace_factor(NPLPSolveWorkspace *workspace,
                                        int p)
{
  size_t gram_elements, i;
  int info = 0;

  if(workspace != NULL){
    workspace->factor_ready = 0;
    workspace->factor_p = 0;
  }
  if(!np_lp_solve_workspace_shape(workspace, p, 1,
                                  &gram_elements, NULL))
    return 0;

  memcpy(workspace->gram_work,
         workspace->gram_source,
         gram_elements*sizeof(double));
  if(p == 1){
    if(!R_FINITE(workspace->gram_work[0]) ||
       (workspace->gram_work[0] == 0.0))
      return 0;
    workspace->factor_ready = 1;
    workspace->factor_p = 1;
    return 1;
  }

  F77_CALL(dgetrf)(&p, &p,
                   workspace->gram_work, &p,
                   workspace->ipiv,
                   &info);
  if(info != 0)
    return 0;
  for(i = 0; i < gram_elements; i++)
    if(!R_FINITE(workspace->gram_work[i]))
      return 0;
  workspace->factor_ready = 1;
  workspace->factor_p = p;
  return 1;
}

static NPLPSolvePolicyStatus np_lp_solve_workspace_prepare_policy_factor(
  NPLPSolveWorkspace *workspace,
  int p,
  int nrhs,
  double ridge_increment,
  NPLPSolvePolicyDiagnostics *diagnostics)
{
  double ridge_total = 0.0;
  int ridge_steps = 0;
  int i;

  if(diagnostics != NULL){
    diagnostics->ridge_steps = 0;
    diagnostics->ridge_total = 0.0;
  }

  if(!np_lp_solve_workspace_shape(workspace, p, nrhs, NULL, NULL) ||
     !R_FINITE(ridge_increment) || (ridge_increment <= 0.0))
    return NP_LP_SOLVE_POLICY_INVALID;
  if(!np_lp_solve_workspace_sources_finite(workspace, p, nrhs))
    return NP_LP_SOLVE_POLICY_NONFINITE;

  while(!np_lp_solve_workspace_factor(workspace, p)){
    if(ridge_steps >= NP_LP_SOLVE_MAX_RIDGE_STEPS)
      return NP_LP_SOLVE_POLICY_RIDGE_EXHAUSTED;
    if(!np_lp_solve_workspace_sources_finite(workspace, p, nrhs))
      return NP_LP_SOLVE_POLICY_NONFINITE;

    for(i = 0; i < p; i++)
      workspace->gram_source[i + i*p] += ridge_increment;
    ridge_total += ridge_increment;
    ridge_steps++;
    if(diagnostics != NULL){
      diagnostics->ridge_steps = ridge_steps;
      diagnostics->ridge_total = ridge_total;
    }
  }

  return NP_LP_SOLVE_POLICY_OK;
}

NPLPSolvePolicyStatus np_lp_solve_workspace_solve_response(
  NPLPSolveWorkspace *workspace,
  int p,
  int nrhs,
  double ridge_increment,
  NPLPSolvePolicyDiagnostics *diagnostics)
{
  NPLPSolvePolicyDiagnostics local_diagnostics;
  NPLPSolvePolicyDiagnostics * const policy_diagnostics =
    (diagnostics != NULL) ? diagnostics : &local_diagnostics;
  const NPLPSolvePolicyStatus factor_status =
    np_lp_solve_workspace_prepare_policy_factor(
      workspace, p, nrhs, ridge_increment, policy_diagnostics);

  if(factor_status != NP_LP_SOLVE_POLICY_OK)
    return factor_status;

  if(policy_diagnostics->ridge_total > 0.0){
    const double denominator =
      (workspace->gram_source[0] > DBL_EPSILON) ?
      workspace->gram_source[0] : DBL_EPSILON;
    int rhs;

    for(rhs = 0; rhs < nrhs; rhs++){
      double * const intercept = workspace->rhs_source + (size_t)rhs*(size_t)p;
      *intercept += policy_diagnostics->ridge_total*(*intercept)/denominator;
    }
  }

  if(!np_lp_solve_workspace_solve_factored(workspace, p, nrhs)){
    if(!np_lp_solve_workspace_sources_finite(workspace, p, nrhs))
      return NP_LP_SOLVE_POLICY_NONFINITE;
    return NP_LP_SOLVE_POLICY_FINAL_FAILED;
  }

  return NP_LP_SOLVE_POLICY_OK;
}

NPLPSolvePolicyStatus np_lp_solve_workspace_solve_adjoint(
  NPLPSolveWorkspace *workspace,
  int p,
  int nrhs,
  double ridge_increment,
  NPLPSolvePolicyDiagnostics *diagnostics)
{
  NPLPSolvePolicyDiagnostics local_diagnostics;
  NPLPSolvePolicyDiagnostics * const policy_diagnostics =
    (diagnostics != NULL) ? diagnostics : &local_diagnostics;
  const NPLPSolvePolicyStatus factor_status =
    np_lp_solve_workspace_prepare_policy_factor(
      workspace, p, nrhs, ridge_increment, policy_diagnostics);

  if(factor_status != NP_LP_SOLVE_POLICY_OK)
    return factor_status;
  if(!np_lp_solve_workspace_solve_factored_with_trans(
       workspace, p, nrhs, 'T')){
    if(!np_lp_solve_workspace_sources_finite(workspace, p, nrhs))
      return NP_LP_SOLVE_POLICY_NONFINITE;
    return NP_LP_SOLVE_POLICY_FINAL_FAILED;
  }

  if(policy_diagnostics->ridge_total > 0.0){
    const double denominator =
      (workspace->gram_source[0] > DBL_EPSILON) ?
      workspace->gram_source[0] : DBL_EPSILON;
    const double intercept_scale =
      1.0 + policy_diagnostics->ridge_total/denominator;
    int rhs;

    for(rhs = 0; rhs < nrhs; rhs++){
      double * const intercept = workspace->rhs_work + (size_t)rhs*(size_t)p;
      *intercept *= intercept_scale;
      if(!R_FINITE(*intercept))
        return NP_LP_SOLVE_POLICY_FINAL_FAILED;
    }
  }

  return NP_LP_SOLVE_POLICY_OK;
}

/*
 * A one-column weighted design has the exact influence row
 *
 *   w_i z_i z_eval / sum_j(w_j z_j^2).
 *
 * Compute it directly in training-row order.  This covers the implicit unit
 * basis used by LP0 without assuming that the sole basis column is constant,
 * retains signed higher-order kernel weights, and prevents width one from
 * allocating solver storage or entering LAPACK.
 */
NPLPWidthOneStatus np_lp_width_one_influence_row(
  const double *basis_train,
  int n,
  const double *kw,
  double basis_eval,
  double *row_out,
  size_t output_stride)
{
  double denominator = 0.0;
  double projection;
  double ridge_total = 0.0;
  int i;

  if((basis_train == NULL) || (kw == NULL) || (row_out == NULL) ||
     (n <= 0) || (output_stride == 0U))
    return NP_LP_WIDTH_ONE_INVALID;

  for(i = 0; i < n; i++){
    const double zi = basis_train[i];
    denominator += kw[i]*zi*zi;
  }

  if(!R_FINITE(denominator) || !R_FINITE(basis_eval))
    return NP_LP_WIDTH_ONE_NONFINITE;

  projection = basis_eval/denominator;
  if(!R_FINITE(projection)){
    const double ridge_increment = 1.0/(double)n;
    int solved = 0;

    for(i = 0; i < NP_LP_SOLVE_MAX_RIDGE_STEPS; i++){
      denominator += ridge_increment;
      ridge_total += ridge_increment;
      projection = basis_eval/denominator;
      if(R_FINITE(projection)){
        solved = 1;
        break;
      }
    }
    if(!solved)
      return NP_LP_WIDTH_ONE_RIDGE_FAILED;

    {
      double correction_denominator = denominator;
      if(fabs(correction_denominator) < DBL_MIN)
        correction_denominator = DBL_MIN;
      projection *= 1.0 + ridge_total/correction_denominator;
    }
    if(!R_FINITE(projection))
      return NP_LP_WIDTH_ONE_RIDGE_FAILED;
  }

  for(i = 0; i < n; i++){
    const double value = kw[i]*basis_train[i]*projection;
    if(!R_FINITE(value))
      return NP_LP_WIDTH_ONE_NONFINITE;
    row_out[(size_t)i*output_stride] = value;
  }

  return NP_LP_WIDTH_ONE_OK;
}

void np_lp_full_row_workspace_init(NPLPFullRowWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  memset(workspace, 0, sizeof(*workspace));
}

void np_lp_full_row_workspace_clear(NPLPFullRowWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  free(workspace->gram);
  free(workspace->rhs);
  free(workspace->ipiv);
  free(workspace->matrix_copy);
  free(workspace->rcond_values);
  free(workspace->rcond_work);
  free(workspace->inverse_work);
  np_lp_full_row_workspace_init(workspace);
}

int np_lp_full_row_workspace_reserve(NPLPFullRowWorkspace *workspace,
                                     int p,
                                     int nrhs)
{
  size_t gram_elements, rhs_elements;
  size_t gram_bytes, rhs_bytes, pivot_bytes, values_bytes;
  double *gram = NULL, *rhs = NULL;
  int *ipiv = NULL;
  double *matrix_copy = NULL, *rcond_values = NULL;

  if((workspace == NULL) || (p <= 0) || (nrhs <= 0))
    return 0;
  if((workspace->p_capacity >= p) &&
     (workspace->nrhs_capacity >= nrhs) &&
     (workspace->gram != NULL) &&
     (workspace->rhs != NULL) &&
     (workspace->ipiv != NULL) &&
     (workspace->matrix_copy != NULL) &&
     (workspace->rcond_values != NULL))
    return 1;
  if(!np_lp_size_product((size_t)p, (size_t)p, &gram_elements) ||
     !np_lp_size_product((size_t)p, (size_t)nrhs, &rhs_elements) ||
     !np_lp_double_bytes(gram_elements, &gram_bytes) ||
     !np_lp_double_bytes(rhs_elements, &rhs_bytes) ||
     !np_lp_size_product((size_t)p, sizeof(int), &pivot_bytes) ||
     !np_lp_double_bytes((size_t)p, &values_bytes))
    return 0;

  gram = (double *)malloc(gram_bytes);
  rhs = (double *)malloc(rhs_bytes);
  ipiv = (int *)malloc(pivot_bytes);
  matrix_copy = (double *)malloc(gram_bytes);
  rcond_values = (double *)malloc(values_bytes);
  if((gram == NULL) || (rhs == NULL) || (ipiv == NULL) ||
     (matrix_copy == NULL) || (rcond_values == NULL)){
    free(gram);
    free(rhs);
    free(ipiv);
    free(matrix_copy);
    free(rcond_values);
    return 0;
  }

  np_lp_full_row_workspace_clear(workspace);
  workspace->p_capacity = p;
  workspace->nrhs_capacity = nrhs;
  workspace->gram_capacity = gram_elements;
  workspace->rhs_capacity = rhs_elements;
  workspace->gram = gram;
  workspace->rhs = rhs;
  workspace->ipiv = ipiv;
  workspace->matrix_copy = matrix_copy;
  workspace->rcond_values = rcond_values;
  return 1;
}

static int np_lp_full_row_bad_rcond(NPLPFullRowWorkspace *workspace,
                                    int p,
                                    double min_rcond)
{
  char jobz = 'N';
  char uplo = 'U';
  int i, j;
  int info = 0;
  double max_eval = 0.0, min_eval = DBL_MAX;

  if((workspace == NULL) || (p <= 0) ||
     (workspace->p_capacity < p) ||
     (workspace->matrix_copy == NULL) ||
     (workspace->rcond_values == NULL) ||
     (workspace->gram == NULL))
    return 1;

  if(p == 1){
    const double gram = workspace->gram[0];
    return (!R_FINITE(gram) || !(fabs(gram) > 0.0) ||
            (1.0 < min_rcond));
  }

  for(j = 0; j < p; j++)
    for(i = 0; i < p; i++)
      workspace->matrix_copy[i + j*p] =
        workspace->gram[i + j*p];

  if((workspace->rcond_work == NULL) ||
     (workspace->rcond_lwork_capacity <= 0)){
    int lwork_query = -1;
    double work_query = 0.0;
    int requested_lwork;
    double *rcond_work;

    F77_CALL(dsyev)(&jobz, &uplo, &p,
                    workspace->matrix_copy, &p,
                    workspace->rcond_values,
                    &work_query, &lwork_query, &info FCONE FCONE);
    if(info != 0)
      return 1;
    requested_lwork = ((int)work_query > 1) ? (int)work_query : 1;
    rcond_work = (double *)malloc((size_t)requested_lwork*sizeof(double));
    if(rcond_work == NULL)
      return 1;
    workspace->rcond_work = rcond_work;
    workspace->rcond_lwork_capacity = requested_lwork;
  }

  F77_CALL(dsyev)(&jobz, &uplo, &p,
                  workspace->matrix_copy, &p,
                  workspace->rcond_values,
                  workspace->rcond_work,
                  &workspace->rcond_lwork_capacity, &info FCONE FCONE);
  if(info != 0)
    return 1;

  for(i = 0; i < p; i++){
    if(!R_FINITE(workspace->rcond_values[i]))
      return 1;
    {
      const double abs_eval = fabs(workspace->rcond_values[i]);
      if(abs_eval > max_eval)
        max_eval = abs_eval;
      if(abs_eval < min_eval)
        min_eval = abs_eval;
    }
  }
  if(!(max_eval > 0.0))
    return 1;
  return ((min_eval / max_eval) < min_rcond);
}

int np_lp_full_row_workspace_solve(NPLPFullRowWorkspace *workspace,
                                   int p,
                                   int nrhs,
                                   double min_rcond)
{
  size_t rhs_elements;
  size_t i;
  int info = 0;

  if(!np_lp_full_row_workspace_reserve(workspace, p, nrhs) ||
     np_lp_full_row_bad_rcond(workspace, p, min_rcond) ||
     !np_lp_size_product((size_t)p, (size_t)nrhs, &rhs_elements) ||
     (rhs_elements > workspace->rhs_capacity))
    return 0;
  if(p == 1)
    return np_lp_solve_width_one(workspace->gram[0],
                                 workspace->rhs,
                                 workspace->rhs,
                                 nrhs);
  F77_CALL(dgesv)(&p, &nrhs,
                  workspace->gram, &p,
                  workspace->ipiv,
                  workspace->rhs, &p,
                  &info);
  if(info != 0)
    return 0;
  for(i = 0; i < rhs_elements; i++)
    if(!R_FINITE(workspace->rhs[i]))
      return 0;
  return 1;
}

static int np_lp_full_row_workspace_factor_invert(
  NPLPFullRowWorkspace *workspace,
  int p)
{
  size_t gram_elements;
  size_t i;
  int info = 0;

  if((workspace == NULL) || (p <= 0) ||
     (workspace->p_capacity < p) ||
     (workspace->gram == NULL) ||
     (workspace->ipiv == NULL) ||
     !np_lp_size_product((size_t)p, (size_t)p, &gram_elements) ||
     (gram_elements > workspace->gram_capacity))
    return 0;

  F77_CALL(dgetrf)(&p, &p,
                   workspace->gram, &p,
                   workspace->ipiv,
                   &info);
  if(info != 0)
    return 0;

  if((workspace->inverse_work == NULL) ||
     (workspace->inverse_lwork_capacity <= 0)){
    int lwork_query = -1;
    double work_query = 0.0;
    int requested_lwork;
    double *inverse_work;

    F77_CALL(dgetri)(&p,
                     workspace->gram, &p,
                     workspace->ipiv,
                     &work_query, &lwork_query,
                     &info);
    if(info != 0)
      return 0;
    requested_lwork = ((int)work_query > p) ? (int)work_query : p;
    inverse_work = (double *)malloc((size_t)requested_lwork*sizeof(double));
    if(inverse_work == NULL)
      return 0;
    workspace->inverse_work = inverse_work;
    workspace->inverse_lwork_capacity = requested_lwork;
  }

  F77_CALL(dgetri)(&p,
                   workspace->gram, &p,
                   workspace->ipiv,
                   workspace->inverse_work,
                   &workspace->inverse_lwork_capacity,
                   &info);
  if(info != 0)
    return 0;
  for(i = 0; i < gram_elements; i++)
    if(!R_FINITE(workspace->gram[i]))
      return 0;
  return 1;
}

int np_lp_full_row_workspace_invert(NPLPFullRowWorkspace *workspace,
                                    int p,
                                    double min_rcond)
{
  if(!np_lp_full_row_workspace_reserve(workspace, p, 1) ||
     np_lp_full_row_bad_rcond(workspace, p, min_rcond))
    return 0;
  if(p == 1){
    const double inverse = 1.0/workspace->gram[0];
    if(!R_FINITE(inverse))
      return 0;
    workspace->gram[0] = inverse;
    return 1;
  }
  return np_lp_full_row_workspace_factor_invert(workspace, p);
}

int np_lp_full_row_workspace_invert_retryable(
  NPLPFullRowWorkspace *workspace,
  int p,
  double ridge_increment,
  int max_ridge_steps)
{
  size_t gram_elements;
  int attempt, i;

  if(!np_lp_full_row_workspace_reserve(workspace, p, 1) ||
     (workspace->matrix_copy == NULL) ||
     (!R_FINITE(ridge_increment)) ||
     (ridge_increment < 0.0) ||
     (max_ridge_steps < 0) ||
     !np_lp_size_product((size_t)p, (size_t)p, &gram_elements) ||
     (gram_elements > workspace->gram_capacity))
    return 0;

  for(attempt = 0; attempt <= max_ridge_steps; attempt++){
    memcpy(workspace->gram,
           workspace->matrix_copy,
           gram_elements*sizeof(double));
    if(p == 1){
      const double inverse = 1.0/workspace->gram[0];
      if(R_FINITE(inverse)){
        workspace->gram[0] = inverse;
        return 1;
      }
    } else {
      if(np_lp_full_row_workspace_factor_invert(workspace, p))
        return 1;
    }
    for(i = 0; i < p; i++)
      workspace->matrix_copy[i + i*p] += ridge_increment;
  }

  return 0;
}

int np_lp_full_row_workspace_pack_inverse_rows(
  NPLPFullRowWorkspace *workspace,
  int p)
{
  int i, j;

  if((workspace == NULL) || (p <= 0) ||
     (workspace->p_capacity < p) ||
     (workspace->gram == NULL) ||
     (workspace->matrix_copy == NULL))
    return 0;

  for(i = 0; i < p; i++)
    for(j = 0; j < p; j++)
      workspace->matrix_copy[i*p + j] =
        workspace->gram[i + j*p];

  return 1;
}
