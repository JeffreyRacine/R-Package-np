/* Reusable contiguous LAPACK workspace for canonical local-polynomial solves. */

#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <R_ext/Arith.h>
#include <R_ext/Applic.h>
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
  return 1;
}

void np_glp_qr_drop_workspace_init(NPGLPQRDropWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  memset(workspace, 0, sizeof(*workspace));
}

void np_glp_qr_drop_workspace_clear(NPGLPQRDropWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  free(workspace->xqr);
  free(workspace->qraux);
  free(workspace->work);
  free(workspace->y);
  free(workspace->qy);
  free(workspace->pivot);
  np_glp_qr_drop_workspace_init(workspace);
}

int np_glp_qr_drop_workspace_reserve(NPGLPQRDropWorkspace *workspace,
                                     int n,
                                     int p)
{
  size_t xqr_elements, xqr_bytes, n_bytes, p_bytes;
  size_t work_elements, work_bytes, pivot_bytes;
  double *xqr = NULL, *qraux = NULL, *work = NULL;
  double *y = NULL, *qy = NULL;
  int *pivot = NULL;

  if((workspace == NULL) || (n <= 0) || (p <= 0))
    return 0;
  if((workspace->n_capacity >= n) && (workspace->p_capacity >= p))
    return 1;
  if(!np_lp_size_product((size_t)n, (size_t)p, &xqr_elements) ||
     !np_lp_double_bytes(xqr_elements, &xqr_bytes) ||
     !np_lp_double_bytes((size_t)n, &n_bytes) ||
     !np_lp_double_bytes((size_t)p, &p_bytes) ||
     !np_lp_size_product((size_t)2, (size_t)p, &work_elements) ||
     !np_lp_double_bytes(work_elements, &work_bytes) ||
     !np_lp_size_product((size_t)p, sizeof(int), &pivot_bytes))
    return 0;

  xqr = (double *)malloc(xqr_bytes);
  qraux = (double *)malloc(p_bytes);
  work = (double *)malloc(work_bytes);
  y = (double *)malloc(n_bytes);
  qy = (double *)malloc(n_bytes);
  pivot = (int *)malloc(pivot_bytes);
  if((xqr == NULL) || (qraux == NULL) || (work == NULL) ||
     (y == NULL) || (qy == NULL) || (pivot == NULL)){
    free(xqr);
    free(qraux);
    free(work);
    free(y);
    free(qy);
    free(pivot);
    return 0;
  }

  np_glp_qr_drop_workspace_clear(workspace);
  workspace->n_capacity = n;
  workspace->p_capacity = p;
  workspace->xqr_capacity = xqr_elements;
  workspace->xqr = xqr;
  workspace->qraux = qraux;
  workspace->work = work;
  workspace->y = y;
  workspace->qy = qy;
  workspace->pivot = pivot;
  return 1;
}

int np_glp_qr_drop_workspace_apply(NPGLPQRDropWorkspace *workspace,
                                   double **basis,
                                   int n,
                                   int p,
                                   const double *kw,
                                   int eval_pos,
                                   double *row_out)
{
  const double tol = 1.0e-7;
  int ldx = n, rank = 0, ny = 1;
  int i, j;

  if((workspace == NULL) || (basis == NULL) || (kw == NULL) ||
     (row_out == NULL) || (n <= 0) || (p <= 0) ||
     (eval_pos < 0) || (eval_pos >= n) ||
     !np_glp_qr_drop_workspace_reserve(workspace, n, p))
    return 1;

  memset(workspace->y, 0, (size_t)n*sizeof(double));
  for(j = 0; j < p; j++){
    workspace->pivot[j] = j + 1;
    for(i = 0; i < n; i++){
      const double w = kw[i];
      workspace->xqr[i + j*n] =
        ((w > 0.0) ? sqrt(w) : 0.0) * basis[j][i];
    }
  }

  F77_NAME(dqrdc2)(workspace->xqr, &ldx, &n, &p, (double *)&tol,
                   &rank, workspace->qraux, workspace->pivot,
                   workspace->work);
  if((rank < 0) || (rank > p) || (rank < p))
    return 1;

  for(i = 0; i < rank; i++){
    double s = basis[i][eval_pos];
    for(j = 0; j < i; j++)
      s -= workspace->xqr[j + i*n]*workspace->y[j];
    if(fabs(workspace->xqr[i + i*n]) <= DBL_MIN)
      return 1;
    workspace->y[i] = s/workspace->xqr[i + i*n];
  }

  F77_NAME(dqrqy)(workspace->xqr, &n, &p, workspace->qraux,
                  workspace->y, &ny, workspace->qy);

  for(i = 0; i < n; i++){
    const double w = kw[i];
    row_out[i] = ((w > 0.0) ? sqrt(w) : 0.0) * workspace->qy[i];
  }

  return 0;
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
  free(workspace->rcond_matrix);
  free(workspace->rcond_values);
  free(workspace->rcond_work);
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
  double *rcond_matrix = NULL, *rcond_values = NULL;

  if((workspace == NULL) || (p <= 0) || (nrhs <= 0))
    return 0;
  if((workspace->p_capacity >= p) &&
     (workspace->nrhs_capacity >= nrhs) &&
     (workspace->gram != NULL) &&
     (workspace->rhs != NULL) &&
     (workspace->ipiv != NULL) &&
     (workspace->rcond_matrix != NULL) &&
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
  rcond_matrix = (double *)malloc(gram_bytes);
  rcond_values = (double *)malloc(values_bytes);
  if((gram == NULL) || (rhs == NULL) || (ipiv == NULL) ||
     (rcond_matrix == NULL) || (rcond_values == NULL)){
    free(gram);
    free(rhs);
    free(ipiv);
    free(rcond_matrix);
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
  workspace->rcond_matrix = rcond_matrix;
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
     (workspace->rcond_matrix == NULL) ||
     (workspace->rcond_values == NULL) ||
     (workspace->gram == NULL))
    return 1;

  for(j = 0; j < p; j++)
    for(i = 0; i < p; i++)
      workspace->rcond_matrix[i + j*p] =
        workspace->gram[i + j*p];

  if((workspace->rcond_work == NULL) ||
     (workspace->rcond_lwork_capacity <= 0)){
    int lwork_query = -1;
    double work_query = 0.0;
    int requested_lwork;
    double *rcond_work;

    F77_CALL(dsyev)(&jobz, &uplo, &p,
                    workspace->rcond_matrix, &p,
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
                  workspace->rcond_matrix, &p,
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
