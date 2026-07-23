/* Reusable contiguous LAPACK workspace for canonical local-polynomial solves. */

#include <limits.h>
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
