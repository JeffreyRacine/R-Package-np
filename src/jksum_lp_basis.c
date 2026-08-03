#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "jksum_lp_basis.h"

static int np_lp_lwork_from_query(const double query, int *lwork)
{
  double rounded;

  if((lwork == NULL) || !isfinite(query) || !(query >= 1.0) ||
     (query > (double)INT_MAX))
    return 0;
  rounded = ceil(query);
  if(rounded > (double)INT_MAX)
    return 0;
  *lwork = (int)rounded;
  return 1;
}

NPLPBasisStatus np_lp_basis_requires_conditioning(
  double * const *source,
  int n,
  int p,
  int leading_dimension,
  double min_rcond,
  int *required)
{
  const char jobz = 'N';
  const char trans = 'T';
  const char uplo = 'U';
  const double alpha = 1.0;
  const double beta = 0.0;
  double *gram = NULL;
  double *eigenvalues = NULL;
  double *work = NULL;
  double work_query = 0.0;
  double maximum_eigenvalue = 0.0;
  double minimum_eigenvalue = DBL_MAX;
  size_t gram_elements;
  int lwork = -1;
  int info = 0;
  int i;
  int j;
  int row;
  NPLPBasisStatus status = NP_LP_BASIS_INVALID;

  if((source == NULL) || (required == NULL) || (n <= 0) ||
     (p <= 0) || (p > n) || (leading_dimension < n) ||
     !isfinite(min_rcond) ||
     (min_rcond < 0.0))
    return status;
  *required = 0;
  if(((size_t)p > SIZE_MAX/(size_t)p) ||
     ((gram_elements = (size_t)p*(size_t)p) > SIZE_MAX/sizeof(double)) ||
     ((size_t)p > SIZE_MAX/sizeof(double)))
    return NP_LP_BASIS_MEMORY;

  gram = (double *)calloc(gram_elements, sizeof(double));
  eigenvalues = (double *)malloc((size_t)p*sizeof(double));
  if((gram == NULL) || (eigenvalues == NULL)){
    status = NP_LP_BASIS_MEMORY;
    goto cleanup;
  }

  for(j = 0; j < p; j++){
    if(source[j] == NULL)
      goto cleanup;
    if(source[j] != source[0] + (size_t)j*(size_t)leading_dimension)
      goto cleanup;
    for(row = 0; row < n; row++)
      if(!isfinite(source[j][row])){
        *required = 1;
        status = NP_LP_BASIS_OK;
        goto cleanup;
      }
  }
  F77_CALL(dsyrk)(&uplo, &trans, &p, &n, &alpha,
                  source[0], &leading_dimension, &beta, gram, &p
                  FCONE FCONE);
  for(j = 0; j < p; j++)
    for(i = 0; i <= j; i++)
      if(!isfinite(gram[i + (size_t)j*(size_t)p])){
        *required = 1;
        status = NP_LP_BASIS_OK;
        goto cleanup;
      }

  F77_CALL(dsyev)(&jobz, &uplo, &p, gram, &p, eigenvalues,
                  &work_query, &lwork, &info FCONE FCONE);
  if((info != 0) || !np_lp_lwork_from_query(work_query, &lwork)){
    *required = 1;
    status = NP_LP_BASIS_OK;
    goto cleanup;
  }
  if((size_t)lwork > SIZE_MAX/sizeof(double)){
    status = NP_LP_BASIS_MEMORY;
    goto cleanup;
  }
  work = (double *)malloc((size_t)lwork*sizeof(double));
  if(work == NULL){
    status = NP_LP_BASIS_MEMORY;
    goto cleanup;
  }
  F77_CALL(dsyev)(&jobz, &uplo, &p, gram, &p, eigenvalues,
                  work, &lwork, &info FCONE FCONE);
  if(info != 0){
    *required = 1;
    status = NP_LP_BASIS_OK;
    goto cleanup;
  }

  for(i = 0; i < p; i++){
    const double eigenvalue = eigenvalues[i];

    if(!isfinite(eigenvalue) || !(eigenvalue > 0.0)){
      *required = 1;
      status = NP_LP_BASIS_OK;
      goto cleanup;
    }
    if(eigenvalue > maximum_eigenvalue)
      maximum_eigenvalue = eigenvalue;
    if(eigenvalue < minimum_eigenvalue)
      minimum_eigenvalue = eigenvalue;
  }
  *required = !(maximum_eigenvalue > 0.0) ||
    ((minimum_eigenvalue/maximum_eigenvalue) < min_rcond);
  status = NP_LP_BASIS_OK;

cleanup:
  free(work);
  free(eigenvalues);
  free(gram);
  return status;
}

NPLPBasisStatus np_lp_conditioned_basis_prepare(
  double * const *source,
  int n,
  int p,
  int leading_dimension,
  double *destination)
{
  int *pivot = NULL;
  double *tau = NULL;
  double *work = NULL;
  double work_query = 0.0;
  double maximum_diagonal = 0.0;
  double minimum_diagonal = DBL_MAX;
  double rank_floor;
  size_t p_size;
  int lwork = -1;
  int info = 0;
  int i;
  int j;
  NPLPBasisStatus status = NP_LP_BASIS_INVALID;

  if((source == NULL) || (destination == NULL) || (n <= 0) ||
     (p <= 0) || (p > n) || (leading_dimension < n))
    return status;

  p_size = (size_t)p;
  if((p_size > SIZE_MAX/sizeof(int)) ||
     (p_size > SIZE_MAX/sizeof(double)))
    return NP_LP_BASIS_MEMORY;

  pivot = (int *)calloc(p_size, sizeof(int));
  tau = (double *)malloc(p_size*sizeof(double));
  if((pivot == NULL) || (tau == NULL)){
    status = NP_LP_BASIS_MEMORY;
    goto cleanup;
  }

  for(j = 0; j < p; j++){
    double column_scale = 0.0;

    if(source[j] == NULL)
      goto cleanup;
    for(i = 0; i < n; i++){
      const double value = source[j][i];
      const double magnitude = fabs(value);

      if(!isfinite(value)){
        status = NP_LP_BASIS_NONFINITE;
        goto cleanup;
      }
      if(magnitude > column_scale)
        column_scale = magnitude;
    }
    if(!(column_scale > 0.0)){
      status = NP_LP_BASIS_RANK_DEFICIENT;
      goto cleanup;
    }
    for(i = 0; i < n; i++)
      destination[i + (size_t)j*(size_t)leading_dimension] =
        source[j][i]/column_scale;
  }

  F77_CALL(dgeqp3)(&n, &p, destination, &leading_dimension,
                   pivot, tau, &work_query, &lwork, &info);
  if((info != 0) || !np_lp_lwork_from_query(work_query, &lwork)){
    status = NP_LP_BASIS_LAPACK;
    goto cleanup;
  }
  if((size_t)lwork > SIZE_MAX/sizeof(double)){
    status = NP_LP_BASIS_MEMORY;
    goto cleanup;
  }
  work = (double *)malloc((size_t)lwork*sizeof(double));
  if(work == NULL){
    status = NP_LP_BASIS_MEMORY;
    goto cleanup;
  }
  F77_CALL(dgeqp3)(&n, &p, destination, &leading_dimension,
                   pivot, tau, work, &lwork, &info);
  if(info != 0){
    status = NP_LP_BASIS_LAPACK;
    goto cleanup;
  }

  for(j = 0; j < p; j++){
    const double diagonal =
      fabs(destination[j + (size_t)j*(size_t)leading_dimension]);

    if(!isfinite(diagonal)){
      status = NP_LP_BASIS_NONFINITE;
      goto cleanup;
    }
    if(diagonal > maximum_diagonal)
      maximum_diagonal = diagonal;
    if(diagonal < minimum_diagonal)
      minimum_diagonal = diagonal;
  }
  if(!(maximum_diagonal > 0.0)){
    status = NP_LP_BASIS_RANK_DEFICIENT;
    goto cleanup;
  }
  rank_floor = DBL_EPSILON*(double)((n > p) ? n : p)*maximum_diagonal;
  if(!(minimum_diagonal > rank_floor)){
    status = NP_LP_BASIS_RANK_DEFICIENT;
    goto cleanup;
  }

  lwork = -1;
  work_query = 0.0;
  F77_CALL(dorgqr)(&n, &p, &p, destination, &leading_dimension,
                   tau, &work_query, &lwork, &info);
  if((info != 0) || !np_lp_lwork_from_query(work_query, &lwork)){
    status = NP_LP_BASIS_LAPACK;
    goto cleanup;
  }
  if((size_t)lwork > SIZE_MAX/sizeof(double)){
    status = NP_LP_BASIS_MEMORY;
    goto cleanup;
  }
  {
    double *larger_work =
      (double *)realloc(work, (size_t)lwork*sizeof(double));
    if(larger_work == NULL){
      status = NP_LP_BASIS_MEMORY;
      goto cleanup;
    }
    work = larger_work;
  }
  F77_CALL(dorgqr)(&n, &p, &p, destination, &leading_dimension,
                   tau, work, &lwork, &info);
  if(info != 0){
    status = NP_LP_BASIS_LAPACK;
    goto cleanup;
  }
  for(j = 0; j < p; j++)
    for(i = 0; i < n; i++)
      if(!isfinite(destination[i + (size_t)j*(size_t)leading_dimension])){
        status = NP_LP_BASIS_NONFINITE;
        goto cleanup;
      }

  status = NP_LP_BASIS_OK;

cleanup:
  free(work);
  free(tau);
  free(pivot);
  return status;
}
