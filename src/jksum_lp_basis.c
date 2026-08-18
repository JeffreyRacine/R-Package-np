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
  int contiguous = 1;
  int i;
  int j;
  int row;
  uintptr_t base_address = 0;
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
    size_t offset_elements;
    uintptr_t offset_bytes;

    if(source[j] == NULL)
      goto cleanup;
    if(j == 0){
      base_address = (uintptr_t)(const void *)source[j];
    } else if(((size_t)j > SIZE_MAX/(size_t)leading_dimension) ||
              ((offset_elements = (size_t)j*(size_t)leading_dimension) >
               UINTPTR_MAX/sizeof(double)) ||
              ((offset_bytes = (uintptr_t)offset_elements*sizeof(double)) >
               UINTPTR_MAX - base_address) ||
              ((uintptr_t)(const void *)source[j] !=
               base_address + offset_bytes)){
      contiguous = 0;
    }
    for(row = 0; row < n; row++)
      if(!isfinite(source[j][row])){
        *required = 1;
        status = NP_LP_BASIS_OK;
        goto cleanup;
      }
  }
  if(contiguous){
    F77_CALL(dsyrk)(&uplo, &trans, &p, &n, &alpha,
                    source[0], &leading_dimension, &beta, gram, &p
                    FCONE FCONE);
  } else {
    for(j = 0; j < p; j++)
      for(i = 0; i <= j; i++){
        double sum = 0.0;

        for(row = 0; row < n; row++)
          sum += source[i][row]*source[j][row];
        gram[i + (size_t)j*(size_t)p] = sum;
      }
  }
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
