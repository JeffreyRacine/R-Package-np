#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include "headers.h"

static void check_nonnegative_dim(int value, const char *name, const char *type)
{
  if(value < 0) {
    error("\nFATAL ERROR: Negative allocation dimension %s=%d (type %s). Program terminated.\n",
          name, value, type);
  }
}

static size_t checked_size_mul(size_t a, size_t b, const char *type)
{
  if(a != 0 && b > ((size_t)-1) / a) {
    error("\nFATAL ERROR: Allocation size overflow (type %s). Program terminated.\n",
          type);
  }
  return a * b;
}

/*
 * This function allocates an n by k array of double precision floating  point numbers.
 */

/* allocates a matrix as a contiguous chunk of memory */
double **alloc_tmatd(int nrows, int ncols)
{
  int i;

  double **m, *f;
  size_t ptr_bytes, cell_count, data_bytes;

  check_nonnegative_dim(nrows, "nrows", "DBL_MATRIX");
  check_nonnegative_dim(ncols, "ncols", "DBL_MATRIX");

	/* malloc() on 64 bit systems seems to barf on malloc(0) */

	/*	if(ncols == 0) ncols++;*/

  if(((size_t)ncols * (size_t)nrows) != 0) {
    ptr_bytes = checked_size_mul((size_t)ncols, sizeof(double*), "DBL_MATRIX");
    cell_count = checked_size_mul((size_t)nrows, (size_t)ncols, "DBL_MATRIX");
    data_bytes = checked_size_mul(cell_count, sizeof(double), "DBL_MATRIX");

    if((m=(double**)malloc(ptr_bytes))==NULL){
      error("\nFATAL ERROR: Memory allocation failure (type DBL_MATRIX). Program terminated.\n");
    }

    if ((m[0]=(double*)malloc(data_bytes))==NULL){
      free(m);
      error("\nFATAL ERROR: Memory allocation failure (type DBL_MATRIX). Program terminated.\n");
    }

    f = m[0];
    for(i=1; i<ncols; i++){
      f += nrows;
      m[i] = f;
    }

    return(m);
  
  } else {
    return(NULL);
  }

}


double **alloc_matd(int nrows, int ncols)
{
  int i;

  double **m;
  size_t ptr_bytes, row_bytes;

  check_nonnegative_dim(nrows, "nrows", "DBL_MATRIX");
  check_nonnegative_dim(ncols, "ncols", "DBL_MATRIX");

  /* malloc() on 64 bit systems seems to barf on malloc(0) */
  
  /*	if(ncols == 0) ncols++;*/
  
  if(ncols != 0) {

    ptr_bytes = checked_size_mul((size_t)ncols, sizeof(double*), "DBL_MATRIX");
    row_bytes = checked_size_mul((size_t)nrows, sizeof(double), "DBL_MATRIX");

    if((m=(double**)malloc(ptr_bytes))==NULL) {
      error("\nFATAL ERROR: Memory allocation failure (type DBL_MATRIX). Program terminated.\n");
    }

    for(i=0;i<ncols;i++) {
      if((m[i]=(double*)malloc(row_bytes))==NULL) {
        free_mat(m, i);
        error("\nFATAL ERROR: Memory allocation failure (type DBL_MATRIX). Program terminated.\n");
      }
    }

    return(m);

  } else {
    return(NULL);
  }

}


/*
 * This function frees an n by k matrix which has been dynamically allocated.
 * It assumes that x has been declared as NULL.
 * Therefore it only frees an allocated array.
 */

void free_tmat(double **x){
  if(x != NULL){
    free(x[0]);
    free(x);
  }
  return;
}


void free_mat(double **x, int n)
{
  int i;

  if(x != NULL)
  {
    for(i=0; i<n ; i++)
    {
      free(x[i]);
    }

    free(x);
  }

  return;
}



#include <stdlib.h>
#include <stdio.h>

/* This function allocates a double precision floating  point vector of observations.
   It returns pointer to memory on success and NULL on failure.    */

/* This function allocates a double precision double precision floating  point vector of observations.
   It returns pointer to memory on success and NULL on failure.    */

double *alloc_vecd(int nobs)
{
  double *a;
  size_t bytes;

  check_nonnegative_dim(nobs, "nobs", "DBL_VECTOR");

	/* malloc() on 64 bit systems seems to barf on malloc(0) */

	/*	if(nobs == 0) nobs++;*/

	if(nobs != 0) {

  bytes = checked_size_mul((size_t)nobs, sizeof(double), "DBL_VECTOR");

  if ((a=(double *)malloc(bytes))==NULL)
  {
    error("\nFATAL ERROR: Memory allocation failure (type DBL_VECTOR). Program terminated.\n");
  }

  return(a);
	} else {
		return(NULL);
	}

}


int *alloc_vecu(int nobs)
{
  int *a;
  size_t bytes;

  check_nonnegative_dim(nobs, "nobs", "INT_VECTOR");

	/* malloc() on 64 bit systems seems to barf on malloc(0) */

	/*	if(nobs == 0) nobs++;*/

	if(nobs != 0) {

  bytes = checked_size_mul((size_t)nobs, sizeof(int), "INT_VECTOR");

  if ((a=(int *)malloc(bytes))==NULL)
  {
    error("\nFATAL ERROR: Memory allocation failure (type INT_VECTOR). Program terminated.\n");
  }

  return(a);
	} else {
		return(NULL);
	}

}
