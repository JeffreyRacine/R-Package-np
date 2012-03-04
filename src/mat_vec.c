#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include "headers.h"

#ifdef RCSID
static char rcsid[] = "$Id: mat_vec.c,v 1.5 2006/11/02 16:56:49 tristen Exp $";
#endif

/*
 * This function allocates an n by k array of double precision floating  point numbers.
 */

/* allocates a matrix as a contiguous chunk of memory */
double **alloc_tmatd(int nrows, int ncols)
{
  int i;

  double **m, *f;

	/* malloc() on 64 bit systems seems to barf on malloc(0) */

	/*	if(ncols == 0) ncols++;*/

  if(ncols*nrows != 0) {
    if((m=(double**)malloc(sizeof(double*)*ncols))==NULL){
      error("\nFATAL ERROR: Memory allocation failure (type DBL_MATRIX). Program terminated.\n");
    }

    if ((m[0]=(double*)malloc(sizeof(double)*nrows*ncols))==NULL){
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

  /* malloc() on 64 bit systems seems to barf on malloc(0) */
  
  /*	if(ncols == 0) ncols++;*/
  
  if(ncols != 0) {

    if((m=(double**)malloc(sizeof(double*)*ncols))==NULL) {
      error("\nFATAL ERROR: Memory allocation failure (type DBL_MATRIX). Program terminated.\n");
    }

    for(i=0;i<ncols;i++) {
      if((m[i]=(double*)malloc(sizeof(double)*nrows))==NULL) {
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

	/* malloc() on 64 bit systems seems to barf on malloc(0) */

	/*	if(nobs == 0) nobs++;*/

	if(nobs != 0) {

  if ((a=(double *)malloc(sizeof(double)*nobs))==NULL)
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

	/* malloc() on 64 bit systems seems to barf on malloc(0) */

	/*	if(nobs == 0) nobs++;*/

	if(nobs != 0) {

  if ((a=(int *)malloc(sizeof(unsigned)*nobs))==NULL)
  {
    error("\nFATAL ERROR: Memory allocation failure (type INT_VECTOR). Program terminated.\n");
  }

  return(a);
	} else {
		return(NULL);
	}

}

/* For pgplot */

