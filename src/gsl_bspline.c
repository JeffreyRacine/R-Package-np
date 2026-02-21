/* We would like to acknowledge the contributions of the GNU GSL
   authors. In particular, we adapt the GNU GSL B-spline routine
   gsl_bspline.c adding automated support for quantile knots (in
   addition to uniform knots), providing missing functionality for
   derivatives, and for extending the splines beyond their
   endpoints. The source files were downloaded from
   http://www.gnu.org/software/gsl/ version 1.14.*, distributed under
   the terms of the GPL, version 2 or later. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl_bspline.h"
#include <R.h>
#include <Rinternals.h>

/* Code to replicate bs() in splines package. Note that feeding in
   x_min and x_max is necessary if you want to replicate predict.bs()
   - use x_min/x_max for your training data and let x be the
   evaluation data. 

   knots_int is an integer (0=uniform knots, 1=quantile knots) and
   quantile_vector a vector of knots.*/

int gsl_bspline(double *x,
                int *n,
                int *degree,
                int *nbreak,
                double *x_min,
                double *x_max,
                double *quantile_vector,
                int *knots_int,
                double *Bx)
{

  int k = *degree + 1; /* k in gsl */
  int ncoeffs;
  int i, j;

  gsl_bspline_workspace *bw = gsl_bspline_alloc(k, *nbreak);
  ncoeffs = (int)gsl_bspline_ncoeffs(bw);  /* *nbreak+k-2 */
  gsl_vector *B = gsl_vector_alloc(ncoeffs);
  gsl_vector *quantile_vec = gsl_vector_alloc(*nbreak);

  /* 7/12/10 added support for quantile knots */

  if(*knots_int == 0) {
    gsl_bspline_knots_uniform(*x_min, *x_max, bw);
  } else {
    for(i = 0; i < *nbreak; i++) gsl_vector_set(quantile_vec, i, quantile_vector[i]);
    gsl_bspline_knots(quantile_vec, bw);
  }

  for (i = 0; i < *n; ++i)
    {

      /* compute B_j(xi) for all j */
      gsl_bspline_eval(x[i], B, bw);
      
      /* fill in row i of Bx */
      for (j = 0; j < ncoeffs; ++j)
        {
          double Bj = gsl_vector_get(B, j);
          Bx[i*ncoeffs+j] = Bj; /* Bx:*n-by-(*nbreak+*degree-1) */
        }
    }
  
  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(quantile_vec);

  return(0);

} /* main() */

/* Provide missing functionality derivative bs() function in package
   splines */

int gsl_bspline_deriv(double *x,
                      int *n,
                      int *degree,
                      int *nbreak,
                      int *order,
											int *order_max, 
                      double *x_min,
                      double *x_max,
                      double *quantile_vector,
                      int *knots_int,
                      double *Bx)
{

  int k = *degree + 1; /* k in gsl */
  int ncoeffs;
  size_t i, j;

  gsl_bspline_workspace *bw = gsl_bspline_alloc(k, *nbreak);
  ncoeffs = (int)gsl_bspline_ncoeffs(bw);
  gsl_vector *dBorder = gsl_vector_alloc(ncoeffs);
  gsl_bspline_deriv_workspace *derivWS = gsl_bspline_deriv_alloc(k);
	gsl_matrix *dB = gsl_matrix_alloc(ncoeffs, *order_max+1);

	gsl_vector *quantile_vec = gsl_vector_alloc(*nbreak);

	/* 7/12/10 added support for quantile knots */

	if(*knots_int == 0) {
			gsl_bspline_knots_uniform(*x_min, *x_max, bw);
	} else {
			for(i = 0; i < *nbreak; i++) gsl_vector_set(quantile_vec, i, quantile_vector[i]);
			gsl_bspline_knots(quantile_vec, bw);
	}

	for (i = 0; i < *n; ++i)
	{

			/* compute B_j(xi) for all j */
			gsl_bspline_deriv_eval(x[i], order[i], dB, bw, derivWS);

			/* fill in row i of Bx */
			gsl_matrix_get_col(dBorder, dB, order[i]);

			for (j = 0; j < ncoeffs; ++j)
			{
					double Bj = gsl_vector_get(dBorder, j);
					Bx[i*ncoeffs+j] = Bj;
			}
	}

	gsl_bspline_free(bw);
	gsl_vector_free(dBorder);
	gsl_matrix_free(dB);
	/*  gsl_vector_free(quantile_vec);*/
	gsl_bspline_deriv_free(derivWS);

	return(0);

} /* main() */

SEXP C_gsl_bspline(SEXP x,
                   SEXP degree,
                   SEXP nbreak,
                   SEXP x_min,
                   SEXP x_max,
                   SEXP knots,
                   SEXP knots_int)
{
  SEXP x_r = R_NilValue;
  SEXP knots_r = R_NilValue;
  SEXP bx = R_NilValue;
  int n;
  int degree_i;
  int nbreak_i;
  int knots_int_i;
  double x_min_d;
  double x_max_d;
  int ncol;

  PROTECT(x_r = coerceVector(x, REALSXP));
  PROTECT(knots_r = coerceVector(knots, REALSXP));

  n = (int) XLENGTH(x_r);
  degree_i = asInteger(degree);
  nbreak_i = asInteger(nbreak);
  knots_int_i = asInteger(knots_int);
  x_min_d = asReal(x_min);
  x_max_d = asReal(x_max);

  ncol = nbreak_i + degree_i - 1;
  PROTECT(bx = allocMatrix(REALSXP, n, ncol));

  gsl_bspline(REAL(x_r),
              &n,
              &degree_i,
              &nbreak_i,
              &x_min_d,
              &x_max_d,
              REAL(knots_r),
              &knots_int_i,
              REAL(bx));

  UNPROTECT(3);
  return bx;
}

SEXP C_gsl_bspline_deriv(SEXP x,
                         SEXP degree,
                         SEXP nbreak,
                         SEXP deriv,
                         SEXP x_min,
                         SEXP x_max,
                         SEXP knots,
                         SEXP knots_int)
{
  SEXP x_r = R_NilValue;
  SEXP deriv_i = R_NilValue;
  SEXP knots_r = R_NilValue;
  SEXP bx = R_NilValue;
  int n;
  int degree_i;
  int nbreak_i;
  int knots_int_i;
  int order_max;
  double x_min_d;
  double x_max_d;
  int ncol;

  PROTECT(x_r = coerceVector(x, REALSXP));
  PROTECT(deriv_i = coerceVector(deriv, INTSXP));
  PROTECT(knots_r = coerceVector(knots, REALSXP));

  n = (int) XLENGTH(x_r);
  degree_i = asInteger(degree);
  nbreak_i = asInteger(nbreak);
  knots_int_i = asInteger(knots_int);
  x_min_d = asReal(x_min);
  x_max_d = asReal(x_max);

  if ((int)XLENGTH(deriv_i) != n)
    error("C_gsl_bspline_deriv: deriv length must match x length");

  order_max = 0;
  if (n > 0) {
    int i;
    for (i = 0; i < n; i++)
      if (INTEGER(deriv_i)[i] > order_max)
        order_max = INTEGER(deriv_i)[i];
  }

  ncol = nbreak_i + degree_i - 1;
  PROTECT(bx = allocMatrix(REALSXP, n, ncol));

  gsl_bspline_deriv(REAL(x_r),
                    &n,
                    &degree_i,
                    &nbreak_i,
                    INTEGER(deriv_i),
                    &order_max,
                    &x_min_d,
                    &x_max_d,
                    REAL(knots_r),
                    &knots_int_i,
                    REAL(bx));

  UNPROTECT(4);
  return bx;
}
