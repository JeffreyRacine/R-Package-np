/* Note:
 *  This file includs the original gsl_bspline.h of gsl/bspline. We
 *  merge all other necessary structures and definitions of gsl routines 
 *  to this file,  so now we can compile and link it without linking to the gsl library.
 *  
 *  The source files are downloaded from http://www.gnu.org/software/gsl/,  
 *  and the current version is 1.14.*
 *
 *
 *
 * */



/* bspline/gsl_bspline.h
 *
 * Copyright (C) 2006 Patrick Alken
 * Copyright (C) 2008 Rhys Ulerich
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_BSPLINE_H__
#define __GSL_BSPLINE_H__

#ifdef WINDOWS
#include <windows.h>   /*   For easy self window porting,   without neading everything R.h needs */
#include <stdarg.h>
#else
#include <R.h>
#endif
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "RuniqueCombs.h"


/* *********************************************************** */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define NULL_VECTOR_VIEW {{0,  0,  0,  0,  0}}
#define NULL_VECTOR {0,  0,  0,  0,  0}

#define GSL_ERROR(reason, gsl_errno) ErrorMessage(reason, gsl_errno)
#define GSL_RANGE_COND(x) (x)
#define RETURN_IF_NULL(x) if (!x) { return ; }
#define GSL_MAX_INT(a, b) ((a) > (b) ? (a) : (b))
#define GSL_MIN_INT(a, b) ((a) < (b) ? (a) : (b))



/* GSL_ERROR suitable for out-of-memory conditions */

enum { 
  GSL_SUCCESS  = 0, 
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;


/*  void ErrorMessage(char *msg,int fatal);
*/
/* ************vector********* */



struct gsl_block_struct
{
  size_t size;
  double *data;
};

typedef struct gsl_block_struct gsl_block;

typedef struct 
{
  size_t size;
  size_t stride;
  double *data;
  gsl_block *block;
  int owner;
} 
gsl_vector;
typedef struct
{
  gsl_vector vector;
} _gsl_vector_view;

typedef _gsl_vector_view gsl_vector_view;

typedef struct
{
  gsl_vector vector;
} _gsl_vector_const_view;

typedef const _gsl_vector_const_view gsl_vector_const_view;



gsl_block *gsl_block_alloc (const size_t n);

gsl_vector *gsl_vector_alloc (const size_t n);
double gsl_vector_get (const gsl_vector * v, const size_t i);
void gsl_vector_set (gsl_vector * v, const size_t i, double x);
void gsl_vector_free (gsl_vector * v);
void gsl_block_free(gsl_block * b);



/* ***********************matrix************* */
typedef struct 
{
  size_t size1;
  size_t size2;
  size_t tda;
  double * data;
  gsl_block * block;
  int owner;
} gsl_matrix;


/* Allocation */

gsl_matrix *gsl_matrix_alloc (const size_t n1, const size_t n2);
_gsl_vector_view gsl_matrix_column (gsl_matrix * m, const size_t j);
int gsl_matrix_get_col(gsl_vector * v, const gsl_matrix * m, const size_t j);
void gsl_matrix_free (gsl_matrix * m);
double   gsl_matrix_get(const gsl_matrix * m, const size_t i, const size_t j);
void    gsl_matrix_set(gsl_matrix * m, const size_t i, const size_t j, const double x);


/* *******bspline******** */
typedef struct
{
    size_t k; /* spline order */
    size_t km1; /* k - 1 (polynomial order) */
    size_t l; /* number of polynomial pieces on interval */
    size_t nbreak; /* number of breakpoints (l + 1) */
    size_t n; /* number of bspline basis functions (l + k - 1) */

    gsl_vector *knots; /* knots vector */
    gsl_vector *deltal; /* left delta */
    gsl_vector *deltar; /* right delta */
    gsl_vector *B; /* temporary spline results */
} gsl_bspline_workspace;

typedef struct
{
    size_t k; /* spline order */
    gsl_matrix *A; /* work matrix */
    gsl_matrix *dB; /* temporary derivative results */
} gsl_bspline_deriv_workspace;


/* *********************************************************** */


gsl_bspline_workspace *gsl_bspline_alloc(const size_t k, const size_t nbreak);

void gsl_bspline_free(gsl_bspline_workspace *w);

size_t gsl_bspline_ncoeffs(gsl_bspline_workspace * w);
size_t gsl_bspline_order(gsl_bspline_workspace * w);
size_t gsl_bspline_nbreak(gsl_bspline_workspace * w);
double gsl_bspline_breakpoint(size_t i, gsl_bspline_workspace * w);
double gsl_bspline_greville_abscissa(size_t i, gsl_bspline_workspace *w);

int gsl_bspline_knots(const gsl_vector *breakpts, gsl_bspline_workspace *w);

int gsl_bspline_knots_uniform(const double a, const double b,
                              gsl_bspline_workspace *w);

int
gsl_bspline_eval(const double x, gsl_vector *B, 
                 gsl_bspline_workspace *w);

int
gsl_bspline_eval_nonzero(const double x,
                         gsl_vector *Bk,
                         size_t *istart,
                         size_t *iend,
                         gsl_bspline_workspace *w);

gsl_bspline_deriv_workspace *
gsl_bspline_deriv_alloc(const size_t k);

void
gsl_bspline_deriv_free(gsl_bspline_deriv_workspace *w);

int
gsl_bspline_deriv_eval(const double x,
                       const size_t nderiv,
                       gsl_matrix *dB,
                       gsl_bspline_workspace *w,
                       gsl_bspline_deriv_workspace *dw);

int
gsl_bspline_deriv_eval_nonzero(const double x,
                               const size_t nderiv,
                               gsl_matrix *dB,
                               size_t *istart,
                               size_t *iend,
                               gsl_bspline_workspace *w,
                               gsl_bspline_deriv_workspace *dw);


#endif /* __GSL_BSPLINE_H__ */
