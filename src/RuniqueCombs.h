/*
 * =====================================================================================
 *
 *       Filename:  RuniqueCombs.h
 *
 *    Description:  The routines are copied from the R package mgcv 
 *    which is downloaded from
 *    http://cran.r-project.org/web/packages/mgcv/
 *
 *    Copyright (C) 2000-2005 Simon N. Wood  simon.wood@r-project.org
 *
 *        Created:  30/04/2011 06:22:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * (www.gnu.org/copyleft/gpl.html)

 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
 * USA. 
 *              
 *
 *
 *        
 * =====================================================================================
 */

#ifndef MATRIX_HEADER_IN
#define MATRIX_HEADER_IN
/* The basic matrix structure */
#define TOL 1e-10
#define VEC M[0]

typedef struct
{ int vec;long r,c,mem,original_r,original_c;double **M,*V;} matrix;

extern matrix null_mat;

extern long matrallocd;


void ErrorMessage(char *msg,int fatal);
void RuniqueCombs(double *X,int *ind,int *r, int *c);
matrix Rmatrix(double *A,long r,long c);
matrix initmat(long rows,long cols);
void mcopy(matrix *A,matrix *B);
void freemat(matrix A);
void RArrayFromMatrix(double *a,long r,matrix *M);
int *Xd_strip(matrix *Xd);
int Xd_row_comp(double *a,double *b,int k);
int real_elemcmp(const void *a,const void *b,int el);
void msort(matrix a);
int melemcmp(const void *a,const void *b);


#endif
