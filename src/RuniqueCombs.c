/*
 * =====================================================================================
 *
 *       Filename:  RuniqueCombs.c
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
#ifdef WINDOWS
#include <windows.h>   /*  For easy self window porting,  without neading everything R.h needs */
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

#define RANGECHECK
#define PAD 1L

#define ROUND(a) ((a)-(int)floor(a)>0.5) ? ((int)floor(a)+1):((int)floor(a))


matrix null_mat; /* matrix for passing when you don't actually need to */
#define PADCON (-1.234565433647588392902028934e270)

/* counter for memory used */


long memused=0L,matrallocd=0L;

/* the routines */
  
struct mrec
{ matrix mat;
  struct mrec *fp,*bp;
};
typedef struct mrec MREC;


void ErrorMessage(char *msg,int fatal)
{ 
#ifdef WINDOWS
  MessageBox(HWND_DESKTOP,msg,"Info!",MB_ICONEXCLAMATION|MB_OK); 
#else  
if (fatal) error("%s",msg);else warning("%s",msg);
#endif
}

matrix null_mat;
MREC *top,*bottom;


matrix Rmatrix(double *A,long r,long c)

/* produces a matrix from the array containing a (default) R matrix stored:
   A[0,0], A[1,0], A[2,0] .... etc */

{ int i,j;
  matrix M;
  M=initmat(r,c);
  for (i=0;i<r;i++) for (j=0;j<c;j++) M.M[i][j]=A[i+j*r];
  return(M);
}

matrix initmat(rows,cols) long rows,cols;

/* Don't alter this without altering freemat() as well !! */

{ matrix A;long i,j,pad;
#ifdef RANGECHECK
  pad=PAD;
#else
  pad=0L;
#endif
  A.vec=0;
  A.M=(double **)calloc((size_t)(rows+2*pad),sizeof(double *));
  if ((cols==1L)||(rows==1L))
  { if (A.M)
    A.M[0]=(double *)calloc((size_t)(cols*rows+2*pad),sizeof(double));
    for (i=1L;i<rows+2*pad;i++)
    A.M[i]=A.M[0]+i*cols;A.vec=1;
  } else
  { if (A.M)
    for (i=0L;i<rows+2*pad;i++)
    A.M[i]=(double *)calloc((size_t)(cols+2*pad),sizeof(double));
  }
  A.mem=rows*cols*sizeof(double);
  memused+=A.mem;matrallocd++;
  A.original_r=A.r=rows;A.original_c=A.c=cols;
  if (((!A.M)||(!A.M[rows-1+2*pad]))&&(rows*cols>0L))
  { ErrorMessage("Failed to initialize memory for matrix.",1);}
  if (pad)  /* This lot is debugging code that checks out matrix errors
		       on allocation and release */
  { if (A.vec)
    { A.V=A.M[0];for (i=0;i<pad;i++) { A.V[i]=PADCON;A.V[i+pad+A.r*A.c]=PADCON;}
    } else
    { for (i=0;i<A.r+2*pad;i++)
      { for (j=0;j<pad;j++) A.M[i][j]=PADCON;
	     for (j=A.c+pad;j<A.c+2*pad;j++) A.M[i][j]=PADCON;
      }
      for (i=0;i<A.c+2*pad;i++)
      { for (j=0;j<pad;j++) A.M[j][i]=PADCON;
	     for (j=A.r+pad;j<A.r+2*pad;j++) A.M[j][i]=PADCON;
      }
    }
    for (i=0;i<A.r+2*pad;i++)
    for (j=0;j<pad;j++) A.M[i]++;  /* shifting pointers forward past padding */
    if (!A.vec) for (j=0;j<pad;j++) A.M++;
    A.V=A.M[0];
  /* putting a record of the matrix on the linked list of all extant matrices */
    if (matrallocd==1) /*new list*/
    { top=bottom=(MREC *)calloc(1,sizeof(MREC));
      bottom->mat=top->mat=A;top->bp=bottom;bottom->fp=top;
    } else  /* expanding the linked list by one */
    { top->fp=(MREC *)calloc(1,sizeof(MREC));
      top->fp->mat=A;top->fp->bp=top;top=top->fp; /* crystal clear, no? */
    }
  }
  A.V=A.M[0];/* This allows vectors to be accessed using A.V[i] */
  return(A);
}


void mcopy(matrix *A,matrix *B)

/* copies A into B */

{ long Ac;
  double *pA,*pB,**AM,**BM;
  if (A->r>B->r||A->c>B->c) ErrorMessage("Target matrix too small in mcopy",1);
  BM=B->M;Ac=A->c;
  for (AM=A->M;AM<A->M+A->r;AM++)
  { pB= *BM;
    for (pA= *AM;pA< *AM+Ac; pA++) *(pB++) = *pA;
    BM++;
  }
}

void freemat(A) matrix A;

{ long i,j,pad;int ok=1;
  MREC *delet;
#ifdef RANGECHECK
  pad=PAD;
#else
  pad=0L;
#endif
/*  if (A.original_r*A.original_c!=0L) */
  { if (pad)
    { if (A.vec)
      { for (i=-pad;i<0;i++)
	     if ((A.V[i]!=PADCON)||(A.V[i+A.original_r*A.original_c+pad]!=PADCON))
	     ok=0;
      } else
      { for (i=-pad;i<A.original_r+pad;i++)
	     { for (j=A.original_c;j<A.original_c+pad;j++) if (A.M[i][j]!=PADCON) ok=0;
	       for (j=-pad;j<0;j++) if (A.M[i][j]!=PADCON) ok=0;
	     }
	     for (i=-pad;i<A.original_c+pad;i++)
	     { for (j=A.original_r;j<A.original_r+pad;j++) if (A.M[j][i]!=PADCON) ok=0;
	       for (j=-pad;j<0;j++) if (A.M[j][i]!=PADCON) ok=0;
	     }
      }
      if (!ok)
      { ErrorMessage("An out of bound write to matrix has occurred!",1);
      }
      /* find the matrix being deleted in the linked list of extant matrices */
      i=0L;delet=bottom;
      while ((i<matrallocd)&&(delet->mat.M!=A.M)) { i++;delet=delet->fp;}
      if (i==matrallocd)
      { ErrorMessage("INTEGRITY PROBLEM in the extant matrix list.",1);
      } else
      { if (i)
	     delet->bp->fp=delet->fp;
	     else bottom=delet->fp;
	     if (i!=matrallocd-1)
	     delet->fp->bp=delet->bp;
	     else top=delet->bp;
	     free(delet);
      }
      /* repositioning pointers so that what was allocated gets freed */
      if (!A.vec) for (i=0;i<pad;i++) A.M--;
      for (i=0;i<A.original_r+2*pad;i++)
      for (j=0;j<pad;j++) A.M[i]--;
    }
    if (A.vec) free(A.M[0]); else
    for (i=0;i<A.original_r+2*pad;i++) if (A.M[i]) free(A.M[i]);
    if (A.M) free(A.M);
    memused -= A.mem;matrallocd--;
  }
}


void RArrayFromMatrix(double *a,long r,matrix *M)

/* copies matrix *M into R array a where r is the number of rows of A treated as
  a matrix by R */

{ int i,j;
  for (i=0;i<M->r;i++) for (j=0;j<M->c;j++) a[i+r*j]=M->M[i][j];
}


int *Xd_strip(matrix *Xd)

/* The rows of Xd (excluding last col) contain the covariate values for
   a set of data to which a thin plate spline is to be fitted. The purpose
   of this routine is to locate co-incident points, and strip out redundant
   copies of these points. At the same time a record is kept of what has 
   been done, so that the function returns an array yxindex, such that 
   yxindex[i] contains the row of the stripped down Xd that corresponds to 
   the ith response datum. Note that the identification of ties involves 
   sorting Xd - even if there are no ties.
    
   Note that this routine assumes that the final column of Xd consists of the 
   integers 0 to Xd->r-1. These are vital for constructing the index.

   On exit Xd->r will contain the number of unique covariate points.
*/

{ int *yxindex,start,stop,ok,i;
  /*  long Xdor;*/
  double xi,**dum;
  yxindex = (int *)calloc((size_t)Xd->r,sizeof(int));
  dum = (double **)calloc((size_t)Xd->r,sizeof(double *));
  msort(*Xd);
  /*  Xdor=Xd->r;  keep record of original length of Xd */
  start=stop=0;ok=1;
  while(ok)
  { /* look for start of run of equal rows ..... */
    while(start<Xd->r-1&&!Xd_row_comp(Xd->M[start],Xd->M[start+1],(int)Xd->c-1)) 
    { /* Xd->M[start] not tied with anything, nothing to erase.... */
      xi=Xd->M[start][Xd->c-1];
      yxindex[ROUND(xi)]=start;
      start++;
    }
    if (start==Xd->r-1) 
    { ok=0; /* reached end with no more ties */
      xi=Xd->M[start][Xd->c-1];
      yxindex[ROUND(xi)]=start; /* final index entry needed */
    }
    if (ok) /* search for end of run */
    { stop=start+1;
      while(stop<Xd->r-1&&Xd_row_comp(Xd->M[stop],Xd->M[stop+1],(int)Xd->c-1)) stop++;
      for (i=start;i<=stop;i++) /* fill out the index array */
      { xi=Xd->M[i][Xd->c-1];
        yxindex[ROUND(xi)]=start;
        dum[i-start]=Xd->M[i]; /* Rows stored to copy back onto end, so matrix can be freed properly */
      }
      for (i=stop+1;i<Xd->r;i++)
      { Xd->M[i-stop+start]=Xd->M[i];}
      Xd->r -= stop-start;
      for (i=1;i<=stop-start;i++)
      { Xd->M[Xd->r-1+i]=dum[i];}
    } 
  }
  free(dum); 
  return(yxindex);
}

int Xd_row_comp(double *a,double *b,int k)

/* service routine for Xd_strip(), compares k elements of two rows for equality */

{ int i;
  for (i=0;i<k;i++) if (a[i]!=b[i]) return(0);
  return(1);
}

int real_elemcmp(const void *a,const void *b,int el)

{ static int k=0;
  int i;
  double *na,*nb;
  if (el>=0) { k=el;return(0);}
  na=(*(double **)a);nb=(*(double **)b);
  for (i=0;i<k;i++) 
  { if (na[i]<nb[i]) return(-1);
    if (na[i]>nb[i]) return(1);
  }
  return(0);
}

int melemcmp(const void *a,const void *b)

{ return(real_elemcmp(a,b,-1));
}


void msort(matrix a)

/* sorts a matrix, in situ, using standard routine qsort so 
   that its first col is in ascending order, its second col
   is in ascending order for any ties in the first col, and 
   so on.....
*/

{ double z=0.0;
  real_elemcmp(&z,&z,(int)a.c); 
  qsort(a.M,(size_t)a.r,sizeof(a.M[0]),melemcmp);
}

void RuniqueCombs(double *X,int *ind,int *r, int *c)

/* X is a matrix. This routine finds its unique rows and strips out the 
   duplicates. This is useful for finding out the number of unique covariate
   combinations present in a set of data. */

{ matrix B,Xd;
  int i,*ind1;
  B=Rmatrix(X,(long)(*r),(long)(*c));
  Xd=initmat(B.r,B.c+1);
  Xd.c--;mcopy(&B,&Xd);freemat(B);Xd.c++;
  for (i=0;i<Xd.r;i++) Xd.M[i][Xd.c-1]=(double)i;
  ind1=Xd_strip(&Xd);
  for (i=0;i<*r;i++) ind[i] = ind1[i]; /* copy index for return */
  Xd.c--; /* hide index array  */
  RArrayFromMatrix(X,Xd.r,&Xd);  /* NOTE: not sure about rows here!!!! */
  *r = (int)Xd.r; 
  freemat(Xd);free(ind1);
#ifdef MEM_CHECK
  dmalloc_log_unfreed();  dmalloc_verify(NULL);
#endif 
}


