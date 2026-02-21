#include <stdlib.h>
#include <math.h>
#include "matrix.h"

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#ifdef RCSID
static char rcsid[] = "$Id: matrix.c,v 1.4 2006/11/02 19:50:13 tristen Exp $";
#endif



/*
 *-----------------------------------------------------------------------------
 * file: matadd.c
 * desc: matrix addition
 * by: ko shu pui, Patrick
 * date: 24 nov 91 v0.1
 * revi: Nov 25 1994 by J. Racine - added error checking, passing matrix
 *       rather than allocating
 * ref:
 *       [1] Mary L.Boas, "Mathematical Methods in the Physical Science,"
 * John Wiley & Sons, 2nd Ed., 1983. Chap 3.
 *
 *-----------------------------------------------------------------------------
 */


/* Major bug in mat_free - memory not freed!!! (missing parentheses)
   Found 7/6/96, jracine */



/*
 *-----------------------------------------------------------------------------
 *       file:   matcreat.c
 *       desc:   matrix mathematics - object creation
 *       by:     ko shu pui, Patrick
 *       date:   24 nov 91 v0.1
 *       revi:   14 may 92 v0.2
 *               21 may 92 v0.3
 *       ref:
 *       [1] Mary L.Boas, "Mathematical Methods in the Physical Science,"
 *       John Wiley & Sons, 2nd Ed., 1983. Chap 3.
 *
 *       [2] Kendall E.Atkinson, "An Introduction to Numerical Analysis,"
 *       John Wiley & Sons, 1978.
 *
 *-----------------------------------------------------------------------------
 */

MATRIX  _mat_creat( int row, int col )
{
	MATBODY *mat;
	int i;

	if ((mat = (MATBODY *)malloc( sizeof(MATHEAD) + sizeof(double *) * row)) == NULL)
	{
          error("mat: malloc error\n" );
	}

	for (i=0; i<row; i++)
	{
		if ((*((double **)(&mat->matrix) + i) = (double *)malloc(sizeof(double) * col)) == NULL)
		{
			error("mat: malloc error\n" );
		}
	}

	mat->head.row = row;
	mat->head.col = col;

	return (&(mat->matrix));
}


/*
 *-----------------------------------------------------------------------------
 *       funct:  mat_creat
 *       desct:  create a matrix
 *       given:  row, col = dimension, type = which kind of matrix
 *       retrn:  allocated matrix (use mat_free() to free memory)
 *-----------------------------------------------------------------------------
 */

MATRIX  mat_creat( int row, int col, int type )
{
	MATRIX  A;

	if ((A =_mat_creat( row, col )) != NULL)
	{
		return (mat_fill(A, type));
	}
	else
	{
		return (NULL);
	}
}


/*
 *-----------------------------------------------------------------------------
 *       funct:  mat_fill
 *       desct:  form a special matrix
 *       given:  A = matrix, type = which kind of matrix
 *       retrn:  A
 *-----------------------------------------------------------------------------
 */

MATRIX mat_fill( MATRIX A, int type )
{
	int    i, j;

	switch (type)
	{
		case UNDEFINED:
			break;
		case ZERO_MATRIX:
		case UNIT_MATRIX:
			for (i=0; i<MatRow(A); i++)
			{
				for (j=0; j<MatCol(A); j++)
				{
					if (type == UNIT_MATRIX)
					{
						if (i==j)
						{
							A[i][j] = 1.0;
							continue;
						}
					}
					A[i][j] = 0.0;
				}
			}
			break;
	}
	return (A);
}


/*
 *-----------------------------------------------------------------------------
 *       funct:  mat_free
 *       desct:  free an allocated matrix
 *       given:  A = matrix
 *       retrn:  nothing <actually 0 = NULL A passed, 1 = normal exit>
 *-----------------------------------------------------------------------------
 */
int mat_free( MATRIX A )
{
	int i;

	if (A==NULL)
	{
		REprintf("\nAttempting to free a non-existent matrix in mat_free()\n");
		return (0);
	}

	for (i=0; i<MatRow(A); i++)
	{
		if(A[i]==NULL)
		{
			REprintf("\nAttempting to free a non-existent matrix row in mat_free()\n");
			return (0);
		}
		else
		{
			free(A[i]);
		}
	}

	free(Mathead(A));

	return (1);
}





/*
 *-----------------------------------------------------------------------------
 * file: materr.c
 * desc: matrix error handler
 * by: ko shu pui, Patrick
 * date: 24 nov 91 v0.1
 * revi:
 * ref:
 *       [1] Mary L.Boas, "Mathematical Methods in the Physical Science,"
 * John Wiley & Sons, 2nd Ed., 1983. Chap 3.
 *
 * [2] Kendall E.Atkinson, "An Introduction to Numerical Analysis,"
 * John Wiley & Sons, 1978.
 *
 *-----------------------------------------------------------------------------
 */

/*MATRIX mat_error( int errno )
{
	switch( errno )
	{
		case MAT_MALLOC:
			error("mat: malloc error\n" );
			break;
		case MAT_FNOTOPEN:
			error("mat: fileopen error\n" );
			break;
		case MAT_FNOTGETMAT:
			error("fgetmat: matrix read error\n");
			break;
	}

	return ((MATRIX)NULL);
  }*/


/*
 *-----------------------------------------------------------------------------
 * file: matinv.c
 * desc: matrix inversion
 * by: ko shu pui, Patrick
 * date: 24 nov 91 v0.1
 * revi: 14 may 92 v0.2
 * revi: Nov 25 1994 by J. Racine - added error checking, passing matrix
 *       rather than allocating
 * ref:
 *       [1] Mary L.Boas, "Mathematical Methods in the Physical Science,"
 * John Wiley & Sons, 2nd Ed., 1983. Chap 3.
 *
 * [2] Kendall E.Atkinson, "An Introduction to Numerical Analysis,"
 * John Wiley & Sons, 1978.
 *
 *-----------------------------------------------------------------------------
 */



/*
 *-----------------------------------------------------------------------------
 * funct:  mat_inv
 * desct:  find inverse of a matrix
 * given:  a = square matrix a, and C, return for inv(a)
 * retrn:  square matrix Inverse(A), C
 *   NULL = fails, singular matrix
 *-----------------------------------------------------------------------------
 */
MATRIX mat_inv( MATRIX a , MATRIX C)
{
	int i, j, n;
	int info = 0;
	int lwork = -1;
	double work_query = 0.0;
	double *A = NULL;
	int *ipiv = NULL;
	double *work = NULL;

#ifdef CONFORM_CHECK

	if(MatCol(a)!=MatCol(C))
	{
		error("\nUnconformable matrices in routine mat_inv(): Col(A)!=Col(B)\n");
	}

	if(MatRow(a)!=MatRow(C))
	{
		error("\nUnconformable matrices in routine mat_inv(): Row(A)!=Row(B)\n");
	}
#endif

	n = MatCol(a);

	A = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
	ipiv = (int *)malloc((size_t)n * sizeof(int));
	if ((A == NULL) || (ipiv == NULL))
		error("mat_inv: malloc error\n");

	/* Copy to LAPACK column-major dense storage. */
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
			A[i + j*n] = a[i][j];

	F77_CALL(dgetrf)(&n, &n, A, &n, ipiv, &info);
	if (info != 0)
	{
		free(A);
		free(ipiv);
		return NULL;
	}

	F77_CALL(dgetri)(&n, A, &n, ipiv, &work_query, &lwork, &info);
	if (info != 0)
	{
		free(A);
		free(ipiv);
		return NULL;
	}

	lwork = (int)work_query;
	if (lwork < n)
		lwork = n;
	work = (double *)malloc((size_t)lwork * sizeof(double));
	if (work == NULL)
		error("mat_inv: malloc error\n");

	F77_CALL(dgetri)(&n, A, &n, ipiv, work, &lwork, &info);
	free(work);
	free(ipiv);
	if (info != 0)
	{
		free(A);
		return NULL;
	}

	/* Copy back from LAPACK column-major dense storage. */
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
			C[i][j] = A[i + j*n];

	free(A);

	if(!isFiniteMatrix(C))
		return NULL;

	return (C);

}


/*
 *-----------------------------------------------------------------------------
 * file: matmul.c
 * desc: matrix multiplication
 * by: ko shu pui, Patrick
 * date: 24 nov 91 v0.1
 * revi: Nov 25 1994 by J. Racine - added error checking, passing matrix
 *       rather than allocating
 * ref:
 *       [1] Mary L.Boas, "Mathematical Methods in the Physical Science,"
 * John Wiley & Sons, 2nd Ed., 1983. Chap 3.
 *
 *-----------------------------------------------------------------------------
 */


/*
 *-----------------------------------------------------------------------------
 * funct:  mat_solve
 * desct:  solve A * X = B directly (prefer SPD via Cholesky, fallback LU)
 * given:  A square, B with matching rows
 * retrn:  X on success, NULL on singular/failed solve
 *-----------------------------------------------------------------------------
 */
MATRIX mat_solve( MATRIX A, MATRIX B, MATRIX X)
{
	int i, j;
	int info = 0;
	const int n = MatRow(A);
	const int nrhs = MatCol(B);
	double *Ac = NULL, *Bc = NULL;
	int *ipiv = NULL;

#ifdef CONFORM_CHECK
	if (MatCol(A) != n)
		error("\nUnconformable matrices in routine mat_solve(): A must be square\n");
	if (MatRow(B) != n)
		error("\nUnconformable matrices in routine mat_solve(): Row(A)!=Row(B)\n");
	if (MatRow(X) != n || MatCol(X) != nrhs)
		error("\nUnconformable matrices in routine mat_solve(): X dims mismatch\n");
#endif

	Ac = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
	Bc = (double *)malloc((size_t)n * (size_t)nrhs * sizeof(double));
	ipiv = (int *)malloc((size_t)n * sizeof(int));
	if ((Ac == NULL) || (Bc == NULL) || (ipiv == NULL))
		error("mat_solve: malloc error\n");

	/* Copy A,B to LAPACK column-major dense storage. */
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
			Ac[i + j*n] = A[i][j];
	for (j = 0; j < nrhs; j++)
		for (i = 0; i < n; i++)
			Bc[i + j*n] = B[i][j];

	F77_CALL(dgesv)(&n, &nrhs, Ac, &n, ipiv, Bc, &n, &info);
	free(Ac);
	free(ipiv);
	if (info != 0) {
		free(Bc);
		return NULL;
	}

	/* Copy solution back from LAPACK column-major dense storage. */
	for (j = 0; j < nrhs; j++)
		for (i = 0; i < n; i++)
			X[i][j] = Bc[i + j*n];
	free(Bc);

	if(!isFiniteMatrix(X))
		return NULL;

	return X;
}

int mat_is_nonsingular( MATRIX A )
{
	int i, j, info = 0;
	const int n = MatRow(A);
	double *Ac = NULL;
	int *ipiv = NULL;

#ifdef CONFORM_CHECK
	if (MatCol(A) != n)
		error("\nUnconformable matrices in routine mat_is_nonsingular(): A must be square\n");
#endif

	Ac = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
	ipiv = (int *)malloc((size_t)n * sizeof(int));
	if ((Ac == NULL) || (ipiv == NULL))
		error("mat_is_nonsingular: malloc error\n");

	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
			Ac[i + j*n] = A[i][j];

	F77_CALL(dgetrf)(&n, &n, Ac, &n, ipiv, &info);
	free(Ac);
	free(ipiv);

	return (info == 0);
}

double mat_inv00( MATRIX A, int *ok )
{
	int i, j, info = 0;
	int nrhs = 1;
	const int n = MatRow(A);
	double *Ac = NULL, *bc = NULL;
	int *ipiv = NULL;
	double v = 0.0;

#ifdef CONFORM_CHECK
	if (MatCol(A) != n)
		error("\nUnconformable matrices in routine mat_inv00(): A must be square\n");
#endif

	Ac = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
	bc = (double *)calloc((size_t)n, sizeof(double));
	ipiv = (int *)malloc((size_t)n * sizeof(int));
	if ((Ac == NULL) || (bc == NULL) || (ipiv == NULL))
		error("mat_inv00: malloc error\n");

	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
			Ac[i + j*n] = A[i][j];
	bc[0] = 1.0;

	F77_CALL(dgesv)(&n, &nrhs, Ac, &n, ipiv, bc, &n, &info);
	if (info != 0) {
		*ok = 0;
		free(Ac);
		free(bc);
		free(ipiv);
		return 0.0;
	}

	v = bc[0];
	*ok = isfinite(v) ? 1 : 0;

	free(Ac);
	free(bc);
	free(ipiv);

	return v;
}


/*
 *-----------------------------------------------------------------------------
 * file: matsolve.c
 * desc: solve linear equations
 * by: ko shu pui, Patrick
 * date: 24 nov 91 v0.1
 * revi: 14 may 92 v0.2
 * ref:
 *       [1] Mary L.Boas, "Mathematical Methods in the Physical Science,"
 * John Wiley & Sons, 2nd Ed., 1983. Chap 3.
 *
 * [2] Kendall E.Atkinson, "An Introduction to Numerical Analysis,"
 * John Wiley & Sons, 1978.
 *
 *-----------------------------------------------------------------------------
 */


int isFiniteMatrix(MATRIX A){
	int i, j;

	const int nc = MatCol(A);
	const int nr = MatRow(A);

  for(i = 0; i < nr; i++){
    for(j = 0; j < nc; j++){
      if(!isfinite(A[i][j]))
        return 0;
    }
  }
  return 1;
}
