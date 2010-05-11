#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

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
          printf("mat: malloc error: exit(EXIT_SUCCESS)\n" );
          exit(EXIT_SUCCESS);
	}

	for (i=0; i<row; i++)
	{
		if ((*((double **)(&mat->matrix) + i) = (double *)malloc(sizeof(double) * col)) == NULL)
		{
			printf("mat: malloc error: exit(EXIT_SUCCESS)\n" );
			exit(EXIT_SUCCESS);
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
		printf("\nAttempting to free a non-existent matrix in mat_free()\n");
		return (0);
	}

	for (i=0; i<MatRow(A); i++)
	{
		if((A[i]==NULL))
		{
			printf("\nAttempting to free a non-existent matrix row in mat_free()\n");
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
 *       funct:  mat_copy
 *       desct:  duplicate a matrix
 *       given:  A = matrix to duplicated
 *       retrn:  C = A
 *       comen:
 *-----------------------------------------------------------------------------
 */
MATRIX mat_copy( MATRIX A, MATRIX C )
{
	int i, j;

	for (i=0; i<MatRow(A); i++)
	{
		for (j=0; j<MatCol(A); j++)
		{
			C[i][j] = A[i][j];
		}
	}
	return (C);
}






/*
 *-----------------------------------------------------------------------------
 * file: matdet.c
 * desc: determinant calculations
 * by: ko shu pui, Patrick
 * date: 21 may 92 v0.3
 * revi:
 * ref:
 *       [1] Mary L.Boas, "Mathematical Methods in the Physical Science,"
 * John Wiley & Sons, 2nd Ed., 1983. Chap 3.
 *
 *-----------------------------------------------------------------------------
 */

static double signa[2] = {1.0, -1.0};


/*
 *-----------------------------------------------------------------------------
 * funct:  mat_det
 * desct:  find determinant
 * given:  A = matrix
 * retrn:  the determinant of A
 * comen:
 *-----------------------------------------------------------------------------
 */
double mat_det( MATRIX a )
{
	MATRIX  A, P;
	int  j;
	int i, n;
	double  result;

	n = MatRow(a);
	A = mat_creat(MatRow(a), MatCol(a), UNDEFINED);
	A = mat_copy(a, A);
	P = mat_creat(n, 1, UNDEFINED);



/*
 * take a LUP-decomposition
 */
	i = mat_lu(A, P);
	switch (i)
	{


/*
 * case for singular matrix
 */
		case -1:
			result = 0.0;
			break;



/*
 * normal case: |A| = |L||U||P|
 * |L| = 1,
 * |U| = multiplication of the diagonal
 * |P| = +-1
 */
		default:
			result = 1.0;
			for (j=0; j<MatRow(A); j++)
			{
				result *= A[(int)P[j][0]][j];
			}
			result *= signa[i%2];
			break;
	}

	mat_free(A);
	mat_free(P);
	return (result);
}


/*
 *-----------------------------------------------------------------------------
 * file: matdump.c
 * desc: matrix mathematics - object dump
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



/*
 *-----------------------------------------------------------------------------
 * funct:  mat_dumpf
 *   desct:  dump a matrix with format string to standard output
 * given:  A = matrix to dumped
 * retrn:  nothing
 * comen:  matrix a dumped to standard output
 *-----------------------------------------------------------------------------
 */
MATRIX mat_dumpf(MATRIX A, char *s)
{
	return (mat_fdumpf(A, s, stdout));
}


MATRIX mat_fdumpf( MATRIX A, char *s, FILE *fp )
{
	int  i, j;

	for (i=0; i<MatRow(A); i++)
	{
		for (j=0; j<MatCol(A); j++)
		{
			fprintf( fp, s, A[i][j] );
		}
		fprintf( fp, "\n" );
	}

	return (A);
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

MATRIX mat_error( int errno )
{
	switch( errno )
	{
		case MAT_MALLOC:
			printf("mat: malloc error: exit(EXIT_SUCCESS)\n" );
			exit(EXIT_SUCCESS);
			break;
		case MAT_FNOTOPEN:
			printf("mat: fileopen error\n" );
			break;
		case MAT_FNOTGETMAT:
			printf("fgetmat: matrix read error\n");
			break;
	}

	return ((MATRIX)NULL);
}


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
	MATRIX  A, B, P;
	int i, n;

#ifdef CONFORM_CHECK

	if(MatCol(a)!=MatCol(C))
	{
		printf("\nUnconformable matrices in routine mat_inv(): Col(A)!=Col(B): exit(EXIT_SUCCESS)\n");
		exit(EXIT_SUCCESS);
	}

	if(MatRow(a)!=MatRow(C))
	{
		printf("\nUnconformable matrices in routine mat_inv(): Row(A)!=Row(B): exit(EXIT_SUCCESS)\n");
		exit(EXIT_SUCCESS);
	}
#endif

	n = MatCol(a);
	A = mat_creat(MatRow(a), MatCol(a), UNDEFINED);
	A = mat_copy(a, A);
	B = mat_creat( n, 1, UNDEFINED );
	P = mat_creat( n, 1, UNDEFINED );



/*
 * - LU-decomposition -
 * also check for singular matrix
 */

	if (mat_lu(A, P) == -1)
	{

		mat_free(A);
		mat_free(B);
		mat_free(P);

		return (NULL);

	}
	else
	{

/* Bug??? was still mat_backsubs1 even when singular??? */

		for (i=0; i<n; i++)
		{
			mat_fill(B, ZERO_MATRIX);
			B[i][0] = 1.0;
			mat_backsubs1( A, B, C, P, i );
		}

		mat_free(A);
		mat_free(B);
		mat_free(P);

		return (C);

	}

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
 * funct:  mat_mul
 * desct:  multiplication of two matrices
 * given:  A, B = compatible matrices to be multiplied
 * retrn:  NULL if malloc() fails
 *   else allocated matrix of A * B
 * comen:
 *-----------------------------------------------------------------------------
 */
MATRIX mat_mul( MATRIX A, MATRIX B , MATRIX C)
{
	int  i, j, k;

#ifdef CONFORM_CHECK

	if(MatCol(A)!=MatRow(B))
	{
		printf("\nUnconformable matrices in routine mat_mul(): Col(A)!=Row(B) (%d/%d): exit(EXIT_SUCCESS)\n", MatCol(A),MatRow(B));
		exit(EXIT_SUCCESS);
	}
	if(MatRow(A)!=MatRow(C))
	{
		printf("\nUnconformable matrices in routine mat_mul(): Row(A)!=Row(C) (%d/%d): exit(EXIT_SUCCESS)\n", MatRow(A), MatRow(C));
		exit(EXIT_SUCCESS);
	}
	if(MatCol(B)!=MatCol(C))
	{
		printf("\nUnconformable matrices in routine mat_mul(): Col(B)!=Col(C) (%d/%d): exit(EXIT_SUCCESS)\n", MatCol(B), MatCol(C));
		exit(EXIT_SUCCESS);
	}
#endif

	for (i=0; i<MatRow(A); i++)
	{
		for (j=0; j<MatCol(B); j++)
		{
			for (k=0, C[i][j]=0.0; k<MatCol(A); k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return (C);
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



/*
 *-----------------------------------------------------------------------------
 * funct:  mat_lu
 * desct:  in-place LU decomposition with partial pivoting
 * given:  !! A = square matrix (n x n) !ATTENTION! see comment
 *   P = permutation vector (n x 1)
 * retrn:  number of permutation performed
 *   -1 means suspected singular matrix
 * comen:  A will be overwritten to be a LU-composite matrix
 *
 * note: the LU decomposed may NOT be equal to the LU of
 *   the original matrix a. But equal to the LU of the
 *   rows interchanged matrix.
 *-----------------------------------------------------------------------------
 */
int mat_lu( MATRIX A, MATRIX P )
{
	int i, j, k, n;
	int maxi;
	double tmp;
	double  c, c1;
	int p;

	n = MatCol(A);

	for (p=0,i=0; i<n; i++)
	{
		P[i][0] = i;
	}

	for (k=0; k<n; k++)
	{


/*
 * --- partial pivoting ---
 */
		for (i=k, maxi=k, c=0.0; i<n; i++)
		{
			c1 = fabs( A[(int)P[i][0]][k] );
			if (c1 > c)
			{
				c = c1;
				maxi = i;
			}
		}



/*
 * row exchange, update permutation vector
 */
		if (k != maxi)
		{
			p++;
			tmp = P[k][0];
			P[k][0] = P[maxi][0];
			P[maxi][0] = tmp;
		}



/*
 * Test for singular matrix (J. Racine, Nov 24 2000)
 */
		if ( fabs(A[(int)P[k][0]][k]) > 0.0 )
		{

			for (i=k+1; i<n; i++)
			{


/*
 * --- calculate m(i,j) ---
 */
				A[(int)P[i][0]][k] = A[(int)P[i][0]][k] / A[(int)P[k][0]][k];



/*
 * --- elimination ---
 */
				for (j=k+1; j<n; j++)
				{
					A[(int)P[i][0]][j] -= A[(int)P[i][0]][k] * A[(int)P[k][0]][j];
				}
			}

		}
		else
		{

			return (-1);

		}

	}

	return (p);

}


/*
 *-----------------------------------------------------------------------------
 * funct:  mat_backsubs1
 * desct:  back substitution
 * given:  A = square matrix A (LU composite)
 *   !! B = column matrix B (attention!, see comen)
 *   !! X = place to put the result of X
 *   P = Permutation vector (after calling mat_lu)
 *   xcol = column of x to put the result
 * retrn:  column matrix X (of AX = B)
 * comen:  B will be overwritten
 *-----------------------------------------------------------------------------
 */
MATRIX mat_backsubs1( MATRIX A, MATRIX B, MATRIX X, MATRIX P, int xcol )
{
	int i, j, k, n;
	double  sum;

	n = MatCol(A);

	for (k=0; k<n; k++)
	{
		for (i=k+1; i<n; i++)
			B[(int)P[i][0]][0] -= A[(int)P[i][0]][k] * B[(int)P[k][0]][0];
	}

	X[n-1][xcol] = B[(int)P[n-1][0]][0] / A[(int)P[n-1][0]][n-1];
	for (k=n-2; k>=0; k--)
	{
		sum = 0.0;
		for (j=k+1; j<n; j++)
		{
			sum += A[(int)P[k][0]][j] * X[j][xcol];
		}
		X[k][xcol] = (B[(int)P[k][0]][0] - sum) / A[(int)P[k][0]][k];
	}

	return (X);
}




/*
 *-----------------------------------------------------------------------------
 * file: matsub.c
 * desc: matrix subtraction
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
 * file: matsubx.c
 * desc: find submatrix
 * by: ko shu pui, Patrick
 * date: 24 may 92 v0.4
 * revi:
 * ref:
 *       [1] Mary L.Boas, "Mathematical Methods in the Physical Science,"
 * John Wiley & Sons, 2nd Ed., 1983. Chap 3.
 *
 *-----------------------------------------------------------------------------
 */
/*
 *-----------------------------------------------------------------------------
 * file: mattran.c
 * desc: matrix mathematics
 * by: ko shu pui, Patrick
 * date: v0.1 - 24 nov 91
 * revi: v0.2 - 14 may 92
 * ref:
 *       [1] Mary L.Boas, "Mathematical Methods in the Physical Science,"
 * John Wiley & Sons, 2nd Ed., 1983. Chap 3.
 *
 *-----------------------------------------------------------------------------
 */




