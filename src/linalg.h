/*
*-----------------------------------------------------------------------------
* file: linalg.h
* desc: matrix mathematics header file
* by: ko shu pui, Patrick
* date: 24 nov 91 v0.1b
* revi: Nov 25 1994 by J. Racine. Changed allocation, added bounds checking
* ref:
* [1] Mary L.Boas, "Mathematical Methods in the Physical Science,"
* John Wiley & Sons, 2nd Ed., 1983. Chap 3.
*
*-----------------------------------------------------------------------------
*/

/* For conformability checking, define CONFORM_CHECK. For run time
 for most apps, do not define to reduce overhead: Added by Racine */

/*#define CONFORM_CHECK*/

/*
*-----------------------------------------------------------------------------
* internal matrix structure
*-----------------------------------------------------------------------------
*/
typedef struct {
 int row;
 int col;
 } MATHEAD;

typedef struct {
 MATHEAD head;
 /*
 * only the starting address of the following will be
 * returned to the C programmer, like malloc() concept
 */
 double *matrix;
 } MATBODY;

typedef double **MATRIX;

#define Mathead(a) ((MATHEAD *)((MATHEAD *)(a) - 1))
#define MatRow(a) (Mathead(a)->row)
#define MatCol(a) (Mathead(a)->col)

/*
*----------------------------------------------------------------------------
* mat_errors definitions
*----------------------------------------------------------------------------
*/
#define MAT_MALLOC 1
#define MAT_FNOTOPEN 2
#define MAT_FNOTGETMAT 3

/*
*----------------------------------------------------------------------------
* matrix types
*----------------------------------------------------------------------------
*/
#define UNDEFINED -1
#define ZERO_MATRIX 0
#define UNIT_MATRIX 1

/*MATRIX mat_error( int errno );*/
MATRIX _mat_creat( int row, int col );
MATRIX mat_creat( int row, int col, int type );
MATRIX mat_fill( MATRIX A, int type );
int mat_free( MATRIX A );
MATRIX mat_solve( MATRIX A, MATRIX B, MATRIX X );
double mat_inv00( MATRIX A, int *ok );
int mat_is_nonsingular( MATRIX A );
MATRIX mat_inv( MATRIX a , MATRIX C );
int isFiniteMatrix(MATRIX A);
