#include <stdint.h>
#include <mpi.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Random.h>
#include <Rinternals.h>

#define CHAR2(x) ((char *)CHAR(x))

/* #define COMM_MAXSIZE 10 */

int mpi_errhandler(int errcode);
int erreturn(int errcode);

MPI_Datatype mpitype(SEXP sexp_type); 

// void mystrcpy(char *new_str, char *old_str, int size);

SEXP AsInt (int n);

struct Dblint {
	double x;
	int rank;

};

// char *acopy_string2(const char *in);
// char *charsxp2char(SEXP x);

// #define CallocCharBuf2(n) 	(char *) R_chk_calloc2((size_t) ((n)+1), sizeof(char))
// extern void *R_chk_calloc2(size_t, size_t);

