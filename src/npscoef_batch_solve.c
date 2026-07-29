#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include <R.h>
#include <R_ext/Arith.h>
#include <R_ext/Lapack.h>
#include <Rinternals.h>

#include "headers.h"

/*
 * Solve the common zero-ridge npscoef moment systems in one native entry.
 *
 * The return value is a p x neval solution matrix when every row is finite,
 * nonsingular, and passes solve.default's reciprocal-condition threshold.
 * Any inadmissible row returns NULL so the R caller can execute its complete
 * historical ridge loop.  The input arrays are never modified.
 */
SEXP C_np_npscoef_batch_zero_solve(SEXP tww_r, SEXP tyw_r)
{
  SEXP tww_dim;
  SEXP tyw_dim;
  SEXP answer;
  const double *tww;
  const double *tyw;
  double *solution;
  double *gram;
  double *rhs;
  double *condition_work;
  int *ipiv;
  size_t gram_count;
  size_t solution_count;
  int p;
  int neval;
  int one = 1;
  int row;

  if((TYPEOF(tww_r) != REALSXP) || (TYPEOF(tyw_r) != REALSXP))
    error("internal npscoef batch solve requires double arrays");

  tww_dim = getAttrib(tww_r, R_DimSymbol);
  tyw_dim = getAttrib(tyw_r, R_DimSymbol);
  if((TYPEOF(tww_dim) != INTSXP) || (XLENGTH(tww_dim) != 3) ||
     (TYPEOF(tyw_dim) != INTSXP) || (XLENGTH(tyw_dim) != 2))
    error("internal npscoef batch solve received invalid dimensions");

  p = INTEGER(tww_dim)[0];
  neval = INTEGER(tww_dim)[2];
  if((p <= 2) || (neval <= 0) ||
     (INTEGER(tww_dim)[1] != p) ||
     (INTEGER(tyw_dim)[0] != p) ||
     (INTEGER(tyw_dim)[1] != neval))
    error("internal npscoef batch solve received incompatible dimensions");

  if(((size_t)p > SIZE_MAX/(size_t)p) ||
     ((size_t)p > SIZE_MAX/(size_t)neval) ||
     ((size_t)p > SIZE_MAX/(size_t)4))
    error("internal npscoef batch solve size overflow");
  gram_count = (size_t)p*(size_t)p;
  solution_count = (size_t)p*(size_t)neval;
  if((R_xlen_t)gram_count > R_XLEN_T_MAX ||
     (R_xlen_t)solution_count > R_XLEN_T_MAX)
    error("internal npscoef batch solve R vector overflow");

  answer = PROTECT(allocMatrix(REALSXP, p, neval));
  solution = REAL(answer);
  tww = REAL(tww_r);
  tyw = REAL(tyw_r);
  gram = (double *)R_alloc(gram_count, sizeof(double));
  rhs = (double *)R_alloc((size_t)p, sizeof(double));
  ipiv = (int *)R_alloc((size_t)p, sizeof(int));
  condition_work = (double *)R_alloc((size_t)4*(size_t)p, sizeof(double));

  for(row = 0; row < neval; row++){
    const char norm = '1';
    const double *gram_source = tww + (size_t)row*gram_count;
    const double *rhs_source = tyw + (size_t)row*(size_t)p;
    double anorm;
    double rcond = 0.0;
    int info = 0;
    size_t j;

    if((row & 63) == 0)
      np_check_user_interrupt();

    for(j = 0; j < gram_count; j++)
      if(!R_FINITE(gram_source[j])){
        UNPROTECT(1);
        return R_NilValue;
      }
    for(j = 0; j < (size_t)p; j++)
      if(!R_FINITE(rhs_source[j])){
        UNPROTECT(1);
        return R_NilValue;
      }

    anorm = F77_CALL(dlange)(&norm, &p, &p, gram_source, &p,
                             (double *)NULL FCONE);
    if(!R_FINITE(anorm) || (anorm <= 0.0)){
      UNPROTECT(1);
      return R_NilValue;
    }

    memcpy(gram, gram_source, gram_count*sizeof(double));
    memcpy(rhs, rhs_source, (size_t)p*sizeof(double));
    F77_CALL(dgesv)(&p, &one, gram, &p, ipiv, rhs, &p, &info);
    if(info != 0){
      UNPROTECT(1);
      return R_NilValue;
    }

    for(j = 0; j < gram_count; j++)
      if(!R_FINITE(gram[j])){
        UNPROTECT(1);
        return R_NilValue;
      }

    F77_CALL(dgecon)(&norm, &p, gram, &p, &anorm, &rcond,
                     condition_work, ipiv, &info FCONE);
    if((info != 0) || !R_FINITE(rcond) || (rcond < DBL_EPSILON)){
      UNPROTECT(1);
      return R_NilValue;
    }

    for(j = 0; j < (size_t)p; j++){
      if(!R_FINITE(rhs[j])){
        UNPROTECT(1);
        return R_NilValue;
      }
      solution[j + (size_t)row*(size_t)p] = rhs[j];
    }
  }

  UNPROTECT(1);
  return answer;
}

/*
 * Project each stable batch-solution column through its corresponding
 * hoisted basis row.  Basis terms remain in ascending order while the
 * register-local loop avoids both R slice/reshape work and tiny-BLAS dispatch.
 */
SEXP C_np_npscoef_batch_project(SEXP theta_r, SEXP wz_r)
{
  SEXP theta_dim;
  SEXP wz_dim;
  SEXP answer;
  const double *theta;
  const double *wz;
  double *coef;
  int ncoef;
  int neval;
  int nbasis;
  int pcoef;
  int row;
  int col;
  int k;

  if((TYPEOF(theta_r) != REALSXP) || (TYPEOF(wz_r) != REALSXP))
    error("internal npscoef batch projection requires double matrices");

  theta_dim = getAttrib(theta_r, R_DimSymbol);
  wz_dim = getAttrib(wz_r, R_DimSymbol);
  if((TYPEOF(theta_dim) != INTSXP) || (XLENGTH(theta_dim) != 2) ||
     (TYPEOF(wz_dim) != INTSXP) || (XLENGTH(wz_dim) != 2))
    error("internal npscoef batch projection received invalid dimensions");

  ncoef = INTEGER(theta_dim)[0];
  neval = INTEGER(theta_dim)[1];
  nbasis = INTEGER(wz_dim)[1];
  if((ncoef <= 0) || (neval <= 0) || (nbasis <= 0) ||
     (INTEGER(wz_dim)[0] != neval) || ((ncoef % nbasis) != 0))
    error("internal npscoef batch projection received incompatible dimensions");
  pcoef = ncoef/nbasis;

  answer = PROTECT(allocMatrix(REALSXP, pcoef, neval));
  theta = REAL(theta_r);
  wz = REAL(wz_r);
  coef = REAL(answer);

  for(row = 0; row < neval; row++){
    if((row > 0) && ((row & 4095) == 0))
      np_check_user_interrupt();
    for(col = 0; col < pcoef; col++){
      const double *theta_col =
        theta + (size_t)row*(size_t)ncoef + (size_t)col*(size_t)nbasis;
      double value = 0.0;
      for(k = 0; k < nbasis; k++)
        value += wz[row + (size_t)k*(size_t)neval] * theta_col[k];
      coef[col + (size_t)row*(size_t)pcoef] = value;
    }
  }

  UNPROTECT(1);
  return answer;
}
