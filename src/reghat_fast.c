#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

/*
 * Exact compiled counterpart of the R loop in
 * .npreghat_exact_lp_matrix_from_kernel_weights().  The caller limits this
 * route to R's BLAS-backed matprod modes.  Matrix products, LAPACK condition
 * checks, ridge increments, and final row construction deliberately retain
 * the incumbent operation order.
 */

static int np_matrix_dims(SEXP x, int *nr, int *nc)
{
  SEXP dim = getAttrib(x, R_DimSymbol);
  if((TYPEOF(dim) != INTSXP) || (XLENGTH(dim) != 2))
    return 0;
  *nr = INTEGER(dim)[0];
  *nc = INTEGER(dim)[1];
  return (*nr >= 0) && (*nc >= 0);
}

/*
 * Exact compiled counterpart of
 *
 *   sweep(t(kw), 1L, denominator, "/")
 *
 * Callers retain ownership of the denominator calculation and any zero-floor
 * policy.  This helper changes only the transpose/allocation traversal.  The
 * scalar division for every output element is therefore identical to the R
 * route, while small tiles keep both the strided input and contiguous output
 * writes resident in cache.
 */
SEXP C_np_lc_hat_normalize(SEXP kw, SEXP denominator)
{
  const int tile = 32;
  int ntrain = 0, neval = 0;
  const double *kw_p = NULL;
  const double *denominator_p = NULL;
  double *out_p = NULL;
  SEXP out = R_NilValue;

  if((TYPEOF(kw) != REALSXP) || (TYPEOF(denominator) != REALSXP) ||
     !np_matrix_dims(kw, &ntrain, &neval) ||
     (XLENGTH(denominator) != (R_xlen_t)neval))
    error("invalid LC hat normalization input");

  kw_p = REAL(kw);
  denominator_p = REAL(denominator);
  out = PROTECT(allocMatrix(REALSXP, neval, ntrain));
  out_p = REAL(out);

  for(int jb = 0; jb < neval; jb += tile){
    const int jend = (jb + tile < neval) ? jb + tile : neval;
    for(int ib = 0; ib < ntrain; ib += tile){
      const int iend = (ib + tile < ntrain) ? ib + tile : ntrain;
      for(int i = ib; i < iend; i++)
        for(int j = jb; j < jend; j++)
          out_p[j + (size_t)neval*(size_t)i] =
            kw_p[i + (size_t)ntrain*(size_t)j]/denominator_p[j];
    }
    R_CheckUserInterrupt();
  }

  UNPROTECT(1);
  return out;
}

static int np_reghat_solve_system(const int nterms,
                                  const double * const matrix,
                                  const double * const rhs,
                                  double * const matrix_work,
                                  double * const solution,
                                  int * const pivot,
                                  double * const condition_work,
                                  int * const condition_iwork)
{
  const char norm = '1';
  const int nrhs = 1;
  int info = 0;
  double anorm = 0.0;
  double rcond = 0.0;

  memcpy(matrix_work, matrix,
         (size_t)nterms*(size_t)nterms*sizeof(double));
  memcpy(solution, rhs, (size_t)nterms*sizeof(double));
  anorm = F77_CALL(dlange)(&norm, &nterms, &nterms, matrix_work, &nterms,
                           condition_work FCONE);
  F77_CALL(dgesv)(&nterms, &nrhs, matrix_work, &nterms, pivot,
                  solution, &nterms, &info);
  if(info != 0)
    return 0;

  F77_CALL(dgecon)(&norm, &nterms, matrix_work, &nterms, &anorm, &rcond,
                   condition_work, condition_iwork, &info FCONE);
  if((info != 0) || !isfinite(rcond) || (rcond < DBL_EPSILON))
    return 0;
  for(int term = 0; term < nterms; term++)
    if(!isfinite(solution[term]))
      return 0;
  return 1;
}

SEXP C_np_reghat_lp_matrix_fast(SEXP kw, SEXP wtrain, SEXP weval)
{
  int ntrain = 0, neval = 0, kw_neval = 0;
  int wtrain_n = 0, nterms = 0, weval_n = 0, weval_p = 0;
  double *weighted_design = NULL;
  double *gram = NULL;
  double *gram_work = NULL;
  double *rhs = NULL;
  double *solution = NULL;
  double *prediction = NULL;
  double *condition_work = NULL;
  int *pivot = NULL;
  int *condition_iwork = NULL;
  SEXP out = R_NilValue;

  if((TYPEOF(kw) != REALSXP) || (TYPEOF(wtrain) != REALSXP) ||
     (TYPEOF(weval) != REALSXP) ||
     !np_matrix_dims(kw, &ntrain, &kw_neval) ||
     !np_matrix_dims(wtrain, &wtrain_n, &nterms) ||
     !np_matrix_dims(weval, &weval_n, &weval_p) ||
     (ntrain <= 0) || (kw_neval <= 0) || (nterms <= 0) ||
     (wtrain_n != ntrain) || (weval_n != kw_neval) ||
     (weval_p != nterms))
    return R_NilValue;

  neval = kw_neval;
  if(((size_t)nterms > SIZE_MAX/sizeof(double)) ||
     ((size_t)ntrain > SIZE_MAX/((size_t)nterms*sizeof(double))) ||
     ((size_t)nterms > SIZE_MAX/((size_t)nterms*sizeof(double))) ||
     ((size_t)nterms > SIZE_MAX/(4*sizeof(double))) ||
     ((size_t)neval > SIZE_MAX/(size_t)ntrain))
    return R_NilValue;

  weighted_design = (double *)R_alloc((size_t)ntrain*(size_t)nterms,
                                      sizeof(double));
  gram = (double *)R_alloc((size_t)nterms*(size_t)nterms, sizeof(double));
  gram_work = (double *)R_alloc((size_t)nterms*(size_t)nterms, sizeof(double));
  rhs = (double *)R_alloc((size_t)nterms, sizeof(double));
  solution = (double *)R_alloc((size_t)nterms, sizeof(double));
  prediction = (double *)R_alloc((size_t)ntrain, sizeof(double));
  condition_work = (double *)R_alloc((size_t)4*(size_t)nterms, sizeof(double));
  pivot = (int *)R_alloc((size_t)nterms, sizeof(int));
  condition_iwork = (int *)R_alloc((size_t)nterms, sizeof(int));

  out = PROTECT(allocMatrix(REALSXP, neval, ntrain));
  for(int j = 0; j < neval; j++){
    const double * const weights = REAL(kw) + (size_t)j*(size_t)ntrain;
    const double alpha = 1.0;
    const double beta = 0.0;
    const int one = 1;
    const double epsilon = 1.0/(double)ntrain;
    double nepsilon = 0.0;
    const char trans_t = 'T';
    const char trans_n = 'N';

    if((j == 0) || (j + 1 == neval) || ((j % 32) == 0))
      R_CheckUserInterrupt();

    for(int term = 0; term < nterms; term++){
      const double * const src = REAL(wtrain) + (size_t)term*(size_t)ntrain;
      double * const dst = weighted_design + (size_t)term*(size_t)ntrain;
      for(int i = 0; i < ntrain; i++)
        dst[i] = src[i]*weights[i];
    }

    F77_CALL(dgemm)(&trans_t, &trans_n,
                    &nterms, &nterms, &ntrain,
                    &alpha, REAL(wtrain), &ntrain,
                    weighted_design, &ntrain,
                    &beta, gram, &nterms FCONE FCONE);

    for(int term = 0; term < nterms; term++)
      rhs[term] = REAL(weval)[j + (size_t)neval*(size_t)term];

    if(!np_reghat_solve_system(nterms, gram, rhs, gram_work, solution,
                               pivot, condition_work, condition_iwork)){
      do {
        for(int term = 0; term < nterms; term++)
          gram[term + (size_t)nterms*(size_t)term] += epsilon;
        nepsilon += epsilon;
        R_CheckUserInterrupt();
      } while(!np_reghat_solve_system(nterms, gram, rhs, gram_work, solution,
                                      pivot, condition_work, condition_iwork));

      {
        double denom = gram[0];
        if(!isfinite(denom) || (fabs(denom) < DBL_MIN))
          denom = DBL_MIN;
        solution[0] *= 1.0 + nepsilon/denom;
      }
    }

    F77_CALL(dgemv)(&trans_n, &ntrain, &nterms, &alpha,
                    REAL(wtrain), &ntrain, solution, &one,
                    &beta, prediction, &one FCONE);
    for(int i = 0; i < ntrain; i++)
      REAL(out)[j + (size_t)neval*(size_t)i] = weights[i]*prediction[i];
  }

  UNPROTECT(1);
  return out;
}
