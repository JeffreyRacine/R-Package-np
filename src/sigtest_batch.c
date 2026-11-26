#include <R.h>
#include <Rinternals.h>
#include "regression_contrast.h"

extern void np_progress_fit_heartbeat(void);
extern int np_sigtest_prepare_contrast(int, const double *, const double *,
                                      NPContrastNumber *, int *);

/* Each kernel column belongs to one training row. Residual preparation
 * uses centered response differences before any accumulation, including
 * when the response units have a large common offset. */
SEXP C_np_sigtest_residual_block(SEXP kernel, SEXP response, SEXP rows)
{
  const int n = nrows(response), draws = ncols(response), m = ncols(kernel);
  if(TYPEOF(kernel) != REALSXP || TYPEOF(response) != REALSXP ||
     !isMatrix(kernel) || !isMatrix(response) || nrows(kernel) != n ||
     n < 1 || draws < 1 || draws > 64 || m < 1 || TYPEOF(rows) != INTSXP || XLENGTH(rows) != m)
    error("invalid categorical significance residual block");
  long double *weight = R_allocLD((size_t)n);
  SEXP out = PROTECT(allocMatrix(REALSXP, m, draws));
  for(int i=0;i<m;++i) {
    const int row = INTEGER(rows)[i];
    if(row < 1 || row > n) error("invalid training row in significance test");
    const int self = row - 1;
    const double *w = REAL(kernel)+(size_t)i*n;
    NPResidualOffDiagonal geometry = {0};
    for(int j=0;j<n;++j)
      if(j != self) np_residual_offdiag_add(&geometry, w[j], 0.0L);
    NPResidualNormalization normalization;
    if(np_residual_offdiag_normalization(&geometry,1.0L,&normalization)
       != NP_RESIDUAL_IDENTIFIED) { UNPROTECT(1); return R_NilValue; }
    for(int j=0;j<n;++j)
      weight[j] = j == self ? 0.0L : np_residual_normalized_weight(&normalization,w[j]);
    for(int b=0;b<draws;++b) {
      const double *y = REAL(response)+(size_t)b*n;
      long double numerator = 0.0L;
      for(int j=0;j<n;++j)
        numerator += weight[j]*((long double)y[self]-y[j]);
      const double z = (double)numerator;
      if(!R_FINITE(z)) { UNPROTECT(1); return R_NilValue; }
      REAL(out)[i+(size_t)m*b] = z;
    }
    if((i & 7)==0) { R_CheckUserInterrupt(); np_progress_fit_heartbeat(); }
  }
  UNPROTECT(1);
  return out;
}

/* Reuse the canonical two-double coefficient, anchor and scaled variance
 * arithmetic across a bounded response tile. Never square raw residual units. */
SEXP C_np_sigtest_contrast_block(SEXP left, SEXP right, SEXP response,
                                  SEXP standardized, SEXP scale)
{
  const int n=nrows(response), draws=ncols(response), m=ncols(left);
  if(TYPEOF(left)!=REALSXP || TYPEOF(right)!=REALSXP || TYPEOF(response)!=REALSXP ||
     TYPEOF(standardized)!=REALSXP || TYPEOF(scale)!=REALSXP ||
     !isMatrix(left) || !isMatrix(right) || !isMatrix(response) || !isMatrix(standardized) ||
     nrows(left)!=n || nrows(right)!=n || ncols(right)!=m ||
     nrows(standardized)!=n || ncols(standardized)!=draws || XLENGTH(scale)!=draws ||
     n<1 || draws<1 || draws>64 || m<1) error("invalid categorical significance contrast block");
  NPContrastNumber *w=(NPContrastNumber *)R_alloc(n,sizeof(*w));
  SEXP out=PROTECT(allocVector(VECSXP,3));
  SEXP g=PROTECT(allocMatrix(REALSXP,m,draws));
  SEXP se=PROTECT(allocMatrix(REALSXP,m,draws));
  SEXP structural=PROTECT(allocVector(LGLSXP,m));
  for(int i=0;i<m;++i) {
    int anchor;
    if(!np_sigtest_prepare_contrast(n,REAL(left)+(size_t)n*i,
        REAL(right)+(size_t)n*i,w,&anchor)) {UNPROTECT(4);return R_NilValue;}
    LOGICAL(structural)[i]=anchor<0;
    for(int b=0;b<draws;++b) {
      const size_t pos=i+(size_t)b*m;
      if(anchor<0) { REAL(g)[pos]=REAL(se)[pos]=0.0;continue; }
      const double *y=REAL(response)+(size_t)n*b;
      const double *z=REAL(standardized)+(size_t)n*b;
      NPRegressionContrastAccumulator sum=np_regression_contrast_begin(y[anchor],1);
      for(int j=0;j<n;++j)
        if(j!=anchor) np_regression_contrast_add(&sum,w[j],y[j],0,NP_RESIDUAL_IDENTIFIED,z[j]);
      const NPRegressionContrastResult contrast=np_regression_contrast_finish(&sum,
        NP_RESIDUAL_IDENTIFIED,z[anchor],0);
      if(contrast.information!=NP_RESIDUAL_IDENTIFIED) {UNPROTECT(4);return R_NilValue;}
      REAL(g)[pos]=contrast.contrast;
      REAL(se)[pos]=contrast.standard_error*REAL(scale)[b];
    }
    if((i & 7)==0) {R_CheckUserInterrupt();np_progress_fit_heartbeat();}
  }
  SET_VECTOR_ELT(out,0,g); SET_VECTOR_ELT(out,1,se); SET_VECTOR_ELT(out,2,structural);
  UNPROTECT(4);
  return out;
}
