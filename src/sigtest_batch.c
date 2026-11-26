#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "sigtest_lp_batch.h"
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

/* The two producers hand us the ordinary fit's actual off-diagonal
 * influences and compensated categorical coefficients. Neither the
 * bandwidth, the local solve nor its rank/ridge policy is reconstructed. */
NPSigtestLPBatch *np_sigtest_lp_begin(SEXP responses, int n)
{
  if (TYPEOF(responses) != REALSXP || !isMatrix(responses) ||
      n < 1 || nrows(responses) != n || ncols(responses) < 1 ||
      ncols(responses) > 64)
    error("invalid local-polynomial significance response tile");
  NPSigtestLPBatch *tile = (NPSigtestLPBatch *)R_alloc(1, sizeof(*tile));
  *tile = (NPSigtestLPBatch){.n = n, .draws = ncols(responses),
                             .response = REAL(responses)};
  const size_t cells = (size_t)XLENGTH(responses);
  for (size_t i = 0; i < cells; ++i)
    if (!R_FINITE(tile->response[i])) error("non-finite significance response tile");
  tile->weight = R_allocLD(n);
  tile->rescale = R_allocLD(n);
  tile->coefficient = (NPContrastNumber *)R_alloc(n, sizeof(*tile->coefficient));
  tile->residual = (double *)R_alloc(cells, sizeof(double));
  tile->scale = (double *)R_alloc(tile->draws, sizeof(double));
  tile->gradient = (double *)R_alloc(cells, sizeof(double));
  tile->stderr = (double *)R_alloc(cells, sizeof(double));
  for (int b = 0; b < tile->draws; ++b) tile->scale[b] = 0.0;
  return tile;
}

void np_sigtest_lp_residual(NPSigtestLPBatch *tile, int self,
                            const NPResidualOffDiagonal *geometry)
{
  if (tile->unsupported) return;
  const int n = tile->n;
  long double scale = 0.0L;
  tile->weight[self] = 0.0L;
  /* Replay the ordinary normalization's sequence of scale changes. This
   * preserves its summation order even for tiny weights or shifted Y. */
  for (int j = 0; j < n; ++j) {
    const long double w = tile->weight[j], a = fabsl(w);
    tile->rescale[j] = 1.0L;
    if (a > scale) {
      tile->rescale[j] = scale / a;
      scale = a;
    }
    if (a != 0.0L) tile->weight[j] = w / scale;
  }
  const long double norm = hypotl(geometry->sum, sqrtl(geometry->sumsq));
  for (int b = 0; b < tile->draws; ++b) {
    const double *y = tile->response + (size_t)b*n;
    long double numerator = 0.0L;
    for (int j = 0; j < n; ++j) {
      numerator *= tile->rescale[j];
      if (tile->weight[j] != 0.0L)
        numerator += tile->weight[j]*((long double)y[self] - y[j]);
    }
    const double value = (double)(numerator / norm);
    if (!isfinite(value)) { tile->unsupported = 1; return; }
    tile->residual[(size_t)b*n + self] = value;
    tile->scale[b] = fmax(tile->scale[b], fabs(value));
  }
  ++tile->residual_rows;
}

void np_sigtest_lp_standardize(NPSigtestLPBatch *tile)
{
  /* Precision and interpolation owners may not supply an identifiable
   * ordinary polynomial row. Their complete per-draw fit keeps ownership. */
  if (tile->unsupported || tile->residual_rows != tile->n) return;
  for (int b = 0; b < tile->draws; ++b) {
    double *r = tile->residual + (size_t)b*tile->n;
    if (tile->scale[b] > 0.0)
      for (int j = 0; j < tile->n; ++j)
        r[j] = (double)((long double)r[j] / tile->scale[b]);
  }
  tile->ready = 1;
}

void np_sigtest_lp_contrast(NPSigtestLPBatch *tile, int row, int anchor, double sign)
{
  if (!tile->ready || tile->unsupported) return;
  for (int b = 0; b < tile->draws; ++b) {
    const size_t offset = (size_t)b*tile->n;
    double gradient = 0.0, stderr = 0.0;
    if (anchor >= 0) {
      const double *y = tile->response + offset;
      const double *r = tile->residual + offset;
      NPRegressionContrastAccumulator sum = np_regression_contrast_begin(y[anchor], 1);
      for (int j = 0; j < tile->n; ++j)
        if (j != anchor)
          np_regression_contrast_add(&sum, tile->coefficient[j], y[j], 0,
                                     NP_RESIDUAL_IDENTIFIED, r[j]);
      const NPRegressionContrastResult applied = np_regression_contrast_finish(
        &sum, NP_RESIDUAL_IDENTIFIED, r[anchor], 0);
      gradient = sign * applied.contrast;
      stderr = (double)((long double)tile->scale[b] * applied.standard_error);
      if (applied.information != NP_RESIDUAL_IDENTIFIED ||
          !isfinite(gradient) || !isfinite(stderr)) {
        tile->unsupported = 1;
        return;
      }
    }
    tile->gradient[offset + row] = gradient;
    tile->stderr[offset + row] = stderr;
  }
  ++tile->contrast_rows;
}

void np_sigtest_lp_attach(NPSigtestLPBatch *tile, SEXP result)
{
  if (!tile->ready || tile->unsupported || tile->contrast_rows != tile->n) return;
  SEXP values = PROTECT(allocVector(VECSXP, 2));
  SEXP gradient = PROTECT(allocMatrix(REALSXP, tile->n, tile->draws));
  SEXP stderr = PROTECT(allocMatrix(REALSXP, tile->n, tile->draws));
  const size_t bytes = (size_t)tile->n * tile->draws * sizeof(double);
  memcpy(REAL(gradient), tile->gradient, bytes);
  memcpy(REAL(stderr), tile->stderr, bytes);
  SET_VECTOR_ELT(values, 0, gradient);
  SET_VECTOR_ELT(values, 1, stderr);
  setAttrib(result, install(".np.sigtest.batch"), values);
  UNPROTECT(3);
}
