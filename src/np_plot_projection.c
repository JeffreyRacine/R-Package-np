/* Replay an accepted fit projection and kernel row. Precision-tagged rows
 * refactor the accepted weighted design; rank and ridge are never reselected. */
#include <R.h>
#include <Rinternals.h>
#include "headers.h"
#include "regression_precision.h"
#include "regression_contrast.h"

typedef struct { SEXP basis,projection,kernel,policy; NPRegressionPrecision q; } NPPrecisionReplay;
static void np_precision_replay_clear(void *data,Rboolean jump)
{ np_regression_precision_clear(&((NPPrecisionReplay *)data)->q); }
static SEXP np_precision_replay(void *data)
{
  NPPrecisionReplay *c=data;
  const int n=nrows(c->basis),p=ncols(c->basis),g=nrows(c->projection);
  if(TYPEOF(c->policy)!=REALSXP||!isMatrix(c->policy)||nrows(c->policy)!=g||ncols(c->policy)!=2||c->kernel==R_NilValue)
    error("invalid precision replay policy");
  SEXP out=PROTECT(allocMatrix(REALSXP,g,n));
  double **b=(double **)R_alloc(p,sizeof(double *));
  double *v=(double *)R_alloc(p,sizeof(double));
  for(int k=0;k<p;++k)b[k]=REAL(c->basis)+(size_t)n*k;
  for(int j=0;j<g;++j){
    R_CheckUserInterrupt();
    double ridge=REAL(c->policy)[j],anchor=REAL(c->policy)[g+j];
    if(!R_FINITE(ridge)||ridge<0||!R_FINITE(anchor)||anchor<0)error("invalid precision replay policy");
    for(int k=0;k<p;++k)v[k]=REAL(c->projection)[j+(size_t)g*k];
    const double *w=REAL(c->kernel)+(size_t)n*j;
    if(anchor>0){
      if(!np_regression_precision_factor(&c->q,b,w,n,p,ridge,anchor)||!np_regression_precision_row(&c->q,v))
        error("accepted precision row replay failed");
      for(int i=0;i<n;++i)REAL(out)[j+(size_t)g*i]=c->q.row[i];
    } else {
      for(int i=0;i<n;++i){
        NPContrastNumber sum={0};
        for(int k=0;k<p;++k)sum=np_contrast_add(sum,np_contrast_two_product(b[k][i],v[k]));
        double value=np_contrast_value(np_contrast_scale(sum,w[i]));
        if(!R_FINITE(value))error("nonfinite precision replay neighbor");
        REAL(out)[j+(size_t)g*i]=value;
      }
    }
  }
  UNPROTECT(1);return out;
}

SEXP C_np_plot_projection_apply(SEXP basis, SEXP projection, SEXP kernel)
{
  if(TYPEOF(basis) != REALSXP || !isMatrix(basis) ||
     TYPEOF(projection) != REALSXP || !isMatrix(projection))
    error("plot projection requires numeric matrices");
  const int n = nrows(basis), p = ncols(basis), g = nrows(projection);
  if(n < 1 || p < 1 || g < 1 || ncols(projection) != p ||
     (kernel != R_NilValue && (TYPEOF(kernel) != REALSXP || !isMatrix(kernel) ||
                              nrows(kernel) != n || ncols(kernel) != g)))
    error("plot projection dimensions do not agree");
  SEXP policy=getAttrib(projection,install(".np.fit.precision"));
  if(policy!=R_NilValue) {
    NPPrecisionReplay c={.basis=basis,.projection=projection,.kernel=kernel,.policy=policy};
    np_regression_precision_init(&c.q);
    return R_UnwindProtect(np_precision_replay,&c,np_precision_replay_clear,&c,NULL);
  }
  SEXP out = PROTECT(allocMatrix(REALSXP, g, n));
  SEXP accepted=getAttrib(projection,install(".np.fit.accepted"));
  if(accepted != R_NilValue) {
    if(TYPEOF(accepted)!=LGLSXP || XLENGTH(accepted)!=1 ||
       LOGICAL(accepted)[0]!=TRUE || kernel==R_NilValue)
      error("invalid accepted regression projection request");
    np_regression_inference_projection_apply(n,p,g,REAL(basis),
      REAL(projection),REAL(kernel),REAL(out));
    UNPROTECT(1);
    return out;
  }
  for(int j = 0; j < g; ++j) {
    R_CheckUserInterrupt();
    for(int i = 0; i < n; ++i) {
      double value = 0.0;
      for(int k = 0; k < p; ++k)
        value += REAL(basis)[i + (size_t)k*n] *
                 REAL(projection)[j + (size_t)k*g];
      if(kernel != R_NilValue) value *= REAL(kernel)[i + (size_t)j*n];
      if(!R_FINITE(value)) error("plot projection has nonfinite influence");
      REAL(out)[j + (size_t)i*g] = value;
    }
  }
  UNPROTECT(1);
  return out;
}
