/* Plot-only complete D8 positive-pilot covariance. No estimator/RNG state.
 * Inputs are borrowed double matrices/vectors; result is R-owned and protected.
 * One column uses O(m) output and constant arithmetic workspace. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>

static double cross_polynomial(int q, double A, double B, double D) {
  if (q == 1) return A*(A*(A*(B*(B*((1.0/2304.0)*B - 1.0/144.0) + 17.0/576.0) - 43.0/1296.0) + B*(B*(-1.0/144.0*B + (1.0/192.0)*D + 137.0/1152.0) - 1.0/18.0*D - 29.0/54.0) + (17.0/144.0)*D + 1631.0/2592.0) + B*(B*((17.0/576.0)*B - 1.0/18.0*D - 29.0/54.0) + D*((1.0/72.0)*D + 5.0/8.0) + 6517.0/2592.0) + D*(-2.0/27.0*D - 112.0/81.0) - 35.0/12.0) + B*(B*(-43.0/1296.0*B + (17.0/144.0)*D + 1631.0/2592.0) + D*(-2.0/27.0*D - 112.0/81.0) - 35.0/12.0) + D*(D*((1.0/162.0)*D + 133.0/324.0) + 665.0/216.0) + 11725.0/3888.0;
  if (q == 2) return A*(A*(A*(B*(B*((1.0/2304.0)*B - 1.0/128.0) + 11.0/288.0) - 65.0/1296.0) + B*(B*(-1.0/128.0*B + (1.0/192.0)*D + 43.0/288.0) - 1.0/16.0*D - 37.0/48.0) + (11.0/72.0)*D + 19.0/18.0) + B*(B*((11.0/288.0)*B - 1.0/16.0*D - 37.0/48.0) + D*((1.0/72.0)*D + 85.0/108.0) + 223.0/54.0) + D*(-1.0/12.0*D - 2) - 461.0/81.0) + B*(B*(-65.0/1296.0*B + (11.0/72.0)*D + 19.0/18.0) + D*(-1.0/12.0*D - 2) - 461.0/81.0) + D*(D*((1.0/162.0)*D + 14.0/27.0) + 418.0/81.0) + 5396.0/729.0;
  if (q == 3) return A*(A*(A*(B*(B*((1.0/2304.0)*B - 5.0/576.0) + 83.0/1728.0) - 47.0/648.0) + B*(B*(-5.0/576.0*B + (1.0/192.0)*D + 211.0/1152.0) - 5.0/72.0*D - 115.0/108.0) + (83.0/432.0)*D + 1447.0/864.0) + B*(B*((83.0/1728.0)*B - 5.0/72.0*D - 115.0/108.0) + D*((1.0/72.0)*D + 209.0/216.0) + 1849.0/288.0) + D*(-5.0/54.0*D - 25.0/9.0) - 6643.0/648.0) + B*(B*(-47.0/648.0*B + (83.0/432.0)*D + 1447.0/864.0) + D*(-5.0/54.0*D - 25.0/9.0) - 6643.0/648.0) + D*(D*((1.0/162.0)*D + 23.0/36.0) + 1757.0/216.0) + 61915.0/3888.0;
  return R_NaN;
}

static void finite_values(SEXP x) {
  for (R_xlen_t i=0; i<XLENGTH(x); ++i) {
    if (!R_FINITE(REAL(x)[i])) Rf_error("nonfinite D8 covariance input");
    if (i % 65536 == 0) R_CheckUserInterrupt();
  }
}

SEXP C_np_plot_udens_d8_column(SEXP xs, SEXP gs, SEXP hs, SEXP mus, SEXP js) {
  SEXP xd=Rf_getAttrib(xs,R_DimSymbol), gd=Rf_getAttrib(gs,R_DimSymbol);
  if (TYPEOF(xs)!=REALSXP || TYPEOF(gs)!=REALSXP ||
      TYPEOF(xd)!=INTSXP || XLENGTH(xd)!=2 ||
      TYPEOF(gd)!=INTSXP || XLENGTH(gd)!=2)
    Rf_error("D8 covariance requires double matrices");
  int n=INTEGER(xd)[0], q=INTEGER(xd)[1], m=INTEGER(gd)[0];
  if (n<2 || m<1 || q<1 || q>3 || INTEGER(gd)[1]!=q ||
      XLENGTH(xs)!=(R_xlen_t)n*q || XLENGTH(gs)!=(R_xlen_t)m*q ||
      TYPEOF(hs)!=REALSXP || XLENGTH(hs)!=q ||
      TYPEOF(mus)!=REALSXP || XLENGTH(mus)!=m ||
      TYPEOF(js)!=INTSXP || XLENGTH(js)!=1 || INTEGER(js)[0]==NA_INTEGER ||
      INTEGER(js)[0]<1 || INTEGER(js)[0]>m)
    Rf_error("invalid D8 covariance dimensions or column");
  finite_values(xs); finite_values(gs); finite_values(hs); finite_values(mus);
  double log_normal=-q*log(2*acos(-1.0))-q*log(3.0)/2;
  for (int j=0;j<q;++j) {
    if (REAL(hs)[j]<=0) Rf_error("D8 covariance bandwidth must be positive");
    log_normal-=2*log(REAL(hs)[j]);
  }
  const double *x=REAL(xs), *g=REAL(gs), *h=REAL(hs), *mu=REAL(mus);
  int b=INTEGER(js)[0]-1;
  SEXP ans=PROTECT(Rf_allocVector(REALSXP,m));
  for (int a=0;a<m;++a) {
    long double total=0;
    for (int i=0;i<n;++i) {
      double A=0,B=0,D=0,exponent=0;
      for (int j=0;j<q;++j) {
        double left=(g[a+(R_xlen_t)m*j]-x[i+(R_xlen_t)n*j])/h[j];
        double right=(g[b+(R_xlen_t)m*j]-x[i+(R_xlen_t)n*j])/h[j];
        double u=(2*left-right)/3, v=(2*right-left)/3;
        A+=u*u; B+=v*v; D+=u*v;
        exponent+=((left-right)*(left-right)+left*left+right*right)/6;
      }
      double product=cross_polynomial(q,A,B,D)*exp(log_normal-exponent);
      if (!R_FINITE(product)) Rf_error("nonfinite complete D8 covariance product");
      total+=product;
      if (i % 65536 == 0) R_CheckUserInterrupt();
    }
    REAL(ans)[a]=((double)(total/n)-mu[a]*mu[b])/n;
    if (!R_FINITE(REAL(ans)[a])) Rf_error("nonfinite complete D8 covariance column");
    if (a % 32 == 0) R_CheckUserInterrupt();
  }
  UNPROTECT(1); return ans;
}
