/* Fixed-bound Gaussian corrected plot operator only. No estimator, RNG or MPI. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <limits.h>

SEXP C_np_plot_fgd_evaluate(SEXP target, SEXP donor, SEXP bandwidth,
                           SEXP coefficients, SEXP logmass)
{
    if (TYPEOF(target)!=REALSXP || TYPEOF(donor)!=REALSXP ||
        TYPEOF(bandwidth)!=REALSXP || XLENGTH(bandwidth)!=1 ||
        TYPEOF(coefficients)!=REALSXP || TYPEOF(logmass)!=REALSXP)
        error("invalid Gaussian plot operator types");
    R_xlen_t m=XLENGTH(target), n=XLENGTH(donor);
    SEXP dim=getAttrib(coefficients,R_DimSymbol);
    if (m<1 || n<1 || m>INT_MAX || n>INT_MAX || m>R_XLEN_T_MAX/n ||
        TYPEOF(dim)!=INTSXP || XLENGTH(dim)!=2 || INTEGER(dim)[0]!=m ||
        INTEGER(dim)[1]!=4 || XLENGTH(logmass)!=m ||
        (double)m*(double)n*8>32.0*1024*1024)
        error("Gaussian plot operator dimensions/workspace");
    double h=REAL(bandwidth)[0], scale=fmin(h,1.0);
    if (!R_FINITE(h) || h<=0) error("invalid Gaussian plot bandwidth");
    for (R_xlen_t i=0;i<m;++i) {
        if (!R_FINITE(REAL(target)[i]) || REAL(target)[i]<0 || REAL(target)[i]>1 ||
            !R_FINITE(REAL(logmass)[i])) error("invalid Gaussian plot target/mass");
        for (int k=0;k<4;++k)
            if (!R_FINITE(REAL(coefficients)[i+m*k])) error("invalid Gaussian plot coefficient");
    }
    for (R_xlen_t j=0;j<n;++j)
        if (!R_FINITE(REAL(donor)[j]) || REAL(donor)[j]<0 || REAL(donor)[j]>1)
            error("invalid Gaussian plot donor");
    SEXP ans=PROTECT(allocMatrix(REALSXP,(int)m,(int)n));
    for (R_xlen_t j=0;j<n;++j) {
        if (!(j%64)) R_CheckUserInterrupt();
        for (R_xlen_t i=0;i<m;++i) {
            double d=REAL(donor)[j]-REAL(target)[i], z=d/scale;
            double p=REAL(coefficients)[i]+z*(REAL(coefficients)[i+m]+
                z*(REAL(coefficients)[i+2*m]+z*REAL(coefficients)[i+3*m]));
            double value=expm1(-(d/h)*(d/h)*0.5-REAL(logmass)[i])*p+(p-1.0);
            if (!R_FINITE(value)) error("nonfinite Gaussian plot operator");
            REAL(ans)[i+m*j]=value;
        }
    }
    UNPROTECT(1);return ans;
}
