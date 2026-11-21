/* Polynomial reproduction in a complete, accepted univariate Gaussian pilot.
   Both operations retain the original zero-ridge admission. */
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include "regression_contrast.h"

typedef NPContrastNumber Wide;
#define number np_contrast_number
#define add np_contrast_add
#define sub np_contrast_subtract
#define mul np_contrast_multiply
#define divw np_contrast_divide
#define scale np_contrast_scale
#define value np_contrast_value

/* Reproduce the accepted weighted polynomial in orthogonal coordinates.
   The derivative of a Gaussian weight is its linear score times the weight.
   Multiplying the reproducing polynomial of degree q by that score leaves
   one degree q+1 term after projecting out the degree q space. */
SEXP C_np_bias_gaussian_extrapolate(SEXP xs, SEXP target, SEXP bandwidth,
                                   SEXP degree, SEXP weights)
{
    if (!isReal(xs) || !isReal(weights) || XLENGTH(xs)!=XLENGTH(weights) ||
        !isReal(target) || XLENGTH(target)!=1 || !isReal(bandwidth) ||
        XLENGTH(bandwidth)!=1 || !isInteger(degree) || XLENGTH(degree)!=1)
        error("invalid Gaussian pilot inputs");
    const int n=LENGTH(xs), q=INTEGER(degree)[0];
    const double *x=REAL(xs), *weight=REAL(weights);
    const double e=REAL(target)[0], bw=REAL(bandwidth)[0];
    if (q<1 || q>=n || !R_FINITE(e) || !R_FINITE(bw) || bw<=0)
        error("invalid Gaussian pilot geometry");
    int anchor=0;
    double span=0;
    for (int i=0; i<n; ++i) {
        if (!R_FINITE(x[i]) || !R_FINITE(weight[i]) || weight[i]<=0)
            error("complete Gaussian pilot support required");
        if (weight[i]>weight[anchor]) anchor=i;
        span=fmax(span,fabs(x[i]-e));
    }
    Wide *t=(Wide *)R_alloc(n,sizeof(Wide)), *omega=(Wide *)R_alloc(n,sizeof(Wide));
    /* Keep the polynomial columns used for reorthogonalization. */
    Wide *basis=(Wide *)R_alloc((size_t)n*(q+1),sizeof(Wide));
    Wide *at=(Wide *)R_alloc(q+1,sizeof(Wide)), *slope=(Wide *)R_alloc(q+1,sizeof(Wide));
    Wide *norm=(Wide *)R_alloc(q+1,sizeof(Wide));
    Wide *next=(Wide *)R_alloc(n,sizeof(Wide)), *row=(Wide *)R_alloc(n,sizeof(Wide));
    Wide *gradient=(Wide *)R_alloc(n,sizeof(Wide));
    const Wide inv_span=divw(number(1),number(span)), inv_bw=divw(number(1),number(bw));
    Wide offset=np_contrast_two_sum(e,-x[anchor]);
    Wide where=mul(offset,inv_span), offset_bw=mul(offset,inv_bw);
    const double kernel_scale=1.0/bw, anchor_u=(x[anchor]-e)*kernel_scale;
    const double anchor_exponent=-0.5*anchor_u*anchor_u;
    Wide mass=number(0);
    for (int i=0; i<n; ++i) {
        Wide delta=np_contrast_two_sum(x[i],-x[anchor]);
        t[i]=mul(delta,inv_span);
        /* Retain the actual kernel weight ratio. Correct rounding of the
           standardized distance in the Gaussian log ratio about its mode. */
        Wide u=mul(delta,inv_bw);
        Wide exact_log=mul(u,sub(offset_bw,scale(u,0.5)));
        double z=(x[i]-e)*kernel_scale;
        Wide kernel_log=np_contrast_two_sum(-0.5*z*z,-anchor_exponent);
        Wide ratio=divw(number(weight[i]),number(weight[anchor]));
        omega[i]=add(ratio,scale(ratio,expm1(value(sub(exact_log,kernel_log)))));
        mass=add(mass,omega[i]);
        basis[i]=number(1); row[i]=gradient[i]=number(0);
    }
    Wide inv_mass=divw(number(1),mass);
    for (int i=0; i<n; ++i) omega[i]=mul(omega[i],inv_mass);
    at[0]=number(1); slope[0]=number(0);
    const Wide score_scale=mul(scale(inv_bw,span),inv_bw);
    for (int j=0; j<=q; ++j) {
        const Wide *poly=basis+(size_t)n*j;
        norm[j]=number(0);
        for (int i=0; i<n; ++i) {
            norm[j]=add(norm[j],mul(mul(omega[i],poly[i]),poly[i]));
            next[i]=mul(t[i],poly[i]);
        }
        Wide coefficient=divw(at[j],norm[j]);
        Wide change=divw(mul(slope[j],inv_span),norm[j]);
        Wide new_at=mul(where,at[j]), new_slope=add(at[j],mul(where,slope[j]));
        /* Orthogonalize t*p[j] twice. A three-term recurrence alone can
           lose orthogonality near saturation on a discrete training grid. */
        for (int pass=0; pass<2; ++pass) {
            for (int k=0; k<=j; ++k) {
                const Wide *previous=basis+(size_t)n*k;
                Wide projection=number(0);
                for (int i=0; i<n; ++i)
                    projection=add(projection,mul(mul(omega[i],previous[i]),next[i]));
                projection=divw(projection,norm[k]);
                for (int i=0; i<n; ++i) next[i]=sub(next[i],mul(projection,previous[i]));
                new_at=sub(new_at,mul(projection,at[k]));
                new_slope=sub(new_slope,mul(projection,slope[k]));
            }
        }
        Wide next_norm=number(0);
        for (int i=0; i<n; ++i) {
            Wide direction=mul(omega[i],poly[i]);
            row[i]=add(row[i],mul(coefficient,direction));
            gradient[i]=add(gradient[i],mul(change,direction));
            if (j==q) gradient[i]=add(gradient[i],mul(coefficient,mul(omega[i],mul(score_scale,next[i]))));
            else next_norm=add(next_norm,mul(mul(omega[i],next[i]),next[i]));
        }
        if (j<q) {
            double length=sqrt(value(next_norm));
            if (!(length>0) || !R_FINITE(length))
                error("nonfinite orthogonal Gaussian pilot on an accepted polynomial design");
            Wide inv_length=divw(number(1),number(length));
            at[j+1]=mul(new_at,inv_length); slope[j+1]=mul(new_slope,inv_length);
            for (int i=0; i<n; ++i) basis[(size_t)n*(j+1)+i]=mul(next[i],inv_length);
        }
    }
    SEXP result=PROTECT(allocMatrix(REALSXP,n,2));
    for(int i=0; i<n; ++i) {
        REAL(result)[i]=value(row[i]); REAL(result)[n+i]=value(gradient[i]);
    }
    UNPROTECT(1);
    return result;
}

/* Clenshaw recurrence for y minus its polynomial control in physical units. */
SEXP C_np_bias_gaussian_center(SEXP xs, SEXP ys, SEXP coefficients,
                              SEXP limits, SEXP magnitude)
{
    if (!isReal(xs) || !isReal(ys) || XLENGTH(xs)!=XLENGTH(ys) || !isReal(coefficients) ||
        !XLENGTH(coefficients) || !isReal(limits) || XLENGTH(limits)!=2 ||
        !isReal(magnitude) || XLENGTH(magnitude)!=1 || !(REAL(limits)[0]<REAL(limits)[1]))
        error("invalid Gaussian residual centering inputs");
    int n=LENGTH(xs), p=LENGTH(coefficients);
    const double *co=REAL(coefficients), lo=REAL(limits)[0];
    Wide inv_width=divw(number(2),np_contrast_two_sum(REAL(limits)[1],-lo));
    SEXP out=PROTECT(allocVector(REALSXP,n));
    for(int i=0;i<n;++i) {
        Wide t=sub(mul(np_contrast_two_sum(REAL(xs)[i],-lo),inv_width),number(1));
        Wide b1=number(0), b2=number(0);
        for(int j=p-1;j>0;--j) {
            Wide b=add(number(co[j]),sub(scale(mul(t,b1),2),b2));
            b2=b1; b1=b;
        }
        Wide fitted=add(number(co[0]),sub(mul(t,b1),b2));
        REAL(out)[i]=value(sub(number(REAL(ys)[i]),scale(fitted,REAL(magnitude)[0])));
    }
    UNPROTECT(1);
    return out;
}
