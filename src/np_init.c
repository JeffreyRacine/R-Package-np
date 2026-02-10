#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void gsl_bspline(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gsl_bspline_deriv(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density_conditional(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density_conditional_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_distribution_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_distribution_conditional_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_kernelsum(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_quantile_conditional(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_regression(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_regression_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_set_seed(void *);
extern void np_set_tgauss2(void *);

static const R_CMethodDef CEntries[] = {
    {"gsl_bspline",                    (DL_FUNC) &gsl_bspline,                     9},
    {"gsl_bspline_deriv",              (DL_FUNC) &gsl_bspline_deriv,              11},
    {"np_density",                     (DL_FUNC) &np_density,                     16},
    {"np_density_bw",                  (DL_FUNC) &np_density_bw,                  14},
    {"np_density_conditional",         (DL_FUNC) &np_density_conditional,         26},
    {"np_density_conditional_bw",      (DL_FUNC) &np_density_conditional_bw,      17},
    {"np_distribution_bw",             (DL_FUNC) &np_distribution_bw,             17},
    {"np_distribution_conditional_bw", (DL_FUNC) &np_distribution_conditional_bw, 20},
    {"np_kernelsum",                   (DL_FUNC) &np_kernelsum,                   17},
    {"np_quantile_conditional",        (DL_FUNC) &np_quantile_conditional,        19},
    {"np_regression",                  (DL_FUNC) &np_regression,                  20},
    {"np_regression_bw",               (DL_FUNC) &np_regression_bw,               15},
    {"np_set_seed",                    (DL_FUNC) &np_set_seed,                     1},
    {"np_set_tgauss2",                 (DL_FUNC) &np_set_tgauss2,                  1},
    {NULL, NULL, 0}
};

void R_init_npRmpi(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    /* npRmpi relies on dynamic symbol lookup for Rmpi bindings */
    R_useDynamicSymbols(dll, TRUE);
}
