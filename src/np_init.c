#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void gsl_bspline(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gsl_bspline_deriv(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density_conditional(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density_conditional_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_distribution_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_distribution_conditional_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_kernelsum(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_quantile_conditional(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_regression(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_regression_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_release_static_buffers(void *);
extern void np_set_seed(void *);
extern void np_set_tgauss2(void *);

static const R_CMethodDef CEntries[] = {
    {"gsl_bspline",                    (DL_FUNC) &gsl_bspline,                     9},
    {"gsl_bspline_deriv",              (DL_FUNC) &gsl_bspline_deriv,              11},
    {"np_density",                     (DL_FUNC) &np_density,                     18},
    {"np_density_bw",                  (DL_FUNC) &np_density_bw,                  16},
    {"np_density_conditional",         (DL_FUNC) &np_density_conditional,         30},
    {"np_density_conditional_bw",      (DL_FUNC) &np_density_conditional_bw,      21},
    {"np_distribution_bw",             (DL_FUNC) &np_distribution_bw,             19},
    {"np_distribution_conditional_bw", (DL_FUNC) &np_distribution_conditional_bw, 24},
    {"np_kernelsum",                   (DL_FUNC) &np_kernelsum,                   19},
    {"np_quantile_conditional",        (DL_FUNC) &np_quantile_conditional,        19},
    {"np_regression",                  (DL_FUNC) &np_regression,                  22},
    {"np_regression_bw",               (DL_FUNC) &np_regression_bw,               17},
    {"np_release_static_buffers",      (DL_FUNC) &np_release_static_buffers,       1},
    {"np_set_seed",                    (DL_FUNC) &np_set_seed,                     1},
    {"np_set_tgauss2",                 (DL_FUNC) &np_set_tgauss2,                  1},
    {NULL, NULL, 0}
};

void R_init_np(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
