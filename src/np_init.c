#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void gsl_bspline(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gsl_bspline_deriv(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density_conditional(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_density_conditional_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_dim_basis(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_distribution_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_distribution_conditional_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_kernelsum(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_quantile_conditional(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_regression(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_regression_bw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void np_release_static_buffers(void *);
extern void np_set_seed(void *);
extern void np_set_tgauss2(void *);

/* .Call calls */
extern SEXP C_gsl_bspline(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsl_bspline_deriv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_dim_basis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_set_seed(SEXP);
extern SEXP C_np_set_tgauss2(SEXP);
extern SEXP C_np_release_static_buffers(void);

static const R_CMethodDef CEntries[] = {
    {"np_density",                     (DL_FUNC) &np_density,                     18},
    {"np_density_bw",                  (DL_FUNC) &np_density_bw,                  16},
    {"np_density_conditional",         (DL_FUNC) &np_density_conditional,         30},
    {"np_density_conditional_bw",      (DL_FUNC) &np_density_conditional_bw,      23},
    {"np_distribution_bw",             (DL_FUNC) &np_distribution_bw,             19},
    {"np_distribution_conditional_bw", (DL_FUNC) &np_distribution_conditional_bw, 26},
    {"np_kernelsum",                   (DL_FUNC) &np_kernelsum,                   19},
    {"np_quantile_conditional",        (DL_FUNC) &np_quantile_conditional,        19},
    {"np_regression",                  (DL_FUNC) &np_regression,                  25},
    {"np_regression_bw",               (DL_FUNC) &np_regression_bw,               22},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"C_gsl_bspline",                  (DL_FUNC) &C_gsl_bspline,                   7},
    {"C_gsl_bspline_deriv",            (DL_FUNC) &C_gsl_bspline_deriv,             8},
    {"C_np_dim_basis",                 (DL_FUNC) &C_np_dim_basis,                  6},
    {"C_np_set_seed",                  (DL_FUNC) &C_np_set_seed,                   1},
    {"C_np_set_tgauss2",               (DL_FUNC) &C_np_set_tgauss2,                1},
    {"C_np_release_static_buffers",    (DL_FUNC) &C_np_release_static_buffers,     0},
    {NULL, NULL, 0}
};

void R_init_np(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
