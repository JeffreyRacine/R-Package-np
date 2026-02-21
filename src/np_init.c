#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP C_gsl_bspline(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsl_bspline_deriv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_dim_basis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_conditional_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_kernelsum(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_quantile_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_set_seed(SEXP);
extern SEXP C_np_set_tgauss2(SEXP);
extern SEXP C_np_release_static_buffers(void);

static const R_CallMethodDef CallEntries[] = {
    {"C_gsl_bspline",                  (DL_FUNC) &C_gsl_bspline,                   7},
    {"C_gsl_bspline_deriv",            (DL_FUNC) &C_gsl_bspline_deriv,             8},
    {"C_np_dim_basis",                 (DL_FUNC) &C_np_dim_basis,                  6},
    {"C_np_density",                   (DL_FUNC) &C_np_density,                   16},
    {"C_np_density_bw",                (DL_FUNC) &C_np_density_bw,                12},
    {"C_np_density_conditional",       (DL_FUNC) &C_np_density_conditional,       27},
    {"C_np_density_conditional_bw",    (DL_FUNC) &C_np_density_conditional_bw,    17},
    {"C_np_distribution_bw",           (DL_FUNC) &C_np_distribution_bw,           15},
    {"C_np_distribution_conditional_bw",(DL_FUNC) &C_np_distribution_conditional_bw,20},
    {"C_np_kernelsum",                 (DL_FUNC) &C_np_kernelsum,                 19},
    {"C_np_quantile_conditional",      (DL_FUNC) &C_np_quantile_conditional,      19},
    {"C_np_regression_bw",             (DL_FUNC) &C_np_regression_bw,             16},
    {"C_np_regression",                (DL_FUNC) &C_np_regression,                23},
    {"C_np_set_seed",                  (DL_FUNC) &C_np_set_seed,                   1},
    {"C_np_set_tgauss2",               (DL_FUNC) &C_np_set_tgauss2,                1},
    {"C_np_release_static_buffers",    (DL_FUNC) &C_np_release_static_buffers,     0},
    {NULL, NULL, 0}
};

void R_init_np(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
