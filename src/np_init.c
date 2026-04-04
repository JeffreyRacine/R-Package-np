#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* Routine registration for the serial np shared library. */

/* .Call calls */
extern SEXP C_gsl_bspline(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsl_bspline_deriv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_dim_basis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_shadow_cv_density_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_shadow_cv_distribution_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_shadow_cv_xweights_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_shadow_cv_xweights_full_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_shadow_cv_yrow_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_lp_apply_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_conditional_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_conditional_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_kernelsum(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_quantile_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_progress_signal(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_progress_fit_begin(SEXP);
extern SEXP C_np_progress_fit_end(void);
extern SEXP C_np_shadow_reset_state(void);
extern SEXP C_np_set_seed(SEXP);
extern SEXP C_np_set_tgauss2(SEXP);
extern SEXP C_np_release_static_buffers(void);

static const R_CallMethodDef CallEntries[] = {
    {"C_gsl_bspline",                  (DL_FUNC) &C_gsl_bspline,                   7},
    {"C_gsl_bspline_deriv",            (DL_FUNC) &C_gsl_bspline_deriv,             8},
    {"C_np_dim_basis",                 (DL_FUNC) &C_np_dim_basis,                  6},
    {"C_np_density",                   (DL_FUNC) &C_np_density,                   16},
    {"C_np_density_bw",                (DL_FUNC) &C_np_density_bw,                12},
    {"C_np_density_conditional",       (DL_FUNC) &C_np_density_conditional,       31},
    {"C_np_shadow_cv_density_conditional",(DL_FUNC) &C_np_shadow_cv_density_conditional,21},
    {"C_np_shadow_cv_distribution_conditional",(DL_FUNC) &C_np_shadow_cv_distribution_conditional,24},
    {"C_np_shadow_cv_xweights_conditional",(DL_FUNC) &C_np_shadow_cv_xweights_conditional,17},
    {"C_np_shadow_cv_xweights_full_conditional",(DL_FUNC) &C_np_shadow_cv_xweights_full_conditional,17},
    {"C_np_shadow_cv_yrow_conditional",(DL_FUNC) &C_np_shadow_cv_yrow_conditional,16},
    {"C_np_regression_lp_apply_conditional",(DL_FUNC) &C_np_regression_lp_apply_conditional,17},
    {"C_np_density_conditional_bw",    (DL_FUNC) &C_np_density_conditional_bw,    21},
    {"C_np_density_conditional_bw_eval",(DL_FUNC) &C_np_density_conditional_bw_eval,21},
    {"C_np_distribution_bw",           (DL_FUNC) &C_np_distribution_bw,           15},
    {"C_np_distribution_conditional_bw",(DL_FUNC) &C_np_distribution_conditional_bw,24},
    {"C_np_distribution_conditional_bw_eval",(DL_FUNC) &C_np_distribution_conditional_bw_eval,24},
    {"C_np_kernelsum",                 (DL_FUNC) &C_np_kernelsum,                 19},
    {"C_np_progress_fit_begin",        (DL_FUNC) &C_np_progress_fit_begin,         1},
    {"C_np_progress_fit_end",          (DL_FUNC) &C_np_progress_fit_end,           0},
    {"C_np_progress_signal",           (DL_FUNC) &C_np_progress_signal,            4},
    {"C_np_shadow_reset_state",        (DL_FUNC) &C_np_shadow_reset_state,         0},
    {"C_np_quantile_conditional",      (DL_FUNC) &C_np_quantile_conditional,      19},
    {"C_np_regression_bw",             (DL_FUNC) &C_np_regression_bw,             16},
    {"C_np_regression_bw_eval",        (DL_FUNC) &C_np_regression_bw_eval,        16},
    {"C_np_regression",                (DL_FUNC) &C_np_regression,                24},
    {"C_np_set_seed",                  (DL_FUNC) &C_np_set_seed,                   1},
    {"C_np_set_tgauss2",               (DL_FUNC) &C_np_set_tgauss2,                1},
    {"C_np_release_static_buffers",    (DL_FUNC) &C_np_release_static_buffers,     0},
    {NULL, NULL, 0}
};

void R_init_np(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    /* Serial package uses registered symbols only. */
    R_useDynamicSymbols(dll, FALSE);
}
