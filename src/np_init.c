#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* Routine registration for the serial np shared library. */

/* .Call calls */
extern SEXP C_np_dim_basis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_nomad_native_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_nomad_shadow_prepare(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_nomad_shadow_eval(SEXP, SEXP);
extern SEXP C_np_density_conditional_nomad_shadow_native_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_nomad_shadow_fixed_native_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_nomad_shadow_clear(void);
extern SEXP C_np_density_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_conditional_count_levels(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_lp_apply_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_lc_hat_normalize(SEXP, SEXP);
extern SEXP C_np_reghat_lp_matrix_fast(SEXP, SEXP, SEXP);
extern SEXP C_np_npscoef_batch_zero_solve(SEXP, SEXP);
extern SEXP C_np_npscoef_batch_project(SEXP, SEXP);
extern SEXP C_np_entropy_gaussian_integrand(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_categorical_profile_kernel_tile(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_nomad_native_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_conditional_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_conditional_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_conditional_nomad_native_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_kernelsum(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_kernelsum_power12(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_quantile_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_nomad_native_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_nomad_r_callback_native_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_lsqregression_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_lsqregression_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_progress_signal(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_progress_fit_begin(SEXP);
extern SEXP C_np_progress_fit_end(void);
extern SEXP C_np_reset_native_estimator_state(void);
extern SEXP C_np_set_seed(SEXP);
extern SEXP C_np_release_static_buffers(void);

static const R_CallMethodDef CallEntries[] = {
    {"C_np_dim_basis",                 (DL_FUNC) &C_np_dim_basis,                  6},
    {"C_np_density",                   (DL_FUNC) &C_np_density,                   16},
    {"C_np_density_bw",                (DL_FUNC) &C_np_density_bw,                12},
    {"C_np_density_bw_eval",           (DL_FUNC) &C_np_density_bw_eval,           12},
    {"C_np_density_nomad_native_search",(DL_FUNC) &C_np_density_nomad_native_search,19},
    {"C_np_density_conditional_nomad_shadow_prepare",(DL_FUNC) &C_np_density_conditional_nomad_shadow_prepare,20},
    {"C_np_density_conditional_nomad_shadow_eval",(DL_FUNC) &C_np_density_conditional_nomad_shadow_eval,2},
    {"C_np_density_conditional_nomad_shadow_native_search",(DL_FUNC) &C_np_density_conditional_nomad_shadow_native_search,12},
    {"C_np_density_conditional_nomad_shadow_fixed_native_search",(DL_FUNC) &C_np_density_conditional_nomad_shadow_fixed_native_search,12},
    {"C_np_density_conditional_nomad_shadow_clear",(DL_FUNC) &C_np_density_conditional_nomad_shadow_clear,0},
    {"C_np_density_conditional",       (DL_FUNC) &C_np_density_conditional,       31},
    {"C_np_conditional_count_levels",  (DL_FUNC) &C_np_conditional_count_levels,  22},
    {"C_np_regression_lp_apply_conditional",(DL_FUNC) &C_np_regression_lp_apply_conditional,23},
    {"C_np_lc_hat_normalize",             (DL_FUNC) &C_np_lc_hat_normalize,              2},
    {"C_np_reghat_lp_matrix_fast",       (DL_FUNC) &C_np_reghat_lp_matrix_fast,        3},
    {"C_np_npscoef_batch_zero_solve",    (DL_FUNC) &C_np_npscoef_batch_zero_solve,     2},
    {"C_np_npscoef_batch_project",       (DL_FUNC) &C_np_npscoef_batch_project,        2},
    {"C_np_entropy_gaussian_integrand",  (DL_FUNC) &C_np_entropy_gaussian_integrand,   4},
    {"C_np_categorical_profile_kernel_tile",(DL_FUNC) &C_np_categorical_profile_kernel_tile,12},
    {"C_np_density_conditional_bw",    (DL_FUNC) &C_np_density_conditional_bw,    21},
    {"C_np_density_conditional_bw_eval",(DL_FUNC) &C_np_density_conditional_bw_eval,21},
    {"C_np_distribution_bw",           (DL_FUNC) &C_np_distribution_bw,           15},
    {"C_np_distribution_bw_eval",      (DL_FUNC) &C_np_distribution_bw_eval,      15},
    {"C_np_distribution_nomad_native_search",(DL_FUNC) &C_np_distribution_nomad_native_search,22},
    {"C_np_distribution_conditional_bw",(DL_FUNC) &C_np_distribution_conditional_bw,24},
    {"C_np_distribution_conditional_bw_eval",(DL_FUNC) &C_np_distribution_conditional_bw_eval,24},
    {"C_np_distribution_conditional_nomad_native_search",(DL_FUNC) &C_np_distribution_conditional_nomad_native_search,31},
    {"C_np_kernelsum",                 (DL_FUNC) &C_np_kernelsum,                 19},
    {"C_np_kernelsum_power12",         (DL_FUNC) &C_np_kernelsum_power12,         19},
    {"C_np_progress_fit_begin",        (DL_FUNC) &C_np_progress_fit_begin,         1},
    {"C_np_progress_fit_end",          (DL_FUNC) &C_np_progress_fit_end,           0},
    {"C_np_progress_signal",           (DL_FUNC) &C_np_progress_signal,            4},
    {"C_np_reset_native_estimator_state",(DL_FUNC) &C_np_reset_native_estimator_state,0},
    {"C_np_quantile_conditional",      (DL_FUNC) &C_np_quantile_conditional,      19},
    {"C_np_regression_bw",             (DL_FUNC) &C_np_regression_bw,             16},
    {"C_np_regression_bw_eval",        (DL_FUNC) &C_np_regression_bw_eval,        16},
    {"C_np_regression_nomad_native_search",(DL_FUNC) &C_np_regression_nomad_native_search,24},
    {"C_np_nomad_r_callback_native_search",(DL_FUNC) &C_np_nomad_r_callback_native_search,11},
    {"C_np_lsqregression_bw",          (DL_FUNC) &C_np_lsqregression_bw,          20},
    {"C_np_lsqregression_bw_eval",     (DL_FUNC) &C_np_lsqregression_bw_eval,     20},
    {"C_np_regression",                (DL_FUNC) &C_np_regression,                24},
    {"C_np_set_seed",                  (DL_FUNC) &C_np_set_seed,                   1},
    {"C_np_release_static_buffers",    (DL_FUNC) &C_np_release_static_buffers,     0},
    {NULL, NULL, 0}
};

void R_init_np(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    /* Serial package uses registered symbols only. */
    R_useDynamicSymbols(dll, FALSE);
}
