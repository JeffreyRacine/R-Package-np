#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* Routine registration for the npRmpi shared library. */

/* .Call calls */
extern SEXP C_gsl_bspline(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsl_bspline_deriv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_dim_basis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_nomad_shadow_prepare(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_density_conditional_nomad_shadow_eval(SEXP, SEXP);
extern SEXP C_np_density_conditional_nomad_shadow_clear(void);
extern SEXP C_np_shadow_cv_xweights_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_conditional_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_distribution_conditional_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_kernelsum(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_progress_fit_begin(SEXP);
extern SEXP C_np_progress_fit_end(void);
extern SEXP C_np_quantile_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_bw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_bw_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_nomad_shadow_prepare(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_nomad_shadow_eval(SEXP, SEXP);
extern SEXP C_np_regression_nomad_shadow_clear(void);
extern SEXP C_np_regression(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_regression_lp_apply_conditional(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_progress_signal(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_np_shadow_reset_state(void);
extern SEXP C_np_set_seed(SEXP);
extern SEXP C_np_set_tgauss2(SEXP);
extern SEXP C_np_set_local_regression_mode(SEXP);
extern SEXP C_np_release_static_buffers(void);
extern SEXP C_np_mpi_init(void);
extern SEXP mpidist(void);
extern SEXP mkstr(SEXP);
extern SEXP mpi_initialize(void);
extern SEXP mpi_finalize(void);
extern SEXP mpi_get_version(void);
extern SEXP mpi_get_processor_name(void);
extern SEXP mpi_universe_size(void);
extern SEXP mpi_any_source(void);
extern SEXP mpi_any_tag(void);
extern SEXP mpi_undefined(void);
extern SEXP mpi_proc_null(void);
extern SEXP mpi_info_create(SEXP);
extern SEXP mpi_info_set(SEXP, SEXP, SEXP);
extern SEXP mpi_info_get(SEXP, SEXP, SEXP);
extern SEXP mpi_info_free(SEXP);
extern SEXP mpi_realloc_comm(SEXP);
extern SEXP mpi_comm_maxsize(void);
extern SEXP mpi_realloc_status(SEXP);
extern SEXP mpi_status_maxsize(void);
extern SEXP mpi_realloc_request(SEXP);
extern SEXP mpi_request_maxsize(void);
extern SEXP mpi_realloc_datatype(SEXP);
extern SEXP mpi_gather(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_gatherv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_scatter(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_scatterv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_allgather(SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_allgatherv(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_bcast(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_send(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_recv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_reduce(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_allreduce(SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_iprobe(SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_probe(SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_get_count(SEXP, SEXP);
extern SEXP mpi_get_sourcetag(SEXP);
extern SEXP mpi_barrier(SEXP);
extern SEXP mpi_comm_is_null(SEXP);
extern SEXP mpi_comm_size(SEXP);
extern SEXP mpi_comm_rank(SEXP);
extern SEXP mpi_comm_dup(SEXP, SEXP);
extern SEXP mpi_comm_c2f(SEXP);
extern SEXP mpi_comm_free(SEXP);
extern SEXP mpi_abort(SEXP);
extern SEXP mpi_comm_set_errhandler(SEXP);
extern SEXP mpi_comm_test_inter(SEXP);
#ifdef MPI2
extern SEXP mpi_comm_spawn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
#endif
extern SEXP mpi_comm_get_parent(SEXP);
extern SEXP mpi_is_master(void);
extern SEXP mpi_comm_disconnect(SEXP);
extern SEXP mpi_intercomm_merge(SEXP, SEXP, SEXP);
extern SEXP mpi_comm_remote_size(SEXP);
extern SEXP mpi_sendrecv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_sendrecv_replace(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_cart_create(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_dims_create(SEXP, SEXP, SEXP);
extern SEXP mpi_cartdim_get(SEXP);
extern SEXP mpi_cart_get(SEXP, SEXP);
extern SEXP mpi_cart_rank(SEXP, SEXP);
extern SEXP mpi_cart_coords(SEXP, SEXP, SEXP);
extern SEXP mpi_cart_shift(SEXP, SEXP, SEXP);
extern SEXP mpi_isend(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_irecv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpi_wait(SEXP, SEXP);
extern SEXP mpi_test(SEXP, SEXP);
extern SEXP mpi_cancel(SEXP);
extern SEXP mpi_test_cancelled(SEXP);
extern SEXP mpi_waitany(SEXP, SEXP);
extern SEXP mpi_testany(SEXP, SEXP);
extern SEXP mpi_waitall(SEXP);
extern SEXP mpi_testall(SEXP);
extern SEXP mpi_testsome(SEXP);
extern SEXP mpi_waitsome(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"C_gsl_bspline",                  (DL_FUNC) &C_gsl_bspline,                   7},
    {"C_gsl_bspline_deriv",            (DL_FUNC) &C_gsl_bspline_deriv,             8},
    {"C_np_dim_basis",                 (DL_FUNC) &C_np_dim_basis,                  6},
    {"C_np_density",                   (DL_FUNC) &C_np_density,                   16},
    {"C_np_density_bw",                (DL_FUNC) &C_np_density_bw,                12},
    {"C_np_density_conditional",       (DL_FUNC) &C_np_density_conditional,       31},
    {"C_np_density_conditional_bw",    (DL_FUNC) &C_np_density_conditional_bw,    21},
    {"C_np_density_conditional_bw_eval",(DL_FUNC) &C_np_density_conditional_bw_eval,21},
    {"C_np_density_conditional_nomad_shadow_prepare",(DL_FUNC) &C_np_density_conditional_nomad_shadow_prepare,20},
    {"C_np_density_conditional_nomad_shadow_eval",(DL_FUNC) &C_np_density_conditional_nomad_shadow_eval,2},
    {"C_np_density_conditional_nomad_shadow_clear",(DL_FUNC) &C_np_density_conditional_nomad_shadow_clear,0},
    {"C_np_shadow_cv_xweights_conditional",(DL_FUNC) &C_np_shadow_cv_xweights_conditional,17},
    {"C_np_distribution_bw",           (DL_FUNC) &C_np_distribution_bw,           15},
    {"C_np_distribution_conditional_bw",(DL_FUNC) &C_np_distribution_conditional_bw,24},
    {"C_np_distribution_conditional_bw_eval",(DL_FUNC) &C_np_distribution_conditional_bw_eval,24},
    {"C_np_kernelsum",                 (DL_FUNC) &C_np_kernelsum,                 19},
    {"C_np_progress_fit_begin",        (DL_FUNC) &C_np_progress_fit_begin,         1},
    {"C_np_progress_fit_end",          (DL_FUNC) &C_np_progress_fit_end,           0},
    {"C_np_progress_signal",           (DL_FUNC) &C_np_progress_signal,            4},
    {"C_np_quantile_conditional",      (DL_FUNC) &C_np_quantile_conditional,      19},
    {"C_np_regression_bw",             (DL_FUNC) &C_np_regression_bw,             16},
    {"C_np_regression_bw_eval",        (DL_FUNC) &C_np_regression_bw_eval,        16},
    {"C_np_regression_nomad_shadow_prepare",(DL_FUNC) &C_np_regression_nomad_shadow_prepare,15},
    {"C_np_regression_nomad_shadow_eval",(DL_FUNC) &C_np_regression_nomad_shadow_eval,2},
    {"C_np_regression_nomad_shadow_clear",(DL_FUNC) &C_np_regression_nomad_shadow_clear,0},
    {"C_np_regression",                (DL_FUNC) &C_np_regression,                24},
    {"C_np_regression_lp_apply_conditional",(DL_FUNC) &C_np_regression_lp_apply_conditional,17},
    {"C_np_set_seed",                  (DL_FUNC) &C_np_set_seed,                   1},
    {"C_np_shadow_reset_state",        (DL_FUNC) &C_np_shadow_reset_state,         0},
    {"C_np_set_tgauss2",               (DL_FUNC) &C_np_set_tgauss2,                1},
    {"C_np_set_local_regression_mode", (DL_FUNC) &C_np_set_local_regression_mode,  1},
    {"C_np_release_static_buffers",    (DL_FUNC) &C_np_release_static_buffers,     0},
    {"C_np_mpi_init",                  (DL_FUNC) &C_np_mpi_init,                   0},
    {"mpidist",                       (DL_FUNC) &mpidist,                         0},
    {"mkstr",                         (DL_FUNC) &mkstr,                           1},
    {"mpi_initialize",                (DL_FUNC) &mpi_initialize,                  0},
    {"mpi_finalize",                  (DL_FUNC) &mpi_finalize,                    0},
    {"mpi_get_version",               (DL_FUNC) &mpi_get_version,                 0},
    {"mpi_get_processor_name",        (DL_FUNC) &mpi_get_processor_name,          0},
    {"mpi_universe_size",             (DL_FUNC) &mpi_universe_size,               0},
    {"mpi_any_source",                (DL_FUNC) &mpi_any_source,                  0},
    {"mpi_any_tag",                   (DL_FUNC) &mpi_any_tag,                     0},
    {"mpi_undefined",                 (DL_FUNC) &mpi_undefined,                   0},
    {"mpi_proc_null",                 (DL_FUNC) &mpi_proc_null,                   0},
    {"mpi_info_create",               (DL_FUNC) &mpi_info_create,                 1},
    {"mpi_info_set",                  (DL_FUNC) &mpi_info_set,                    3},
    {"mpi_info_get",                  (DL_FUNC) &mpi_info_get,                    3},
    {"mpi_info_free",                 (DL_FUNC) &mpi_info_free,                   1},
    {"mpi_realloc_comm",              (DL_FUNC) &mpi_realloc_comm,                1},
    {"mpi_comm_maxsize",              (DL_FUNC) &mpi_comm_maxsize,                0},
    {"mpi_realloc_status",            (DL_FUNC) &mpi_realloc_status,              1},
    {"mpi_status_maxsize",            (DL_FUNC) &mpi_status_maxsize,              0},
    {"mpi_realloc_request",           (DL_FUNC) &mpi_realloc_request,             1},
    {"mpi_request_maxsize",           (DL_FUNC) &mpi_request_maxsize,             0},
    {"mpi_realloc_datatype",          (DL_FUNC) &mpi_realloc_datatype,            1},
    {"mpi_gather",                    (DL_FUNC) &mpi_gather,                      5},
    {"mpi_gatherv",                   (DL_FUNC) &mpi_gatherv,                     6},
    {"mpi_scatter",                   (DL_FUNC) &mpi_scatter,                     5},
    {"mpi_scatterv",                  (DL_FUNC) &mpi_scatterv,                    6},
    {"mpi_allgather",                 (DL_FUNC) &mpi_allgather,                   4},
    {"mpi_allgatherv",                (DL_FUNC) &mpi_allgatherv,                  5},
    {"mpi_bcast",                     (DL_FUNC) &mpi_bcast,                       5},
    {"mpi_send",                      (DL_FUNC) &mpi_send,                        5},
    {"mpi_recv",                      (DL_FUNC) &mpi_recv,                        6},
    {"mpi_reduce",                    (DL_FUNC) &mpi_reduce,                      5},
    {"mpi_allreduce",                 (DL_FUNC) &mpi_allreduce,                   4},
    {"mpi_iprobe",                    (DL_FUNC) &mpi_iprobe,                      4},
    {"mpi_probe",                     (DL_FUNC) &mpi_probe,                       4},
    {"mpi_get_count",                 (DL_FUNC) &mpi_get_count,                   2},
    {"mpi_get_sourcetag",             (DL_FUNC) &mpi_get_sourcetag,               1},
    {"mpi_barrier",                   (DL_FUNC) &mpi_barrier,                     1},
    {"mpi_comm_is_null",              (DL_FUNC) &mpi_comm_is_null,                1},
    {"mpi_comm_size",                 (DL_FUNC) &mpi_comm_size,                   1},
    {"mpi_comm_rank",                 (DL_FUNC) &mpi_comm_rank,                   1},
    {"mpi_comm_dup",                  (DL_FUNC) &mpi_comm_dup,                    2},
    {"mpi_comm_c2f",                  (DL_FUNC) &mpi_comm_c2f,                    1},
    {"mpi_comm_free",                 (DL_FUNC) &mpi_comm_free,                   1},
    {"mpi_abort",                     (DL_FUNC) &mpi_abort,                       1},
    {"mpi_comm_set_errhandler",       (DL_FUNC) &mpi_comm_set_errhandler,         1},
    {"mpi_comm_test_inter",           (DL_FUNC) &mpi_comm_test_inter,             1},
#ifdef MPI2
    {"mpi_comm_spawn",                (DL_FUNC) &mpi_comm_spawn,                  7},
#endif
    {"mpi_comm_get_parent",           (DL_FUNC) &mpi_comm_get_parent,             1},
    {"mpi_is_master",                 (DL_FUNC) &mpi_is_master,                   0},
    {"mpi_comm_disconnect",           (DL_FUNC) &mpi_comm_disconnect,             1},
    {"mpi_intercomm_merge",           (DL_FUNC) &mpi_intercomm_merge,             3},
    {"mpi_comm_remote_size",          (DL_FUNC) &mpi_comm_remote_size,            1},
    {"mpi_sendrecv",                  (DL_FUNC) &mpi_sendrecv,                   10},
    {"mpi_sendrecv_replace",          (DL_FUNC) &mpi_sendrecv_replace,            8},
    {"mpi_cart_create",               (DL_FUNC) &mpi_cart_create,                 5},
    {"mpi_dims_create",               (DL_FUNC) &mpi_dims_create,                 3},
    {"mpi_cartdim_get",               (DL_FUNC) &mpi_cartdim_get,                 1},
    {"mpi_cart_get",                  (DL_FUNC) &mpi_cart_get,                    2},
    {"mpi_cart_rank",                 (DL_FUNC) &mpi_cart_rank,                   2},
    {"mpi_cart_coords",               (DL_FUNC) &mpi_cart_coords,                 3},
    {"mpi_cart_shift",                (DL_FUNC) &mpi_cart_shift,                  3},
    {"mpi_isend",                     (DL_FUNC) &mpi_isend,                       6},
    {"mpi_irecv",                     (DL_FUNC) &mpi_irecv,                       6},
    {"mpi_wait",                      (DL_FUNC) &mpi_wait,                        2},
    {"mpi_test",                      (DL_FUNC) &mpi_test,                        2},
    {"mpi_cancel",                    (DL_FUNC) &mpi_cancel,                      1},
    {"mpi_test_cancelled",            (DL_FUNC) &mpi_test_cancelled,              1},
    {"mpi_waitany",                   (DL_FUNC) &mpi_waitany,                     2},
    {"mpi_testany",                   (DL_FUNC) &mpi_testany,                     2},
    {"mpi_waitall",                   (DL_FUNC) &mpi_waitall,                     1},
    {"mpi_testall",                   (DL_FUNC) &mpi_testall,                     1},
    {"mpi_testsome",                  (DL_FUNC) &mpi_testsome,                    1},
    {"mpi_waitsome",                  (DL_FUNC) &mpi_waitsome,                    1},
    {NULL, NULL, 0}
};

void R_init_npRmpi(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
