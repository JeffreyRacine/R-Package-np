/* Copyright (C) Jeff Racine, 1995-2004 */

#include <float.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <R.h>
#include <R_ext/Arith.h>
#include <Rinternals.h>
#include <stdint.h>

#ifdef MPI2
#include "mpi.h"
int my_rank;
int iNum_Processors;
int source;
int dest;
int tag;
int iSeed_my_rank;
MPI_Status status;
extern MPI_Comm	*comm;
#endif

/*
-f to enable Fast resampling (memory intensive)
-n to use Nonparametric measures of central tendency and dispersion
-o to use Ordered categorical gradients
-p to compute univariate categorical conditional Predictions
-r to enable a different Random seed with each program invocation

*/

/* headers.h has all definitions of routines used by main() and related modules */

#include "headers.h"
#include "tree.h"

// categorical hashing
#include "hash.h"


int int_DEBUG;
int int_VERBOSE;
int int_TAYLOR;
int int_WEIGHTS = 0;
int int_LARGE_SF;
int int_NOKEYPRESS;
int int_DISPLAY_CV;
int int_RANDOM_SEED = 42;
int int_MINIMIZE_IO;
int int_ORDERED_CATEGORICAL_GRADIENT;
int int_PREDICT;
int int_ROBUST;
int int_SIMULATION;

int int_RESTART_FROM_MIN;

int int_TREE_X;
int int_TREE_Y;
int int_TREE_XY;
int int_nn_k_min_extern = 1;

static int np_mpi_local_regression_mode = 0;
#ifdef MPI2
static int np_mpi_local_regression_saved_rank = 0;
static int np_mpi_local_regression_saved_nproc = 1;
static MPI_Comm np_mpi_local_regression_saved_comm1 = MPI_COMM_NULL;
#endif
static int fit_progress_active = 0;
static int fit_progress_total = 0;
static int fit_progress_offset = 0;
static clock_t fit_progress_last_signal_clock = 0;
static int fit_progress_last_signal_eval = 0;

static void np_progress_signal(const char *event, const char *surface, const int current, const int total)
{
  SEXP ns = R_NilValue;
  SEXP fn = R_NilValue;
  SEXP event_s = R_NilValue;
  SEXP surface_s = R_NilValue;
  SEXP current_s = R_NilValue;
  SEXP total_s = R_NilValue;
  SEXP call = R_NilValue;
  int err = 0;

  if (event == NULL || surface == NULL)
    return;

  PROTECT(ns = R_FindNamespace(Rf_ScalarString(Rf_mkChar("npRmpi"))));
  if (ns == R_NilValue) {
    UNPROTECT(1);
    return;
  }

  if (!R_existsVarInFrame(ns, Rf_install(".np_progress_signal_from_c"))) {
    UNPROTECT(1);
    return;
  }

  PROTECT(fn = Rf_findFun(Rf_install(".np_progress_signal_from_c"), ns));

  PROTECT(event_s = Rf_mkString(event));
  PROTECT(surface_s = Rf_mkString(surface));
  PROTECT(current_s = Rf_ScalarInteger(current));
  PROTECT(total_s = Rf_ScalarInteger(total));
  PROTECT(call = Rf_lang5(fn, event_s, surface_s, current_s, total_s));
  R_tryEval(call, ns, &err);
  UNPROTECT(7);
}

SEXP C_np_progress_signal(SEXP event, SEXP surface, SEXP current, SEXP total)
{
  const char *event_c = NULL;
  const char *surface_c = NULL;
  int current_i = 0;
  int total_i = 0;

  if (TYPEOF(event) != STRSXP || XLENGTH(event) < 1 ||
      TYPEOF(surface) != STRSXP || XLENGTH(surface) < 1)
    error("C_np_progress_signal: event and surface must be character scalars");

  event_c = CHAR(STRING_ELT(event, 0));
  surface_c = CHAR(STRING_ELT(surface, 0));
  current_i = Rf_asInteger(current);
  total_i = Rf_asInteger(total);

  np_progress_signal(event_c, surface_c, current_i, total_i);

  return R_NilValue;
}

static void np_progress_bandwidth_multistart_step(const int done, const int total)
{
  if (done < 1 || total <= 1)
    return;

  np_progress_signal("bandwidth_multistart_step", "bandwidth", done, total);
}

static void np_progress_bandwidth_activity_step(const int done)
{
  if (done < 1)
    return;

  np_progress_signal("bandwidth_activity_step", "bandwidth", done, 0);
}

static void np_progress_fit_clear_state(void)
{
  fit_progress_active = 0;
  fit_progress_total = 0;
  fit_progress_offset = 0;
  fit_progress_last_signal_clock = 0;
  fit_progress_last_signal_eval = 0;
}

static void np_progress_fit_activate(const int total)
{
  np_progress_fit_clear_state();

  if (total < 1)
    return;

  fit_progress_active = 1;
  fit_progress_total = total;
}

void np_progress_fit_set_offset(const int offset)
{
  if (!fit_progress_active)
    return;

  fit_progress_offset = MAX(0, offset);
  fit_progress_last_signal_clock = 0;
  fit_progress_last_signal_eval = fit_progress_offset;
}

static void np_progress_fit_maybe_signal(const int global_done)
{
  const clock_t now = clock();
  const double signal_after_sec = 0.5;
  const int signal_every_evals = 64;
  double since_last = 0.0;
  int bounded_done = global_done;

  if (!fit_progress_active || fit_progress_total < 1 || global_done < 1)
    return;

  if (bounded_done > fit_progress_total)
    bounded_done = fit_progress_total;

  if (bounded_done < fit_progress_total) {
    if (bounded_done <= fit_progress_last_signal_eval)
      return;

    if (fit_progress_last_signal_eval > 0) {
      if ((fit_progress_last_signal_clock > 0) && (now > fit_progress_last_signal_clock)) {
        since_last = ((double)(now - fit_progress_last_signal_clock)) / (double)CLOCKS_PER_SEC;

        if ((bounded_done - fit_progress_last_signal_eval) < signal_every_evals &&
            since_last < signal_after_sec)
          return;
      }
    }
  }

  np_progress_signal("fit_step", "bandwidth", bounded_done, fit_progress_total);
  fit_progress_last_signal_eval = bounded_done;
  fit_progress_last_signal_clock = now;
}

void np_progress_fit_step(const int done)
{
  np_progress_fit_maybe_signal(done);
}

void np_progress_fit_loop_step(const int done, const int natural_total)
{
  if (!fit_progress_active || natural_total <= 1)
    return;

  if ((fit_progress_offset + natural_total) > fit_progress_total)
    return;

  np_progress_fit_maybe_signal(fit_progress_offset + done);
}

int np_progress_fit_is_active(void)
{
  return fit_progress_active;
}

SEXP C_np_progress_fit_begin(SEXP total)
{
  np_progress_fit_activate(Rf_asInteger(total));
  return R_NilValue;
}

SEXP C_np_progress_fit_end(void)
{
  np_progress_fit_clear_state();
  return R_NilValue;
}

/* Some externals for numerical routines */
/* Some externals for numerical routines */

int num_obs_train_extern=0;
int num_obs_eval_extern=0;
int num_var_continuous_extern=0;
int num_var_unordered_extern=0;
int num_var_ordered_extern=0;
int num_reg_continuous_extern=0;
int num_reg_unordered_extern=0;
int num_reg_ordered_extern=0;

int int_cker_bound_extern=0;
double *vector_ckerlb_extern=NULL;
double *vector_ckerub_extern=NULL;

int np_mpi_local_regression_active(void)
{
  return np_mpi_local_regression_mode;
}
int int_cxker_bound_extern=0;
int int_cyker_bound_extern=0;
int int_cxyker_bound_extern=0;
double *vector_cxkerlb_extern=NULL;
double *vector_cxkerub_extern=NULL;
double *vector_cykerlb_extern=NULL;
double *vector_cykerub_extern=NULL;
double *vector_cxykerlb_extern=NULL;
double *vector_cxykerub_extern=NULL;


int *num_categories_extern;
double **matrix_categorical_vals_extern;

double **matrix_X_continuous_train_extern;
double **matrix_X_unordered_train_extern;
double **matrix_X_ordered_train_extern;
double **matrix_X_continuous_eval_extern;
double **matrix_X_unordered_eval_extern;
double **matrix_X_ordered_eval_extern;

double **matrix_Y_continuous_train_extern;
double **matrix_Y_unordered_train_extern;
double **matrix_Y_ordered_train_extern;

double **matrix_Y_continuous_eval_extern;
double **matrix_Y_unordered_eval_extern;
double **matrix_Y_ordered_eval_extern;

/* these are data which are sorted into an 'alternate' order */
/* this allows us to support 2 trees simultaneously !*/

int * num_categories_extern_XY;
int * num_categories_extern_X;
int * num_categories_extern_Y;

double **matrix_categorical_vals_extern_XY;
double **matrix_categorical_vals_extern_X;
double **matrix_categorical_vals_extern_Y;

double **matrix_XY_continuous_train_extern;
double **matrix_XY_unordered_train_extern;
double **matrix_XY_ordered_train_extern;
double **matrix_XY_continuous_eval_extern;
double **matrix_XY_unordered_eval_extern;
double **matrix_XY_ordered_eval_extern;

int * ipt_extern_X;
int * ipt_extern_Y;
int * ipt_extern_XY;

static double (*bwmfunc_raw)(double *) = NULL;
static double bwm_eval_count = 0.0;
static double bwm_invalid_count = 0.0;
static double bwm_fast_eval_count = 0.0;
static clock_t bwm_progress_started_clock = 0;
static clock_t bwm_progress_last_signal_clock = 0;
static int bwm_progress_last_signal_eval = 0;

typedef struct {
  int active;
  int owned;
  int num_var;
  int num_reg_continuous;
  int num_reg_unordered;
  int num_reg_ordered;
  int *ipt;
  int *glp_degree;
  double *ckerlb;
  double *ckerub;
  double *vector_scale_factor;
} NPRegressionNomadShadowCtx;

static NPRegressionNomadShadowCtx np_regression_nomad_shadow =
  {0, 0, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL};

typedef struct {
  int active;
  int owned;
  int num_all_var;
  int num_reg_continuous;
  int num_var_continuous;
  int num_reg_unordered;
  int num_var_unordered;
  int num_reg_ordered;
  int num_var_ordered;
  int need_y_side;
  int old_cdens;
  int penalty_mode;
  double penalty_multiplier;
  int *glp_degree;
  int *ipt_x;
  int *ipt_y;
  int *ipt_xy;
  int *ipt_lookup_x;
  int *ipt_lookup_y;
  int *ipt_lookup_xy;
  double *vector_scale_factor;
  double *cxkerlb;
  double *cxkerub;
  double *cykerlb;
  double *cykerub;
  double *cxykerlb;
  double *cxykerub;
} NPConditionalDensityNomadShadowCtx;

static NPConditionalDensityNomadShadowCtx np_conditional_density_nomad_shadow =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, NULL, NULL, NULL, NULL, NULL, NULL,
   NULL, NULL, NULL, NULL, NULL, NULL, NULL};

static void np_regression_nomad_shadow_clear_internal(void);
static void np_conditional_density_nomad_shadow_clear_internal(void);
static void np_conditional_density_nomad_shadow_refresh_penalty(void);

static int np_mpi_comm1_all_ok(const int local_ok)
{
#ifdef MPI2
  int global_ok = local_ok ? 1 : 0;
  if (comm[1] != MPI_COMM_NULL)
    MPI_Allreduce(MPI_IN_PLACE, &global_ok, 1, MPI_INT, MPI_MIN, comm[1]);
  return global_ok;
#else
  return local_ok ? 1 : 0;
#endif
}

static int np_mpi_rank_failure_injected(const char *env_name)
{
#ifdef MPI2
  const char *value = getenv(env_name);
  char *end = NULL;
  long rank;

  if ((value == NULL) || (value[0] == '\0') || strcmp(value, "0") == 0)
    return 0;
  if (strcmp(value, "all") == 0)
    return 1;

  errno = 0;
  rank = strtol(value, &end, 10);
  if ((errno != 0) || (end == value) || (*end != '\0'))
    return 0;

  return ((int)rank == my_rank);
#else
  (void)env_name;
  return 0;
#endif
}

static void resolve_bounds_or_default(SEXP lb_r, SEXP ub_r, int ncon, double ** lb_p, double ** ub_p)
{
  int i;
  *lb_p = REAL(lb_r);
  *ub_p = REAL(ub_r);
  if ((XLENGTH(lb_r) == 0 || XLENGTH(ub_r) == 0) && ncon > 0) {
    *lb_p = (double *)R_alloc((size_t)ncon, sizeof(double));
    *ub_p = (double *)R_alloc((size_t)ncon, sizeof(double));
    for (i = 0; i < ncon; i++) {
      (*lb_p)[i] = -INFINITY;
      (*ub_p)[i] = INFINITY;
    }
  }
}

static void np_reset_y_side_extern(void)
{
  num_var_unordered_extern = 0;
  num_var_ordered_extern = 0;
  num_var_continuous_extern = 0;

  matrix_Y_unordered_train_extern = NULL;
  matrix_Y_ordered_train_extern = NULL;
  matrix_Y_continuous_train_extern = NULL;
  matrix_Y_unordered_eval_extern = NULL;
  matrix_Y_ordered_eval_extern = NULL;
  matrix_Y_continuous_eval_extern = NULL;

  num_categories_extern_Y = NULL;
  matrix_categorical_vals_extern_Y = NULL;
}

static void np_shadow_matrix_dims(SEXP x, int *nrow, int *ncol)
{
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(dim == R_NilValue || XLENGTH(dim) != 2){
    *nrow = 0;
    *ncol = 0;
    return;
  }
  *nrow = INTEGER(dim)[0];
  *ncol = INTEGER(dim)[1];
}

static void np_shadow_fill_matrix(double **dest, const double *src, int nrow, int ncol)
{
  int i, j;
  for(j = 0; j < ncol; j++)
    for(i = 0; i < nrow; i++)
      dest[j][i] = src[(size_t)j * (size_t)nrow + (size_t)i];
}

static int bwm_use_transform = 0;
static int bwm_num_reg_continuous = 0;
static int bwm_num_reg_unordered = 0;
static int bwm_num_reg_ordered = 0;
static int bwm_kernel_unordered = 0;
static int *bwm_num_categories = NULL;
static double *bwm_transform_buf = NULL;
static int bwm_transform_buf_len = 0;

static void bwm_reserve_transform_buf(int needed_len)
{
  double *tmp = NULL;
  if (needed_len <= bwm_transform_buf_len)
    return;
  tmp = (double *) realloc(bwm_transform_buf, (size_t) needed_len * sizeof(double));
  if (tmp == NULL)
    error("bwm_reserve_transform_buf: memory allocation failed");
  bwm_transform_buf = tmp;
  bwm_transform_buf_len = needed_len;
}

void np_release_static_buffers(int *unused)
{
  (void)unused;
  if (bwm_transform_buf != NULL) {
    free(bwm_transform_buf);
    bwm_transform_buf = NULL;
  }
  bwm_transform_buf_len = 0;
}
static int bwm_penalty_mode = 0;
static double bwm_penalty_value = DBL_MAX;
static int *bwm_kernel_unordered_vec = NULL;
static int bwm_kernel_unordered_len = 0;
static double bwm_scale_factor_lower_bound = 0.1;
static const char *bwm_deferred_error = NULL;

void np_bwm_set_deferred_error(const char *msg)
{
  bwm_deferred_error = msg;
}

const char *np_bwm_get_deferred_error(void)
{
  return bwm_deferred_error;
}

void np_bwm_clear_deferred_error(void)
{
  bwm_deferred_error = NULL;
}

static void bwm_set_scale_factor_lower_bound(double value)
{
  bwm_scale_factor_lower_bound =
    (R_FINITE(value) && (value >= 0.0)) ? value : 0.1;
}

static int np_has_finite_cker_bounds(const double *lb, const double *ub, const int n)
{
  int i;
  const double big = 0.5*DBL_MAX;
  if(lb == NULL || ub == NULL || n <= 0)
    return 0;
  for(i = 0; i < n; i++) {
    const int lb_fin = isfinite(lb[i]) && (fabs(lb[i]) < big);
    const int ub_fin = isfinite(ub[i]) && (fabs(ub[i]) < big);
    if(lb_fin || ub_fin)
      return 1;
  }
  return 0;
}

static int np_cmp_desc_int(const void *a, const void *b)
{
  const int ia = *((const int *)a);
  const int ib = *((const int *)b);
  return (ib > ia) - (ib < ia);
}

static void np_dim_basis_two_dimen(const int d1, const int d2, double *nd1, double *pd12)
{
  int i, j, low;
  double d12 = *pd12;
  double s;
  double *nd2 = NULL;

  if ((d1 <= 0) || (d2 <= 0))
    return;

  if (d2 == 1){
    *pd12 = d12;
    return;
  }

  d12 = (double)d2;
  if (d1 > d2){
    for (i = 1; i <= (d1 - d2); i++)
      d12 += ((double)d2) * nd1[i];
  }

  for (i = 2; i <= d2; i++)
    d12 += ((double)i) * nd1[d1 - i + 1];

  d12 += nd1[d1];

  nd2 = (double *)malloc((size_t)(d1 + 1) * sizeof(double));
  if (nd2 == NULL)
    error("np_dim_basis_two_dimen: memory allocation failed");

  for (i = 0; i <= d1; i++)
    nd2[i] = nd1[i];

  if (d1 > 1){
    for (j = 1; j <= (d1 - 1); j++){
      s = 0.0;
      low = MAX(0, j - d2 + 1);
      for (i = j; i >= low; i--)
        s += (i > 0) ? nd1[i] : 1.0;
      nd2[j] = s;
    }
  }

  nd2[d1] = nd1[d1];
  for (i = MAX(1, d1 - d2 + 1); i <= (d1 - 1); i++)
    nd2[d1] += nd1[i];

  for (i = 0; i <= d1; i++)
    nd1[i] = nd2[i];

  free(nd2);
  *pd12 = d12;
}

void np_dim_basis(int *basis_code,
                  int *kernel,
                  int *degree,
                  int *segments,
                  int *k,
                  int *include,
                  int *categories,
                  int *ninclude,
                  double *result)
{
  int i, m, dsum, rcount = 0;
  int basis;
  int use_kernel;
  int ndeg;
  int ninc;
  int *rows = NULL;
  int *dims = NULL;
  double ncol_bs = 0.0;

  if ((basis_code == NULL) || (kernel == NULL) || (k == NULL) || (result == NULL)){
    if (result != NULL) *result = NA_REAL;
    return;
  }

  basis = *basis_code;      /* 0=additive, 1=glp, 2=tensor */
  use_kernel = *kernel;
  ndeg = MAX(0, *k);
  ninc = MAX(0, (ninclude == NULL) ? 0 : *ninclude);

  if ((basis < 0) || (basis > 2)){
    *result = NA_REAL;
    return;
  }

  if (ndeg > 0){
    rows = (int *)malloc((size_t)ndeg * sizeof(int));
    if (rows == NULL)
      error("np_dim_basis: memory allocation failed");
    for (i = 0; i < ndeg; i++)
      rows[i] = degree[i] + segments[i];
  }

  if (use_kernel){
    if (basis == 0){ /* additive */
      for (i = 0; i < ndeg; i++)
        if (degree[i] > 0)
          ncol_bs += (double)(rows[i] - 1);
    } else if (basis == 2){ /* tensor */
      ncol_bs = 1.0;
      m = 0;
      for (i = 0; i < ndeg; i++){
        if (degree[i] > 0){
          ncol_bs *= (double)rows[i];
          m++;
        }
      }
      if (m == 0)
        ncol_bs = 0.0;
    } else { /* glp */
      if (ndeg > 0){
        dims = (int *)malloc((size_t)ndeg * sizeof(int));
        if (dims == NULL)
          error("np_dim_basis: memory allocation failed");
      }
      for (i = 0; i < ndeg; i++){
        if (degree[i] > 0){
          dsum = rows[i] - 1;
          if (dsum > 0)
            dims[rcount++] = dsum;
        }
      }
    }
  } else {
    if (basis == 0){ /* additive */
      for (i = 0; i < ndeg; i++)
        if (degree[i] > 0)
          ncol_bs += (double)(rows[i] - 1);
      for (i = 0; i < ninc; i++)
        ncol_bs += (double)(include[i] * categories[i] - 1);
    } else if (basis == 2){ /* tensor */
      ncol_bs = 1.0;
      m = 0;
      for (i = 0; i < ndeg; i++){
        if (degree[i] > 0){
          ncol_bs *= (double)rows[i];
          m++;
        }
      }
      for (i = 0; i < ninc; i++){
        ncol_bs *= (double)(include[i] * categories[i] - 1);
        m++;
      }
      if (m == 0)
        ncol_bs = 0.0;
    } else { /* glp */
      if ((ndeg + ninc) > 0){
        dims = (int *)malloc((size_t)(ndeg + ninc) * sizeof(int));
        if (dims == NULL)
          error("np_dim_basis: memory allocation failed");
      }
      for (i = 0; i < ndeg; i++){
        if (degree[i] > 0){
          dsum = rows[i] - 1;
          if (dsum > 0)
            dims[rcount++] = dsum;
        }
      }
      for (i = 0; i < ninc; i++){
        dsum = include[i] * categories[i] - 1;
        if (dsum > 0)
          dims[rcount++] = dsum;
      }
    }
  }

  if (basis == 1){ /* glp recurrence */
    if (rcount == 0){
      ncol_bs = 0.0;
    } else {
      int dmax, idx;
      double *nd1 = NULL;
      qsort(dims, (size_t)rcount, sizeof(int), np_cmp_desc_int);
      dmax = dims[0];
      nd1 = (double *)malloc((size_t)(dmax + 1) * sizeof(double));
      if (nd1 == NULL)
        error("np_dim_basis: memory allocation failed");
      for (idx = 0; idx <= dmax; idx++)
        nd1[idx] = (idx == dmax) ? 0.0 : 1.0;
      ncol_bs = (double)dmax;
      if (rcount > 1){
        for (idx = 1; idx < rcount; idx++)
          np_dim_basis_two_dimen(dmax, dims[idx], nd1, &ncol_bs);
        ncol_bs += (double)(rcount - 1);
      }
      free(nd1);
    }
  }

  if (rows != NULL) free(rows);
  if (dims != NULL) free(dims);
  *result = ncol_bs;
}

static void bwm_reset_counters(void)
{
  bwm_eval_count = 0.0;
  bwm_invalid_count = 0.0;
  bwm_fast_eval_count = 0.0;
  bwm_progress_started_clock = clock();
  bwm_progress_last_signal_clock = bwm_progress_started_clock;
  bwm_progress_last_signal_eval = 0;
  np_fastcv_alllarge_hits_reset();
}

static void bwm_maybe_signal_activity(void)
{
  const int current_eval = (int) bwm_eval_count;
  const clock_t now = clock();
  const double signal_after_sec = 0.5;
  const int signal_every_evals = 64;
  double since_last = 0.0;

  if (current_eval < 1)
    return;

  if ((bwm_progress_last_signal_clock > 0) && (now > bwm_progress_last_signal_clock))
    since_last = ((double) (now - bwm_progress_last_signal_clock)) / (double) CLOCKS_PER_SEC;

  if ((current_eval - bwm_progress_last_signal_eval) < signal_every_evals &&
      since_last < signal_after_sec)
    return;

  np_progress_bandwidth_activity_step(current_eval);
  bwm_progress_last_signal_eval = current_eval;
  bwm_progress_last_signal_clock = now;
}

static inline void bwm_snapshot_fast_counters(void)
{
  bwm_fast_eval_count = np_fastcv_alllarge_hits_get();
}

static double bwm_sigmoid(double x)
{
  if (x >= 0.0) {
    double z = exp(-x);
    return 1.0/(1.0+z);
  } else {
    double z = exp(x);
    return z/(1.0+z);
  }
}

static double bwm_safe_exp(double x)
{
  if (x > 700.0) return DBL_MAX/2.0;
  if (x < -700.0) return 0.0;
  return exp(x);
}

static double bwm_logit(double p)
{
  const double eps = 1e-12;
  if (p < eps) p = eps;
  if (p > 1.0 - eps) p = 1.0 - eps;
  return log(p/(1.0-p));
}

static void bwm_apply_transform(const double *p, double *out, int n)
{
  int i;
  int idx;

  /* copy */
  for (i = 1; i <= n; i++)
    out[i] = p[i];

  if (!bwm_use_transform)
    return;

  /* continuous (positive) */
  for (i = 1; i <= bwm_num_reg_continuous; i++)
    out[i] = bwm_safe_exp(p[i]);

  /* unordered categorical */
  for (i = 0; i < bwm_num_reg_unordered; i++) {
    idx = bwm_num_reg_continuous + 1 + i;
    if (bwm_num_categories != NULL) {
    int kern = (bwm_kernel_unordered_vec != NULL && i < bwm_kernel_unordered_len) ?
      bwm_kernel_unordered_vec[i] : bwm_kernel_unordered;
    double maxbw = max_unordered_bw(bwm_num_categories[i], kern);
      out[idx] = bwm_sigmoid(p[idx]) * maxbw;
    } else {
      out[idx] = bwm_sigmoid(p[idx]);
    }
  }

  /* ordered categorical (0..1) */
  for (i = 0; i < bwm_num_reg_ordered; i++) {
    idx = bwm_num_reg_continuous + bwm_num_reg_unordered + 1 + i;
    out[idx] = bwm_sigmoid(p[idx]);
  }
}

static void bwm_to_unconstrained(double *p, int n)
{
  int i;
  int idx;
  const double eps = 1e-12;

  if (!bwm_use_transform)
    return;

  for (i = 1; i <= n; i++)
    bwm_transform_buf[i] = p[i];

  /* continuous */
  for (i = 1; i <= bwm_num_reg_continuous; i++) {
    double v = bwm_transform_buf[i];
    if (v <= 0.0) v = eps;
    p[i] = log(v);
  }

  /* unordered */
  for (i = 0; i < bwm_num_reg_unordered; i++) {
    idx = bwm_num_reg_continuous + 1 + i;
    int kern = (bwm_kernel_unordered_vec != NULL && i < bwm_kernel_unordered_len) ?
      bwm_kernel_unordered_vec[i] : bwm_kernel_unordered;
    double maxbw = (bwm_num_categories != NULL) ? max_unordered_bw(bwm_num_categories[i], kern) : 1.0;
    double v = bwm_transform_buf[idx];
    if (maxbw <= 0.0) {
      p[idx] = 0.0;
    } else {
      double frac = v / maxbw;
      p[idx] = bwm_logit(frac);
    }
  }

  /* ordered */
  for (i = 0; i < bwm_num_reg_ordered; i++) {
    idx = bwm_num_reg_continuous + bwm_num_reg_unordered + 1 + i;
    p[idx] = bwm_logit(bwm_transform_buf[idx]);
  }
}

static void bwm_to_constrained(double *p, int n)
{
  int i;
  if (!bwm_use_transform)
    return;

  bwm_apply_transform(p, bwm_transform_buf, n);
  for (i = 1; i <= n; i++)
    p[i] = bwm_transform_buf[i];
}

static int np_bw_candidate_is_admissible_with_floor(
  int n,
  int use_transform,
  int KERNEL,
  int KERNEL_unordered_liracine,
  int BANDWIDTH,
  int BANDWIDTH_den_ml,
  int REGRESSION_ML,
  int num_obs,
  int num_var_continuous,
  int num_var_unordered,
  int num_var_ordered,
  int num_reg_continuous,
  int num_reg_unordered,
  int num_reg_ordered,
  int *num_categories,
  double *vector_scale_factor,
  double floor_coeff);

static int bwm_floor_context_active = 0;
static int bwm_floor_context_n = 0;
static int bwm_floor_context_use_transform = 0;
static int bwm_floor_context_kernel = 0;
static int bwm_floor_context_kernel_unordered = 0;
static int bwm_floor_context_bandwidth = 0;
static int bwm_floor_context_bandwidth_den_ml = 0;
static int bwm_floor_context_regression_ml = 0;
static int bwm_floor_context_num_obs = 0;
static int bwm_floor_context_num_var_continuous = 0;
static int bwm_floor_context_num_var_unordered = 0;
static int bwm_floor_context_num_var_ordered = 0;
static int bwm_floor_context_num_reg_continuous = 0;
static int bwm_floor_context_num_reg_unordered = 0;
static int bwm_floor_context_num_reg_ordered = 0;
static int *bwm_floor_context_num_categories = NULL;
static double bwm_floor_context_coeff = 0.1;

static void bwm_clear_floor_context(void)
{
  bwm_floor_context_active = 0;
  bwm_floor_context_num_categories = NULL;
}

static void bwm_set_floor_context(
  int active,
  int n,
  int use_transform,
  int KERNEL,
  int KERNEL_unordered_liracine,
  int BANDWIDTH,
  int BANDWIDTH_den_ml,
  int REGRESSION_ML,
  int num_obs,
  int num_var_continuous,
  int num_var_unordered,
  int num_var_ordered,
  int num_reg_continuous,
  int num_reg_unordered,
  int num_reg_ordered,
  int *num_categories,
  double floor_coeff)
{
  bwm_floor_context_active = active;
  bwm_floor_context_n = n;
  bwm_floor_context_use_transform = use_transform;
  bwm_floor_context_kernel = KERNEL;
  bwm_floor_context_kernel_unordered = KERNEL_unordered_liracine;
  bwm_floor_context_bandwidth = BANDWIDTH;
  bwm_floor_context_bandwidth_den_ml = BANDWIDTH_den_ml;
  bwm_floor_context_regression_ml = REGRESSION_ML;
  bwm_floor_context_num_obs = num_obs;
  bwm_floor_context_num_var_continuous = num_var_continuous;
  bwm_floor_context_num_var_unordered = num_var_unordered;
  bwm_floor_context_num_var_ordered = num_var_ordered;
  bwm_floor_context_num_reg_continuous = num_reg_continuous;
  bwm_floor_context_num_reg_unordered = num_reg_unordered;
  bwm_floor_context_num_reg_ordered = num_reg_ordered;
  bwm_floor_context_num_categories = num_categories;
  bwm_floor_context_coeff = floor_coeff;
}

static int bwm_active_floor_candidate_ok(double *p)
{
  if (!bwm_floor_context_active)
    return 1;

  return np_bw_candidate_is_admissible_with_floor(
    bwm_floor_context_n,
    bwm_floor_context_use_transform,
    bwm_floor_context_kernel,
    bwm_floor_context_kernel_unordered,
    bwm_floor_context_bandwidth,
    bwm_floor_context_bandwidth_den_ml,
    bwm_floor_context_regression_ml,
    bwm_floor_context_num_obs,
    bwm_floor_context_num_var_continuous,
    bwm_floor_context_num_var_unordered,
    bwm_floor_context_num_var_ordered,
    bwm_floor_context_num_reg_continuous,
    bwm_floor_context_num_reg_unordered,
    bwm_floor_context_num_reg_ordered,
    bwm_floor_context_num_categories,
    p,
    bwm_floor_context_coeff);
}

static double bwmfunc_wrapper(double *p)
{
  double val;
  double *use_p = p;

  bwm_eval_count += 1.0;
  if (!bwm_active_floor_candidate_ok(p)) {
    bwm_invalid_count += 1.0;
    if (bwm_penalty_mode == 1 && R_FINITE(bwm_penalty_value))
      return bwm_penalty_value;
    return DBL_MAX;
  }

  if (bwm_use_transform) {
    int n = bwm_num_reg_continuous + bwm_num_reg_unordered + bwm_num_reg_ordered;
    bwm_reserve_transform_buf(n + 1);
    bwm_apply_transform(p, bwm_transform_buf, n);
    use_p = bwm_transform_buf;
  }

  val = bwmfunc_raw(use_p);
  bwm_maybe_signal_activity();

  if (!R_FINITE(val) || val == DBL_MAX) {
    bwm_invalid_count += 1.0;
    if (bwm_penalty_mode == 1 && R_FINITE(bwm_penalty_value))
      return bwm_penalty_value;
  }

  return val;
}

static void np_copy_scale_factor(double *dest, const double *src, int n)
{
  int i;

  for (i = 1; i <= n; i++)
    dest[i] = src[i];
}

static void np_copy_scale_factor_for_raw(double *dest, const double *src, int n)
{
  np_copy_scale_factor(dest, src, n);
  if (bwm_use_transform)
    bwm_to_constrained(dest, n);
}

static double bwmfunc_raw_current_scale(double *vector_scale_factor, int n)
{
  double val;

  if (!bwm_use_transform)
    return bwmfunc_raw(vector_scale_factor);

  double *tmp = alloc_vecd(n + 1);
  np_copy_scale_factor_for_raw(tmp, vector_scale_factor, n);
  val = bwmfunc_raw(tmp);
  safe_free(tmp);
  return val;
}

extern double *vector_continuous_stddev_extern;
extern double nconfac_extern;

static double np_fixed_continuous_floor_cv_with_coeff(const int continuous_index,
                                                      const double temp_pow,
                                                      const double floor_coeff)
{
  const double scale_pow =
    (R_FINITE(nconfac_extern) && (nconfac_extern > 0.0)) ? nconfac_extern : temp_pow;

  if (int_LARGE_SF == SF_NORMAL)
    return floor_coeff;

  if ((vector_continuous_stddev_extern != NULL) &&
      R_FINITE(vector_continuous_stddev_extern[continuous_index]) &&
      (vector_continuous_stddev_extern[continuous_index] > 0.0)) {
    return floor_coeff * vector_continuous_stddev_extern[continuous_index] * scale_pow;
  }

  return floor_coeff * scale_pow;
}

static double np_fixed_continuous_temp_pow_cv(
  int KERNEL,
  int num_obs,
  int num_var_continuous,
  int num_reg_continuous)
{
  const double dim = (double) num_reg_continuous + num_var_continuous;

  switch (KERNEL) {
    case 0:
    case 4:
    case 8:
      return 1.0 / pow((double) num_obs, 2.0 / (4.0 + dim));
    case 1:
    case 5:
      return 1.0 / pow((double) num_obs, 2.0 / (8.0 + dim));
    case 2:
    case 6:
      return 1.0 / pow((double) num_obs, 2.0 / (12.0 + dim));
    case 3:
    case 7:
      return 1.0 / pow((double) num_obs, 2.0 / (16.0 + dim));
    default:
      return DBL_MAX;
  }
}

static int np_fixed_continuous_below_floor_cv(double candidate, double floor)
{
  /* NOMAD/Powell handoffs can land exactly on this floor up to roundoff. */
  const double scale = fmax(1.0, fmax(fabs(candidate), fabs(floor)));
  const double tol = 64.0 * DBL_EPSILON * scale;

  return candidate < (floor - tol);
}

static int np_fixed_continuous_floor_ok_cv_with_coeff(
  int KERNEL,
  int num_obs,
  int num_var_continuous,
  int num_reg_continuous,
  double *candidate,
  double floor_coeff)
{
  int i;
  const double temp_pow = np_fixed_continuous_temp_pow_cv(
    KERNEL,
    num_obs,
    num_var_continuous,
    num_reg_continuous);
  const int total_continuous = num_reg_continuous + num_var_continuous;

  for (i = 1; i <= total_continuous; i++) {
    const double bw_floor = np_fixed_continuous_floor_cv_with_coeff(i - 1, temp_pow, floor_coeff);
    if (np_fixed_continuous_below_floor_cv(candidate[i], bw_floor))
      return 0;
  }

  return 1;
}

static int np_bw_candidate_is_admissible_with_floor(
  int n,
  int use_transform,
  int KERNEL,
  int KERNEL_unordered_liracine,
  int BANDWIDTH,
  int BANDWIDTH_den_ml,
  int REGRESSION_ML,
  int num_obs,
  int num_var_continuous,
  int num_var_unordered,
  int num_var_ordered,
  int num_reg_continuous,
  int num_reg_unordered,
  int num_reg_ordered,
  int *num_categories,
  double *vector_scale_factor,
  double floor_coeff)
{
  double *candidate = vector_scale_factor;
  double *tmp = NULL;
  int invalid = 0;

  if (BANDWIDTH != BW_FIXED)
    return 1;

  if (use_transform) {
    tmp = alloc_vecd(n + 1);
    np_copy_scale_factor(tmp, vector_scale_factor, n);
    bwm_to_constrained(tmp, n);
    candidate = tmp;
  }

  invalid = check_valid_scale_factor_cv(
    KERNEL,
    KERNEL_unordered_liracine,
    BANDWIDTH,
    BANDWIDTH_den_ml,
    REGRESSION_ML,
    num_obs,
    num_var_continuous,
    num_var_unordered,
    num_var_ordered,
    num_reg_continuous,
    num_reg_unordered,
    num_reg_ordered,
    num_categories,
    candidate);

  if ((invalid == 0) &&
      !np_fixed_continuous_floor_ok_cv_with_coeff(
        KERNEL,
        num_obs,
        num_var_continuous,
        num_reg_continuous,
        candidate,
        floor_coeff))
    invalid = 1;

  if (tmp != NULL)
    safe_free(tmp);

  return (invalid == 0);
}

static int np_bw_candidate_is_admissible(
  int n,
  int use_transform,
  int KERNEL,
  int KERNEL_unordered_liracine,
  int BANDWIDTH,
  int BANDWIDTH_den_ml,
  int REGRESSION_ML,
  int num_obs,
  int num_var_continuous,
  int num_var_unordered,
  int num_var_ordered,
  int num_reg_continuous,
  int num_reg_unordered,
  int num_reg_ordered,
  int *num_categories,
  double *vector_scale_factor)
{
  return np_bw_candidate_is_admissible_with_floor(
    n,
    use_transform,
    KERNEL,
    KERNEL_unordered_liracine,
    BANDWIDTH,
    BANDWIDTH_den_ml,
    REGRESSION_ML,
    num_obs,
    num_var_continuous,
    num_var_unordered,
    num_var_ordered,
    num_reg_continuous,
    num_reg_unordered,
    num_reg_ordered,
    num_categories,
    vector_scale_factor,
    bwm_scale_factor_lower_bound);
}

int * ipt_lookup_extern_X;
int * ipt_lookup_extern_Y;
int * ipt_lookup_extern_XY;

/* Quantile - no Y ordered or unordered used, but defined anyways */

double **matrix_Y_continuous_quantile_extern;
double **matrix_Y_unordered_quantile_extern;
double **matrix_Y_ordered_quantile_extern;
double **matrix_X_continuous_quantile_extern;
double **matrix_X_unordered_quantile_extern;
double **matrix_X_ordered_quantile_extern;

double *vector_Y_extern;
double *vector_T_extern;
double *vector_T_resample;
double *vector_Y_eval_extern;
double *vector_Y_null;


int int_ll_extern=0;
int *vector_glp_degree_extern=NULL;
int *vector_glp_gradient_order_extern=NULL;
int int_glp_bernstein_extern=0;
int int_glp_basis_extern=1;
int int_bounded_cvls_quadrature_grid_extern=1;
int int_bounded_cvls_quadrature_points_extern=0;
double double_bounded_cvls_quadrature_extend_factor_extern=1.0;
double double_bounded_cvls_quadrature_ratios_extern[3]={
  0.20, 0.55, 0.25
};

int KERNEL_reg_extern=0;
int KERNEL_reg_unordered_extern=0;
int KERNEL_reg_ordered_extern=0;
int KERNEL_den_extern=0;
int KERNEL_den_unordered_extern=0;
int KERNEL_den_ordered_extern=0;
int BANDWIDTH_reg_extern;
int BANDWIDTH_den_extern;

// cdf algorithm extern
double dbl_memfac_ccdf_extern = 1.0;
double dbl_memfac_dls_extern = 1.0;
int cdfontrain_extern = 0;

/* Statics for dependence metric */

int num_lag_extern;
int int_lag_extern;
int int_iter_extern;

int itmax_extern;
double small_extern;

double *vector_scale_factor_dep_met_bivar_extern;
double *vector_scale_factor_dep_met_univar_extern;
double *vector_scale_factor_dep_met_univar_lag_extern;

double *vector_scale_factor_extern;
double gamma_extern = 0.5;
double y_min_extern;
double y_max_extern;

int imsnum = 0;
int imstot = 0;

KDT * kdt_extern_X = NULL;
KDT * kdt_extern_Y = NULL;
KDT * kdt_extern_XY = NULL;

// to facilitate bandwidth->scale-factor conversions
double * vector_continuous_stddev_extern = NULL;
double nconfac_extern = 0.0;
double ncatfac_extern = 0.0;
static int np_shadow_state_active = 0;

extern int iff;

double np_tgauss2_b = 3.0, np_tgauss2_alpha = 1.030174731161562;
double np_tgauss2_c0 = .004565578246317041;

double np_tgauss2_a0 = 0.2993759121518507, np_tgauss2_a1 = 2.0844504723243343E-5;
double np_tgauss2_a2 = -2.0*0.002351671671248367;

double np_tgauss2_k = 2.90113075268188e-01, np_tgauss2_k2 = 9.17819591566274e-01;
double np_tgauss2_k22 = 1.40866160472795e-01, np_tgauss2_km = 2.23983611906613e-01;

extern double cksup[OP_NCFUN][2];

double timing_extern  = -1.0;

void np_set_tgauss2(double * coefficients){
  np_tgauss2_b = coefficients[TG2_B];
  np_tgauss2_alpha = coefficients[TG2_ALPHA];
  np_tgauss2_c0 = coefficients[TG2_C0];

  np_tgauss2_a0 = coefficients[TG2_A0];
  np_tgauss2_a1 = coefficients[TG2_A1];
  np_tgauss2_a2 = coefficients[TG2_A2];

  np_tgauss2_k = coefficients[TG2_K];
  np_tgauss2_k2 = coefficients[TG2_K2];
  np_tgauss2_k22 = coefficients[TG2_K22];
  np_tgauss2_km = coefficients[TG2_KM];

  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_NORMAL]][0] = -np_tgauss2_b;
  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_NORMAL]][1] = np_tgauss2_b;

  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_CONVOLUTION]][0] = -2.0*np_tgauss2_b;
  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_CONVOLUTION]][1] = 2.0*np_tgauss2_b;

  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_DERIVATIVE]][0] = -np_tgauss2_b;
  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_DERIVATIVE]][1] = np_tgauss2_b;

  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_INTEGRAL]][0] = -np_tgauss2_b;
  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_INTEGRAL]][1] = DBL_MAX;
}

void spinner(int num) {
  (void) num;
}

void np_set_seed(int * num){
  int_RANDOM_SEED = *num;
  iff = 0;
}

SEXP C_np_dim_basis(SEXP basis_code,
                    SEXP kernel,
                    SEXP degree,
                    SEXP segments,
                    SEXP include,
                    SEXP categories)
{
  SEXP degree_i = R_NilValue;
  SEXP segments_i = R_NilValue;
  SEXP include_i = R_NilValue;
  SEXP categories_i = R_NilValue;
  int basis_code_i = asInteger(basis_code);
  int kernel_i = asInteger(kernel);
  int k;
  int ninclude;
  double result = NA_REAL;

  PROTECT(degree_i = coerceVector(degree, INTSXP));
  PROTECT(segments_i = coerceVector(segments, INTSXP));
  PROTECT(include_i = coerceVector(include, INTSXP));
  PROTECT(categories_i = coerceVector(categories, INTSXP));

  k = (int) XLENGTH(degree_i);
  ninclude = (int) XLENGTH(include_i);

  np_dim_basis(&basis_code_i,
               &kernel_i,
               INTEGER(degree_i),
               INTEGER(segments_i),
               &k,
               INTEGER(include_i),
               INTEGER(categories_i),
               &ninclude,
               &result);

  UNPROTECT(4);
  return ScalarReal(result);
}

SEXP C_np_set_seed(SEXP seed)
{
  int num = 0;

  if (XLENGTH(seed) != 1)
    error("C_np_set_seed: seed must have length 1");

  if (TYPEOF(seed) == INTSXP) {
    const int raw = INTEGER(seed)[0];
    if (raw == NA_INTEGER)
      error("C_np_set_seed: seed must be finite");
    if (raw == INT_MIN)
      error("C_np_set_seed: seed must be representable after abs()");
    num = abs(raw);
  } else if (TYPEOF(seed) == REALSXP) {
    const double raw = REAL(seed)[0];
    const double normalized = fabs(raw);
    if (!R_finite(raw))
      error("C_np_set_seed: seed must be finite");
    if (normalized > (double)INT_MAX || normalized != floor(normalized))
      error("C_np_set_seed: seed must be representable as a non-negative integer after abs()");
    num = (int) normalized;
  } else {
    error("C_np_set_seed: seed must be numeric");
  }

  np_set_seed(&num);
  return R_NilValue;
}

SEXP C_np_set_tgauss2(SEXP coefficients)
{
  SEXP coef_r = R_NilValue;

  PROTECT(coef_r = coerceVector(coefficients, REALSXP));
  if (XLENGTH(coef_r) != 10)
    error("C_np_set_tgauss2: coefficients must have length 10");

  np_set_tgauss2(REAL(coef_r));

  UNPROTECT(1);
  return R_NilValue;
}

SEXP C_np_set_local_regression_mode(SEXP active)
{
  const int requested = asLogical(active);
  const int previous = np_mpi_local_regression_mode;

  if (requested == NA_LOGICAL)
    error("C_np_set_local_regression_mode: 'active' must be TRUE or FALSE");

#ifdef MPI2
  if (requested && !np_mpi_local_regression_mode) {
    np_mpi_local_regression_saved_rank = my_rank;
    np_mpi_local_regression_saved_nproc = iNum_Processors;
    np_mpi_local_regression_saved_comm1 = comm[1];
    comm[1] = MPI_COMM_SELF;
    my_rank = 0;
    iNum_Processors = 1;
    np_mpi_local_regression_mode = 1;
  } else if (!requested && np_mpi_local_regression_mode) {
    if (np_mpi_local_regression_saved_comm1 != MPI_COMM_NULL)
      comm[1] = np_mpi_local_regression_saved_comm1;
    my_rank = np_mpi_local_regression_saved_rank;
    iNum_Processors = np_mpi_local_regression_saved_nproc;
    np_mpi_local_regression_saved_comm1 = MPI_COMM_NULL;
    np_mpi_local_regression_mode = 0;
  }
#else
  np_mpi_local_regression_mode = requested ? 1 : 0;
#endif

  return ScalarLogical(previous);
}

SEXP C_np_release_static_buffers(void)
{
  int unused = 0;
  np_release_static_buffers(&unused);
  return R_NilValue;
}

static void np_shadow_reset_state_internal(void)
{
  double **yuno_train = matrix_Y_unordered_train_extern;
  double **yord_train = matrix_Y_ordered_train_extern;
  double **ycon_train = matrix_Y_continuous_train_extern;
  double **yuno_eval = matrix_Y_unordered_eval_extern;
  double **yord_eval = matrix_Y_ordered_eval_extern;
  double **ycon_eval = matrix_Y_continuous_eval_extern;
  double **xuno_train = matrix_X_unordered_train_extern;
  double **xord_train = matrix_X_ordered_train_extern;
  double **xcon_train = matrix_X_continuous_train_extern;
  double **xuno_eval = matrix_X_unordered_eval_extern;
  double **xord_eval = matrix_X_ordered_eval_extern;
  double **xcon_eval = matrix_X_continuous_eval_extern;
  double **xyuno_train = matrix_XY_unordered_train_extern;
  double **xyord_train = matrix_XY_ordered_train_extern;
  double **xycon_train = matrix_XY_continuous_train_extern;
  double **xyuno_eval = matrix_XY_unordered_eval_extern;
  double **xyord_eval = matrix_XY_ordered_eval_extern;
  double **xycon_eval = matrix_XY_continuous_eval_extern;
  int nuno = num_var_unordered_extern;
  int nord = num_var_ordered_extern;
  int ncon = num_var_continuous_extern;
  int runo = num_reg_unordered_extern;
  int rord = num_reg_ordered_extern;
  int rcon = num_reg_continuous_extern;

  np_regression_nomad_shadow_clear_internal();
  np_conditional_density_nomad_shadow_clear_internal();

  if (!np_shadow_state_active)
    return;

  if (kdt_extern_X != NULL) free_kdtree(&kdt_extern_X);
  if (kdt_extern_Y != NULL) free_kdtree(&kdt_extern_Y);
  if (kdt_extern_XY != NULL) free_kdtree(&kdt_extern_XY);

  safe_free(ipt_extern_X);
  safe_free(ipt_extern_Y);
  safe_free(ipt_extern_XY);
  safe_free(ipt_lookup_extern_X);
  safe_free(ipt_lookup_extern_Y);
  safe_free(ipt_lookup_extern_XY);

  if (yuno_eval != NULL && yuno_eval != yuno_train) free_mat(yuno_eval, nuno);
  if (yord_eval != NULL && yord_eval != yord_train) free_mat(yord_eval, nord);
  if (ycon_eval != NULL && ycon_eval != ycon_train) free_mat(ycon_eval, ncon);
  if (xuno_eval != NULL && xuno_eval != xuno_train) free_mat(xuno_eval, runo);
  if (xord_eval != NULL && xord_eval != xord_train) free_mat(xord_eval, rord);
  if (xcon_eval != NULL && xcon_eval != xcon_train) free_mat(xcon_eval, rcon);
  if (xyuno_eval != NULL && xyuno_eval != xyuno_train) free_mat(xyuno_eval, nuno + runo);
  if (xyord_eval != NULL && xyord_eval != xyord_train) free_mat(xyord_eval, nord + rord);
  if (xycon_eval != NULL && xycon_eval != xycon_train) free_mat(xycon_eval, ncon + rcon);

  if (yuno_train != NULL) free_mat(yuno_train, nuno);
  if (yord_train != NULL) free_mat(yord_train, nord);
  if (ycon_train != NULL) free_mat(ycon_train, ncon);
  if (xuno_train != NULL) free_mat(xuno_train, runo);
  if (xord_train != NULL) free_mat(xord_train, rord);
  if (xcon_train != NULL) free_mat(xcon_train, rcon);
  if (xyuno_train != NULL) free_mat(xyuno_train, nuno + runo);
  if (xyord_train != NULL) free_mat(xyord_train, nord + rord);
  if (xycon_train != NULL) free_mat(xycon_train, ncon + rcon);

  if (matrix_categorical_vals_extern != NULL) free_mat(matrix_categorical_vals_extern, nuno + nord + runo + rord);
  if (matrix_categorical_vals_extern_X != NULL) free_mat(matrix_categorical_vals_extern_X, runo + rord);
  if (matrix_categorical_vals_extern_Y != NULL) free_mat(matrix_categorical_vals_extern_Y, nuno + nord);
  if (matrix_categorical_vals_extern_XY != NULL) free_mat(matrix_categorical_vals_extern_XY, nuno + nord + runo + rord);

  safe_free(num_categories_extern);
  safe_free(num_categories_extern_X);
  safe_free(num_categories_extern_Y);
  safe_free(num_categories_extern_XY);
  safe_free(vector_continuous_stddev_extern);

  matrix_Y_unordered_train_extern = NULL;
  matrix_Y_ordered_train_extern = NULL;
  matrix_Y_continuous_train_extern = NULL;
  matrix_Y_unordered_eval_extern = NULL;
  matrix_Y_ordered_eval_extern = NULL;
  matrix_Y_continuous_eval_extern = NULL;
  matrix_X_unordered_train_extern = NULL;
  matrix_X_ordered_train_extern = NULL;
  matrix_X_continuous_train_extern = NULL;
  matrix_X_unordered_eval_extern = NULL;
  matrix_X_ordered_eval_extern = NULL;
  matrix_X_continuous_eval_extern = NULL;
  matrix_XY_unordered_train_extern = NULL;
  matrix_XY_ordered_train_extern = NULL;
  matrix_XY_continuous_train_extern = NULL;
  matrix_XY_unordered_eval_extern = NULL;
  matrix_XY_ordered_eval_extern = NULL;
  matrix_XY_continuous_eval_extern = NULL;
  matrix_categorical_vals_extern = NULL;
  matrix_categorical_vals_extern_X = NULL;
  matrix_categorical_vals_extern_Y = NULL;
  matrix_categorical_vals_extern_XY = NULL;
  num_categories_extern = NULL;
  num_categories_extern_X = NULL;
  num_categories_extern_Y = NULL;
  num_categories_extern_XY = NULL;
  vector_continuous_stddev_extern = NULL;

  num_obs_train_extern = 0;
  num_obs_eval_extern = 0;
  num_var_unordered_extern = 0;
  num_var_ordered_extern = 0;
  num_var_continuous_extern = 0;
  num_reg_unordered_extern = 0;
  num_reg_ordered_extern = 0;
  num_reg_continuous_extern = 0;
  int_ll_extern = LL_LC;
  vector_glp_degree_extern = NULL;
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = 0;
  int_glp_basis_extern = 1;
  int_TREE_X = NP_TREE_FALSE;
  int_TREE_Y = NP_TREE_FALSE;
  int_TREE_XY = NP_TREE_FALSE;
  int_LARGE_SF = 0;
  nconfac_extern = 0.0;
  ncatfac_extern = 0.0;
  KERNEL_reg_extern = 0;
  KERNEL_reg_unordered_extern = 0;
  KERNEL_reg_ordered_extern = 0;
  KERNEL_den_extern = 0;
  KERNEL_den_unordered_extern = 0;
  KERNEL_den_ordered_extern = 0;
  BANDWIDTH_den_extern = 0;
  cdfontrain_extern = 0;

  np_glp_cv_clear_extern();
  np_reg_cv_core_clear_extern();
  np_shadow_state_active = 0;
}

static void np_clear_estimator_extern_aliases(void)
{
  matrix_Y_unordered_train_extern = NULL;
  matrix_Y_ordered_train_extern = NULL;
  matrix_Y_continuous_train_extern = NULL;
  matrix_Y_unordered_eval_extern = NULL;
  matrix_Y_ordered_eval_extern = NULL;
  matrix_Y_continuous_eval_extern = NULL;
  matrix_X_unordered_train_extern = NULL;
  matrix_X_ordered_train_extern = NULL;
  matrix_X_continuous_train_extern = NULL;
  matrix_X_unordered_eval_extern = NULL;
  matrix_X_ordered_eval_extern = NULL;
  matrix_X_continuous_eval_extern = NULL;
  matrix_XY_unordered_train_extern = NULL;
  matrix_XY_ordered_train_extern = NULL;
  matrix_XY_continuous_train_extern = NULL;
  matrix_XY_unordered_eval_extern = NULL;
  matrix_XY_ordered_eval_extern = NULL;
  matrix_XY_continuous_eval_extern = NULL;
  matrix_categorical_vals_extern = NULL;
  matrix_categorical_vals_extern_X = NULL;
  matrix_categorical_vals_extern_Y = NULL;
  matrix_categorical_vals_extern_XY = NULL;
  num_categories_extern = NULL;
  num_categories_extern_X = NULL;
  num_categories_extern_Y = NULL;
  num_categories_extern_XY = NULL;
  vector_continuous_stddev_extern = NULL;
  ipt_extern_X = NULL;
  ipt_extern_Y = NULL;
  ipt_extern_XY = NULL;
  ipt_lookup_extern_X = NULL;
  ipt_lookup_extern_Y = NULL;
  ipt_lookup_extern_XY = NULL;
  kdt_extern_X = NULL;
  kdt_extern_Y = NULL;
  kdt_extern_XY = NULL;
}

SEXP C_np_shadow_reset_state(void)
{
  np_shadow_reset_state_internal();
  return R_NilValue;
}

static void np_regression_nomad_shadow_clear_internal(void)
{
  if (!np_regression_nomad_shadow.active && !np_regression_nomad_shadow.owned)
    return;

  if (kdt_extern_X != NULL)
    free_kdtree(&kdt_extern_X);
  int_TREE_X = NP_TREE_FALSE;

  if (matrix_X_unordered_train_extern != NULL)
    free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  if (matrix_X_ordered_train_extern != NULL)
    free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  if (matrix_X_continuous_train_extern != NULL)
    free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  if (matrix_categorical_vals_extern != NULL)
    free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  safe_free(vector_Y_extern);
  safe_free(num_categories_extern);
  safe_free(vector_continuous_stddev_extern);
  safe_free(np_regression_nomad_shadow.ipt);
  safe_free(np_regression_nomad_shadow.glp_degree);
  safe_free(np_regression_nomad_shadow.ckerlb);
  safe_free(np_regression_nomad_shadow.ckerub);
  safe_free(np_regression_nomad_shadow.vector_scale_factor);
  safe_free(bwm_kernel_unordered_vec);

  matrix_X_unordered_train_extern = NULL;
  matrix_X_ordered_train_extern = NULL;
  matrix_X_continuous_train_extern = NULL;
  matrix_X_unordered_eval_extern = NULL;
  matrix_X_ordered_eval_extern = NULL;
  matrix_X_continuous_eval_extern = NULL;
  matrix_categorical_vals_extern = NULL;
  vector_Y_extern = NULL;
  num_categories_extern = NULL;
  vector_continuous_stddev_extern = NULL;
  bwm_num_categories = NULL;
  bwm_kernel_unordered_vec = NULL;
  bwm_kernel_unordered_len = 0;
  bwmfunc_raw = NULL;
  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  bwm_use_transform = 0;
  bwm_num_reg_continuous = 0;
  bwm_num_reg_unordered = 0;
  bwm_num_reg_ordered = 0;
  bwm_kernel_unordered = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;
  int_cker_bound_extern = 0;
  vector_glp_degree_extern = NULL;
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = 0;
  int_glp_basis_extern = 1;
  int_nn_k_min_extern = 1;
  BANDWIDTH_reg_extern = 0;
  BANDWIDTH_den_extern = 0;
  KERNEL_reg_extern = 0;
  KERNEL_reg_unordered_extern = 0;
  KERNEL_reg_ordered_extern = 0;
  num_obs_train_extern = 0;
  num_reg_continuous_extern = 0;
  num_reg_unordered_extern = 0;
  num_reg_ordered_extern = 0;
  np_reset_y_side_extern();
  np_glp_cv_clear_extern();
  np_reg_cv_core_clear_extern();

  np_regression_nomad_shadow.active = 0;
  np_regression_nomad_shadow.owned = 0;
  np_regression_nomad_shadow.num_var = 0;
  np_regression_nomad_shadow.num_reg_continuous = 0;
  np_regression_nomad_shadow.num_reg_unordered = 0;
  np_regression_nomad_shadow.num_reg_ordered = 0;
}

static int np_regression_nomad_shadow_refresh_degree(const int *degree)
{
  int i;
  int changed = 0;

  if (np_regression_nomad_shadow.num_reg_continuous <= 0)
    return 1;
  if (degree == NULL || np_regression_nomad_shadow.glp_degree == NULL)
    return 0;

  for (i = 0; i < np_regression_nomad_shadow.num_reg_continuous; i++) {
    if (np_regression_nomad_shadow.glp_degree[i] != degree[i]) {
      changed = 1;
      break;
    }
  }

  if (!changed)
    return 1;

  for (i = 0; i < np_regression_nomad_shadow.num_reg_continuous; i++)
    np_regression_nomad_shadow.glp_degree[i] = degree[i];

  vector_glp_degree_extern = np_regression_nomad_shadow.glp_degree;
  np_glp_cv_clear_extern();

  if ((int_ll_extern == LL_LP) &&
      (!np_glp_cv_prepare_extern(int_ll_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 matrix_X_continuous_train_extern)))
    return 0;

  return 1;
}

static int np_regression_nomad_shadow_prepare_internal(double *runo,
                                                       double *rord,
                                                       double *rcon,
                                                       double *y,
                                                       double *mysd,
                                                       int *myopti,
                                                       double *myoptd,
                                                       double *rbw,
                                                       int *penalty_mode,
                                                       double *penalty_mult,
                                                       int *glp_degree,
                                                       int *glp_bernstein,
                                                       int *glp_basis,
                                                       double *ckerlb,
                                                       double *ckerub)
{
  int i, j;
  int num_var;
  int scale_cat;
  int int_use_starting_values;
  int *ipt = NULL;
  int dfc_dir;
  double lbc_dir, c_dir, initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init;
  double lbd_init, hbd_init, d_init;
  double **matrix_y = NULL;
  double *vsfh = NULL;

  np_regression_nomad_shadow_clear_internal();
  np_regression_nomad_shadow.owned = 1;
  num_reg_continuous_extern = myopti[RBW_NCONI];
  num_reg_unordered_extern = myopti[RBW_NUNOI];
  num_reg_ordered_extern = myopti[RBW_NORDI];
  np_reset_y_side_extern();

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[RBW_NOBSI];
  KERNEL_reg_extern = myopti[RBW_CKRNEVI];
  KERNEL_reg_unordered_extern = myopti[RBW_UKRNEVI];
  KERNEL_reg_ordered_extern = myopti[RBW_OKRNEVI];

  np_regression_nomad_shadow.num_var = num_var;
  np_regression_nomad_shadow.num_reg_continuous = num_reg_continuous_extern;
  np_regression_nomad_shadow.num_reg_unordered = num_reg_unordered_extern;
  np_regression_nomad_shadow.num_reg_ordered = num_reg_ordered_extern;
  if (num_reg_continuous_extern > 0) {
    np_regression_nomad_shadow.ckerlb = alloc_vecd(num_reg_continuous_extern);
    np_regression_nomad_shadow.ckerub = alloc_vecd(num_reg_continuous_extern);
    if (np_regression_nomad_shadow.ckerlb == NULL ||
        np_regression_nomad_shadow.ckerub == NULL)
      goto fail;
    for (i = 0; i < num_reg_continuous_extern; i++) {
      np_regression_nomad_shadow.ckerlb[i] = ckerlb[i];
      np_regression_nomad_shadow.ckerub[i] = ckerub[i];
    }
    vector_ckerlb_extern = np_regression_nomad_shadow.ckerlb;
    vector_ckerub_extern = np_regression_nomad_shadow.ckerub;
  } else {
    vector_ckerlb_extern = NULL;
    vector_ckerub_extern = NULL;
  }
  int_cker_bound_extern = np_has_finite_cker_bounds(vector_ckerlb_extern,
                                                    vector_ckerub_extern,
                                                    num_reg_continuous_extern);

  int_use_starting_values = myopti[RBW_USTARTI];
  int_LARGE_SF = myopti[RBW_LSFI];

  BANDWIDTH_reg_extern = myopti[RBW_REGI];
  BANDWIDTH_den_extern = 0;

  int_RESTART_FROM_MIN = RE_MIN_FALSE;
  int_MINIMIZE_IO = IO_MIN_TRUE;

  int_ll_extern = myopti[RBW_LL];
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = *glp_bernstein;
  int_glp_basis_extern = *glp_basis;
  int_nn_k_min_extern =
    ((BANDWIDTH_reg_extern != BW_FIXED) && (num_reg_continuous_extern > 0)) ? 2 : 1;

  if (num_reg_continuous_extern > 0) {
    np_regression_nomad_shadow.glp_degree = alloc_vecu(num_reg_continuous_extern);
    if (np_regression_nomad_shadow.glp_degree == NULL)
      goto fail;
    for (i = 0; i < num_reg_continuous_extern; i++)
      np_regression_nomad_shadow.glp_degree[i] = glp_degree[i];
    vector_glp_degree_extern = np_regression_nomad_shadow.glp_degree;
  } else {
    vector_glp_degree_extern = NULL;
  }

  int_TREE_X = myopti[RBW_DOTREEI];
  scale_cat = myopti[RBW_SCATI];
  bwm_use_transform = myopti[RBW_TBNDI];
  if (BANDWIDTH_reg_extern != BW_FIXED)
    bwm_use_transform = 0;

  nconfac_extern = myoptd[RBW_NCONFD];
  ncatfac_extern = myoptd[RBW_NCATFD];
  bwm_set_scale_factor_lower_bound(myoptd[RBW_SFLOORD]);
  dfc_dir = myopti[RBW_DFC_DIRI];
  lbc_dir = myoptd[RBW_LBC_DIRD];
  c_dir = myoptd[RBW_C_DIRD];
  initc_dir = myoptd[RBW_INITC_DIRD];
  lbd_dir = myoptd[RBW_LBD_DIRD];
  hbd_dir = myoptd[RBW_HBD_DIRD];
  d_dir = myoptd[RBW_D_DIRD];
  initd_dir = myoptd[RBW_INITD_DIRD];
  lbc_init = myoptd[RBW_LBC_INITD];
  hbc_init = myoptd[RBW_HBC_INITD];
  c_init = myoptd[RBW_C_INITD];
  lbd_init = myoptd[RBW_LBD_INITD];
  hbd_init = myoptd[RBW_HBD_INITD];
  d_init = myoptd[RBW_D_INITD];

  imsnum = 0;
  imstot = 0;

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);
  vector_Y_extern = alloc_vecd(num_obs_train_extern);
  num_categories_extern = alloc_vecu(num_reg_unordered_extern + num_reg_ordered_extern);
  matrix_categorical_vals_extern = alloc_matd(num_obs_train_extern,
                                              num_reg_unordered_extern + num_reg_ordered_extern);
  vector_continuous_stddev_extern = alloc_vecd(num_reg_continuous_extern);
  np_regression_nomad_shadow.vector_scale_factor = alloc_vecd(num_var + 1);
  matrix_y = alloc_matd(num_var + 1, num_var + 1);
  vsfh = alloc_vecd(num_var + 1);

  if (((num_reg_unordered_extern > 0) && (matrix_X_unordered_train_extern == NULL)) ||
      ((num_reg_ordered_extern > 0) && (matrix_X_ordered_train_extern == NULL)) ||
      ((num_reg_continuous_extern > 0) && (matrix_X_continuous_train_extern == NULL)) ||
      (vector_Y_extern == NULL) ||
      (((num_reg_unordered_extern + num_reg_ordered_extern) > 0) && (num_categories_extern == NULL)) ||
      (((num_reg_unordered_extern + num_reg_ordered_extern) > 0) && (matrix_categorical_vals_extern == NULL)) ||
      ((num_reg_continuous_extern > 0) && (vector_continuous_stddev_extern == NULL)) ||
      (np_regression_nomad_shadow.vector_scale_factor == NULL) ||
      (matrix_y == NULL) ||
      (vsfh == NULL))
    goto fail;

  for (j = 0; j < num_reg_continuous_extern; j++)
    vector_continuous_stddev_extern[j] = mysd[j];

  if (int_use_starting_values) {
    for (i = 0; i < num_var; i++)
      np_regression_nomad_shadow.vector_scale_factor[i + 1] = rbw[i];
  }

  for (j = 0; j < num_reg_unordered_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_X_unordered_train_extern[j][i] = runo[j * num_obs_train_extern + i];

  for (j = 0; j < num_reg_ordered_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_X_ordered_train_extern[j][i] = rord[j * num_obs_train_extern + i];

  for (j = 0; j < num_reg_continuous_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_X_continuous_train_extern[j][i] = rcon[j * num_obs_train_extern + i];

  for (i = 0; i < num_obs_train_extern; i++)
    vector_Y_extern[i] = y[i];
  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  if (penalty_mode[0] == 1) {
    double pmult = penalty_mult[0];
    double y_mean = 0.0;
    double mse0 = 0.0;
    if (pmult < 1.0) pmult = 1.0;
    for (i = 0; i < num_obs_train_extern; i++)
      y_mean += vector_Y_extern[i];
    y_mean /= (double) num_obs_train_extern;
    for (i = 0; i < num_obs_train_extern; i++) {
      const double dy = vector_Y_extern[i] - y_mean;
      mse0 += dy * dy;
    }
    mse0 /= (double) num_obs_train_extern;
    if (mse0 <= 0.0) mse0 = DBL_MIN;
    if (myopti[RBW_MI] == RBWM_CVAIC) {
      const double denom = 1.0 - 2.0 / ((double) num_obs_train_extern);
      if (denom > 0.0) {
        const double base_aic = log(mse0) + (1.0 / denom);
        bwm_penalty_value = base_aic + log(pmult);
      }
    } else {
      bwm_penalty_value = pmult * mse0;
    }
    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }

  ipt = (int *)malloc((size_t)num_obs_train_extern * sizeof(int));
  if (ipt == NULL)
    goto fail;
  for (i = 0; i < num_obs_train_extern; i++)
    ipt[i] = i;

  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);
  if (int_TREE_X == NP_TREE_TRUE) {
    build_kdtree(matrix_X_continuous_train_extern,
                 num_obs_train_extern,
                 num_reg_continuous_extern,
                 4 * num_reg_continuous_extern,
                 ipt,
                 &kdt_extern_X);

    for (j = 0; j < num_reg_unordered_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_X_unordered_train_extern[j][i] = runo[j * num_obs_train_extern + ipt[i]];

    for (j = 0; j < num_reg_ordered_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_X_ordered_train_extern[j][i] = rord[j * num_obs_train_extern + ipt[i]];

    for (j = 0; j < num_reg_continuous_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_X_continuous_train_extern[j][i] = rcon[j * num_obs_train_extern + ipt[i]];

    for (i = 0; i < num_obs_train_extern; i++)
      vector_Y_extern[i] = y[ipt[i]];
  }

  determine_categorical_vals(num_obs_train_extern,
                             0,
                             0,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             matrix_Y_unordered_train_extern,
                             matrix_Y_ordered_train_extern,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern,
                             matrix_categorical_vals_extern);

  if ((int_ll_extern == LL_LP) &&
      (!np_glp_cv_prepare_extern(int_ll_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 matrix_X_continuous_train_extern)))
    goto fail;

  initialize_nr_vector_scale_factor(BANDWIDTH_reg_extern,
                                    0,
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0,
                                    KERNEL_reg_unordered_extern,
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev_extern,
                                    np_regression_nomad_shadow.vector_scale_factor,
                                    lbc_init, hbc_init, c_init,
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_reg_extern,
                                    0,
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0,
                                    KERNEL_reg_unordered_extern,
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev_extern,
                                    vsfh,
                                    lbc_init, hbc_init, c_init,
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_reg_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           0,
                           0,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED,
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);

  switch (myopti[RBW_MI]) {
    case RBWM_CVAIC:
      bwmfunc_raw = cv_func_regression_categorical_aic_c;
      break;
    case RBWM_CVLS:
      bwmfunc_raw = cv_func_regression_categorical_ls;
      break;
    default:
      error("np.c: invalid bandwidth selection method.");
  }

  bwm_num_reg_continuous = num_reg_continuous_extern;
  bwm_num_reg_unordered = num_reg_unordered_extern;
  bwm_num_reg_ordered = num_reg_ordered_extern;
  bwm_kernel_unordered = KERNEL_reg_unordered_extern;
  bwm_num_categories = num_categories_extern;
  bwm_kernel_unordered_vec = NULL;
  bwm_kernel_unordered_len = 0;
  if (bwm_use_transform) {
    int n = bwm_num_reg_continuous + bwm_num_reg_unordered + bwm_num_reg_ordered;
    bwm_reserve_transform_buf(n + 1);
  }
  bwm_reset_counters();

  if (np_mpi_rank_failure_injected("NP_RMPI_INJECT_REG_SHADOW_PREPARE_FAIL_RANK"))
    goto fail;

  np_regression_nomad_shadow.ipt = ipt;
  ipt = NULL;
  free_mat(matrix_y, num_var + 1);
  safe_free(vsfh);
  np_regression_nomad_shadow.active = 1;
  return 1;

fail:
  if (matrix_y != NULL)
    free_mat(matrix_y, num_var + 1);
  safe_free(vsfh);
  safe_free(ipt);
  np_regression_nomad_shadow_clear_internal();
  return 0;
}

SEXP C_np_regression_nomad_shadow_prepare(SEXP runo,
                                          SEXP rord,
                                          SEXP rcon,
                                          SEXP y,
                                          SEXP mysd,
                                          SEXP myopti,
                                          SEXP myoptd,
                                          SEXP rbw,
                                          SEXP penalty_mode,
                                          SEXP penalty_mult,
                                          SEXP glp_degree,
                                          SEXP glp_bernstein,
                                          SEXP glp_basis,
                                          SEXP ckerlb,
                                          SEXP ckerub)
{
  SEXP runo_r = R_NilValue, rord_r = R_NilValue, rcon_r = R_NilValue, y_r = R_NilValue;
  SEXP mysd_r = R_NilValue, myopti_i = R_NilValue, myoptd_r = R_NilValue;
  SEXP rbw_r = R_NilValue;
  SEXP penalty_mode_i = R_NilValue, penalty_mult_r = R_NilValue;
  SEXP glp_degree_i = R_NilValue, glp_bernstein_i = R_NilValue, glp_basis_i = R_NilValue;
  SEXP ckerlb_r = R_NilValue, ckerub_r = R_NilValue;
  int ok;

  PROTECT(runo_r = coerceVector(runo, REALSXP));
  PROTECT(rord_r = coerceVector(rord, REALSXP));
  PROTECT(rcon_r = coerceVector(rcon, REALSXP));
  PROTECT(y_r = coerceVector(y, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(myoptd_r = coerceVector(myoptd, REALSXP));
  PROTECT(rbw_r = coerceVector(rbw, REALSXP));
  PROTECT(penalty_mode_i = coerceVector(penalty_mode, INTSXP));
  PROTECT(penalty_mult_r = coerceVector(penalty_mult, REALSXP));
  PROTECT(glp_degree_i = coerceVector(glp_degree, INTSXP));
  PROTECT(glp_bernstein_i = coerceVector(glp_bernstein, INTSXP));
  PROTECT(glp_basis_i = coerceVector(glp_basis, INTSXP));
  PROTECT(ckerlb_r = coerceVector(ckerlb, REALSXP));
  PROTECT(ckerub_r = coerceVector(ckerub, REALSXP));

  if (XLENGTH(myoptd_r) <= RBW_SFLOORD) {
    ok = 0;
  } else {
    ok = np_regression_nomad_shadow_prepare_internal(REAL(runo_r),
                                                     REAL(rord_r),
                                                     REAL(rcon_r),
                                                     REAL(y_r),
                                                     REAL(mysd_r),
                                                     INTEGER(myopti_i),
                                                     REAL(myoptd_r),
                                                     REAL(rbw_r),
                                                     INTEGER(penalty_mode_i),
                                                     REAL(penalty_mult_r),
                                                     INTEGER(glp_degree_i),
                                                     INTEGER(glp_bernstein_i),
                                                     INTEGER(glp_basis_i),
                                                     REAL(ckerlb_r),
                                                     REAL(ckerub_r));
  }

  ok = np_mpi_comm1_all_ok(ok);
  if (!ok)
    np_regression_nomad_shadow_clear_internal();

  UNPROTECT(15);

  return ScalarLogical(ok ? 1 : 0);
}

SEXP C_np_regression_nomad_shadow_eval(SEXP rbw, SEXP glp_degree)
{
  SEXP rbw_r = R_NilValue, degree_i = R_NilValue, out = R_NilValue;
  int i;
  double val, fast = 0.0, fast_before = 0.0;
  int degree_ok;
  int bw_ok;

  if (!np_regression_nomad_shadow.active)
    error("resident npreg NOMAD shadow state is not active");

  PROTECT(rbw_r = coerceVector(rbw, REALSXP));
  PROTECT(degree_i = coerceVector(glp_degree, INTSXP));

  if (XLENGTH(rbw_r) != np_regression_nomad_shadow.num_var) {
    UNPROTECT(2);
    error("resident npreg NOMAD shadow received bandwidth vector of unexpected length");
  }
  if (XLENGTH(degree_i) != np_regression_nomad_shadow.num_reg_continuous) {
    UNPROTECT(2);
    error("resident npreg NOMAD shadow received degree vector of unexpected length");
  }

  degree_ok = np_regression_nomad_shadow_refresh_degree(INTEGER(degree_i)) ? 1 : 0;
  if (comm[1] != MPI_COMM_NULL)
    MPI_Allreduce(MPI_IN_PLACE, &degree_ok, 1, MPI_INT, MPI_MIN, comm[1]);

  if (!degree_ok) {
    bwm_eval_count += 1.0;
    bwm_invalid_count += 1.0;
    bwm_maybe_signal_activity();
    val = (bwm_penalty_mode == 1 && R_FINITE(bwm_penalty_value)) ? bwm_penalty_value : DBL_MAX;
    PROTECT(out = allocVector(REALSXP, 2));
    REAL(out)[0] = val;
    REAL(out)[1] = 0.0;
    UNPROTECT(3);
    return out;
  }

  for (i = 0; i < np_regression_nomad_shadow.num_var; i++)
    np_regression_nomad_shadow.vector_scale_factor[i + 1] = REAL(rbw_r)[i];

  bw_ok = np_bw_candidate_is_admissible(
    np_regression_nomad_shadow.num_var,
    bwm_use_transform,
    KERNEL_reg_extern,
    KERNEL_reg_unordered_extern,
    BANDWIDTH_reg_extern,
    BANDWIDTH_reg_extern,
    0,
    num_obs_train_extern,
    0,
    0,
    0,
    num_reg_continuous_extern,
    num_reg_unordered_extern,
    num_reg_ordered_extern,
    num_categories_extern,
    np_regression_nomad_shadow.vector_scale_factor
  ) ? 1 : 0;
  if (comm[1] != MPI_COMM_NULL)
    MPI_Allreduce(MPI_IN_PLACE, &bw_ok, 1, MPI_INT, MPI_MIN, comm[1]);

  if (!bw_ok) {
    bwm_eval_count += 1.0;
    bwm_invalid_count += 1.0;
    bwm_maybe_signal_activity();
    val = (bwm_penalty_mode == 1 && R_FINITE(bwm_penalty_value)) ? bwm_penalty_value : DBL_MAX;
    PROTECT(out = allocVector(REALSXP, 2));
    REAL(out)[0] = val;
    REAL(out)[1] = 0.0;
    UNPROTECT(3);
    return out;
  }

  fast_before = np_fastcv_alllarge_hits_get();
  val = bwmfunc_wrapper(np_regression_nomad_shadow.vector_scale_factor);
  fast = np_fastcv_alllarge_hits_get() - fast_before;
  if (fast < 0.0)
    fast = 0.0;
  if (comm[1] != MPI_COMM_NULL)
    MPI_Bcast(&val, 1, MPI_DOUBLE, 0, comm[1]);
  if (comm[1] != MPI_COMM_NULL)
    MPI_Bcast(&fast, 1, MPI_DOUBLE, 0, comm[1]);

  PROTECT(out = allocVector(REALSXP, 2));
  REAL(out)[0] = val;
  REAL(out)[1] = fast;
  UNPROTECT(3);
  return out;
}

SEXP C_np_regression_nomad_shadow_clear(void)
{
  np_regression_nomad_shadow_clear_internal();
  return R_NilValue;
}

static void np_conditional_density_nomad_shadow_clear_internal(void)
{
  const int num_all_cvar = num_reg_continuous_extern + num_var_continuous_extern;
  const int num_all_uvar = num_reg_unordered_extern + num_var_unordered_extern;
  const int num_all_ovar = num_reg_ordered_extern + num_var_ordered_extern;

  np_bounded_cvls_conditional_quad_context_clear_extern();

  if (!np_conditional_density_nomad_shadow.active && !np_conditional_density_nomad_shadow.owned)
    return;

  if (kdt_extern_X != NULL)
    free_kdtree(&kdt_extern_X);
  if (kdt_extern_Y != NULL)
    free_kdtree(&kdt_extern_Y);
  if (kdt_extern_XY != NULL)
    free_kdtree(&kdt_extern_XY);
  int_TREE_X = NP_TREE_FALSE;
  int_TREE_Y = NP_TREE_FALSE;
  int_TREE_XY = NP_TREE_FALSE;

  if (matrix_Y_unordered_train_extern != NULL)
    free_mat(matrix_Y_unordered_train_extern, num_var_unordered_extern);
  if (matrix_Y_ordered_train_extern != NULL)
    free_mat(matrix_Y_ordered_train_extern, num_var_ordered_extern);
  if (matrix_Y_continuous_train_extern != NULL)
    free_mat(matrix_Y_continuous_train_extern, num_var_continuous_extern);

  if (matrix_X_unordered_train_extern != NULL)
    free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  if (matrix_X_ordered_train_extern != NULL)
    free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  if (matrix_X_continuous_train_extern != NULL)
    free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  if (matrix_XY_continuous_train_extern != NULL)
    free_mat(matrix_XY_continuous_train_extern, num_all_cvar);
  if (matrix_XY_unordered_train_extern != NULL)
    free_mat(matrix_XY_unordered_train_extern, num_all_uvar);
  if (matrix_XY_ordered_train_extern != NULL)
    free_mat(matrix_XY_ordered_train_extern, num_all_ovar);

  if (matrix_categorical_vals_extern != NULL)
    free_mat(matrix_categorical_vals_extern,
             num_var_unordered_extern + num_var_ordered_extern +
             num_reg_unordered_extern + num_reg_ordered_extern);
  if (matrix_categorical_vals_extern_X != NULL)
    free_mat(matrix_categorical_vals_extern_X, num_reg_unordered_extern + num_reg_ordered_extern);
  if (matrix_categorical_vals_extern_Y != NULL)
    free_mat(matrix_categorical_vals_extern_Y, num_var_unordered_extern + num_var_ordered_extern);
  if (matrix_categorical_vals_extern_XY != NULL)
    free_mat(matrix_categorical_vals_extern_XY,
             num_var_unordered_extern + num_var_ordered_extern +
             num_reg_unordered_extern + num_reg_ordered_extern);

  safe_free(num_categories_extern);
  safe_free(num_categories_extern_X);
  safe_free(num_categories_extern_Y);
  safe_free(num_categories_extern_XY);
  safe_free(vector_continuous_stddev_extern);
  safe_free(np_conditional_density_nomad_shadow.glp_degree);
  safe_free(np_conditional_density_nomad_shadow.vector_scale_factor);
  safe_free(np_conditional_density_nomad_shadow.cxkerlb);
  safe_free(np_conditional_density_nomad_shadow.cxkerub);
  safe_free(np_conditional_density_nomad_shadow.cykerlb);
  safe_free(np_conditional_density_nomad_shadow.cykerub);
  safe_free(np_conditional_density_nomad_shadow.cxykerlb);
  safe_free(np_conditional_density_nomad_shadow.cxykerub);
  safe_free(bwm_kernel_unordered_vec);

  safe_free(np_conditional_density_nomad_shadow.ipt_x);
  safe_free(np_conditional_density_nomad_shadow.ipt_lookup_x);
  if (np_conditional_density_nomad_shadow.need_y_side) {
    safe_free(np_conditional_density_nomad_shadow.ipt_y);
    safe_free(np_conditional_density_nomad_shadow.ipt_lookup_y);
  }
  safe_free(np_conditional_density_nomad_shadow.ipt_xy);
  safe_free(np_conditional_density_nomad_shadow.ipt_lookup_xy);

  matrix_Y_unordered_train_extern = NULL;
  matrix_Y_ordered_train_extern = NULL;
  matrix_Y_continuous_train_extern = NULL;
  matrix_Y_unordered_eval_extern = NULL;
  matrix_Y_ordered_eval_extern = NULL;
  matrix_Y_continuous_eval_extern = NULL;
  matrix_X_unordered_train_extern = NULL;
  matrix_X_ordered_train_extern = NULL;
  matrix_X_continuous_train_extern = NULL;
  matrix_X_unordered_eval_extern = NULL;
  matrix_X_ordered_eval_extern = NULL;
  matrix_X_continuous_eval_extern = NULL;
  matrix_XY_continuous_train_extern = NULL;
  matrix_XY_unordered_train_extern = NULL;
  matrix_XY_ordered_train_extern = NULL;
  matrix_XY_continuous_eval_extern = NULL;
  matrix_XY_unordered_eval_extern = NULL;
  matrix_XY_ordered_eval_extern = NULL;
  matrix_categorical_vals_extern = NULL;
  matrix_categorical_vals_extern_X = NULL;
  matrix_categorical_vals_extern_Y = NULL;
  matrix_categorical_vals_extern_XY = NULL;
  ipt_extern_X = NULL;
  ipt_extern_Y = NULL;
  ipt_extern_XY = NULL;
  ipt_lookup_extern_X = NULL;
  ipt_lookup_extern_Y = NULL;
  ipt_lookup_extern_XY = NULL;

  bwmfunc_raw = NULL;
  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  bwm_use_transform = 0;
  bwm_num_reg_continuous = 0;
  bwm_num_reg_unordered = 0;
  bwm_num_reg_ordered = 0;
  bwm_kernel_unordered = 0;
  bwm_num_categories = NULL;
  bwm_kernel_unordered_vec = NULL;
  bwm_kernel_unordered_len = 0;

  int_cxker_bound_extern = 0;
  int_cyker_bound_extern = 0;
  int_cxyker_bound_extern = 0;
  int_cker_bound_extern = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;
  vector_cxkerlb_extern = NULL;
  vector_cxkerub_extern = NULL;
  vector_cykerlb_extern = NULL;
  vector_cykerub_extern = NULL;
  vector_cxykerlb_extern = NULL;
  vector_cxykerub_extern = NULL;
  vector_glp_degree_extern = NULL;
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = 0;
  int_glp_basis_extern = 1;
  int_ll_extern = LL_LC;
  int_nn_k_min_extern = 1;
  BANDWIDTH_den_extern = 0;
  BANDWIDTH_reg_extern = 0;
  KERNEL_reg_extern = 0;
  KERNEL_den_extern = 0;
  KERNEL_reg_unordered_extern = 0;
  KERNEL_den_unordered_extern = 0;
  KERNEL_reg_ordered_extern = 0;
  KERNEL_den_ordered_extern = 0;
  cdfontrain_extern = 0;
  int_WEIGHTS = 0;

  np_glp_cv_clear_extern();
  np_reset_y_side_extern();

  np_conditional_density_nomad_shadow.active = 0;
  np_conditional_density_nomad_shadow.owned = 0;
  np_conditional_density_nomad_shadow.num_all_var = 0;
  np_conditional_density_nomad_shadow.num_reg_continuous = 0;
  np_conditional_density_nomad_shadow.num_var_continuous = 0;
  np_conditional_density_nomad_shadow.num_reg_unordered = 0;
  np_conditional_density_nomad_shadow.num_var_unordered = 0;
  np_conditional_density_nomad_shadow.num_reg_ordered = 0;
  np_conditional_density_nomad_shadow.num_var_ordered = 0;
  np_conditional_density_nomad_shadow.need_y_side = 0;
  np_conditional_density_nomad_shadow.old_cdens = 0;
  np_conditional_density_nomad_shadow.penalty_mode = 0;
  np_conditional_density_nomad_shadow.penalty_multiplier = 0.0;
}

static void np_conditional_density_nomad_shadow_refresh_penalty(void)
{
  int i;

  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;

  if (np_conditional_density_nomad_shadow.penalty_mode != 1)
    return;

  {
    double pmult = np_conditional_density_nomad_shadow.penalty_multiplier;
    double baseline;

    if (pmult < 1.0)
      pmult = 1.0;

    baseline = bwmfunc_raw_current_scale(
      np_conditional_density_nomad_shadow.vector_scale_factor,
      np_conditional_density_nomad_shadow.num_all_var);
    bwm_eval_count += 1.0;
    if (!R_FINITE(baseline) || baseline == DBL_MAX)
      bwm_invalid_count += 1.0;
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      double *tmp = alloc_vecd(np_conditional_density_nomad_shadow.num_all_var + 1);
      if (tmp == NULL)
        return;
      np_copy_scale_factor_for_raw(
        tmp,
        np_conditional_density_nomad_shadow.vector_scale_factor,
        np_conditional_density_nomad_shadow.num_all_var);
      for (i = 1; i <= (bwm_num_reg_continuous); i++)
        tmp[i] *= 2.0;
      for (i = 0; i < bwm_num_reg_unordered; i++) {
        int idx = bwm_num_reg_continuous + 1 + i;
        double maxbw = max_unordered_bw(num_categories_extern[i], KERNEL_den_unordered_extern);
        tmp[idx] = 0.5 * maxbw;
      }
      for (i = 0; i < bwm_num_reg_ordered; i++) {
        int idx = bwm_num_reg_continuous + bwm_num_reg_unordered + 1 + i;
        tmp[idx] = 0.5;
      }
      baseline = bwmfunc_raw(tmp);
      bwm_eval_count += 1.0;
      if (!R_FINITE(baseline) || baseline == DBL_MAX)
        bwm_invalid_count += 1.0;
      safe_free(tmp);
    }

    if (!R_FINITE(baseline) || baseline == DBL_MAX)
      bwm_penalty_value = pmult * 1.0e6;
    else
      bwm_penalty_value = baseline + (fabs(baseline) + 1.0) * pmult;

    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }
}

static int np_conditional_density_nomad_shadow_refresh_degree(const int *degree)
{
  int i;
  int changed = 0;

  if (int_ll_extern != LL_LP || np_conditional_density_nomad_shadow.num_reg_continuous <= 0)
    return 1;
  if (degree == NULL || np_conditional_density_nomad_shadow.glp_degree == NULL)
    return 0;

  for (i = 0; i < np_conditional_density_nomad_shadow.num_reg_continuous; i++) {
    if (np_conditional_density_nomad_shadow.glp_degree[i] != degree[i]) {
      changed = 1;
      break;
    }
  }

  if (!changed)
    return 1;

  for (i = 0; i < np_conditional_density_nomad_shadow.num_reg_continuous; i++)
    np_conditional_density_nomad_shadow.glp_degree[i] = degree[i];

  vector_glp_degree_extern = np_conditional_density_nomad_shadow.glp_degree;
  np_glp_cv_clear_extern();

  if (!np_glp_cv_prepare_extern(int_ll_extern,
                                num_obs_train_extern,
                                num_reg_continuous_extern,
                                matrix_X_continuous_train_extern))
    return 0;

  return 1;
}

static int np_conditional_density_nomad_shadow_prepare_internal(double *c_uno,
                                                                double *c_ord,
                                                                double *c_con,
                                                                double *u_uno,
                                                                double *u_ord,
                                                                double *u_con,
                                                                double *mysd,
                                                                int *myopti,
                                                                double *myoptd,
                                                                double *rbw,
                                                                int *penalty_mode,
                                                                double *penalty_mult,
                                                                int *glp_degree,
                                                                int *glp_bernstein,
                                                                int *glp_basis,
                                                                int *regtype,
                                                                double *cxkerlb,
                                                                double *cxkerub,
                                                                double *cykerlb,
                                                                double *cykerub)
{
  int i, j;
  int num_all_var, num_all_cvar, num_all_uvar, num_all_ovar;
  int int_use_starting_values, ibwmfunc, scale_cat, need_y_side;
  int dfc_dir;
  int *ipt_X = NULL, *ipt_Y = NULL, *ipt_XY = NULL;
  int *ipt_lookup_X = NULL, *ipt_lookup_Y = NULL, *ipt_lookup_XY = NULL;
  double **matrix_y = NULL;
  double *vsfh = NULL;
  double *vector_continuous_stddev = NULL;
  double lbc_dir, c_dir, initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init;
  double lbd_init, hbd_init, d_init;

  np_conditional_density_nomad_shadow_clear_internal();
  np_conditional_density_nomad_shadow.owned = 1;

  num_var_unordered_extern = myopti[CBW_CNUNOI];
  num_var_ordered_extern = myopti[CBW_CNORDI];
  num_var_continuous_extern = myopti[CBW_CNCONI];
  num_reg_unordered_extern = myopti[CBW_UNUNOI];
  num_reg_ordered_extern = myopti[CBW_UNORDI];
  num_reg_continuous_extern = myopti[CBW_UNCONI];
  num_obs_train_extern = myopti[CBW_NOBSI];

  num_all_var = num_reg_continuous_extern + num_var_continuous_extern +
    num_var_unordered_extern + num_var_ordered_extern +
    num_reg_unordered_extern + num_reg_ordered_extern;
  num_all_cvar = num_reg_continuous_extern + num_var_continuous_extern;
  num_all_uvar = num_reg_unordered_extern + num_var_unordered_extern;
  num_all_ovar = num_reg_ordered_extern + num_var_ordered_extern;

  KERNEL_reg_extern = myopti[CBW_CXKRNEVI];
  KERNEL_den_extern = myopti[CBW_CYKRNEVI];
  KERNEL_reg_unordered_extern = myopti[CBW_UXKRNEVI];
  KERNEL_den_unordered_extern = myopti[CBW_UYKRNEVI];
  KERNEL_reg_ordered_extern = myopti[CBW_OXKRNEVI];
  KERNEL_den_ordered_extern = myopti[CBW_OYKRNEVI];

  np_conditional_density_nomad_shadow.num_all_var = num_all_var;
  np_conditional_density_nomad_shadow.num_reg_continuous = num_reg_continuous_extern;
  np_conditional_density_nomad_shadow.num_var_continuous = num_var_continuous_extern;
  np_conditional_density_nomad_shadow.num_reg_unordered = num_reg_unordered_extern;
  np_conditional_density_nomad_shadow.num_var_unordered = num_var_unordered_extern;
  np_conditional_density_nomad_shadow.num_reg_ordered = num_reg_ordered_extern;
  np_conditional_density_nomad_shadow.num_var_ordered = num_var_ordered_extern;

  if (num_reg_continuous_extern > 0) {
    np_conditional_density_nomad_shadow.cxkerlb = alloc_vecd(num_reg_continuous_extern);
    np_conditional_density_nomad_shadow.cxkerub = alloc_vecd(num_reg_continuous_extern);
    if (np_conditional_density_nomad_shadow.cxkerlb == NULL ||
        np_conditional_density_nomad_shadow.cxkerub == NULL)
      goto fail;
    for (i = 0; i < num_reg_continuous_extern; i++) {
      np_conditional_density_nomad_shadow.cxkerlb[i] = cxkerlb[i];
      np_conditional_density_nomad_shadow.cxkerub[i] = cxkerub[i];
    }
    vector_cxkerlb_extern = np_conditional_density_nomad_shadow.cxkerlb;
    vector_cxkerub_extern = np_conditional_density_nomad_shadow.cxkerub;
  } else {
    vector_cxkerlb_extern = NULL;
    vector_cxkerub_extern = NULL;
  }
  int_cxker_bound_extern = np_has_finite_cker_bounds(vector_cxkerlb_extern,
                                                     vector_cxkerub_extern,
                                                     num_reg_continuous_extern);

  if (num_var_continuous_extern > 0) {
    np_conditional_density_nomad_shadow.cykerlb = alloc_vecd(num_var_continuous_extern);
    np_conditional_density_nomad_shadow.cykerub = alloc_vecd(num_var_continuous_extern);
    if (np_conditional_density_nomad_shadow.cykerlb == NULL ||
        np_conditional_density_nomad_shadow.cykerub == NULL)
      goto fail;
    for (i = 0; i < num_var_continuous_extern; i++) {
      np_conditional_density_nomad_shadow.cykerlb[i] = cykerlb[i];
      np_conditional_density_nomad_shadow.cykerub[i] = cykerub[i];
    }
    vector_cykerlb_extern = np_conditional_density_nomad_shadow.cykerlb;
    vector_cykerub_extern = np_conditional_density_nomad_shadow.cykerub;
  } else {
    vector_cykerlb_extern = NULL;
    vector_cykerub_extern = NULL;
  }
  int_cyker_bound_extern = np_has_finite_cker_bounds(vector_cykerlb_extern,
                                                     vector_cykerub_extern,
                                                     num_var_continuous_extern);

  if (num_all_cvar > 0) {
    np_conditional_density_nomad_shadow.cxykerlb = alloc_vecd(num_all_cvar);
    np_conditional_density_nomad_shadow.cxykerub = alloc_vecd(num_all_cvar);
    if (np_conditional_density_nomad_shadow.cxykerlb == NULL ||
        np_conditional_density_nomad_shadow.cxykerub == NULL)
      goto fail;
    for (i = 0; i < num_reg_continuous_extern; i++) {
      np_conditional_density_nomad_shadow.cxykerlb[i] = vector_cxkerlb_extern[i];
      np_conditional_density_nomad_shadow.cxykerub[i] = vector_cxkerub_extern[i];
    }
    for (i = 0; i < num_var_continuous_extern; i++) {
      np_conditional_density_nomad_shadow.cxykerlb[num_reg_continuous_extern + i] = vector_cykerlb_extern[i];
      np_conditional_density_nomad_shadow.cxykerub[num_reg_continuous_extern + i] = vector_cykerub_extern[i];
    }
    vector_cxykerlb_extern = np_conditional_density_nomad_shadow.cxykerlb;
    vector_cxykerub_extern = np_conditional_density_nomad_shadow.cxykerub;
    int_cxyker_bound_extern = np_has_finite_cker_bounds(vector_cxykerlb_extern,
                                                        vector_cxykerub_extern,
                                                        num_all_cvar);
  } else {
    vector_cxykerlb_extern = NULL;
    vector_cxykerub_extern = NULL;
    int_cxyker_bound_extern = 0;
  }

  int_use_starting_values = myopti[CBW_USTARTI];
  int_LARGE_SF = myopti[CBW_LSFI];
  BANDWIDTH_den_extern = myopti[CBW_DENI];
  int_RESTART_FROM_MIN = RE_MIN_FALSE;
  int_MINIMIZE_IO = IO_MIN_TRUE;
  int_WEIGHTS = 0;
  np_conditional_density_nomad_shadow.old_cdens = myopti[CBW_OLDI];

  ibwmfunc = myopti[CBW_MI];
  int_ll_extern = ((ibwmfunc == CBWM_CVML) || (ibwmfunc == CBWM_CVLS)) ? *regtype : LL_LC;
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = ((ibwmfunc == CBWM_CVML) || (ibwmfunc == CBWM_CVLS)) && (int_ll_extern == LL_LP) ? *glp_bernstein : 0;
  int_glp_basis_extern = ((ibwmfunc == CBWM_CVML) || (ibwmfunc == CBWM_CVLS)) && (int_ll_extern == LL_LP) ? *glp_basis : 1;
  int_bounded_cvls_quadrature_grid_extern = myopti[CBW_CVLS_QUAD_GRIDI];
  if ((int_bounded_cvls_quadrature_grid_extern < 0) ||
      (int_bounded_cvls_quadrature_grid_extern > 2))
    goto fail;
  int_bounded_cvls_quadrature_points_extern = myopti[CBW_CVLS_QUAD_POINTSI];
  if (int_bounded_cvls_quadrature_points_extern < 2)
    int_bounded_cvls_quadrature_points_extern = 0;
  int_nn_k_min_extern = 1;

  if ((int_ll_extern == LL_LP) && (num_reg_continuous_extern > 0)) {
    np_conditional_density_nomad_shadow.glp_degree = alloc_vecu(num_reg_continuous_extern);
    if (np_conditional_density_nomad_shadow.glp_degree == NULL)
      goto fail;
    for (i = 0; i < num_reg_continuous_extern; i++)
      np_conditional_density_nomad_shadow.glp_degree[i] = glp_degree[i];
    vector_glp_degree_extern = np_conditional_density_nomad_shadow.glp_degree;
  } else {
    vector_glp_degree_extern = NULL;
  }

  need_y_side = (ibwmfunc == CBWM_CVLS) || ((ibwmfunc == CBWM_CVML) && (int_ll_extern == LL_LP));
  np_conditional_density_nomad_shadow.need_y_side = need_y_side;
  np_conditional_density_nomad_shadow.penalty_mode = penalty_mode[0];
  np_conditional_density_nomad_shadow.penalty_multiplier = penalty_mult[0];

  int_TREE_XY = int_TREE_Y = int_TREE_X = myopti[CBW_TREEI];
  if(int_ll_extern == LL_LP){
    int_TREE_Y = NP_TREE_FALSE;
    int_TREE_XY = NP_TREE_FALSE;
  }
  scale_cat = myopti[CBW_SCATI];
  bwm_use_transform = myopti[CBW_TBNDI];
  if (BANDWIDTH_den_extern != BW_FIXED)
    bwm_use_transform = 0;
  if (bwm_use_transform)
    bwm_reserve_transform_buf(num_all_var + 1);

  nconfac_extern = myoptd[CBW_NCONFD];
  ncatfac_extern = myoptd[CBW_NCATFD];
  dbl_memfac_ccdf_extern = myoptd[CBW_MEMFACD];
  bwm_set_scale_factor_lower_bound(myoptd[CBW_SFLOORD]);
  double_bounded_cvls_quadrature_extend_factor_extern = myoptd[CBW_QUAD_EXTD];
  if (!R_FINITE(double_bounded_cvls_quadrature_extend_factor_extern) ||
      double_bounded_cvls_quadrature_extend_factor_extern <= 0.0)
    goto fail;
  dfc_dir = myopti[CBW_DFC_DIRI];
  lbc_dir = myoptd[CBW_LBC_DIRD];
  c_dir = myoptd[CBW_C_DIRD];
  initc_dir = myoptd[CBW_INITC_DIRD];
  lbd_dir = myoptd[CBW_LBD_DIRD];
  hbd_dir = myoptd[CBW_HBD_DIRD];
  d_dir = myoptd[CBW_D_DIRD];
  initd_dir = myoptd[CBW_INITD_DIRD];
  lbc_init = myoptd[CBW_LBC_INITD];
  hbc_init = myoptd[CBW_HBC_INITD];
  c_init = myoptd[CBW_C_INITD];
  lbd_init = myoptd[CBW_LBD_INITD];
  hbd_init = myoptd[CBW_HBD_INITD];
  d_init = myoptd[CBW_D_INITD];

  imsnum = 0;
  imstot = 0;

  matrix_Y_unordered_train_extern = alloc_matd(num_obs_train_extern, num_var_unordered_extern);
  matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);
  matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);
  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);
  num_categories_extern = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                     num_reg_unordered_extern + num_reg_ordered_extern);
  num_categories_extern_X = alloc_vecu(num_reg_unordered_extern + num_reg_ordered_extern);
  num_categories_extern_Y = need_y_side ? alloc_vecu(num_var_unordered_extern + num_var_ordered_extern) : NULL;
  num_categories_extern_XY = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                        num_reg_unordered_extern + num_reg_ordered_extern);
  matrix_y = alloc_matd(num_all_var + 1, num_all_var + 1);
  np_conditional_density_nomad_shadow.vector_scale_factor = alloc_vecd(num_all_var + 1);
  vsfh = alloc_vecd(num_all_var + 1);
  matrix_categorical_vals_extern =
    alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern +
               num_reg_unordered_extern + num_reg_ordered_extern);
  matrix_categorical_vals_extern_X =
    alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);
  matrix_categorical_vals_extern_Y =
    need_y_side ? alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern) : NULL;
  matrix_categorical_vals_extern_XY =
    alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern +
               num_reg_unordered_extern + num_reg_ordered_extern);
  matrix_XY_continuous_train_extern = alloc_matd(num_obs_train_extern, num_all_cvar);
  matrix_XY_unordered_train_extern = alloc_matd(num_obs_train_extern, num_all_uvar);
  matrix_XY_ordered_train_extern = alloc_matd(num_obs_train_extern, num_all_ovar);
  vector_continuous_stddev = alloc_vecd(num_all_cvar);
  vector_continuous_stddev_extern = vector_continuous_stddev;

  if ((((num_var_unordered_extern > 0) || (num_reg_unordered_extern > 0)) && (matrix_Y_unordered_train_extern == NULL && num_var_unordered_extern > 0)) ||
      (((num_var_ordered_extern > 0) || (num_reg_ordered_extern > 0)) && (matrix_Y_ordered_train_extern == NULL && num_var_ordered_extern > 0)) ||
      ((num_var_continuous_extern > 0) && (matrix_Y_continuous_train_extern == NULL)) ||
      ((num_reg_unordered_extern > 0) && (matrix_X_unordered_train_extern == NULL)) ||
      ((num_reg_ordered_extern > 0) && (matrix_X_ordered_train_extern == NULL)) ||
      ((num_reg_continuous_extern > 0) && (matrix_X_continuous_train_extern == NULL)) ||
      (((num_var_unordered_extern + num_var_ordered_extern + num_reg_unordered_extern + num_reg_ordered_extern) > 0) && (num_categories_extern == NULL)) ||
      (((num_reg_unordered_extern + num_reg_ordered_extern) > 0) && (num_categories_extern_X == NULL)) ||
      (need_y_side && ((num_var_unordered_extern + num_var_ordered_extern) > 0) && (num_categories_extern_Y == NULL)) ||
      (((num_var_unordered_extern + num_var_ordered_extern + num_reg_unordered_extern + num_reg_ordered_extern) > 0) && (num_categories_extern_XY == NULL)) ||
      (matrix_y == NULL) ||
      (np_conditional_density_nomad_shadow.vector_scale_factor == NULL) ||
      (vsfh == NULL) ||
      (((num_var_unordered_extern + num_var_ordered_extern + num_reg_unordered_extern + num_reg_ordered_extern) > 0) && (matrix_categorical_vals_extern == NULL)) ||
      (((num_reg_unordered_extern + num_reg_ordered_extern) > 0) && (matrix_categorical_vals_extern_X == NULL)) ||
      (need_y_side && ((num_var_unordered_extern + num_var_ordered_extern) > 0) && (matrix_categorical_vals_extern_Y == NULL)) ||
      (((num_var_unordered_extern + num_var_ordered_extern + num_reg_unordered_extern + num_reg_ordered_extern) > 0) && (matrix_categorical_vals_extern_XY == NULL)) ||
      ((num_all_cvar > 0) && (matrix_XY_continuous_train_extern == NULL)) ||
      ((num_all_uvar > 0) && (matrix_XY_unordered_train_extern == NULL)) ||
      ((num_all_ovar > 0) && (matrix_XY_ordered_train_extern == NULL)) ||
      ((num_all_cvar > 0) && (vector_continuous_stddev == NULL)))
    goto fail;

  if (int_use_starting_values) {
    for (i = 0; i < num_all_var; i++)
      np_conditional_density_nomad_shadow.vector_scale_factor[i + 1] = rbw[i];
  }

  for (j = 0; j < num_var_unordered_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_Y_unordered_train_extern[j][i] = c_uno[j * num_obs_train_extern + i];
  for (j = 0; j < num_var_ordered_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_Y_ordered_train_extern[j][i] = c_ord[j * num_obs_train_extern + i];
  for (j = 0; j < num_var_continuous_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_Y_continuous_train_extern[j][i] = c_con[j * num_obs_train_extern + i];
  for (j = 0; j < num_reg_unordered_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_X_unordered_train_extern[j][i] = u_uno[j * num_obs_train_extern + i];
  for (j = 0; j < num_reg_ordered_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_X_ordered_train_extern[j][i] = u_ord[j * num_obs_train_extern + i];
  for (j = 0; j < num_reg_continuous_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_X_continuous_train_extern[j][i] = u_con[j * num_obs_train_extern + i];

  ipt_X = (int *)malloc((size_t)num_obs_train_extern * sizeof(int));
  ipt_lookup_X = (int *)malloc((size_t)num_obs_train_extern * sizeof(int));
  if (ipt_X == NULL || ipt_lookup_X == NULL)
    goto fail;
  for (i = 0; i < num_obs_train_extern; i++) {
    ipt_X[i] = i;
    ipt_lookup_X[i] = i;
  }
  ipt_extern_X = ipt_X;
  ipt_lookup_extern_X = ipt_lookup_X;

  if (need_y_side) {
    ipt_Y = (int *)malloc((size_t)num_obs_train_extern * sizeof(int));
    ipt_lookup_Y = (int *)malloc((size_t)num_obs_train_extern * sizeof(int));
    if (ipt_Y == NULL || ipt_lookup_Y == NULL)
      goto fail;
    for (i = 0; i < num_obs_train_extern; i++) {
      ipt_Y[i] = i;
      ipt_lookup_Y[i] = i;
    }
  } else {
    ipt_Y = ipt_X;
    ipt_lookup_Y = ipt_lookup_X;
  }
  ipt_extern_Y = ipt_Y;
  ipt_lookup_extern_Y = ipt_lookup_Y;

  ipt_XY = (int *)malloc((size_t)num_obs_train_extern * sizeof(int));
  ipt_lookup_XY = (int *)malloc((size_t)num_obs_train_extern * sizeof(int));
  if (ipt_XY == NULL || ipt_lookup_XY == NULL)
    goto fail;
  for (i = 0; i < num_obs_train_extern; i++) {
    ipt_XY[i] = i;
    ipt_lookup_XY[i] = i;
  }
  ipt_extern_XY = ipt_XY;
  ipt_lookup_extern_XY = ipt_lookup_XY;

  int_TREE_XY = int_TREE_XY && ((num_all_cvar != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);
  int_TREE_Y = int_TREE_Y && need_y_side && ((num_var_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);

  if (int_TREE_X == NP_TREE_TRUE) {
    build_kdtree(matrix_X_continuous_train_extern,
                 num_obs_train_extern,
                 num_reg_continuous_extern,
                 4 * num_reg_continuous_extern,
                 ipt_X,
                 &kdt_extern_X);

    for (j = 0; j < num_reg_unordered_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_X_unordered_train_extern[j][i] = u_uno[j * num_obs_train_extern + ipt_X[i]];
    for (j = 0; j < num_reg_ordered_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_X_ordered_train_extern[j][i] = u_ord[j * num_obs_train_extern + ipt_X[i]];
    for (j = 0; j < num_reg_continuous_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_X_continuous_train_extern[j][i] = u_con[j * num_obs_train_extern + ipt_X[i]];
    for (i = 0; i < num_obs_train_extern; i++)
      ipt_lookup_X[ipt_X[i]] = i;
  }

  if (int_TREE_Y == NP_TREE_TRUE) {
    build_kdtree(matrix_Y_continuous_train_extern,
                 num_obs_train_extern,
                 num_var_continuous_extern,
                 4 * num_var_continuous_extern,
                 ipt_Y,
                 &kdt_extern_Y);

    for (j = 0; j < num_var_unordered_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_Y_unordered_train_extern[j][i] = c_uno[j * num_obs_train_extern + ipt_Y[i]];
    for (j = 0; j < num_var_ordered_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_Y_ordered_train_extern[j][i] = c_ord[j * num_obs_train_extern + ipt_Y[i]];
    for (j = 0; j < num_var_continuous_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_Y_continuous_train_extern[j][i] = c_con[j * num_obs_train_extern + ipt_Y[i]];
    for (i = 0; i < num_obs_train_extern; i++)
      ipt_lookup_Y[ipt_Y[i]] = i;
  }

  for (j = 0; j < num_reg_unordered_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_XY_unordered_train_extern[j][i] = u_uno[j * num_obs_train_extern + i];
  for (j = num_reg_unordered_extern; j < num_all_uvar; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_XY_unordered_train_extern[j][i] = c_uno[(j - num_reg_unordered_extern) * num_obs_train_extern + i];
  for (j = 0; j < num_reg_ordered_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_XY_ordered_train_extern[j][i] = u_ord[j * num_obs_train_extern + i];
  for (j = num_reg_ordered_extern; j < num_all_ovar; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_XY_ordered_train_extern[j][i] = c_ord[(j - num_reg_ordered_extern) * num_obs_train_extern + i];
  for (j = 0; j < num_reg_continuous_extern; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_XY_continuous_train_extern[j][i] = u_con[j * num_obs_train_extern + i];
  for (j = num_reg_continuous_extern; j < num_all_cvar; j++)
    for (i = 0; i < num_obs_train_extern; i++)
      matrix_XY_continuous_train_extern[j][i] = c_con[(j - num_reg_continuous_extern) * num_obs_train_extern + i];

  if (int_TREE_XY == NP_TREE_TRUE) {
    build_kdtree(matrix_XY_continuous_train_extern,
                 num_obs_train_extern,
                 num_all_cvar,
                 4 * num_all_cvar,
                 ipt_XY,
                 &kdt_extern_XY);

    for (j = 0; j < num_reg_unordered_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i] = u_uno[j * num_obs_train_extern + ipt_XY[i]];
    for (j = num_reg_unordered_extern; j < num_all_uvar; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i] = c_uno[(j - num_reg_unordered_extern) * num_obs_train_extern + ipt_XY[i]];
    for (j = 0; j < num_reg_ordered_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i] = u_ord[j * num_obs_train_extern + ipt_XY[i]];
    for (j = num_reg_ordered_extern; j < num_all_ovar; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i] = c_ord[(j - num_reg_ordered_extern) * num_obs_train_extern + ipt_XY[i]];
    for (j = 0; j < num_reg_continuous_extern; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i] = u_con[j * num_obs_train_extern + ipt_XY[i]];
    for (j = num_reg_continuous_extern; j < num_all_cvar; j++)
      for (i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i] = c_con[(j - num_reg_continuous_extern) * num_obs_train_extern + ipt_XY[i]];
    for (i = 0; i < num_obs_train_extern; i++)
      ipt_lookup_XY[ipt_XY[i]] = i;
  }

  determine_categorical_vals(num_obs_train_extern,
                             num_var_unordered_extern,
                             num_var_ordered_extern,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             matrix_Y_unordered_train_extern,
                             matrix_Y_ordered_train_extern,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern,
                             matrix_categorical_vals_extern);

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern, num_var_ordered_extern, num_var_continuous_extern,
                        num_reg_unordered_extern, num_reg_ordered_extern, num_reg_continuous_extern,
                        np_conditional_density_nomad_shadow.vector_scale_factor + 1,
                        num_categories_extern,
                        matrix_categorical_vals_extern,
                        NULL, NULL, NULL,
                        num_categories_extern_X,
                        need_y_side ? num_categories_extern_Y : NULL,
                        num_categories_extern_XY,
                        matrix_categorical_vals_extern_X,
                        need_y_side ? matrix_categorical_vals_extern_Y : NULL,
                        matrix_categorical_vals_extern_XY);

  for (j = 0; j < num_all_cvar; j++)
    vector_continuous_stddev_extern[j] = mysd[j];

  if ((int_ll_extern == LL_LP) &&
      (!np_glp_cv_prepare_extern(int_ll_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 matrix_X_continuous_train_extern)))
    goto fail;

  if ((ibwmfunc == CBWM_CVLS) && (int_ll_extern == LL_LP) &&
      (np_bounded_cvls_conditional_quad_context_prepare_extern() != 0))
    goto fail;

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    num_var_continuous_extern,
                                    num_var_unordered_extern,
                                    num_var_ordered_extern,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    KERNEL_den_unordered_extern,
                                    KERNEL_reg_unordered_extern,
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0, 0.2),
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev_extern,
                                    np_conditional_density_nomad_shadow.vector_scale_factor,
                                    lbc_init, hbc_init, c_init,
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    num_var_continuous_extern,
                                    num_var_unordered_extern,
                                    num_var_ordered_extern,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    KERNEL_den_unordered_extern,
                                    KERNEL_reg_unordered_extern,
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0, 0.2),
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev_extern,
                                    vsfh,
                                    lbc_init, hbc_init, c_init,
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           num_var_continuous_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED,
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);

  if (np_conditional_density_nomad_shadow.old_cdens) {
    switch (ibwmfunc) {
      case CBWM_CVML:
        bwmfunc_raw = cv_func_con_density_categorical_ml;
        break;
      case CBWM_CVLS:
        bwmfunc_raw = cv_func_con_density_categorical_ls;
        break;
      case CBWM_NPLS:
        bwmfunc_raw = np_cv_func_con_density_categorical_ls;
        break;
      default:
        error("np.c: invalid bandwidth selection method.");
    }
  } else {
    switch (ibwmfunc) {
      case CBWM_CVML:
        bwmfunc_raw = np_cv_func_con_density_categorical_ml;
        break;
      case CBWM_CVLS:
        bwmfunc_raw = np_cv_func_con_density_categorical_ls_npksum;
        break;
      default:
        error("np.c: invalid bandwidth selection method.");
    }
  }

  bwm_num_reg_continuous = num_all_cvar;
  bwm_num_reg_unordered = num_all_uvar;
  bwm_num_reg_ordered = num_all_ovar;
  bwm_kernel_unordered = KERNEL_den_unordered_extern;
  bwm_kernel_unordered_len = bwm_num_reg_unordered;
  if (bwm_kernel_unordered_len > 0) {
    bwm_kernel_unordered_vec = alloc_vecu(bwm_kernel_unordered_len);
    if (bwm_kernel_unordered_vec == NULL)
      goto fail;
    for (i = 0; i < num_var_unordered_extern; i++)
      bwm_kernel_unordered_vec[i] = KERNEL_den_unordered_extern;
    for (i = 0; i < num_reg_unordered_extern; i++)
      bwm_kernel_unordered_vec[num_var_unordered_extern + i] = KERNEL_reg_unordered_extern;
  } else {
    bwm_kernel_unordered_vec = NULL;
  }
  bwm_num_categories = num_categories_extern;
  bwm_reset_counters();

  np_conditional_density_nomad_shadow_refresh_penalty();

  if (np_mpi_rank_failure_injected("NP_RMPI_INJECT_CDEN_SHADOW_PREPARE_FAIL_RANK"))
    goto fail;

  np_conditional_density_nomad_shadow.ipt_x = ipt_X;
  np_conditional_density_nomad_shadow.ipt_y = need_y_side ? ipt_Y : NULL;
  np_conditional_density_nomad_shadow.ipt_xy = ipt_XY;
  np_conditional_density_nomad_shadow.ipt_lookup_x = ipt_lookup_X;
  np_conditional_density_nomad_shadow.ipt_lookup_y = need_y_side ? ipt_lookup_Y : NULL;
  np_conditional_density_nomad_shadow.ipt_lookup_xy = ipt_lookup_XY;
  ipt_X = ipt_Y = ipt_XY = NULL;
  ipt_lookup_X = ipt_lookup_Y = ipt_lookup_XY = NULL;

  free_mat(matrix_y, num_all_var + 1);
  safe_free(vsfh);
  np_conditional_density_nomad_shadow.active = 1;
  return 1;

fail:
  if (matrix_y != NULL)
    free_mat(matrix_y, num_all_var + 1);
  safe_free(vsfh);
  safe_free(ipt_X);
  if (ipt_Y != ipt_X)
    safe_free(ipt_Y);
  safe_free(ipt_XY);
  safe_free(ipt_lookup_X);
  if (ipt_lookup_Y != ipt_lookup_X)
    safe_free(ipt_lookup_Y);
  safe_free(ipt_lookup_XY);
  np_conditional_density_nomad_shadow_clear_internal();
  return 0;
}

SEXP C_np_density_conditional_nomad_shadow_prepare(SEXP c_uno,
                                                   SEXP c_ord,
                                                   SEXP c_con,
                                                   SEXP u_uno,
                                                   SEXP u_ord,
                                                   SEXP u_con,
                                                   SEXP mysd,
                                                   SEXP myopti,
                                                   SEXP myoptd,
                                                   SEXP rbw,
                                                   SEXP penalty_mode,
                                                   SEXP penalty_mult,
                                                   SEXP glp_degree,
                                                   SEXP glp_bernstein,
                                                   SEXP glp_basis,
                                                   SEXP regtype,
                                                   SEXP cxkerlb,
                                                   SEXP cxkerub,
                                                   SEXP cykerlb,
                                                   SEXP cykerub)
{
  SEXP c_uno_r = R_NilValue, c_ord_r = R_NilValue, c_con_r = R_NilValue;
  SEXP u_uno_r = R_NilValue, u_ord_r = R_NilValue, u_con_r = R_NilValue;
  SEXP mysd_r = R_NilValue, myopti_i = R_NilValue, myoptd_r = R_NilValue, rbw_r = R_NilValue;
  SEXP penalty_mode_i = R_NilValue, penalty_mult_r = R_NilValue;
  SEXP glp_degree_i = R_NilValue, glp_bernstein_i = R_NilValue, glp_basis_i = R_NilValue, regtype_i = R_NilValue;
  SEXP cxkerlb_r = R_NilValue, cxkerub_r = R_NilValue, cykerlb_r = R_NilValue, cykerub_r = R_NilValue;
  int ok;

  PROTECT(c_uno_r = coerceVector(c_uno, REALSXP));
  PROTECT(c_ord_r = coerceVector(c_ord, REALSXP));
  PROTECT(c_con_r = coerceVector(c_con, REALSXP));
  PROTECT(u_uno_r = coerceVector(u_uno, REALSXP));
  PROTECT(u_ord_r = coerceVector(u_ord, REALSXP));
  PROTECT(u_con_r = coerceVector(u_con, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(myoptd_r = coerceVector(myoptd, REALSXP));
  PROTECT(rbw_r = coerceVector(rbw, REALSXP));
  PROTECT(penalty_mode_i = coerceVector(penalty_mode, INTSXP));
  PROTECT(penalty_mult_r = coerceVector(penalty_mult, REALSXP));
  PROTECT(glp_degree_i = coerceVector(glp_degree, INTSXP));
  PROTECT(glp_bernstein_i = coerceVector(glp_bernstein, INTSXP));
  PROTECT(glp_basis_i = coerceVector(glp_basis, INTSXP));
  PROTECT(regtype_i = coerceVector(regtype, INTSXP));
  PROTECT(cxkerlb_r = coerceVector(cxkerlb, REALSXP));
  PROTECT(cxkerub_r = coerceVector(cxkerub, REALSXP));
  PROTECT(cykerlb_r = coerceVector(cykerlb, REALSXP));
  PROTECT(cykerub_r = coerceVector(cykerub, REALSXP));

  if (XLENGTH(myopti_i) <= CBW_CVLS_QUAD_POINTSI ||
      XLENGTH(myoptd_r) <= CBW_QUAD_EXTD) {
    ok = 0;
  } else {
    ok = np_conditional_density_nomad_shadow_prepare_internal(REAL(c_uno_r),
                                                              REAL(c_ord_r),
                                                              REAL(c_con_r),
                                                              REAL(u_uno_r),
                                                              REAL(u_ord_r),
                                                              REAL(u_con_r),
                                                              REAL(mysd_r),
                                                              INTEGER(myopti_i),
                                                              REAL(myoptd_r),
                                                              REAL(rbw_r),
                                                              INTEGER(penalty_mode_i),
                                                              REAL(penalty_mult_r),
                                                              INTEGER(glp_degree_i),
                                                              INTEGER(glp_bernstein_i),
                                                              INTEGER(glp_basis_i),
                                                              INTEGER(regtype_i),
                                                              REAL(cxkerlb_r),
                                                              REAL(cxkerub_r),
                                                              REAL(cykerlb_r),
                                                              REAL(cykerub_r));
  }

  ok = np_mpi_comm1_all_ok(ok);
  if (!ok)
    np_conditional_density_nomad_shadow_clear_internal();

  UNPROTECT(20);

  return ScalarLogical(ok ? 1 : 0);
}

SEXP C_np_density_conditional_nomad_shadow_eval(SEXP rbw, SEXP glp_degree)
{
  SEXP rbw_r = R_NilValue, degree_i = R_NilValue, out = R_NilValue;
  int i;
  double val, fast = 0.0, fast_before = 0.0, evals = 0.0, eval_before = 0.0;

  if (!np_conditional_density_nomad_shadow.active)
    error("resident npcdens NOMAD shadow state is not active");

  PROTECT(rbw_r = coerceVector(rbw, REALSXP));
  PROTECT(degree_i = coerceVector(glp_degree, INTSXP));

  if (XLENGTH(rbw_r) != np_conditional_density_nomad_shadow.num_all_var) {
    UNPROTECT(2);
    error("resident npcdens NOMAD shadow received bandwidth vector of unexpected length");
  }
  if (XLENGTH(degree_i) != np_conditional_density_nomad_shadow.num_reg_continuous) {
    UNPROTECT(2);
    error("resident npcdens NOMAD shadow received degree vector of unexpected length");
  }

  if (!np_conditional_density_nomad_shadow_refresh_degree(INTEGER(degree_i))) {
    bwm_eval_count += 1.0;
    bwm_invalid_count += 1.0;
    bwm_maybe_signal_activity();
    val = (bwm_penalty_mode == 1 && R_FINITE(bwm_penalty_value)) ? bwm_penalty_value : DBL_MAX;
    PROTECT(out = allocVector(REALSXP, 3));
    REAL(out)[0] = -val;
    REAL(out)[1] = 1.0;
    REAL(out)[2] = 0.0;
    UNPROTECT(3);
    return out;
  }

  for (i = 0; i < np_conditional_density_nomad_shadow.num_all_var; i++)
    np_conditional_density_nomad_shadow.vector_scale_factor[i + 1] = REAL(rbw_r)[i];

  if (np_conditional_density_nomad_shadow.penalty_mode == 1) {
    bwm_penalty_mode = 0;
    bwm_penalty_value = DBL_MAX;
  }
  eval_before = bwm_eval_count;
  fast_before = np_fastcv_alllarge_hits_get();
  val = bwmfunc_wrapper(np_conditional_density_nomad_shadow.vector_scale_factor);
  fast = np_fastcv_alllarge_hits_get() - fast_before;
  evals = bwm_eval_count - eval_before;
  if (fast < 0.0)
    fast = 0.0;
  if (evals < 0.0)
    evals = 0.0;

  if ((!R_FINITE(val) || val == DBL_MAX) &&
      np_conditional_density_nomad_shadow.penalty_mode == 1) {
    np_conditional_density_nomad_shadow_refresh_penalty();
    val = (bwm_penalty_mode == 1 && R_FINITE(bwm_penalty_value)) ? bwm_penalty_value : DBL_MAX;
    evals = bwm_eval_count - eval_before;
    if (evals < 0.0)
      evals = 0.0;
  }

  if ((BANDWIDTH_den_extern == BW_FIXED) &&
      (!np_bw_candidate_is_admissible(
          np_conditional_density_nomad_shadow.num_all_var,
          bwm_use_transform,
          KERNEL_den_extern,
          KERNEL_reg_unordered_extern,
          BANDWIDTH_den_extern,
          BANDWIDTH_den_extern,
          0,
          num_obs_train_extern,
          num_var_continuous_extern,
          num_var_unordered_extern,
          num_var_ordered_extern,
          num_reg_continuous_extern,
          num_reg_unordered_extern,
          num_reg_ordered_extern,
          num_categories_extern,
          np_conditional_density_nomad_shadow.vector_scale_factor))) {
    val = DBL_MAX;
    fast = 0.0;
  }

  PROTECT(out = allocVector(REALSXP, 3));
  REAL(out)[0] = -val;
  REAL(out)[1] = evals;
  REAL(out)[2] = fast;
  UNPROTECT(3);
  return out;
}

SEXP C_np_density_conditional_nomad_shadow_clear(void)
{
  np_conditional_density_nomad_shadow_clear_internal();
  return R_NilValue;
}

static void np_regression_bw_mode(double * runo, double * rord, double * rcon, double * y,
                                  double * mysd, int * myopti, double * myoptd, double * rbw, double * fval,
                                  double * objective_function_values, double * objective_function_evals,
                                  double * objective_function_invalid, double * timing,
                                  double * objective_function_fast,
                                  int * penalty_mode, double * penalty_mult,
                                  int * glp_degree,
                                  int * glp_bernstein,
                                  int * glp_basis,
                                  double * ckerlb, double * ckerub,
                                  const int eval_only);

void np_regression(double * tuno, double * tord, double * tcon, double * ty,
                   double * euno, double * eord, double * econ, double * ey,
                   double * rbw,
                   double * mcv, double * padnum,
                   double * nconfac, double * ncatfac, double * mysd,
                   int * myopti,
                   int * glp_degree,
                   int * glp_gradient_order,
                   int * glp_bernstein,
                   int * glp_basis,
                   double * cm, double * cmerr, double * g, double *gerr,
                   double * xtra,
                   double * ckerlb, double * ckerub);

void np_density(double * tuno, double * tord, double * tcon,
                double * euno, double * eord, double * econ,
                double * rbw,
                double * mcv, double * padnum,
                double * nconfac, double * ncatfac, double * mysd,
                int * myopti,
                double * mydens, double * myderr, double * ll,
                double * ckerlb, double * ckerub);

void np_density_conditional(double * tyuno, double * tyord, double * tycon,
                            double * txuno, double * txord, double * txcon,
                            double * eyuno, double * eyord, double * eycon,
                            double * exuno, double * exord, double * excon,
                            double * rbw,
                            double * ymcv, double * ypadnum, double * xmcv, double * xpadnum,
                            double * nconfac, double * ncatfac, double * mysd,
                            int * myopti,
                            double * cmean, double * cmean_stderr,
                            double * gradients, double * gradients_stderr, double * ll,
                            double * cxkerlb, double * cxkerub,
                            double * cykerlb, double * cykerub);
void np_density_bw(double * myuno, double * myord, double * mycon,
                   double * mysd, int * myopti, double * myoptd, double * myans, double * fval,
                   double * objective_function_values, double * objective_function_evals,
                   double * objective_function_invalid, double * timing,
                   double * objective_function_fast,
                   int * penalty_mode, double * penalty_mult,
                   double * ckerlb, double * ckerub);
void np_distribution_bw(double * myuno, double * myord, double * mycon,
                        double * myeuno, double * myeord, double * myecon, double * mysd,
                        int * myopti, double * myoptd, double * myans, double * fval,
                        double * objective_function_values, double * objective_function_evals,
                        double * objective_function_invalid, double * timing,
                        double * objective_function_fast,
                        int * penalty_mode, double * penalty_mult,
                        double * ckerlb, double * ckerub);
void np_density_conditional_bw(double * c_uno, double * c_ord, double * c_con,
                               double * u_uno, double * u_ord, double * u_con,
                               double * mysd,
                               int * myopti, double * myoptd, double * myans, double * fval,
                               double * objective_function_values, double * objective_function_evals,
                               double * objective_function_invalid, double * timing,
                               double * objective_function_fast,
                               int * penalty_mode, double * penalty_mult,
                               int * glp_degree,
                               int * glp_bernstein,
                               int * glp_basis,
                               int * regtype,
                               double * cxkerlb, double * cxkerub,
                               double * cykerlb, double * cykerub,
                               const int eval_only);
void np_distribution_conditional_bw(double * c_uno, double * c_ord, double * c_con,
                                    double * u_uno, double * u_ord, double * u_con,
                                    double * cg_uno, double * cg_ord, double * cg_con, double * mysd,
                                    int * myopti, double * myoptd, double * myans, double * fval,
                                    double * objective_function_values, double * objective_function_evals,
                                    double * objective_function_invalid, double * timing,
                                    double * objective_function_fast,
                                    int * penalty_mode, double * penalty_mult,
                                    int * glp_degree,
                                    int * glp_bernstein,
                                    int * glp_basis,
                                    int * regtype,
                                    double * cxkerlb, double * cxkerub,
                                    double * cykerlb, double * cykerub,
                                    const int eval_only);
void np_kernelsum(double * tuno, double * tord, double * tcon,
                  double * ty, double * weights,
                  double * euno, double * eord, double * econ,
                  double * bw,
                  double * mcv, double * padnum,
                  int * operator,
                  int * myopti, double * kpow,
                  double * weighted_sum, double * weighted_p_sum,
                  double * kernel_weights,
                  double * permutation_kernel_weights,
                  double * ckerlb, double * ckerub);
void np_quantile_conditional(double * tc_con,
                             double * tu_uno, double * tu_ord, double * tu_con,
                             double * eu_uno, double * eu_ord, double * eu_con,
                             double * quantile,
                             double * mybw,
                             double * mcv, double *padnum,
                             double * nconfac, double * ncatfac, double * mysd,
                             int * myopti, double * myoptd,
                             double * yq, double * yqerr, double *yg);

static SEXP C_np_regression_bw_common(SEXP runo,
                                      SEXP rord,
                                      SEXP rcon,
                                      SEXP y,
                                      SEXP mysd,
                                      SEXP myopti,
                                      SEXP myoptd,
                                      SEXP rbw,
                                      SEXP hist_len,
                                      SEXP penalty_mode,
                                      SEXP penalty_mult,
                                      SEXP glp_degree,
                                      SEXP glp_bernstein,
                                      SEXP glp_basis,
                                      SEXP ckerlb,
                                      SEXP ckerub,
                                      const int eval_only)
{
  SEXP runo_r = R_NilValue, rord_r = R_NilValue, rcon_r = R_NilValue;
  SEXP y_r = R_NilValue, mysd_r = R_NilValue, myopti_i = R_NilValue, myoptd_r = R_NilValue;
  SEXP rbw_r = R_NilValue, degree_i = R_NilValue, ckerlb_r = R_NilValue, ckerub_r = R_NilValue;
  SEXP out = R_NilValue, out_names = R_NilValue;
  SEXP out_bw = R_NilValue, out_fval = R_NilValue, out_fval_hist = R_NilValue;
  SEXP out_eval_hist = R_NilValue, out_invalid_hist = R_NilValue, out_timing = R_NilValue;
  SEXP out_fast = R_NilValue;
  int hlen = asInteger(hist_len);
  int pmode = asInteger(penalty_mode);
  double pmult = asReal(penalty_mult);
  int bern = asInteger(glp_bernstein);
  int basis = asInteger(glp_basis);
  int ncon = 0;
  double * ckerlb_p = NULL;
  double * ckerub_p = NULL;

  if (hlen < 1)
    hlen = 1;

  PROTECT(runo_r = coerceVector(runo, REALSXP));
  PROTECT(rord_r = coerceVector(rord, REALSXP));
  PROTECT(rcon_r = coerceVector(rcon, REALSXP));
  PROTECT(y_r = coerceVector(y, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(myoptd_r = coerceVector(myoptd, REALSXP));
  PROTECT(rbw_r = coerceVector(rbw, REALSXP));
  PROTECT(degree_i = coerceVector(glp_degree, INTSXP));
  PROTECT(ckerlb_r = coerceVector(ckerlb, REALSXP));
  PROTECT(ckerub_r = coerceVector(ckerub, REALSXP));

  if (XLENGTH(myoptd_r) <= RBW_SFLOORD)
    error("C_np_regression_bw: myoptd is missing scale.factor.lower.bound");

  ncon = (int)INTEGER(myopti_i)[REG_NCONI];
  resolve_bounds_or_default(ckerlb_r, ckerub_r, ncon, &ckerlb_p, &ckerub_p);

  PROTECT(out_bw = allocVector(REALSXP, XLENGTH(rbw_r)));
  PROTECT(out_fval = allocVector(REALSXP, 3));
  PROTECT(out_fval_hist = allocVector(REALSXP, hlen));
  PROTECT(out_eval_hist = allocVector(REALSXP, hlen));
  PROTECT(out_invalid_hist = allocVector(REALSXP, hlen));
  PROTECT(out_timing = allocVector(REALSXP, 1));
  PROTECT(out_fast = allocVector(REALSXP, 1));

  memcpy(REAL(out_bw), REAL(rbw_r), (size_t)XLENGTH(rbw_r) * sizeof(double));

  np_regression_bw_mode(REAL(runo_r), REAL(rord_r), REAL(rcon_r), REAL(y_r),
                        REAL(mysd_r), INTEGER(myopti_i), REAL(myoptd_r), REAL(out_bw), REAL(out_fval),
                        REAL(out_fval_hist), REAL(out_eval_hist), REAL(out_invalid_hist), REAL(out_timing),
                        REAL(out_fast),
                        &pmode, &pmult, INTEGER(degree_i), &bern, &basis, ckerlb_p, ckerub_p,
                        eval_only);

  PROTECT(out = allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, out_bw);
  SET_VECTOR_ELT(out, 1, out_fval);
  SET_VECTOR_ELT(out, 2, out_fval_hist);
  SET_VECTOR_ELT(out, 3, out_eval_hist);
  SET_VECTOR_ELT(out, 4, out_invalid_hist);
  SET_VECTOR_ELT(out, 5, out_timing);
  SET_VECTOR_ELT(out, 6, out_fast);

  PROTECT(out_names = allocVector(STRSXP, 7));
  SET_STRING_ELT(out_names, 0, mkChar("bw"));
  SET_STRING_ELT(out_names, 1, mkChar("fval"));
  SET_STRING_ELT(out_names, 2, mkChar("fval.history"));
  SET_STRING_ELT(out_names, 3, mkChar("eval.history"));
  SET_STRING_ELT(out_names, 4, mkChar("invalid.history"));
  SET_STRING_ELT(out_names, 5, mkChar("timing"));
  SET_STRING_ELT(out_names, 6, mkChar("fast.history"));
  setAttrib(out, R_NamesSymbol, out_names);

  UNPROTECT(20);
  return out;
}

SEXP C_np_regression_bw(SEXP runo,
                        SEXP rord,
                        SEXP rcon,
                        SEXP y,
                        SEXP mysd,
                        SEXP myopti,
                        SEXP myoptd,
                        SEXP rbw,
                        SEXP hist_len,
                        SEXP penalty_mode,
                        SEXP penalty_mult,
                        SEXP glp_degree,
                        SEXP glp_bernstein,
                        SEXP glp_basis,
                        SEXP ckerlb,
                        SEXP ckerub)
{
  return C_np_regression_bw_common(runo, rord, rcon, y, mysd, myopti, myoptd,
                                   rbw, hist_len, penalty_mode, penalty_mult,
                                   glp_degree, glp_bernstein, glp_basis,
                                   ckerlb, ckerub, 0);
}

SEXP C_np_regression_bw_eval(SEXP runo,
                             SEXP rord,
                             SEXP rcon,
                             SEXP y,
                             SEXP mysd,
                             SEXP myopti,
                             SEXP myoptd,
                             SEXP rbw,
                             SEXP hist_len,
                             SEXP penalty_mode,
                             SEXP penalty_mult,
                             SEXP glp_degree,
                             SEXP glp_bernstein,
                             SEXP glp_basis,
                             SEXP ckerlb,
                             SEXP ckerub)
{
  return C_np_regression_bw_common(runo, rord, rcon, y, mysd, myopti, myoptd,
                                   rbw, hist_len, penalty_mode, penalty_mult,
                                   glp_degree, glp_bernstein, glp_basis,
                                   ckerlb, ckerub, 1);
}

SEXP C_np_regression(SEXP tuno,
                     SEXP tord,
                     SEXP tcon,
                     SEXP ty,
                     SEXP euno,
                     SEXP eord,
                     SEXP econ,
                     SEXP ey,
                     SEXP rbw,
                     SEXP mcv,
                     SEXP padnum,
                     SEXP nconfac,
                     SEXP ncatfac,
                     SEXP mysd,
                     SEXP myopti,
                     SEXP glp_degree,
                     SEXP glp_gradient_order,
                     SEXP glp_bernstein,
                     SEXP glp_basis,
                     SEXP enrow,
                     SEXP ncol,
                     SEXP gradients,
                     SEXP ckerlb,
                     SEXP ckerub)
{
  SEXP tuno_r = R_NilValue, tord_r = R_NilValue, tcon_r = R_NilValue, ty_r = R_NilValue;
  SEXP euno_r = R_NilValue, eord_r = R_NilValue, econ_r = R_NilValue, ey_r = R_NilValue;
  SEXP rbw_r = R_NilValue, mcv_r = R_NilValue, padnum_r = R_NilValue;
  SEXP nconfac_r = R_NilValue, ncatfac_r = R_NilValue, mysd_r = R_NilValue, myopti_i = R_NilValue;
  SEXP degree_i = R_NilValue, gradient_order_i = R_NilValue, ckerlb_r = R_NilValue, ckerub_r = R_NilValue;
  SEXP out = R_NilValue, out_names = R_NilValue;
  SEXP out_mean = R_NilValue, out_merr = R_NilValue, out_g = R_NilValue, out_gerr = R_NilValue, out_xtra = R_NilValue;
  int bern = asInteger(glp_bernstein);
  int basis = asInteger(glp_basis);
  int en = asInteger(enrow);
  int nc = asInteger(ncol);
  int do_grad = asLogical(gradients);
  int ncon = 0;
  double * ckerlb_p = NULL;
  double * ckerub_p = NULL;
  R_xlen_t gsize;

  if (en < 0) en = 0;
  if (nc < 0) nc = 0;
  if (do_grad == NA_LOGICAL) do_grad = 0;
  gsize = (do_grad ? ((R_xlen_t)en * (R_xlen_t)nc) : 1);

  PROTECT(tuno_r = coerceVector(tuno, REALSXP));
  PROTECT(tord_r = coerceVector(tord, REALSXP));
  PROTECT(tcon_r = coerceVector(tcon, REALSXP));
  PROTECT(ty_r = coerceVector(ty, REALSXP));
  PROTECT(euno_r = coerceVector(euno, REALSXP));
  PROTECT(eord_r = coerceVector(eord, REALSXP));
  PROTECT(econ_r = coerceVector(econ, REALSXP));
  PROTECT(ey_r = coerceVector(ey, REALSXP));
  PROTECT(rbw_r = coerceVector(rbw, REALSXP));
  PROTECT(mcv_r = coerceVector(mcv, REALSXP));
  PROTECT(padnum_r = coerceVector(padnum, REALSXP));
  PROTECT(nconfac_r = coerceVector(nconfac, REALSXP));
  PROTECT(ncatfac_r = coerceVector(ncatfac, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(degree_i = coerceVector(glp_degree, INTSXP));
  PROTECT(gradient_order_i = coerceVector(glp_gradient_order, INTSXP));
  PROTECT(ckerlb_r = coerceVector(ckerlb, REALSXP));
  PROTECT(ckerub_r = coerceVector(ckerub, REALSXP));

  ncon = (int)INTEGER(myopti_i)[REG_NCONI];
  resolve_bounds_or_default(ckerlb_r, ckerub_r, ncon, &ckerlb_p, &ckerub_p);

  PROTECT(out_mean = allocVector(REALSXP, en));
  PROTECT(out_merr = allocVector(REALSXP, en));
  PROTECT(out_g = allocVector(REALSXP, gsize));
  PROTECT(out_gerr = allocVector(REALSXP, gsize));
  PROTECT(out_xtra = allocVector(REALSXP, 6));

  np_regression(REAL(tuno_r), REAL(tord_r), REAL(tcon_r), REAL(ty_r),
                REAL(euno_r), REAL(eord_r), REAL(econ_r), REAL(ey_r),
                REAL(rbw_r), REAL(mcv_r), REAL(padnum_r),
                REAL(nconfac_r), REAL(ncatfac_r), REAL(mysd_r),
                INTEGER(myopti_i),
                INTEGER(degree_i), INTEGER(gradient_order_i), &bern, &basis,
                REAL(out_mean), REAL(out_merr), REAL(out_g), REAL(out_gerr), REAL(out_xtra),
                ckerlb_p, ckerub_p);

  PROTECT(out = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, out_mean);
  SET_VECTOR_ELT(out, 1, out_merr);
  SET_VECTOR_ELT(out, 2, out_g);
  SET_VECTOR_ELT(out, 3, out_gerr);
  SET_VECTOR_ELT(out, 4, out_xtra);

  PROTECT(out_names = allocVector(STRSXP, 5));
  SET_STRING_ELT(out_names, 0, mkChar("mean"));
  SET_STRING_ELT(out_names, 1, mkChar("merr"));
  SET_STRING_ELT(out_names, 2, mkChar("g"));
  SET_STRING_ELT(out_names, 3, mkChar("gerr"));
  SET_STRING_ELT(out_names, 4, mkChar("xtra"));
  setAttrib(out, R_NamesSymbol, out_names);

  UNPROTECT(26);
  return out;
}

SEXP C_np_density(SEXP tuno,
                  SEXP tord,
                  SEXP tcon,
                  SEXP euno,
                  SEXP eord,
                  SEXP econ,
                  SEXP rbw,
                  SEXP mcv,
                  SEXP padnum,
                  SEXP nconfac,
                  SEXP ncatfac,
                  SEXP mysd,
                  SEXP myopti,
                  SEXP enrow,
                  SEXP ckerlb,
                  SEXP ckerub)
{
  SEXP tuno_r=R_NilValue, tord_r=R_NilValue, tcon_r=R_NilValue;
  SEXP euno_r=R_NilValue, eord_r=R_NilValue, econ_r=R_NilValue;
  SEXP rbw_r=R_NilValue, mcv_r=R_NilValue, padnum_r=R_NilValue;
  SEXP nconfac_r=R_NilValue, ncatfac_r=R_NilValue, mysd_r=R_NilValue, myopti_i=R_NilValue;
  SEXP ckerlb_r=R_NilValue, ckerub_r=R_NilValue;
  SEXP out=R_NilValue, out_names=R_NilValue, out_dens=R_NilValue, out_derr=R_NilValue, out_ll=R_NilValue;
  int en = asInteger(enrow);
  int ncon = 0;
  double * ckerlb_p = NULL;
  double * ckerub_p = NULL;

  if (en < 0) en = 0;

  PROTECT(tuno_r = coerceVector(tuno, REALSXP));
  PROTECT(tord_r = coerceVector(tord, REALSXP));
  PROTECT(tcon_r = coerceVector(tcon, REALSXP));
  PROTECT(euno_r = coerceVector(euno, REALSXP));
  PROTECT(eord_r = coerceVector(eord, REALSXP));
  PROTECT(econ_r = coerceVector(econ, REALSXP));
  PROTECT(rbw_r = coerceVector(rbw, REALSXP));
  PROTECT(mcv_r = coerceVector(mcv, REALSXP));
  PROTECT(padnum_r = coerceVector(padnum, REALSXP));
  PROTECT(nconfac_r = coerceVector(nconfac, REALSXP));
  PROTECT(ncatfac_r = coerceVector(ncatfac, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(ckerlb_r = coerceVector(ckerlb, REALSXP));
  PROTECT(ckerub_r = coerceVector(ckerub, REALSXP));

  ncon = (int)INTEGER(myopti_i)[DEN_NCONI];
  resolve_bounds_or_default(ckerlb_r, ckerub_r, ncon, &ckerlb_p, &ckerub_p);

  PROTECT(out_dens = allocVector(REALSXP, en));
  PROTECT(out_derr = allocVector(REALSXP, en));
  PROTECT(out_ll = allocVector(REALSXP, 1));

  np_density(REAL(tuno_r), REAL(tord_r), REAL(tcon_r),
             REAL(euno_r), REAL(eord_r), REAL(econ_r),
             REAL(rbw_r),
             REAL(mcv_r), REAL(padnum_r),
             REAL(nconfac_r), REAL(ncatfac_r), REAL(mysd_r),
             INTEGER(myopti_i),
             REAL(out_dens), REAL(out_derr), REAL(out_ll),
             ckerlb_p, ckerub_p);

  PROTECT(out = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, out_dens);
  SET_VECTOR_ELT(out, 1, out_derr);
  SET_VECTOR_ELT(out, 2, out_ll);

  PROTECT(out_names = allocVector(STRSXP, 3));
  SET_STRING_ELT(out_names, 0, mkChar("dens"));
  SET_STRING_ELT(out_names, 1, mkChar("derr"));
  SET_STRING_ELT(out_names, 2, mkChar("log_likelihood"));
  setAttrib(out, R_NamesSymbol, out_names);

  UNPROTECT(20);
  return out;
}

SEXP C_np_density_conditional(SEXP tyuno,
                              SEXP tyord,
                              SEXP tycon,
                              SEXP txuno,
                              SEXP txord,
                              SEXP txcon,
                              SEXP eyuno,
                              SEXP eyord,
                              SEXP eycon,
                              SEXP exuno,
                              SEXP exord,
                              SEXP excon,
                              SEXP rbw,
                              SEXP ymcv,
                              SEXP ypadnum,
                              SEXP xmcv,
                              SEXP xpadnum,
                              SEXP nconfac,
                              SEXP ncatfac,
                              SEXP mysd,
                              SEXP myopti,
                              SEXP enrow,
                              SEXP xndim,
                              SEXP ckerlbx,
                              SEXP ckerubx,
                              SEXP ckerlby,
                              SEXP ckeruby,
                              SEXP regtype,
                              SEXP glp_degree,
                              SEXP glp_bernstein,
                              SEXP glp_basis)
{
  SEXP tyuno_r=R_NilValue, tyord_r=R_NilValue, tycon_r=R_NilValue;
  SEXP txuno_r=R_NilValue, txord_r=R_NilValue, txcon_r=R_NilValue;
  SEXP eyuno_r=R_NilValue, eyord_r=R_NilValue, eycon_r=R_NilValue;
  SEXP exuno_r=R_NilValue, exord_r=R_NilValue, excon_r=R_NilValue;
  SEXP rbw_r=R_NilValue, ymcv_r=R_NilValue, ypadnum_r=R_NilValue, xmcv_r=R_NilValue, xpadnum_r=R_NilValue;
  SEXP nconfac_r=R_NilValue, ncatfac_r=R_NilValue, mysd_r=R_NilValue, myopti_i=R_NilValue;
  SEXP ckerlbx_r=R_NilValue, ckerubx_r=R_NilValue, ckerlby_r=R_NilValue, ckeruby_r=R_NilValue;
  SEXP regtype_i=R_NilValue, glp_degree_i=R_NilValue, glp_bernstein_i=R_NilValue, glp_basis_i=R_NilValue;
  SEXP out=R_NilValue, out_names=R_NilValue;
  SEXP out_cond=R_NilValue, out_cderr=R_NilValue, out_grad=R_NilValue, out_gerr=R_NilValue, out_ll=R_NilValue;
  int en = asInteger(enrow);
  int xd = asInteger(xndim);
  int ncon_x = 0;
  int ncon_y = 0;
  double * cxkerlb_p = NULL;
  double * cxkerub_p = NULL;
  double * cykerlb_p = NULL;
  double * cykerub_p = NULL;
  R_xlen_t gsize;

  if (en < 0) en = 0;
  if (xd < 0) xd = 0;
  gsize = (R_xlen_t)en * (R_xlen_t)xd;

  PROTECT(tyuno_r = coerceVector(tyuno, REALSXP));
  PROTECT(tyord_r = coerceVector(tyord, REALSXP));
  PROTECT(tycon_r = coerceVector(tycon, REALSXP));
  PROTECT(txuno_r = coerceVector(txuno, REALSXP));
  PROTECT(txord_r = coerceVector(txord, REALSXP));
  PROTECT(txcon_r = coerceVector(txcon, REALSXP));
  PROTECT(eyuno_r = coerceVector(eyuno, REALSXP));
  PROTECT(eyord_r = coerceVector(eyord, REALSXP));
  PROTECT(eycon_r = coerceVector(eycon, REALSXP));
  PROTECT(exuno_r = coerceVector(exuno, REALSXP));
  PROTECT(exord_r = coerceVector(exord, REALSXP));
  PROTECT(excon_r = coerceVector(excon, REALSXP));
  PROTECT(rbw_r = coerceVector(rbw, REALSXP));
  PROTECT(ymcv_r = coerceVector(ymcv, REALSXP));
  PROTECT(ypadnum_r = coerceVector(ypadnum, REALSXP));
  PROTECT(xmcv_r = coerceVector(xmcv, REALSXP));
  PROTECT(xpadnum_r = coerceVector(xpadnum, REALSXP));
  PROTECT(nconfac_r = coerceVector(nconfac, REALSXP));
  PROTECT(ncatfac_r = coerceVector(ncatfac, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(ckerlbx_r = coerceVector(ckerlbx, REALSXP));
  PROTECT(ckerubx_r = coerceVector(ckerubx, REALSXP));
  PROTECT(ckerlby_r = coerceVector(ckerlby, REALSXP));
  PROTECT(ckeruby_r = coerceVector(ckeruby, REALSXP));
  PROTECT(regtype_i = coerceVector(regtype, INTSXP));
  PROTECT(glp_degree_i = coerceVector(glp_degree, INTSXP));
  PROTECT(glp_bernstein_i = coerceVector(glp_bernstein, INTSXP));
  PROTECT(glp_basis_i = coerceVector(glp_basis, INTSXP));

  ncon_x = (int)INTEGER(myopti_i)[CD_UNCONI];
  ncon_y = (int)INTEGER(myopti_i)[CD_CNCONI];
  resolve_bounds_or_default(ckerlbx_r, ckerubx_r, ncon_x, &cxkerlb_p, &cxkerub_p);
  resolve_bounds_or_default(ckerlby_r, ckeruby_r, ncon_y, &cykerlb_p, &cykerub_p);

  int_ll_extern = asInteger(regtype_i);
  if ((int_ll_extern == LL_LP) && (ncon_x > 0)) {
    if ((int)XLENGTH(glp_degree_i) != ncon_x)
      error("C_np_density_conditional: length(glp_degree) must equal number of continuous x variables");
    vector_glp_degree_extern = INTEGER(glp_degree_i);
    int_glp_bernstein_extern = asInteger(glp_bernstein_i);
    int_glp_basis_extern = asInteger(glp_basis_i);
  } else {
    vector_glp_degree_extern = NULL;
    int_glp_bernstein_extern = 0;
    int_glp_basis_extern = 1;
  }

  PROTECT(out_cond = allocVector(REALSXP, en));
  PROTECT(out_cderr = allocVector(REALSXP, en));
  PROTECT(out_grad = allocVector(REALSXP, gsize));
  PROTECT(out_gerr = allocVector(REALSXP, gsize));
  PROTECT(out_ll = allocVector(REALSXP, 1));

  np_density_conditional(REAL(tyuno_r), REAL(tyord_r), REAL(tycon_r),
                         REAL(txuno_r), REAL(txord_r), REAL(txcon_r),
                         REAL(eyuno_r), REAL(eyord_r), REAL(eycon_r),
                         REAL(exuno_r), REAL(exord_r), REAL(excon_r),
                         REAL(rbw_r),
                         REAL(ymcv_r), REAL(ypadnum_r), REAL(xmcv_r), REAL(xpadnum_r),
                         REAL(nconfac_r), REAL(ncatfac_r), REAL(mysd_r),
                         INTEGER(myopti_i),
                         REAL(out_cond), REAL(out_cderr), REAL(out_grad), REAL(out_gerr), REAL(out_ll),
                         cxkerlb_p, cxkerub_p, cykerlb_p, cykerub_p);

  PROTECT(out = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, out_cond);
  SET_VECTOR_ELT(out, 1, out_cderr);
  SET_VECTOR_ELT(out, 2, out_grad);
  SET_VECTOR_ELT(out, 3, out_gerr);
  SET_VECTOR_ELT(out, 4, out_ll);

  PROTECT(out_names = allocVector(STRSXP, 5));
  SET_STRING_ELT(out_names, 0, mkChar("condens"));
  SET_STRING_ELT(out_names, 1, mkChar("conderr"));
  SET_STRING_ELT(out_names, 2, mkChar("congrad"));
  SET_STRING_ELT(out_names, 3, mkChar("congerr"));
  SET_STRING_ELT(out_names, 4, mkChar("log_likelihood"));
  setAttrib(out, R_NamesSymbol, out_names);
  vector_glp_degree_extern = NULL;
  int_glp_bernstein_extern = 0;
  int_glp_basis_extern = 1;
  int_ll_extern = LL_LC;

  UNPROTECT(36);
  return out;
}

SEXP C_np_density_bw(SEXP myuno,
                     SEXP myord,
                     SEXP mycon,
                     SEXP mysd,
                     SEXP myopti,
                     SEXP myoptd,
                     SEXP bw,
                     SEXP hist_len,
                     SEXP penalty_mode,
                     SEXP penalty_mult,
                     SEXP ckerlb,
                     SEXP ckerub)
{
  SEXP myuno_r=R_NilValue, myord_r=R_NilValue, mycon_r=R_NilValue, mysd_r=R_NilValue;
  SEXP myopti_i=R_NilValue, myoptd_r=R_NilValue, bw_r=R_NilValue, ckerlb_r=R_NilValue, ckerub_r=R_NilValue;
  SEXP out=R_NilValue, out_names=R_NilValue;
  SEXP out_bw=R_NilValue, out_fval=R_NilValue, out_fval_hist=R_NilValue, out_eval_hist=R_NilValue;
  SEXP out_invalid_hist=R_NilValue, out_timing=R_NilValue, out_fast=R_NilValue;
  int hlen = asInteger(hist_len);
  int pmode = asInteger(penalty_mode);
  double pmult = asReal(penalty_mult);
  int ncon = 0;
  double * ckerlb_p = NULL;
  double * ckerub_p = NULL;

  if(hlen < 1) hlen = 1;

  PROTECT(myuno_r = coerceVector(myuno, REALSXP));
  PROTECT(myord_r = coerceVector(myord, REALSXP));
  PROTECT(mycon_r = coerceVector(mycon, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(myoptd_r = coerceVector(myoptd, REALSXP));
  PROTECT(bw_r = coerceVector(bw, REALSXP));
  PROTECT(ckerlb_r = coerceVector(ckerlb, REALSXP));
  PROTECT(ckerub_r = coerceVector(ckerub, REALSXP));

  if (XLENGTH(myoptd_r) <= BW_SFLOORD)
    error("C_np_density_bw: myoptd is missing scale.factor.lower.bound");

  ncon = (int)INTEGER(myopti_i)[BW_NCONI];
  resolve_bounds_or_default(ckerlb_r, ckerub_r, ncon, &ckerlb_p, &ckerub_p);

  PROTECT(out_bw = allocVector(REALSXP, XLENGTH(bw_r)));
  PROTECT(out_fval = allocVector(REALSXP, 2));
  PROTECT(out_fval_hist = allocVector(REALSXP, hlen));
  PROTECT(out_eval_hist = allocVector(REALSXP, hlen));
  PROTECT(out_invalid_hist = allocVector(REALSXP, hlen));
  PROTECT(out_timing = allocVector(REALSXP, 1));
  PROTECT(out_fast = allocVector(REALSXP, 1));

  memcpy(REAL(out_bw), REAL(bw_r), (size_t)XLENGTH(bw_r) * sizeof(double));
  np_density_bw(REAL(myuno_r), REAL(myord_r), REAL(mycon_r),
                REAL(mysd_r), INTEGER(myopti_i), REAL(myoptd_r), REAL(out_bw), REAL(out_fval),
                REAL(out_fval_hist), REAL(out_eval_hist), REAL(out_invalid_hist), REAL(out_timing),
                REAL(out_fast),
                &pmode, &pmult, ckerlb_p, ckerub_p);

  PROTECT(out = allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, out_bw);
  SET_VECTOR_ELT(out, 1, out_fval);
  SET_VECTOR_ELT(out, 2, out_fval_hist);
  SET_VECTOR_ELT(out, 3, out_eval_hist);
  SET_VECTOR_ELT(out, 4, out_invalid_hist);
  SET_VECTOR_ELT(out, 5, out_timing);
  SET_VECTOR_ELT(out, 6, out_fast);

  PROTECT(out_names = allocVector(STRSXP, 7));
  SET_STRING_ELT(out_names, 0, mkChar("bw"));
  SET_STRING_ELT(out_names, 1, mkChar("fval"));
  SET_STRING_ELT(out_names, 2, mkChar("fval.history"));
  SET_STRING_ELT(out_names, 3, mkChar("eval.history"));
  SET_STRING_ELT(out_names, 4, mkChar("invalid.history"));
  SET_STRING_ELT(out_names, 5, mkChar("timing"));
  SET_STRING_ELT(out_names, 6, mkChar("fast.history"));
  setAttrib(out, R_NamesSymbol, out_names);

  UNPROTECT(18);
  return out;
}

SEXP C_np_distribution_bw(SEXP myuno,
                          SEXP myord,
                          SEXP mycon,
                          SEXP myeuno,
                          SEXP myeord,
                          SEXP myecon,
                          SEXP mysd,
                          SEXP myopti,
                          SEXP myoptd,
                          SEXP bw,
                          SEXP hist_len,
                          SEXP penalty_mode,
                          SEXP penalty_mult,
                          SEXP ckerlb,
                          SEXP ckerub)
{
  SEXP myuno_r=R_NilValue, myord_r=R_NilValue, mycon_r=R_NilValue;
  SEXP myeuno_r=R_NilValue, myeord_r=R_NilValue, myecon_r=R_NilValue, mysd_r=R_NilValue;
  SEXP myopti_i=R_NilValue, myoptd_r=R_NilValue, bw_r=R_NilValue, ckerlb_r=R_NilValue, ckerub_r=R_NilValue;
  SEXP out=R_NilValue, out_names=R_NilValue;
  SEXP out_bw=R_NilValue, out_fval=R_NilValue, out_fval_hist=R_NilValue, out_eval_hist=R_NilValue;
  SEXP out_invalid_hist=R_NilValue, out_timing=R_NilValue, out_fast=R_NilValue;
  int hlen = asInteger(hist_len);
  int pmode = asInteger(penalty_mode);
  double pmult = asReal(penalty_mult);
  int ncon = 0;
  double * ckerlb_p = NULL;
  double * ckerub_p = NULL;

  if(hlen < 1) hlen = 1;

  PROTECT(myuno_r = coerceVector(myuno, REALSXP));
  PROTECT(myord_r = coerceVector(myord, REALSXP));
  PROTECT(mycon_r = coerceVector(mycon, REALSXP));
  PROTECT(myeuno_r = coerceVector(myeuno, REALSXP));
  PROTECT(myeord_r = coerceVector(myeord, REALSXP));
  PROTECT(myecon_r = coerceVector(myecon, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(myoptd_r = coerceVector(myoptd, REALSXP));
  PROTECT(bw_r = coerceVector(bw, REALSXP));
  PROTECT(ckerlb_r = coerceVector(ckerlb, REALSXP));
  PROTECT(ckerub_r = coerceVector(ckerub, REALSXP));

  ncon = (int)INTEGER(myopti_i)[DBW_NCONI];
  resolve_bounds_or_default(ckerlb_r, ckerub_r, ncon, &ckerlb_p, &ckerub_p);

  PROTECT(out_bw = allocVector(REALSXP, XLENGTH(bw_r)));
  PROTECT(out_fval = allocVector(REALSXP, 3));
  PROTECT(out_fval_hist = allocVector(REALSXP, hlen));
  PROTECT(out_eval_hist = allocVector(REALSXP, hlen));
  PROTECT(out_invalid_hist = allocVector(REALSXP, hlen));
  PROTECT(out_timing = allocVector(REALSXP, 1));
  PROTECT(out_fast = allocVector(REALSXP, 1));

  memcpy(REAL(out_bw), REAL(bw_r), (size_t)XLENGTH(bw_r) * sizeof(double));
  np_distribution_bw(REAL(myuno_r), REAL(myord_r), REAL(mycon_r),
                     REAL(myeuno_r), REAL(myeord_r), REAL(myecon_r), REAL(mysd_r),
                     INTEGER(myopti_i), REAL(myoptd_r), REAL(out_bw), REAL(out_fval),
                     REAL(out_fval_hist), REAL(out_eval_hist), REAL(out_invalid_hist), REAL(out_timing),
                     REAL(out_fast),
                     &pmode, &pmult, ckerlb_p, ckerub_p);

  PROTECT(out = allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, out_bw);
  SET_VECTOR_ELT(out, 1, out_fval);
  SET_VECTOR_ELT(out, 2, out_fval_hist);
  SET_VECTOR_ELT(out, 3, out_eval_hist);
  SET_VECTOR_ELT(out, 4, out_invalid_hist);
  SET_VECTOR_ELT(out, 5, out_timing);
  SET_VECTOR_ELT(out, 6, out_fast);

  PROTECT(out_names = allocVector(STRSXP, 7));
  SET_STRING_ELT(out_names, 0, mkChar("bw"));
  SET_STRING_ELT(out_names, 1, mkChar("fval"));
  SET_STRING_ELT(out_names, 2, mkChar("fval.history"));
  SET_STRING_ELT(out_names, 3, mkChar("eval.history"));
  SET_STRING_ELT(out_names, 4, mkChar("invalid.history"));
  SET_STRING_ELT(out_names, 5, mkChar("timing"));
  SET_STRING_ELT(out_names, 6, mkChar("fast.history"));
  setAttrib(out, R_NamesSymbol, out_names);

  UNPROTECT(21);
  return out;
}

static SEXP C_np_density_conditional_bw_common(SEXP c_uno,
                                               SEXP c_ord,
                                               SEXP c_con,
                                               SEXP u_uno,
                                               SEXP u_ord,
                                               SEXP u_con,
                                               SEXP mysd,
                                               SEXP myopti,
                                               SEXP myoptd,
                                               SEXP bw,
                                               SEXP hist_len,
                                               SEXP penalty_mode,
                                               SEXP penalty_mult,
                                               SEXP glp_degree,
                                               SEXP glp_bernstein,
                                               SEXP glp_basis,
                                               SEXP regtype,
                                               SEXP cxkerlb,
                                               SEXP cxkerub,
                                               SEXP cykerlb,
                                               SEXP cykerub,
                                               const int eval_only)
{
  SEXP c_uno_r=R_NilValue, c_ord_r=R_NilValue, c_con_r=R_NilValue, u_uno_r=R_NilValue, u_ord_r=R_NilValue, u_con_r=R_NilValue;
  SEXP mysd_r=R_NilValue, myopti_i=R_NilValue, myoptd_r=R_NilValue, bw_r=R_NilValue;
  SEXP degree_i=R_NilValue;
  SEXP cxkerlb_r=R_NilValue, cxkerub_r=R_NilValue, cykerlb_r=R_NilValue, cykerub_r=R_NilValue;
  SEXP out=R_NilValue, out_names=R_NilValue;
  SEXP out_bw=R_NilValue, out_fval=R_NilValue, out_fval_hist=R_NilValue, out_eval_hist=R_NilValue;
  SEXP out_invalid_hist=R_NilValue, out_timing=R_NilValue, out_fast=R_NilValue;
  int hlen = asInteger(hist_len);
  int pmode = asInteger(penalty_mode);
  double pmult = asReal(penalty_mult);
  int bern = asInteger(glp_bernstein);
  int basis = asInteger(glp_basis);
  int ll_mode = asInteger(regtype);
  int ncon_x = 0;
  int ncon_y = 0;
  double * cxkerlb_p = NULL;
  double * cxkerub_p = NULL;
  double * cykerlb_p = NULL;
  double * cykerub_p = NULL;

  if(hlen < 1) hlen = 1;

  PROTECT(c_uno_r = coerceVector(c_uno, REALSXP));
  PROTECT(c_ord_r = coerceVector(c_ord, REALSXP));
  PROTECT(c_con_r = coerceVector(c_con, REALSXP));
  PROTECT(u_uno_r = coerceVector(u_uno, REALSXP));
  PROTECT(u_ord_r = coerceVector(u_ord, REALSXP));
  PROTECT(u_con_r = coerceVector(u_con, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(myoptd_r = coerceVector(myoptd, REALSXP));
  PROTECT(bw_r = coerceVector(bw, REALSXP));
  PROTECT(degree_i = coerceVector(glp_degree, INTSXP));
  PROTECT(cxkerlb_r = coerceVector(cxkerlb, REALSXP));
  PROTECT(cxkerub_r = coerceVector(cxkerub, REALSXP));
  PROTECT(cykerlb_r = coerceVector(cykerlb, REALSXP));
  PROTECT(cykerub_r = coerceVector(cykerub, REALSXP));

  if (XLENGTH(myopti_i) <= CBW_CVLS_QUAD_POINTSI)
    error("C_np_density_conditional_bw: myopti is missing cvls.quadrature.points");
  if (XLENGTH(myoptd_r) <= CBW_QUAD_EXTD)
    error("C_np_density_conditional_bw: myoptd is missing cvls.quadrature.extend.factor");

  ncon_x = (int)INTEGER(myopti_i)[CDBW_UNCONI];
  ncon_y = (int)INTEGER(myopti_i)[CDBW_CNCONI];
  resolve_bounds_or_default(cxkerlb_r, cxkerub_r, ncon_x, &cxkerlb_p, &cxkerub_p);
  resolve_bounds_or_default(cykerlb_r, cykerub_r, ncon_y, &cykerlb_p, &cykerub_p);

  PROTECT(out_bw = allocVector(REALSXP, XLENGTH(bw_r)));
  PROTECT(out_fval = allocVector(REALSXP, 3));
  PROTECT(out_fval_hist = allocVector(REALSXP, hlen));
  PROTECT(out_eval_hist = allocVector(REALSXP, hlen));
  PROTECT(out_invalid_hist = allocVector(REALSXP, hlen));
  PROTECT(out_timing = allocVector(REALSXP, 1));
  PROTECT(out_fast = allocVector(REALSXP, 1));

  memcpy(REAL(out_bw), REAL(bw_r), (size_t)XLENGTH(bw_r) * sizeof(double));
  np_density_conditional_bw(REAL(c_uno_r), REAL(c_ord_r), REAL(c_con_r),
                            REAL(u_uno_r), REAL(u_ord_r), REAL(u_con_r), REAL(mysd_r),
                            INTEGER(myopti_i), REAL(myoptd_r), REAL(out_bw), REAL(out_fval),
                            REAL(out_fval_hist), REAL(out_eval_hist), REAL(out_invalid_hist), REAL(out_timing),
                            REAL(out_fast), &pmode, &pmult,
                            INTEGER(degree_i), &bern, &basis, &ll_mode,
                            cxkerlb_p, cxkerub_p, cykerlb_p, cykerub_p,
                            eval_only);

  PROTECT(out = allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, out_bw);
  SET_VECTOR_ELT(out, 1, out_fval);
  SET_VECTOR_ELT(out, 2, out_fval_hist);
  SET_VECTOR_ELT(out, 3, out_eval_hist);
  SET_VECTOR_ELT(out, 4, out_invalid_hist);
  SET_VECTOR_ELT(out, 5, out_timing);
  SET_VECTOR_ELT(out, 6, out_fast);

  PROTECT(out_names = allocVector(STRSXP, 7));
  SET_STRING_ELT(out_names, 0, mkChar("bw"));
  SET_STRING_ELT(out_names, 1, mkChar("fval"));
  SET_STRING_ELT(out_names, 2, mkChar("fval.history"));
  SET_STRING_ELT(out_names, 3, mkChar("eval.history"));
  SET_STRING_ELT(out_names, 4, mkChar("invalid.history"));
  SET_STRING_ELT(out_names, 5, mkChar("timing"));
  SET_STRING_ELT(out_names, 6, mkChar("fast.history"));
  setAttrib(out, R_NamesSymbol, out_names);

  UNPROTECT(24);
  return out;
}

static SEXP C_np_distribution_conditional_bw_common(SEXP c_uno,
                                                    SEXP c_ord,
                                                    SEXP c_con,
                                                    SEXP u_uno,
                                                    SEXP u_ord,
                                                    SEXP u_con,
                                                    SEXP cg_uno,
                                                    SEXP cg_ord,
                                                    SEXP cg_con,
                                                    SEXP mysd,
                                                    SEXP myopti,
                                                    SEXP myoptd,
                                                    SEXP bw,
                                                    SEXP hist_len,
                                                    SEXP penalty_mode,
                                                    SEXP penalty_mult,
                                                    SEXP glp_degree,
                                                    SEXP glp_bernstein,
                                                    SEXP glp_basis,
                                                    SEXP regtype,
                                                    SEXP cxkerlb,
                                                    SEXP cxkerub,
                                                    SEXP cykerlb,
                                                    SEXP cykerub,
                                                    const int eval_only)
{
  SEXP c_uno_r=R_NilValue, c_ord_r=R_NilValue, c_con_r=R_NilValue, u_uno_r=R_NilValue, u_ord_r=R_NilValue, u_con_r=R_NilValue;
  SEXP cg_uno_r=R_NilValue, cg_ord_r=R_NilValue, cg_con_r=R_NilValue, mysd_r=R_NilValue;
  SEXP myopti_i=R_NilValue, myoptd_r=R_NilValue, bw_r=R_NilValue, degree_i=R_NilValue, cxkerlb_r=R_NilValue, cxkerub_r=R_NilValue, cykerlb_r=R_NilValue, cykerub_r=R_NilValue;
  SEXP out=R_NilValue, out_names=R_NilValue;
  SEXP out_bw=R_NilValue, out_fval=R_NilValue, out_fval_hist=R_NilValue, out_eval_hist=R_NilValue;
  SEXP out_invalid_hist=R_NilValue, out_timing=R_NilValue, out_fast=R_NilValue;
  int hlen = asInteger(hist_len);
  int pmode = asInteger(penalty_mode);
  double pmult = asReal(penalty_mult);
  int bern = asInteger(glp_bernstein);
  int basis = asInteger(glp_basis);
  int ll_mode = asInteger(regtype);
  int ncon_x = 0;
  int ncon_y = 0;
  double * cxkerlb_p = NULL;
  double * cxkerub_p = NULL;
  double * cykerlb_p = NULL;
  double * cykerub_p = NULL;

  if(hlen < 1) hlen = 1;

  PROTECT(c_uno_r = coerceVector(c_uno, REALSXP));
  PROTECT(c_ord_r = coerceVector(c_ord, REALSXP));
  PROTECT(c_con_r = coerceVector(c_con, REALSXP));
  PROTECT(u_uno_r = coerceVector(u_uno, REALSXP));
  PROTECT(u_ord_r = coerceVector(u_ord, REALSXP));
  PROTECT(u_con_r = coerceVector(u_con, REALSXP));
  PROTECT(cg_uno_r = coerceVector(cg_uno, REALSXP));
  PROTECT(cg_ord_r = coerceVector(cg_ord, REALSXP));
  PROTECT(cg_con_r = coerceVector(cg_con, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(myoptd_r = coerceVector(myoptd, REALSXP));
  PROTECT(bw_r = coerceVector(bw, REALSXP));
  PROTECT(degree_i = coerceVector(glp_degree, INTSXP));
  PROTECT(cxkerlb_r = coerceVector(cxkerlb, REALSXP));
  PROTECT(cxkerub_r = coerceVector(cxkerub, REALSXP));
  PROTECT(cykerlb_r = coerceVector(cykerlb, REALSXP));
  PROTECT(cykerub_r = coerceVector(cykerub, REALSXP));

  ncon_x = (int)INTEGER(myopti_i)[CBW_UNCONI];
  ncon_y = (int)INTEGER(myopti_i)[CBW_CNCONI];
  resolve_bounds_or_default(cxkerlb_r, cxkerub_r, ncon_x, &cxkerlb_p, &cxkerub_p);
  resolve_bounds_or_default(cykerlb_r, cykerub_r, ncon_y, &cykerlb_p, &cykerub_p);

  PROTECT(out_bw = allocVector(REALSXP, XLENGTH(bw_r)));
  PROTECT(out_fval = allocVector(REALSXP, 3));
  PROTECT(out_fval_hist = allocVector(REALSXP, hlen));
  PROTECT(out_eval_hist = allocVector(REALSXP, hlen));
  PROTECT(out_invalid_hist = allocVector(REALSXP, hlen));
  PROTECT(out_timing = allocVector(REALSXP, 1));
  PROTECT(out_fast = allocVector(REALSXP, 1));

  memcpy(REAL(out_bw), REAL(bw_r), (size_t)XLENGTH(bw_r) * sizeof(double));
  np_distribution_conditional_bw(REAL(c_uno_r), REAL(c_ord_r), REAL(c_con_r),
                                 REAL(u_uno_r), REAL(u_ord_r), REAL(u_con_r),
                                 REAL(cg_uno_r), REAL(cg_ord_r), REAL(cg_con_r), REAL(mysd_r),
                                 INTEGER(myopti_i), REAL(myoptd_r), REAL(out_bw), REAL(out_fval),
                                 REAL(out_fval_hist), REAL(out_eval_hist), REAL(out_invalid_hist), REAL(out_timing),
                                 REAL(out_fast), &pmode, &pmult,
                                 INTEGER(degree_i), &bern, &basis, &ll_mode,
                                 cxkerlb_p, cxkerub_p, cykerlb_p, cykerub_p,
                                 eval_only);

  PROTECT(out = allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, out_bw);
  SET_VECTOR_ELT(out, 1, out_fval);
  SET_VECTOR_ELT(out, 2, out_fval_hist);
  SET_VECTOR_ELT(out, 3, out_eval_hist);
  SET_VECTOR_ELT(out, 4, out_invalid_hist);
  SET_VECTOR_ELT(out, 5, out_timing);
  SET_VECTOR_ELT(out, 6, out_fast);

  PROTECT(out_names = allocVector(STRSXP, 7));
  SET_STRING_ELT(out_names, 0, mkChar("bw"));
  SET_STRING_ELT(out_names, 1, mkChar("fval"));
  SET_STRING_ELT(out_names, 2, mkChar("fval.history"));
  SET_STRING_ELT(out_names, 3, mkChar("eval.history"));
  SET_STRING_ELT(out_names, 4, mkChar("invalid.history"));
  SET_STRING_ELT(out_names, 5, mkChar("timing"));
  SET_STRING_ELT(out_names, 6, mkChar("fast.history"));
  setAttrib(out, R_NamesSymbol, out_names);

  UNPROTECT(27);
  return out;
}

SEXP C_np_distribution_conditional_bw(SEXP c_uno,
                                      SEXP c_ord,
                                      SEXP c_con,
                                      SEXP u_uno,
                                      SEXP u_ord,
                                      SEXP u_con,
                                      SEXP cg_uno,
                                      SEXP cg_ord,
                                      SEXP cg_con,
                                      SEXP mysd,
                                      SEXP myopti,
                                      SEXP myoptd,
                                      SEXP bw,
                                      SEXP hist_len,
                                      SEXP penalty_mode,
                                      SEXP penalty_mult,
                                      SEXP glp_degree,
                                      SEXP glp_bernstein,
                                      SEXP glp_basis,
                                      SEXP regtype,
                                      SEXP cxkerlb,
                                      SEXP cxkerub,
                                      SEXP cykerlb,
                                      SEXP cykerub)
{
  return C_np_distribution_conditional_bw_common(c_uno, c_ord, c_con, u_uno, u_ord, u_con,
                                                 cg_uno, cg_ord, cg_con, mysd, myopti, myoptd,
                                                 bw, hist_len, penalty_mode, penalty_mult,
                                                 glp_degree, glp_bernstein, glp_basis, regtype,
                                                 cxkerlb, cxkerub, cykerlb, cykerub, 0);
}

SEXP C_np_distribution_conditional_bw_eval(SEXP c_uno,
                                           SEXP c_ord,
                                           SEXP c_con,
                                           SEXP u_uno,
                                           SEXP u_ord,
                                           SEXP u_con,
                                           SEXP cg_uno,
                                           SEXP cg_ord,
                                           SEXP cg_con,
                                           SEXP mysd,
                                           SEXP myopti,
                                           SEXP myoptd,
                                           SEXP bw,
                                           SEXP hist_len,
                                           SEXP penalty_mode,
                                           SEXP penalty_mult,
                                           SEXP glp_degree,
                                           SEXP glp_bernstein,
                                           SEXP glp_basis,
                                           SEXP regtype,
                                           SEXP cxkerlb,
                                           SEXP cxkerub,
                                           SEXP cykerlb,
                                           SEXP cykerub)
{
  return C_np_distribution_conditional_bw_common(c_uno, c_ord, c_con, u_uno, u_ord, u_con,
                                                 cg_uno, cg_ord, cg_con, mysd, myopti, myoptd,
                                                 bw, hist_len, penalty_mode, penalty_mult,
                                                 glp_degree, glp_bernstein, glp_basis, regtype,
                                                 cxkerlb, cxkerub, cykerlb, cykerub, 1);
}

SEXP C_np_density_conditional_bw(SEXP c_uno,
                                 SEXP c_ord,
                                 SEXP c_con,
                                 SEXP u_uno,
                                 SEXP u_ord,
                                 SEXP u_con,
                                 SEXP mysd,
                                 SEXP myopti,
                                 SEXP myoptd,
                                 SEXP bw,
                                 SEXP hist_len,
                                 SEXP penalty_mode,
                                 SEXP penalty_mult,
                                 SEXP glp_degree,
                                 SEXP glp_bernstein,
                                 SEXP glp_basis,
                                 SEXP regtype,
                                 SEXP cxkerlb,
                                 SEXP cxkerub,
                                 SEXP cykerlb,
                                 SEXP cykerub)
{
  return C_np_density_conditional_bw_common(c_uno, c_ord, c_con, u_uno, u_ord, u_con,
                                            mysd, myopti, myoptd, bw, hist_len,
                                            penalty_mode, penalty_mult,
                                            glp_degree, glp_bernstein, glp_basis, regtype,
                                            cxkerlb, cxkerub, cykerlb, cykerub, 0);
}

SEXP C_np_density_conditional_bw_eval(SEXP c_uno,
                                      SEXP c_ord,
                                      SEXP c_con,
                                      SEXP u_uno,
                                      SEXP u_ord,
                                      SEXP u_con,
                                      SEXP mysd,
                                      SEXP myopti,
                                      SEXP myoptd,
                                      SEXP bw,
                                      SEXP hist_len,
                                      SEXP penalty_mode,
                                      SEXP penalty_mult,
                                      SEXP glp_degree,
                                      SEXP glp_bernstein,
                                      SEXP glp_basis,
                                      SEXP regtype,
                                      SEXP cxkerlb,
                                      SEXP cxkerub,
                                      SEXP cykerlb,
                                      SEXP cykerub)
{
  return C_np_density_conditional_bw_common(c_uno, c_ord, c_con, u_uno, u_ord, u_con,
                                            mysd, myopti, myoptd, bw, hist_len,
                                            penalty_mode, penalty_mult,
                                            glp_degree, glp_bernstein, glp_basis, regtype,
                                            cxkerlb, cxkerub, cykerlb, cykerub, 1);
}

SEXP C_np_shadow_cv_xweights_conditional(SEXP tyuno,
                                         SEXP tyord,
                                         SEXP tycon,
                                         SEXP txuno,
                                         SEXP txord,
                                         SEXP txcon,
                                         SEXP rbw,
                                         SEXP bwtype,
                                         SEXP kernel_x,
                                         SEXP kernel_xu,
                                         SEXP kernel_xo,
                                         SEXP use_tree,
                                         SEXP regtype,
                                         SEXP glp_degree,
                                         SEXP glp_bernstein,
                                         SEXP glp_basis,
                                         SEXP row_index)
{
  SEXP tycon_r = R_NilValue;
  SEXP txuno_r = R_NilValue, txord_r = R_NilValue, txcon_r = R_NilValue;
  SEXP rbw_r = R_NilValue, degree_i = R_NilValue;
  SEXP out = R_NilValue, out_names = R_NilValue, out_streamed = R_NilValue;
  int nrow_yuno = 0, ncol_yuno = 0, nrow_yord = 0, ncol_yord = 0, nrow_ycon = 0, ncol_ycon = 0;
  int nrow_xuno = 0, ncol_xuno = 0, nrow_xord = 0, ncol_xord = 0, nrow_xcon = 0, ncol_xcon = 0;
  int num_obs = 0;
  int tree_flag = asLogical(use_tree);
  int int_large_sf_save = int_LARGE_SF;
  double nconfac_save = nconfac_extern;
  double ncatfac_save = ncatfac_extern;
  double *vector_continuous_stddev_save = vector_continuous_stddev_extern;
  double *shadow_continuous_stddev = NULL;
  double *streamed_row = NULL;
  int row_idx = asInteger(row_index) - 1;
  int i;

  tycon_r = PROTECT(coerceVector(tycon, REALSXP));
  txuno_r = PROTECT(coerceVector(txuno, REALSXP));
  txord_r = PROTECT(coerceVector(txord, REALSXP));
  txcon_r = PROTECT(coerceVector(txcon, REALSXP));
  rbw_r = PROTECT(coerceVector(rbw, REALSXP));
  degree_i = PROTECT(coerceVector(glp_degree, INTSXP));

  np_shadow_matrix_dims(tyuno, &nrow_yuno, &ncol_yuno);
  np_shadow_matrix_dims(tyord, &nrow_yord, &ncol_yord);
  np_shadow_matrix_dims(tycon, &nrow_ycon, &ncol_ycon);
  np_shadow_matrix_dims(txuno, &nrow_xuno, &ncol_xuno);
  np_shadow_matrix_dims(txord, &nrow_xord, &ncol_xord);
  np_shadow_matrix_dims(txcon, &nrow_xcon, &ncol_xcon);

  num_obs = MAX(MAX(nrow_yuno, nrow_yord), MAX(MAX(nrow_ycon, nrow_xuno), MAX(nrow_xord, nrow_xcon)));
  if(num_obs <= 0)
    error("C_np_shadow_cv_xweights_conditional: zero-row inputs are not supported");
  if((nrow_yuno > 0 && nrow_yuno != num_obs) || (nrow_yord > 0 && nrow_yord != num_obs) ||
     (nrow_ycon > 0 && nrow_ycon != num_obs) || (nrow_xuno > 0 && nrow_xuno != num_obs) ||
     (nrow_xord > 0 && nrow_xord != num_obs) || (nrow_xcon > 0 && nrow_xcon != num_obs))
    error("C_np_shadow_cv_xweights_conditional: all inputs must share the same row count");
  if((row_idx < 0) || (row_idx >= num_obs))
    error("C_np_shadow_cv_xweights_conditional: row_index out of range");
  if((asInteger(bwtype) != BW_FIXED) &&
     (asInteger(bwtype) != BW_GEN_NN))
    error("C_np_shadow_cv_xweights_conditional: fixed/generalized-nn bandwidths only");

  np_shadow_state_active = 1;
  num_obs_train_extern = num_obs_eval_extern = num_obs;
  num_var_unordered_extern = ncol_yuno;
  num_var_ordered_extern = ncol_yord;
  num_var_continuous_extern = ncol_ycon;
  num_reg_unordered_extern = ncol_xuno;
  num_reg_ordered_extern = ncol_xord;
  num_reg_continuous_extern = ncol_xcon;
  int_LARGE_SF = 1;
  nconfac_extern = 1.0;
  ncatfac_extern = 1.0;

  KERNEL_den_extern = CK_GAUSS2;
  KERNEL_den_unordered_extern = UKERNEL_UAA;
  KERNEL_den_ordered_extern = 0;
  KERNEL_reg_extern = asInteger(kernel_x);
  KERNEL_reg_unordered_extern = asInteger(kernel_xu);
  KERNEL_reg_ordered_extern = asInteger(kernel_xo);
  BANDWIDTH_den_extern = asInteger(bwtype);

  matrix_X_unordered_train_extern = alloc_matd(num_obs, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs, num_reg_continuous_extern);
  np_shadow_fill_matrix(matrix_X_unordered_train_extern, REAL(txuno_r), num_obs, num_reg_unordered_extern);
  np_shadow_fill_matrix(matrix_X_ordered_train_extern, REAL(txord_r), num_obs, num_reg_ordered_extern);
  np_shadow_fill_matrix(matrix_X_continuous_train_extern, REAL(txcon_r), num_obs, num_reg_continuous_extern);

  if((num_reg_continuous_extern + num_var_continuous_extern) > 0){
    shadow_continuous_stddev =
      (double *)malloc((size_t)(num_reg_continuous_extern + num_var_continuous_extern) * sizeof(double));
    if(shadow_continuous_stddev == NULL)
      error("C_np_shadow_cv_xweights_conditional: stddev allocation failed");
    for(i = 0; i < num_reg_continuous_extern; i++)
      shadow_continuous_stddev[i] =
        standerrd(num_obs, matrix_X_continuous_train_extern[i]);
    for(i = 0; i < num_var_continuous_extern; i++)
      shadow_continuous_stddev[num_reg_continuous_extern + i] =
        standerrd(num_obs, REAL(tycon_r) + i*num_obs);
    vector_continuous_stddev_extern = shadow_continuous_stddev;
  } else {
    vector_continuous_stddev_extern = NULL;
  }

  if(tree_flag){
    int_TREE_X = (num_reg_continuous_extern > 0) ? NP_TREE_TRUE : NP_TREE_FALSE;
    int_TREE_Y = int_TREE_XY = NP_TREE_FALSE;
    if(int_TREE_X == NP_TREE_TRUE){
      ipt_extern_X = (int *)malloc((size_t)num_obs*sizeof(int));
      ipt_lookup_extern_X = (int *)malloc((size_t)num_obs*sizeof(int));
      if((ipt_extern_X == NULL) || (ipt_lookup_extern_X == NULL))
        error("C_np_shadow_cv_xweights_conditional: x-tree allocation failed");
      for(i = 0; i < num_obs; i++) ipt_extern_X[i] = i;
      build_kdtree(matrix_X_continuous_train_extern, num_obs, num_reg_continuous_extern,
                   4*num_reg_continuous_extern, ipt_extern_X, &kdt_extern_X);
      for(i = 0; i < num_obs; i++) ipt_lookup_extern_X[ipt_extern_X[i]] = i;
      for(int j = 0; j < num_reg_unordered_extern; j++)
        for(i = 0; i < num_obs; i++)
          matrix_X_unordered_train_extern[j][i] = REAL(txuno_r)[(size_t)j * (size_t)num_obs + (size_t)ipt_extern_X[i]];
      for(int j = 0; j < num_reg_ordered_extern; j++)
        for(i = 0; i < num_obs; i++)
          matrix_X_ordered_train_extern[j][i] = REAL(txord_r)[(size_t)j * (size_t)num_obs + (size_t)ipt_extern_X[i]];
      for(int j = 0; j < num_reg_continuous_extern; j++)
        for(i = 0; i < num_obs; i++)
          matrix_X_continuous_train_extern[j][i] = REAL(txcon_r)[(size_t)j * (size_t)num_obs + (size_t)ipt_extern_X[i]];
    } else {
      ipt_extern_X = NULL;
      ipt_lookup_extern_X = NULL;
      kdt_extern_X = NULL;
    }
  } else {
    int_TREE_X = int_TREE_Y = int_TREE_XY = NP_TREE_FALSE;
    ipt_extern_X = ipt_extern_Y = ipt_extern_XY = NULL;
    ipt_lookup_extern_X = ipt_lookup_extern_Y = ipt_lookup_extern_XY = NULL;
    kdt_extern_X = kdt_extern_Y = kdt_extern_XY = NULL;
  }

  num_categories_extern_X = alloc_vecu(num_reg_unordered_extern + num_reg_ordered_extern);
  matrix_categorical_vals_extern_X = alloc_matd(num_obs, num_reg_unordered_extern + num_reg_ordered_extern);
  determine_categorical_vals(num_obs,
                             0,
                             0,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             NULL,
                             NULL,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern_X,
                             matrix_categorical_vals_extern_X);

  int_ll_extern = asInteger(regtype);
  if((int_ll_extern == LL_LP) && (num_reg_continuous_extern > 0)){
    if((int)XLENGTH(degree_i) != num_reg_continuous_extern)
      error("C_np_shadow_cv_xweights_conditional: glp_degree length mismatch");
    vector_glp_degree_extern = INTEGER(degree_i);
    int_glp_bernstein_extern = asInteger(glp_bernstein);
    int_glp_basis_extern = asInteger(glp_basis);
  } else {
    vector_glp_degree_extern = NULL;
    int_glp_bernstein_extern = 0;
    int_glp_basis_extern = 1;
  }

  streamed_row = (double *)malloc((size_t)num_obs*sizeof(double));
  if(streamed_row == NULL)
    error("C_np_shadow_cv_xweights_conditional: weight allocation failed");

  if(np_shadow_proof_conditional_x_weight_row_stream(REAL(rbw_r), row_idx, streamed_row) != 0)
    error("C_np_shadow_cv_xweights_conditional: streamed row helper failed");

  out_streamed = PROTECT(allocVector(REALSXP, num_obs));
  for(i = 0; i < num_obs; i++)
    REAL(out_streamed)[i] = streamed_row[i];
  out = PROTECT(allocVector(VECSXP, 1));
  SET_VECTOR_ELT(out, 0, out_streamed);
  out_names = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(out_names, 0, mkChar("streamed"));
  setAttrib(out, R_NamesSymbol, out_names);

  if(kdt_extern_X != NULL) free_kdtree(&kdt_extern_X);
  safe_free(ipt_extern_X); safe_free(ipt_lookup_extern_X);
  if(matrix_X_unordered_train_extern != NULL) free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  if(matrix_X_ordered_train_extern != NULL) free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  if(matrix_X_continuous_train_extern != NULL) free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  if(num_categories_extern_X != NULL) safe_free(num_categories_extern_X);
  if(matrix_categorical_vals_extern_X != NULL) free_mat(matrix_categorical_vals_extern_X, num_reg_unordered_extern + num_reg_ordered_extern);
  safe_free(streamed_row);
  safe_free(shadow_continuous_stddev);
  matrix_X_unordered_train_extern = matrix_X_ordered_train_extern = matrix_X_continuous_train_extern = NULL;
  num_categories_extern_X = NULL;
  matrix_categorical_vals_extern_X = NULL;
  ipt_extern_X = NULL;
  ipt_lookup_extern_X = NULL;
  int_ll_extern = LL_LC;
  vector_glp_degree_extern = NULL;
  int_glp_bernstein_extern = 0;
  int_glp_basis_extern = 1;
  int_TREE_X = int_TREE_Y = int_TREE_XY = NP_TREE_FALSE;
  int_LARGE_SF = int_large_sf_save;
  nconfac_extern = nconfac_save;
  ncatfac_extern = ncatfac_save;
  vector_continuous_stddev_extern = vector_continuous_stddev_save;
  np_glp_cv_clear_extern();
  np_shadow_state_active = 0;

  UNPROTECT(9);
  return out;
}

SEXP C_np_regression_lp_apply_conditional(SEXP txuno,
                                          SEXP txord,
                                          SEXP txcon,
                                          SEXP exuno,
                                          SEXP exord,
                                          SEXP excon,
                                          SEXP rhs,
                                          SEXP rbw,
                                          SEXP bwtype,
                                          SEXP kernel_x,
                                          SEXP kernel_xu,
                                          SEXP kernel_xo,
                                          SEXP use_tree,
                                          SEXP glp_degree,
                                          SEXP glp_gradient_order,
                                          SEXP glp_bernstein,
                                          SEXP glp_basis)
{
  SEXP txuno_r = R_NilValue, txord_r = R_NilValue, txcon_r = R_NilValue;
  SEXP exuno_r = R_NilValue, exord_r = R_NilValue, excon_r = R_NilValue;
  SEXP rhs_r = R_NilValue, rbw_r = R_NilValue, degree_i = R_NilValue, grad_i = R_NilValue, out = R_NilValue;
  int nrow_txuno = 0, ncol_txuno = 0, nrow_txord = 0, ncol_txord = 0, nrow_txcon = 0, ncol_txcon = 0;
  int nrow_exuno = 0, ncol_exuno = 0, nrow_exord = 0, ncol_exord = 0, nrow_excon = 0, ncol_excon = 0;
  int nrow_rhs = 0, ncol_rhs = 0;
  int num_obs_train = 0, num_obs_eval = 0;
  int tree_flag = asLogical(use_tree);
  int int_large_sf_save = int_LARGE_SF;
  double nconfac_save = nconfac_extern;
  double ncatfac_save = ncatfac_extern;
  double *vector_continuous_stddev_save = vector_continuous_stddev_extern;
  double *shadow_continuous_stddev = NULL;
  double **rhs_cols = NULL;
  double *rhs_sorted = NULL;
  int helper_status = 1;
  int i;

  txuno_r = PROTECT(coerceVector(txuno, REALSXP));
  txord_r = PROTECT(coerceVector(txord, REALSXP));
  txcon_r = PROTECT(coerceVector(txcon, REALSXP));
  exuno_r = PROTECT(coerceVector(exuno, REALSXP));
  exord_r = PROTECT(coerceVector(exord, REALSXP));
  excon_r = PROTECT(coerceVector(excon, REALSXP));
  rhs_r = PROTECT(coerceVector(rhs, REALSXP));
  rbw_r = PROTECT(coerceVector(rbw, REALSXP));
  degree_i = PROTECT(coerceVector(glp_degree, INTSXP));
  grad_i = PROTECT(coerceVector(glp_gradient_order, INTSXP));

  np_shadow_matrix_dims(txuno, &nrow_txuno, &ncol_txuno);
  np_shadow_matrix_dims(txord, &nrow_txord, &ncol_txord);
  np_shadow_matrix_dims(txcon, &nrow_txcon, &ncol_txcon);
  np_shadow_matrix_dims(exuno, &nrow_exuno, &ncol_exuno);
  np_shadow_matrix_dims(exord, &nrow_exord, &ncol_exord);
  np_shadow_matrix_dims(excon, &nrow_excon, &ncol_excon);
  np_shadow_matrix_dims(rhs, &nrow_rhs, &ncol_rhs);

  num_obs_train = MAX(nrow_txuno, MAX(nrow_txord, nrow_txcon));
  num_obs_eval = MAX(nrow_exuno, MAX(nrow_exord, nrow_excon));
  if((num_obs_train <= 0) || (num_obs_eval <= 0))
    error("C_np_regression_lp_apply_conditional: zero-row inputs are not supported");
  if((nrow_txuno > 0 && nrow_txuno != num_obs_train) ||
     (nrow_txord > 0 && nrow_txord != num_obs_train) ||
     (nrow_txcon > 0 && nrow_txcon != num_obs_train))
    error("C_np_regression_lp_apply_conditional: training inputs must share the same row count");
  if((nrow_exuno > 0 && nrow_exuno != num_obs_eval) ||
     (nrow_exord > 0 && nrow_exord != num_obs_eval) ||
     (nrow_excon > 0 && nrow_excon != num_obs_eval))
    error("C_np_regression_lp_apply_conditional: evaluation inputs must share the same row count");
  if((nrow_rhs <= 0) || (ncol_rhs <= 0))
    error("C_np_regression_lp_apply_conditional: rhs must be a non-empty numeric matrix");
  if(nrow_rhs != num_obs_train)
    error("C_np_regression_lp_apply_conditional: rhs row count must match training row count");
  if((asInteger(bwtype) != BW_FIXED) &&
     (asInteger(bwtype) != BW_GEN_NN))
    error("C_np_regression_lp_apply_conditional: fixed/generalized-nn bandwidths only");

  np_shadow_state_active = 1;
  num_obs_train_extern = num_obs_train;
  num_obs_eval_extern = num_obs_eval;
  num_var_unordered_extern = 0;
  num_var_ordered_extern = 0;
  num_var_continuous_extern = 0;
  num_reg_unordered_extern = ncol_txuno;
  num_reg_ordered_extern = ncol_txord;
  num_reg_continuous_extern = ncol_txcon;
  int_LARGE_SF = 1;
  nconfac_extern = 1.0;
  ncatfac_extern = 1.0;

  KERNEL_den_extern = CK_GAUSS2;
  KERNEL_den_unordered_extern = UKERNEL_UAA;
  KERNEL_den_ordered_extern = 0;
  KERNEL_reg_extern = asInteger(kernel_x);
  KERNEL_reg_unordered_extern = asInteger(kernel_xu);
  KERNEL_reg_ordered_extern = asInteger(kernel_xo);
  BANDWIDTH_den_extern = asInteger(bwtype);

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train, num_reg_continuous_extern);
  matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval, num_reg_unordered_extern);
  matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval, num_reg_ordered_extern);
  matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval, num_reg_continuous_extern);
  np_shadow_fill_matrix(matrix_X_unordered_train_extern, REAL(txuno_r), num_obs_train, num_reg_unordered_extern);
  np_shadow_fill_matrix(matrix_X_ordered_train_extern, REAL(txord_r), num_obs_train, num_reg_ordered_extern);
  np_shadow_fill_matrix(matrix_X_continuous_train_extern, REAL(txcon_r), num_obs_train, num_reg_continuous_extern);
  np_shadow_fill_matrix(matrix_X_unordered_eval_extern, REAL(exuno_r), num_obs_eval, num_reg_unordered_extern);
  np_shadow_fill_matrix(matrix_X_ordered_eval_extern, REAL(exord_r), num_obs_eval, num_reg_ordered_extern);
  np_shadow_fill_matrix(matrix_X_continuous_eval_extern, REAL(excon_r), num_obs_eval, num_reg_continuous_extern);

  if(num_reg_continuous_extern > 0){
    shadow_continuous_stddev =
      (double *)malloc((size_t)num_reg_continuous_extern*sizeof(double));
    if(shadow_continuous_stddev == NULL)
      error("C_np_regression_lp_apply_conditional: stddev allocation failed");
    for(i = 0; i < num_reg_continuous_extern; i++)
      shadow_continuous_stddev[i] =
        standerrd(num_obs_train, matrix_X_continuous_train_extern[i]);
    vector_continuous_stddev_extern = shadow_continuous_stddev;
  } else {
    vector_continuous_stddev_extern = NULL;
  }

  if(tree_flag){
    int_TREE_X = (num_reg_continuous_extern > 0) ? NP_TREE_TRUE : NP_TREE_FALSE;
    int_TREE_Y = int_TREE_XY = NP_TREE_FALSE;
    if(int_TREE_X == NP_TREE_TRUE){
      ipt_extern_X = (int *)malloc((size_t)num_obs_train*sizeof(int));
      ipt_lookup_extern_X = (int *)malloc((size_t)num_obs_train*sizeof(int));
      if((ipt_extern_X == NULL) || (ipt_lookup_extern_X == NULL))
        error("C_np_regression_lp_apply_conditional: x-tree allocation failed");
      for(i = 0; i < num_obs_train; i++) ipt_extern_X[i] = i;
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train, num_reg_continuous_extern,
                   4*num_reg_continuous_extern, ipt_extern_X, &kdt_extern_X);
      for(i = 0; i < num_obs_train; i++) ipt_lookup_extern_X[ipt_extern_X[i]] = i;
      for(int j = 0; j < num_reg_unordered_extern; j++)
        for(i = 0; i < num_obs_train; i++)
          matrix_X_unordered_train_extern[j][i] = REAL(txuno_r)[(size_t)j * (size_t)num_obs_train + (size_t)ipt_extern_X[i]];
      for(int j = 0; j < num_reg_ordered_extern; j++)
        for(i = 0; i < num_obs_train; i++)
          matrix_X_ordered_train_extern[j][i] = REAL(txord_r)[(size_t)j * (size_t)num_obs_train + (size_t)ipt_extern_X[i]];
      for(int j = 0; j < num_reg_continuous_extern; j++)
        for(i = 0; i < num_obs_train; i++)
          matrix_X_continuous_train_extern[j][i] = REAL(txcon_r)[(size_t)j * (size_t)num_obs_train + (size_t)ipt_extern_X[i]];
    } else {
      ipt_extern_X = NULL;
      ipt_lookup_extern_X = NULL;
      kdt_extern_X = NULL;
    }
  } else {
    int_TREE_X = int_TREE_Y = int_TREE_XY = NP_TREE_FALSE;
    ipt_extern_X = ipt_extern_Y = ipt_extern_XY = NULL;
    ipt_lookup_extern_X = ipt_lookup_extern_Y = ipt_lookup_extern_XY = NULL;
    kdt_extern_X = kdt_extern_Y = kdt_extern_XY = NULL;
  }

  num_categories_extern_X = alloc_vecu(num_reg_unordered_extern + num_reg_ordered_extern);
  matrix_categorical_vals_extern_X = alloc_matd(num_obs_train, num_reg_unordered_extern + num_reg_ordered_extern);
  determine_categorical_vals(num_obs_train,
                             0,
                             0,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             NULL,
                             NULL,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern_X,
                             matrix_categorical_vals_extern_X);

  int_ll_extern = LL_LP;
  if((int)XLENGTH(degree_i) != num_reg_continuous_extern)
    error("C_np_regression_lp_apply_conditional: glp_degree length mismatch");
  if((XLENGTH(grad_i) != 0) && ((int)XLENGTH(grad_i) != num_reg_continuous_extern))
    error("C_np_regression_lp_apply_conditional: glp_gradient_order length mismatch");
  vector_glp_degree_extern = INTEGER(degree_i);
  vector_glp_gradient_order_extern = (XLENGTH(grad_i) > 0) ? INTEGER(grad_i) : NULL;
  int_glp_bernstein_extern = asInteger(glp_bernstein);
  int_glp_basis_extern = asInteger(glp_basis);

  rhs_cols = (double **)malloc((size_t)ncol_rhs*sizeof(double *));
  if(rhs_cols == NULL)
    error("C_np_regression_lp_apply_conditional: rhs column allocation failed");

  if(int_TREE_X == NP_TREE_TRUE){
    rhs_sorted = (double *)malloc((size_t)num_obs_train*(size_t)ncol_rhs*sizeof(double));
    if(rhs_sorted == NULL)
      error("C_np_regression_lp_apply_conditional: rhs reorder allocation failed");
    for(int j = 0; j < ncol_rhs; j++){
      double *rhs_col = rhs_sorted + (size_t)j*(size_t)num_obs_train;
      for(i = 0; i < num_obs_train; i++)
        rhs_col[i] = REAL(rhs_r)[(size_t)ipt_extern_X[i] + (size_t)j*(size_t)num_obs_train];
      rhs_cols[j] = rhs_col;
    }
  } else {
    for(i = 0; i < ncol_rhs; i++)
      rhs_cols[i] = REAL(rhs_r) + ((size_t)i * (size_t)nrow_rhs);
  }

  out = PROTECT(allocMatrix(REALSXP, num_obs_eval, ncol_rhs));
  helper_status = np_regression_lp_apply_matrix(REAL(rbw_r), rhs_cols, ncol_rhs, REAL(out));

  if(kdt_extern_X != NULL) free_kdtree(&kdt_extern_X);
  safe_free(ipt_extern_X); safe_free(ipt_lookup_extern_X);
  if(matrix_X_unordered_train_extern != NULL) free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  if(matrix_X_ordered_train_extern != NULL) free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  if(matrix_X_continuous_train_extern != NULL) free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  if(matrix_X_unordered_eval_extern != NULL) free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
  if(matrix_X_ordered_eval_extern != NULL) free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
  if(matrix_X_continuous_eval_extern != NULL) free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  if(num_categories_extern_X != NULL) safe_free(num_categories_extern_X);
  if(matrix_categorical_vals_extern_X != NULL) free_mat(matrix_categorical_vals_extern_X, num_reg_unordered_extern + num_reg_ordered_extern);
  safe_free(shadow_continuous_stddev);
  safe_free(rhs_sorted);
  safe_free(rhs_cols);
  matrix_X_unordered_train_extern = matrix_X_ordered_train_extern = matrix_X_continuous_train_extern = NULL;
  matrix_X_unordered_eval_extern = matrix_X_ordered_eval_extern = matrix_X_continuous_eval_extern = NULL;
  num_categories_extern_X = NULL;
  matrix_categorical_vals_extern_X = NULL;
  ipt_extern_X = NULL;
  ipt_lookup_extern_X = NULL;
  int_ll_extern = LL_LC;
  vector_glp_degree_extern = NULL;
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = 0;
  int_glp_basis_extern = 1;
  int_TREE_X = int_TREE_Y = int_TREE_XY = NP_TREE_FALSE;
  vector_continuous_stddev_extern = vector_continuous_stddev_save;
  int_LARGE_SF = int_large_sf_save;
  nconfac_extern = nconfac_save;
  ncatfac_extern = ncatfac_save;
  np_glp_cv_clear_extern();
  np_shadow_state_active = 0;

  if(helper_status != 0)
    error("C_np_regression_lp_apply_conditional: lp apply helper failed");

  UNPROTECT(11);
  return out;
}

SEXP C_np_kernelsum(SEXP tuno,
                    SEXP tord,
                    SEXP tcon,
                    SEXP ty,
                    SEXP weights,
                    SEXP euno,
                    SEXP eord,
                    SEXP econ,
                    SEXP bw,
                    SEXP mcv,
                    SEXP padnum,
                    SEXP op,
                    SEXP myopti,
                    SEXP kpow,
                    SEXP ksum_len,
                    SEXP pksum_len,
                    SEXP kw_len,
                    SEXP ckerlb,
                    SEXP ckerub)
{
  SEXP tuno_r=R_NilValue, tord_r=R_NilValue, tcon_r=R_NilValue, ty_r=R_NilValue, weights_r=R_NilValue;
  SEXP euno_r=R_NilValue, eord_r=R_NilValue, econ_r=R_NilValue, bw_r=R_NilValue, mcv_r=R_NilValue;
  SEXP padnum_r=R_NilValue, op_i=R_NilValue, myopti_i=R_NilValue, kpow_r=R_NilValue, ckerlb_r=R_NilValue, ckerub_r=R_NilValue;
  SEXP out=R_NilValue, out_names=R_NilValue, out_ksum=R_NilValue, out_pksum=R_NilValue, out_kw=R_NilValue, out_pkw=R_NilValue;
  int n_ksum = asInteger(ksum_len);
  int n_pksum = asInteger(pksum_len);
  int n_kw = asInteger(kw_len);
  int ncon = 0, nuno = 0, nord = 0, p_operator = 0, do_score = 0, do_ocg = 0, return_kernel_weights = 0, p_nvar = 0;
  R_xlen_t n_pkw = 0;
  int * myopti_p = NULL;
  double * ckerlb_p = NULL;
  double * ckerub_p = NULL;

  if(n_ksum < 0) n_ksum = 0;
  if(n_pksum < 0) n_pksum = 0;
  if(n_kw < 0) n_kw = 0;

  PROTECT(tuno_r = coerceVector(tuno, REALSXP));
  PROTECT(tord_r = coerceVector(tord, REALSXP));
  PROTECT(tcon_r = coerceVector(tcon, REALSXP));
  PROTECT(ty_r = coerceVector(ty, REALSXP));
  PROTECT(weights_r = coerceVector(weights, REALSXP));
  PROTECT(euno_r = coerceVector(euno, REALSXP));
  PROTECT(eord_r = coerceVector(eord, REALSXP));
  PROTECT(econ_r = coerceVector(econ, REALSXP));
  PROTECT(bw_r = coerceVector(bw, REALSXP));
  PROTECT(mcv_r = coerceVector(mcv, REALSXP));
  PROTECT(padnum_r = coerceVector(padnum, REALSXP));
  PROTECT(op_i = coerceVector(op, INTSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(kpow_r = coerceVector(kpow, REALSXP));
  PROTECT(ckerlb_r = coerceVector(ckerlb, REALSXP));
  PROTECT(ckerub_r = coerceVector(ckerub, REALSXP));

  ncon = (int)INTEGER(myopti_i)[KWS_NCONI];
  nuno = (int)INTEGER(myopti_i)[KWS_NUNOI];
  nord = (int)INTEGER(myopti_i)[KWS_NORDI];
  p_operator = (int)INTEGER(myopti_i)[KWS_POPI];
  do_score = (int)INTEGER(myopti_i)[KWS_PSCOREI];
  do_ocg = (int)INTEGER(myopti_i)[KWS_POCGI];
  return_kernel_weights = (int)INTEGER(myopti_i)[KWS_RKWI];
  p_nvar = ((p_operator != OP_NOOP) ? ncon : 0) + ((do_score || do_ocg) ? (nuno + nord) : 0);
  n_pkw = (return_kernel_weights && (p_nvar > 0)) ? ((R_xlen_t)n_kw * (R_xlen_t)p_nvar) : 0;
  myopti_p = INTEGER(myopti_i);
  if(XLENGTH(myopti_i) <= KWS_SPARI){
    int i;
    int * myopti_tmp = (int *)R_alloc((size_t)(KWS_SPARI + 1), sizeof(int));
    for(i = 0; i <= KWS_SPARI; i++)
      myopti_tmp[i] = 0;
    for(i = 0; i < XLENGTH(myopti_i); i++)
      myopti_tmp[i] = myopti_p[i];
    myopti_p = myopti_tmp;
  }
  resolve_bounds_or_default(ckerlb_r, ckerub_r, ncon, &ckerlb_p, &ckerub_p);

  PROTECT(out_ksum = allocVector(REALSXP, n_ksum));
  PROTECT(out_pksum = allocVector(REALSXP, n_pksum));
  PROTECT(out_kw = allocVector(REALSXP, n_kw));
  PROTECT(out_pkw = allocVector(REALSXP, n_pkw));

  np_kernelsum(REAL(tuno_r), REAL(tord_r), REAL(tcon_r),
               REAL(ty_r), REAL(weights_r),
               REAL(euno_r), REAL(eord_r), REAL(econ_r),
               REAL(bw_r),
               REAL(mcv_r), REAL(padnum_r),
               INTEGER(op_i), myopti_p, REAL(kpow_r),
               REAL(out_ksum), REAL(out_pksum), REAL(out_kw), REAL(out_pkw),
               ckerlb_p, ckerub_p);

  PROTECT(out = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, out_ksum);
  SET_VECTOR_ELT(out, 1, out_pksum);
  SET_VECTOR_ELT(out, 2, out_kw);
  SET_VECTOR_ELT(out, 3, out_pkw);

  PROTECT(out_names = allocVector(STRSXP, 4));
  SET_STRING_ELT(out_names, 0, mkChar("ksum"));
  SET_STRING_ELT(out_names, 1, mkChar("p.ksum"));
  SET_STRING_ELT(out_names, 2, mkChar("kernel.weights"));
  SET_STRING_ELT(out_names, 3, mkChar("p.kernel.weights"));
  setAttrib(out, R_NamesSymbol, out_names);

  UNPROTECT(22);
  return out;
}

SEXP C_np_quantile_conditional(SEXP tc_con,
                               SEXP tu_uno,
                               SEXP tu_ord,
                               SEXP tu_con,
                               SEXP eu_uno,
                               SEXP eu_ord,
                               SEXP eu_con,
                               SEXP quantile,
                               SEXP mybw,
                               SEXP mcv,
                               SEXP padnum,
                               SEXP nconfac,
                               SEXP ncatfac,
                               SEXP mysd,
                               SEXP myopti,
                               SEXP myoptd,
                               SEXP enrow,
                               SEXP xndim,
                               SEXP gradients)
{
  SEXP tc_con_r=R_NilValue, tu_uno_r=R_NilValue, tu_ord_r=R_NilValue, tu_con_r=R_NilValue;
  SEXP eu_uno_r=R_NilValue, eu_ord_r=R_NilValue, eu_con_r=R_NilValue, quantile_r=R_NilValue;
  SEXP mybw_r=R_NilValue, mcv_r=R_NilValue, padnum_r=R_NilValue, nconfac_r=R_NilValue, ncatfac_r=R_NilValue, mysd_r=R_NilValue;
  SEXP myopti_i=R_NilValue, myoptd_r=R_NilValue;
  SEXP out=R_NilValue, out_names=R_NilValue, out_yq=R_NilValue, out_yqerr=R_NilValue, out_yg=R_NilValue;
  int en = asInteger(enrow);
  int xd = asInteger(xndim);
  int do_grad = asLogical(gradients);
  R_xlen_t gsize;

  if(en < 0) en = 0;
  if(xd < 0) xd = 0;
  if(do_grad == NA_LOGICAL) do_grad = 0;
  gsize = (R_xlen_t)en * (R_xlen_t)xd * (do_grad ? 1 : 0);

  PROTECT(tc_con_r = coerceVector(tc_con, REALSXP));
  PROTECT(tu_uno_r = coerceVector(tu_uno, REALSXP));
  PROTECT(tu_ord_r = coerceVector(tu_ord, REALSXP));
  PROTECT(tu_con_r = coerceVector(tu_con, REALSXP));
  PROTECT(eu_uno_r = coerceVector(eu_uno, REALSXP));
  PROTECT(eu_ord_r = coerceVector(eu_ord, REALSXP));
  PROTECT(eu_con_r = coerceVector(eu_con, REALSXP));
  PROTECT(quantile_r = coerceVector(quantile, REALSXP));
  PROTECT(mybw_r = coerceVector(mybw, REALSXP));
  PROTECT(mcv_r = coerceVector(mcv, REALSXP));
  PROTECT(padnum_r = coerceVector(padnum, REALSXP));
  PROTECT(nconfac_r = coerceVector(nconfac, REALSXP));
  PROTECT(ncatfac_r = coerceVector(ncatfac, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(myoptd_r = coerceVector(myoptd, REALSXP));

  PROTECT(out_yq = allocVector(REALSXP, en));
  PROTECT(out_yqerr = allocVector(REALSXP, en));
  PROTECT(out_yg = allocVector(REALSXP, gsize));

  np_quantile_conditional(REAL(tc_con_r),
                          REAL(tu_uno_r), REAL(tu_ord_r), REAL(tu_con_r),
                          REAL(eu_uno_r), REAL(eu_ord_r), REAL(eu_con_r),
                          REAL(quantile_r),
                          REAL(mybw_r),
                          REAL(mcv_r), REAL(padnum_r),
                          REAL(nconfac_r), REAL(ncatfac_r), REAL(mysd_r),
                          INTEGER(myopti_i), REAL(myoptd_r),
                          REAL(out_yq), REAL(out_yqerr), REAL(out_yg));

  PROTECT(out = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, out_yq);
  SET_VECTOR_ELT(out, 1, out_yqerr);
  SET_VECTOR_ELT(out, 2, out_yg);

  PROTECT(out_names = allocVector(STRSXP, 3));
  SET_STRING_ELT(out_names, 0, mkChar("yq"));
  SET_STRING_ELT(out_names, 1, mkChar("yqerr"));
  SET_STRING_ELT(out_names, 2, mkChar("yqgrad"));
  setAttrib(out, R_NamesSymbol, out_names);

  UNPROTECT(21);
  return out;
}

void np_mpi_init(int * mpi_status){
#ifdef MPI2 
  MPI_Comm_rank(comm[1], &my_rank);
  MPI_Comm_size(comm[1], &iNum_Processors);
  mpi_status[MPI_RANKI] = my_rank;
  mpi_status[MPI_NUMPI] = iNum_Processors;
#else
  mpi_status[MPI_RANKI] = -1;
  mpi_status[MPI_NUMPI] = -1;
#endif
}

SEXP C_np_mpi_init(void)
{
  SEXP ans = PROTECT(allocVector(INTSXP, 2));
  np_mpi_init(INTEGER(ans));
  UNPROTECT(1);
  return ans;
}


void np_density_bw(double * myuno, double * myord, double * mycon, 
                   double * mysd, int * myopti, double * myoptd, double * myans, double * fval,
                   double * objective_function_values, double * objective_function_evals,
                   double * objective_function_invalid, double * timing,
                   double * objective_function_fast,
                   int * penalty_mode, double * penalty_mult,
                   double * ckerlb, double * ckerub){
  int_nn_k_min_extern = 1;
  /* Likelihood bandwidth selection for density estimation */

  double **matrix_y;

  double *vector_continuous_stddev;
  double *vsfh, *vector_scale_factor, *vector_scale_factor_multistart;
  double *vector_scale_factor_startbest;

  double fret, fret_best, fret_start_best, fret_initial;
  double ftol, tol;
  double (* bwmfunc)(double *) = NULL;

  double small, lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init; 
  double lbd_init, hbd_init, d_init;
  int dfc_dir;
  
  int i,j;
  double fast_eval_total = 0.0;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, iImproved;
  int enforce_fixed_feasibility;
  int have_start_best, have_multistart_best;
  int itmax, iter;
  int int_use_starting_values;
  int scale_cat;
  const char *bw_error_msg = NULL;

  int * ipt = NULL;  // point permutation, see tree.c
  int old_bw;


  num_reg_unordered_extern = myopti[BW_NUNOI];
  num_reg_ordered_extern = myopti[BW_NORDI];
  num_reg_continuous_extern = myopti[BW_NCONI];
  np_reset_y_side_extern();
  bwm_clear_floor_context();

  vector_ckerlb_extern = ckerlb;
  vector_ckerub_extern = ckerub;
  int_cker_bound_extern = np_has_finite_cker_bounds(ckerlb, ckerub, num_reg_continuous_extern);

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[BW_NOBSI];
  iMultistart = myopti[BW_IMULTII];
  iNum_Multistart = myopti[BW_NMULTII];
  if (iNum_Multistart < 1)
    error("C_np_density_bw: nmulti must be a positive integer");

  KERNEL_den_extern = myopti[BW_CKRNEVI];
  KERNEL_den_unordered_extern = myopti[BW_UKRNEVI];
  KERNEL_den_ordered_extern = myopti[BW_OKRNEVI];

  int_use_starting_values= myopti[BW_USTARTI];
  int_LARGE_SF=myopti[BW_LSFI];
  BANDWIDTH_den_extern=myopti[BW_DENI];
  enforce_fixed_feasibility = (BANDWIDTH_den_extern == BW_FIXED);
  int_RESTART_FROM_MIN = myopti[BW_REMINI];
  int_MINIMIZE_IO = myopti[BW_MINIOI];

  itmax=myopti[BW_ITMAXI];
  old_bw=myopti[BW_OLDBW];
  int_TREE_X = myopti[BW_DOTREEI];
  scale_cat = myopti[BW_SCATI];
  bwm_use_transform = myopti[BW_TBNDI];
  if (BANDWIDTH_den_extern != BW_FIXED)
    bwm_use_transform = 0;
  if (bwm_use_transform) {
    int n = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
    bwm_reserve_transform_buf(n + 1);
  }

  ftol=myoptd[BW_FTOLD];
  tol=myoptd[BW_TOLD];
  small=myoptd[BW_SMALLD];

  dfc_dir = myopti[BW_DFC_DIRI];

  lbc_dir = myoptd[BW_LBC_DIRD];
  c_dir = myoptd[BW_C_DIRD];
  initc_dir = myoptd[BW_INITC_DIRD]; 
  lbd_dir = myoptd[BW_LBD_DIRD]; 
  hbd_dir = myoptd[BW_HBD_DIRD]; 
  d_dir = myoptd[BW_D_DIRD]; 
  initd_dir = myoptd[BW_INITD_DIRD]; 

  lbc_init = myoptd[BW_LBC_INITD]; 
  hbc_init = myoptd[BW_HBC_INITD]; 
  c_init = myoptd[BW_C_INITD]; 

  lbd_init = myoptd[BW_LBD_INITD]; 
  hbd_init = myoptd[BW_HBD_INITD]; 
  d_init = myoptd[BW_D_INITD]; 

  nconfac_extern = myoptd[BW_NCONFD];
  ncatfac_extern = myoptd[BW_NCATFD];
  bwm_set_scale_factor_lower_bound(myoptd[BW_SFLOORD]);

/* Allocate memory for objects */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);


  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_y = alloc_matd(num_var + 1, num_var +1);
  vector_scale_factor = alloc_vecd(num_var + 1);
  vector_scale_factor_startbest = alloc_vecd(num_var + 1);
  vsfh = alloc_vecd(num_var + 1);

  matrix_categorical_vals_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  
  if (int_use_starting_values)
    for( i=0;i<num_var; i++ )
      vector_scale_factor[i+1] = myans[i];

/* Parse data */

  for( j=0;j<num_reg_unordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_unordered_train_extern[j][i]=myuno[j*num_obs_train_extern+i];
    

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=myord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=mycon[j*num_obs_train_extern+i];

  ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  // attempt tree build, if enabled 
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                 4*num_reg_continuous_extern, ipt, &kdt_extern_X);
  

    //put training data into tree-order using the index array

    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_unordered_train_extern[j][i]=myuno[j*num_obs_train_extern+ipt[i]];
    
    
    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_ordered_train_extern[j][i]=myord[j*num_obs_train_extern+ipt[i]];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_continuous_train_extern[j][i]=mycon[j*num_obs_train_extern+ipt[i]];

  }


  determine_categorical_vals(
                             num_obs_train_extern,
                             0,
                             0,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             matrix_Y_unordered_train_extern,
                             matrix_Y_ordered_train_extern,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern,
                             matrix_categorical_vals_extern);

  vector_continuous_stddev = alloc_vecd(num_reg_continuous_extern);

  for (j = 0; j < num_reg_continuous_extern; j++)
    vector_continuous_stddev[j] = mysd[j];

  vector_continuous_stddev_extern = vector_continuous_stddev;

  /* Initialize scale factors and Hessian for NR modules */

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0, 
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0, 
                                    KERNEL_den_unordered_extern,                                    
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0, 
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0, 
                                    KERNEL_den_unordered_extern,                                    
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vsfh,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           0,
                           0,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED, 
                           lbc_dir, dfc_dir, c_dir,initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);

  /* When multistarting, set counter */

  imsnum = iMs_counter = 0;
  imstot = iNum_Multistart;
  
  /* Conduct direction set search */

  /* assign the function to be optimized */
  if(old_bw){
    switch(myopti[BW_MI]){
    case BWM_CVML : bwmfunc = cv_func_density_categorical_ml; break;
    case BWM_CVLS : bwmfunc = cv_func_density_categorical_ls; break;
      //case BWM_CVML_NP : bwmfunc = cv_func_np_density_categorical_ml; break;
    default : REprintf("np.c: invalid bandwidth selection method.");
      error("np.c: invalid bandwidth selection method."); break;
    }
  } else {
    switch(myopti[BW_MI]){
    case BWM_CVML : bwmfunc = np_cv_func_density_categorical_ml; break;
    case BWM_CVLS : bwmfunc = np_cv_func_density_categorical_ls; break;
    default : REprintf("np.c: invalid bandwidth selection method.");
      error("np.c: invalid bandwidth selection method."); break;
    }

  }

  if (bwm_use_transform)
    bwm_to_unconstrained(vector_scale_factor, num_var);

  spinner(0);

  np_bwm_clear_deferred_error();
  bwmfunc_raw = bwmfunc;
  bwm_num_reg_continuous = num_reg_continuous_extern;
  bwm_num_reg_unordered = num_reg_unordered_extern;
  bwm_num_reg_ordered = num_reg_ordered_extern;
  bwm_kernel_unordered = KERNEL_den_unordered_extern;
  bwm_kernel_unordered_vec = NULL;
  bwm_kernel_unordered_len = 0;
  bwm_num_categories = num_categories_extern;
  bwm_set_floor_context(
    enforce_fixed_feasibility,
    num_var,
    bwm_use_transform,
    KERNEL_den_extern,
    KERNEL_den_unordered_extern,
    BANDWIDTH_den_extern,
    BANDWIDTH_den_extern,
    0,
    num_obs_train_extern,
    0,
    0,
    0,
    num_reg_continuous_extern,
    num_reg_unordered_extern,
    num_reg_ordered_extern,
    num_categories_extern,
    bwm_scale_factor_lower_bound);
  bwm_reset_counters();
  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  if (penalty_mode[0] == 1) {
    double pmult = penalty_mult[0];
    double baseline;
    if (pmult < 1.0) pmult = 1.0;
    baseline = bwmfunc_raw_current_scale(vector_scale_factor, num_var);
    bwm_eval_count += 1.0;
    if (!R_FINITE(baseline) || baseline == DBL_MAX)
      bwm_invalid_count += 1.0;
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      double *tmp = alloc_vecd(num_var + 1);
      np_copy_scale_factor_for_raw(tmp, vector_scale_factor, num_var);
      for (i = 1; i <= num_reg_continuous_extern; i++)
        tmp[i] *= 2.0;
      for (i = 0; i < num_reg_unordered_extern; i++) {
        int idx = num_reg_continuous_extern + 1 + i;
        double maxbw = max_unordered_bw(num_categories_extern[i], KERNEL_den_unordered_extern);
        tmp[idx] = 0.5*maxbw;
      }
      for (i = 0; i < num_reg_ordered_extern; i++) {
        int idx = num_reg_continuous_extern + num_reg_unordered_extern + 1 + i;
        tmp[idx] = 0.5;
      }
      baseline = bwmfunc_raw(tmp);
      bwm_eval_count += 1.0;
      if (!R_FINITE(baseline) || baseline == DBL_MAX)
        bwm_invalid_count += 1.0;
      safe_free(tmp);
    }
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      bwm_penalty_value = pmult * 1.0e6;
    } else {
      bwm_penalty_value = baseline + (fabs(baseline) + 1.0) * pmult;
    }
    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }

  fret_initial = fret_best = bwmfunc_wrapper(vector_scale_factor);
  iImproved = 0;
  have_start_best = 0;
  have_multistart_best = 0;
  fret_start_best = DBL_MAX;
  if (enforce_fixed_feasibility &&
      np_bw_candidate_is_admissible(
        num_var,
        bwm_use_transform,
        KERNEL_den_extern,
        KERNEL_den_unordered_extern,
        BANDWIDTH_den_extern,
        BANDWIDTH_den_extern,
        0,
        num_obs_train_extern,
        0,
        0,
        0,
        num_reg_continuous_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_categories_extern,
        vector_scale_factor)) {
    have_start_best = 1;
    fret_start_best = fret_initial;
    np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor, num_var);
  }

  powell(0,
         0,
         vector_scale_factor,
         vector_scale_factor,
         matrix_y,
         num_var,
         ftol,
         tol,
         small,
         itmax,
         &iter,
         &fret,
         bwmfunc_wrapper);

  /* int_RESTART_FROM_MIN needs to be set */

  if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

    initialize_nr_directions(BANDWIDTH_den_extern,
                             num_obs_train_extern,
                             num_reg_continuous_extern,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             0,
                             0,
                             0,
                             vsfh,
                             num_categories_extern,
                             matrix_y,
                             0, int_RANDOM_SEED, 
                             lbc_dir, dfc_dir, c_dir, initc_dir,
                             lbd_dir, hbd_dir, d_dir, initd_dir,
                             matrix_X_continuous_train_extern,
                             matrix_Y_continuous_train_extern);



    powell(0,
           0,
           vector_scale_factor,
           vector_scale_factor,
           matrix_y,
           num_var,
           ftol,
           tol,
           small,
           itmax,
           &iter,
           &fret,
           bwmfunc_wrapper);
  }

  if (enforce_fixed_feasibility &&
      np_bw_candidate_is_admissible(
        num_var,
        bwm_use_transform,
        KERNEL_den_extern,
        KERNEL_den_unordered_extern,
        BANDWIDTH_den_extern,
        BANDWIDTH_den_extern,
        0,
        num_obs_train_extern,
        0,
        0,
        0,
        num_reg_continuous_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_categories_extern,
        vector_scale_factor) &&
      ((!have_start_best) || (fret < fret_start_best))) {
    have_start_best = 1;
    fret_start_best = fret;
    np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor, num_var);
  }

  if (enforce_fixed_feasibility) {
    if (have_start_best) {
      fret = fret_start_best;
      np_copy_scale_factor(vector_scale_factor, vector_scale_factor_startbest, num_var);
    } else {
      fret = DBL_MAX;
    }
  }

  iImproved = (enforce_fixed_feasibility && have_start_best) ? (fret_start_best < fret_initial) : (fret < fret_best);
  *timing = timing_extern;

  /* When multistarting save initial minimum of objective function and scale factors */
  objective_function_values[0]=-fret;
  objective_function_evals[0]=bwm_eval_count;
  objective_function_invalid[0]=bwm_invalid_count;
  bwm_snapshot_fast_counters();
  fast_eval_total += bwm_fast_eval_count;

  if(iMultistart == IMULTI_TRUE){
    if (enforce_fixed_feasibility) {
      if (have_start_best) {
        have_multistart_best = 1;
        fret_best = fret;
      } else {
        have_multistart_best = 0;
        fret_best = DBL_MAX;
      }
    } else {
      fret_best = fret;
    }
    vector_scale_factor_multistart = alloc_vecd(num_var + 1);
    for(i = 1; i <= num_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
    np_progress_bandwidth_multistart_step(1, iNum_Multistart);
    		
    /* Conduct search from new random values of the search parameters */
       	
    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and directions for NR modules */

      initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                        1,                /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        0, 
                                        0,
                                        0,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        0, 
                                        KERNEL_den_unordered_extern,                                    
                                        0,
                                        scale_cat,
                                        pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                        nconfac_extern, ncatfac_extern,
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor,
                                        lbc_init, hbc_init, c_init, 
                                        lbd_init, hbd_init, d_init,
                                        matrix_X_continuous_train_extern,
                                        matrix_Y_continuous_train_extern);


      initialize_nr_directions(BANDWIDTH_den_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               0,
                               0,
                               0,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               1, int_RANDOM_SEED, 
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);



      if (bwm_use_transform)
        bwm_to_unconstrained(vector_scale_factor, num_var);

      /* Conduct direction set search */

      bwm_reset_counters();

      powell(0,
             0,
             vector_scale_factor,
             vector_scale_factor,
             matrix_y,
             num_var,
             ftol,
             tol,
             small,
             itmax,
             &iter,
             &fret,
             bwmfunc_wrapper);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE){
        initialize_nr_directions(BANDWIDTH_den_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 num_reg_unordered_extern,
                                 num_reg_ordered_extern,
                                 0,
                                 0,
                                 0,
                                 vsfh,
                                 num_categories_extern,
                                 matrix_y,
                                 0, int_RANDOM_SEED, 
                                 lbc_dir, dfc_dir, c_dir, initc_dir,
                                 lbd_dir, hbd_dir, d_dir, initd_dir,
                                 matrix_X_continuous_train_extern,
                                 matrix_Y_continuous_train_extern);



        powell(0,
               0,
               vector_scale_factor,
               vector_scale_factor,
               matrix_y,
               num_var,
               ftol,
               tol,
               small,
               itmax,
               &iter,
               &fret,
               bwmfunc_wrapper);
      }
       			
      /* If this run resulted in an improved minimum save information */
       		
      if (enforce_fixed_feasibility) {
        if (np_bw_candidate_is_admissible(
              num_var,
              bwm_use_transform,
              KERNEL_den_extern,
              KERNEL_den_unordered_extern,
              BANDWIDTH_den_extern,
              BANDWIDTH_den_extern,
              0,
              num_obs_train_extern,
              0,
              0,
              0,
              num_reg_continuous_extern,
              num_reg_unordered_extern,
              num_reg_ordered_extern,
              num_categories_extern,
              vector_scale_factor) &&
            ((!have_multistart_best) || (fret < fret_best))) {
          fret_best = fret;
          have_multistart_best = 1;
          iImproved = iMs_counter+1;
          *timing = timing_extern;
          np_copy_scale_factor(vector_scale_factor_multistart, vector_scale_factor, num_var);
        }
      } else if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        *timing = timing_extern;
       
        for(i = 1; i <= num_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }
      objective_function_values[iMs_counter]=-fret;
      objective_function_evals[iMs_counter]=bwm_eval_count;
      objective_function_invalid[iMs_counter]=bwm_invalid_count;
      bwm_snapshot_fast_counters();
      fast_eval_total += bwm_fast_eval_count;
      np_progress_bandwidth_multistart_step(iMs_counter+1, iNum_Multistart);
    }

    /* Save best for estimation */

    if (enforce_fixed_feasibility) {
      if (have_multistart_best) {
        fret = fret_best;
        np_copy_scale_factor(vector_scale_factor, vector_scale_factor_multistart, num_var);
        have_start_best = 1;
        fret_start_best = fret_best;
        np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor_multistart, num_var);
      } else {
        have_start_best = 0;
        fret = DBL_MAX;
      }
    } else {
      fret = fret_best;

      for(i = 1; i <= num_var; i++)
        vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];
    }

    free(vector_scale_factor_multistart);

  }

  if (enforce_fixed_feasibility) {
    double final_raw;
    if (!have_start_best) {
      bw_error_msg = "C_np_density_bw: optimizer failed to produce a feasible fixed-bandwidth candidate";
      goto cleanup_np_density_bw;
    }
    if (!np_bw_candidate_is_admissible(
          num_var,
          bwm_use_transform,
          KERNEL_den_extern,
          KERNEL_den_unordered_extern,
          BANDWIDTH_den_extern,
          BANDWIDTH_den_extern,
          0,
          num_obs_train_extern,
          0,
          0,
          0,
          num_reg_continuous_extern,
          num_reg_unordered_extern,
          num_reg_ordered_extern,
          num_categories_extern,
          vector_scale_factor)) {
      bw_error_msg = "C_np_density_bw: optimizer returned an infeasible fixed-bandwidth candidate";
      goto cleanup_np_density_bw;
    }
    final_raw = bwmfunc_raw_current_scale(vector_scale_factor, num_var);
    if (!R_FINITE(final_raw) || final_raw == DBL_MAX) {
      bw_error_msg = "C_np_density_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective";
      goto cleanup_np_density_bw;
    }
    fret = final_raw;
  }

  if (bwm_use_transform)
    bwm_to_constrained(vector_scale_factor, num_var);

  /* return data to R */
  if (BANDWIDTH_den_extern == BW_GEN_NN || 
      BANDWIDTH_den_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern; i++ )
      vector_scale_factor[i+1]=np_fround(vector_scale_factor[i+1]);
  }

  for( i=0; i<num_var; i++ )
    myans[i]=vector_scale_factor[i+1];

  fval[0] = -fret;
  fval[1] = iImproved;
  objective_function_fast[0] = fast_eval_total;

  /* end return data */

cleanup_np_density_bw:
  /* Free data objects */
  bwm_clear_floor_context();

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  free_mat(matrix_y, num_var + 1);
  free(vector_scale_factor);
  free(vector_scale_factor_startbest);
  free(vsfh);
  free(num_categories_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  free(vector_continuous_stddev);

  free(ipt);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  int_cker_bound_extern = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;
  np_reset_y_side_extern();
  np_clear_estimator_extern_aliases();

  if (bw_error_msg != NULL)
    error("%s", bw_error_msg);

  return ;
  
}

 
// For distributions the bandwidth selection involves evaluating the CDF at a number of points.
// We allow one to specify those points, passed in by mye{uno,ord,con}.

void np_distribution_bw(double * myuno, double * myord, double * mycon, 
                        double * myeuno, double * myeord, double * myecon, double * mysd,
                        int * myopti, double * myoptd, double * myans, double * fval,
                        double * objective_function_values, double * objective_function_evals,
                        double * objective_function_invalid, double * timing,
                        double * objective_function_fast,
                        int * penalty_mode, double * penalty_mult,
                        double * ckerlb, double * ckerub){
  int_nn_k_min_extern = 1;
  /* Likelihood bandwidth selection for density estimation */

  double **matrix_y;

  double *vector_continuous_stddev;
  double *vsfh, *vector_scale_factor, *vector_scale_factor_multistart;
  double *vector_scale_factor_startbest;

  double fret, fret_best, fret_start_best, fret_initial;
  double fast_eval_total = 0.0;
  double ftol, tol;
  double (* bwmfunc)(double *) = NULL;

  double small, lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init; 
  double lbd_init, hbd_init, d_init;
  int dfc_dir;
  
  int i,j;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, iImproved;
  int enforce_fixed_feasibility;
  int have_start_best, have_multistart_best;
  int itmax, iter;
  int int_use_starting_values, cdfontrain;

  int scale_cat;
  const char *bw_error_msg = NULL;

  int * ipt = NULL, * ipe = NULL;

  cdfontrain_extern = cdfontrain =  myopti[DBW_CDFONTRAIN];

  num_reg_unordered_extern = myopti[DBW_NUNOI];
  num_reg_ordered_extern = myopti[DBW_NORDI];
  num_reg_continuous_extern = myopti[DBW_NCONI];
  np_reset_y_side_extern();
  bwm_clear_floor_context();

  vector_ckerlb_extern = ckerlb;
  vector_ckerub_extern = ckerub;
  int_cker_bound_extern = np_has_finite_cker_bounds(ckerlb, ckerub, num_reg_continuous_extern);

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[DBW_NOBSI];
  
  num_obs_eval_extern = cdfontrain ? num_obs_train_extern : myopti[DBW_NEVALI];

  iMultistart = myopti[DBW_IMULTII];
  iNum_Multistart = myopti[DBW_NMULTII];
  if (iNum_Multistart < 1)
    error("C_np_distribution_bw: nmulti must be a positive integer");

  KERNEL_den_extern = myopti[DBW_CKRNEVI];
  KERNEL_den_unordered_extern = myopti[DBW_UKRNEVI];
  KERNEL_den_ordered_extern = myopti[DBW_OKRNEVI];

  int_use_starting_values= myopti[DBW_USTARTI];
  int_LARGE_SF=myopti[DBW_LSFI];
  BANDWIDTH_den_extern=myopti[DBW_DENI];
  enforce_fixed_feasibility = (BANDWIDTH_den_extern == BW_FIXED);
  int_RESTART_FROM_MIN = myopti[DBW_REMINI];
  int_MINIMIZE_IO = myopti[DBW_MINIOI];

  itmax=myopti[DBW_ITMAXI];

  int_TREE_X = myopti[DBW_DOTREEI];
  scale_cat = myopti[DBW_SCATI];
  bwm_use_transform = myopti[DBW_TBNDI];
  if (BANDWIDTH_den_extern != BW_FIXED)
    bwm_use_transform = 0;
  if (bwm_use_transform) {
    int n = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
    bwm_reserve_transform_buf(n + 1);
  }

  ftol=myoptd[DBW_FTOLD];
  tol=myoptd[DBW_TOLD];
  small=myoptd[DBW_SMALLD];

  dfc_dir = myopti[DBW_DFC_DIRI];
  lbc_dir = myoptd[DBW_LBC_DIRD];
  c_dir = myoptd[DBW_C_DIRD];
  initc_dir = myoptd[DBW_INITC_DIRD]; 

  lbd_dir = myoptd[DBW_LBD_DIRD]; 
  hbd_dir = myoptd[DBW_HBD_DIRD]; 
  d_dir = myoptd[DBW_D_DIRD]; 
  initd_dir = myoptd[DBW_INITD_DIRD]; 

  lbc_init = myoptd[DBW_LBC_INITD]; 
  hbc_init = myoptd[DBW_HBC_INITD]; 
  c_init = myoptd[DBW_C_INITD]; 

  lbd_init = myoptd[DBW_LBD_INITD]; 
  hbd_init = myoptd[DBW_HBD_INITD]; 
  d_init = myoptd[DBW_D_INITD]; 

  nconfac_extern = myoptd[DBW_NCONFD];
  ncatfac_extern = myoptd[DBW_NCATFD];

  dbl_memfac_dls_extern = myoptd[DBW_MEMORYD];
  bwm_set_scale_factor_lower_bound(myoptd[DBW_SFLOORD]);

/* Allocate memory for objects */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  if(cdfontrain){
    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;
  } else {
    matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
    matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
    matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);
  }

  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_y = alloc_matd(num_var + 1, num_var +1);
  vector_scale_factor = alloc_vecd(num_var + 1);
  vector_scale_factor_startbest = alloc_vecd(num_var + 1);
  vsfh = alloc_vecd(num_var + 1);
  // nb check vals
  matrix_categorical_vals_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  if(num_reg_unordered_extern > 0)
    error("np.c: distribution bw selection only works on ordered and continuous data."); 

  if (int_use_starting_values)
    for( i=0;i<num_var; i++ )
      vector_scale_factor[i+1] = myans[i];

/* Parse data */

  for( j=0;j<num_reg_unordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_unordered_train_extern[j][i]=myuno[j*num_obs_train_extern+i];
    

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=myord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=mycon[j*num_obs_train_extern+i];


  // points for evaluating the cdf
  if(!cdfontrain){
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_unordered_eval_extern[j][i]=myeuno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=myeord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=myecon[j*num_obs_eval_extern+i];
  }

  ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  if(!cdfontrain) {
    ipe = (int *)malloc(num_obs_eval_extern*sizeof(int));
    if(!(ipe != NULL))
      error("!(ipe != NULL)");

    for(i = 0; i < num_obs_eval_extern; i++){
      ipe[i] = i;
    }
  } else {
    ipe = ipt;
  }

  // attempt tree build, if enabled 
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    if((BANDWIDTH_den_extern != BW_ADAP_NN) || ((BANDWIDTH_den_extern == BW_ADAP_NN) && cdfontrain)){
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipt, &kdt_extern_X);

      //put training data into tree-order using the index array

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_unordered_train_extern[j][i]=myuno[j*num_obs_train_extern+ipt[i]];  
    
      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_ordered_train_extern[j][i]=myord[j*num_obs_train_extern+ipt[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_continuous_train_extern[j][i]=mycon[j*num_obs_train_extern+ipt[i]];

    } else {
      build_kdtree(matrix_X_continuous_eval_extern, num_obs_eval_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipe, &kdt_extern_X);      

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_unordered_eval_extern[j][i]=myeuno[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_ordered_eval_extern[j][i]=myeord[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_continuous_eval_extern[j][i]=myecon[j*num_obs_eval_extern+ipe[i]];

    }

  }

  determine_categorical_vals(
                             num_obs_train_extern,
                             0,
                             0,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             matrix_Y_unordered_train_extern,
                             matrix_Y_ordered_train_extern,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern,
                             matrix_categorical_vals_extern);

  vector_continuous_stddev = alloc_vecd(num_reg_continuous_extern);

  for (j = 0; j < num_reg_continuous_extern; j++)
    vector_continuous_stddev[j] = mysd[j];

  vector_continuous_stddev_extern = vector_continuous_stddev;


  /* Initialize scale factors and Directions for NR modules */

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0,
                                    KERNEL_den_unordered_extern,
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_eval_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0,
                                    KERNEL_den_unordered_extern,
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vsfh,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_eval_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           0,
                           0,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED, 
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);


  /* When multistarting, set counter */

  imsnum = iMs_counter = 0;
  imstot = iNum_Multistart;
  
  /* Conduct direction set search */

  /* assign the function to be optimized */
  switch(myopti[DBW_MI]){
  case DBWM_CVLS : bwmfunc = cv_func_distribution_categorical_ls; break;
  default : REprintf("np.c: invalid bandwidth selection method.");
    error("np.c: invalid bandwidth selection method."); break;
  }

  if (bwm_use_transform)
    bwm_to_unconstrained(vector_scale_factor, num_var);

  spinner(0);

  np_bwm_clear_deferred_error();
  bwmfunc_raw = bwmfunc;
  bwm_num_reg_continuous = num_reg_continuous_extern;
  bwm_num_reg_unordered = num_reg_unordered_extern;
  bwm_num_reg_ordered = num_reg_ordered_extern;
  bwm_kernel_unordered = KERNEL_den_unordered_extern;
  bwm_kernel_unordered_vec = NULL;
  bwm_kernel_unordered_len = 0;
  bwm_num_categories = num_categories_extern;
  bwm_set_floor_context(
    enforce_fixed_feasibility,
    num_var,
    bwm_use_transform,
    KERNEL_den_extern,
    KERNEL_den_unordered_extern,
    BANDWIDTH_den_extern,
    BANDWIDTH_den_extern,
    0,
    num_obs_train_extern,
    0,
    0,
    0,
    num_reg_continuous_extern,
    num_reg_unordered_extern,
    num_reg_ordered_extern,
    num_categories_extern,
    bwm_scale_factor_lower_bound);
  bwm_reset_counters();
  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  if (penalty_mode[0] == 1) {
    double pmult = penalty_mult[0];
    double baseline;
    if (pmult < 1.0) pmult = 1.0;
    baseline = bwmfunc_raw_current_scale(vector_scale_factor, num_var);
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      double *tmp = alloc_vecd(num_var + 1);
      np_copy_scale_factor_for_raw(tmp, vector_scale_factor, num_var);
      for (i = 1; i <= num_reg_continuous_extern; i++)
        tmp[i] *= 2.0;
      for (i = 0; i < num_reg_unordered_extern; i++) {
        int idx = num_reg_continuous_extern + 1 + i;
        double maxbw = max_unordered_bw(num_categories_extern[i], KERNEL_den_unordered_extern);
        tmp[idx] = 0.5*maxbw;
      }
      for (i = 0; i < num_reg_ordered_extern; i++) {
        int idx = num_reg_continuous_extern + num_reg_unordered_extern + 1 + i;
        tmp[idx] = 0.5;
      }
      baseline = bwmfunc_raw(tmp);
      safe_free(tmp);
    }
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      bwm_penalty_value = pmult * 1.0e6;
    } else {
      bwm_penalty_value = baseline + (fabs(baseline) + 1.0) * pmult;
    }
    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }

  fret_initial = fret_best = bwmfunc_wrapper(vector_scale_factor);
  iImproved = 0;
  have_start_best = 0;
  have_multistart_best = 0;
  fret_start_best = DBL_MAX;
  if (enforce_fixed_feasibility &&
      np_bw_candidate_is_admissible(
        num_var,
        bwm_use_transform,
        KERNEL_den_extern,
        KERNEL_den_unordered_extern,
        BANDWIDTH_den_extern,
        BANDWIDTH_den_extern,
        0,
        num_obs_train_extern,
        0,
        0,
        0,
        num_reg_continuous_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_categories_extern,
        vector_scale_factor)) {
    have_start_best = 1;
    fret_start_best = fret_initial;
    np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor, num_var);
  }

  powell(0,
         0,
         vector_scale_factor,
         vector_scale_factor,
         matrix_y,
         num_var,
         ftol,
         tol,
         small,
         itmax,
         &iter,
         &fret,
         bwmfunc_wrapper);

  /* int_RESTART_FROM_MIN needs to be set */

  if(int_RESTART_FROM_MIN == RE_MIN_TRUE){
    initialize_nr_directions(BANDWIDTH_den_extern,
                             num_obs_train_extern,
                             num_reg_continuous_extern,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             0,
                             0,
                             0,
                             vsfh,
                             num_categories_extern,
                             matrix_y,
                             0, int_RANDOM_SEED,  
                             lbc_dir, dfc_dir, c_dir, initc_dir,
                             lbd_dir, hbd_dir, d_dir, initd_dir,
                             matrix_X_continuous_train_extern,
                             matrix_Y_continuous_train_extern);



    powell(0,
           0,
           vector_scale_factor,
           vector_scale_factor,
           matrix_y,
           num_var,
           ftol,
           tol,
           small,
           itmax,
           &iter,
           &fret,
           bwmfunc_wrapper);
  }

  if (enforce_fixed_feasibility &&
      np_bw_candidate_is_admissible(
        num_var,
        bwm_use_transform,
        KERNEL_den_extern,
        KERNEL_den_unordered_extern,
        BANDWIDTH_den_extern,
        BANDWIDTH_den_extern,
        0,
        num_obs_train_extern,
        0,
        0,
        0,
        num_reg_continuous_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_categories_extern,
        vector_scale_factor) &&
      ((!have_start_best) || (fret < fret_start_best))) {
    have_start_best = 1;
    fret_start_best = fret;
    np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor, num_var);
  }

  if (enforce_fixed_feasibility) {
    if (have_start_best) {
      fret = fret_start_best;
      np_copy_scale_factor(vector_scale_factor, vector_scale_factor_startbest, num_var);
    } else {
      fret = DBL_MAX;
    }
  }

  iImproved = (enforce_fixed_feasibility && have_start_best) ? (fret_start_best < fret_initial) : (fret < fret_best);
  *timing = timing_extern;

  objective_function_values[0]=fret;
  objective_function_evals[0]=bwm_eval_count;
  objective_function_invalid[0]=bwm_invalid_count;
  bwm_snapshot_fast_counters();
  fast_eval_total += bwm_fast_eval_count;
  /* When multistarting save initial minimum of objective function and scale factors */

  if(iMultistart == IMULTI_TRUE){
    if (enforce_fixed_feasibility) {
      if (have_start_best) {
        have_multistart_best = 1;
        fret_best = fret;
      } else {
        have_multistart_best = 0;
        fret_best = DBL_MAX;
      }
    } else {
      fret_best = fret;
    }
    vector_scale_factor_multistart = alloc_vecd(num_var + 1);
    for(i = 1; i <= num_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
    np_progress_bandwidth_multistart_step(1, iNum_Multistart);
    		
    /* Conduct search from new random values of the search parameters */
       	
    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and directions for NR modules */

      initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                        1,        /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        0,
                                        0,
                                        0,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        0,
                                        KERNEL_den_unordered_extern,
                                        0,
                                        scale_cat,
                                        pow((double)4.0/(double)3.0,0.2),     /* Init for continuous vars */
                                        nconfac_extern, ncatfac_extern,
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor,
                                        lbc_init, hbc_init, c_init, 
                                        lbd_init, hbd_init, d_init,
                                        matrix_X_continuous_eval_extern,
                                        matrix_Y_continuous_train_extern);

      initialize_nr_directions(BANDWIDTH_den_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               0,
                               0,
                               0,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               1, int_RANDOM_SEED,  
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);

      /* Conduct direction set search */

      if (bwm_use_transform)
        bwm_to_unconstrained(vector_scale_factor, num_var);

      bwm_reset_counters();

      powell(0,
             0,
             vector_scale_factor,
             vector_scale_factor,
             matrix_y,
             num_var,
             ftol,
             tol,
             small,
             itmax,
             &iter,
             &fret,
             bwmfunc_wrapper);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

        initialize_nr_directions(BANDWIDTH_den_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 num_reg_unordered_extern,
                                 num_reg_ordered_extern,
                                 0,
                                 0,
                                 0,
                                 vsfh,
                                 num_categories_extern,
                                 matrix_y,
                                 0, int_RANDOM_SEED,  
                                 lbc_dir, dfc_dir, c_dir, initc_dir,
                                 lbd_dir, hbd_dir, d_dir, initd_dir,
                                 matrix_X_continuous_train_extern,
                                 matrix_Y_continuous_train_extern);

        powell(0,
               0,
               vector_scale_factor,
               vector_scale_factor,
               matrix_y,
               num_var,
               ftol,
               tol,
               small,
               itmax,
               &iter,
               &fret,
               bwmfunc_wrapper);
      }
       			
      /* If this run resulted in an improved minimum save information */
       		
      if (enforce_fixed_feasibility) {
        if (np_bw_candidate_is_admissible(
              num_var,
              bwm_use_transform,
              KERNEL_den_extern,
              KERNEL_den_unordered_extern,
              BANDWIDTH_den_extern,
              BANDWIDTH_den_extern,
              0,
              num_obs_train_extern,
              0,
              0,
              0,
              num_reg_continuous_extern,
              num_reg_unordered_extern,
              num_reg_ordered_extern,
              num_categories_extern,
              vector_scale_factor) &&
            ((!have_multistart_best) || (fret < fret_best))) {
          fret_best = fret;
          have_multistart_best = 1;
          iImproved = iMs_counter+1;
          *timing = timing_extern;
          np_copy_scale_factor(vector_scale_factor_multistart, vector_scale_factor, num_var);
        }
      } else if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        *timing = timing_extern;
       
        for(i = 1; i <= num_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }
      objective_function_values[iMs_counter]=fret;
      objective_function_evals[iMs_counter]=bwm_eval_count;
      objective_function_invalid[iMs_counter]=bwm_invalid_count;
      bwm_snapshot_fast_counters();
      fast_eval_total += bwm_fast_eval_count;
      np_progress_bandwidth_multistart_step(iMs_counter+1, iNum_Multistart);
    }

    /* Save best for estimation */

    if (enforce_fixed_feasibility) {
      if (have_multistart_best) {
        fret = fret_best;
        np_copy_scale_factor(vector_scale_factor, vector_scale_factor_multistart, num_var);
        have_start_best = 1;
        fret_start_best = fret_best;
        np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor_multistart, num_var);
      } else {
        have_start_best = 0;
        fret = DBL_MAX;
      }
    } else {
      fret = fret_best;

      for(i = 1; i <= num_var; i++)
        vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];
    }

    free(vector_scale_factor_multistart);

  }

  if (enforce_fixed_feasibility) {
    double final_raw;
    if (!have_start_best) {
      bw_error_msg = "C_np_distribution_bw: optimizer failed to produce a feasible fixed-bandwidth candidate";
      goto cleanup_np_distribution_bw;
    }
    if (!np_bw_candidate_is_admissible(
          num_var,
          bwm_use_transform,
          KERNEL_den_extern,
          KERNEL_den_unordered_extern,
          BANDWIDTH_den_extern,
          BANDWIDTH_den_extern,
          0,
          num_obs_train_extern,
          0,
          0,
          0,
          num_reg_continuous_extern,
          num_reg_unordered_extern,
          num_reg_ordered_extern,
          num_categories_extern,
          vector_scale_factor)) {
      bw_error_msg = "C_np_distribution_bw: optimizer returned an infeasible fixed-bandwidth candidate";
      goto cleanup_np_distribution_bw;
    }
    final_raw = bwmfunc_raw_current_scale(vector_scale_factor, num_var);
    if (!R_FINITE(final_raw) || final_raw == DBL_MAX) {
      bw_error_msg = "C_np_distribution_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective";
      goto cleanup_np_distribution_bw;
    }
    fret = final_raw;
  }

  if (bwm_use_transform)
    bwm_to_constrained(vector_scale_factor, num_var);

  /* return data to R */
  if (BANDWIDTH_den_extern == BW_GEN_NN || 
      BANDWIDTH_den_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern; i++ )
      vector_scale_factor[i+1]=np_fround(vector_scale_factor[i+1]);
  }

  for( i=0; i<num_var; i++ )
    myans[i]=vector_scale_factor[i+1];

  fval[0] = fret;
  fval[1] = iImproved;
  fval[2] = fret_initial;
  objective_function_fast[0] = fast_eval_total;

  /* end return data */

cleanup_np_distribution_bw:
  /* Free data objects */
  bwm_clear_floor_context();

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  if(!cdfontrain){
    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  free_mat(matrix_y, num_var + 1);
  free(vector_scale_factor);
  free(vector_scale_factor_startbest);
  free(vsfh);
  free(num_categories_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  free(vector_continuous_stddev);

  free(ipt);
  if(!cdfontrain)
    free(ipe);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  int_cker_bound_extern = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;
  np_reset_y_side_extern();
  np_clear_estimator_extern_aliases();

  if (bw_error_msg != NULL)
    error("%s", bw_error_msg);

  return ;
  
}


void np_density_conditional_bw(double * c_uno, double * c_ord, double * c_con, 
                               double * u_uno, double * u_ord, double * u_con,
                               double * mysd,
                               int * myopti, double * myoptd, double * myans, double * fval,
                               double * objective_function_values, double * objective_function_evals,
                               double * objective_function_invalid, double * timing,
                               double * objective_function_fast,
                               int * penalty_mode, double * penalty_mult,
                               int * glp_degree,
                               int * glp_bernstein,
                               int * glp_basis,
                               int * regtype,
                               double * cxkerlb, double * cxkerub,
                               double * cykerlb, double * cykerub,
                               const int eval_only){
  int_nn_k_min_extern = 1;
/* Likelihood bandwidth selection for density estimation */

  double **matrix_y = NULL;

  double *vector_continuous_stddev;
  double *vsfh, *vector_scale_factor, *vector_scale_factor_multistart;
  double *vector_scale_factor_startbest;
  double *cxylb = NULL, *cxyub = NULL;

  double fret, fret_best, fret_start_best, fret_initial;
  double ftol, tol;
  double (* bwmfunc)(double *) = NULL;

  double small, lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init; 
  double lbd_init, hbd_init, d_init;
  double scale_factor_lower_bound;
  int dfc_dir;
  
  int i,j;
  int need_y_side = 0;
  double fast_eval_total = 0.0;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, num_all_var, num_var_var, iImproved;
  int enforce_fixed_feasibility;
  int have_start_best, have_multistart_best;
  int itmax, iter;
  int int_use_starting_values, ibwmfunc, old_cdens, scale_cat;
  const char *bw_error_msg = NULL;

  int num_all_cvar, num_all_uvar, num_all_ovar;

  int * ipt_X = NULL, * ipt_XY = NULL, * ipt_Y = NULL; 
  int * ipt_lookup_XY = NULL, * ipt_lookup_Y = NULL, * ipt_lookup_X = NULL;

  /* Ensure optional Y-only categorical arrays are reset each call */
  num_categories_extern_Y = NULL;
  matrix_categorical_vals_extern_Y = NULL;
  np_bounded_cvls_conditional_quad_context_clear_extern();
  bwm_clear_floor_context();

  num_var_unordered_extern = myopti[CBW_CNUNOI];
  num_var_ordered_extern = myopti[CBW_CNORDI];
  num_var_continuous_extern = myopti[CBW_CNCONI];

  num_reg_unordered_extern = myopti[CBW_UNUNOI];
  num_reg_ordered_extern = myopti[CBW_UNORDI];
  num_reg_continuous_extern = myopti[CBW_UNCONI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;
  num_var_var = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  num_all_var = num_var+num_var_var;

  num_obs_train_extern = myopti[CBW_NOBSI];
  iMultistart = myopti[CBW_IMULTII];
  iNum_Multistart = myopti[CBW_NMULTII];
  if (!eval_only && iNum_Multistart < 1)
    error("C_np_density_conditional_bw: nmulti must be a positive integer");

  KERNEL_reg_extern = myopti[CBW_CXKRNEVI];
  KERNEL_den_extern = myopti[CBW_CYKRNEVI];

  KERNEL_reg_unordered_extern = myopti[CBW_UXKRNEVI];
  KERNEL_den_unordered_extern = myopti[CBW_UYKRNEVI];

  KERNEL_reg_ordered_extern = myopti[CBW_OXKRNEVI];
  KERNEL_den_ordered_extern = myopti[CBW_OYKRNEVI];

  vector_cxkerlb_extern = cxkerlb;
  vector_cxkerub_extern = cxkerub;
  int_cxker_bound_extern = np_has_finite_cker_bounds(cxkerlb, cxkerub, num_reg_continuous_extern);

  vector_cykerlb_extern = cykerlb;
  vector_cykerub_extern = cykerub;
  int_cyker_bound_extern = np_has_finite_cker_bounds(cykerlb, cykerub, num_var_continuous_extern);

  if((num_reg_continuous_extern + num_var_continuous_extern) > 0){
    cxylb = alloc_vecd(num_reg_continuous_extern + num_var_continuous_extern);
    cxyub = alloc_vecd(num_reg_continuous_extern + num_var_continuous_extern);
    for(i = 0; i < num_reg_continuous_extern; i++){
      cxylb[i] = (cxkerlb != NULL) ? cxkerlb[i] : DBL_MAX;
      cxyub[i] = (cxkerub != NULL) ? cxkerub[i] : DBL_MAX;
    }
    for(i = 0; i < num_var_continuous_extern; i++){
      cxylb[num_reg_continuous_extern + i] = (cykerlb != NULL) ? cykerlb[i] : DBL_MAX;
      cxyub[num_reg_continuous_extern + i] = (cykerub != NULL) ? cykerub[i] : DBL_MAX;
    }
    vector_cxykerlb_extern = cxylb;
    vector_cxykerub_extern = cxyub;
    int_cxyker_bound_extern = np_has_finite_cker_bounds(cxylb, cxyub, num_reg_continuous_extern + num_var_continuous_extern);
  } else {
    vector_cxykerlb_extern = NULL;
    vector_cxykerub_extern = NULL;
    int_cxyker_bound_extern = 0;
  }

  int_use_starting_values= myopti[CBW_USTARTI];
  int_LARGE_SF=myopti[CBW_LSFI];
  BANDWIDTH_den_extern=myopti[CBW_DENI];
  enforce_fixed_feasibility = ((BANDWIDTH_den_extern == BW_FIXED) && (!eval_only));
  int_RESTART_FROM_MIN = myopti[CBW_REMINI];
  int_MINIMIZE_IO = myopti[CBW_MINIOI];

  itmax=myopti[CBW_ITMAXI];
  int_WEIGHTS = 0;
  old_cdens = myopti[CBW_OLDI];
  int_TREE_XY = int_TREE_Y = int_TREE_X = myopti[CBW_TREEI];
  scale_cat = myopti[CBW_SCATI];
  bwm_use_transform = 0;
  
  ibwmfunc = myopti[CBW_MI];
  int_ll_extern = ((ibwmfunc == CBWM_CVML) || (ibwmfunc == CBWM_CVLS)) ? *regtype : LL_LC;
  vector_glp_degree_extern = (((ibwmfunc == CBWM_CVML) || (ibwmfunc == CBWM_CVLS)) && (int_ll_extern == LL_LP)) ? glp_degree : NULL;
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = (((ibwmfunc == CBWM_CVML) || (ibwmfunc == CBWM_CVLS)) && (int_ll_extern == LL_LP)) ? *glp_bernstein : 0;
  int_glp_basis_extern = (((ibwmfunc == CBWM_CVML) || (ibwmfunc == CBWM_CVLS)) && (int_ll_extern == LL_LP)) ? *glp_basis : 1;
  int_bounded_cvls_quadrature_grid_extern = myopti[CBW_CVLS_QUAD_GRIDI];
  if ((int_bounded_cvls_quadrature_grid_extern < 0) ||
      (int_bounded_cvls_quadrature_grid_extern > 2))
    error("C_np_density_conditional_bw: cvls.quadrature.grid is invalid");
  int_bounded_cvls_quadrature_points_extern = myopti[CBW_CVLS_QUAD_POINTSI];
  if (int_bounded_cvls_quadrature_points_extern < 2)
    int_bounded_cvls_quadrature_points_extern = 0;
  need_y_side = (ibwmfunc == CBWM_CVLS) || ((ibwmfunc == CBWM_CVML) && (int_ll_extern == LL_LP));
  bwm_use_transform = myopti[CBW_TBNDI];
  if (BANDWIDTH_den_extern != BW_FIXED)
    bwm_use_transform = 0;
  if (bwm_use_transform) {
    int n = num_var_continuous_extern + num_reg_continuous_extern +
      num_var_unordered_extern + num_reg_unordered_extern +
      num_var_ordered_extern + num_reg_ordered_extern;
    bwm_reserve_transform_buf(n + 1);
  }

  ftol=myoptd[CBW_FTOLD];
  tol=myoptd[CBW_TOLD];
  small=myoptd[CBW_SMALLD];
  dbl_memfac_ccdf_extern = myoptd[CBW_MEMFACD];
  scale_factor_lower_bound = myoptd[CBW_SFLOORD];
  if (!R_FINITE(scale_factor_lower_bound) || scale_factor_lower_bound < 0.0)
    scale_factor_lower_bound = 0.1;
  double_bounded_cvls_quadrature_extend_factor_extern = myoptd[CBW_QUAD_EXTD];
  if (!R_FINITE(double_bounded_cvls_quadrature_extend_factor_extern) ||
      double_bounded_cvls_quadrature_extend_factor_extern <= 0.0)
    error("C_np_density_conditional_bw: cvls.quadrature.extend.factor must be positive and finite");
  double_bounded_cvls_quadrature_ratios_extern[0] = myoptd[CBW_QUAD_RATIO_UNIFORMD];
  double_bounded_cvls_quadrature_ratios_extern[1] = myoptd[CBW_QUAD_RATIO_SAMPLED];
  double_bounded_cvls_quadrature_ratios_extern[2] = myoptd[CBW_QUAD_RATIO_GLD];
  if (!R_FINITE(double_bounded_cvls_quadrature_ratios_extern[0]) ||
      !R_FINITE(double_bounded_cvls_quadrature_ratios_extern[1]) ||
      !R_FINITE(double_bounded_cvls_quadrature_ratios_extern[2]) ||
      (double_bounded_cvls_quadrature_ratios_extern[0] < 0.0) ||
      (double_bounded_cvls_quadrature_ratios_extern[1] < 0.0) ||
      (double_bounded_cvls_quadrature_ratios_extern[2] < 0.0) ||
      (fabs(double_bounded_cvls_quadrature_ratios_extern[0] +
            double_bounded_cvls_quadrature_ratios_extern[1] +
            double_bounded_cvls_quadrature_ratios_extern[2] - 1.0) > 1.0e-8))
    error("C_np_density_conditional_bw: cvls.quadrature.ratios must be non-negative and sum to one");

  dfc_dir = myopti[CBW_DFC_DIRI];
  lbc_dir = myoptd[CBW_LBC_DIRD];
  c_dir = myoptd[CBW_C_DIRD];
  initc_dir = myoptd[CBW_INITC_DIRD]; 

  lbd_dir = myoptd[CBW_LBD_DIRD]; 
  hbd_dir = myoptd[CBW_HBD_DIRD]; 
  d_dir = myoptd[CBW_D_DIRD]; 
  initd_dir = myoptd[CBW_INITD_DIRD]; 

  lbc_init = myoptd[CBW_LBC_INITD]; 
  hbc_init = myoptd[CBW_HBC_INITD]; 
  c_init = myoptd[CBW_C_INITD]; 

  lbd_init = myoptd[CBW_LBD_INITD]; 
  hbd_init = myoptd[CBW_HBD_INITD]; 
  d_init = myoptd[CBW_D_INITD]; 

  nconfac_extern = myoptd[CBW_NCONFD];
  ncatfac_extern = myoptd[CBW_NCATFD];

/* Allocate memory for objects */

  //if((BANDWIDTH_den_extern != BW_FIXED) && (ibwmfunc == CBWM_CVLS))
  //old_cdens = 1;

  matrix_Y_unordered_train_extern = alloc_matd(num_obs_train_extern, num_var_unordered_extern);
  matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);
  matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);
	
  num_categories_extern = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                     num_reg_unordered_extern + num_reg_ordered_extern);

  num_categories_extern_XY = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                     num_reg_unordered_extern + num_reg_ordered_extern);

  num_categories_extern_X = alloc_vecu(num_reg_unordered_extern + num_reg_ordered_extern);
  
  if(need_y_side){
    num_categories_extern_Y = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern);
  }

  matrix_y = alloc_matd(num_all_var + 1, num_all_var + 1);
  vector_scale_factor = alloc_vecd(num_all_var + 1);
  vector_scale_factor_startbest = alloc_vecd(num_all_var + 1);
  vsfh = alloc_vecd(num_all_var + 1);
  
  matrix_categorical_vals_extern = 
    alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern + 
               num_reg_unordered_extern + num_reg_ordered_extern);

  matrix_categorical_vals_extern_X = 
    alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  if(need_y_side){
    matrix_categorical_vals_extern_Y = 
      alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern);
  }

  matrix_categorical_vals_extern_XY = 
    alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern + 
               num_reg_unordered_extern + num_reg_ordered_extern);


  /* in v_s_f order is creg, cvar, uvar, ovar, ureg, oreg  */

  if (int_use_starting_values)
    for( i=0;i<num_all_var; i++ )
      vector_scale_factor[i+1] = myans[i];

/* Parse data */

  for(j=0;j<num_var_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_unordered_train_extern[j][i]=c_uno[j*num_obs_train_extern+i];

  for(j=0;j<num_var_ordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_ordered_train_extern[j][i]=c_ord[j*num_obs_train_extern+i];

  for(j=0;j<num_var_continuous_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_continuous_train_extern[j][i]=c_con[j*num_obs_train_extern+i];


  for(j=0;j<num_reg_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_X_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+i];

  /* data has been parsed, make the joint xy dataset */

  // the main tree is x, the alt tree is xy

  num_all_cvar = num_reg_continuous_extern + num_var_continuous_extern;
  num_all_uvar = num_reg_unordered_extern + num_var_unordered_extern;
  num_all_ovar = num_reg_ordered_extern + num_var_ordered_extern;

  // we use 3 trees for cpdf ls, and 2 for cpdf ml

  ipt_X = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_X != NULL))
    error("!(ipt_X != NULL)");

  ipt_lookup_X = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_lookup_X != NULL)){
    safe_free(ipt_X);
    error("!(ipt_lookup_X != NULL)");
  }

  for(i = 0; i < num_obs_train_extern; i++){
    ipt_lookup_X[i] = ipt_X[i] = i;
  }

  ipt_extern_X = ipt_X;
  ipt_lookup_extern_X = ipt_lookup_X;

  if(need_y_side){
    ipt_Y = (int *)malloc(num_obs_train_extern*sizeof(int));
    if(!(ipt_Y != NULL)){
      safe_free(ipt_X);
      safe_free(ipt_lookup_X);
      error("!(ipt_Y != NULL)");
    }

    ipt_lookup_Y = (int *)malloc(num_obs_train_extern*sizeof(int));
    if(!(ipt_lookup_Y != NULL)){
      safe_free(ipt_X);
      safe_free(ipt_lookup_X);
      safe_free(ipt_Y);
      error("!(ipt_lookup_Y != NULL)");
    }

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_Y[i] = ipt_Y[i] = i;
    }

  } else {
    ipt_Y = ipt_X;
    ipt_lookup_Y = ipt_lookup_X;
  }

  ipt_extern_Y = ipt_Y;
  ipt_lookup_extern_Y = ipt_lookup_Y;

  ipt_XY = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_XY != NULL)){
    safe_free(ipt_X);
    safe_free(ipt_lookup_X);
    if(need_y_side){
      safe_free(ipt_Y);
      safe_free(ipt_lookup_Y);
    }
    error("!(ipt_XY != NULL)");
  }

  ipt_lookup_XY = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_lookup_XY != NULL)){
    safe_free(ipt_X);
    safe_free(ipt_lookup_X);
    if(need_y_side){
      safe_free(ipt_Y);
      safe_free(ipt_lookup_Y);
    }
    safe_free(ipt_XY);
    error("!(ipt_lookup_XY != NULL)");
  }

  for(i = 0; i < num_obs_train_extern; i++){
    ipt_lookup_XY[i] = ipt_XY[i] = i;
  }

  ipt_extern_XY = ipt_XY;
  ipt_lookup_extern_XY = ipt_lookup_XY;

  int_TREE_XY = int_TREE_XY && (((num_all_cvar) != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);

  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);

  int_TREE_Y = int_TREE_Y && need_y_side && ((num_var_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);


  if(int_TREE_X == NP_TREE_TRUE){
    build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                 4*num_reg_continuous_extern, ipt_X, &kdt_extern_X);
  

    //put training data into tree-order using the index array

    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+ipt_X[i]];
    
    
    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+ipt_X[i]];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+ipt_X[i]];

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_X[ipt_X[i]] = i;
    }
  }

  if(int_TREE_Y == NP_TREE_TRUE){
    build_kdtree(matrix_Y_continuous_train_extern, num_obs_train_extern, num_var_continuous_extern, 
                 4*num_var_continuous_extern, ipt_Y, &kdt_extern_Y);

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_Y[ipt_Y[i]] = i;
    }

  }


  if((int_TREE_X == NP_TREE_TRUE) || (int_TREE_Y == NP_TREE_TRUE)){
    for(j=0;j<num_var_unordered_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_unordered_train_extern[j][i]=c_uno[j*num_obs_train_extern+ipt_Y[i]];

    for(j=0;j<num_var_ordered_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_ordered_train_extern[j][i]=c_ord[j*num_obs_train_extern+ipt_Y[i]];

    for(j=0;j<num_var_continuous_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_continuous_train_extern[j][i]=c_con[j*num_obs_train_extern+ipt_Y[i]];

  }


  matrix_XY_continuous_train_extern = alloc_matd(num_obs_train_extern, num_all_cvar);
  matrix_XY_unordered_train_extern = alloc_matd(num_obs_train_extern, num_all_uvar);
  matrix_XY_ordered_train_extern = alloc_matd(num_obs_train_extern, num_all_ovar);


  for(j = 0; j < num_reg_unordered_extern; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+i];

  for(j = num_reg_unordered_extern; j < num_all_uvar; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_unordered_train_extern[j][i]=c_uno[(j-num_reg_unordered_extern)*num_obs_train_extern+i];


  for(j = 0; j < num_reg_ordered_extern; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+i];

  for(j = num_reg_ordered_extern; j < num_all_ovar; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_ordered_train_extern[j][i]=c_ord[(j-num_reg_ordered_extern)*num_obs_train_extern+i];


  for(j = 0; j < num_reg_continuous_extern; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+i];

  for(j = num_reg_continuous_extern; j < num_all_cvar; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_continuous_train_extern[j][i]=c_con[(j-num_reg_continuous_extern)*num_obs_train_extern+i];

  if(int_TREE_XY == NP_TREE_TRUE){

    build_kdtree(matrix_XY_continuous_train_extern, num_obs_train_extern, num_all_cvar, 
                 4*num_all_cvar, ipt_XY, &kdt_extern_XY);

    // put data into tree-order
    for(j = 0; j < num_reg_unordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_unordered_extern; j < num_all_uvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=c_uno[(j-num_reg_unordered_extern)*num_obs_train_extern+ipt_XY[i]];


    for(j = 0; j < num_reg_ordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_ordered_extern; j < num_all_ovar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=c_ord[(j-num_reg_ordered_extern)*num_obs_train_extern+ipt_XY[i]];


    for(j = 0; j < num_reg_continuous_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_continuous_extern; j < num_all_cvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=c_con[(j-num_reg_continuous_extern)*num_obs_train_extern+ipt_XY[i]];

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_XY[ipt_XY[i]] = i;
    }

  }

  determine_categorical_vals(num_obs_train_extern,
                             num_var_unordered_extern,
                             num_var_ordered_extern,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             matrix_Y_unordered_train_extern,
                             matrix_Y_ordered_train_extern,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern,
                             matrix_categorical_vals_extern);

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern, num_var_ordered_extern, num_var_continuous_extern,
                        num_reg_unordered_extern, num_reg_ordered_extern, num_reg_continuous_extern,
                        vector_scale_factor+1,
                        num_categories_extern,
                        matrix_categorical_vals_extern,
                        NULL, NULL, NULL,
                        num_categories_extern_X, need_y_side ? num_categories_extern_Y : NULL, num_categories_extern_XY,
                        matrix_categorical_vals_extern_X, need_y_side ? matrix_categorical_vals_extern_Y : NULL, matrix_categorical_vals_extern_XY);
  

  vector_continuous_stddev = alloc_vecd(num_var_continuous_extern + num_reg_continuous_extern);

  for(j = 0; j < (num_var_continuous_extern + num_reg_continuous_extern); j++)
    vector_continuous_stddev[j] = mysd[j];

  vector_continuous_stddev_extern = vector_continuous_stddev;

  if((ibwmfunc == CBWM_CVLS) && (int_ll_extern == LL_LP)){
    if(np_bounded_cvls_conditional_quad_context_prepare_extern() != 0){
      bw_error_msg = "C_np_density_conditional_bw: failed to prepare bounded cv.ls quadrature context";
      goto cleanup_np_density_conditional_bw;
    }
  }

  /* Initialize scale factors and Directions for NR modules */

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    num_var_continuous_extern,
                                    num_var_unordered_extern,
                                    num_var_ordered_extern,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    KERNEL_den_unordered_extern,
                                    KERNEL_reg_unordered_extern,
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    num_var_continuous_extern,
                                    num_var_unordered_extern,
                                    num_var_ordered_extern,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    KERNEL_den_unordered_extern,
                                    KERNEL_reg_unordered_extern,
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vsfh,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           num_var_continuous_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED,  
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);

  /* When multistarting, set counter */

  imsnum = iMs_counter = 0;
  imstot = iNum_Multistart;

  /* Conduct direction set search */

  /* assign the function to be optimized */

  /* 7/2/2010 */  
  if(old_cdens){
    switch(ibwmfunc){
    case CBWM_CVML : bwmfunc = cv_func_con_density_categorical_ml; break;
    case CBWM_CVLS : bwmfunc = cv_func_con_density_categorical_ls; break;
    case CBWM_NPLS : bwmfunc = np_cv_func_con_density_categorical_ls;break;
    case CBWM_CCDF : bwmfunc = cv_func_con_distribution_categorical_ccdf; break;
    default : REprintf("np.c: invalid bandwidth selection method.");
      error("np.c: invalid bandwidth selection method."); break;
    }
  } else {
    switch(ibwmfunc){
    case CBWM_CVML : bwmfunc = np_cv_func_con_density_categorical_ml; break;
    case CBWM_CVLS : bwmfunc = np_cv_func_con_density_categorical_ls_npksum; break;
    default : REprintf("np.c: invalid bandwidth selection method.");
      error("np.c: invalid bandwidth selection method."); break;
    }
  }

  if (bwm_use_transform)
    bwm_to_unconstrained(vector_scale_factor, num_all_var);

  spinner(0);

  np_bwm_clear_deferred_error();
  bwmfunc_raw = bwmfunc;
  bwm_num_reg_continuous = num_var_continuous_extern + num_reg_continuous_extern;
  bwm_num_reg_unordered = num_var_unordered_extern + num_reg_unordered_extern;
  bwm_num_reg_ordered = num_var_ordered_extern + num_reg_ordered_extern;
  bwm_kernel_unordered = KERNEL_den_unordered_extern;
  bwm_kernel_unordered_len = bwm_num_reg_unordered;
  if (bwm_kernel_unordered_len > 0) {
    bwm_kernel_unordered_vec = alloc_vecu(bwm_kernel_unordered_len);
    for (i = 0; i < num_var_unordered_extern; i++)
      bwm_kernel_unordered_vec[i] = KERNEL_den_unordered_extern;
    for (i = 0; i < num_reg_unordered_extern; i++)
      bwm_kernel_unordered_vec[num_var_unordered_extern + i] = KERNEL_reg_unordered_extern;
  } else {
    bwm_kernel_unordered_vec = NULL;
  }
  bwm_num_categories = num_categories_extern;
  bwm_set_floor_context(
    enforce_fixed_feasibility,
    num_all_var,
    bwm_use_transform,
    KERNEL_den_extern,
    KERNEL_reg_unordered_extern,
    BANDWIDTH_den_extern,
    BANDWIDTH_den_extern,
    0,
    num_obs_train_extern,
    num_var_continuous_extern,
    num_var_unordered_extern,
    num_var_ordered_extern,
    num_reg_continuous_extern,
    num_reg_unordered_extern,
    num_reg_ordered_extern,
    num_categories_extern,
    scale_factor_lower_bound);
  bwm_reset_counters();
  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  if (penalty_mode[0] == 1) {
    double pmult = penalty_mult[0];
    double baseline;
    if (pmult < 1.0) pmult = 1.0;
    baseline = bwmfunc_raw_current_scale(vector_scale_factor, num_all_var);
    /* Penalty-baseline probes are real objective evaluations and must be
       included in the evaluation/invalid accounting. */
    bwm_eval_count += 1.0;
    if (!R_FINITE(baseline) || baseline == DBL_MAX)
      bwm_invalid_count += 1.0;
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      double *tmp = alloc_vecd(num_all_var + 1);
      np_copy_scale_factor_for_raw(tmp, vector_scale_factor, num_all_var);
      for (i = 1; i <= (num_var_continuous_extern + num_reg_continuous_extern); i++)
        tmp[i] *= 2.0;
      for (i = 0; i < (num_var_unordered_extern + num_reg_unordered_extern); i++) {
        int idx = num_var_continuous_extern + num_reg_continuous_extern + 1 + i;
        double maxbw = max_unordered_bw(num_categories_extern[i], KERNEL_den_unordered_extern);
        tmp[idx] = 0.5*maxbw;
      }
      for (i = 0; i < (num_var_ordered_extern + num_reg_ordered_extern); i++) {
        int idx = num_var_continuous_extern + num_reg_continuous_extern +
          num_var_unordered_extern + num_reg_unordered_extern + 1 + i;
        tmp[idx] = 0.5;
      }
      baseline = bwmfunc_raw(tmp);
      bwm_eval_count += 1.0;
      if (!R_FINITE(baseline) || baseline == DBL_MAX)
        bwm_invalid_count += 1.0;
      safe_free(tmp);
    }
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      bwm_penalty_value = pmult * 1.0e6;
    } else {
      bwm_penalty_value = baseline + (fabs(baseline) + 1.0) * pmult;
    }
    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }

  fret_initial = fret_best = bwmfunc_wrapper(vector_scale_factor);
  iImproved = 0;
  have_start_best = 0;
  have_multistart_best = 0;
  fret_start_best = DBL_MAX;
  if (enforce_fixed_feasibility &&
      np_bw_candidate_is_admissible_with_floor(
        num_all_var,
        bwm_use_transform,
        KERNEL_den_extern,
        KERNEL_reg_unordered_extern,
        BANDWIDTH_den_extern,
        BANDWIDTH_den_extern,
        0,
        num_obs_train_extern,
        num_var_continuous_extern,
        num_var_unordered_extern,
        num_var_ordered_extern,
        num_reg_continuous_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_categories_extern,
        vector_scale_factor,
        scale_factor_lower_bound)) {
    have_start_best = 1;
    fret_start_best = fret_initial;
    np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor, num_all_var);
  }

  if(!eval_only){
    powell(0,
           0,
           vector_scale_factor,
           vector_scale_factor,
           matrix_y,
           num_all_var,
           ftol,
           tol,
           small,
           itmax,
           &iter,
           &fret,
           bwmfunc_wrapper);

    if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

      initialize_nr_directions(BANDWIDTH_den_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               num_var_continuous_extern,
                               num_var_unordered_extern,
                               num_var_ordered_extern,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               0, int_RANDOM_SEED,  
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);


      powell(0,
             0,
             vector_scale_factor,
             vector_scale_factor,
             matrix_y,
             num_all_var,
             ftol,
             tol,
             small,
             itmax,
             &iter,
             &fret,
             bwmfunc_wrapper);

    }
  } else {
    fret = fret_best;
  }

  if (enforce_fixed_feasibility &&
      np_bw_candidate_is_admissible_with_floor(
        num_all_var,
        bwm_use_transform,
        KERNEL_den_extern,
        KERNEL_reg_unordered_extern,
        BANDWIDTH_den_extern,
        BANDWIDTH_den_extern,
        0,
        num_obs_train_extern,
        num_var_continuous_extern,
        num_var_unordered_extern,
        num_var_ordered_extern,
        num_reg_continuous_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_categories_extern,
        vector_scale_factor,
        scale_factor_lower_bound) &&
      ((!have_start_best) || (fret < fret_start_best))) {
    have_start_best = 1;
    fret_start_best = fret;
    np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor, num_all_var);
  }

  if (np_bwm_get_deferred_error() != NULL) {
    bw_error_msg = np_bwm_get_deferred_error();
    goto cleanup_np_density_conditional_bw;
  }

  if (enforce_fixed_feasibility) {
    if (have_start_best) {
      fret = fret_start_best;
      np_copy_scale_factor(vector_scale_factor, vector_scale_factor_startbest, num_all_var);
    } else {
      fret = DBL_MAX;
    }
  }

  iImproved = (enforce_fixed_feasibility && have_start_best) ? (fret_start_best < fret_initial) : (fret < fret_best);
  *timing = timing_extern;

  objective_function_values[0]=-fret;
  objective_function_evals[0]=bwm_eval_count;
  objective_function_invalid[0]=bwm_invalid_count;
  bwm_snapshot_fast_counters();
  fast_eval_total += bwm_fast_eval_count;
  /* When multistarting save initial minimum of objective function and scale factors */


  if((!eval_only) && (iMultistart == IMULTI_TRUE)){
    if (enforce_fixed_feasibility) {
      if (have_start_best) {
        have_multistart_best = 1;
        fret_best = fret;
      } else {
        have_multistart_best = 0;
        fret_best = DBL_MAX;
      }
    } else {
      fret_best = fret;
    }
    vector_scale_factor_multistart = alloc_vecd(num_all_var + 1);
    for(i = 1; i <= num_all_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
    np_progress_bandwidth_multistart_step(1, iNum_Multistart);
			

    /* Conduct search from new random values of the search parameters */
		
    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and directions for NR modules */
      initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                        1,                /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        num_var_continuous_extern,
                                        num_var_unordered_extern,
                                        num_var_ordered_extern,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        KERNEL_den_unordered_extern,
                                        KERNEL_reg_unordered_extern,
                                        0,
                                        scale_cat,
                                        pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                        nconfac_extern, ncatfac_extern,
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor,
                                        lbc_init, hbc_init, c_init, 
                                        lbd_init, hbd_init, d_init,
                                        matrix_X_continuous_train_extern,
                                        matrix_Y_continuous_train_extern);

      initialize_nr_directions(BANDWIDTH_den_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               num_var_continuous_extern,
                               num_var_unordered_extern,
                               num_var_ordered_extern,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               1, int_RANDOM_SEED,  
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);

      /* Conduct direction set search */

      if (bwm_use_transform)
        bwm_to_unconstrained(vector_scale_factor, num_all_var);

      bwm_reset_counters();
      
      powell(0,
             0,
             vector_scale_factor,
             vector_scale_factor,
             matrix_y,
             num_all_var,
             ftol,
             tol,
             small,
             itmax,
             &iter,
             &fret,
             bwmfunc_wrapper);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

        initialize_nr_directions(BANDWIDTH_den_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 num_reg_unordered_extern,
                                 num_reg_ordered_extern,
                                 num_var_continuous_extern,
                                 num_var_unordered_extern,
                                 num_var_ordered_extern,
                                 vsfh,
                                 num_categories_extern,
                                 matrix_y,
                                 0, int_RANDOM_SEED,  
                                 lbc_dir, dfc_dir, c_dir, initc_dir,
                                 lbd_dir, hbd_dir, d_dir, initd_dir,
                                 matrix_X_continuous_train_extern,
                                 matrix_Y_continuous_train_extern);

        powell(0,
               0,
               vector_scale_factor,
               vector_scale_factor,
               matrix_y,
               num_all_var,
               ftol,
               tol,
               small,
               itmax,
               &iter,
               &fret,
               bwmfunc_wrapper);
      }
				
      /* If this run resulted in an improved minimum save information */
      
      if (enforce_fixed_feasibility) {
        if (np_bw_candidate_is_admissible_with_floor(
              num_all_var,
              bwm_use_transform,
              KERNEL_den_extern,
              KERNEL_reg_unordered_extern,
              BANDWIDTH_den_extern,
              BANDWIDTH_den_extern,
              0,
              num_obs_train_extern,
              num_var_continuous_extern,
              num_var_unordered_extern,
              num_var_ordered_extern,
              num_reg_continuous_extern,
              num_reg_unordered_extern,
              num_reg_ordered_extern,
              num_categories_extern,
              vector_scale_factor,
              scale_factor_lower_bound) &&
            ((!have_multistart_best) || (fret < fret_best))) {
          fret_best = fret;
          have_multistart_best = 1;
          iImproved = iMs_counter+1;
          *timing = timing_extern;
          np_copy_scale_factor(vector_scale_factor_multistart, vector_scale_factor, num_all_var);
        }
      } else if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        *timing = timing_extern;
        
        for(i = 1; i <= num_all_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }
      objective_function_values[iMs_counter]=-fret;
      objective_function_evals[iMs_counter]=bwm_eval_count;
      objective_function_invalid[iMs_counter]=bwm_invalid_count;
      bwm_snapshot_fast_counters();
      fast_eval_total += bwm_fast_eval_count;
      np_progress_bandwidth_multistart_step(iMs_counter+1, iNum_Multistart);
    }

    /* Save best for estimation */

    if (enforce_fixed_feasibility) {
      if (have_multistart_best) {
        fret = fret_best;
        np_copy_scale_factor(vector_scale_factor, vector_scale_factor_multistart, num_all_var);
        have_start_best = 1;
        fret_start_best = fret_best;
        np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor_multistart, num_all_var);
      } else {
        have_start_best = 0;
        fret = DBL_MAX;
      }
    } else {
      fret = fret_best;
      for(i = 1; i <= num_all_var; i++)
        vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];
    }
    free(vector_scale_factor_multistart);
  }

  if (np_bwm_get_deferred_error() != NULL) {
    bw_error_msg = np_bwm_get_deferred_error();
    goto cleanup_np_density_conditional_bw;
  }

  if (enforce_fixed_feasibility) {
    double final_raw;
    if (!have_start_best) {
      bw_error_msg = "C_np_density_conditional_bw: optimizer failed to produce a feasible fixed-bandwidth candidate";
      goto cleanup_np_density_conditional_bw;
    }
    if (!np_bw_candidate_is_admissible_with_floor(
          num_all_var,
          bwm_use_transform,
          KERNEL_den_extern,
          KERNEL_reg_unordered_extern,
          BANDWIDTH_den_extern,
          BANDWIDTH_den_extern,
          0,
          num_obs_train_extern,
          num_var_continuous_extern,
          num_var_unordered_extern,
          num_var_ordered_extern,
          num_reg_continuous_extern,
          num_reg_unordered_extern,
          num_reg_ordered_extern,
          num_categories_extern,
          vector_scale_factor,
          scale_factor_lower_bound)) {
      bw_error_msg = "C_np_density_conditional_bw: optimizer returned an infeasible fixed-bandwidth candidate";
      goto cleanup_np_density_conditional_bw;
    }
    final_raw = bwmfunc_raw_current_scale(vector_scale_factor, num_all_var);
    if (!R_FINITE(final_raw) || final_raw == DBL_MAX) {
      bw_error_msg = "C_np_density_conditional_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective";
      goto cleanup_np_density_conditional_bw;
    }
    fret = final_raw;
  }

  if (bwm_use_transform)
    bwm_to_constrained(vector_scale_factor, num_all_var);

  /* return data to R */
  if (BANDWIDTH_den_extern == BW_GEN_NN || 
      BANDWIDTH_den_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern+num_var_continuous_extern; i++ )
      vector_scale_factor[i+1]=np_fround(vector_scale_factor[i+1]);
  }

  for( i=0; i<num_all_var; i++ )
    myans[i]=vector_scale_factor[i+1];

  fval[0] = -fret;
  fval[1] = iImproved;
  fval[2] = -fret_initial;
  objective_function_fast[0] = fast_eval_total;
  /* end return data */

cleanup_np_density_conditional_bw:
  /* Free data objects */

  np_bounded_cvls_conditional_quad_context_clear_extern();
  bwm_clear_floor_context();

  free_mat(matrix_Y_unordered_train_extern, num_var_unordered_extern);
  free_mat(matrix_Y_ordered_train_extern, num_var_ordered_extern);
  free_mat(matrix_Y_continuous_train_extern, num_var_continuous_extern);

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  free_mat(matrix_y, num_all_var + 1);
  safe_free(vector_scale_factor);
  safe_free(vector_scale_factor_startbest);
  safe_free(vsfh);
  safe_free(num_categories_extern);
  safe_free(num_categories_extern_XY);
  safe_free(num_categories_extern_X);
  safe_free(num_categories_extern_Y);
  safe_free(bwm_kernel_unordered_vec);
  bwm_kernel_unordered_vec = NULL;
  bwm_kernel_unordered_len = 0;

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);

  free_mat(matrix_categorical_vals_extern_X, num_reg_unordered_extern + num_reg_ordered_extern);

  free_mat(matrix_categorical_vals_extern_XY, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);

  if(need_y_side)
    free_mat(matrix_categorical_vals_extern_Y, num_var_unordered_extern + num_var_ordered_extern);
  matrix_categorical_vals_extern_Y = NULL;

  safe_free(vector_continuous_stddev);

  safe_free(ipt_X);
  safe_free(ipt_XY);

  safe_free(ipt_lookup_X);
  safe_free(ipt_lookup_XY);

  if(need_y_side){
    safe_free(ipt_Y);
    safe_free(ipt_lookup_Y);
  }
  num_categories_extern_Y = NULL;

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  if(int_TREE_Y == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_Y);
    int_TREE_Y = NP_TREE_FALSE;
  }

  free_mat(matrix_XY_continuous_train_extern, num_all_cvar);
  free_mat(matrix_XY_unordered_train_extern, num_all_uvar);
  free_mat(matrix_XY_ordered_train_extern, num_all_ovar);

  if(int_TREE_XY == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_XY);
    int_TREE_XY = NP_TREE_FALSE;
  }


  int_WEIGHTS = 0;

  int_cxker_bound_extern = 0;
  int_cyker_bound_extern = 0;
  int_cxyker_bound_extern = 0;
  int_ll_extern = LL_LC;
  vector_glp_degree_extern = NULL;
  int_glp_bernstein_extern = 0;
  int_glp_basis_extern = 1;
  vector_cxkerlb_extern = NULL;
  vector_cxkerub_extern = NULL;
  vector_cykerlb_extern = NULL;
  vector_cykerub_extern = NULL;
  vector_cxykerlb_extern = NULL;
  vector_cxykerub_extern = NULL;
  int_cker_bound_extern = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;
  safe_free(cxylb);
  safe_free(cxyub);
  np_clear_estimator_extern_aliases();

  if (bw_error_msg != NULL) {
    np_bwm_clear_deferred_error();
    error("%s", bw_error_msg);
  }

  return ;
}

void np_distribution_conditional_bw(double * c_uno, double * c_ord, double * c_con, 
                                    double * u_uno, double * u_ord, double * u_con,
                                    double * cg_uno, double * cg_ord, double * cg_con, double * mysd,
                                    int * myopti, double * myoptd, double * myans, double * fval,
                                    double * objective_function_values, double * objective_function_evals,
                                    double * objective_function_invalid, double * timing,
                                    double * objective_function_fast,
                                    int * penalty_mode, double * penalty_mult,
                                    int * glp_degree,
                                    int * glp_bernstein,
                                    int * glp_basis,
                                    int * regtype,
                                    double * cxkerlb, double * cxkerub,
                                    double * cykerlb, double * cykerub,
                                    const int eval_only){
  int_nn_k_min_extern = 1;
/* Likelihood bandwidth selection for density estimation */

  double **matrix_y;

  double *vector_continuous_stddev;
  double *vsfh, *vector_scale_factor, *vector_scale_factor_multistart;
  double *vector_scale_factor_startbest;
  double *cxylb = NULL, *cxyub = NULL;

  double fret, fret_best, fret_start_best, fret_initial;
  double ftol, tol;
  double (* bwmfunc)(double *) = NULL;

  double small, lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init; 
  double lbd_init, hbd_init, d_init;
  int dfc_dir;
  
  int i,j;
  double fast_eval_total = 0.0;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, num_all_var, num_var_var, iImproved;
  int enforce_fixed_feasibility;
  int have_start_best, have_multistart_best;
  int itmax, iter;
  int int_use_starting_values, ibwmfunc;
  int cdfontrain;
  int scale_cat;
  const char *bw_error_msg = NULL;

  int num_all_cvar, num_all_uvar, num_all_ovar;

  int * ipt_X = NULL, * ipt_XY = NULL, * ipt_Y = NULL; 
  int * ipt_lookup_XY = NULL, * ipt_lookup_Y = NULL, * ipt_lookup_X = NULL;
  int num_obs_alt;

  cdfontrain_extern = cdfontrain =  myopti[CDBW_CDFONTRAIN];

  num_var_unordered_extern = myopti[CDBW_CNUNOI];
  num_var_ordered_extern = myopti[CDBW_CNORDI];
  num_var_continuous_extern = myopti[CDBW_CNCONI];

  num_reg_unordered_extern = myopti[CDBW_UNUNOI];
  num_reg_ordered_extern = myopti[CDBW_UNORDI];
  num_reg_continuous_extern = myopti[CDBW_UNCONI];
  bwm_clear_floor_context();

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;
  num_var_var = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  num_all_var = num_var+num_var_var;

  num_obs_train_extern = myopti[CDBW_NOBSI];
  num_obs_eval_extern  = cdfontrain ? num_obs_train_extern : myopti[CDBW_NEVALI];

  iMultistart = myopti[CDBW_IMULTII];
  iNum_Multistart = myopti[CDBW_NMULTII];
  if (!eval_only && iNum_Multistart < 1)
    error("C_np_distribution_conditional_bw: nmulti must be a positive integer");

  KERNEL_reg_extern = myopti[CDBW_CXKRNEVI];
  KERNEL_den_extern = myopti[CDBW_CYKRNEVI];

  KERNEL_reg_unordered_extern = myopti[CDBW_UXKRNEVI];
  KERNEL_den_unordered_extern = myopti[CDBW_UYKRNEVI];

  KERNEL_reg_ordered_extern = myopti[CDBW_OXKRNEVI];
  KERNEL_den_ordered_extern = myopti[CDBW_OYKRNEVI];

  vector_cxkerlb_extern = cxkerlb;
  vector_cxkerub_extern = cxkerub;
  int_cxker_bound_extern = np_has_finite_cker_bounds(cxkerlb, cxkerub, num_reg_continuous_extern);

  vector_cykerlb_extern = cykerlb;
  vector_cykerub_extern = cykerub;
  int_cyker_bound_extern = np_has_finite_cker_bounds(cykerlb, cykerub, num_var_continuous_extern);

  if((num_reg_continuous_extern + num_var_continuous_extern) > 0){
    cxylb = alloc_vecd(num_reg_continuous_extern + num_var_continuous_extern);
    cxyub = alloc_vecd(num_reg_continuous_extern + num_var_continuous_extern);
    for(i = 0; i < num_reg_continuous_extern; i++){
      cxylb[i] = (cxkerlb != NULL) ? cxkerlb[i] : DBL_MAX;
      cxyub[i] = (cxkerub != NULL) ? cxkerub[i] : DBL_MAX;
    }
    for(i = 0; i < num_var_continuous_extern; i++){
      cxylb[num_reg_continuous_extern + i] = (cykerlb != NULL) ? cykerlb[i] : DBL_MAX;
      cxyub[num_reg_continuous_extern + i] = (cykerub != NULL) ? cykerub[i] : DBL_MAX;
    }
    vector_cxykerlb_extern = cxylb;
    vector_cxykerub_extern = cxyub;
    int_cxyker_bound_extern = np_has_finite_cker_bounds(cxylb, cxyub, num_reg_continuous_extern + num_var_continuous_extern);
  } else {
    vector_cxykerlb_extern = NULL;
    vector_cxykerub_extern = NULL;
    int_cxyker_bound_extern = 0;
  }

  int_use_starting_values= myopti[CDBW_USTARTI];
  int_LARGE_SF=myopti[CDBW_LSFI];
  BANDWIDTH_den_extern=myopti[CDBW_DENI];
  enforce_fixed_feasibility = ((BANDWIDTH_den_extern == BW_FIXED) && (!eval_only));
  int_RESTART_FROM_MIN = myopti[CDBW_REMINI];
  int_MINIMIZE_IO = myopti[CDBW_MINIOI];

  itmax=myopti[CDBW_ITMAXI];
  ibwmfunc = myopti[CDBW_MI];
  int_ll_extern = (ibwmfunc == CDBWM_CVLS) ? *regtype : LL_LC;
  vector_glp_degree_extern = ((ibwmfunc == CDBWM_CVLS) && (int_ll_extern == LL_LP)) ? glp_degree : NULL;
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = ((ibwmfunc == CDBWM_CVLS) && (int_ll_extern == LL_LP)) ? *glp_bernstein : 0;
  int_glp_basis_extern = ((ibwmfunc == CDBWM_CVLS) && (int_ll_extern == LL_LP)) ? *glp_basis : 1;

  int_TREE_XY = int_TREE_Y = int_TREE_X = myopti[CDBW_TREEI];
  if(int_ll_extern == LL_LP){
    int_TREE_Y = NP_TREE_FALSE;
    int_TREE_XY = NP_TREE_FALSE;
  }

  scale_cat = myopti[CDBW_SCATI];
  bwm_use_transform = myopti[CDBW_TBNDI];
  if (BANDWIDTH_den_extern != BW_FIXED)
    bwm_use_transform = 0;
  if (bwm_use_transform) {
    int n = num_var_continuous_extern + num_reg_continuous_extern +
      num_var_unordered_extern + num_reg_unordered_extern +
      num_var_ordered_extern + num_reg_ordered_extern;
    bwm_reserve_transform_buf(n + 1);
  }

  ftol=myoptd[CDBW_FTOLD];
  tol=myoptd[CDBW_TOLD];
  small=myoptd[CDBW_SMALLD];
  dbl_memfac_ccdf_extern = myoptd[CDBW_MEMFACD];

  dfc_dir = myopti[CDBW_DFC_DIRI];
  lbc_dir = myoptd[CDBW_LBC_DIRD];
  c_dir = myoptd[CDBW_C_DIRD];
  initc_dir = myoptd[CDBW_INITC_DIRD]; 

  lbd_dir = myoptd[CDBW_LBD_DIRD]; 
  hbd_dir = myoptd[CDBW_HBD_DIRD]; 
  d_dir = myoptd[CDBW_D_DIRD]; 
  initd_dir = myoptd[CDBW_INITD_DIRD]; 

  lbc_init = myoptd[CDBW_LBC_INITD]; 
  hbc_init = myoptd[CDBW_HBC_INITD]; 
  c_init = myoptd[CDBW_C_INITD]; 

  lbd_init = myoptd[CDBW_LBD_INITD]; 
  hbd_init = myoptd[CDBW_HBD_INITD]; 
  d_init = myoptd[CDBW_D_INITD]; 

  nconfac_extern = myoptd[CDBW_NCONFD];
  ncatfac_extern = myoptd[CDBW_NCATFD];
  bwm_set_scale_factor_lower_bound(myoptd[CDBW_SFLOORD]);

/* Allocate memory for objects */

  matrix_Y_unordered_train_extern = alloc_matd(num_obs_train_extern, num_var_unordered_extern);
  matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);
  matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  if(cdfontrain){
    matrix_Y_unordered_eval_extern = matrix_Y_unordered_train_extern;
    matrix_Y_ordered_eval_extern = matrix_Y_ordered_train_extern;
    matrix_Y_continuous_eval_extern = matrix_Y_continuous_train_extern;
  } else {
    matrix_Y_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_var_unordered_extern);
    matrix_Y_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_var_ordered_extern);
    matrix_Y_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_var_continuous_extern);
  }

  num_categories_extern = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                     num_reg_unordered_extern + num_reg_ordered_extern);

  num_categories_extern_X = alloc_vecu(num_reg_unordered_extern + num_reg_ordered_extern);
  num_categories_extern_Y = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern);

  num_categories_extern_XY = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                        num_reg_unordered_extern + num_reg_ordered_extern);
  
  matrix_y = alloc_matd(num_all_var + 1, num_all_var + 1);
  vector_scale_factor = alloc_vecd(num_all_var + 1);
  vector_scale_factor_startbest = alloc_vecd(num_all_var + 1);
  vsfh = alloc_vecd(num_all_var + 1);
  
  matrix_categorical_vals_extern = 
    alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern + 
               num_reg_unordered_extern + num_reg_ordered_extern);

  matrix_categorical_vals_extern_X = 
    alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  matrix_categorical_vals_extern_Y = 
    alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern);

  matrix_categorical_vals_extern_XY = 
    alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern + 
               num_reg_unordered_extern + num_reg_ordered_extern);

  /* in v_s_f order is creg, cvar, uvar, ovar, ureg, oreg  */

  if (int_use_starting_values)
    for( i=0;i<num_all_var; i++ )
      vector_scale_factor[i+1] = myans[i];

/* Parse data */

  for(j=0;j<num_var_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_unordered_train_extern[j][i]=c_uno[j*num_obs_train_extern+i];

  for(j=0;j<num_var_ordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_ordered_train_extern[j][i]=c_ord[j*num_obs_train_extern+i];

  for(j=0;j<num_var_continuous_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_continuous_train_extern[j][i]=c_con[j*num_obs_train_extern+i];


  for(j=0;j<num_reg_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_X_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+i];

  if(!cdfontrain){
    for(j=0;j<num_var_unordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_Y_unordered_eval_extern[j][i]=cg_uno[j*num_obs_eval_extern+i];

    for(j=0;j<num_var_ordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_Y_ordered_eval_extern[j][i]=cg_ord[j*num_obs_eval_extern+i];

    for(j=0;j<num_var_continuous_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_Y_continuous_eval_extern[j][i]=cg_con[j*num_obs_eval_extern+i];

  }

  num_all_cvar = num_reg_continuous_extern + num_var_continuous_extern;
  num_all_uvar = num_reg_unordered_extern + num_var_unordered_extern;
  num_all_ovar = num_reg_ordered_extern + num_var_ordered_extern;

  // we need 3 trees to accelerate cg_concdf :)

  ipt_X = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_X != NULL))
    error("!(ipt_X != NULL)");

  ipt_lookup_X = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_lookup_X != NULL)){
    safe_free(ipt_X);
    error("!(ipt_lookup_X != NULL)");
  }

  for(i = 0; i < num_obs_train_extern; i++){
    ipt_lookup_X[i] = ipt_X[i] = i;
  }

  ipt_extern_X = ipt_X;
  ipt_lookup_extern_X = ipt_lookup_X;


  ipt_Y = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_Y != NULL)){
    safe_free(ipt_X);
    safe_free(ipt_lookup_X);
    error("!(ipt_Y != NULL)");
  }

  ipt_lookup_Y = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_lookup_Y != NULL)){
    safe_free(ipt_X);
    safe_free(ipt_lookup_X);
    safe_free(ipt_Y);
    error("!(ipt_lookup_Y != NULL)");
  }

  for(i = 0; i < num_obs_train_extern; i++){
    ipt_lookup_Y[i] = ipt_Y[i] = i;
  }

  ipt_extern_Y = ipt_Y;
  ipt_lookup_extern_Y = ipt_lookup_Y;

  num_obs_alt = (BANDWIDTH_den_extern != BW_ADAP_NN) ? num_obs_train_extern : 0;

  ipt_XY = (int *)malloc(num_obs_alt*sizeof(int));
  if(!(ipt_XY != NULL)){
    safe_free(ipt_X);
    safe_free(ipt_lookup_X);
    safe_free(ipt_Y);
    safe_free(ipt_lookup_Y);
    error("!(ipt_XY != NULL)");
  }

  ipt_lookup_XY = (int *)malloc(num_obs_alt*sizeof(int));
  if(!(ipt_lookup_XY != NULL)){
    safe_free(ipt_X);
    safe_free(ipt_lookup_X);
    safe_free(ipt_Y);
    safe_free(ipt_lookup_Y);
    safe_free(ipt_XY);
    error("!(ipt_lookup_XY != NULL)");
  }

  for(i = 0; i < num_obs_alt; i++){
    ipt_lookup_XY[i] = ipt_XY[i] = i;
  }

  ipt_extern_XY = ipt_XY;
  ipt_lookup_extern_XY = ipt_lookup_XY;

  int_TREE_XY = int_TREE_XY && (((num_all_cvar) != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);

  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);

  int_TREE_Y = int_TREE_Y && ((num_var_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);

  if(int_TREE_X == NP_TREE_TRUE){
    build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                 4*num_reg_continuous_extern, ipt_X, &kdt_extern_X);
  
    // put x data into x-tree order
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+ipt_X[i]];
    
    
    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+ipt_X[i]];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+ipt_X[i]];

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_X[ipt_X[i]] = i;
    }

  }

  if(int_TREE_Y == NP_TREE_TRUE){
    build_kdtree(matrix_Y_continuous_train_extern, num_obs_train_extern, num_var_continuous_extern, 
                 4*num_var_continuous_extern, ipt_Y, &kdt_extern_Y);
  
    // put y data into y-tree order

    for(j=0;j<num_var_unordered_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_unordered_train_extern[j][i]=c_uno[j*num_obs_train_extern+ipt_Y[i]];

    for(j=0;j<num_var_ordered_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_ordered_train_extern[j][i]=c_ord[j*num_obs_train_extern+ipt_Y[i]];

    for(j=0;j<num_var_continuous_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_continuous_train_extern[j][i]=c_con[j*num_obs_train_extern+ipt_Y[i]];

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_Y[ipt_Y[i]] = i;
    }
    
  }

  matrix_XY_continuous_train_extern = alloc_matd(num_obs_alt, num_all_cvar);
  matrix_XY_unordered_train_extern = alloc_matd(num_obs_alt, num_all_uvar);
  matrix_XY_ordered_train_extern = alloc_matd(num_obs_alt, num_all_ovar);

  if(int_TREE_XY == NP_TREE_TRUE){

    for(j = 0; j < num_reg_unordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+i];

    for(j = num_reg_unordered_extern; j < num_all_uvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=c_uno[(j-num_reg_unordered_extern)*num_obs_train_extern+i];


    for(j = 0; j < num_reg_ordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+i];

    for(j = num_reg_ordered_extern; j < num_all_ovar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=c_ord[(j-num_reg_ordered_extern)*num_obs_train_extern+i];


    for(j = 0; j < num_reg_continuous_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+i];

    for(j = num_reg_continuous_extern; j < num_all_cvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=c_con[(j-num_reg_continuous_extern)*num_obs_train_extern+i];

  // XY tree!

    build_kdtree(matrix_XY_continuous_train_extern, num_obs_train_extern, num_all_cvar, 
                 4*num_all_cvar, ipt_XY, &kdt_extern_XY);

    // put data into xy-tree order
    for(j = 0; j < num_reg_unordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_unordered_extern; j < num_all_uvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=c_uno[(j-num_reg_unordered_extern)*num_obs_train_extern+ipt_XY[i]];


    for(j = 0; j < num_reg_ordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_ordered_extern; j < num_all_ovar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=c_ord[(j-num_reg_ordered_extern)*num_obs_train_extern+ipt_XY[i]];


    for(j = 0; j < num_reg_continuous_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_continuous_extern; j < num_all_cvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=c_con[(j-num_reg_continuous_extern)*num_obs_train_extern+ipt_XY[i]];

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_XY[ipt_XY[i]] = i;
    }
    
  }

  determine_categorical_vals(
                             num_obs_train_extern,
                             num_var_unordered_extern,
                             num_var_ordered_extern,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             matrix_Y_unordered_train_extern,
                             matrix_Y_ordered_train_extern,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern,
                             matrix_categorical_vals_extern);

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern, num_var_ordered_extern, num_var_continuous_extern,
                        num_reg_unordered_extern, num_reg_ordered_extern, num_reg_continuous_extern,
                        vector_scale_factor+1,
                        num_categories_extern,
                        matrix_categorical_vals_extern,
                        NULL, NULL, NULL,
                        num_categories_extern_X, num_categories_extern_Y, num_categories_extern_XY,
                        matrix_categorical_vals_extern_X, matrix_categorical_vals_extern_Y, matrix_categorical_vals_extern_XY);


  vector_continuous_stddev = vector_continuous_stddev_extern = mysd;


  /* Initialize scale factors and Directions for NR modules */

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    num_var_continuous_extern,
                                    num_var_unordered_extern,
                                    num_var_ordered_extern,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    KERNEL_den_unordered_extern,
                                    KERNEL_reg_unordered_extern,
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    num_var_continuous_extern,
                                    num_var_unordered_extern,
                                    num_var_ordered_extern,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    KERNEL_den_unordered_extern,
                                    KERNEL_reg_unordered_extern,
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vsfh,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           num_var_continuous_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED,  
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);


  /* When multistarting, set counter */

  imsnum = iMs_counter = 0;
  imstot = iNum_Multistart;

  /* Conduct direction set search */

  /* assign the function to be optimized */

  ibwmfunc = myopti[CDBW_MI];

  switch(ibwmfunc){
  case CDBWM_CVLS : bwmfunc = cv_func_con_distribution_categorical_ls; break;
  default : REprintf("np.c: invalid bandwidth selection method.");
    error("np.c: invalid bandwidth selection method."); break;
  }

  if (bwm_use_transform)
    bwm_to_unconstrained(vector_scale_factor, num_all_var);

  spinner(0);

  bwmfunc_raw = bwmfunc;
  bwm_num_reg_continuous = num_var_continuous_extern + num_reg_continuous_extern;
  bwm_num_reg_unordered = num_var_unordered_extern + num_reg_unordered_extern;
  bwm_num_reg_ordered = num_var_ordered_extern + num_reg_ordered_extern;
  bwm_kernel_unordered = KERNEL_den_unordered_extern;
  bwm_kernel_unordered_len = bwm_num_reg_unordered;
  if (bwm_kernel_unordered_len > 0) {
    bwm_kernel_unordered_vec = alloc_vecu(bwm_kernel_unordered_len);
    for (i = 0; i < num_var_unordered_extern; i++)
      bwm_kernel_unordered_vec[i] = KERNEL_den_unordered_extern;
    for (i = 0; i < num_reg_unordered_extern; i++)
      bwm_kernel_unordered_vec[num_var_unordered_extern + i] = KERNEL_reg_unordered_extern;
  } else {
    bwm_kernel_unordered_vec = NULL;
  }
  bwm_num_categories = num_categories_extern;
  bwm_set_floor_context(
    enforce_fixed_feasibility,
    num_all_var,
    bwm_use_transform,
    KERNEL_den_extern,
    KERNEL_reg_unordered_extern,
    BANDWIDTH_den_extern,
    BANDWIDTH_den_extern,
    0,
    num_obs_train_extern,
    num_var_continuous_extern,
    num_var_unordered_extern,
    num_var_ordered_extern,
    num_reg_continuous_extern,
    num_reg_unordered_extern,
    num_reg_ordered_extern,
    num_categories_extern,
    bwm_scale_factor_lower_bound);
  bwm_reset_counters();
  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  if (penalty_mode[0] == 1) {
    double pmult = penalty_mult[0];
    double baseline;
    if (pmult < 1.0) pmult = 1.0;
    baseline = bwmfunc_raw_current_scale(vector_scale_factor, num_all_var);
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      double *tmp = alloc_vecd(num_all_var + 1);
      np_copy_scale_factor_for_raw(tmp, vector_scale_factor, num_all_var);
      for (i = 1; i <= (num_var_continuous_extern + num_reg_continuous_extern); i++)
        tmp[i] *= 2.0;
      for (i = 0; i < (num_var_unordered_extern + num_reg_unordered_extern); i++) {
        int idx = num_var_continuous_extern + num_reg_continuous_extern + 1 + i;
        double maxbw = max_unordered_bw(num_categories_extern[i], KERNEL_den_unordered_extern);
        tmp[idx] = 0.5*maxbw;
      }
      for (i = 0; i < (num_var_ordered_extern + num_reg_ordered_extern); i++) {
        int idx = num_var_continuous_extern + num_reg_continuous_extern +
          num_var_unordered_extern + num_reg_unordered_extern + 1 + i;
        tmp[idx] = 0.5;
      }
      baseline = bwmfunc_raw(tmp);
      safe_free(tmp);
    }
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      bwm_penalty_value = pmult * 1.0e6;
    } else {
      bwm_penalty_value = baseline + (fabs(baseline) + 1.0) * pmult;
    }
    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }

  fret_initial = fret_best = bwmfunc_wrapper(vector_scale_factor);
  iImproved = 0;
  have_start_best = 0;
  have_multistart_best = 0;
  fret_start_best = DBL_MAX;
  if (enforce_fixed_feasibility &&
      np_bw_candidate_is_admissible(
        num_all_var,
        bwm_use_transform,
        KERNEL_den_extern,
        KERNEL_reg_unordered_extern,
        BANDWIDTH_den_extern,
        BANDWIDTH_den_extern,
        0,
        num_obs_train_extern,
        num_var_continuous_extern,
        num_var_unordered_extern,
        num_var_ordered_extern,
        num_reg_continuous_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_categories_extern,
        vector_scale_factor)) {
    have_start_best = 1;
    fret_start_best = fret_initial;
    np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor, num_all_var);
  }

  if(!eval_only){
    powell(0,
           0,
           vector_scale_factor,
           vector_scale_factor,
           matrix_y,
           num_all_var,
           ftol,
           tol,
           small,
           itmax,
           &iter,
           &fret,
           bwmfunc_wrapper);

    if(int_RESTART_FROM_MIN == RE_MIN_TRUE){
      initialize_nr_directions(BANDWIDTH_den_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               num_var_continuous_extern,
                               num_var_unordered_extern,
                               num_var_ordered_extern,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               0, int_RANDOM_SEED,  
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);

      powell(0,
             0,
             vector_scale_factor,
             vector_scale_factor,
             matrix_y,
             num_all_var,
             ftol,
             tol,
             small,
             itmax,
             &iter,
             &fret,
             bwmfunc_wrapper);

    }
  } else {
    fret = fret_best;
  }

  if (enforce_fixed_feasibility &&
      np_bw_candidate_is_admissible(
        num_all_var,
        bwm_use_transform,
        KERNEL_den_extern,
        KERNEL_reg_unordered_extern,
        BANDWIDTH_den_extern,
        BANDWIDTH_den_extern,
        0,
        num_obs_train_extern,
        num_var_continuous_extern,
        num_var_unordered_extern,
        num_var_ordered_extern,
        num_reg_continuous_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_categories_extern,
        vector_scale_factor) &&
      ((!have_start_best) || (fret < fret_start_best))) {
    have_start_best = 1;
    fret_start_best = fret;
    np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor, num_all_var);
  }

  if (enforce_fixed_feasibility) {
    if (have_start_best) {
      fret = fret_start_best;
      np_copy_scale_factor(vector_scale_factor, vector_scale_factor_startbest, num_all_var);
    } else {
      fret = DBL_MAX;
    }
  }

  iImproved = (enforce_fixed_feasibility && have_start_best) ? (fret_start_best < fret_initial) : (fret < fret_best);
  *timing = timing_extern;

  objective_function_values[0]=fret;
  objective_function_evals[0]=bwm_eval_count;
  objective_function_invalid[0]=bwm_invalid_count;
  bwm_snapshot_fast_counters();
  fast_eval_total += bwm_fast_eval_count;
  /* When multistarting save initial minimum of objective function and scale factors */


  if((!eval_only) && (iMultistart == IMULTI_TRUE)){
    if (enforce_fixed_feasibility) {
      if (have_start_best) {
        have_multistart_best = 1;
        fret_best = fret;
      } else {
        have_multistart_best = 0;
        fret_best = DBL_MAX;
      }
    } else {
      fret_best = fret;
    }
    vector_scale_factor_multistart = alloc_vecd(num_all_var + 1);
    for(i = 1; i <= num_all_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
    np_progress_bandwidth_multistart_step(1, iNum_Multistart);
			

    /* Conduct search from new random values of the search parameters */
		
    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and directions for NR modules */
      initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                        1,                /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        num_var_continuous_extern,
                                        num_var_unordered_extern,
                                        num_var_ordered_extern,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        KERNEL_den_unordered_extern,
                                        KERNEL_reg_unordered_extern,
                                        0,
                                        scale_cat,
                                        pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                        nconfac_extern, ncatfac_extern,
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor,
                                        lbc_init, hbc_init, c_init, 
                                        lbd_init, hbd_init, d_init,
                                        matrix_X_continuous_train_extern,
                                        matrix_Y_continuous_train_extern);

      initialize_nr_directions(BANDWIDTH_den_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               num_var_continuous_extern,
                               num_var_unordered_extern,
                               num_var_ordered_extern,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               1, int_RANDOM_SEED,  
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);


      /* Conduct direction set search */

      bwm_reset_counters();

      powell(0,
             0,
             vector_scale_factor,
             vector_scale_factor,
             matrix_y,
             num_all_var,
             ftol,
             tol,
             small,
             itmax,
             &iter,
             &fret,
             bwmfunc_wrapper);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

        initialize_nr_directions(BANDWIDTH_den_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 num_reg_unordered_extern,
                                 num_reg_ordered_extern,
                                 num_var_continuous_extern,
                                 num_var_unordered_extern,
                                 num_var_ordered_extern,
                                 vsfh,
                                 num_categories_extern,
                                 matrix_y, 
                                 0, int_RANDOM_SEED,  
                                 lbc_dir, dfc_dir, c_dir, initc_dir,
                                 lbd_dir, hbd_dir, d_dir, initd_dir,
                                 matrix_X_continuous_train_extern,
                                 matrix_Y_continuous_train_extern);


        powell(0,
               0,
               vector_scale_factor,
               vector_scale_factor,
               matrix_y,
               num_all_var,
               ftol,
               tol,
               small,
               itmax,
               &iter,
               &fret,
               bwmfunc_wrapper);
      }
				
      /* If this run resulted in an improved minimum save information */
      
      if (enforce_fixed_feasibility) {
        if (np_bw_candidate_is_admissible(
              num_all_var,
              bwm_use_transform,
              KERNEL_den_extern,
              KERNEL_reg_unordered_extern,
              BANDWIDTH_den_extern,
              BANDWIDTH_den_extern,
              0,
              num_obs_train_extern,
              num_var_continuous_extern,
              num_var_unordered_extern,
              num_var_ordered_extern,
              num_reg_continuous_extern,
              num_reg_unordered_extern,
              num_reg_ordered_extern,
              num_categories_extern,
              vector_scale_factor) &&
            ((!have_multistart_best) || (fret < fret_best))) {
          fret_best = fret;
          have_multistart_best = 1;
          iImproved = iMs_counter+1;
          *timing = timing_extern;
          np_copy_scale_factor(vector_scale_factor_multistart, vector_scale_factor, num_all_var);
        }
      } else if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        *timing = timing_extern;
        
        for(i = 1; i <= num_all_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }
      objective_function_values[iMs_counter]=fret;
      objective_function_evals[iMs_counter]=bwm_eval_count;
      objective_function_invalid[iMs_counter]=bwm_invalid_count;
      bwm_snapshot_fast_counters();
      fast_eval_total += bwm_fast_eval_count;
      np_progress_bandwidth_multistart_step(iMs_counter+1, iNum_Multistart);
    }

    /* Save best for estimation */

    if (enforce_fixed_feasibility) {
      if (have_multistart_best) {
        fret = fret_best;
        np_copy_scale_factor(vector_scale_factor, vector_scale_factor_multistart, num_all_var);
        have_start_best = 1;
        fret_start_best = fret_best;
        np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor_multistart, num_all_var);
      } else {
        have_start_best = 0;
        fret = DBL_MAX;
      }
    } else {
      fret = fret_best;
      for(i = 1; i <= num_all_var; i++)
        vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];
    }
    free(vector_scale_factor_multistart);
  }

  if (enforce_fixed_feasibility) {
    double final_raw;
    if (!have_start_best) {
      bw_error_msg = "C_np_distribution_conditional_bw: optimizer failed to produce a feasible fixed-bandwidth candidate";
      goto cleanup_np_distribution_conditional_bw;
    }
    if (!np_bw_candidate_is_admissible(
          num_all_var,
          bwm_use_transform,
          KERNEL_den_extern,
          KERNEL_reg_unordered_extern,
          BANDWIDTH_den_extern,
          BANDWIDTH_den_extern,
          0,
          num_obs_train_extern,
          num_var_continuous_extern,
          num_var_unordered_extern,
          num_var_ordered_extern,
          num_reg_continuous_extern,
          num_reg_unordered_extern,
          num_reg_ordered_extern,
          num_categories_extern,
          vector_scale_factor)) {
      bw_error_msg = "C_np_distribution_conditional_bw: optimizer returned an infeasible fixed-bandwidth candidate";
      goto cleanup_np_distribution_conditional_bw;
    }
    final_raw = bwmfunc_raw_current_scale(vector_scale_factor, num_all_var);
    if (!R_FINITE(final_raw) || final_raw == DBL_MAX) {
      bw_error_msg = "C_np_distribution_conditional_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective";
      goto cleanup_np_distribution_conditional_bw;
    }
    fret = final_raw;
  }

  if (bwm_use_transform)
    bwm_to_constrained(vector_scale_factor, num_all_var);

  /* return data to R */
  if (BANDWIDTH_den_extern == BW_GEN_NN || 
      BANDWIDTH_den_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern+num_var_continuous_extern; i++ )
      vector_scale_factor[i+1]=np_fround(vector_scale_factor[i+1]);
  }

  for( i=0; i<num_all_var; i++ )
    myans[i]=vector_scale_factor[i+1];

  fval[0] = fret;
  fval[1] = iImproved;
  fval[2] = fret_initial;
  objective_function_fast[0] = fast_eval_total;
  /* end return data */

cleanup_np_distribution_conditional_bw:
  /* Free data objects */
  bwm_clear_floor_context();

  free_mat(matrix_Y_unordered_train_extern, num_var_unordered_extern);
  free_mat(matrix_Y_ordered_train_extern, num_var_ordered_extern);
  free_mat(matrix_Y_continuous_train_extern, num_var_continuous_extern);

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  if(!cdfontrain){
    free_mat(matrix_Y_unordered_eval_extern, num_var_unordered_extern);
    free_mat(matrix_Y_ordered_eval_extern, num_var_ordered_extern);
    free_mat(matrix_Y_continuous_eval_extern, num_var_continuous_extern);
  }

  free_mat(matrix_y, num_all_var + 1);
  safe_free(vector_scale_factor);
  safe_free(vector_scale_factor_startbest);
  safe_free(vsfh);
  safe_free(num_categories_extern);
  safe_free(bwm_kernel_unordered_vec);
  bwm_kernel_unordered_vec = NULL;
  bwm_kernel_unordered_len = 0;
  safe_free(num_categories_extern_X);
  safe_free(num_categories_extern_Y);
  safe_free(num_categories_extern_XY);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);

  free_mat(matrix_categorical_vals_extern_X, num_reg_unordered_extern + num_reg_ordered_extern);

  free_mat(matrix_categorical_vals_extern_Y, num_var_unordered_extern + num_var_ordered_extern);
  free_mat(matrix_categorical_vals_extern_XY, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);

  safe_free(ipt_X);
  safe_free(ipt_Y);
  safe_free(ipt_XY);

  safe_free(ipt_lookup_X);
  safe_free(ipt_lookup_Y);
  safe_free(ipt_lookup_XY);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

 if(int_TREE_Y == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_Y);
    int_TREE_Y = NP_TREE_FALSE;
  }

  free_mat(matrix_XY_continuous_train_extern, num_all_cvar);
  free_mat(matrix_XY_unordered_train_extern, num_all_uvar);
  free_mat(matrix_XY_ordered_train_extern, num_all_ovar);

  if(int_TREE_XY == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_XY);
    int_TREE_XY = NP_TREE_FALSE;
  }


  int_WEIGHTS = 0;

  int_cxker_bound_extern = 0;
  int_cyker_bound_extern = 0;
  int_cxyker_bound_extern = 0;
  int_ll_extern = LL_LC;
  vector_glp_degree_extern = NULL;
  int_glp_bernstein_extern = 0;
  int_glp_basis_extern = 1;
  vector_cxkerlb_extern = NULL;
  vector_cxkerub_extern = NULL;
  vector_cykerlb_extern = NULL;
  vector_cykerub_extern = NULL;
  vector_cxykerlb_extern = NULL;
  vector_cxykerub_extern = NULL;
  int_cker_bound_extern = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;
  safe_free(cxylb);
  safe_free(cxyub);
  np_clear_estimator_extern_aliases();

  if (bw_error_msg != NULL)
    error("%s", bw_error_msg);

  return ;
}


void np_density_conditional(double * tc_uno, double * tc_ord, double * tc_con, 
                            double * tu_uno, double * tu_ord, double * tu_con,
                            double * ec_uno, double * ec_ord, double * ec_con, 
                            double * eu_uno, double * eu_ord, double * eu_con,
                            double * mybw, 
                            double * ymcv, double * ypadnum,
                            double * xmcv, double * xpadnum,
                            double * nconfac, double * ncatfac, double * mysd,
                            int * myopti, 
                            double * cdens, double * cderr, 
                            double * cg, double * cgerr,
                            double * ll,
                            double * cxkerlb, double * cxkerub,
                            double * cykerlb, double * cykerub){
  /* Likelihood bandwidth selection for density estimation */

  double *vector_scale_factor, *pdf, *pdf_stderr, log_likelihood = 0.0;
  double ** pdf_deriv = NULL, ** pdf_deriv_stderr = NULL;
  double *cxylb = NULL, *cxyub = NULL;
  double xpad_num, ypad_num;

  int i,j;
  int num_var;

  int num_all_var, num_var_var, train_is_eval, do_grad, num_obs_eval_alloc;
  int num_all_cvar, num_all_uvar, num_all_ovar, num_all_catvar;
  int xmax_lev, ymax_lev, dens_or_dist, t_num;

  int * ipt_XY = NULL, *ipe_XY = NULL;
  int operator;

  int int_ll_eff;
  int ycat_offset;
  int xcat_offset;
  int saved_cker_bound;
  double * saved_ckerlb = NULL;
  double * saved_ckerub = NULL;


  num_var_unordered_extern = myopti[CD_CNUNOI];
  num_var_ordered_extern = myopti[CD_CNORDI];
  num_var_continuous_extern = myopti[CD_CNCONI];

  num_reg_unordered_extern = myopti[CD_UNUNOI];
  num_reg_ordered_extern = myopti[CD_UNORDI];
  num_reg_continuous_extern = myopti[CD_UNCONI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;
  num_var_var = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  num_all_var = num_var + num_var_var;

  num_all_cvar = num_reg_continuous_extern + num_var_continuous_extern;
  num_all_uvar = num_reg_unordered_extern + num_var_unordered_extern;
  num_all_ovar = num_reg_ordered_extern + num_var_ordered_extern;

  num_all_catvar = num_all_uvar + num_all_ovar;

  num_obs_train_extern = myopti[CD_TNOBSI];
  num_obs_eval_extern = myopti[CD_ENOBSI];

  if((train_is_eval = myopti[CD_TISEI]) && 
     (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_density_conditional): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    error("\n(np_density_conditional): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
  }

  KERNEL_reg_extern = myopti[CD_CXKRNEVI];
  KERNEL_den_extern = myopti[CD_CYKRNEVI];

  KERNEL_reg_unordered_extern = myopti[CD_UXKRNEVI];
  KERNEL_den_unordered_extern = myopti[CD_UYKRNEVI];

  KERNEL_reg_ordered_extern = myopti[CD_OXKRNEVI];
  KERNEL_den_ordered_extern = myopti[CD_OYKRNEVI];

  vector_cxkerlb_extern = cxkerlb;
  vector_cxkerub_extern = cxkerub;
  int_cxker_bound_extern = np_has_finite_cker_bounds(cxkerlb, cxkerub, num_reg_continuous_extern);

  vector_cykerlb_extern = cykerlb;
  vector_cykerub_extern = cykerub;
  int_cyker_bound_extern = np_has_finite_cker_bounds(cykerlb, cykerub, num_var_continuous_extern);

  if((num_reg_continuous_extern + num_var_continuous_extern) > 0){
    cxylb = alloc_vecd(num_reg_continuous_extern + num_var_continuous_extern);
    cxyub = alloc_vecd(num_reg_continuous_extern + num_var_continuous_extern);
    for(i = 0; i < num_reg_continuous_extern; i++){
      cxylb[i] = (cxkerlb != NULL) ? cxkerlb[i] : DBL_MAX;
      cxyub[i] = (cxkerub != NULL) ? cxkerub[i] : DBL_MAX;
    }
    for(i = 0; i < num_var_continuous_extern; i++){
      cxylb[num_reg_continuous_extern + i] = (cykerlb != NULL) ? cykerlb[i] : DBL_MAX;
      cxyub[num_reg_continuous_extern + i] = (cykerub != NULL) ? cykerub[i] : DBL_MAX;
    }
    vector_cxykerlb_extern = cxylb;
    vector_cxykerub_extern = cxyub;
    int_cxyker_bound_extern = np_has_finite_cker_bounds(cxylb, cxyub, num_reg_continuous_extern + num_var_continuous_extern);
  } else {
    vector_cxykerlb_extern = NULL;
    vector_cxykerub_extern = NULL;
    int_cxyker_bound_extern = 0;
  }

  int_LARGE_SF = myopti[CD_LSFI];
  BANDWIDTH_den_extern = myopti[CD_DENI];
  int_MINIMIZE_IO = myopti[CD_MINIOI];
  do_grad = myopti[CD_GRAD];

  ymax_lev = myopti[CD_YMLEVI];
  xmax_lev = myopti[CD_XMLEVI];

  ypad_num = *ypadnum;
  xpad_num = *xpadnum;

  nconfac_extern = *nconfac;
  ncatfac_extern = *ncatfac;

  dens_or_dist = myopti[CD_DODENI];

  int_TREE_XY = myopti[CD_TREEI]; // we just build a single xy tree
  
  int_TREE_X = int_TREE_Y = NP_TREE_FALSE;

  operator = (dens_or_dist == NP_DO_DENS) ? OP_NORMAL : OP_INTEGRAL;

#ifdef MPI2
  num_obs_eval_alloc = MAX(ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif

  // our method of evaluation uses a single joint xy data matrix

  matrix_XY_continuous_train_extern = alloc_matd(num_obs_train_extern, num_all_cvar);
  matrix_XY_unordered_train_extern = alloc_matd(num_obs_train_extern, num_all_uvar);
  matrix_XY_ordered_train_extern = alloc_matd(num_obs_train_extern, num_all_ovar);

  if(train_is_eval) {
    matrix_XY_continuous_eval_extern = matrix_XY_continuous_train_extern;
    matrix_XY_unordered_eval_extern = matrix_XY_unordered_train_extern;
    matrix_XY_ordered_eval_extern = matrix_XY_ordered_train_extern;
  } else {
    matrix_XY_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_all_cvar);
    matrix_XY_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_all_uvar);
    matrix_XY_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_all_ovar);
  }

  num_categories_extern = alloc_vecu(num_all_catvar);
	
  num_categories_extern_XY = alloc_vecu(num_all_catvar);

  vector_scale_factor = alloc_vecd(num_all_var + 1);
  
  matrix_categorical_vals_extern = alloc_matd(MAX(xmax_lev,ymax_lev), num_all_catvar);
  matrix_categorical_vals_extern_XY = alloc_matd(MAX(xmax_lev,ymax_lev), num_all_catvar);


  /* notice use of num_obs_eval_alloc for MPI compatibility */
  pdf = alloc_vecd(num_obs_eval_alloc);
  pdf_stderr = alloc_vecd(num_obs_eval_alloc);

  if (do_grad){
    pdf_deriv = alloc_matd(num_obs_eval_alloc, num_var);
    pdf_deriv_stderr = alloc_matd(num_obs_eval_alloc, num_var);
  }

  /* in v_s_f order is creg, cvar, uvar, ovar, ureg, oreg  */

  for( i=0;i<num_all_var; i++ )
    vector_scale_factor[i+1] = mybw[i];

  vector_continuous_stddev_extern = mysd;

  /* Parse data */

  // xy arrays contain first x, then y data

  for(j=0;j<num_reg_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_XY_unordered_train_extern[j][i]=tu_uno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_XY_ordered_train_extern[j][i]=tu_ord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_XY_continuous_train_extern[j][i]=tu_con[j*num_obs_train_extern+i];

  // y data
  for(j=0;j<num_var_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_XY_unordered_train_extern[j+num_reg_unordered_extern][i]=tc_uno[j*num_obs_train_extern+i];

  for(j=0;j<num_var_ordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_XY_ordered_train_extern[j+num_reg_ordered_extern][i]=tc_ord[j*num_obs_train_extern+i];

  for(j=0;j<num_var_continuous_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_XY_continuous_train_extern[j+num_reg_continuous_extern][i]=tc_con[j*num_obs_train_extern+i];

  /* eval */
  if(!train_is_eval){
    for(j=0;j<num_reg_unordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_XY_unordered_eval_extern[j][i]=eu_uno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_XY_ordered_eval_extern[j][i]=eu_ord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_XY_continuous_eval_extern[j][i]=eu_con[j*num_obs_eval_extern+i];

    // y data
    for(j=0;j<num_var_unordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_XY_unordered_eval_extern[j+num_reg_unordered_extern][i]=ec_uno[j*num_obs_eval_extern+i];

    for(j=0;j<num_var_ordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_XY_ordered_eval_extern[j+num_reg_ordered_extern][i]=ec_ord[j*num_obs_eval_extern+i];

    for(j=0;j<num_var_continuous_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_XY_continuous_eval_extern[j+num_reg_continuous_extern][i]=ec_con[j*num_obs_eval_extern+i];
  }

  /* fix up categories */
  for(j=0; j < (num_var_unordered_extern + num_var_ordered_extern); j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = ymcv[j*ymax_lev+i];
    } while(++i < ymax_lev && ymcv[j*ymax_lev+i] != ypad_num);
    num_categories_extern[j] = i;
  }

  t_num = j;

  for(j=0; j < (num_reg_unordered_extern+num_reg_ordered_extern); j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j+t_num][i] = xmcv[j*xmax_lev+i];
    } while(++i < xmax_lev && xmcv[j*xmax_lev+i] != xpad_num);
    num_categories_extern[j+t_num] = i;
  }

  // properly fill-in xy categorical data
  np_splitxy_vsf_mcv_nc(num_var_unordered_extern, num_var_ordered_extern, num_var_continuous_extern,
                        num_reg_unordered_extern, num_reg_ordered_extern, num_reg_continuous_extern,
                        vector_scale_factor+1,
                        num_categories_extern,
                        matrix_categorical_vals_extern,
                        NULL, NULL, NULL,
                        NULL, NULL, num_categories_extern_XY,
                        NULL, NULL, matrix_categorical_vals_extern_XY);

  
  // set up indexing data for tree
  ipt_XY = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_XY != NULL))
    error("!(ipt_XY != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt_XY[i] = i;
  }

  if(!train_is_eval) {
    ipe_XY = (int *)malloc(num_obs_eval_extern*sizeof(int));
    if(!(ipe_XY != NULL))
      error("!(ipe_XY != NULL)");

    for(i = 0; i < num_obs_eval_extern; i++){
      ipe_XY[i] = i;
    }
  } else {
    ipe_XY = ipt_XY;
  }

  int_TREE_XY = int_TREE_XY && ((num_all_cvar != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_XY == NP_TREE_TRUE){
    if((BANDWIDTH_den_extern != BW_ADAP_NN) || ((BANDWIDTH_den_extern == BW_ADAP_NN) && train_is_eval)){
      build_kdtree(matrix_XY_continuous_train_extern, num_obs_train_extern, num_all_cvar, 
                   4*num_all_cvar, ipt_XY, &kdt_extern_XY);

      // x
      for(j = 0; j < num_reg_unordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_unordered_train_extern[j][i]=tu_uno[j*num_obs_train_extern+ipt_XY[i]];
        
      for(j = 0; j < num_reg_ordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_ordered_train_extern[j][i]=tu_ord[j*num_obs_train_extern+ipt_XY[i]];

      for(j = 0; j < num_reg_continuous_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_continuous_train_extern[j][i]=tu_con[j*num_obs_train_extern+ipt_XY[i]];
    
      // y
      for(j = 0; j < num_var_unordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_unordered_train_extern[j+num_reg_unordered_extern][i]=tc_uno[j*num_obs_train_extern+ipt_XY[i]];
        
      for(j = 0; j < num_var_ordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_ordered_train_extern[j+num_reg_ordered_extern][i]=tc_ord[j*num_obs_train_extern+ipt_XY[i]];

      for(j = 0; j < num_var_continuous_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_continuous_train_extern[j+num_reg_continuous_extern][i]=tc_con[j*num_obs_train_extern+ipt_XY[i]];
    } else {
      build_kdtree(matrix_XY_continuous_eval_extern, num_obs_eval_extern, num_all_cvar, 
                   4*num_all_cvar, ipe_XY, &kdt_extern_XY);

      // x
      for(j = 0; j < num_reg_unordered_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_unordered_eval_extern[j][i]=eu_uno[j*num_obs_eval_extern+ipe_XY[i]];
        
      for(j = 0; j < num_reg_ordered_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_ordered_eval_extern[j][i]=eu_ord[j*num_obs_eval_extern+ipe_XY[i]];

      for(j = 0; j < num_reg_continuous_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_continuous_eval_extern[j][i]=eu_con[j*num_obs_eval_extern+ipe_XY[i]];
    
      // y
      for(j = 0; j < num_var_unordered_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_unordered_eval_extern[j+num_reg_unordered_extern][i]=ec_uno[j*num_obs_eval_extern+ipe_XY[i]];
        
      for(j = 0; j < num_var_ordered_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_ordered_eval_extern[j+num_reg_ordered_extern][i]=ec_ord[j*num_obs_eval_extern+ipe_XY[i]];

      for(j = 0; j < num_var_continuous_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_continuous_eval_extern[j+num_reg_continuous_extern][i]=ec_con[j*num_obs_eval_extern+ipe_XY[i]];

    }
  }

  int_ll_eff = int_ll_extern;
  if((int_ll_eff == LL_LP) && (num_reg_continuous_extern == 0))
    int_ll_eff = LL_LC;

  if((int_ll_eff == LL_LC) ||
     ((int_ll_eff == LL_LP) && (BANDWIDTH_den_extern == BW_ADAP_NN))){
    np_kernel_estimate_con_dens_dist_categorical(KERNEL_den_extern,
                                                 KERNEL_den_unordered_extern,
                                                 KERNEL_den_ordered_extern,
                                                 KERNEL_reg_extern,
                                                 KERNEL_reg_unordered_extern,
                                                 KERNEL_reg_ordered_extern,
                                                 BANDWIDTH_den_extern,
                                                 operator,
                                                 num_obs_train_extern,
                                                 num_obs_eval_extern,
                                                 num_var_unordered_extern,
                                                 num_var_ordered_extern,
                                                 num_var_continuous_extern,
                                                 num_reg_unordered_extern,
                                                 num_reg_ordered_extern,
                                                 num_reg_continuous_extern,
                                                 matrix_XY_unordered_train_extern,
                                                 matrix_XY_ordered_train_extern,
                                                 matrix_XY_continuous_train_extern,
                                                 matrix_XY_unordered_eval_extern,
                                                 matrix_XY_ordered_eval_extern,
                                                 matrix_XY_continuous_eval_extern,
                                                 &vector_scale_factor[1],
                                                 num_categories_extern,
                                                 num_categories_extern_XY,
                                                 matrix_categorical_vals_extern,
                                                 matrix_categorical_vals_extern_XY,
                                                 pdf,
                                                 pdf_stderr,
                                                 pdf_deriv,
                                                 pdf_deriv_stderr,
                                                 &log_likelihood);
  } else {
    int status = 0;
    int lp_eval_alloc = 1;
    double *vsf_x = NULL, *vsf_y = NULL, *ykw = NULL, *y_eval_one = NULL;
    double *mean_one = NULL, *stderr_one = NULL;
    double **xuno_eval_one = NULL, **xord_eval_one = NULL, **xcon_eval_one = NULL;
    double **yuno_eval_one = NULL, **yord_eval_one = NULL, **ycon_eval_one = NULL;
    double **grad_one = NULL, **graderr_one = NULL;
    int *kernel_cy = NULL, *kernel_uy = NULL, *kernel_oy = NULL, *operator_y = NULL;
    double RS = 0.0, MSE = 0.0, MAE = 0.0, MAPE = 0.0, CORR = 0.0, SIGN = 0.0;
    int num_y_vars = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
    int num_x_vars = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;

#ifdef MPI2
    lp_eval_alloc = MAX((int)ceil(1.0 / (double)iNum_Processors), 1) * iNum_Processors;
#endif

    if((int_ll_eff == LL_LP) &&
       ((vector_glp_degree_extern == NULL) || (num_reg_continuous_extern <= 0)))
      error("np_density_conditional: LP conditional path requires continuous x variables and GLP degree metadata");

    ycat_offset = 0;
    xcat_offset = num_var_unordered_extern + num_var_ordered_extern;

    vsf_x = alloc_vecd(MAX(1, num_x_vars));
    vsf_y = alloc_vecd(MAX(1, num_y_vars));
    ykw = alloc_vecd(MAX(1, num_obs_train_extern));
    y_eval_one = alloc_vecd(1);
    mean_one = alloc_vecd(MAX(1, lp_eval_alloc));
    stderr_one = alloc_vecd(MAX(1, lp_eval_alloc));

    if(num_reg_unordered_extern > 0) xuno_eval_one = alloc_matd(1, num_reg_unordered_extern);
    if(num_reg_ordered_extern > 0) xord_eval_one = alloc_matd(1, num_reg_ordered_extern);
    if(num_reg_continuous_extern > 0) xcon_eval_one = alloc_matd(1, num_reg_continuous_extern);

    if(num_var_unordered_extern > 0) yuno_eval_one = alloc_matd(1, num_var_unordered_extern);
    if(num_var_ordered_extern > 0) yord_eval_one = alloc_matd(1, num_var_ordered_extern);
    if(num_var_continuous_extern > 0) ycon_eval_one = alloc_matd(1, num_var_continuous_extern);

    if(do_grad && (num_x_vars > 0)){
      grad_one = alloc_matd(MAX(1, lp_eval_alloc), num_x_vars);
      graderr_one = alloc_matd(MAX(1, lp_eval_alloc), num_x_vars);
    }

    kernel_cy = (int *)calloc((size_t)MAX(1, num_var_continuous_extern), sizeof(int));
    kernel_uy = (int *)calloc((size_t)MAX(1, num_var_unordered_extern), sizeof(int));
    kernel_oy = (int *)calloc((size_t)MAX(1, num_var_ordered_extern), sizeof(int));
    operator_y = (int *)calloc((size_t)MAX(1, num_y_vars), sizeof(int));

    if((vsf_x == NULL) || (vsf_y == NULL) || (ykw == NULL) ||
       (y_eval_one == NULL) || (mean_one == NULL) || (stderr_one == NULL) ||
       (kernel_cy == NULL) || (kernel_uy == NULL) || (kernel_oy == NULL) ||
       (operator_y == NULL) ||
       ((num_reg_unordered_extern > 0) && (xuno_eval_one == NULL)) ||
       ((num_reg_ordered_extern > 0) && (xord_eval_one == NULL)) ||
       ((num_reg_continuous_extern > 0) && (xcon_eval_one == NULL)) ||
       ((num_var_unordered_extern > 0) && (yuno_eval_one == NULL)) ||
       ((num_var_ordered_extern > 0) && (yord_eval_one == NULL)) ||
       ((num_var_continuous_extern > 0) && (ycon_eval_one == NULL)) ||
       (do_grad && (num_x_vars > 0) && ((grad_one == NULL) || (graderr_one == NULL))))
      error("np_density_conditional: memory allocation failed in conditional LP path");

    for(i = 0; i < num_var_continuous_extern; i++) kernel_cy[i] = KERNEL_den_extern;
    for(i = 0; i < num_var_unordered_extern; i++) kernel_uy[i] = KERNEL_den_unordered_extern;
    for(i = 0; i < num_var_ordered_extern; i++) kernel_oy[i] = KERNEL_den_ordered_extern;
    for(i = 0; i < num_y_vars; i++) operator_y[i] = operator;

    for(i = 0; i < num_reg_continuous_extern; i++)
      vsf_x[i] = vector_scale_factor[1 + i];
    for(i = 0; i < num_reg_unordered_extern; i++)
      vsf_x[num_reg_continuous_extern + i] =
        vector_scale_factor[1 + num_reg_continuous_extern + num_var_continuous_extern +
                            num_var_unordered_extern + num_var_ordered_extern + i];
    for(i = 0; i < num_reg_ordered_extern; i++)
      vsf_x[num_reg_continuous_extern + num_reg_unordered_extern + i] =
        vector_scale_factor[1 + num_reg_continuous_extern + num_var_continuous_extern +
                            num_var_unordered_extern + num_var_ordered_extern +
                            num_reg_unordered_extern + i];

    for(i = 0; i < num_var_continuous_extern; i++)
      vsf_y[i] = vector_scale_factor[1 + num_reg_continuous_extern + i];
    for(i = 0; i < num_var_unordered_extern; i++)
      vsf_y[num_var_continuous_extern + i] =
        vector_scale_factor[1 + num_reg_continuous_extern + num_var_continuous_extern + i];
    for(i = 0; i < num_var_ordered_extern; i++)
      vsf_y[num_var_continuous_extern + num_var_unordered_extern + i] =
        vector_scale_factor[1 + num_reg_continuous_extern + num_var_continuous_extern +
                            num_var_unordered_extern + i];

    saved_cker_bound = int_cker_bound_extern;
    saved_ckerlb = vector_ckerlb_extern;
    saved_ckerub = vector_ckerub_extern;
    log_likelihood = 0.0;

    for(j = 0; j < num_obs_eval_extern; j++){
      for(i = 0; i < num_reg_unordered_extern; i++)
        xuno_eval_one[i][0] = matrix_XY_unordered_eval_extern[i][j];
      for(i = 0; i < num_reg_ordered_extern; i++)
        xord_eval_one[i][0] = matrix_XY_ordered_eval_extern[i][j];
      for(i = 0; i < num_reg_continuous_extern; i++)
        xcon_eval_one[i][0] = matrix_XY_continuous_eval_extern[i][j];

      for(i = 0; i < num_var_unordered_extern; i++)
        yuno_eval_one[i][0] = matrix_XY_unordered_eval_extern[num_reg_unordered_extern + i][j];
      for(i = 0; i < num_var_ordered_extern; i++)
        yord_eval_one[i][0] = matrix_XY_ordered_eval_extern[num_reg_ordered_extern + i][j];
      for(i = 0; i < num_var_continuous_extern; i++)
        ycon_eval_one[i][0] = matrix_XY_continuous_eval_extern[num_reg_continuous_extern + i][j];

      int_cker_bound_extern = int_cyker_bound_extern;
      vector_ckerlb_extern = vector_cykerlb_extern;
      vector_ckerub_extern = vector_cykerub_extern;

      status = kernel_weighted_sum_np(kernel_cy,
                                      kernel_uy,
                                      kernel_oy,
                                      BANDWIDTH_den_extern,
                                      num_obs_train_extern,
                                      1,
                                      num_var_unordered_extern,
                                      num_var_ordered_extern,
                                      num_var_continuous_extern,
                                      0,
                                      0,
                                      1,
                                      1,
                                      1,
                                      0,
                                      0,
                                      0,
                                      0,
                                      operator_y,
                                      OP_NOOP,
                                      0,
                                      0,
                                      NULL,
                                      1,
                                      0,
                                      0,
                                      int_TREE_Y,
                                      0,
                                      NULL,
                                      NULL,
                                      NULL,
                                      NULL,
                                      (num_var_unordered_extern > 0) ? matrix_XY_unordered_train_extern + num_reg_unordered_extern : NULL,
                                      (num_var_ordered_extern > 0) ? matrix_XY_ordered_train_extern + num_reg_ordered_extern : NULL,
                                      (num_var_continuous_extern > 0) ? matrix_XY_continuous_train_extern + num_reg_continuous_extern : NULL,
                                      yuno_eval_one,
                                      yord_eval_one,
                                      ycon_eval_one,
                                      NULL,
                                      NULL,
                                      NULL,
                                      vsf_y,
                                      0,
                                      NULL,
                                      NULL,
                                      NULL,
                                      num_categories_extern + ycat_offset,
                                      matrix_categorical_vals_extern + ycat_offset,
                                      NULL,
                                      NULL,
                                      NULL,
                                      ykw,
                                      NULL);
      if(status != 0)
        error("np_density_conditional: y-kernel response construction failed in LP path");

      y_eval_one[0] = ykw[0];

      int_cker_bound_extern = int_cxker_bound_extern;
      vector_ckerlb_extern = vector_cxkerlb_extern;
      vector_ckerub_extern = vector_cxkerub_extern;

      status = kernel_estimate_regression_categorical_tree_np(int_ll_eff,
                                                               KERNEL_reg_extern,
                                                               KERNEL_reg_unordered_extern,
                                                               KERNEL_reg_ordered_extern,
                                                               BANDWIDTH_den_extern,
                                                               num_obs_train_extern,
                                                               1,
                                                               num_reg_unordered_extern,
                                                               num_reg_ordered_extern,
                                                               num_reg_continuous_extern,
                                                               matrix_XY_unordered_train_extern,
                                                               matrix_XY_ordered_train_extern,
                                                               matrix_XY_continuous_train_extern,
                                                               xuno_eval_one,
                                                               xord_eval_one,
                                                               xcon_eval_one,
                                                               ykw,
                                                               y_eval_one,
                                                               vsf_x,
                                                               num_categories_extern + xcat_offset,
                                                               matrix_categorical_vals_extern + xcat_offset,
                                                               mean_one,
                                                               grad_one,
                                                               stderr_one,
                                                               graderr_one,
                                                               &RS,
                                                               &MSE,
                                                               &MAE,
                                                               &MAPE,
                                                               &CORR,
                                                               &SIGN);

      if(status != 0)
        error("np_density_conditional: regression LP solve failed in conditional LP path");

      pdf[j] = mean_one[0];
      pdf_stderr[j] = stderr_one[0];

      if(dens_or_dist == NP_DO_DENS){
        const double val = (pdf[j] < DBL_MIN) ? DBL_MIN : pdf[j];
        log_likelihood += log(val);
      }

      if(do_grad){
        for(i = 0; i < num_x_vars; i++){
          pdf_deriv[i][j] = (grad_one != NULL) ? grad_one[i][0] : 0.0;
          pdf_deriv_stderr[i][j] = (graderr_one != NULL) ? graderr_one[i][0] : 0.0;
        }
      }

      np_progress_fit_step(j + 1);
    }

    int_cker_bound_extern = saved_cker_bound;
    vector_ckerlb_extern = saved_ckerlb;
    vector_ckerub_extern = saved_ckerub;

    safe_free(vsf_x);
    safe_free(vsf_y);
    safe_free(ykw);
    safe_free(y_eval_one);
    safe_free(mean_one);
    safe_free(stderr_one);

    if(xuno_eval_one != NULL) free_mat(xuno_eval_one, num_reg_unordered_extern);
    if(xord_eval_one != NULL) free_mat(xord_eval_one, num_reg_ordered_extern);
    if(xcon_eval_one != NULL) free_mat(xcon_eval_one, num_reg_continuous_extern);

    if(yuno_eval_one != NULL) free_mat(yuno_eval_one, num_var_unordered_extern);
    if(yord_eval_one != NULL) free_mat(yord_eval_one, num_var_ordered_extern);
    if(ycon_eval_one != NULL) free_mat(ycon_eval_one, num_var_continuous_extern);

    if(grad_one != NULL) free_mat(grad_one, num_x_vars);
    if(graderr_one != NULL) free_mat(graderr_one, num_x_vars);

    safe_free(kernel_cy);
    safe_free(kernel_uy);
    safe_free(kernel_oy);
    safe_free(operator_y);
  }


  /* return data to R */
  for( i=0; i<num_obs_eval_extern; i++ )
    cdens[ipe_XY[i]]=pdf[i];

  for( i=0; i<num_obs_eval_extern; i++ )
    cderr[ipe_XY[i]]=pdf_stderr[i];
  
  if (do_grad) {
    for(j=0;j<num_var;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        cgerr[j*num_obs_eval_extern+ipe_XY[i]]=pdf_deriv_stderr[j][i];

    for(j=0;j<num_var;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        cg[j*num_obs_eval_extern+ipe_XY[i]]=pdf_deriv[j][i];
  }



  *ll = log_likelihood;
  /* end return data */

  /* Free data objects */

  free_mat(matrix_XY_unordered_train_extern, num_all_uvar);
  free_mat(matrix_XY_ordered_train_extern, num_all_ovar);
  free_mat(matrix_XY_continuous_train_extern, num_all_cvar);

  if(train_is_eval){
    matrix_XY_unordered_eval_extern = NULL;
    matrix_XY_ordered_eval_extern = NULL;
    matrix_XY_continuous_eval_extern = NULL;
  } else {
    free_mat(matrix_XY_unordered_eval_extern, num_all_uvar);
    free_mat(matrix_XY_ordered_eval_extern, num_all_ovar);
    free_mat(matrix_XY_continuous_eval_extern, num_all_cvar);
  }

  if (do_grad){
    free_mat(pdf_deriv, num_var);
    free_mat(pdf_deriv_stderr, num_var);
  }

  vector_continuous_stddev_extern = NULL;

  safe_free(vector_scale_factor);
  safe_free(num_categories_extern);
  safe_free(num_categories_extern_XY);
  safe_free(pdf);
  safe_free(pdf_stderr);

  free_mat(matrix_categorical_vals_extern, num_all_catvar);
  free_mat(matrix_categorical_vals_extern_XY, num_all_catvar);

  safe_free(ipt_XY);

  if(!train_is_eval)
    safe_free(ipe_XY);

  if(int_TREE_XY == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_XY);
    int_TREE_XY = NP_TREE_FALSE;
  }

  int_cxker_bound_extern = 0;
  int_cyker_bound_extern = 0;
  int_cxyker_bound_extern = 0;
  vector_cxkerlb_extern = NULL;
  vector_cxkerub_extern = NULL;
  vector_cykerlb_extern = NULL;
  vector_cykerub_extern = NULL;
  vector_cxykerlb_extern = NULL;
  vector_cxykerub_extern = NULL;
  int_cker_bound_extern = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;
  safe_free(cxylb);
  safe_free(cxyub);

  return;
}


void np_density(double * tuno, double * tord, double * tcon, 
                double * euno, double * eord, double * econ, 
                double * dbw, 
                double * mcv, double * padnum, 
                double * nconfac, double * ncatfac, double * mysd,
                int * myopti, double * mydens, double * myderr, double * ll,
                double * ckerlb, double * ckerub){


  double small = 1.0e-16;
  double * vector_scale_factor, * pdf, * pdf_stderr, log_likelihood = 0.0;
  double pad_num;

  int itmax = 10000;
  int i,j;
  int num_var, num_obs_eval_alloc, max_lev, train_is_eval, dens_or_dist, old_dens;

  int * ipt = NULL, * ipe = NULL;
  

  /* match integer options with their globals */

  num_reg_continuous_extern = myopti[DEN_NCONI];
  num_reg_unordered_extern = myopti[DEN_NUNOI];
  num_reg_ordered_extern = myopti[DEN_NORDI];

  vector_ckerlb_extern = ckerlb;
  vector_ckerub_extern = ckerub;
  int_cker_bound_extern = np_has_finite_cker_bounds(ckerlb, ckerub, num_reg_continuous_extern);

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[DEN_TNOBSI];
  num_obs_eval_extern = myopti[DEN_ENOBSI];

  KERNEL_den_extern = myopti[DEN_CKRNEVI];
  KERNEL_den_unordered_extern = myopti[DEN_UKRNEVI];
  KERNEL_den_ordered_extern = myopti[DEN_OKRNEVI];

  int_LARGE_SF = myopti[DEN_LSFI];
  int_MINIMIZE_IO = myopti[DEN_MINIOI];
  BANDWIDTH_den_extern = myopti[DEN_DENI];

  train_is_eval = myopti[DEN_TISEI];

  max_lev = myopti[DEN_MLEVI];
  pad_num = *padnum;

  nconfac_extern = *nconfac;
  ncatfac_extern = *ncatfac;

  dens_or_dist = myopti[DEN_DODENI];
  old_dens = myopti[DEN_OLDI];
  int_TREE_X = myopti[DEN_TREEI];

#ifdef MPI2
  num_obs_eval_alloc = MAX(ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif

  /* Allocate memory for objects */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  if(!train_is_eval){
    matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
    matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
    matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);
  } else {
    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;
  }

  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  vector_scale_factor = alloc_vecd(num_var + 1);
  matrix_categorical_vals_extern = alloc_matd(max_lev, num_reg_unordered_extern + num_reg_ordered_extern);

  /* note use of num_obs_eval_alloc */
  pdf = alloc_vecd(num_obs_eval_alloc);
  pdf_stderr = alloc_vecd(num_obs_eval_alloc);
  
  vector_continuous_stddev_extern = mysd;

  /* Parse data */
	
  /* train */

  for( j=0;j<num_reg_unordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+i];

  /* eval */
  if (!train_is_eval) {
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+i];
  }

  /*  bandwidths/scale factors */

  for( i=0; i<num_var; i++ )
    vector_scale_factor[i+1]=dbw[i];

  /* fix up categories */
  
  for(j=0; j < (num_reg_unordered_extern + num_reg_ordered_extern); j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;
  }

  /* data has been copied, now build tree */

  ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  if(!train_is_eval) {
    ipe = (int *)malloc(num_obs_eval_extern*sizeof(int));
    if(!(ipe != NULL))
      error("!(ipe != NULL)");

    for(i = 0; i < num_obs_eval_extern; i++){
      ipe[i] = i;
    }
  } else {
    ipe = ipt;
  }

  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    if((BANDWIDTH_den_extern != BW_ADAP_NN) || ((BANDWIDTH_den_extern == BW_ADAP_NN) && train_is_eval)){
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipt, &kdt_extern_X);

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+ipt[i]];
    
    
      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+ipt[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+ipt[i]];

    } else {
      build_kdtree(matrix_X_continuous_eval_extern, num_obs_eval_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipe, &kdt_extern_X);


      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+ipe[i]];
    }

  }


  /* Conduct estimation */
  
  if(old_dens){
    if (dens_or_dist == NP_DO_DENS){
      /* nb - KERNEL_(|un)ordered_den are set to zero upon declaration 
         - they have only one kernel type each at the moment */
      kernel_estimate_density_categorical(KERNEL_den_extern,
                                          KERNEL_den_unordered_extern,
                                          KERNEL_den_ordered_extern,
                                          BANDWIDTH_den_extern,
                                          num_obs_train_extern,
                                          num_obs_eval_extern,
                                          num_reg_unordered_extern,
                                          num_reg_ordered_extern,
                                          num_reg_continuous_extern,
                                          /* Train */
                                          matrix_X_unordered_train_extern,
                                          matrix_X_ordered_train_extern,
                                          matrix_X_continuous_train_extern,
                                          /* Eval */
                                          matrix_X_unordered_eval_extern,
                                          matrix_X_ordered_eval_extern,
                                          matrix_X_continuous_eval_extern,
                                          &vector_scale_factor[1],
                                          num_categories_extern,
                                          pdf,
                                          pdf_stderr,
                                          &log_likelihood);
    } else if (dens_or_dist == NP_DO_DIST) {
      kernel_estimate_distribution_categorical(KERNEL_den_extern,
                                               KERNEL_den_unordered_extern,
                                               KERNEL_den_ordered_extern,
                                               BANDWIDTH_den_extern,
                                               num_obs_train_extern,
                                               num_obs_eval_extern,
                                               num_reg_unordered_extern,
                                               num_reg_ordered_extern,
                                               num_reg_continuous_extern,
                                               /* Train */
                                               matrix_X_unordered_train_extern,
                                               matrix_X_ordered_train_extern,
                                               matrix_X_continuous_train_extern,
                                               /* Eval */
                                               matrix_X_unordered_eval_extern,
                                               matrix_X_ordered_eval_extern,
                                               matrix_X_continuous_eval_extern,
                                               &vector_scale_factor[1],
                                               num_categories_extern,
                                               matrix_categorical_vals_extern,
                                               pdf,
                                               pdf_stderr,
                                               small, itmax);

    }
  } else {
    const int dop = (dens_or_dist == NP_DO_DENS) ? OP_NORMAL : OP_INTEGRAL;

      np_progress_fit_set_offset(0);
      kernel_estimate_dens_dist_categorical_np(KERNEL_den_extern,
                                               KERNEL_den_unordered_extern,
                                               KERNEL_den_ordered_extern,
                                               BANDWIDTH_den_extern,
                                               num_obs_train_extern,
                                               num_obs_eval_extern,
                                               num_reg_unordered_extern,
                                               num_reg_ordered_extern,
                                               num_reg_continuous_extern,
                                               dop,
                                               /* Train */
                                               matrix_X_unordered_train_extern,
                                               matrix_X_ordered_train_extern,
                                               matrix_X_continuous_train_extern,
                                               /* Eval */
                                               matrix_X_unordered_eval_extern,
                                               matrix_X_ordered_eval_extern,
                                               matrix_X_continuous_eval_extern,
                                               &vector_scale_factor[1],
                                               num_categories_extern,
                                               matrix_categorical_vals_extern,
                                               pdf,
                                               pdf_stderr,
                                               &log_likelihood);
  }
  
  
  /* write the return values */

  for(i=0;i<num_obs_eval_extern;i++){
    mydens[ipe[i]] = pdf[i];
    myderr[ipe[i]] = pdf_stderr[i];
  }
  *ll = log_likelihood;

  /* clean up and wave goodbye */

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  if (!train_is_eval){
    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  vector_continuous_stddev_extern = NULL;

  safe_free(vector_scale_factor);
  safe_free(num_categories_extern);
  safe_free(pdf_stderr);
  safe_free(pdf);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  safe_free(ipt);

  if(!train_is_eval)
    safe_free(ipe);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  int_cker_bound_extern = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;

  return;
}


static void np_regression_bw_mode(double * runo, double * rord, double * rcon, double * y,
                                  double * mysd, int * myopti, double * myoptd, double * rbw, double * fval,
                                  double * objective_function_values, double * objective_function_evals,
                                  double * objective_function_invalid, double * timing,
                                  double * objective_function_fast,
                                  int * penalty_mode, double * penalty_mult,
                                  int * glp_degree,
                                  int * glp_bernstein,
                                  int * glp_basis,
                                  double * ckerlb, double * ckerub,
                                  const int eval_only){
  //KDT * kdt = NULL; // tree structure
  //NL nl = { .node = NULL, .n = 0, .nalloc = 0 };// a node list structure -- used for searching - here for testing
  //double tb[4] = {0.25, 0.5, 0.3, 0.75};
  int * ipt = NULL;  // point permutation, see tree.c

  double **matrix_y;

  double *vector_continuous_stddev;
  double *vector_scale_factor, *vector_scale_factor_multistart, *vsfh;
  double *vector_scale_factor_startbest;

  double fret, fret_best, fret_start_best, fret_initial;
  double ftol, tol, small;
  double (* bwmfunc)(double *) = NULL;

  double lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init; 
  double lbd_init, hbd_init, d_init;
  int dfc_dir;

  int i,j;
  double fast_eval_total = 0.0;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, iImproved;
  int enforce_fixed_feasibility;
  int have_start_best, have_multistart_best;
  int itmax, iter;
  int int_use_starting_values;

  int scale_cat;
  const char *bw_error_msg = NULL;

  num_reg_continuous_extern = myopti[RBW_NCONI];
  num_reg_unordered_extern = myopti[RBW_NUNOI];
  num_reg_ordered_extern = myopti[RBW_NORDI];
  np_reset_y_side_extern();
  bwm_clear_floor_context();

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[RBW_NOBSI];
  iMultistart = myopti[RBW_IMULTII];
  iNum_Multistart = myopti[RBW_NMULTII];
  if (!eval_only && iNum_Multistart < 1)
    error("C_np_regression_bw: nmulti must be a positive integer");

  KERNEL_reg_extern = myopti[RBW_CKRNEVI];
  KERNEL_reg_unordered_extern = myopti[RBW_UKRNEVI];
  KERNEL_reg_ordered_extern = myopti[RBW_OKRNEVI];

  vector_ckerlb_extern = ckerlb;
  vector_ckerub_extern = ckerub;
  int_cker_bound_extern = np_has_finite_cker_bounds(ckerlb, ckerub, num_reg_continuous_extern);

  int_use_starting_values= myopti[RBW_USTARTI];
  int_LARGE_SF=myopti[RBW_LSFI];

  BANDWIDTH_reg_extern=myopti[RBW_REGI];
  BANDWIDTH_den_extern=0;
  enforce_fixed_feasibility = ((BANDWIDTH_reg_extern == BW_FIXED) && (!eval_only));

  itmax=myopti[RBW_ITMAXI];
  int_RESTART_FROM_MIN = myopti[RBW_REMINI];
  int_MINIMIZE_IO = myopti[RBW_MINIOI];

  int_ll_extern = myopti[RBW_LL];
  vector_glp_degree_extern = glp_degree;
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = *glp_bernstein;
  int_glp_basis_extern = *glp_basis;
  int_nn_k_min_extern =
    ((BANDWIDTH_reg_extern != BW_FIXED) && (num_reg_continuous_extern > 0)) ? 2 : 1;

  int_TREE_X = myopti[RBW_DOTREEI];
  scale_cat = myopti[RBW_SCATI];
  bwm_use_transform = myopti[RBW_TBNDI];
  if (BANDWIDTH_reg_extern != BW_FIXED)
    bwm_use_transform = 0;

  ftol=myoptd[RBW_FTOLD];
  tol=myoptd[RBW_TOLD];
  small=myoptd[RBW_SMALLD];

  dfc_dir = myopti[RBW_DFC_DIRI];
  lbc_dir = myoptd[RBW_LBC_DIRD];
  c_dir = myoptd[RBW_C_DIRD];
  initc_dir = myoptd[RBW_INITC_DIRD]; 

  lbd_dir = myoptd[RBW_LBD_DIRD]; 
  hbd_dir = myoptd[RBW_HBD_DIRD]; 
  d_dir = myoptd[RBW_D_DIRD]; 
  initd_dir = myoptd[RBW_INITD_DIRD]; 

  lbc_init = myoptd[RBW_LBC_INITD]; 
  hbc_init = myoptd[RBW_HBC_INITD]; 
  c_init = myoptd[RBW_C_INITD]; 

  lbd_init = myoptd[RBW_LBD_INITD]; 
  hbd_init = myoptd[RBW_HBD_INITD]; 
  d_init = myoptd[RBW_D_INITD]; 

  nconfac_extern = myoptd[RBW_NCONFD];
  ncatfac_extern = myoptd[RBW_NCATFD];
  bwm_set_scale_factor_lower_bound(myoptd[RBW_SFLOORD]);

  imsnum = 0;
  imstot = iNum_Multistart;

  /* Allocate memory for objects */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  vector_Y_extern = alloc_vecd(num_obs_train_extern);
	
  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_y = alloc_matd(num_var + 1, num_var +1);
  vector_scale_factor = alloc_vecd(num_var + 1);
  vector_scale_factor_startbest = alloc_vecd(num_var + 1);
  vsfh = alloc_vecd(num_var + 1);
  matrix_categorical_vals_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  vector_continuous_stddev = alloc_vecd(num_reg_continuous_extern);

  for(j = 0; j < num_reg_continuous_extern; j++)
    vector_continuous_stddev[j] = mysd[j];

  vector_continuous_stddev_extern = vector_continuous_stddev;

  /* Request starting values for optimization if values already exist */

  /* bandwidths */

  if (int_use_starting_values)
    for( i=0;i<num_var; i++ )
      vector_scale_factor[i+1] = rbw[i];

  /* regressors */

  for( j=0;j<num_reg_unordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_unordered_train_extern[j][i]=runo[j*num_obs_train_extern+i];
    

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=rord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=rcon[j*num_obs_train_extern+i];

  /* response variable */
  for( i=0;i<num_obs_train_extern;i++ )
    vector_Y_extern[i] = y[i];

  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  if (penalty_mode[0] == 1) {
    double pmult = penalty_mult[0];
    double y_mean = 0.0;
    double mse0 = 0.0;
    if (pmult < 1.0) pmult = 1.0;
    for (i = 0; i < num_obs_train_extern; i++)
      y_mean += vector_Y_extern[i];
    y_mean /= (double) num_obs_train_extern;
    for (i = 0; i < num_obs_train_extern; i++) {
      const double dy = vector_Y_extern[i] - y_mean;
      mse0 += dy*dy;
    }
    mse0 /= (double) num_obs_train_extern;
    if (mse0 <= 0.0) mse0 = DBL_MIN;
    if (myopti[RBW_MI] == RBWM_CVAIC) {
      const double denom = 1.0 - 2.0/((double) num_obs_train_extern);
      if (denom > 0.0) {
        const double base_aic = log(mse0) + (1.0/denom);
        bwm_penalty_value = base_aic + log(pmult);
      }
    } else {
      bwm_penalty_value = pmult * mse0;
    }
    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }

  // initialize permutation arrays
  ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  // attempt tree build, if enabled 
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                 4*num_reg_continuous_extern, ipt, &kdt_extern_X);

    //put training data into tree-order using the index array

    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_unordered_train_extern[j][i]=runo[j*num_obs_train_extern+ipt[i]];
    
    
    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_ordered_train_extern[j][i]=rord[j*num_obs_train_extern+ipt[i]];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_continuous_train_extern[j][i]=rcon[j*num_obs_train_extern+ipt[i]];

    /* response variable */
    for( i=0;i<num_obs_train_extern;i++ )
      vector_Y_extern[i] = y[ipt[i]];
    
    //boxSearch(kdt_extern, 0, tb, &nl);
  }

  determine_categorical_vals(
                             num_obs_train_extern,
                             0,
                             0,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             matrix_Y_unordered_train_extern,
                             matrix_Y_ordered_train_extern,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern,
                             matrix_categorical_vals_extern);

  if((int_ll_extern == LL_LP) &&
     (!np_glp_cv_prepare_extern(int_ll_extern,
                                num_obs_train_extern,
                                num_reg_continuous_extern,
                                matrix_X_continuous_train_extern))){
    error("failed to prepare LP CV basis cache");
  }


  /* Initialize scale factors and Directions for NR modules */

  initialize_nr_vector_scale_factor(BANDWIDTH_reg_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0,
                                    KERNEL_reg_unordered_extern,
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_reg_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0,
                                    KERNEL_reg_unordered_extern,
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vsfh,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_reg_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           0,
                           0,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED, 
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);


  /* When multistarting, set counter */

  iMs_counter = 0;

  /* assign the function to be optimized */
  switch(myopti[RBW_MI]){
  case RBWM_CVAIC : bwmfunc = cv_func_regression_categorical_aic_c; break;
  case RBWM_CVLS : bwmfunc = cv_func_regression_categorical_ls; break;
  default : REprintf("np.c: invalid bandwidth selection method.");
    error("np.c: invalid bandwidth selection method.");break;
  }

  spinner(0);

  bwmfunc_raw = bwmfunc;
  bwm_num_reg_continuous = num_reg_continuous_extern;
  bwm_num_reg_unordered = num_reg_unordered_extern;
  bwm_num_reg_ordered = num_reg_ordered_extern;
  bwm_kernel_unordered = KERNEL_reg_unordered_extern;
  bwm_num_categories = num_categories_extern;
  bwm_set_floor_context(
    enforce_fixed_feasibility,
    num_var,
    bwm_use_transform,
    KERNEL_reg_extern,
    KERNEL_reg_unordered_extern,
    BANDWIDTH_reg_extern,
    BANDWIDTH_reg_extern,
    0,
    num_obs_train_extern,
    0,
    0,
    0,
    num_reg_continuous_extern,
    num_reg_unordered_extern,
    num_reg_ordered_extern,
    num_categories_extern,
    bwm_scale_factor_lower_bound);
  if (bwm_use_transform) {
    int n = bwm_num_reg_continuous + bwm_num_reg_unordered + bwm_num_reg_ordered;
    bwm_reserve_transform_buf(n + 1);
  }
  bwm_reset_counters();

  fret_initial = fret_best = bwmfunc_wrapper(vector_scale_factor);
  iImproved = 0;
  have_start_best = 0;
  have_multistart_best = 0;
  fret_start_best = DBL_MAX;
  if (enforce_fixed_feasibility &&
      np_bw_candidate_is_admissible(
        num_var,
        bwm_use_transform,
        KERNEL_reg_extern,
        KERNEL_reg_unordered_extern,
        BANDWIDTH_reg_extern,
        BANDWIDTH_reg_extern,
        0,
        num_obs_train_extern,
        0,
        0,
        0,
        num_reg_continuous_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_categories_extern,
        vector_scale_factor)) {
    have_start_best = 1;
    fret_start_best = fret_initial;
    np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor, num_var);
  }

  if(!eval_only){
    powell(0,
           0,
           vector_scale_factor,
           vector_scale_factor,
           matrix_y,
           num_var,
           ftol,
           tol,
           small,
           itmax,
           &iter,
           &fret,
           bwmfunc_wrapper);


    if(int_RESTART_FROM_MIN == RE_MIN_TRUE){
      initialize_nr_directions(BANDWIDTH_reg_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               0,
                               0,
                               0,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               0, int_RANDOM_SEED, 
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);

      powell(0,
             0,
             vector_scale_factor,
             vector_scale_factor,
             matrix_y,
             num_var,
             ftol,
             tol,
             small,
             itmax,
             &iter,
             &fret,
             bwmfunc_wrapper);

    }
  } else {
    fret = fret_best;
  }

  if (enforce_fixed_feasibility &&
      np_bw_candidate_is_admissible(
        num_var,
        bwm_use_transform,
        KERNEL_reg_extern,
        KERNEL_reg_unordered_extern,
        BANDWIDTH_reg_extern,
        BANDWIDTH_reg_extern,
        0,
        num_obs_train_extern,
        0,
        0,
        0,
        num_reg_continuous_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_categories_extern,
        vector_scale_factor) &&
      ((!have_start_best) || (fret < fret_start_best))) {
    have_start_best = 1;
    fret_start_best = fret;
    np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor, num_var);
  }

  if (enforce_fixed_feasibility) {
    if (have_start_best) {
      fret = fret_start_best;
      np_copy_scale_factor(vector_scale_factor, vector_scale_factor_startbest, num_var);
    } else {
      fret = DBL_MAX;
    }
  }

  iImproved = (enforce_fixed_feasibility && have_start_best) ? (fret_start_best < fret_initial) : (fret < fret_best);
  *timing = timing_extern;

  objective_function_values[0]=fret;
  objective_function_evals[0]=bwm_eval_count;
  objective_function_invalid[0]=bwm_invalid_count;
  bwm_snapshot_fast_counters();
  fast_eval_total += bwm_fast_eval_count;
  /* When multistarting save initial minimum of objective function and scale factors */


  if((!eval_only) && (iMultistart == IMULTI_TRUE)){
    if (enforce_fixed_feasibility) {
      if (have_start_best) {
        have_multistart_best = 1;
        fret_best = fret;
      } else {
        have_multistart_best = 0;
        fret_best = DBL_MAX;
      }
    } else {
      fret_best = fret;
    }
    vector_scale_factor_multistart = alloc_vecd(num_var + 1);

    for(i = 1; i <= num_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
    np_progress_bandwidth_multistart_step(1, iNum_Multistart);

    /* Conduct search from new random values of the search parameters */

    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and directions for NR modules */
				
      initialize_nr_vector_scale_factor(BANDWIDTH_reg_extern,
                                        1,        /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        0,
                                        0,
                                        0,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        0,
                                        KERNEL_reg_unordered_extern,
                                        int_use_starting_values,
                                        scale_cat,
                                        pow((double)4.0/(double)3.0,0.2),     /* Init for continuous vars */
                                        nconfac_extern, ncatfac_extern,
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor,
                                        lbc_init, hbc_init, c_init, 
                                        lbd_init, hbd_init, d_init,
                                        matrix_X_continuous_train_extern,
                                        matrix_Y_continuous_train_extern);

      initialize_nr_directions(BANDWIDTH_reg_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               0,
                               0,
                               0,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               1, int_RANDOM_SEED, 
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);


      /* Conduct direction set search */

      bwm_reset_counters();

      powell(0,
             0,
             vector_scale_factor,
             vector_scale_factor,
             matrix_y,
             num_var,
             ftol,
             tol,
             small,
             itmax,
             &iter,
             &fret,
             bwmfunc_wrapper);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE)	{
						
        initialize_nr_directions(BANDWIDTH_reg_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 num_reg_unordered_extern,
                                 num_reg_ordered_extern,
                                 0,
                                 0,
                                 0,
                                 vsfh,
                                 num_categories_extern,
                                 matrix_y,
                                 0, int_RANDOM_SEED, 
                                 lbc_dir, dfc_dir, c_dir, initc_dir,
                                 lbd_dir, hbd_dir, d_dir, initd_dir,
                                 matrix_X_continuous_train_extern,
                                 matrix_Y_continuous_train_extern);

						
        powell(0,
               0,
               vector_scale_factor,
               vector_scale_factor,
               matrix_y,
               num_var,
               ftol,
               tol,
               small,
               itmax,
               &iter,
               &fret,
               bwmfunc_wrapper);

      }

      /* If this run resulted in an improved minimum save information */

      if (enforce_fixed_feasibility) {
        if (np_bw_candidate_is_admissible(
              num_var,
              bwm_use_transform,
              KERNEL_reg_extern,
              KERNEL_reg_unordered_extern,
              BANDWIDTH_reg_extern,
              BANDWIDTH_reg_extern,
              0,
              num_obs_train_extern,
              0,
              0,
              0,
              num_reg_continuous_extern,
              num_reg_unordered_extern,
              num_reg_ordered_extern,
              num_categories_extern,
              vector_scale_factor) &&
            ((!have_multistart_best) || (fret < fret_best))) {
          fret_best = fret;
          have_multistart_best = 1;
          iImproved = iMs_counter+1;
          *timing = timing_extern;
          np_copy_scale_factor(vector_scale_factor_multistart, vector_scale_factor, num_var);
        }
      } else if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        *timing = timing_extern;
        
        for(i = 1; i <= num_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }
      objective_function_values[iMs_counter]=fret;
      objective_function_evals[iMs_counter]=bwm_eval_count;
      objective_function_invalid[iMs_counter]=bwm_invalid_count;
      bwm_snapshot_fast_counters();
      fast_eval_total += bwm_fast_eval_count;
      np_progress_bandwidth_multistart_step(iMs_counter+1, iNum_Multistart);

    }

    /* Save best for estimation */

    if (enforce_fixed_feasibility) {
      if (have_multistart_best) {
        fret = fret_best;
        np_copy_scale_factor(vector_scale_factor, vector_scale_factor_multistart, num_var);
        have_start_best = 1;
        fret_start_best = fret_best;
        np_copy_scale_factor(vector_scale_factor_startbest, vector_scale_factor_multistart, num_var);
      } else {
        have_start_best = 0;
        fret = DBL_MAX;
      }
    } else {
      fret = fret_best;

      for(i = 1; i <= num_var; i++)
        vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];
    }

    free(vector_scale_factor_multistart);

  }

  if (enforce_fixed_feasibility) {
    double final_raw;
    if (!have_start_best) {
      bw_error_msg = "C_np_regression_bw: optimizer failed to produce a feasible fixed-bandwidth candidate";
      goto cleanup_np_regression_bw_mode;
    }
    if (!np_bw_candidate_is_admissible(
          num_var,
          bwm_use_transform,
          KERNEL_reg_extern,
          KERNEL_reg_unordered_extern,
          BANDWIDTH_reg_extern,
          BANDWIDTH_reg_extern,
          0,
          num_obs_train_extern,
          0,
          0,
          0,
          num_reg_continuous_extern,
          num_reg_unordered_extern,
          num_reg_ordered_extern,
          num_categories_extern,
          vector_scale_factor)) {
      bw_error_msg = "C_np_regression_bw: optimizer returned an infeasible fixed-bandwidth candidate";
      goto cleanup_np_regression_bw_mode;
    }
    final_raw = bwmfunc_raw_current_scale(vector_scale_factor, num_var);
    if (!R_FINITE(final_raw) || final_raw == DBL_MAX) {
      bw_error_msg = "C_np_regression_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective";
      goto cleanup_np_regression_bw_mode;
    }
    fret = final_raw;
  }

  if (bwm_use_transform)
    bwm_to_constrained(vector_scale_factor, num_var);

  /* return data to R */
  if (BANDWIDTH_reg_extern == BW_GEN_NN || 
      BANDWIDTH_reg_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern; i++ )
      vector_scale_factor[i+1]=np_fround(vector_scale_factor[i+1]);
  }
  for( i=0; i<num_var; i++ )
    rbw[i]=vector_scale_factor[i+1];

  fval[0] = fret;
  fval[1] = iImproved;
  fval[2] = fret_initial;
  objective_function_fast[0] = fast_eval_total;
  /* end return data */

cleanup_np_regression_bw_mode:
  /* Free data objects */
  bwm_clear_floor_context();

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  safe_free(vector_Y_extern);

  free_mat(matrix_y, num_var + 1);
  safe_free(vector_scale_factor);
  safe_free(vector_scale_factor_startbest);
  safe_free(vsfh);
  safe_free(num_categories_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  free(vector_continuous_stddev);

  safe_free(ipt);
  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  np_glp_cv_clear_extern();
  np_reg_cv_core_clear_extern();

  int_cker_bound_extern = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;
  np_reset_y_side_extern();
  vector_glp_degree_extern = NULL;
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = 0;
  int_glp_basis_extern = 1;
  np_clear_estimator_extern_aliases();
  int_nn_k_min_extern = 1;

  if (bw_error_msg != NULL)
    error("%s", bw_error_msg);

  //fprintf(stderr,"\nNP TOASTY\n");
  return ;
  
}

void np_regression_bw(double * runo, double * rord, double * rcon, double * y,
                      double * mysd, int * myopti, double * myoptd, double * rbw, double * fval,
                      double * objective_function_values, double * objective_function_evals,
                      double * objective_function_invalid, double * timing,
                      double * objective_function_fast,
                      int * penalty_mode, double * penalty_mult,
                      int * glp_degree,
                      int * glp_bernstein,
                      int * glp_basis,
                      double * ckerlb, double * ckerub){
  int_nn_k_min_extern = 1;
  np_regression_bw_mode(runo, rord, rcon, y,
                        mysd, myopti, myoptd, rbw, fval,
                        objective_function_values, objective_function_evals,
                        objective_function_invalid, timing,
                        objective_function_fast,
                        penalty_mode, penalty_mult,
                        glp_degree, glp_bernstein, glp_basis,
                        ckerlb, ckerub, 0);
}


void np_regression(double * tuno, double * tord, double * tcon, double * ty,
                   double * euno, double * eord, double * econ, double * ey,
                   double * rbw, 
                   double * mcv, double * padnum, 
                   double * nconfac, double * ncatfac, double * mysd,
                   int * myopti, 
                   int * glp_degree,
                   int * glp_gradient_order,
                   int * glp_bernstein,
                   int * glp_basis,
                   double * cm, double * cmerr, double * g, double *gerr, 
                   double * xtra,
                   double * ckerlb, double * ckerub){

  double * vector_scale_factor, * ecm = NULL, * ecmerr = NULL, ** eg = NULL, **egerr = NULL;
  double * lambda, ** matrix_bandwidth;
  double RS, MSE, MAE, MAPE, CORR, SIGN, pad_num;

  int i,j, num_var;
  int ey_is_ty, do_grad, train_is_eval, num_obs_eval_alloc, max_lev, old_reg;

  int * ipt = NULL, * ipe = NULL;  // point permutation, see tree.c
  /* match integer options with their globals */

  num_reg_continuous_extern = myopti[REG_NCONI];
  num_reg_unordered_extern = myopti[REG_NUNOI];
  num_reg_ordered_extern = myopti[REG_NORDI];
  np_reset_y_side_extern();

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  train_is_eval = myopti[REG_TISEI];
  ey_is_ty = myopti[REG_EY];

  num_obs_train_extern = myopti[REG_TNOBSI];
  num_obs_eval_extern = myopti[REG_ENOBSI];

  if(train_is_eval && (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_regression): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    error("\n(np_regression): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
  }

  KERNEL_reg_extern = myopti[REG_CKRNEVI];
  KERNEL_reg_unordered_extern = myopti[REG_UKRNEVI];
  KERNEL_reg_ordered_extern = myopti[REG_OKRNEVI];

  vector_ckerlb_extern = ckerlb;
  vector_ckerub_extern = ckerub;
  int_cker_bound_extern = np_has_finite_cker_bounds(ckerlb, ckerub, num_reg_continuous_extern);

  int_LARGE_SF = myopti[REG_LSFI];
  int_MINIMIZE_IO = myopti[REG_MINIOI];
  BANDWIDTH_reg_extern = myopti[REG_BWI];

  do_grad = myopti[REG_GRAD];
  int_ll_extern = myopti[REG_LL];
  vector_glp_degree_extern = glp_degree;
  vector_glp_gradient_order_extern = glp_gradient_order;
  int_glp_bernstein_extern = *glp_bernstein;
  int_glp_basis_extern = *glp_basis;

  max_lev = myopti[REG_MLEVI];
  pad_num = *padnum;

  nconfac_extern = *nconfac;
  ncatfac_extern = *ncatfac;

  int_TREE_X = myopti[REG_DOTREEI];
  old_reg = myopti[REG_OLDREGI];

#ifdef MPI2
  num_obs_eval_alloc = MAX((int)ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif


  /* Allocate memory for objects */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  vector_Y_extern = alloc_vecd(num_obs_train_extern);

  if(!train_is_eval){
    matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
    matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
    matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);

    if(!ey_is_ty)
      vector_Y_eval_extern = alloc_vecd(num_obs_eval_extern);
    else
      vector_Y_eval_extern = NULL;

  } else {
    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;

    if(!ey_is_ty)
      vector_Y_eval_extern = alloc_vecd(num_obs_eval_extern);
    else
      vector_Y_eval_extern = vector_Y_extern;

  }

  ecm = alloc_vecd(num_obs_eval_alloc);
  ecmerr = alloc_vecd(num_obs_eval_alloc);
  

  eg = alloc_matd(num_obs_eval_alloc, num_var);
  egerr = alloc_matd(num_obs_eval_alloc, num_var);
  
  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  vector_scale_factor = alloc_vecd(num_var + 1);
  matrix_categorical_vals_extern = alloc_matd(max_lev, num_reg_unordered_extern + num_reg_ordered_extern);

  lambda =  alloc_vecd(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_bandwidth = alloc_matd((BANDWIDTH_reg_extern==BW_GEN_NN)?num_obs_eval_extern:
                                ((BANDWIDTH_reg_extern==BW_ADAP_NN)?num_obs_train_extern:1),num_reg_continuous_extern);  

  vector_continuous_stddev_extern = mysd;
  /* train */

  for( j=0;j<num_reg_unordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+i];

  for( i=0;i<num_obs_train_extern;i++ )
    vector_Y_extern[i] = ty[i];

  /* eval */
  if(!train_is_eval){
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+i];
  }

  if (!ey_is_ty)
    for(i=0;i<num_obs_eval_extern;i++)
      vector_Y_eval_extern[i] = ey[i];

  /*  bandwidths/scale factors */

  for( i=0; i<num_var; i++ )
    vector_scale_factor[i+1] = rbw[i];

  /* fix up categories */

  for(j=0; j < (num_reg_unordered_extern + num_reg_ordered_extern); j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;
  }

  ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  if(!train_is_eval) {
    ipe = (int *)malloc(num_obs_eval_extern*sizeof(int));
    if(!(ipe != NULL))
      error("!(ipe != NULL)");

    for(i = 0; i < num_obs_eval_extern; i++){
      ipe[i] = i;
    }
  } else {
    ipe = ipt;
  }

  // attempt tree build, if enabled 
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    if((BANDWIDTH_reg_extern != BW_ADAP_NN) || ((BANDWIDTH_reg_extern == BW_ADAP_NN) && train_is_eval)){
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipt, &kdt_extern_X);

      //put training data into tree-order using the index array

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+ipt[i]];
    
    
      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+ipt[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+ipt[i]];

      /* response variable */
      for( i=0;i<num_obs_train_extern;i++ )
        vector_Y_extern[i] = ty[ipt[i]];

    } else {
      build_kdtree(matrix_X_continuous_eval_extern, num_obs_eval_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipe, &kdt_extern_X);

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+ipe[i]];

      if(!ey_is_ty)
        for(i=0;i<num_obs_eval_extern;i++)
          vector_Y_eval_extern[i] = ey[ipe[i]];

    }
  }


  /* Conduct estimation */
	
  /* 
     nb - KERNEL_(|un)ordered_den are set to zero upon declaration 
     - they have only one kernel type each at the moment 
  */

  if(old_reg){
    kernel_estimate_regression_categorical(int_ll_extern,
                                           KERNEL_reg_extern,
                                           KERNEL_reg_unordered_extern,
                                           KERNEL_reg_ordered_extern,
                                           BANDWIDTH_reg_extern,
                                           num_obs_train_extern,
                                           num_obs_eval_extern,
                                           num_reg_unordered_extern,
                                           num_reg_ordered_extern,
                                           num_reg_continuous_extern,
                                           /* Train */
                                           matrix_X_unordered_train_extern,
                                           matrix_X_ordered_train_extern,
                                           matrix_X_continuous_train_extern,
                                           /* Eval */
                                           matrix_X_unordered_eval_extern,
                                           matrix_X_ordered_eval_extern,
                                           matrix_X_continuous_eval_extern,
                                           /* Bandwidth */
                                           matrix_X_continuous_train_extern,
                                           vector_Y_extern,
                                           vector_Y_eval_extern,
                                           &vector_scale_factor[1],
                                           num_categories_extern,
                                           ecm,
                                           eg,
                                           ecmerr,
                                           egerr,
                                           &RS,
                                           &MSE,
                                           &MAE,
                                           &MAPE,
                                           &CORR,
                                           &SIGN);

    if (do_grad){
      kernel_bandwidth_mean(KERNEL_reg_extern,
                            BANDWIDTH_reg_extern,
                            num_obs_train_extern,
                            num_obs_eval_extern,
                            0,
                            0,
                            0,
                            num_reg_continuous_extern,
                            num_reg_unordered_extern,
                            num_reg_ordered_extern,
                            np_mpi_local_regression_active() ? 1 : 0,
                            &vector_scale_factor[1],
                            /* Not used */
                            matrix_Y_continuous_train_extern,
                            /* Not used */
                            matrix_Y_continuous_train_extern,
                            matrix_X_continuous_train_extern,
                            matrix_X_continuous_eval_extern,
                            matrix_bandwidth,/* Not used */
                            matrix_bandwidth,
                            lambda);
      kernel_estimate_categorical_gradient_ocg_fast(1,
                                                    NULL,
                                                    0,
                                                    KERNEL_reg_extern,
                                                    KERNEL_reg_unordered_extern,
                                                    KERNEL_reg_ordered_extern,
                                                    BANDWIDTH_reg_extern,
                                                    int_ll_extern,
                                                    0,
                                                    num_obs_train_extern,
                                                    num_obs_eval_extern,
                                                    num_reg_unordered_extern,
                                                    num_reg_ordered_extern,
                                                    num_reg_continuous_extern,
                                                    vector_Y_extern,
                                                    matrix_X_unordered_train_extern,
                                                    matrix_X_ordered_train_extern,
                                                    matrix_X_continuous_train_extern,
                                                    matrix_X_unordered_eval_extern,
                                                    matrix_X_ordered_eval_extern,
                                                    matrix_X_continuous_eval_extern,
                                                    matrix_bandwidth,
                                                    NULL,
                                                    lambda,
                                                    num_categories_extern,
                                                    matrix_categorical_vals_extern,
                                                    ecm,
                                                    &eg[num_reg_continuous_extern]);

    }
  } else {

    kernel_estimate_regression_categorical_tree_np(int_ll_extern,
                                                   KERNEL_reg_extern,
                                                   KERNEL_reg_unordered_extern,
                                                   KERNEL_reg_ordered_extern,
                                                   BANDWIDTH_reg_extern,
                                                   num_obs_train_extern,
                                                   num_obs_eval_extern,
                                                   num_reg_unordered_extern,
                                                   num_reg_ordered_extern,
                                                   num_reg_continuous_extern,
                                                   /* Train */
                                                   matrix_X_unordered_train_extern,
                                                   matrix_X_ordered_train_extern,
                                                   matrix_X_continuous_train_extern,
                                                   /* Eval */
                                                   matrix_X_unordered_eval_extern,
                                                   matrix_X_ordered_eval_extern,
                                                   matrix_X_continuous_eval_extern,
                                                   vector_Y_extern,
                                                   vector_Y_eval_extern,
                                                   &vector_scale_factor[1],
                                                   num_categories_extern,
                                                   matrix_categorical_vals_extern,
                                                   ecm,
                                                   do_grad ? eg : NULL,
                                                   ecmerr,
                                                   do_grad ? egerr : NULL,
                                                   &RS,
                                                   &MSE,
                                                   &MAE,
                                                   &MAPE,
                                                   &CORR,
                                                   &SIGN);


  }

  for(i=0;i<num_obs_eval_extern;i++)
    cm[ipe[i]] = ecm[i];
      
  for(i=0;i<num_obs_eval_extern;i++)
    cmerr[ipe[i]] = ecmerr[i];


  if(do_grad){
    for(j=0;j<num_var;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        g[j*num_obs_eval_extern+ipe[i]]=eg[j][i];

    for(j=0;j<num_reg_continuous_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        gerr[j*num_obs_eval_extern+ipe[i]]=egerr[j][i];
  }


  /* write the return values */


  xtra[0] = RS;
  xtra[1] = MSE;
  xtra[2] = MAE;
  xtra[3] = MAPE;
  xtra[4] = CORR;
  xtra[5] = SIGN;

  /* clean up and wave goodbye */

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  if(!train_is_eval){
    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  safe_free(ipt);

  if(!train_is_eval)
    safe_free(ipe);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  free_mat(eg, num_var);
  free_mat(egerr, num_var);
  
  free_mat(matrix_bandwidth, num_reg_continuous_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  safe_free(vector_Y_extern);
  if(!ey_is_ty)
    safe_free(vector_Y_eval_extern);

  safe_free(ecm);
  safe_free(ecmerr);

  safe_free(num_categories_extern);
  safe_free(vector_scale_factor);
  vector_continuous_stddev_extern = NULL;

  safe_free(lambda);

  int_cker_bound_extern = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;
  np_reset_y_side_extern();
  vector_glp_degree_extern = NULL;
  vector_glp_gradient_order_extern = NULL;
  int_glp_bernstein_extern = 0;
  int_glp_basis_extern = 1;

  return;
}

void np_kernelsum(double * tuno, double * tord, double * tcon, 
                  double * ty, double * weights,
                  double * euno, double * eord, double * econ, 
                  double * bw,
                  double * mcv, double * padnum, 
                  int * operator,
                  int * myopti, double * kpow, 
                  double * weighted_sum, double * weighted_p_sum,
                  double * kernel_weights,
                  double * permutation_kernel_weights,
                  double * ckerlb, double * ckerub){

  int * ipt = NULL, * ipe = NULL;  // point permutation, see tree.c
      
  /* the ys are the weights */

  double * vector_scale_factor, * ksum, * p_ksum = NULL, pad_num, * kw = NULL, * pkw = NULL;
  int i,j,k, num_var, num_obs_eval_alloc;
  int no_y, leave_one_out, train_is_eval, do_divide_bw;
  int max_lev, no_weights, sum_element_length, return_kernel_weights;
  int p_operator, do_score, do_ocg, p_nvar = 0;
  int suppress_parallel_ksum = 0;

  int use_tree = 0;
  int allocated_X_train = 1, allocated_X_eval = 1;
  int allocated_Y = 1, allocated_W = 1;
  int ksum_is_output = 0, pksum_is_output = 0, kw_is_output = 0, pkw_is_output = 0;

  struct th_table * otabs = NULL;
  struct th_entry * ret = NULL;
  int ** matrix_ordered_indices = NULL;

  int ncol_Y, ncol_W;

  int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;

  int npks_err = 0;

  /* match integer options with their globals */

  num_reg_continuous_extern = myopti[KWS_NCONI];
  num_reg_unordered_extern = myopti[KWS_NUNOI];
  num_reg_ordered_extern = myopti[KWS_NORDI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[KWS_TNOBSI];
  num_obs_eval_extern = myopti[KWS_ENOBSI];

  KERNEL_reg_extern = myopti[KWS_CKRNEVI];
  KERNEL_reg_unordered_extern = myopti[KWS_UKRNEVI];
  KERNEL_reg_ordered_extern = myopti[KWS_OKRNEVI];
  vector_ckerlb_extern = ckerlb;
  vector_ckerub_extern = ckerub;
  int_cker_bound_extern = np_has_finite_cker_bounds(ckerlb, ckerub, num_reg_continuous_extern);

  int_LARGE_SF = myopti[KWS_LSFI];
  int_MINIMIZE_IO = myopti[KWS_MINIOI];
  BANDWIDTH_reg_extern = myopti[KWS_BWI];

  train_is_eval = myopti[KWS_TISEI];
  // no_y = myopti[KWS_NOYI];
  leave_one_out = myopti[KWS_LOOI];
  do_divide_bw = myopti[KWS_BDIVI];
  
  max_lev = myopti[KWS_MLEVI];
  pad_num = *padnum;

  /* the y and weight matrices will be contained in these variables */
  ncol_Y = myopti[KWS_YNCOLI];
  ncol_W = myopti[KWS_WNCOLI];

  int_TREE_X = myopti[KWS_DOTREEI];
  return_kernel_weights = myopti[KWS_RKWI];
  p_operator = myopti[KWS_POPI];
  do_score = myopti[KWS_PSCOREI];
  do_ocg = myopti[KWS_POCGI];
  suppress_parallel_ksum = myopti[KWS_SPARI];

  nconfac_extern = ncatfac_extern = 0.0;

  no_y = (ncol_Y == 0);
  no_weights = (ncol_W == 0);

  sum_element_length = (no_y ? 1 : ncol_Y)*(no_weights ? 1 : ncol_W);

  use_tree = (int_TREE_X == NP_TREE_TRUE);

#ifdef MPI2
  num_obs_eval_alloc = MAX(ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif

  if(train_is_eval && (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_kernelsum): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    error("\n(np_kernelsum): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
  }

  /* allocate */

  if(use_tree){
    matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
    matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
    matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);
  } else {
    allocated_X_train = 0;
    matrix_X_unordered_train_extern = (num_reg_unordered_extern > 0) ?
      (double **)R_alloc((size_t)num_reg_unordered_extern, sizeof(double*)) : NULL;
    matrix_X_ordered_train_extern = (num_reg_ordered_extern > 0) ?
      (double **)R_alloc((size_t)num_reg_ordered_extern, sizeof(double*)) : NULL;
    matrix_X_continuous_train_extern = (num_reg_continuous_extern > 0) ?
      (double **)R_alloc((size_t)num_reg_continuous_extern, sizeof(double*)) : NULL;

    for(j = 0; j < num_reg_unordered_extern; j++)
      matrix_X_unordered_train_extern[j] = tuno + j*num_obs_train_extern;
    for(j = 0; j < num_reg_ordered_extern; j++)
      matrix_X_ordered_train_extern[j] = tord + j*num_obs_train_extern;
    for(j = 0; j < num_reg_continuous_extern; j++)
      matrix_X_continuous_train_extern[j] = tcon + j*num_obs_train_extern;
  }
  
  /* for the moment we will just allocate a vector of ones */
  /* vector_Y_extern = (no_y)?NULL:alloc_vecd(num_obs_train_extern); */

  if(use_tree){
    matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, ncol_Y);
    matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, ncol_W);
  } else {
    allocated_Y = 0;
    allocated_W = 0;
    matrix_Y_continuous_train_extern = (ncol_Y > 0) ?
      (double **)R_alloc((size_t)ncol_Y, sizeof(double*)) : NULL;
    matrix_Y_ordered_train_extern = (ncol_W > 0) ?
      (double **)R_alloc((size_t)ncol_W, sizeof(double*)) : NULL;

    for(j = 0; j < ncol_Y; j++)
      matrix_Y_continuous_train_extern[j] = ty + j*num_obs_train_extern;
    for(j = 0; j < ncol_W; j++)
      matrix_Y_ordered_train_extern[j] = weights + j*num_obs_train_extern;
  }

  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_categorical_vals_extern = alloc_matd(max_lev, num_reg_unordered_extern + num_reg_ordered_extern);

  vector_scale_factor = alloc_vecd(num_var + 1);
  if(!use_tree && (num_obs_eval_alloc == num_obs_eval_extern)){
    ksum = weighted_sum;
    ksum_is_output = 1;
  } else {
    ksum = alloc_vecd(num_obs_eval_alloc*sum_element_length);
  }

  if((p_operator != OP_NOOP) || do_ocg){
    p_nvar = ((p_operator != OP_NOOP) ? num_reg_continuous_extern : 0) + ((do_score || do_ocg) ? num_reg_unordered_extern + num_reg_ordered_extern : 0);
    if(!use_tree && (num_obs_eval_alloc == num_obs_eval_extern)){
      p_ksum = weighted_p_sum;
      pksum_is_output = 1;
    } else {
      p_ksum = alloc_vecd(num_obs_eval_alloc*sum_element_length*p_nvar);
    }
  }

  if(!train_is_eval){
    if(use_tree){
      matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
      matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
      matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);
    } else {
      allocated_X_eval = 0;
      matrix_X_unordered_eval_extern = (num_reg_unordered_extern > 0) ?
        (double **)R_alloc((size_t)num_reg_unordered_extern, sizeof(double*)) : NULL;
      matrix_X_ordered_eval_extern = (num_reg_ordered_extern > 0) ?
        (double **)R_alloc((size_t)num_reg_ordered_extern, sizeof(double*)) : NULL;
      matrix_X_continuous_eval_extern = (num_reg_continuous_extern > 0) ?
        (double **)R_alloc((size_t)num_reg_continuous_extern, sizeof(double*)) : NULL;

      for(j = 0; j < num_reg_unordered_extern; j++)
        matrix_X_unordered_eval_extern[j] = euno + j*num_obs_eval_extern;
      for(j = 0; j < num_reg_ordered_extern; j++)
        matrix_X_ordered_eval_extern[j] = eord + j*num_obs_eval_extern;
      for(j = 0; j < num_reg_continuous_extern; j++)
        matrix_X_continuous_eval_extern[j] = econ + j*num_obs_eval_extern;
    }
  } else {
    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;
  }

  /* train */

  if(use_tree){
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+i];

    for( j = 0; j < ncol_Y; j++ )
      for( i = 0; i < num_obs_train_extern; i++ )
        matrix_Y_continuous_train_extern[j][i] = ty[j*num_obs_train_extern+i];

    for( j = 0; j < ncol_W; j++ )
      for( i = 0; i < num_obs_train_extern; i++ )
        matrix_Y_ordered_train_extern[j][i] = weights[j*num_obs_train_extern+i];
  }

  if(use_tree && !train_is_eval){
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+i];
  }

  if(use_tree){
    ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
    if(!(ipt != NULL))
      error("!(ipt != NULL)");

    for(i = 0; i < num_obs_train_extern; i++){
      ipt[i] = i;
    }

    if(!train_is_eval) {
      ipe = (int *)malloc(num_obs_eval_extern*sizeof(int));
      if(!(ipe != NULL))
        error("!(ipe != NULL)");

      for(i = 0; i < num_obs_eval_extern; i++){
        ipe[i] = i;
      }
    } else {
      ipe = ipt;
    }
  } else {
    ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
    if(!(ipt != NULL))
      error("!(ipt != NULL)");

    for(i = 0; i < num_obs_train_extern; i++){
      ipt[i] = i;
    }

    if(!train_is_eval){
      ipe = (int *)malloc(num_obs_eval_extern*sizeof(int));
      if(!(ipe != NULL))
        error("!(ipe != NULL)");

      for(i = 0; i < num_obs_eval_extern; i++){
        ipe[i] = i;
      }
    } else {
      ipe = ipt;
    }
  }

  // attempt tree build, if enabled 
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    if((BANDWIDTH_reg_extern != BW_ADAP_NN) || ((BANDWIDTH_reg_extern == BW_ADAP_NN) && train_is_eval)){
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipt, &kdt_extern_X);

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+ipt[i]];
    
    
      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+ipt[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+ipt[i]];

      for( j = 0; j < ncol_Y; j++ )
        for( i = 0; i < num_obs_train_extern; i++ )
          matrix_Y_continuous_train_extern[j][i] = ty[j*num_obs_train_extern+ipt[i]];

      for( j = 0; j < ncol_W; j++ )
        for( i = 0; i < num_obs_train_extern; i++ )
          matrix_Y_ordered_train_extern[j][i] = weights[j*num_obs_train_extern+ipt[i]];

    } else {
      build_kdtree(matrix_X_continuous_eval_extern, num_obs_eval_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipe, &kdt_extern_X);


      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+ipe[i]];

    }
  }



  /* bandwidths */
  for( i=0; i<num_var; i++ )
    vector_scale_factor[i+1] = bw[i];
  
  /* fix up categories */

  for(j=0; j<num_reg_unordered_extern; j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;
  }

  if(do_ocg && (num_reg_ordered_extern > 0)){
    otabs = (struct th_table *)malloc(num_reg_ordered_extern*sizeof(struct th_table));
    matrix_ordered_indices = (int **)malloc(num_reg_ordered_extern*sizeof(int *));
    int * tc = (int *)malloc(num_reg_ordered_extern*num_obs_eval_extern*sizeof(int));
    for(i = 0; i < num_reg_ordered_extern; i++)
      matrix_ordered_indices[i] = tc + i*num_obs_eval_extern;
  }

  for(j=num_reg_unordered_extern, k = 0; j < (num_reg_unordered_extern+num_reg_ordered_extern); j++, k++){
    i = 0;

    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;

    if(do_ocg){
      if(thcreate_r((size_t)ceil(1.6*num_categories_extern[j]), otabs + k) == TH_ERROR)
        error("hash table creation failed");

      for(i = 0; i < num_categories_extern[j]; i++){
        struct th_entry centry;
        centry.key.dkey = matrix_categorical_vals_extern[j][i];
        centry.data = i;

        if(thsearch_r(&centry, TH_ENTER, &ret, otabs+k) == TH_FAILURE)
          error("insertion into hash table failed");
      }
      
      // now do lookups
      struct th_entry te;
      te.key.dkey = pad_num;
      te.data = -1;

      ret = &te;
      for(i = 0; i < num_obs_eval_extern; i++){
        if(ret->key.dkey != matrix_X_ordered_eval_extern[k][i]){
          te.key.dkey = matrix_X_ordered_eval_extern[k][i];
          if(thsearch_r(&te, TH_SEARCH, &ret, otabs+k) == TH_FAILURE)
            error("hash table lookup failed (which should be impossible)");

        } 

        matrix_ordered_indices[k][i] = ret->data;
      }
    }
  }

  if(return_kernel_weights){
    if(!use_tree && (BANDWIDTH_reg_extern != BW_ADAP_NN) && (num_obs_eval_alloc == num_obs_eval_extern)){
      kw = kernel_weights;
      kw_is_output = 1;
    } else {
      kw = alloc_vecd(num_obs_train_extern*num_obs_eval_extern);
    }
  }

  if(return_kernel_weights && (p_nvar > 0)){
    pkw = alloc_vecd(num_obs_train_extern*num_obs_eval_extern*p_nvar);
  }
  
  kernel_c = (int *)malloc(sizeof(int)*num_reg_continuous_extern);

  for(i = 0; i < num_reg_continuous_extern; i++)
    kernel_c[i] = KERNEL_reg_extern;

  kernel_u = (int *)malloc(sizeof(int)*num_reg_unordered_extern);

  for(i = 0; i < num_reg_unordered_extern; i++)
    kernel_u[i] = KERNEL_reg_unordered_extern;

  kernel_o = (int *)malloc(sizeof(int)*num_reg_ordered_extern);

  for(i = 0; i < num_reg_ordered_extern; i++)
    kernel_o[i] = KERNEL_reg_ordered_extern;

  int *bpso = NULL;
  {
    const int bpso_len = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
    if(bpso_len > 0){
      const int do_perm = (p_operator != OP_NOOP);
      const int doscoreocg = (do_score || do_ocg);

      bpso = (int *)R_alloc((size_t)bpso_len, sizeof(int));
      for(i = 0; i < num_reg_continuous_extern; i++)
        bpso[i] = do_perm;
      for(i = num_reg_continuous_extern; i < bpso_len; i++)
        bpso[i] = doscoreocg;
    }
  }
  
  
  npks_err = kernel_weighted_sum_np(kernel_c,
                                      kernel_u,
                                      kernel_o,
                                      BANDWIDTH_reg_extern,
                                      num_obs_train_extern,
                                      num_obs_eval_extern,
                                      num_reg_unordered_extern,
                                      num_reg_ordered_extern,
                                      num_reg_continuous_extern,
                                      leave_one_out,
                                      0,
                                      (int)(*kpow),
                                      do_divide_bw,
                                      (BANDWIDTH_reg_extern == BW_ADAP_NN) ? do_divide_bw : 0,
                                      0, //not symmetric
                                      0, //disable 'twisting'
                                      0, // do not drop train
                                      0, // do not drop train
                                      operator,
                                      p_operator,
                                      do_score,
                                      do_ocg, // no ocg (for now)
                                      bpso,
                                      suppress_parallel_ksum,
                                      ncol_Y,
                                      ncol_W,
                                      int_TREE_X,
                                      0,
                                      kdt_extern_X,
                                      NULL, NULL, NULL,
                                      matrix_X_unordered_train_extern,
                                      matrix_X_ordered_train_extern,
                                      matrix_X_continuous_train_extern,
                                      matrix_X_unordered_eval_extern,
                                      matrix_X_ordered_eval_extern,
                                      matrix_X_continuous_eval_extern,
                                      /* ys matrix */
                                      matrix_Y_continuous_train_extern,
                                      /* weights matrix */
                                      matrix_Y_ordered_train_extern,
                                      NULL,
                                      &vector_scale_factor[1],
                                      0,NULL,NULL,NULL,
                                      num_categories_extern,
                                      matrix_categorical_vals_extern,
                                      matrix_ordered_indices,
                                      ksum,
                                      p_ksum,
                                      kw,
                                      pkw);
  if(npks_err != 0){
    error("kernel_weighted_sum_np failed with code %d", npks_err);
  }


  if(!npks_err){
    if(use_tree || !ksum_is_output){
      for(j = 0; j < num_obs_eval_extern; j++)
        for(i = 0; i < sum_element_length; i++)
          weighted_sum[ipe[j]*sum_element_length + i] = ksum[j*sum_element_length+i];
    }

    if(return_kernel_weights){
      if(use_tree || (BANDWIDTH_reg_extern == BW_ADAP_NN) || !kw_is_output){
        if(BANDWIDTH_reg_extern != BW_ADAP_NN){ // adaptive weights are currently returned transposed...
          for(j = 0; j < num_obs_eval_extern; j++)
            for(i = 0; i < num_obs_train_extern; i++)
              kernel_weights[ipe[j]*num_obs_train_extern + ipt[i]] = kw[j*num_obs_train_extern + i];
        } else {
          for(j = 0; j < num_obs_train_extern; j++)
            for(i = 0; i < num_obs_eval_extern; i++)
              kernel_weights[ipe[i]*num_obs_train_extern + ipt[j]] = kw[j*num_obs_eval_extern + i];      
        }
      }
    }

    if(return_kernel_weights && (p_nvar > 0) && (pkw != NULL)){
      if(BANDWIDTH_reg_extern != BW_ADAP_NN){
        for(k = 0; k < p_nvar; k++){
          const int koff = k*num_obs_train_extern*num_obs_eval_extern;
          for(j = 0; j < num_obs_eval_extern; j++)
            for(i = 0; i < num_obs_train_extern; i++)
              permutation_kernel_weights[koff + ipe[j]*num_obs_train_extern + ipt[i]] =
                pkw[koff + j*num_obs_train_extern + i];
        }
      } else {
        for(k = 0; k < p_nvar; k++){
          const int koff = k*num_obs_train_extern*num_obs_eval_extern;
          for(j = 0; j < num_obs_train_extern; j++)
            for(i = 0; i < num_obs_eval_extern; i++)
              permutation_kernel_weights[koff + ipe[i]*num_obs_train_extern + ipt[j]] =
                pkw[koff + j*num_obs_eval_extern + i];
        }
      }
    }

    if(p_nvar > 0){
      if(use_tree || !pksum_is_output){
        for(k = 0; k < p_nvar; k++){
          const int kidx = k*num_obs_eval_extern*sum_element_length;
          for(j = 0; j < num_obs_eval_extern; j++)
            for(i = 0; i < sum_element_length; i++)
              weighted_p_sum[kidx + ipe[j]*sum_element_length + i] = p_ksum[kidx + j*sum_element_length + i];
        }
      }
    }

  }
  /* clean up */

  if(allocated_X_train){
    free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  }
  
  if(!train_is_eval && allocated_X_eval){
    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);
  if(allocated_Y)
    free_mat(matrix_Y_continuous_train_extern, ncol_Y);
  if(allocated_W)
    free_mat(matrix_Y_ordered_train_extern, ncol_W);

  safe_free(num_categories_extern);
  safe_free(vector_scale_factor);
  if(!ksum_is_output)
    safe_free(ksum);

  if(!kw_is_output)
    safe_free(kw);
  if(!pkw_is_output)
    safe_free(pkw);

  safe_free(ipt);

  free(kernel_c);
  free(kernel_u);
  free(kernel_o);

  if(!train_is_eval)
    safe_free(ipe);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  if(p_nvar > 0){
    if(!pksum_is_output)
      safe_free(p_ksum);
  }

  if(do_ocg && (num_reg_ordered_extern > 0)){

    for(k = 0; k < num_reg_ordered_extern; k++){
      thdestroy_r(otabs+k);
    }
    free(otabs);
    free(matrix_ordered_indices[0]);
    free(matrix_ordered_indices);
  }

  int_cker_bound_extern = 0;
  vector_ckerlb_extern = NULL;
  vector_ckerub_extern = NULL;

  return;
}


void np_quantile_conditional(double * tc_con,
                             double * tu_uno, double * tu_ord, double * tu_con,
                             double * eu_uno, double * eu_ord, double * eu_con,
                             double * quantile,
                             double * mybw, 
                             double * mcv, double *padnum,
                             double * nconfac, double * ncatfac, double * mysd,
                             int * myopti, double * myoptd,
                             double * yq, double * yqerr, double *yg){
  /* Likelihood bandwidth selection for density estimation */

  double **g = NULL, * eq, * eqerr;
  double ftol, small, tol;
  double pad_num;

  int i,j, max_lev;
  int num_var, num_obs_eval_alloc;
  int num_all_var, num_var_var, train_is_eval, do_gradients;
  int itmax;

  int iNum_Multistart;

  double lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  int dfc_dir;

  iNum_Multistart = myopti[CQ_NMULTII];
  imsnum = 0;
  imstot = myopti[CQ_NMULTII]; /* iNum_Multistart */

  num_var_unordered_extern = 0;
  num_var_ordered_extern = 0;
  num_var_continuous_extern = 1;

  num_reg_unordered_extern = myopti[CQ_UNUNOI];
  num_reg_ordered_extern = myopti[CQ_UNORDI];
  num_reg_continuous_extern = myopti[CQ_UNCONI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;
  num_var_var = 1;
  num_all_var = num_var + num_var_var;

  num_obs_train_extern = myopti[CQ_TNOBSI];
  num_obs_eval_extern = myopti[CQ_ENOBSI];

  if((train_is_eval = myopti[CQ_TISEI]) && 
     (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_quantile_conditional): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    error("\n(np_quantile_conditional): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
  }

  KERNEL_reg_extern = myopti[CQ_CXKRNEVI];
  KERNEL_den_extern = myopti[CQ_CYKRNEVI];

  KERNEL_reg_unordered_extern = myopti[CQ_UXKRNEVI];
  KERNEL_den_unordered_extern = myopti[CQ_UYKRNEVI];

  KERNEL_reg_ordered_extern = myopti[CQ_OXKRNEVI];
  KERNEL_den_ordered_extern = myopti[CQ_OYKRNEVI];

  int_LARGE_SF = myopti[CQ_LSFI];
  BANDWIDTH_den_extern = myopti[CQ_DENI];
  int_MINIMIZE_IO = myopti[CQ_MINIOI];
  do_gradients = myopti[CQ_GRADI];
  itmax = myopti[CQ_ITMAXI];

  max_lev = myopti[CQ_MLEVI];
  pad_num = *padnum;

  ftol = myoptd[CQ_FTOLD];
  tol = myoptd[CQ_TOLD];
  small = myoptd[CQ_SMALLD];

  dfc_dir = myopti[CQ_DFC_DIRI];
  lbc_dir = myoptd[CQ_LBC_DIRD];
  c_dir = myoptd[CQ_C_DIRD];
  initc_dir = myoptd[CQ_INITC_DIRD]; 

  lbd_dir = myoptd[CQ_LBD_DIRD]; 
  hbd_dir = myoptd[CQ_HBD_DIRD]; 
  d_dir = myoptd[CQ_D_DIRD]; 
  initd_dir = myoptd[CQ_INITD_DIRD]; 

  gamma_extern = *quantile;
  nconfac_extern = *nconfac;
  ncatfac_extern = *ncatfac;
  vector_continuous_stddev_extern = mysd;

#ifdef MPI2
  num_obs_eval_alloc = MAX(ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif


  /* Allocate memory for objects */
  matrix_Y_continuous_quantile_extern = alloc_matd(1, num_var_continuous_extern);
  matrix_X_unordered_quantile_extern = alloc_matd(1, num_reg_unordered_extern);
  matrix_X_ordered_quantile_extern = alloc_matd(1, num_reg_ordered_extern);
  matrix_X_continuous_quantile_extern = alloc_matd(1, num_reg_continuous_extern);
  /* */

  matrix_Y_unordered_train_extern = alloc_matd(num_obs_train_extern, num_var_unordered_extern);
  matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);
  matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  if(train_is_eval) {
    matrix_Y_unordered_eval_extern = alloc_matd(num_obs_train_extern, num_var_unordered_extern);
    matrix_Y_ordered_eval_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);
    matrix_Y_continuous_eval_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);

    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;
  } else {
    matrix_Y_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_var_unordered_extern);
    matrix_Y_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_var_ordered_extern);
    matrix_Y_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_var_continuous_extern);

    matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
    matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
    matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);
  }
	
  num_categories_extern = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                     num_reg_unordered_extern + num_reg_ordered_extern);
  vector_scale_factor_extern = alloc_vecd(num_all_var + 1);
  
  matrix_categorical_vals_extern = 
    alloc_matd(max_lev, num_var_unordered_extern + num_var_ordered_extern + 
               num_reg_unordered_extern + num_reg_ordered_extern);

  eq = alloc_vecd(num_obs_eval_alloc);
  eqerr = alloc_vecd(num_obs_eval_alloc);

  if(do_gradients) {
    g = alloc_matd(num_obs_eval_alloc, num_var);
    for(j = 0; j < num_var; j++)
      for(i = 0; i < num_obs_eval_alloc; i++)
        g[j][i] = 0.0;
  }

  /* in v_s_f order is creg, cvar, uvar, ovar, ureg, oreg  */

  for( i=0;i<num_all_var; i++ )
    vector_scale_factor_extern[i+1] = mybw[i];

  /* Parse data */

  /* train */

  for(j=0;j<num_var_continuous_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_continuous_train_extern[j][i]=tc_con[j*num_obs_train_extern+i];

  for(j=0;j<num_reg_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_X_unordered_train_extern[j][i]=tu_uno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=tu_ord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=tu_con[j*num_obs_train_extern+i];

  /* eval */
  if(!train_is_eval){
    for(j=0;j<num_reg_unordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_X_unordered_eval_extern[j][i]=eu_uno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=eu_ord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=eu_con[j*num_obs_eval_extern+i];
  }

  /* fix up categories */
  
  for(j=0; j < (num_reg_unordered_extern + num_reg_ordered_extern); j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;
  }

  kernel_estimate_quantile(do_gradients,
                           KERNEL_den_extern,
                           KERNEL_den_unordered_extern,
                           KERNEL_den_ordered_extern,
                           BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_obs_eval_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           num_var_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           num_reg_continuous_extern,
                           matrix_Y_unordered_train_extern,
                           matrix_Y_ordered_train_extern,
                           matrix_Y_continuous_train_extern,
                           matrix_Y_unordered_eval_extern,
                           matrix_Y_ordered_eval_extern,
                           matrix_Y_continuous_eval_extern,
                           matrix_X_unordered_train_extern,
                           matrix_X_ordered_train_extern,
                           matrix_X_continuous_train_extern,
                           matrix_X_unordered_eval_extern,
                           matrix_X_ordered_eval_extern,
                           matrix_X_continuous_eval_extern,
                           &vector_scale_factor_extern[1],
                           eq,
                           eqerr,
                           g,
                           int_RANDOM_SEED,
                           ftol,
                           tol,
                           small,
                           itmax,
                           iNum_Multistart,            /* Maximum number of multistarts */
                           1.0e-10,         /* Zero for all intents and purposes */
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir);

  /* return data to R */

  for(i=0; i < num_obs_eval_extern; i++)
    yq[i] = eq[i];

  for(i=0; i < num_obs_eval_extern; i++)
    yqerr[i] = eqerr[i];

  if(do_gradients)
    for(j=0; j < num_var; j++)
      for(i=0; i < num_obs_eval_extern; i++)
        yg[j*num_obs_eval_extern+i] = g[j][i];
  
  /* end return data */

  /* Free data objects */

  free_mat(matrix_Y_unordered_train_extern, num_var_unordered_extern);
  free_mat(matrix_Y_ordered_train_extern, num_var_ordered_extern);
  free_mat(matrix_Y_continuous_train_extern, num_var_continuous_extern);

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  free_mat(matrix_Y_continuous_quantile_extern, num_var_continuous_extern); 
  free_mat(matrix_X_unordered_quantile_extern, num_reg_unordered_extern); 
  free_mat(matrix_X_ordered_quantile_extern, num_reg_ordered_extern); 
  free_mat(matrix_X_continuous_quantile_extern, num_reg_continuous_extern); 

  free_mat(g, num_var);

  free_mat(matrix_Y_unordered_eval_extern, num_var_unordered_extern);
  free_mat(matrix_Y_ordered_eval_extern, num_var_ordered_extern);
  free_mat(matrix_Y_continuous_eval_extern, num_var_continuous_extern);


  if(train_is_eval){
    matrix_X_unordered_eval_extern = NULL;
    matrix_X_ordered_eval_extern = NULL;
    matrix_X_continuous_eval_extern = NULL;
  } else {
    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  safe_free(vector_scale_factor_extern);
  safe_free(num_categories_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);

  safe_free(eq);
  safe_free(eqerr);
  np_clear_estimator_extern_aliases();

  return ;
}
