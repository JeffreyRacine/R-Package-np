/* This module contains the functions for the kernel bandwidth function. */

/* Copyright (C) J. Racine, 1995-2001 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <errno.h>

// timing tests
#include <time.h>

#include "headers.h"
#include "continuous_kernel_row.h"
#include "jksum_gaussian_density.h"
#include "kernel_registry.h"

#include <R.h>
#include <Rmath.h>

#if defined(__GNUC__) || defined(__clang__)
#define NP_KERNELCV_NOINLINE __attribute__((noinline))
#else
#define NP_KERNELCV_NOINLINE
#endif

#ifdef MPI2

#include "mpi.h"

extern  int my_rank;
extern  int source;
extern  int dest;
extern  int tag;
extern  int iNum_Processors;
extern  int iSeed_my_rank;
extern  MPI_Status status;
#endif

/*
int int_LARGE_SF;
int int_DEBUG;
int int_VERBOSE;
int int_NOKEYPRESS;
int int_DISPLAY_CV;
int int_RANDOM_SEED;
int int_MINIMIZE_IO;
int int_ORDERED_CATEGORICAL_GRADIENT;
int int_PREDICT;
int int_ROBUST;
int int_SIMULATION;
int int_TAYLOR;
int int_WEIGHTS;
*/

/* Some externals for numerical routines */

extern int num_obs_train_extern;
extern int num_obs_eval_extern;
extern int num_var_continuous_extern;
extern int num_var_unordered_extern;
extern int num_var_ordered_extern;
extern int num_reg_continuous_extern;
extern int num_reg_unordered_extern;
extern int num_reg_ordered_extern;
extern int *num_categories_extern;
extern double **matrix_categorical_vals_extern;

extern double **matrix_X_continuous_train_extern;
extern double **matrix_X_unordered_train_extern;
extern double **matrix_X_ordered_train_extern;
extern double **matrix_X_continuous_eval_extern;
extern double **matrix_X_unordered_eval_extern;
extern double **matrix_X_ordered_eval_extern;

extern double **matrix_Y_continuous_train_extern;
extern double **matrix_Y_unordered_train_extern;
extern double **matrix_Y_ordered_train_extern;
extern double **matrix_Y_continuous_eval_extern;
extern double **matrix_Y_unordered_eval_extern;
extern double **matrix_Y_ordered_eval_extern;

extern double *vector_Y_extern;
extern double *vector_lsq_scale_extern;
extern double *vector_lsq_loss_extern;
extern double *vector_lsq_q_extern;
extern double np_lsq_tau_extern;
extern double np_lsq_delta_lower_extern;
extern double np_lsq_delta_upper_extern;
extern double *vector_T_extern;
extern double *vector_Y_eval_extern;

/* Quantile - no Y ordered or unordered used, but defined anyways */

extern double **matrix_Y_continuous_quantile_extern;
extern double **matrix_Y_unordered_quantile_extern;
extern double **matrix_Y_ordered_quantile_extern;
extern double **matrix_X_continuous_quantile_extern;
extern double **matrix_X_unordered_quantile_extern;
extern double **matrix_X_ordered_quantile_extern;

extern int np_lp_engine_extern;

extern int KERNEL_reg_extern;
extern int KERNEL_reg_unordered_extern;
extern int KERNEL_reg_ordered_extern;
extern int KERNEL_den_extern;
extern int KERNEL_den_unordered_extern;
extern int KERNEL_den_ordered_extern;
extern int BANDWIDTH_reg_extern;
extern int BANDWIDTH_den_extern;

extern int itmax_extern;
extern double small_extern;
extern double gamma_extern;
extern double *vector_scale_factor_extern;

extern double y_min_extern;
extern double y_max_extern;

// cdens + trees extern
extern double **matrix_XY_continuous_train_extern;
extern double **matrix_XY_unordered_train_extern;
extern double **matrix_XY_ordered_train_extern;
extern double **matrix_XY_continuous_eval_extern;
extern double **matrix_XY_unordered_eval_extern;
extern double **matrix_XY_ordered_eval_extern;

// cdf extern
extern double dbl_memfac_ccdf_extern;
extern double dbl_memfac_dls_extern;
extern int cdfontrain_extern;

// timing
extern double timing_extern;
extern int np_beta_bw_order_extern;
extern int np_regression_bw_categorical_compress_extern;
extern int np_density_bw_categorical_compress_extern;
extern int np_distribution_bw_categorical_compress_extern;
extern int np_conditional_density_bw_categorical_compress_extern;
extern int np_conditional_distribution_bw_categorical_compress_extern;
extern int np_beta_cx_bw_order_extern;
extern int np_beta_cy_bw_order_extern;
extern double *vector_ckerlb_extern;
extern double *vector_ckerub_extern;
extern double *vector_cxkerlb_extern;
extern double *vector_cxkerub_extern;
extern double *vector_cykerlb_extern;
extern double *vector_cykerub_extern;
extern int int_bounded_cvls_quadrature_points_extern;

#define NP_LP_ENGINE_SCALAR  0

#define BW_FIXED   0
#define BW_GEN_NN  1
#define BW_ADAP_NN 2

static int np_beta_conditional_bw_active(void)
{
  return KERNEL_reg_extern == NP_CKERNEL_COORDINATE_CODE ||
    KERNEL_den_extern == NP_CKERNEL_COORDINATE_CODE;
}

static void np_beta_conditional_bw_route_or_error(
  NPContinuousKernelRoute * const route,
  NPContinuousKernelDerivativeDiagnostics * const diagnostics,
  const int order,
  const int num_continuous,
  const double * const lower,
  const double * const upper,
  const char * const where)
{
  route->segment_count = 1;
  route->segment[0].descriptor.family = NP_CKERNEL_FAMILY_BETA;
  route->segment[0].descriptor.legacy_code = NP_CKERNEL_COORDINATE_CODE;
  route->segment[0].descriptor.order = order;
  route->segment[0].coordinate_offset = 0;
  route->segment[0].coordinate_count = num_continuous;
  route->segment[0].lower = lower;
  route->segment[0].upper = upper;
  if(np_continuous_kernel_route_validate(route, num_continuous) !=
     NP_CKERNEL_ROUTE_OK)
    error("%s beta route has an invalid layout", where);
  diagnostics->bad_coordinate = -1;
  diagnostics->bad_observation = -1;
  diagnostics->undefined_count = 0;
  diagnostics->beta_status = NP_BETA_OK;
}

/*
 * Keep route construction and the expanded context call outside the shared
 * CVML wrapper.  The beta branch is selected before entry, while the legacy
 * wrapper retains its established stack frame and hot instruction layout.
 */
static NP_KERNELCV_NOINLINE double
np_beta_conditional_density_bw_objective_ctx(double *vector_scale_factor)
{
  double cv = 0.0;
  NPContinuousKernelRoute beta_x_route;
  NPContinuousKernelRoute beta_y_route;
  NPContinuousKernelDerivativeDiagnostics beta_x_diagnostics;
  NPContinuousKernelDerivativeDiagnostics beta_y_diagnostics;
  NPConditionalKernelExecutionContext execution_context = {0};

  if(KERNEL_reg_extern == NP_CKERNEL_COORDINATE_CODE) {
    np_beta_conditional_bw_route_or_error(
      &beta_x_route, &beta_x_diagnostics,
      np_beta_cx_bw_order_extern, num_reg_continuous_extern,
      vector_cxkerlb_extern, vector_cxkerub_extern,
      "conditional density CVML X");
    execution_context.x_route = &beta_x_route;
    execution_context.x_diagnostics = &beta_x_diagnostics;
  }
  if(KERNEL_den_extern == NP_CKERNEL_COORDINATE_CODE) {
    np_beta_conditional_bw_route_or_error(
      &beta_y_route, &beta_y_diagnostics,
      np_beta_cy_bw_order_extern, num_var_continuous_extern,
      vector_cykerlb_extern, vector_cykerub_extern,
      "conditional density CVML Y");
    execution_context.y_route = &beta_y_route;
    execution_context.y_diagnostics = &beta_y_diagnostics;
  }
  execution_context.categorical_compress =
    np_conditional_density_bw_categorical_compress_extern;
  if(np_kernel_estimate_con_density_categorical_leave_one_out_cv_ctx(
       KERNEL_den_extern,
       KERNEL_den_unordered_extern,
       KERNEL_den_ordered_extern,
       KERNEL_reg_extern,
       KERNEL_reg_unordered_extern,
       KERNEL_reg_ordered_extern,
       BANDWIDTH_den_extern,
       num_obs_train_extern,
       num_var_unordered_extern,
       num_var_ordered_extern,
       num_var_continuous_extern,
       num_reg_unordered_extern,
       num_reg_ordered_extern,
       num_reg_continuous_extern,
       matrix_Y_unordered_train_extern,
       matrix_Y_ordered_train_extern,
       matrix_Y_continuous_train_extern,
       matrix_X_unordered_train_extern,
       matrix_X_ordered_train_extern,
       matrix_X_continuous_train_extern,
       matrix_XY_unordered_train_extern,
       matrix_XY_ordered_train_extern,
       matrix_XY_continuous_train_extern,
       &vector_scale_factor[1],
       num_categories_extern,
       &execution_context,
       &cv) == 1)
    return DBL_MAX;
  return cv;
}

/*
 * Conditional-density CVLS uses the same immutable X/Y route descriptors as
 * CVML, but enters the shared quadrature/analytic CVLS adapter.  Keeping this
 * beta-only callback out of line preserves the incumbent legacy callback and
 * makes failure after route selection terminal rather than a sidecar fallback.
 */
static NP_KERNELCV_NOINLINE double
np_beta_conditional_density_bw_objective_ls_ctx(
  double *vector_scale_factor)
{
  double cv = 0.0;
  NPContinuousKernelRoute beta_x_route;
  NPContinuousKernelRoute beta_y_route;
  NPContinuousKernelDerivativeDiagnostics beta_x_diagnostics;
  NPContinuousKernelDerivativeDiagnostics beta_y_diagnostics;
  NPConditionalKernelExecutionContext execution_context = {0};

  if(KERNEL_reg_extern == NP_CKERNEL_COORDINATE_CODE) {
    np_beta_conditional_bw_route_or_error(
      &beta_x_route, &beta_x_diagnostics,
      np_beta_cx_bw_order_extern, num_reg_continuous_extern,
      vector_cxkerlb_extern, vector_cxkerub_extern,
      "conditional density CVLS X");
    execution_context.x_route = &beta_x_route;
    execution_context.x_diagnostics = &beta_x_diagnostics;
  }
  if(KERNEL_den_extern == NP_CKERNEL_COORDINATE_CODE) {
    np_beta_conditional_bw_route_or_error(
      &beta_y_route, &beta_y_diagnostics,
      np_beta_cy_bw_order_extern, num_var_continuous_extern,
      vector_cykerlb_extern, vector_cykerub_extern,
      "conditional density CVLS Y");
    execution_context.y_route = &beta_y_route;
    execution_context.y_diagnostics = &beta_y_diagnostics;
  }
  execution_context.categorical_compress =
    np_conditional_density_bw_categorical_compress_extern;
  if(np_conditional_density_cvls_lp_stream_ctx(
       &vector_scale_factor[1], &execution_context, &cv) != 0)
    return DBL_MAX;
  return cv;
}

/*
 * Conditional-distribution CVLS shares the immutable X/Y descriptors and
 * canonical signed LP row owner with conditional density.  The response
 * provider changes only the operator to OP_INTEGRAL inside jksum.c; the
 * callback owns no estimator algebra and has no beta-sidecar fallback.
 */
static NP_KERNELCV_NOINLINE double
np_beta_conditional_distribution_bw_objective_ls_ctx(
  double *vector_scale_factor)
{
  double cv = 0.0;
  NPContinuousKernelRoute beta_x_route;
  NPContinuousKernelRoute beta_y_route;
  NPContinuousKernelDerivativeDiagnostics beta_x_diagnostics;
  NPContinuousKernelDerivativeDiagnostics beta_y_diagnostics;
  NPConditionalKernelExecutionContext execution_context = {0};

  if(KERNEL_reg_extern == NP_CKERNEL_COORDINATE_CODE) {
    np_beta_conditional_bw_route_or_error(
      &beta_x_route, &beta_x_diagnostics,
      np_beta_cx_bw_order_extern, num_reg_continuous_extern,
      vector_cxkerlb_extern, vector_cxkerub_extern,
      "conditional distribution CVLS X");
    execution_context.x_route = &beta_x_route;
    execution_context.x_diagnostics = &beta_x_diagnostics;
  }
  if(KERNEL_den_extern == NP_CKERNEL_COORDINATE_CODE) {
    np_beta_conditional_bw_route_or_error(
      &beta_y_route, &beta_y_diagnostics,
      np_beta_cy_bw_order_extern, num_var_continuous_extern,
      vector_cykerlb_extern, vector_cykerub_extern,
      "conditional distribution CVLS Y");
    execution_context.y_route = &beta_y_route;
    execution_context.y_diagnostics = &beta_y_diagnostics;
  }
  execution_context.categorical_compress =
    np_conditional_distribution_bw_categorical_compress_extern;
  if(np_conditional_distribution_cvls_lp_stream_ctx(
       &vector_scale_factor[1], &execution_context, &cv) != 0)
    return DBL_MAX;
  return cv;
}

static void np_beta_regression_bw_route_or_error(
  NPContinuousKernelRoute * const route,
  NPContinuousKernelDerivativeDiagnostics * const diagnostics,
  const char * const where)
{
  route->segment_count = 1;
  route->segment[0].descriptor.family = NP_CKERNEL_FAMILY_BETA;
  route->segment[0].descriptor.legacy_code = NP_CKERNEL_COORDINATE_CODE;
  route->segment[0].descriptor.order = np_beta_bw_order_extern;
  route->segment[0].coordinate_offset = 0;
  route->segment[0].coordinate_count = num_reg_continuous_extern;
  route->segment[0].lower = vector_ckerlb_extern;
  route->segment[0].upper = vector_ckerub_extern;
  if(np_continuous_kernel_route_validate(
       route, num_reg_continuous_extern) != NP_CKERNEL_ROUTE_OK)
    error("%s beta route has an invalid layout", where);
  diagnostics->bad_coordinate = -1;
  diagnostics->bad_observation = -1;
  diagnostics->undefined_count = 0;
  diagnostics->beta_status = NP_BETA_OK;
}


double cv_func_regression_categorical_ls(double *vector_scale_factor){
  double cv = 0.0;
  clock_t start, diff;
  NPContinuousKernelRoute beta_route;
  NPContinuousKernelDerivativeDiagnostics beta_diagnostics;

  if(KERNEL_reg_extern == NP_CKERNEL_COORDINATE_CODE) {
    np_beta_regression_bw_route_or_error(
      &beta_route, &beta_diagnostics, "regression CVLS");
    start = clock();
    cv = np_kernel_estimate_regression_categorical_ls_aic_ctx(
      np_lp_engine_extern, RBWM_CVLS,
      KERNEL_reg_extern, KERNEL_reg_unordered_extern,
      KERNEL_reg_ordered_extern, BANDWIDTH_reg_extern,
      num_obs_train_extern, num_reg_unordered_extern,
      num_reg_ordered_extern, num_reg_continuous_extern,
      matrix_X_unordered_train_extern,
      matrix_X_ordered_train_extern,
      matrix_X_continuous_train_extern, vector_Y_extern,
      &vector_scale_factor[1], num_categories_extern,
      &beta_route, &beta_diagnostics,
      np_regression_bw_categorical_compress_extern);
    diff = clock() - start;
    timing_extern = ((double)diff) / ((double)CLOCKS_PER_SEC);
    return cv;
  }

  if(check_valid_scale_factor_cv(
                                 KERNEL_reg_extern,
                                 KERNEL_reg_unordered_extern,
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
                                 vector_scale_factor) == 1)
    {
      //Rprintf("toasty!\n");
      //for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern; ii++)
      //Rprintf("%3.15g ", vector_scale_factor[ii]);
      //Rprintf("\n");

      return(DBL_MAX);
    }
    start = clock();

    cv = (np_kernel_estimate_regression_categorical_ls_aic(
                                                            np_lp_engine_extern,
                                                            RBWM_CVLS,
                                                            KERNEL_reg_extern,
                                                            KERNEL_reg_unordered_extern,
                                                            KERNEL_reg_ordered_extern,
                                                            BANDWIDTH_reg_extern,
                                                            num_obs_train_extern,
                                                            num_reg_unordered_extern,
                                                            num_reg_ordered_extern,
                                                            num_reg_continuous_extern,
                                                            matrix_X_unordered_train_extern,
                                                            matrix_X_ordered_train_extern,
                                                            matrix_X_continuous_train_extern,
                                                            vector_Y_extern,
                                                            &vector_scale_factor[1],
                                                            num_categories_extern));
    diff = clock() - start;
    timing_extern = ((double)diff)/((double)CLOCKS_PER_SEC);

    return(cv);

}

double cv_func_regression_categorical_ks(double *vector_scale_factor){
  double cv = 0.0;
  clock_t start, diff;

  if(check_valid_scale_factor_cv(
                                 KERNEL_reg_extern,
                                 KERNEL_reg_unordered_extern,
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
                                 vector_scale_factor) == 1)
    {
      return(DBL_MAX);
    }

    start = clock();

    cv = (np_kernel_estimate_regression_categorical_ls_aic(
                                                            np_lp_engine_extern,
                                                            RBWM_CVKS,
                                                            KERNEL_reg_extern,
                                                            KERNEL_reg_unordered_extern,
                                                            KERNEL_reg_ordered_extern,
                                                            BANDWIDTH_reg_extern,
                                                            num_obs_train_extern,
                                                            num_reg_unordered_extern,
                                                            num_reg_ordered_extern,
                                                            num_reg_continuous_extern,
                                                            matrix_X_unordered_train_extern,
                                                            matrix_X_ordered_train_extern,
                                                            matrix_X_continuous_train_extern,
                                                            vector_Y_extern,
                                                            &vector_scale_factor[1],
                                                            num_categories_extern));
    diff = clock() - start;
    timing_extern = ((double)diff)/((double)CLOCKS_PER_SEC);

    return(cv);

}

double cv_func_lsqregression_categorical_check(double *vector_scale_factor){
  double cv = 0.0;
  clock_t start, diff;
  const int nvar =
    num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
  const double delta = vector_scale_factor[nvar + 1];
  double zdelta;
  int i;

  if((delta <= np_lsq_delta_lower_extern) ||
     (delta >= np_lsq_delta_upper_extern) ||
     (!R_FINITE(delta)))
    return(DBL_MAX);

  if(check_valid_scale_factor_cv(
                                 KERNEL_reg_extern,
                                 KERNEL_reg_unordered_extern,
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
                                 vector_scale_factor) == 1)
    return(DBL_MAX);

  if((vector_lsq_scale_extern == NULL) ||
     (vector_lsq_loss_extern == NULL) ||
     (vector_lsq_q_extern == NULL))
    return(DBL_MAX);

  zdelta = qnorm(delta, 0.0, 1.0, 1, 0);
  if(!R_FINITE(zdelta))
    return(DBL_MAX);

  for(i = 0; i < num_obs_train_extern; i++)
    vector_lsq_q_extern[i] =
      vector_lsq_loss_extern[i] + vector_lsq_scale_extern[i]*zdelta;

  start = clock();
  cv = (np_kernel_estimate_regression_categorical_ls_aic(
                                                          np_lp_engine_extern,
                                                          RBWM_CVCHECK,
                                                          KERNEL_reg_extern,
                                                          KERNEL_reg_unordered_extern,
                                                          KERNEL_reg_ordered_extern,
                                                          BANDWIDTH_reg_extern,
                                                          num_obs_train_extern,
                                                          num_reg_unordered_extern,
                                                          num_reg_ordered_extern,
                                                          num_reg_continuous_extern,
                                                          matrix_X_unordered_train_extern,
                                                          matrix_X_ordered_train_extern,
                                                          matrix_X_continuous_train_extern,
                                                          vector_lsq_q_extern,
                                                          &vector_scale_factor[1],
                                                          num_categories_extern));
  diff = clock() - start;
  timing_extern = ((double)diff)/((double)CLOCKS_PER_SEC);

  return(cv);
}

double np_cv_func_density_categorical_ml(double *vector_scale_factor)
{

/* Numerical recipes wrapper function for likelihood density
                    cross-validation */

/* Declarations */

    double cv = 0.0;
    clock_t start, diff;
    NPContinuousKernelRoute beta_route;
    NPContinuousKernelDerivativeDiagnostics beta_diagnostics;
    const NPContinuousKernelRoute *active_route = NULL;
    NPContinuousKernelDerivativeDiagnostics *active_diagnostics = NULL;

    if(KERNEL_den_extern == NP_CKERNEL_COORDINATE_CODE) {
      beta_route.segment_count = 1;
      beta_route.segment[0].descriptor.family = NP_CKERNEL_FAMILY_BETA;
      beta_route.segment[0].descriptor.legacy_code =
        NP_CKERNEL_COORDINATE_CODE;
      beta_route.segment[0].descriptor.order = np_beta_bw_order_extern;
      beta_route.segment[0].coordinate_offset = 0;
      beta_route.segment[0].coordinate_count = num_reg_continuous_extern;
      beta_route.segment[0].lower = vector_ckerlb_extern;
      beta_route.segment[0].upper = vector_ckerub_extern;
      if(np_continuous_kernel_route_validate(
           &beta_route, num_reg_continuous_extern) != NP_CKERNEL_ROUTE_OK)
        error("density CVML beta route has an invalid layout");
      beta_diagnostics.bad_coordinate = -1;
      beta_diagnostics.bad_observation = -1;
      beta_diagnostics.undefined_count = 0;
      beta_diagnostics.beta_status = NP_BETA_OK;
      active_route = &beta_route;
      active_diagnostics = &beta_diagnostics;
    }

    if(active_route == NULL && check_valid_scale_factor_cv(
        KERNEL_den_extern,
        KERNEL_den_unordered_extern,
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
        vector_scale_factor) == 1)
    {
        return(DBL_MAX);
    }

/* Compute the cross-validation function */
    start = clock();
    
    if(np_kernel_estimate_density_categorical_leave_one_out_cv(KERNEL_den_extern,
        KERNEL_den_unordered_extern,
        KERNEL_den_ordered_extern,
        BANDWIDTH_den_extern,
        num_obs_train_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_reg_continuous_extern,
        matrix_X_unordered_train_extern,
        matrix_X_ordered_train_extern,
        matrix_X_continuous_train_extern,
        &vector_scale_factor[1],
        num_categories_extern,
        active_route,
        active_diagnostics,
        np_density_bw_categorical_compress_extern,
        &cv)==1)
    {
        return(DBL_MAX);
    }

    diff = clock() - start;
    timing_extern = ((double)diff)/((double)CLOCKS_PER_SEC);

    return(cv);

}

double cv_func_con_distribution_categorical_ls(double *vector_scale_factor)
{

/* Numerical recipes wrapper function for likelihood density
                    cross-validation */

/* Declarations */

    double cv = 0.0;
    clock_t start, diff;

    if(np_beta_conditional_bw_active()) {
      start = clock();
      cv = np_beta_conditional_distribution_bw_objective_ls_ctx(
        vector_scale_factor);
      diff = clock() - start;
      timing_extern = ((double)diff) / ((double)CLOCKS_PER_SEC);
      return cv;
    }

    if(check_valid_scale_factor_cv(
        KERNEL_den_extern,
        KERNEL_den_unordered_extern,
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
        vector_scale_factor) == 1) {

      //                        Rprintf("toasty!\n");
      //                  for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
      //                    Rprintf("%3.15g ", vector_scale_factor[ii]);
      //                  Rprintf("\n");

      return(DBL_MAX);
    }

/* Compute the cross-validation function */
    start = clock();

    if(np_kernel_estimate_con_distribution_categorical_leave_one_out_ls_cv(KERNEL_den_extern,
                                                                           KERNEL_den_unordered_extern,
                                                                           KERNEL_den_ordered_extern,
                                                                           KERNEL_reg_extern,
                                                                           KERNEL_reg_unordered_extern,
                                                                           KERNEL_reg_ordered_extern,
                                                                           BANDWIDTH_den_extern,
                                                                           num_obs_train_extern,
                                                                           num_obs_eval_extern,
                                                                           num_var_unordered_extern,
                                                                           num_var_ordered_extern,
                                                                           num_var_continuous_extern,
                                                                           num_reg_unordered_extern,
                                                                           num_reg_ordered_extern,
                                                                           num_reg_continuous_extern,
                                                                           cdfontrain_extern,
                                                                           dbl_memfac_ccdf_extern,
                                                                           matrix_Y_unordered_train_extern,
                                                                           matrix_Y_ordered_train_extern,
                                                                           matrix_Y_continuous_train_extern,
                                                                           matrix_X_unordered_train_extern,
                                                                           matrix_X_ordered_train_extern,
                                                                           matrix_X_continuous_train_extern,
                                                                           matrix_XY_unordered_train_extern, 
                                                                           matrix_XY_ordered_train_extern, 
                                                                           matrix_XY_continuous_train_extern,
                                                                           matrix_Y_unordered_eval_extern,
                                                                           matrix_Y_ordered_eval_extern,
                                                                           matrix_Y_continuous_eval_extern,
                                                                           &vector_scale_factor[1],
                                                                           num_categories_extern,
                                                                           matrix_categorical_vals_extern,
                                                                           &cv)==1)
      {
        //                        Rprintf("toaster!\n");
        //                        for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
        //                          Rprintf("%3.15g ", vector_scale_factor[ii]);
        //                        Rprintf("\n");

        return(DBL_MAX);
      }
    diff = clock() - start;
    timing_extern = ((double)diff)/((double)CLOCKS_PER_SEC);


    //        for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
    //          Rprintf("%3.15g ", vector_scale_factor[ii]);
    //                Rprintf("%3.15g ", cv);
    //                Rprintf("\n");

    return(cv);

}

double np_cv_func_con_density_categorical_ml(double *vector_scale_factor){

/* Numerical recipes wrapper function for likelihood density
                    cross-validation */

/* Declarations */

    double cv = 0.0;
    clock_t start, diff;

    if(np_beta_conditional_bw_active()) {
      start = clock();
      cv = np_beta_conditional_density_bw_objective_ctx(
        vector_scale_factor);
      diff = clock() - start;
      timing_extern = ((double)diff) / ((double)CLOCKS_PER_SEC);
      return cv;
    }

    if(check_valid_scale_factor_cv(
        KERNEL_den_extern,
        KERNEL_den_unordered_extern,
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
        vector_scale_factor) == 1) {

      //                  Rprintf("toasty!\n");
      //            for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
      //              Rprintf("%3.15g ", vector_scale_factor[ii]);
      //            Rprintf("\n");

      return(DBL_MAX);
    }
/* Compute the cross-validation function */
    start = clock();

    if(np_kernel_estimate_con_density_categorical_leave_one_out_cv(KERNEL_den_extern,
        KERNEL_den_unordered_extern,
        KERNEL_den_ordered_extern,
				KERNEL_reg_extern,
        KERNEL_reg_unordered_extern,
        KERNEL_reg_ordered_extern,
        BANDWIDTH_den_extern,
        num_obs_train_extern,
        num_var_unordered_extern,
        num_var_ordered_extern,
        num_var_continuous_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_reg_continuous_extern,
        matrix_Y_unordered_train_extern,
        matrix_Y_ordered_train_extern,
        matrix_Y_continuous_train_extern,
        matrix_X_unordered_train_extern,
        matrix_X_ordered_train_extern,
        matrix_X_continuous_train_extern,
        matrix_XY_unordered_train_extern,
        matrix_XY_ordered_train_extern,
        matrix_XY_continuous_train_extern,
        &vector_scale_factor[1],
        num_categories_extern,
        &cv)==1)
    {
      //                  Rprintf("toaster!\n");
      //                  for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
      //                    Rprintf("%3.15g ", vector_scale_factor[ii]);
      //                  Rprintf("\n");

        return(DBL_MAX);
    }
    diff = clock() - start;
    timing_extern = ((double)diff)/((double)CLOCKS_PER_SEC);


    //    for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
    //      Rprintf("%3.15g ", vector_scale_factor[ii]);
    //            Rprintf("%3.15g ", cv);
    //            Rprintf("\n");


    return(cv);

}

double np_cv_func_con_density_categorical_ls_npksum(double *vector_scale_factor){

/* Numerical recipes wrapper function for least squares conditional density
                    cross-validation */

/* Declarations */

  double cv = 0.0;
  clock_t start, diff;

  if(np_beta_conditional_bw_active()) {
    start = clock();
    cv = np_beta_conditional_density_bw_objective_ls_ctx(
      vector_scale_factor);
    diff = clock() - start;
    timing_extern = ((double)diff) / ((double)CLOCKS_PER_SEC);
    return cv;
  }

  if(check_valid_scale_factor_cv(KERNEL_den_extern,
                                 KERNEL_den_unordered_extern,
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
                                 vector_scale_factor) == 1) {
    //        Rprintf("toasty\n");
    //        for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
    //          Rprintf("%3.15g ", vector_scale_factor[ii]);
    //        Rprintf("\n");

    return(DBL_MAX);
  }
  /* Compute the cross-validation function */
    start = clock();

    if(np_kernel_estimate_con_density_categorical_leave_one_out_ls_cv(KERNEL_den_extern,
                                                                      KERNEL_den_unordered_extern,
                                                                      KERNEL_den_ordered_extern,
                                                                      KERNEL_reg_extern,
                                                                      KERNEL_reg_unordered_extern,
                                                                      KERNEL_reg_ordered_extern,
                                                                      BANDWIDTH_den_extern,
                                                                      num_obs_train_extern,
                                                                      num_var_unordered_extern,
                                                                      num_var_ordered_extern,
                                                                      num_var_continuous_extern,
                                                                      num_reg_unordered_extern,
                                                                      num_reg_ordered_extern,
                                                                      num_reg_continuous_extern,
                                                                      dbl_memfac_ccdf_extern,
                                                                      matrix_Y_unordered_train_extern,
                                                                      matrix_Y_ordered_train_extern,
                                                                      matrix_Y_continuous_train_extern,
                                                                      matrix_X_unordered_train_extern,
                                                                      matrix_X_ordered_train_extern,
                                                                      matrix_X_continuous_train_extern,
                                                                      matrix_XY_unordered_train_extern, 
                                                                      matrix_XY_ordered_train_extern, 
                                                                      matrix_XY_continuous_train_extern,
                                                                      &vector_scale_factor[1],
                                                                      num_categories_extern,
                                                                      matrix_categorical_vals_extern,
                                                                      &cv)==1)
      {
        //                Rprintf("toaster!!\n");
        //                for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
        //                  Rprintf("%3.15g ", vector_scale_factor[ii]);
        //                Rprintf("\n");

        return(DBL_MAX);
      }
    diff = clock() - start;
    timing_extern = ((double)diff)/((double)CLOCKS_PER_SEC);

    //        for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
    //          Rprintf("%3.15g ", vector_scale_factor[ii]);
    //          Rprintf("%3.15g ", cv);
    //        Rprintf("\n");

  return(cv);

}

double np_cv_func_density_categorical_ls(double *vector_scale_factor){

/* Numerical recipes wrapper function for likelihood density
                    cross-validation */

/* Declarations */

    double cv = 0.0;
    clock_t start, diff;
    NPContinuousKernelRoute beta_route;
    NPContinuousKernelDerivativeDiagnostics beta_diagnostics;
    const NPContinuousKernelRoute *active_route = NULL;
    NPContinuousKernelDerivativeDiagnostics *active_diagnostics = NULL;

    if(KERNEL_den_extern == NP_CKERNEL_COORDINATE_CODE) {
      beta_route.segment_count = 1;
      beta_route.segment[0].descriptor.family = NP_CKERNEL_FAMILY_BETA;
      beta_route.segment[0].descriptor.legacy_code =
        NP_CKERNEL_COORDINATE_CODE;
      beta_route.segment[0].descriptor.order = np_beta_bw_order_extern;
      beta_route.segment[0].coordinate_offset = 0;
      beta_route.segment[0].coordinate_count = num_reg_continuous_extern;
      beta_route.segment[0].lower = vector_ckerlb_extern;
      beta_route.segment[0].upper = vector_ckerub_extern;
      if(np_continuous_kernel_route_validate(
           &beta_route, num_reg_continuous_extern) != NP_CKERNEL_ROUTE_OK)
        error("density CVLS beta route has an invalid layout");
      beta_diagnostics.bad_coordinate = -1;
      beta_diagnostics.bad_observation = -1;
      beta_diagnostics.undefined_count = 0;
      beta_diagnostics.beta_status = NP_BETA_OK;
      active_route = &beta_route;
      active_diagnostics = &beta_diagnostics;
    }

    if(active_route == NULL && check_valid_scale_factor_cv(
        KERNEL_den_extern,
        KERNEL_den_unordered_extern,
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
        vector_scale_factor) == 1) return(DBL_MAX);

/* Compute the cross-validation function */
    start = clock();

    if(np_fixed_gaussian_density_cvls_pair_dispatch_try(
        KERNEL_den_extern,
        BANDWIDTH_den_extern,
        num_obs_train_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_reg_continuous_extern,
        matrix_X_continuous_train_extern,
        &vector_scale_factor[1],
        &cv))
    {
        diff = clock() - start;
        timing_extern = ((double)diff)/((double)CLOCKS_PER_SEC);
        return(cv);
    }

    if(np_kernel_estimate_density_categorical_convolution_cv(KERNEL_den_extern,
        KERNEL_den_unordered_extern,
        KERNEL_den_ordered_extern,
        BANDWIDTH_den_extern,
        num_obs_train_extern,
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_reg_continuous_extern,
        matrix_X_unordered_train_extern,
        matrix_X_ordered_train_extern,
        matrix_X_continuous_train_extern,
        &vector_scale_factor[1],
        num_categories_extern,
        matrix_categorical_vals_extern,
        active_route,
        active_diagnostics,
        np_density_bw_categorical_compress_extern,
        &cv)==1)
    {
        return(DBL_MAX);
    }

    diff = clock() - start;
    timing_extern = ((double)diff)/((double)CLOCKS_PER_SEC);

    return(cv);

}

double cv_func_distribution_categorical_ls(double *vector_scale_factor)
{

/* Numerical recipes wrapper function for likelihood density
                    cross-validation */

/* Declarations */

    double cv = 0.0;
    clock_t start, diff;
    NPContinuousKernelRoute beta_route;
    NPContinuousKernelDerivativeDiagnostics beta_diagnostics;
    const NPContinuousKernelRoute *active_route = NULL;
    NPContinuousKernelDerivativeDiagnostics *active_diagnostics = NULL;

    if(KERNEL_den_extern == NP_CKERNEL_COORDINATE_CODE) {
      beta_route.segment_count = 1;
      beta_route.segment[0].descriptor.family = NP_CKERNEL_FAMILY_BETA;
      beta_route.segment[0].descriptor.legacy_code =
        NP_CKERNEL_COORDINATE_CODE;
      beta_route.segment[0].descriptor.order = np_beta_bw_order_extern;
      beta_route.segment[0].coordinate_offset = 0;
      beta_route.segment[0].coordinate_count = num_reg_continuous_extern;
      beta_route.segment[0].lower = vector_ckerlb_extern;
      beta_route.segment[0].upper = vector_ckerub_extern;
      if(np_continuous_kernel_route_validate(
           &beta_route, num_reg_continuous_extern) != NP_CKERNEL_ROUTE_OK)
        error("distribution CVLS beta route has an invalid layout");
      beta_diagnostics.bad_coordinate = -1;
      beta_diagnostics.bad_observation = -1;
      beta_diagnostics.undefined_count = 0;
      beta_diagnostics.beta_status = NP_BETA_OK;
      active_route = &beta_route;
      active_diagnostics = &beta_diagnostics;
    }

    if(active_route == NULL && check_valid_scale_factor_cv(
        KERNEL_den_extern,
        KERNEL_den_unordered_extern,
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
        vector_scale_factor) == 1) return(DBL_MAX);

/* Compute the cross-validation function */
    start = clock();
    if(np_kernel_estimate_distribution_ls_cv(KERNEL_den_extern,
                                             KERNEL_den_unordered_extern,
                                             KERNEL_den_ordered_extern,
                                             BANDWIDTH_den_extern,
                                             num_obs_train_extern,
                                             num_obs_eval_extern,
                                             num_reg_unordered_extern,
                                             num_reg_ordered_extern,
                                             num_reg_continuous_extern,
                                             cdfontrain_extern,
                                             dbl_memfac_dls_extern,
                                             matrix_X_unordered_train_extern,
                                             matrix_X_ordered_train_extern,
                                             matrix_X_continuous_train_extern,
                                             matrix_X_unordered_eval_extern,
                                             matrix_X_ordered_eval_extern,
                                             matrix_X_continuous_eval_extern,
                                             &vector_scale_factor[1],
                                             num_categories_extern,
                                             matrix_categorical_vals_extern,
                                             active_route,
                                             active_diagnostics,
                                             np_distribution_bw_categorical_compress_extern,
                                             &cv)==1)
    {
        return(DBL_MAX);
    }

    diff = clock() - start;
    timing_extern = ((double)diff)/((double)CLOCKS_PER_SEC);

    return(cv);

}


double func_con_density_quantile(double *quantile)
{

/* Declarations */

    double func = 0.0;
    double cdf[1];
    double cdf_stderr[1];

    if((quantile[1] < y_min_extern)||(quantile[1] > y_max_extern))
    {
        return(DBL_MAX);
    }

    matrix_Y_continuous_quantile_extern[0][0]=quantile[1];

/* Compute the conditional density at y = quantile */

/* Can we disable MPI temporarily if it is on? */

    kernel_estimate_con_distribution_categorical_no_mpi(
        KERNEL_den_extern,
        KERNEL_den_unordered_extern,
        KERNEL_den_ordered_extern,
				KERNEL_reg_extern,
        KERNEL_reg_unordered_extern,
        KERNEL_reg_ordered_extern,
        BANDWIDTH_den_extern,
        num_obs_train_extern,
        1,                                        /* One evaluation observation */
        0,                                        /* Zero discrete Y */
        0,                                        /* Zero discrete Y */
        1,                                        /* One continuous Y */
        num_reg_unordered_extern,
        num_reg_ordered_extern,
        num_reg_continuous_extern,
        matrix_Y_unordered_train_extern,
        matrix_Y_ordered_train_extern,
        matrix_Y_continuous_train_extern,
        matrix_Y_unordered_quantile_extern, /* Not used */
        matrix_Y_ordered_quantile_extern,   /* Not used */
        matrix_Y_continuous_quantile_extern,
        matrix_X_unordered_train_extern,
        matrix_X_ordered_train_extern,
        matrix_X_continuous_train_extern,
        matrix_X_unordered_quantile_extern,
        matrix_X_ordered_quantile_extern,
        matrix_X_continuous_quantile_extern,
        &vector_scale_factor_extern[1],
        num_categories_extern,
        matrix_categorical_vals_extern,
        cdf,
        cdf_stderr,
        small_extern,
        itmax_extern);

    func = ipow(gamma_extern - cdf[0], 2);


    return(func);

}


double cv_func_regression_categorical_aic_c(double *vector_scale_factor)
{

/* Numerical recipes wrapper function for Hurvich/Simonoff/Tsai JRSS B 1998 */

/* Declarations */
  double cv = 0.0;
  clock_t start, diff;
  NPContinuousKernelRoute beta_route;
  NPContinuousKernelDerivativeDiagnostics beta_diagnostics;

  if(KERNEL_reg_extern == NP_CKERNEL_COORDINATE_CODE) {
    np_beta_regression_bw_route_or_error(
      &beta_route, &beta_diagnostics, "regression CVAIC");
    start = clock();
    cv = np_kernel_estimate_regression_categorical_ls_aic_ctx(
      np_lp_engine_extern, RBWM_CVAIC,
      KERNEL_reg_extern, KERNEL_reg_unordered_extern,
      KERNEL_reg_ordered_extern, BANDWIDTH_reg_extern,
      num_obs_train_extern, num_reg_unordered_extern,
      num_reg_ordered_extern, num_reg_continuous_extern,
      matrix_X_unordered_train_extern,
      matrix_X_ordered_train_extern,
      matrix_X_continuous_train_extern, vector_Y_extern,
      &vector_scale_factor[1], num_categories_extern,
      &beta_route, &beta_diagnostics,
      np_regression_bw_categorical_compress_extern);
    diff = clock() - start;
    timing_extern = ((double)diff) / ((double)CLOCKS_PER_SEC);
    return cv;
  }

  if(check_valid_scale_factor_cv(
                                 KERNEL_reg_extern,
                                 KERNEL_reg_unordered_extern,
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
                                 vector_scale_factor) == 1)
    {
      return(DBL_MAX);
    }

    start = clock();

    cv = (np_kernel_estimate_regression_categorical_ls_aic(np_lp_engine_extern,
                                                            RBWM_CVAIC,
                                                            KERNEL_reg_extern,
                                                            KERNEL_reg_unordered_extern,
                                                            KERNEL_reg_ordered_extern,
                                                            BANDWIDTH_reg_extern,
                                                            num_obs_train_extern,
                                                            num_reg_unordered_extern,
                                                            num_reg_ordered_extern,
                                                            num_reg_continuous_extern,
                                                            matrix_X_unordered_train_extern,
                                                            matrix_X_ordered_train_extern,
                                                            matrix_X_continuous_train_extern,
                                                            vector_Y_extern,
                                                            &vector_scale_factor[1],
                                                            num_categories_extern));
    diff = clock() - start;
    timing_extern = ((double)diff)/((double)CLOCKS_PER_SEC);

    return(cv);
}
