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

#include <R.h>

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
extern int *num_categories_extern_X;
extern int *num_categories_extern_Y;
extern double **matrix_categorical_vals_extern_X;
extern double **matrix_categorical_vals_extern_Y;

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
extern KDT * kdt_extern_Y;

extern double *vector_Y_extern;
extern double *vector_T_extern;
extern double *vector_Y_eval_extern;

/* Quantile - no Y ordered or unordered used, but defined anyways */

extern double **matrix_Y_continuous_quantile_extern;
extern double **matrix_Y_unordered_quantile_extern;
extern double **matrix_Y_ordered_quantile_extern;
extern double **matrix_X_continuous_quantile_extern;
extern double **matrix_X_unordered_quantile_extern;
extern double **matrix_X_ordered_quantile_extern;

extern int int_ll_extern;

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

#define RBWM_CVAIC 0
#define RBWM_CVLS 1

#define LL_LC  0
#define LL_LL  1
#define LL_LP  2

#define BW_FIXED   0
#define BW_GEN_NN  1
#define BW_ADAP_NN 2

extern int int_TREE_Y;
extern int int_TREE_X;

#define NP_COND_CV_DENS_ML 1
#define NP_COND_CV_DENS_LS 2
#define NP_COND_CV_DIST_LS 3

static double np_lp_conditional_cv_objective(const int mode,
                                             double *vector_scale_factor){
  const int ntrain = num_obs_train_extern;
  const int neval = (mode == NP_COND_CV_DIST_LS) ? num_obs_eval_extern : num_obs_train_extern;
  const int nx = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
  const int ny = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  const int y_operator_mode = (mode == NP_COND_CV_DIST_LS) ? OP_INTEGRAL : OP_NORMAL;
  double cv = 0.0;
  double RS = 0.0, MSE = 0.0, MAE = 0.0, MAPE = 0.0, CORR = 0.0, SIGN = 0.0;
  int status = 0;
  int i, j, l;
  int saved_tree_x = int_TREE_X;

  int *kernel_cy = NULL, *kernel_uy = NULL, *kernel_oy = NULL, *operator_y = NULL;
  double *vsf_x = NULL, *vsf_y = NULL, *ykw = NULL;
  double *pred_one = NULL, *pred_se_one = NULL, *y_eval_one = NULL;
  double **yuno_eval_one = NULL, **yord_eval_one = NULL, **ycon_eval_one = NULL;
  double **xuno_eval_one = NULL, **xord_eval_one = NULL, **xcon_eval_one = NULL;

  if((ntrain <= 1) || (neval <= 0))
    return DBL_MAX;

  if(check_valid_scale_factor_cv(KERNEL_den_extern,
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
                                 vector_scale_factor) == 1)
    return DBL_MAX;

  vsf_x = alloc_vecd(MAX(1, nx));
  vsf_y = alloc_vecd(MAX(1, ny));
  ykw = alloc_vecd(MAX(1, ntrain));
  pred_one = alloc_vecd(1);
  pred_se_one = alloc_vecd(1);
  y_eval_one = alloc_vecd(1);

  if((num_var_unordered_extern > 0) && ((yuno_eval_one = alloc_matd(1, num_var_unordered_extern)) == NULL))
    goto fail_lp_cond_cv;
  if((num_var_ordered_extern > 0) && ((yord_eval_one = alloc_matd(1, num_var_ordered_extern)) == NULL))
    goto fail_lp_cond_cv;
  if((num_var_continuous_extern > 0) && ((ycon_eval_one = alloc_matd(1, num_var_continuous_extern)) == NULL))
    goto fail_lp_cond_cv;
  if((num_reg_unordered_extern > 0) && ((xuno_eval_one = alloc_matd(1, num_reg_unordered_extern)) == NULL))
    goto fail_lp_cond_cv;
  if((num_reg_ordered_extern > 0) && ((xord_eval_one = alloc_matd(1, num_reg_ordered_extern)) == NULL))
    goto fail_lp_cond_cv;
  if((num_reg_continuous_extern > 0) && ((xcon_eval_one = alloc_matd(1, num_reg_continuous_extern)) == NULL))
    goto fail_lp_cond_cv;

  kernel_cy = (int *)malloc((size_t)MAX(1, num_var_continuous_extern)*sizeof(int));
  kernel_uy = (int *)malloc((size_t)MAX(1, num_var_unordered_extern)*sizeof(int));
  kernel_oy = (int *)malloc((size_t)MAX(1, num_var_ordered_extern)*sizeof(int));
  operator_y = (int *)malloc((size_t)MAX(1, ny)*sizeof(int));

  if((vsf_x == NULL) || (vsf_y == NULL) || (ykw == NULL) ||
     (pred_one == NULL) || (pred_se_one == NULL) || (y_eval_one == NULL) ||
     (kernel_cy == NULL) || (kernel_uy == NULL) || (kernel_oy == NULL) || (operator_y == NULL))
    goto fail_lp_cond_cv;

  for(i = 0; i < ny; i++)
    operator_y[i] = y_operator_mode;
  for(i = 0; i < num_var_continuous_extern; i++)
    kernel_cy[i] = KERNEL_den_extern;
  for(i = 0; i < num_var_unordered_extern; i++)
    kernel_uy[i] = KERNEL_den_unordered_extern;
  for(i = 0; i < num_var_ordered_extern; i++)
    kernel_oy[i] = KERNEL_den_ordered_extern;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        vsf_x,
                        vsf_y,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  int_TREE_X = NP_TREE_FALSE;

  for(j = 0; j < neval; j++){
    if(mode == NP_COND_CV_DIST_LS){
      for(l = 0; l < num_var_unordered_extern; l++)
        yuno_eval_one[l][0] = matrix_Y_unordered_eval_extern[l][j];
      for(l = 0; l < num_var_ordered_extern; l++)
        yord_eval_one[l][0] = matrix_Y_ordered_eval_extern[l][j];
      for(l = 0; l < num_var_continuous_extern; l++)
        ycon_eval_one[l][0] = matrix_Y_continuous_eval_extern[l][j];
    } else {
      for(l = 0; l < num_var_unordered_extern; l++)
        yuno_eval_one[l][0] = matrix_Y_unordered_train_extern[l][j];
      for(l = 0; l < num_var_ordered_extern; l++)
        yord_eval_one[l][0] = matrix_Y_ordered_train_extern[l][j];
      for(l = 0; l < num_var_continuous_extern; l++)
        ycon_eval_one[l][0] = matrix_Y_continuous_train_extern[l][j];
    }

    status = kernel_weighted_sum_np(kernel_cy,
                                    kernel_uy,
                                    kernel_oy,
                                    BANDWIDTH_den_extern,
                                    ntrain,
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
                                    0,
                                    0,
                                    0,
                                    NP_TREE_FALSE,
                                    0,
                                    NULL,
                                    NULL, NULL, NULL,
                                    matrix_Y_unordered_train_extern,
                                    matrix_Y_ordered_train_extern,
                                    matrix_Y_continuous_train_extern,
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
                                    num_categories_extern_Y,
                                    matrix_categorical_vals_extern_Y,
                                    NULL,
                                    NULL,
                                    NULL,
                                    ykw);
    if(status != 0)
      goto fail_lp_cond_cv;

    for(i = 0; i < ntrain; i++){
      for(l = 0; l < num_reg_unordered_extern; l++)
        xuno_eval_one[l][0] = matrix_X_unordered_train_extern[l][i];
      for(l = 0; l < num_reg_ordered_extern; l++)
        xord_eval_one[l][0] = matrix_X_ordered_train_extern[l][i];
      for(l = 0; l < num_reg_continuous_extern; l++)
        xcon_eval_one[l][0] = matrix_X_continuous_train_extern[l][i];

      y_eval_one[0] = ykw[i];

      status = kernel_estimate_regression_categorical_tree_np(int_ll_extern,
                                                              KERNEL_reg_extern,
                                                              KERNEL_reg_unordered_extern,
                                                              KERNEL_reg_ordered_extern,
                                                              BANDWIDTH_den_extern,
                                                              ntrain,
                                                              1,
                                                              num_reg_unordered_extern,
                                                              num_reg_ordered_extern,
                                                              num_reg_continuous_extern,
                                                              matrix_X_unordered_train_extern,
                                                              matrix_X_ordered_train_extern,
                                                              matrix_X_continuous_train_extern,
                                                              xuno_eval_one,
                                                              xord_eval_one,
                                                              xcon_eval_one,
                                                              ykw,
                                                              y_eval_one,
                                                              vsf_x,
                                                              num_categories_extern_X,
                                                              matrix_categorical_vals_extern_X,
                                                              pred_one,
                                                              NULL,
                                                              pred_se_one,
                                                              NULL,
                                                              &RS,
                                                              &MSE,
                                                              &MAE,
                                                              &MAPE,
                                                              &CORR,
                                                              &SIGN,
                                                              0,
                                                              1,
                                                              i);
      if(status != 0)
        goto fail_lp_cond_cv;

      if((mode == NP_COND_CV_DENS_ML) && (i == j)){
        const double p = pred_one[0];
        cv -= (p < DBL_MIN) ? log(DBL_MIN) : log(p);
      } else if(mode == NP_COND_CV_DENS_LS){
        const double d = ykw[i] - pred_one[0];
        cv += d*d;
      } else if(mode == NP_COND_CV_DIST_LS){
        double indy = 1.0;
        for(l = 0; l < num_var_ordered_extern; l++)
          indy *= (matrix_Y_ordered_train_extern[l][i] <= matrix_Y_ordered_eval_extern[l][j]);
        for(l = 0; l < num_var_continuous_extern; l++)
          indy *= (matrix_Y_continuous_train_extern[l][i] <= matrix_Y_continuous_eval_extern[l][j]);
        {
          const double d = indy - pred_one[0];
          cv += d*d;
        }
      }
    }
  }

  if((mode == NP_COND_CV_DENS_LS) || (mode == NP_COND_CV_DIST_LS))
    cv /= ((double)ntrain * (double)neval);

  goto cleanup_lp_cond_cv;

fail_lp_cond_cv:
  cv = DBL_MAX;

cleanup_lp_cond_cv:
  int_TREE_X = saved_tree_x;
  if(yuno_eval_one != NULL) free_mat(yuno_eval_one, num_var_unordered_extern);
  if(yord_eval_one != NULL) free_mat(yord_eval_one, num_var_ordered_extern);
  if(ycon_eval_one != NULL) free_mat(ycon_eval_one, num_var_continuous_extern);
  if(xuno_eval_one != NULL) free_mat(xuno_eval_one, num_reg_unordered_extern);
  if(xord_eval_one != NULL) free_mat(xord_eval_one, num_reg_ordered_extern);
  if(xcon_eval_one != NULL) free_mat(xcon_eval_one, num_reg_continuous_extern);
  if(kernel_cy != NULL) free(kernel_cy);
  if(kernel_uy != NULL) free(kernel_uy);
  if(kernel_oy != NULL) free(kernel_oy);
  if(operator_y != NULL) free(operator_y);
  if(vsf_x != NULL) free(vsf_x);
  if(vsf_y != NULL) free(vsf_y);
  if(ykw != NULL) free(ykw);
  if(pred_one != NULL) free(pred_one);
  if(pred_se_one != NULL) free(pred_se_one);
  if(y_eval_one != NULL) free(y_eval_one);
  return cv;
}


double cv_func_regression_categorical_ls(double *vector_scale_factor){
  double cv = 0.0;
  clock_t start, diff;

  if(check_valid_scale_factor_cv(
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
                                                            int_ll_extern,
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

double cv_func_regression_categorical_ls_nn(double *vector_scale_factor)
{

/* Numerical recipes wrapper function for least squares regression
                    cross-validation */

/* Declarations */

    double cv = 0.0;

    int i;

    double *mean;

    double *py;
    double *pm;

#ifdef MPI2
    int stride;
#endif

/* Allocate memory for objects */

#ifndef MPI2
    mean = alloc_vecd(num_obs_train_extern);
#endif

#ifdef MPI2

    stride = (int)ceil((double) num_obs_train_extern / (double) iNum_Processors);
    if(stride < 1) stride = 1;
    mean = alloc_vecd(stride*iNum_Processors);
#endif

/* Compute the cross-validation function */

    if(kernel_estimate_regression_categorical_leave_one_out(
        int_ll_extern,
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
        num_categories_extern,
        mean)==1)
    {
        free(mean);
        return(DBL_MAX);
    }

    py = &vector_Y_extern[0];
    pm = &mean[0];

    for(i=0;i<num_obs_train_extern;i++)
    {
        cv += ipow((*py++ - *pm++),2);
    }

    cv /= (double) num_obs_train_extern;


    free(mean);

    return(cv);

}



double cv_func_density_categorical_ml(double *vector_scale_factor)
{

/* Numerical recipes wrapper function for likelihood density
                    cross-validation */

/* Declarations */

    double cv = 0.0;

    if(check_valid_scale_factor_cv(
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
        vector_scale_factor) == 1)
    {
        return(DBL_MAX);
    }

/* Compute the cross-validation function */

    if(kernel_estimate_density_categorical_leave_one_out_cv(KERNEL_den_extern,
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
        &cv)==1)
    {
        return(DBL_MAX);
    }
    

    return(cv);

}

double np_cv_func_density_categorical_ml(double *vector_scale_factor)
{

/* Numerical recipes wrapper function for likelihood density
                    cross-validation */

/* Declarations */

    double cv = 0.0;
    clock_t start, diff;

    if(check_valid_scale_factor_cv(
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

    if((int_ll_extern != LL_LC) && (num_reg_continuous_extern > 0))
      return np_lp_conditional_cv_objective(NP_COND_CV_DIST_LS, vector_scale_factor);

    if(check_valid_scale_factor_cv(
        KERNEL_den_extern,
        KERNEL_reg_unordered_extern, /* Only for conditioning vars in conditional den */
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

double cv_func_con_density_categorical_ml(double *vector_scale_factor)
{

/* Numerical recipes wrapper function for likelihood density
                    cross-validation */

/* Declarations */

    double cv = 0.0;

    if((int_ll_extern != LL_LC) && (num_reg_continuous_extern > 0))
      return np_lp_conditional_cv_objective(NP_COND_CV_DENS_ML, vector_scale_factor);

    if(check_valid_scale_factor_cv(
        KERNEL_den_extern,
        KERNEL_reg_unordered_extern, /* Only for conditioning vars in conditional den */
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
        vector_scale_factor) == 1) return(DBL_MAX);

/* Compute the cross-validation function */

    if(kernel_estimate_con_density_categorical_leave_one_out_cv(KERNEL_den_extern,
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
        &vector_scale_factor[1],
        num_categories_extern,
        &cv)==1)
    {
        return(DBL_MAX);
    }


    return(cv);

}

double np_cv_func_con_density_categorical_ml(double *vector_scale_factor){

/* Numerical recipes wrapper function for likelihood density
                    cross-validation */

/* Declarations */

    double cv = 0.0;
    clock_t start, diff;

    if((int_ll_extern != LL_LC) && (num_reg_continuous_extern > 0))
      return np_lp_conditional_cv_objective(NP_COND_CV_DENS_ML, vector_scale_factor);

    if(check_valid_scale_factor_cv(
        KERNEL_den_extern,
        KERNEL_reg_unordered_extern, /* Only for conditioning vars in conditional den */
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

double np_cv_func_con_density_categorical_ls(double *vector_scale_factor){

/* Numerical recipes wrapper function for least squares conditional density
                    cross-validation */

/* Declarations */

  double cv = 0.0;

  if((int_ll_extern != LL_LC) && (num_reg_continuous_extern > 0))
    return np_lp_conditional_cv_objective(NP_COND_CV_DENS_LS, vector_scale_factor);

  if(check_valid_scale_factor_cv(KERNEL_den_extern,
                                 KERNEL_reg_unordered_extern,  /* Only for conditioning vars in conditional den */
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

  if(np_kernel_estimate_con_density_categorical_convolution_cv(KERNEL_den_extern,
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
                                                               &vector_scale_factor[1],
                                                               num_categories_extern,
                                                               matrix_categorical_vals_extern,
                                                               &cv)==1) {
    //                Rprintf("toaster!!\n");
    //                for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
    //                  Rprintf("%3.15g ", vector_scale_factor[ii]);
    //                Rprintf("\n");

    return(DBL_MAX);
  }

  //        for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
  //          Rprintf("%3.15g ", vector_scale_factor[ii]);
  //          Rprintf("%3.15g ", cv);
  //        Rprintf("\n");

  return(cv);

}

double np_cv_func_con_density_categorical_ls_npksum(double *vector_scale_factor){

/* Numerical recipes wrapper function for least squares conditional density
                    cross-validation */

/* Declarations */

  double cv = 0.0;
  clock_t start, diff;

  if((int_ll_extern != LL_LC) && (num_reg_continuous_extern > 0))
    return np_lp_conditional_cv_objective(NP_COND_CV_DENS_LS, vector_scale_factor);

  if(check_valid_scale_factor_cv(KERNEL_den_extern,
                                 KERNEL_reg_unordered_extern,  /* Only for conditioning vars in conditional den */
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

double cv_func_con_density_categorical_ls(double *vector_scale_factor)
{

/* Numerical recipes wrapper function for least squares conditional density
                    cross-validation */

/* Declarations */

    double cv = 0.0;

    if((int_ll_extern != LL_LC) && (num_reg_continuous_extern > 0))
      return np_lp_conditional_cv_objective(NP_COND_CV_DENS_LS, vector_scale_factor);

    if(check_valid_scale_factor_cv(
        KERNEL_den_extern,
        KERNEL_reg_unordered_extern, /* Only for conditioning vars in conditional den */
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
        vector_scale_factor) == 1){

      //      Rprintf("toasty!!\n");
      //      for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
      //        Rprintf("%3.15g ", vector_scale_factor[ii]);
      //      Rprintf("\n");

      return(DBL_MAX);
    }
/* Compute the cross-validation function */

    if(kernel_estimate_con_density_categorical_convolution_cv(KERNEL_den_extern,
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
        &vector_scale_factor[1],
        num_categories_extern,
        matrix_categorical_vals_extern,
        &cv)==1)
    {
      //      Rprintf("toaster!!\n");
      //      for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
      //        Rprintf("%3.15g ", vector_scale_factor[ii]);
      //      Rprintf("\n");

      return(DBL_MAX);
    }

    //    for(int ii = 1; ii <= num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern + num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern; ii++)
    //      Rprintf("%3.15g ", vector_scale_factor[ii]);
    //      Rprintf("%3.15g ", cv);
    //    Rprintf("\n");

    return(cv);

}

/* Feb 7 2010 */

double cv_func_con_distribution_categorical_ccdf(double *vector_scale_factor)
{

/* Numerical recipes wrapper function for conditional distribution
function */

/* Declarations */

    double cv = 0.0;

    if(check_valid_scale_factor_cv(
        KERNEL_den_extern,
        KERNEL_reg_unordered_extern, /* Only for conditioning vars in conditional den */
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
        vector_scale_factor) == 1) return(DBL_MAX);

    if(kernel_estimate_con_distribution_categorical_leave_one_out_ccdf(KERNEL_den_extern,
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
        &vector_scale_factor[1],
        num_categories_extern,
        matrix_categorical_vals_extern,
        &cv,
        small_extern,
        itmax_extern)==1)
    {
        return(DBL_MAX);
    }


    return(cv);

}

double cv_func_density_categorical_ls(double *vector_scale_factor)
{

/* Numerical recipes wrapper function for likelihood density
                    cross-validation */

/* Declarations */

    double cv = 0.0;

    if(check_valid_scale_factor_cv(
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
        vector_scale_factor) == 1) return(DBL_MAX);

/* Compute the cross-validation function */

    if(kernel_estimate_density_categorical_convolution_cv(KERNEL_den_extern,
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
        &cv)==1)
    {
        return(DBL_MAX);
    }


    return(cv);

}

double np_cv_func_density_categorical_ls(double *vector_scale_factor){

/* Numerical recipes wrapper function for likelihood density
                    cross-validation */

/* Declarations */

    double cv = 0.0;
    clock_t start, diff;

    if(check_valid_scale_factor_cv(
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
        vector_scale_factor) == 1) return(DBL_MAX);

/* Compute the cross-validation function */
    start = clock();

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

    if(check_valid_scale_factor_cv(
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

  if(check_valid_scale_factor_cv(
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
                                 vector_scale_factor) == 1)
    {
      return(DBL_MAX);
    }

    start = clock();

    cv = (np_kernel_estimate_regression_categorical_ls_aic(int_ll_extern,
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
