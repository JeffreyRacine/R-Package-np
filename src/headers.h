#ifndef _NP_
#define _NP_

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <stdint.h>

#include "tree.h"

#define MAX_OBS INT_MAX
#define MAX_REG 100                               /* Allows flexibility while trapping error... */
#define VERSION 1.1

#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif


/* SUNOS does not have ANSI difftime(), so create a workaround */

#ifdef SUNOS
int difftime(time_t t1, time_t t0);
int difftime(time_t t1, time_t t0)
{
  t0 = t0 - time(NULL);
  t1 = t1 - time(NULL);
  return (t0 <= t1 ? (int) (t1-t0) : -(int)(t0-t1));
}
#endif

/* Function prototypes */

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
/* No zero divide allowing for negative or positive values */
#define NZD(a) ((a) < 0 ? (MIN(-DBL_EPSILON,a)) : (MAX(DBL_EPSILON,a)))

double **alloc_matd(int nrows, int ncols);
double **alloc_tmatd(int nrows, int ncols);

double *alloc_vecd(int nobs);
double *vector(int nl,int nh);
double brent(double ax, double bx, double cx, double (*f)(double), double tol, double small, int itmax, double *xmin);
double cdf_kernel(int KERNEL, double z);
double cdf_kernel_unordered(int KERNEL, double x, double y, double lambda, int c, double *categorical_vals);
double cdf_kernel_ordered(int KERNEL, double x, double y, double lambda, int c, double *categorical_vals);
double cv_func_binary_choice(double *scale_factor);
double cv_func_con_density_categorical_ls(double *vector_scale_factor);
double cv_func_con_density_categorical_ml(double *vector_scale_factor);
double cv_func_density_categorical_ls(double *vector_scale_factor);
double cv_func_density_categorical_ml(double *vector_scale_factor);

double np_cv_func_density_categorical_ml(double *vector_scale_factor);
double np_cv_func_density_categorical_ls(double *vector_scale_factor);

double cv_func_con_distribution_categorical_ls(double *vector_scale_factor);

double cv_func_distribution_categorical_ls(double *vector_scale_factor);
/* double cv_func_np_density_categorical_ml(double *vector_scale_factor); */
double cv_func_regression_categorical_ls(double *vector_scale_factor);

/*double erfun(double x, double small, int itmax);*/
double erfun(double x);
double f1dim(double x);
double fMSE(int iNum_Obs, double *fvector_Y, double *fkernel_fit);
double fMAE(int iNum_Obs, double *fvector_Y, double *fkernel_fit);
double fMAPE(int iNum_Obs, double *fvector_Y, double *fkernel_fit);
double fSIGN(int iNum_Obs, double *fvector_Y, double *fkernel_fit);
double fCORR(int iNum_Obs, double *fvector_Y, double *fkernel_fit);
double fGoodness_of_Fit(int iNum_Obs, double *fvector_Y, double *fkernel_fit);

double kernel(int KERNEL, double z);
double kernel_unordered(int KERNEL, double x, double y, double lambda, int c);
double kernel_unordered_ratio(int KERNEL, double x, double y, double lambda, int c);
double kernel_ordered(int KERNEL, double x, double y, double lambda);
double kernel_unordered_convolution(int KERNEL, double x, double y, double lambda, int c, double *c_vals);
double kernel_ordered_convolution(int KERNEL, double x, double y, double lambda, int c, double *c_vals);
double kernel_convol(int KERNEL, int BANDWIDTH, double z, double h1, double h2);
double kernel_deriv(int KERNEL, double z);
double meand(int n, double *vector);

double ran3(int *idum);
double gasdev(int *idum);
double chidev(int *idum, int df);

double standerrd(int n, double *vector);

int compute_nn_distance(int num_obs, int suppress_parallel, double *vector_data, int int_k_nn, double *nn_distance);
int compute_nn_distance_train_eval(int num_obs_train, int num_obs_eval, int suppress_parallel, double *vector_data_train, double *vector_data_eval, int int_k_nn, double *nn_distance);

int determine_categorical_vals(int num_obs, int num_var_unordered, int num_var_ordered, int num_reg_unordered, int num_reg_ordered, double **matrix_Y_unordered, double **matrix_Y_ordered, double **matrix_X_unordered, double **matrix_X_ordered, int *num_categories, double **matrix_categorical_vals);

int np_fround(double x);

int initialize_kernel_density_asymptotic_constants(int KERNEL, int num_var_continuous, double *INT_KERNEL_P, double *K_INT_KERNEL_P);  
int initialize_kernel_regression_asymptotic_constants(int KERNEL, int num_reg_continuous, double *INT_KERNEL_P, double *K_INT_KERNEL_P, double *INT_KERNEL_PM_HALF, double *DIFF_KER_PPM); 

int initialize_nr_directions(int BANDWIDTH,int num_obs,int num_reg_continuous,int num_reg_unordered,int num_reg_ordered,int num_var_continuous,int num_var_unordered,int num_var_ordered,double * vector_scale_factor,int * num_categories,double **matrix_y,int random,int seed,double lbc_dir,int dfc_dir,double c_dir,double initc_dir,double lbd_dir,double hbd_dir,double d_dir, double initd_dir, double ** matrix_x_continuous,double ** matrix_y_continuous);

int kernel_bandwidth(int KERNEL, int BANDWIDTH, int num_obs_train, int num_obs_eval, int num_var_cont, int num_var_un, int num_var_or, int num_reg_cont, int num_reg_un, int num_reg_or, double *vector_scale_factor, double **matrix_Y_train, double **matrix_Y_eval, double **matrix_X_train, double **matrix_X_eval, double **matrix_bandwidth_Y, double **matrix_bandwidth_X, double *vector_lambda, double **matrix_bandwidth_deriv);
int kernel_bandwidth_mean(int KERNEL, int BANDWIDTH, int num_obs_train, int num_obs_eval, int num_var_cont, int num_var_un, int num_var_or, int num_reg_cont, int num_reg_un, int num_reg_or, int suppress_parallel, double *vector_scale_factor, double **matrix_Y_train, double **matrix_Y_eval, double **matrix_X_train, double **matrix_X_eval, double **matrix_bandwidth_Y, double **matrix_bandwidth_X, double *vector_lambda);

int kernel_estimate_categorical_gradient_ocg_fast(int int_COMPUTE, int *var_index_int, int num_var_test_int, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_reg, int int_ll, int int_ordered_categorical_gradient, int num_obs_train, int num_obs_eval, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double *vector_Y, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double **matrix_bandwidth, double **matrix_bandwidth_deriv, double *lambda, int *num_categories, double **matrix_categorical_vals, double *mean, double **gradient_categorical);

int kernel_estimate_con_density_categorical(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double *pdf, double *pdf_stderr, double *log_likelihood);

int kernel_estimate_con_density_categorical_convolution_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered, double **matrix_Y_ordered, double **matrix_Y_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cv);

int kernel_estimate_con_density_categorical_gradient(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double *pdf, double *pdf_stderr, double **pdf_deriv, double **pdf_deriv_stderr, double *log_likelihood);

int kernel_estimate_con_density_categorical_gradient_categorical(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, int int_ordered_categorical_gradient, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, double **matrix_categorical_vals, int *num_categories, double *pdf, double **pdf_deriv, double **pdf_deriv_stderr);

int kernel_estimate_con_density_categorical_leave_one_out_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered, double **matrix_Y_ordered, double **matrix_Y_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_scale_factor, int *num_categories, double *cv);

int np_kernel_estimate_con_density_categorical_leave_one_out_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered, double **matrix_Y_ordered, double **matrix_Y_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double **matrix_XY_unordered, double **matrix_XY_ordered, double **matrix_XY_continuous, double *vector_scale_factor, int *num_categories, double *cv);

int kernel_estimate_con_distribution_categorical(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cdf, double *cdf_stderr, double small, int itmax);

int kernel_estimate_con_distribution_categorical_leave_one_out(int KERNEL_den,int KERNEL_unordered_den,int KERNEL_ordered_den,int KERNEL_reg,int KERNEL_unordered_reg,int KERNEL_ordered_reg,int BANDWIDTH_den,int num_obs_train,int num_obs_eval,int num_var_unordered,int num_var_ordered,int num_var_continuous,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double **matrix_Y_unordered_train,double **matrix_Y_ordered_train,double **matrix_Y_continuous_train,double **matrix_Y_unordered_eval,double **matrix_Y_ordered_eval,double **matrix_Y_continuous_eval,double **matrix_X_unordered_train,double **matrix_X_ordered_train,double **matrix_X_continuous_train,double **matrix_X_unordered_eval,double **matrix_X_ordered_eval,double **matrix_X_continuous_eval,double *vector_scale_factor,int *num_categories,double **matrix_categorical_vals,double *cdf,double small,int itmax);

int indfunc(double a);
double cv_func_con_distribution_categorical_ccdf(double *vector_scale_factor);

int kernel_estimate_con_distribution_categorical_leave_one_out_ccdf(int KERNEL_den,int KERNEL_unordered_den,int KERNEL_ordered_den,int KERNEL_reg,int KERNEL_unordered_reg,int KERNEL_ordered_reg,int BANDWIDTH_den,int num_obs_train,int num_var_unordered,int num_var_ordered,int num_var_continuous,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double **matrix_Y_unordered_train,double **matrix_Y_ordered_train,double **matrix_Y_continuous_train,double **matrix_X_unordered_train,double **matrix_X_ordered_train,double **matrix_X_continuous_train,double *vector_scale_factor,int *num_categories,double **matrix_categorical_vals,double *cv,double small,int itmax);

int kernel_estimate_con_distribution_categorical_no_mpi(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cdf, double *cdf_stderr, double small, int itmax);

int kernel_estimate_con_distribution_categorical_gradient(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cdf, double *cdf_stderr, double **cdf_deriv, double **cdf_deriv_stderr, double small, int itmax);

int kernel_estimate_con_distribution_categorical_gradient_categorical(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, int int_ordered_categorical_gradient, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, double **matrix_categorical_vals, int *num_categories, double *cdf, double **cdf_deriv, double **cdf_deriv_stderr, double small, int itmax);

int kernel_estimate_density_categorical(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double *pdf, double *pdf_stderr, double *log_likelihood);

void kernel_estimate_dens_dist_categorical_np(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, int dop, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double ** matrix_categorical_vals, double *pdf, double *pdf_stderr, double *log_likelihood);


int kernel_estimate_density_categorical_convolution_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cv);

int np_kernel_estimate_density_categorical_convolution_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cv);

int kernel_estimate_density_categorical_leave_one_out_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_scale_factor, int *num_categories, double *cv);

int np_kernel_estimate_density_categorical_leave_one_out_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_scale_factor, int *num_categories, double *cv);

int kernel_estimate_distribution_categorical(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cdf, double *cdf_stderr, double small, int itmax);

int kernel_estimate_regression_categorical_leave_one_out(int int_ll, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_reg, int num_obs, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_Y, double *vector_scale_factor, int *num_categories, double *mean);
double np_kernel_estimate_regression_categorical_ls_aic(int int_ll, int bwm, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_reg, int num_obs, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_Y, double *vector_scale_factor, int *num_categories);
double cv_func_regression_categorical_ls_nn(double *vector_scale_factor);

int kernel_estimate_regression_categorical_no_stderr(int int_compute_gradient, int int_ll, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_reg, int int_WEIGHTS, int *var_index_int, int num_var_test_int, double **matrix_weights_K, double ***matrix_weights_K_deriv, int num_obs_train, int num_obs_eval, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double **matrix_bandwidth, double **matrix_bandwidth_deriv, double *vector_Y, double *lambda, int *num_categories, double *mean, double **gradient);

int kernel_weighted_sum(int KERNEL_reg,int KERNEL_unordered_reg,int KERNEL_ordered_reg,int BANDWIDTH_reg,int num_obs_train,int num_obs_eval,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double **matrix_X_unordered_train,double **matrix_X_ordered_train,double **matrix_X_continuous_train,double **matrix_X_unordered_eval,double **matrix_X_ordered_eval,double **matrix_X_continuous_eval,double *vector_Y,double *vector_scale_factor,int *num_categories,double *kernel_sum);

int *alloc_vecu(int nobs);
void amoeba(double **p, double *y, int ndim, double ftol, double (*funk)(double *), int *nfunk);

void cv_func_bandwidth_gradient_regression(double *vector_scale_factor, double *vector_bandwidth_gradient);
void cv_func_bandwidth_gradient_regression_ml(double *vector_scale_factor, double *vector_bandwidth_gradient);
void cv_func_binary_choice_gradient(double *vector_scale_factor, double *vector_parameter_gradient);

void free_mat(double **x, int n);
void free_tmat(double **x);

void free_vector(double *v, int nl);

void linmin(int RESTRICT, int INTEGER, double *p_restrict, double *p, double *xi, int n, double tol, double small, int itmax, double *fret, double (*func)(double *));
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double));
void nrerror(char error_text[]);
void powell(int RESTRICT, int INTEGER, double *p_restrict, double *p, double **xi, int n, double ftol, double tol, double small, int itmax, int *iter, double *fret, double (*func)(double *));

void sort(int n, double ra[]);

int pgplot_xy_errorbars(int int_GENERATE, char *output, int num_obs, int num_var_unordered, int num_var_continuous, double *x_categorical_vec, double *x_continuous_vec, double *y_vec, double *y_std, char *x_label, char *y_label, char *title);

int compute_continuous_stddev(int int_LARGE, int num_obs, int num_var_continuous, int num_reg_continuous, double **matrix_Y_continuous, double **matrix_X_continuous, double *vector_continuous_stddev);

void initialize_nr_vector_scale_factor(int BANDWIDTH,int RANDOM,int seed,int int_large,int num_obs,int num_var_continuous,int num_var_unordered,int num_var_ordered,int num_reg_continuous,int num_reg_unordered,int num_reg_ordered,int kernel_yu,int kernel_xu,int int_use_starting_values,int scale_cat,double init_continuous,double nconfac,double ncatfac,int *num_categories,double *vector_continuous_stddev,double *vector_scale_factor,double lbc_init,double hbc_init,double c_init,double lbd_init,double hbd_init,double d_init,double ** matrix_x_continuous,double ** matrix_y_continuous);

int kernel_weights_conditional_convolution_cv(int int_WEIGHTS, int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered, double **matrix_Y_ordered, double **matrix_Y_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *lambda, double **matrix_bandwidth_var, double **matrix_bandwidth_reg, int *num_categories, double **matrix_categorical_vals, double **matrix_weights_K_x, double **matrix_weights_K_xy, double **matrix_weights_K_convol_y);

int check_valid_scale_factor_cv(int KERNEL, int KERNEL_unordered_liracine, int BANDWIDTH, int BANDWIDTH_den_ml, int REGRESSION_ML, int num_obs, int num_var_continuous, int num_var_unordered, int num_var_ordered, int num_reg_continuous, int num_reg_unordered, int num_reg_ordered, int *num_categories, double *vector_scale_factor);

int kernel_estimate_regression_categorical(int int_ll, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_reg, int num_obs_train, int num_obs_eval, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double **matrix_X_continuous_bandwidth, double *vector_Y, double *vector_Y_eval, double *vector_scale_factor, int *num_categories, double *mean, double **gradient, double *mean_stderr, double **gradient_stderr, double *R_squared, double *MSE, double *MAE, double *MAPE, double *CORR, double *SIGN);

int kernel_estimate_regression_categorical_tree_np(int int_ll,int KERNEL_reg,int KERNEL_unordered_reg,int KERNEL_ordered_reg,int BANDWIDTH_reg,int num_obs_train,int num_obs_eval,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double **matrix_X_unordered_train,double **matrix_X_ordered_train,double **matrix_X_continuous_train,double **matrix_X_unordered_eval,double **matrix_X_ordered_eval,double **matrix_X_continuous_eval,double *vector_Y,double *vector_Y_eval,double *vector_scale_factor,int *num_categories, double ** matrix_categorical_vals, double *mean,double **gradient,double *mean_stderr,double **gradient_stderr,double *R_squared,double *MSE,double *MAE,double *MAPE,double *CORR,double *SIGN);

double func_con_density_quantile(double *quantile);

int kernel_estimate_quantile(int gradient_compute, int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, double *quan, double *quan_stderr, double **quan_gradient, int seed, double ftol, double tol, double small, int itmax, int iMax_Num_Multistart, double zero, double lbc_dir, int dfc_dir, double c_dir,double initc_dir,double lbd_dir,double  hbd_dir,double  d_dir,double  initd_dir);

double ipow(double x, int n);

double kernel_estimate_regression_categorical_aic_c(int int_ll,int KERNEL_reg,int KERNEL_unordered_reg,int KERNEL_ordered_reg,int BANDWIDTH_reg,int num_obs,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double **matrix_X_unordered,double **matrix_X_ordered,double **matrix_X_continuous,double *vector_Y,double *vector_scale_factor,int *num_categories);

double cv_func_regression_categorical_aic_c(double *vector_scale_factor);

int simple_unique(int n, double * vector);
int unique(int num_obs, double *x);
void spinner(int num);

int kernel_weighted_sum_np(int * KERNEL_reg, int * KERNEL_unordered_reg, int * KERNEL_ordered_reg, const int BANDWIDTH_reg, const int num_obs_train, const int num_obs_eval, const int num_reg_unordered, const int num_reg_ordered, const int num_reg_continuous, const int leave_one_out, const int leave_one_out_offset, const int kernel_pow, const int bandwidth_divide, const int bandwidth_divide_weights, const int symmetric, const int gather_scatter, const int drop_one_train, const int drop_which_train, const int * const operator, const int permutation_operator, int do_score, int do_ocg, int * bpso, const int suppress_parallel, const int ncol_Y, const int ncol_W, const int int_TREE, const int do_partial_tree, KDT * const kdt, NL * const inl, int * const nld, int * const idx, double ** matrix_X_unordered_train,double **matrix_X_ordered_train,double **matrix_X_continuous_train,double **matrix_X_unordered_eval,double **matrix_X_ordered_eval,double **matrix_X_continuous_eval,double ** matrix_Y, double ** matrix_W, double * sgn, double *vector_scale_factor,int bandwidth_provided,double ** matrix_bw_train,double ** matrix_bw_eval,double * lambda_pre,int *num_categories,double ** matrix_categorical_vals, int ** matrix_ordered_indices, double * const restrict weighted_sum,  double * const restrict weighted_permutation_sum, double * const restrict kw);

int kernel_convolution_weighted_sum(int KERNEL_reg,int KERNEL_unordered_reg,int KERNEL_ordered_reg,int BANDWIDTH_reg,int num_obs_train,int num_obs_eval,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double **matrix_X_unordered_train,double **matrix_X_ordered_train,double **matrix_X_continuous_train,double **matrix_X_unordered_eval,double **matrix_X_ordered_eval,double **matrix_X_continuous_eval,double *vector_Y,double *vector_scale_factor,int *num_categories,double **matrix_categorical_vals,double *kernel_sum);

/*
int np_cuokernelv_loo_mlcv(int KERNEL, int uKERNEL, int oKERNEL,
                           int BANDWIDTH_den,
                           const int num,
                           const int unum, const int onum, const int cnum, 
                           double * const * const xtu,
                           double * const * const xto,
                           double * const * const xtc,
                           double * vector_scale_factor,
                           int * const ncat,                            
                           double * const cv);
*/

int np_kernel_estimate_con_density_categorical_convolution_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered, double **matrix_Y_ordered, double **matrix_Y_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_scale_factor, int *num_categories, double ** matrix_categorical_vals, double *cv);
double np_cv_func_con_density_categorical_ls(double *vector_scale_factor);
double np_cv_func_con_density_categorical_ml(double *vector_scale_factor);
double np_cv_func_con_density_categorical_ls_npksum(double *vector_scale_factor);

double np_kernel_estimate_distribution_ls_cv(int KERNEL_den,int KERNEL_den_unordered,int KERNEL_den_ordered,int BANDWIDTH_den,int num_obs_train,int num_obs_eval,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double memfac,double ** matrix_X_unordered_train,double ** matrix_X_ordered_train,double ** matrix_X_continuous_train,double ** matrix_X_unordered_eval,double ** matrix_X_ordered_eval,double ** matrix_X_continuous_eval,double * vsf,int * num_categories,double ** matrix_categorical_vals,double * cv);

int np_kernel_estimate_con_distribution_categorical_leave_one_out_ls_cv(int KERNEL_den,int KERNEL_unordered_den,int KERNEL_ordered_den,int KERNEL_reg,int KERNEL_unordered_reg,int KERNEL_ordered_reg,int BANDWIDTH_den,int64_t num_obs_train,int64_t num_obs_eval,int num_var_unordered,int num_var_ordered,int num_var_continuous,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double memfac, double **matrix_Y_unordered_train,double **matrix_Y_ordered_train,double **matrix_Y_continuous_train,double **matrix_X_unordered_train,double **matrix_X_ordered_train,double **matrix_X_continuous_train,double **matrix_XY_unordered_train, double **matrix_XY_ordered_train, double **matrix_XY_continuous_train, double **matrix_Y_unordered_eval,double **matrix_Y_ordered_eval,double **matrix_Y_continuous_eval,double *vector_scale_factor,int *num_categories,double **matrix_categorical_vals,double *cv);

int np_kernel_estimate_con_density_categorical_leave_one_out_ls_cv(int KERNEL_var,int KERNEL_unordered_var,int KERNEL_ordered_var,int KERNEL_reg,int KERNEL_unordered_reg,int KERNEL_ordered_reg,int BANDWIDTH_den,int64_t num_obs_train,int num_var_unordered,int num_var_ordered,int num_var_continuous,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double memfac,double **matrix_Y_unordered_train,double **matrix_Y_ordered_train,double **matrix_Y_continuous_train,double **matrix_X_unordered_train,double **matrix_X_ordered_train,double **matrix_X_continuous_train,double **matrix_XY_unordered_train,double **matrix_XY_ordered_train,double **matrix_XY_continuous_train,double *vector_scale_factor,int *num_categories,double **matrix_categorical_vals,double *cv);


void np_kernel_estimate_con_dens_dist_categorical(int KERNEL_Y,int KERNEL_unordered_Y,int KERNEL_ordered_Y,int KERNEL_X,int KERNEL_unordered_X,int KERNEL_ordered_X,int BANDWIDTH_den,int yop,int num_obs_train,int num_obs_eval,int num_Y_unordered,int num_Y_ordered,int num_Y_continuous,int num_X_unordered,int num_X_ordered,int num_X_continuous,double **matrix_XY_unordered_train, double **matrix_XY_ordered_train, double **matrix_XY_continuous_train, double **matrix_XY_unordered_eval, double **matrix_XY_ordered_eval, double **matrix_XY_continuous_eval, double *vector_scale_factor,int *num_categories,int *num_categories_XY, double ** matrix_categorical_vals, double ** matrix_categorical_vals_XY, double * kdf,double * kdf_stderr,double ** kdf_deriv,double ** kdf_deriv_stderr,double * log_likelihood);

void np_splitxy_vsf_mcv_nc(const int num_var_unordered, const int num_var_ordered, const int num_var_continuous, const int num_reg_unordered, const int num_reg_ordered, const int num_reg_continuous, const double * const vector_scale_factor, const int * const num_categories, double ** matrix_categorical_vals, double * vsf_x, double * vsf_y, double * vsf_xy, int * nc_x, int * nc_y, int * nc_xy, double ** mcv_x, double ** mcv_y, double ** mcv_xy);


void np_kernelop_xy(const int kernel_var_continuous, const int kernel_var_unordered, const int kernel_var_ordered, const int kernel_reg_continuous, const int kernel_reg_unordered, const int kernel_reg_ordered, const int operator_var, const int operator_reg, const int num_var_unordered, const int num_var_ordered, const int num_var_continuous, const int num_reg_unordered, const int num_reg_ordered, const int num_reg_continuous, int * kernel_cx, int * kernel_ux, int * kernel_ox, int * kernel_cy, int * kernel_uy, int * kernel_oy, int * kernel_cxy, int * kernel_uxy, int * kernel_oxy, int * operator_x, int * operator_y, int * operator_xy);


int is_valid_unordered_bw(double lambda,
                          int num_categories,
                          int kernel);
double max_unordered_bw(int num_categories,
                        int kernel);

// some general np and R-c interface related defines
#define safe_free(x) if((x) != NULL) free((x))

#define SF_NORMAL 0
#define SF_ARB    1

#define BW_FIXED   0
#define BW_GEN_NN  1
#define BW_ADAP_NN 2

#define IMULTI_TRUE  1
#define IMULTI_FALSE 0

#define RE_MIN_TRUE  0
#define RE_MIN_FALSE 1

#define IO_MIN_TRUE  1
#define IO_MIN_FALSE 0

#define LL_LC  0
#define LL_LL  1

#define OCG_UNO 0
#define OCG_ORD 1

#define BWM_CVML 0
#define BWM_CVLS 1
#define BWM_CVML_NP 2

#define DBWM_CVLS 0

#define NP_DO_DENS 1
#define NP_DO_DIST 0

#define OP_NOOP       -1
#define OP_NORMAL      0
#define OP_CONVOLUTION 1
#define OP_DERIVATIVE  2
#define OP_INTEGRAL    3

#define NP_NTREE 3

#define NP_TREE_FALSE 0
#define NP_TREE_TRUE  1

// offsets array : ordinary, convolution, derivative, cdfative
static const int OP_CFUN_OFFSETS[4] = { 0, 10, 20, 30 };
#define OP_NCFUN 40
static const int OP_UFUN_OFFSETS[4] = { 0, 2, 4, 0 };
static const int OP_OFUN_OFFSETS[4] = { 0, 3, 6, 9 };
// these defines are to facilitate accessing the continuous kernels in their various arrays

#define CK_GAUSS2 0
#define CK_GAUSS4 1
#define CK_GAUSS6 2
#define CK_GAUSS8 3
#define CK_EPAN2 4
#define CK_EPAN4 5
#define CK_EPAN6 6
#define CK_EPAN8 7
#define CK_UNIF 8
#define CK_TGAUSS2 9

#define ONE_OVER_SQRT_TWO_PI 0.39894228040143267794
#define SQRT_5 2.236067978

#define RBWM_CVAIC 0
#define RBWM_CVLS 1

#define RBW_NOBSI   0
#define RBW_IMULTII 1
#define RBW_NMULTII 2
#define RBW_USTARTI 3
#define RBW_LSFI    4
#define RBW_REGI  5
#define RBW_ITMAXI  6
#define RBW_REMINI  7
#define RBW_MINIOI  8
#define RBW_MI    9
#define RBW_CKRNEVI 10
#define RBW_UKRNEVI 11
#define RBW_OKRNEVI 12
#define RBW_NUNOI 13
#define RBW_NORDI 14
#define RBW_NCONI 15
#define RBW_LL 16
#define RBW_DOTREEI 17
#define RBW_SCATI 18
#define RBW_DFC_DIRI 19

#define RBW_FTOLD  0
#define RBW_TOLD   1
#define RBW_SMALLD 2
#define RBW_LBC_DIRD 3
#define RBW_C_DIRD 4
#define RBW_INITC_DIRD 5
#define RBW_LBD_DIRD 6
#define RBW_HBD_DIRD 7
#define RBW_D_DIRD 8
#define RBW_INITD_DIRD 9
#define RBW_LBC_INITD 10
#define RBW_HBC_INITD 11
#define RBW_C_INITD 12
#define RBW_LBD_INITD 13
#define RBW_HBD_INITD 14
#define RBW_D_INITD 15
#define RBW_NCONFD   16
#define RBW_NCATFD   17

#define MPI_RANKI 0
#define MPI_NUMPI 1

#define BW_NOBSI   0
#define BW_IMULTII 1
#define BW_NMULTII 2
#define BW_USTARTI 3
#define BW_LSFI    4
#define BW_DENI  5
#define BW_ITMAXI  6
#define BW_REMINI  7
#define BW_MINIOI  8
#define BW_MI    9
#define BW_CKRNEVI 10
#define BW_UKRNEVI 11
#define BW_OKRNEVI 12
#define BW_NUNOI 13
#define BW_NORDI 14
#define BW_NCONI 15
#define BW_OLDBW 16
#define BW_DOTREEI 17
#define BW_SCATI 18
#define BW_DFC_DIRI 19

#define BW_FTOLD  0
#define BW_TOLD   1
#define BW_SMALLD 2
#define BW_LBC_DIRD 3
#define BW_C_DIRD 4
#define BW_INITC_DIRD 5
#define BW_LBD_DIRD 6
#define BW_HBD_DIRD 7
#define BW_D_DIRD 8
#define BW_INITD_DIRD 9
#define BW_LBC_INITD 10
#define BW_HBC_INITD 11
#define BW_C_INITD 12
#define BW_LBD_INITD 13
#define BW_HBD_INITD 14
#define BW_D_INITD 15
#define BW_NCONFD     16
#define BW_NCATFD     17

// distribution defines
#define DBW_NOBSI   0
#define DBW_NEVALI  1
#define DBW_IMULTII 2
#define DBW_NMULTII 3
#define DBW_USTARTI 4
#define DBW_LSFI    5
#define DBW_DENI  6
#define DBW_ITMAXI  7
#define DBW_REMINI  8
#define DBW_MINIOI  9
#define DBW_MI    10
#define DBW_CKRNEVI 11
#define DBW_UKRNEVI 12
#define DBW_OKRNEVI 13
#define DBW_CDFONTRAIN 14
#define DBW_NUNOI 15
#define DBW_NORDI 16
#define DBW_NCONI 17
#define DBW_DOTREEI 18
#define DBW_SCATI 19
#define DBW_DFC_DIRI 20

#define DBW_FTOLD  0
#define DBW_TOLD   1
#define DBW_SMALLD 2
#define DBW_LBC_DIRD 3
#define DBW_C_DIRD 4
#define DBW_INITC_DIRD 5
#define DBW_LBD_DIRD 6
#define DBW_HBD_DIRD 7
#define DBW_D_DIRD 8
#define DBW_INITD_DIRD 9
#define DBW_LBC_INITD 10
#define DBW_HBC_INITD 11
#define DBW_C_INITD 12
#define DBW_LBD_INITD 13
#define DBW_HBD_INITD 14
#define DBW_D_INITD 15
#define DBW_NCONFD     16
#define DBW_NCATFD     17
#define DBW_MEMORYD     18


#define CBW_NOBSI   0
#define CBW_IMULTII 1
#define CBW_NMULTII 2
#define CBW_USTARTI 3
#define CBW_LSFI    4
#define CBW_DENI  5
#define CBW_ITMAXI  6
#define CBW_REMINI  7
#define CBW_MINIOI  8
#define CBW_MI    9
#define CBW_CXKRNEVI 10
#define CBW_CYKRNEVI 11
#define CBW_UXKRNEVI 12
#define CBW_UYKRNEVI 13
#define CBW_OXKRNEVI 14
#define CBW_OYKRNEVI 15
#define CBW_CNUNOI 16
#define CBW_CNORDI 17
#define CBW_CNCONI 18
#define CBW_UNUNOI 19
#define CBW_UNORDI 20
#define CBW_UNCONI 21
#define CBW_FASTI 22
#define CBW_OLDI 23
#define CBW_TREEI 24
#define CBW_SCATI 25
#define CBW_DFC_DIRI 26

#define CBW_FTOLD  0
#define CBW_TOLD   1
#define CBW_SMALLD 2
#define CBW_MEMFACD 3
#define CBW_LBC_DIRD 4
#define CBW_C_DIRD 5
#define CBW_INITC_DIRD 6
#define CBW_LBD_DIRD 7
#define CBW_HBD_DIRD 8
#define CBW_D_DIRD 9
#define CBW_INITD_DIRD 10
#define CBW_LBC_INITD 11
#define CBW_HBC_INITD 12
#define CBW_C_INITD 13
#define CBW_LBD_INITD 14
#define CBW_HBD_INITD 15
#define CBW_D_INITD 16
#define CBW_NCONFD     17
#define CBW_NCATFD     18


#define CBWM_CVML 0
#define CBWM_CVLS 1
#define CBWM_NPLS 2
#define CBWM_CCDF 3

#define CBW_MINOBS 1024

// ccdf defines
#define CDBW_NOBSI   0
#define CDBW_NEVALI  1
#define CDBW_IMULTII 2
#define CDBW_NMULTII 3
#define CDBW_USTARTI 4
#define CDBW_LSFI    5
#define CDBW_DENI  6
#define CDBW_ITMAXI  7
#define CDBW_REMINI  8
#define CDBW_MINIOI  9
#define CDBW_MI    10
#define CDBW_CXKRNEVI 11
#define CDBW_CYKRNEVI 12
#define CDBW_UXKRNEVI 13
#define CDBW_UYKRNEVI 14
#define CDBW_OXKRNEVI 15
#define CDBW_OYKRNEVI 16
#define CDBW_CNUNOI 17
#define CDBW_CNORDI 18
#define CDBW_CNCONI 19
#define CDBW_UNUNOI 20
#define CDBW_UNORDI 21
#define CDBW_UNCONI 22
#define CDBW_CDFONTRAIN 23
#define CDBW_TREEI 24
#define CDBW_SCATI 25
#define CDBW_DFC_DIRI 26

#define CDBW_FTOLD  0
#define CDBW_TOLD   1
#define CDBW_SMALLD 2
#define CDBW_MEMFACD 3
#define CDBW_LBC_DIRD 4
#define CDBW_C_DIRD 5
#define CDBW_INITC_DIRD 6
#define CDBW_LBD_DIRD 7
#define CDBW_HBD_DIRD 8
#define CDBW_D_DIRD 9
#define CDBW_INITD_DIRD 10
#define CDBW_LBC_INITD 11
#define CDBW_HBC_INITD 12
#define CDBW_C_INITD 13
#define CDBW_LBD_INITD 14
#define CDBW_HBD_INITD 15
#define CDBW_D_INITD 16
#define CDBW_NCONFD     17
#define CDBW_NCATFD     18

#define CDBWM_CVLS 0

//
#define CD_TNOBSI 0
#define CD_ENOBSI   1
#define CD_LSFI    2
#define CD_DENI  3
#define CD_MINIOI  4
#define CD_CXKRNEVI 5
#define CD_CYKRNEVI 6
#define CD_UXKRNEVI 7
#define CD_UYKRNEVI 8
#define CD_OXKRNEVI 9
#define CD_OYKRNEVI 10
#define CD_CNUNOI 11
#define CD_CNORDI 12
#define CD_CNCONI 13
#define CD_UNUNOI 14
#define CD_UNORDI 15
#define CD_UNCONI 16
#define CD_TISEI 17
#define CD_GRAD 18
#define CD_YMLEVI 19
#define CD_XMLEVI 20
#define CD_DODENI 21
#define CD_TREEI 22

#define DEN_TNOBSI   0
#define DEN_ENOBSI   1
#define DEN_NUNOI 2
#define DEN_NORDI 3
#define DEN_NCONI 4
#define DEN_LSFI 5
#define DEN_DENI 6
#define DEN_MINIOI 7
#define DEN_CKRNEVI    8
#define DEN_UKRNEVI    9
#define DEN_OKRNEVI    10
#define DEN_TISEI 11
#define DEN_MLEVI 12
#define DEN_DODENI 13
#define DEN_OLDI 14
#define DEN_TREEI 15

#define REG_TNOBSI   0
#define REG_ENOBSI   1
#define REG_NUNOI 2
#define REG_NORDI 3
#define REG_NCONI 4
#define REG_LSFI 5
#define REG_BWI 6
#define REG_MINIOI 7
#define REG_CKRNEVI    8
#define REG_UKRNEVI 9
#define REG_OKRNEVI 10
#define REG_EY 11
#define REG_GRAD 12
#define REG_LL 13
#define REG_TISEI 14
#define REG_MLEVI 15
#define REG_DOTREEI 16
#define REG_OLDREGI 17

#define KWS_TNOBSI   0
#define KWS_ENOBSI 1
#define KWS_NUNOI 2
#define KWS_NORDI 3
#define KWS_NCONI 4
#define KWS_LSFI 5
#define KWS_BWI 6
#define KWS_MINIOI 7
#define KWS_CKRNEVI    8
#define KWS_UKRNEVI 9
#define KWS_OKRNEVI 10
#define KWS_TISEI 11
#define KWS_LOOI 12
#define KWS_BDIVI 13
#define KWS_MLEVI 14
#define KWS_WNCOLI 15
#define KWS_YNCOLI 16
#define KWS_DOTREEI 17
#define KWS_RKWI 18
#define KWS_POPI 19
#define KWS_PSCOREI 20
#define KWS_POCGI 21

#define CQ_TNOBSI 0
#define CQ_ENOBSI   1
#define CQ_LSFI    2
#define CQ_DENI  3
#define CQ_MINIOI  4
#define CQ_CXKRNEVI 5
#define CQ_CYKRNEVI 6
#define CQ_UXKRNEVI 7
#define CQ_UYKRNEVI 8
#define CQ_OXKRNEVI 9
#define CQ_OYKRNEVI 10
#define CQ_CNUNOI 11
#define CQ_CNORDI 12
#define CQ_CNCONI 13
#define CQ_UNUNOI 14
#define CQ_UNORDI 15
#define CQ_UNCONI 16
#define CQ_TISEI 17
#define CQ_GRADI 18
#define CQ_ITMAXI 19
#define CQ_MLEVI 20
#define CQ_NMULTII 21
#define CQ_DFC_DIRI 22

#define CQ_FTOLD  0
#define CQ_TOLD   1
#define CQ_SMALLD 2
#define CQ_LBC_DIRD 3
#define CQ_C_DIRD 4
#define CQ_INITC_DIRD 5
#define CQ_LBD_DIRD 6
#define CQ_HBD_DIRD 7
#define CQ_D_DIRD 8
#define CQ_INITD_DIRD 9

#define TG2_B     0
#define TG2_ALPHA 1
#define TG2_C0    2
#define TG2_A0    3
#define TG2_A1    4
#define TG2_A2    5
#define TG2_K     6
#define TG2_K2    7
#define TG2_K22   8
#define TG2_KM    9

#define KWSNP_ERR_NOEVAL 1
#define KWSNP_ERR_BADBW 2
#define KWSNP_ERR_BADINVOC 3

#define UKERNEL_UAA 0
#define UKERNEL_ULIR 1

#endif // _NP__
