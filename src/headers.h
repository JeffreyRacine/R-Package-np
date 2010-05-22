#ifndef _NP_
#define _NP_

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <time.h>

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
#define NZD(a) ((a) < 0 ? (MIN(-DBL_MIN,a)) : (MAX(DBL_MIN,a)))

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
double standerrd(int n, double *vector);

int compute_nn_distance(int num_obs, double *vector_data, int int_k_nn, double *nn_distance);
int compute_nn_distance_train_eval(int num_obs_train, int num_obs_eval, double *vector_data_train, double *vector_data_eval, int int_k_nn, double *nn_distance);

int determine_categorical_vals(int num_obs, int num_var_unordered, int num_var_ordered, int num_reg_unordered, int num_reg_ordered, double **matrix_Y_unordered, double **matrix_Y_ordered, double **matrix_X_unordered, double **matrix_X_ordered, int *num_categories, double **matrix_categorical_vals);

int fround(double x);

int initialize_kernel_density_asymptotic_constants(int KERNEL, int num_var_continuous, double *INT_KERNEL_P, double *K_INT_KERNEL_P);  
int initialize_kernel_regression_asymptotic_constants(int KERNEL, int num_reg_continuous, double *INT_KERNEL_P, double *K_INT_KERNEL_P, double *INT_KERNEL_PM_HALF, double *DIFF_KER_PPM); 
int initialize_nr_hessian(int num_var, double **matrix_y);

int kernel_bandwidth(int KERNEL, int BANDWIDTH, int num_obs_train, int num_obs_eval, int num_var_cont, int num_var_un, int num_var_or, int num_reg_cont, int num_reg_un, int num_reg_or, double *vector_scale_factor, double **matrix_Y_train, double **matrix_Y_eval, double **matrix_X_train, double **matrix_X_eval, double **matrix_bandwidth_Y, double **matrix_bandwidth_X, double *vector_lambda, double **matrix_bandwidth_deriv);
int kernel_bandwidth_mean(int KERNEL, int BANDWIDTH, int num_obs_train, int num_obs_eval, int num_var_cont, int num_var_un, int num_var_or, int num_reg_cont, int num_reg_un, int num_reg_or, double *vector_scale_factor, double **matrix_Y_train, double **matrix_Y_eval, double **matrix_X_train, double **matrix_X_eval, double **matrix_bandwidth_Y, double **matrix_bandwidth_X, double *vector_lambda);

int kernel_estimate_categorical_gradient_ocg_fast(int int_COMPUTE, int *var_index_int, int num_var_test_int, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_reg, int int_ll, int int_ordered_categorical_gradient, int num_obs_train, int num_obs_eval, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double *vector_Y, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double **matrix_bandwidth, double **matrix_bandwidth_deriv, double *lambda, int *num_categories, double **matrix_categorical_vals, double *mean, double **gradient_categorical);

int kernel_estimate_con_density_categorical(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double *pdf, double *pdf_stderr, double *log_likelihood);

int kernel_estimate_con_density_categorical_convolution_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered, double **matrix_Y_ordered, double **matrix_Y_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cv);

int kernel_estimate_con_density_categorical_gradient(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double *pdf, double *pdf_stderr, double **pdf_deriv, double **pdf_deriv_stderr, double *log_likelihood);

int kernel_estimate_con_density_categorical_gradient_categorical(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, int int_ordered_categorical_gradient, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, double **matrix_categorical_vals, int *num_categories, double *pdf, double **pdf_deriv, double **pdf_deriv_stderr);

int kernel_estimate_con_density_categorical_leave_one_out_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered, double **matrix_Y_ordered, double **matrix_Y_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_scale_factor, int *num_categories, double *cv);

int kernel_estimate_con_distribution_categorical(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cdf, double *cdf_stderr, double small, int itmax);

int kernel_estimate_con_distribution_categorical_leave_one_out(int KERNEL_den,int KERNEL_unordered_den,int KERNEL_ordered_den,int KERNEL_reg,int KERNEL_unordered_reg,int KERNEL_ordered_reg,int BANDWIDTH_den,int num_obs_train,int num_obs_eval,int num_var_unordered,int num_var_ordered,int num_var_continuous,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double **matrix_Y_unordered_train,double **matrix_Y_ordered_train,double **matrix_Y_continuous_train,double **matrix_Y_unordered_eval,double **matrix_Y_ordered_eval,double **matrix_Y_continuous_eval,double **matrix_X_unordered_train,double **matrix_X_ordered_train,double **matrix_X_continuous_train,double **matrix_X_unordered_eval,double **matrix_X_ordered_eval,double **matrix_X_continuous_eval,double *vector_scale_factor,int *num_categories,double **matrix_categorical_vals,double *cdf,double small,int itmax);

int indfunc(double a);
double cv_func_con_distribution_categorical_ccdf(double *vector_scale_factor);

int kernel_estimate_con_distribution_categorical_leave_one_out_ccdf(int KERNEL_den,int KERNEL_unordered_den,int KERNEL_ordered_den,int KERNEL_reg,int KERNEL_unordered_reg,int KERNEL_ordered_reg,int BANDWIDTH_den,int num_obs_train,int num_var_unordered,int num_var_ordered,int num_var_continuous,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double **matrix_Y_unordered_train,double **matrix_Y_ordered_train,double **matrix_Y_continuous_train,double **matrix_X_unordered_train,double **matrix_X_ordered_train,double **matrix_X_continuous_train,double *vector_scale_factor,int *num_categories,double **matrix_categorical_vals,double *cv,double small,int itmax);

int kernel_estimate_con_distribution_categorical_no_mpi(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cdf, double *cdf_stderr, double small, int itmax);

int kernel_estimate_con_distribution_categorical_gradient(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cdf, double *cdf_stderr, double **cdf_deriv, double **cdf_deriv_stderr, double small, int itmax);

int kernel_estimate_con_distribution_categorical_gradient_categorical(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, int int_ordered_categorical_gradient, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, double **matrix_categorical_vals, int *num_categories, double *cdf, double **cdf_deriv, double **cdf_deriv_stderr, double small, int itmax);

int kernel_estimate_density_categorical(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, int *num_categories, double *pdf, double *pdf_stderr, double *log_likelihood);

int kernel_estimate_density_categorical_convolution_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_scale_factor, int *num_categories, double **matrix_categorical_vals, double *cv);

int kernel_estimate_density_categorical_leave_one_out_cv(int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *vector_scale_factor, int *num_categories, double *cv);

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

void free_vector(double *v, int nl, int nh);

void linmin(int RESTRICT, int INTEGER, double *p_restrict, double *p, double *xi, int n, double tol, double small, int itmax, double *fret, double (*func)(double *));
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double small, double (*func)(double));
void nrerror(char error_text[]);
void powell(int RESTRICT, int INTEGER, double *p_restrict, double *p, double **xi, int n, double ftol, double tol, double small, int itmax, int *iter, double *fret, double (*func)(double *));

void sort(int n, double ra[]);

int pgplot_xy_errorbars(int int_GENERATE, char *output, int num_obs, int num_var_unordered, int num_var_continuous, double *x_categorical_vec, double *x_continuous_vec, double *y_vec, double *y_std, char *x_label, char *y_label, char *title);

int compute_continuous_stddev(int int_LARGE, int num_obs, int num_var_continuous, int num_reg_continuous, double **matrix_Y_continuous, double **matrix_X_continuous, double *vector_continuous_stddev);

int initialize_nr_vector_scale_factor(int BANDWIDTH, int BANDWIDTH_den, int RANDOM, int seed, int REGRESSION_ML, int int_large, int num_obs, int num_var_continuous, int num_var_unordered, int num_var_ordered, int num_reg_continuous, int num_reg_unordered, int num_reg_ordered, double **matrix_Y_continuous, double **matrix_X_continuous, int int_use_starting_values, double init_continuous, int *num_categories, double *vector_continuous_stddev, double *vector_scale_factor);

int kernel_weights_conditional_convolution_cv(int int_WEIGHTS, int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_den, int num_obs, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered, double **matrix_Y_ordered, double **matrix_Y_continuous, double **matrix_X_unordered, double **matrix_X_ordered, double **matrix_X_continuous, double *lambda, double **matrix_bandwidth_var, double **matrix_bandwidth_reg, int *num_categories, double **matrix_categorical_vals, double **matrix_weights_K_x, double **matrix_weights_K_xy, double **matrix_weights_K_convol_y);

int check_valid_scale_factor_cv(int KERNEL, int KERNEL_unordered_liracine, int BANDWIDTH, int BANDWIDTH_den_ml, int REGRESSION_ML, int num_obs, int num_var_continuous, int num_var_unordered, int num_var_ordered, int num_reg_continuous, int num_reg_unordered, int num_reg_ordered, int *num_categories, double *vector_scale_factor);

int kernel_estimate_regression_categorical(int int_ll, int KERNEL_reg, int KERNEL_unordered_reg, int KERNEL_ordered_reg, int BANDWIDTH_reg, int num_obs_train, int num_obs_eval, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double **matrix_X_continuous_bandwidth, double *vector_Y, double *vector_Y_eval, double *vector_scale_factor, int *num_categories, double *mean, double **gradient, double *mean_stderr, double **gradient_stderr, double *R_squared, double *MSE, double *MAE, double *MAPE, double *CORR, double *SIGN);

double func_con_density_quantile(double *quantile);

int kernel_estimate_quantile(int gradient_compute, int KERNEL_den, int KERNEL_unordered_den, int KERNEL_ordered_den, int BANDWIDTH_den, int num_obs_train, int num_obs_eval, int num_var_unordered, int num_var_ordered, int num_var_continuous, int num_reg_unordered, int num_reg_ordered, int num_reg_continuous, double **matrix_Y_unordered_train, double **matrix_Y_ordered_train, double **matrix_Y_continuous_train, double **matrix_Y_unordered_eval, double **matrix_Y_ordered_eval, double **matrix_Y_continuous_eval, double **matrix_X_unordered_train, double **matrix_X_ordered_train, double **matrix_X_continuous_train, double **matrix_X_unordered_eval, double **matrix_X_ordered_eval, double **matrix_X_continuous_eval, double *vector_scale_factor, double *quan, double *quan_stderr, double **quan_gradient, int seed, double ftol, double tol, double small, int itmax, int iMax_Num_Multistart, double zero);

double ipow(double x, int n);

double kernel_estimate_regression_categorical_aic_c(int int_ll,int KERNEL_reg,int KERNEL_unordered_reg,int KERNEL_ordered_reg,int BANDWIDTH_reg,int num_obs,int num_reg_unordered,int num_reg_ordered,int num_reg_continuous,double **matrix_X_unordered,double **matrix_X_ordered,double **matrix_X_continuous,double *vector_Y,double *vector_scale_factor,int *num_categories);

double cv_func_regression_categorical_aic_c(double *vector_scale_factor);

int unique(int num_obs, double *x);
void spinner(int num);

int kernel_weighted_sum_np(const int KERNEL_reg, const int KERNEL_unordered_reg, const int KERNEL_ordered_reg, const int BANDWIDTH_reg, const int num_obs_train, const int num_obs_eval, const int num_reg_unordered, const int num_reg_ordered, const int num_reg_continuous, const int leave_one_out, const int kernel_pow, const int bandwidth_divide, const int do_smooth_coef_weights, const int symmetric, const int gather_scatter, const int operator, double * const * const matrix_X_unordered_train,double **matrix_X_ordered_train,double **matrix_X_continuous_train,double **matrix_X_unordered_eval,double **matrix_X_ordered_eval,double **matrix_X_continuous_eval,double ** matrix_Y, double ** matrix_W, double * sgn, double *vector_scale_factor,int *num_categories,double ** matrix_categorical_vals, double *weighted_sum);

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

#define NP_DO_DENS 1
#define NP_DO_DIST 0

#define OP_NORMAL      0
#define OP_CONVOLUTION 1
#define OP_DERIVATIVE  2
#define OP_INTEGRAL    3

static const int OP_CFUN_OFFSETS[4] = { 0, 9, 17, 26 };
static const int OP_UFUN_OFFSETS[4] = { 0, 2, 0, 0 };

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

#define RBW_FTOLD  0
#define RBW_TOLD   1
#define RBW_SMALLD 2

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
#define BW_NUNOI 11
#define BW_NORDI 12
#define BW_NCONI 13

#define BW_FTOLD  0
#define BW_TOLD   1
#define BW_SMALLD 2

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
#define CBW_AUTOI 23

#define CBW_FTOLD  0
#define CBW_TOLD   1
#define CBW_SMALLD 2

#define CBWM_CVML 0
#define CBWM_CVLS 1
#define CBWM_NPLS 2
#define CBWM_CCDF 3

#define CBW_MINOBS 1024

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

#define DEN_TNOBSI   0
#define DEN_ENOBSI   1
#define DEN_NUNOI 2
#define DEN_NORDI 3
#define DEN_NCONI 4
#define DEN_LSFI 5
#define DEN_DENI 6
#define DEN_MINIOI 7
#define DEN_CKRNEVI    8
#define DEN_TISEI 9
#define DEN_MLEVI 10
#define DEN_DODENI 11

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
#define KWS_IPOWI 13
#define KWS_BDIVI 14
#define KWS_OPI 15
#define KWS_MLEVI 16
#define KWS_SCOEFI 17
#define KWS_WNCOLI 18
#define KWS_YNCOLI 19

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

#define CQ_FTOLD  0
#define CQ_TOLD   1
#define CQ_SMALLD 2

#endif // _NP_
