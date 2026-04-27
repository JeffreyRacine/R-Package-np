/* Copyright (C) J. Racine, 1995-2001 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <string.h>
#include <time.h>

#include <R.h>
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <Rinternals.h>

#include "headers.h"
#include "linalg.h"
#include "gsl_bspline.h"

#include "hash.h"
#include "tree.h"

#include <assert.h>

#include <inttypes.h>
#ifdef MPI2

#include "mpi.h"

extern  int my_rank;
extern  int source;
extern  int dest;
extern  int tag;
extern  int iNum_Processors;
extern  int iSeed_my_rank;
extern  MPI_Status status;
extern MPI_Comm	*comm;
#endif

static int np_mpi_rank_failure_injected(const char *env_name){
#ifdef MPI2
  const char *value = getenv(env_name);
  char *end = NULL;
  long rank;

  if((value == NULL) || (value[0] == '\0') || strcmp(value, "0") == 0)
    return 0;
  if(strcmp(value, "all") == 0)
    return 1;

  errno = 0;
  rank = strtol(value, &end, 10);
  if((errno != 0) || (end == value) || (*end != '\0'))
    return 0;

  return ((int)rank == my_rank);
#else
  (void)env_name;
  return 0;
#endif
}

extern int int_DEBUG;
extern int int_VERBOSE;
extern int int_TAYLOR;
extern int int_WEIGHTS;
extern int int_LARGE_SF;

extern int int_TREE_X;
extern int int_TREE_Y;
extern int int_TREE_XY;

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
extern int int_cker_bound_extern;
extern double *vector_ckerlb_extern;
extern double *vector_ckerub_extern;
extern int int_cxker_bound_extern;
extern int int_cyker_bound_extern;
extern int int_cxyker_bound_extern;
extern double *vector_cxkerlb_extern;
extern double *vector_cxkerub_extern;
extern double *vector_cykerlb_extern;
extern double *vector_cykerub_extern;
extern double *vector_cxykerlb_extern;
extern double *vector_cxykerub_extern;
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
extern int *vector_glp_degree_extern;
extern int *vector_glp_gradient_order_extern;
extern int int_glp_bernstein_extern;
extern int int_glp_basis_extern;
extern int int_bounded_cvls_quadrature_grid_extern;
extern int int_bounded_cvls_quadrature_points_extern;
extern double double_bounded_cvls_quadrature_extend_factor_extern;
extern double double_bounded_cvls_quadrature_ratios_extern[3];

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

/* Statics for dependence metric */

extern int num_lag_extern;
extern int int_lag_extern;
extern int int_iter_extern;

extern double *vector_scale_factor_dep_met_bivar_extern;
extern double *vector_scale_factor_dep_met_univar_extern;
extern double *vector_scale_factor_dep_met_univar_lag_extern;

extern double y_min_extern;
extern double y_max_extern;

// tree
extern KDT * kdt_extern_X;
extern KDT * kdt_extern_Y;
extern KDT * kdt_extern_XY;

extern int * ipt_extern_X;
extern int * ipt_extern_Y;
extern int * ipt_extern_XY;

static double np_blas_ddot_int(const int n, const double *x, const double *y){
  const int inc = 1;
  if((n <= 0) || (x == NULL) || (y == NULL))
    return 0.0;
  return F77_CALL(ddot)(&n, x, &inc, y, &inc);
}

static void np_blas_dgemm_tn_int(const int m,
                                 const int n,
                                 const int k,
                                 const double *a,
                                 const double *b,
                                 double *c){
  const char transa = 'T';
  const char transb = 'N';
  const double alpha = 1.0;
  const double beta = 0.0;

  if((m <= 0) || (n <= 0) || (k <= 0) || (a == NULL) || (b == NULL) || (c == NULL))
    return;

  F77_CALL(dgemm)(&transa,
                  &transb,
                  &m,
                  &n,
                  &k,
                  &alpha,
                  a,
                  &k,
                  b,
                  &k,
                  &beta,
                  c,
                  &m
                  FCONE FCONE);
}

extern int * ipt_lookup_extern_X;
extern int * ipt_lookup_extern_Y;
extern int * ipt_lookup_extern_XY;

extern int *num_categories_extern_XY;
extern int *num_categories_extern_X;
extern int *num_categories_extern_Y;
extern int cdfontrain_extern;

extern double ** matrix_categorical_vals_extern_X;
extern double ** matrix_categorical_vals_extern_Y;
extern double ** matrix_categorical_vals_extern_XY;

static size_t np_jksum_size_mul_or_die(size_t a, size_t b, const char *what)
{
  if ((a != 0) && (b > (SIZE_MAX / a)))
    error("%s: allocation size overflow", what);
  return a * b;
}

static size_t np_jksum_size_mul3_or_die(size_t a, size_t b, size_t c, const char *what)
{
  return np_jksum_size_mul_or_die(np_jksum_size_mul_or_die(a, b, what), c, what);
}

static void *np_jksum_malloc_bytes_or_die(size_t nbytes, const char *what)
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = malloc(nbytes);
  if (ptr == NULL)
    error("%s: memory allocation failed", what);

  return ptr;
}

static void *np_jksum_malloc_array_or_die(size_t count, size_t elem_size, const char *what)
{
  return np_jksum_malloc_bytes_or_die(np_jksum_size_mul_or_die(count, elem_size, what), what);
}

static void *np_jksum_malloc_array3_or_die(size_t a, size_t b, size_t c, const char *what)
{
  return np_jksum_malloc_bytes_or_die(np_jksum_size_mul3_or_die(a, b, c, what), what);
}

int kernel_convolution_weighted_sum(
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_reg,
int num_obs_train,
int num_obs_eval,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_Y,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *kernel_sum)
{

	/* This function takes a vector Y and returns a convolution kernel
		  weighted sum. By default Y should be a vector of ones (simply
		  compute the kernel sum). This function will allow users to `roll
		  their own' with mixed data convolution kernel sums. */

	/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel;
	double sum_y_ker;

	double *lambda;
	double **matrix_bandwidth = NULL;

#ifndef MPI2
  double * psum;
#endif
	double *py;

#ifdef MPI2
	int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
#endif

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

	if((BANDWIDTH_reg == 0)||(BANDWIDTH_reg == 1))
	{
		matrix_bandwidth = alloc_matd(num_obs_eval,num_reg_continuous);
	}
	else if(BANDWIDTH_reg == 2)
	{
		matrix_bandwidth = alloc_matd(num_obs_train,num_reg_continuous);
	}

	/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

	if(kernel_bandwidth_mean(
		KERNEL_reg,
		BANDWIDTH_reg,
		num_obs_train,
		num_obs_eval,
		0,
		0,
		0,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
    0, // do not suppress parallel
		vector_scale_factor,
		matrix_X_continuous_train,	 /* Not used */
		matrix_X_continuous_eval,		 /* Not used */
		matrix_X_continuous_train,
		matrix_X_continuous_eval,
		matrix_bandwidth,						 /* Not used */
		matrix_bandwidth,
		lambda) == 1)
	{
		error("\n** Error: invalid bandwidth.");
	}

#ifndef MPI2
	if(BANDWIDTH_reg == 0)
	{

		psum = &kernel_sum[0];

		for(j=0; j < num_obs_eval; j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0],
						matrix_bandwidth[l][0],
						matrix_bandwidth[l][0]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			*psum++ = sum_y_ker;

		}

	}
	else if(BANDWIDTH_reg == 1)
	{

		psum = &kernel_sum[0];

		for(j=0; j < num_obs_eval; j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j],
						matrix_bandwidth[l][i],
						matrix_bandwidth[l][j]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			*psum++ = sum_y_ker;

		}

	}
	else
	{

		psum = &kernel_sum[0];

		for(j=0; j < num_obs_eval; j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i],
						matrix_bandwidth[l][j],
						matrix_bandwidth[l][i]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			*psum++ = sum_y_ker;

		}

	}
#endif

#ifdef MPI2

	if(BANDWIDTH_reg == 0)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0],
						matrix_bandwidth[l][0],
						matrix_bandwidth[l][0]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			kernel_sum[j-my_rank*stride] = sum_y_ker;

		}

	}
	else if(BANDWIDTH_reg == 1)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j],
						matrix_bandwidth[l][i],
						matrix_bandwidth[l][j]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			kernel_sum[j-my_rank*stride] = sum_y_ker;

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i],
						matrix_bandwidth[l][j],
						matrix_bandwidth[l][i]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			kernel_sum[j-my_rank*stride] = sum_y_ker;

		}

	}

	if(my_rank == 0)
		MPI_Gather(MPI_IN_PLACE, stride, MPI_DOUBLE, kernel_sum, stride, MPI_DOUBLE, 0, comm[1]);
	else
		MPI_Gather(kernel_sum, stride, MPI_DOUBLE, NULL, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(kernel_sum, num_obs_eval, MPI_DOUBLE, 0, comm[1]);
#endif

	free(lambda);

	free_mat(matrix_bandwidth,num_reg_continuous);

	return(0);

}

extern double np_tgauss2_b, np_tgauss2_alpha, np_tgauss2_c0;
// convolution kernel constants
extern double np_tgauss2_a0, np_tgauss2_a1, np_tgauss2_a2;


double np_tgauss2(const double z){
  return (fabs(z) >= np_tgauss2_b) ? 0.0 : np_tgauss2_alpha*ONE_OVER_SQRT_TWO_PI*exp(-0.5*z*z) - np_tgauss2_c0;
}

double np_gauss2(const double z){
  return ONE_OVER_SQRT_TWO_PI*exp(-0.5*z*z);
}

double np_gauss4(const double z){
  return ONE_OVER_SQRT_TWO_PI*(1.5-0.5*z*z)*exp(-0.5*z*z);
}

double np_gauss6(const double z){
  const double z2 = z*z;
  return ONE_OVER_SQRT_TWO_PI*(1.875+z2*(z2*0.125-1.25))*exp(-0.5*z2);
}

double np_gauss8(const double z){
  const double z2 = z*z;
  return ONE_OVER_SQRT_TWO_PI*(2.1875+z2*(-2.1875+z2*(0.4375-z2*0.02083333333)))*exp(-0.5*z2);
}

double np_epan2(const double z){
  return (z*z < 5.0)?((double)(0.33541019662496845446-0.067082039324993690892*z*z)):0.0;
}

double np_epan4(const double z){
  const double z2 = z*z;
  return (z2 < 5.0)?((double)(0.008385254916*(-15.0+7.0*z2)*(-5.0+z2))):0.0;
}

double np_epan6(const double z){
  const double z2 = z*z;
  return (z2 < 5.0)?((double)(0.33541019662496845446*(2.734375+z2*(-3.28125+0.721875*z2))*(1.0-0.2*z2))):0.0;
}

double np_epan8(const double z){
  const double z2 = z*z;
  return (z2 < 5.0)?((double)(0.33541019662496845446*(3.5888671875+z2*(-7.8955078125+z2*(4.1056640625-0.5865234375*z2)))*(1.0-0.2*z2))):0.0;
}

double np_rect(const double z){
  return (z*z < 1.0)?0.5:0.0;
}

double np_uaa(const int same_cat,const double lambda, const int c){
  if(c < 2)
    return same_cat ? 1.0 : 0.0;
  return (same_cat)?(1.0-lambda):lambda/((double)c-1.0);
}

double np_score_uaa(const int same_cat,const double lambda, const int c){
  if(c < 2)
    return 0.0;
  return (same_cat)?-1.0:(1.0/((double)c-1.0));
}

double np_uli_racine(const int same_cat, const double lambda, const int c){
  return (same_cat)?1.0:lambda;
}

double np_unli_racine(const int same_cat, const double lambda, const int c){
  return ((same_cat)?1.0:lambda)/((c-1.0)*lambda + 1.0);
}

double np_score_uli_racine(const int same_cat, const double lambda, const int c){
  return (same_cat)?0.0:1.0;
}

double np_score_unli_racine(const int same_cat, const double lambda, const int c){
  const double inorm = 1.0/((c-1.0)*lambda + 1.0);
  return (same_cat)?(1.0-c)*inorm*inorm:inorm*(lambda*(1.0-c)*inorm + 1.0);
}

double np_owang_van_ryzin(const double x, const double y, const double lambda, const double cl, const double ch){
  return (x == y)?(1.0-lambda):ipow(lambda, (int)fabs(x-y))*(1.0-lambda)*0.5;
}

double np_score_owang_van_ryzin(const double x, const double y, const double lambda, const double cl, const double ch){
  const int cxy = (int)fabs(x-y);
  if(cxy == 0) return -1.0;
  if(lambda == 0.0) return (cxy == 1) ? 0.5 : 0.0;
  return 0.5*ipow(lambda, cxy - 1)*(cxy - 2.0*lambda);
}

double np_oli_racine(const double x, const double y, const double lambda, const double cl, const double ch){
  return ipow(lambda, (int)fabs(x-y));
}

double np_score_oli_racine(const double x, const double y, const double lambda, const double cl, const double ch){
  const int cxy = (int)fabs(x-y);
  if (cxy == 0) return 0.0;
  return (cxy * ipow(lambda, cxy - 1));
}

double np_onli_racine(const double x, const double y, const double lambda, const double cl, const double ch){
  return ipow(lambda, (int)fabs(x-y))*(1.0 - lambda)/(1.0 + lambda);
}

double np_score_onli_racine(const double x, const double y, const double lambda, const double cl, const double ch){
  const int cxy = (int)fabs(x-y);
  if(cxy == 0) return -2.0;
  return ipow(lambda, cxy - 1)*(cxy*(1.0 - lambda*lambda) - 2.0*lambda);
}

static inline void np_orly_term_deriv(const int d, const double lambda, double *term, double *dterm){
  if(d <= 0){
    *term = 1.0;
    *dterm = 0.0;
    return;
  }

  if(lambda == 0.0){
    *term = 0.0;
    *dterm = (d == 1) ? 1.0 : 0.0;
    return;
  }

  *term = R_pow_di(lambda, d);
  *dterm = d * R_pow_di(lambda, d - 1);
}

static inline double np_orly_denom_support(const double x,
                                           const double lambda,
                                           const double * const cats,
                                           const int ncat,
                                           const double cl,
                                           const double ch){
  double denom = 0.0;
  int z;

  if(cats != NULL && ncat > 0){
    int i;
    for(i = 0; i < ncat; i++)
      denom += R_pow_di(lambda, (int)fabs(x - cats[i]));
    return denom;
  }

  for(z = (int)cl; z <= (int)ch; z++)
    denom += R_pow_di(lambda, (int)fabs(x - (double)z));

  return denom;
}

static inline double np_orly_kernel_support(const double x,
                                            const double y,
                                            const double lambda,
                                            const double * const cats,
                                            const int ncat,
                                            const double cl,
                                            const double ch){
  const double num = R_pow_di(lambda, (int)fabs(x - y));
  const double den = np_orly_denom_support(x, lambda, cats, ncat, cl, ch);
  return (den > 0.0) ? (num / den) : 0.0;
}

static inline double np_orly_score_support(const double x,
                                           const double y,
                                           const double lambda,
                                           const double * const cats,
                                           const int ncat,
                                           const double cl,
                                           const double ch){
  const int dxy = (int)fabs(x - y);
  double num, dnum, den = 0.0, dden = 0.0;
  int z;

  np_orly_term_deriv(dxy, lambda, &num, &dnum);

  if(cats != NULL && ncat > 0){
    int i;
    for(i = 0; i < ncat; i++){
      const int d = (int)fabs(x - cats[i]);
      double t, dt;
      np_orly_term_deriv(d, lambda, &t, &dt);
      den += t;
      dden += dt;
    }
  } else {
    for(z = (int)cl; z <= (int)ch; z++){
      const int d = (int)fabs(x - (double)z);
      double t, dt;
      np_orly_term_deriv(d, lambda, &t, &dt);
      den += t;
      dden += dt;
    }
  }

  if(!(den > 0.0))
    return 0.0;

  return (dnum * den - num * dden)/(den * den);
}

double np_oracine_li_yan(const double x, const double y, const double lambda, const double cl, const double ch){
  return np_orly_kernel_support(x, y, lambda, NULL, 0, cl, ch);
}

double np_score_oracine_li_yan(const double x, const double y, const double lambda, const double cl, const double ch){
  return np_orly_score_support(x, y, lambda, NULL, 0, cl, ch);
}

static inline double np_ordered_eval_cached012(const int kernel,
                                               const double x,
                                               const double y,
                                               const double lambda,
                                               const int max_cxy,
                                               const double * const lpow){
  const int cxy = (int)fabs(x-y);
  const double gee = (lpow != NULL && cxy <= max_cxy) ? lpow[cxy] : R_pow_di(lambda, cxy);
  switch(kernel){
    case 0:
      return (cxy == 0) ? (1.0-lambda) : (0.5*(1.0-lambda)*gee);
    case 1:
      return gee;
    case 2:
      return gee*(1.0-lambda)/(1.0+lambda);
    default:
      return 0.0;
  }
}

static inline double np_ordered_eval_kernel(const int kernel,
                                            const double x,
                                            const double y,
                                            const double lambda,
                                            const int max_cxy,
                                            const double * const lpow,
                                            const double * const cats,
                                            const int ncat,
                                            const double cl,
                                            const double ch){
  if(kernel >= 0 && kernel <= 2)
    return np_ordered_eval_cached012(kernel, x, y, lambda, max_cxy, lpow);
  if(kernel == 3)
    return np_orly_kernel_support(x, y, lambda, cats, ncat, cl, ch);

  return 0.0;
}

static inline uint64_t np_mix_u64(uint64_t x){
  x ^= x >> 30;
  x *= UINT64_C(0xbf58476d1ce4e5b9);
  x ^= x >> 27;
  x *= UINT64_C(0x94d049bb133111eb);
  x ^= x >> 31;
  return x;
}

static inline uint64_t np_hash_discrete_profile_idx(const int idx,
                                                    const int num_reg_unordered,
                                                    const int num_reg_ordered,
                                                    double * const * const xtu,
                                                    double * const * const xto){
  uint64_t h = UINT64_C(0x9e3779b97f4a7c15);
  uint64_t bits;

  for(int u = 0; u < num_reg_unordered; u++){
    memcpy(&bits, &xtu[u][idx], sizeof(bits));
    h ^= np_mix_u64(bits + UINT64_C(0x9e3779b97f4a7c15) + ((uint64_t)u << 1));
  }

  for(int o = 0; o < num_reg_ordered; o++){
    memcpy(&bits, &xto[o][idx], sizeof(bits));
    h ^= np_mix_u64(bits + UINT64_C(0x517cc1b727220a95) + ((uint64_t)o << 1));
  }

  return h;
}

static inline int np_same_discrete_profile_idx(const int ia,
                                               const int ib,
                                               const int num_reg_unordered,
                                               const int num_reg_ordered,
                                               double * const * const xtu,
                                               double * const * const xto){
  for(int u = 0; u < num_reg_unordered; u++)
    if(xtu[u][ia] != xtu[u][ib]) return 0;

  for(int o = 0; o < num_reg_ordered; o++)
    if(xto[o][ia] != xto[o][ib]) return 0;

  return 1;
}

static int np_build_discrete_profile_index(const int num_xt,
                                           const int num_reg_unordered,
                                           const int num_reg_ordered,
                                           double * const * const xtu,
                                           double * const * const xto,
                                           int **out_disc_prof_id,
                                           int **out_disc_prof_rep,
                                           int *out_disc_nprof){
  int i, pos, pid;
  int tsz = 1;
  int nprof = 0;
  int *disc_prof_id = NULL, *disc_prof_rep = NULL, *htable = NULL;
  uint64_t *disc_prof_hash = NULL;

  if(out_disc_prof_id == NULL || out_disc_prof_rep == NULL || out_disc_nprof == NULL)
    return 0;

  *out_disc_prof_id = NULL;
  *out_disc_prof_rep = NULL;
  *out_disc_nprof = 0;

  if(num_xt <= 0 || (num_reg_unordered + num_reg_ordered) <= 1)
    return 0;

  disc_prof_id = (int *)malloc((size_t)num_xt*sizeof(int));
  disc_prof_rep = (int *)malloc((size_t)num_xt*sizeof(int));
  disc_prof_hash = (uint64_t *)malloc((size_t)num_xt*sizeof(uint64_t));

  if(disc_prof_id == NULL || disc_prof_rep == NULL || disc_prof_hash == NULL)
    goto fail;

  while(tsz < (2*num_xt)) tsz <<= 1;
  htable = (int *)malloc((size_t)tsz*sizeof(int));
  if(htable == NULL)
    goto fail;

  for(i = 0; i < tsz; i++)
    htable[i] = -1;

  for(i = 0; i < num_xt; i++){
    const uint64_t h = np_hash_discrete_profile_idx(i, num_reg_unordered, num_reg_ordered, xtu, xto);
    pos = (int)(h & (uint64_t)(tsz - 1));
    pid = -1;

    while(1){
      const int hs = htable[pos];
      if(hs < 0){
        pid = nprof++;
        htable[pos] = pid;
        disc_prof_rep[pid] = i;
        disc_prof_hash[pid] = h;
        break;
      }
      if(disc_prof_hash[hs] == h &&
         np_same_discrete_profile_idx(i, disc_prof_rep[hs], num_reg_unordered, num_reg_ordered, xtu, xto)){
        pid = hs;
        break;
      }
      pos = (pos + 1) & (tsz - 1);
    }

    disc_prof_id[i] = pid;
  }

  free(htable);
  free(disc_prof_hash);

  *out_disc_prof_id = disc_prof_id;
  *out_disc_prof_rep = disc_prof_rep;
  *out_disc_nprof = nprof;
  return 1;

 fail:
  if(htable != NULL) free(htable);
  if(disc_prof_hash != NULL) free(disc_prof_hash);
  if(disc_prof_id != NULL) free(disc_prof_id);
  if(disc_prof_rep != NULL) free(disc_prof_rep);
  return 0;
}

typedef struct {
  int valid;
  int num_xt;
  int num_reg_unordered;
  int num_reg_ordered;
  int nprof;
  const double **xtu_rows;
  const double **xto_rows;
  int *disc_prof_id;
  int *disc_prof_rep;
} NP_DiscProfileCache;

static NP_DiscProfileCache np_disc_profile_cache = {0};

typedef struct {
  int valid;
  int num_obs_train;
  int num_obs_eval;
  int num_reg_continuous;
  double rel_tol;
  const double **x_train_rows;
  const double **x_eval_rows;
  int *kernel_c;
  int *cont_ok;
  double *cont_hmin;
  double *cont_k0;
} NP_ContLargeHCache;

static NP_ContLargeHCache np_cont_largeh_cache = {0};
static uint64_t np_fastcv_alllarge_hits = 0;
/*
  Keep only true all-large shortcuts enabled.
  Partial per-variable gate shortcuts are disabled by default to avoid
  per-evaluation overhead on the normal kernel path.
*/
static const int np_partial_gate_features_enabled = 0;

void np_fastcv_alllarge_hits_reset(void){
  np_fastcv_alllarge_hits = 0;
}

double np_fastcv_alllarge_hits_get(void){
  return (double)np_fastcv_alllarge_hits;
}

static inline void np_disc_profile_cache_clear(void){
  if(np_disc_profile_cache.xtu_rows != NULL) free((void *)np_disc_profile_cache.xtu_rows);
  if(np_disc_profile_cache.xto_rows != NULL) free((void *)np_disc_profile_cache.xto_rows);
  if(np_disc_profile_cache.disc_prof_id != NULL) free(np_disc_profile_cache.disc_prof_id);
  if(np_disc_profile_cache.disc_prof_rep != NULL) free(np_disc_profile_cache.disc_prof_rep);
  memset(&np_disc_profile_cache, 0, sizeof(np_disc_profile_cache));
}

static inline void np_cont_largeh_cache_clear(void){
  if(np_cont_largeh_cache.x_train_rows != NULL) free((void *)np_cont_largeh_cache.x_train_rows);
  if(np_cont_largeh_cache.x_eval_rows != NULL) free((void *)np_cont_largeh_cache.x_eval_rows);
  if(np_cont_largeh_cache.kernel_c != NULL) free(np_cont_largeh_cache.kernel_c);
  if(np_cont_largeh_cache.cont_ok != NULL) free(np_cont_largeh_cache.cont_ok);
  if(np_cont_largeh_cache.cont_hmin != NULL) free(np_cont_largeh_cache.cont_hmin);
  if(np_cont_largeh_cache.cont_k0 != NULL) free(np_cont_largeh_cache.cont_k0);
  memset(&np_cont_largeh_cache, 0, sizeof(np_cont_largeh_cache));
}

static int np_disc_profile_cache_match(const int num_xt,
                                       const int num_reg_unordered,
                                       const int num_reg_ordered,
                                       double * const * const xtu,
                                       double * const * const xto){
  int i;
  if(!np_disc_profile_cache.valid) return 0;
  if(np_disc_profile_cache.num_xt != num_xt) return 0;
  if(np_disc_profile_cache.num_reg_unordered != num_reg_unordered) return 0;
  if(np_disc_profile_cache.num_reg_ordered != num_reg_ordered) return 0;
  if(num_reg_unordered > 0){
    if(np_disc_profile_cache.xtu_rows == NULL || xtu == NULL) return 0;
    for(i = 0; i < num_reg_unordered; i++)
      if(np_disc_profile_cache.xtu_rows[i] != xtu[i]) return 0;
  }
  if(num_reg_ordered > 0){
    if(np_disc_profile_cache.xto_rows == NULL || xto == NULL) return 0;
    for(i = 0; i < num_reg_ordered; i++)
      if(np_disc_profile_cache.xto_rows[i] != xto[i]) return 0;
  }
  return 1;
}

static int np_disc_profile_cache_get_or_build(const int num_xt,
                                              const int num_reg_unordered,
                                              const int num_reg_ordered,
                                              double * const * const xtu,
                                              double * const * const xto,
                                              int **out_disc_prof_id,
                                              int **out_disc_prof_rep,
                                              int *out_disc_nprof){
  int i;
  int *disc_prof_id = NULL, *disc_prof_rep = NULL;
  const double **xtu_rows = NULL, **xto_rows = NULL;
  int nprof = 0;

  if(out_disc_prof_id == NULL || out_disc_prof_rep == NULL || out_disc_nprof == NULL)
    return 0;

  *out_disc_prof_id = NULL;
  *out_disc_prof_rep = NULL;
  *out_disc_nprof = 0;

  if(np_disc_profile_cache_match(num_xt, num_reg_unordered, num_reg_ordered, xtu, xto)){
    *out_disc_prof_id = np_disc_profile_cache.disc_prof_id;
    *out_disc_prof_rep = np_disc_profile_cache.disc_prof_rep;
    *out_disc_nprof = np_disc_profile_cache.nprof;
    return 1;
  }

  if(!np_build_discrete_profile_index(num_xt,
                                      num_reg_unordered,
                                      num_reg_ordered,
                                      xtu,
                                      xto,
                                      &disc_prof_id,
                                      &disc_prof_rep,
                                      &nprof))
    return 0;

  if(num_reg_unordered > 0){
    xtu_rows = (const double **)malloc((size_t)num_reg_unordered*sizeof(double *));
    if(xtu_rows == NULL) goto fail;
    for(i = 0; i < num_reg_unordered; i++) xtu_rows[i] = xtu[i];
  }

  if(num_reg_ordered > 0){
    xto_rows = (const double **)malloc((size_t)num_reg_ordered*sizeof(double *));
    if(xto_rows == NULL) goto fail;
    for(i = 0; i < num_reg_ordered; i++) xto_rows[i] = xto[i];
  }

  np_disc_profile_cache_clear();
  np_disc_profile_cache.valid = 1;
  np_disc_profile_cache.num_xt = num_xt;
  np_disc_profile_cache.num_reg_unordered = num_reg_unordered;
  np_disc_profile_cache.num_reg_ordered = num_reg_ordered;
  np_disc_profile_cache.nprof = nprof;
  np_disc_profile_cache.xtu_rows = xtu_rows;
  np_disc_profile_cache.xto_rows = xto_rows;
  np_disc_profile_cache.disc_prof_id = disc_prof_id;
  np_disc_profile_cache.disc_prof_rep = disc_prof_rep;

  *out_disc_prof_id = np_disc_profile_cache.disc_prof_id;
  *out_disc_prof_rep = np_disc_profile_cache.disc_prof_rep;
  *out_disc_nprof = np_disc_profile_cache.nprof;
  return 1;

 fail:
  if(xtu_rows != NULL) free((void *)xtu_rows);
  if(xto_rows != NULL) free((void *)xto_rows);
  if(disc_prof_id != NULL) free(disc_prof_id);
  if(disc_prof_rep != NULL) free(disc_prof_rep);
  return 0;
}

/*
  Large-bandwidth shortcut for ordinary continuous kernels:
  if max_i |(x - x_i)/h| <= u_tol, replace K((x-x_i)/h) by K(0) to avoid
  repeated kernel evaluations. This is intentionally conservative and only
  enabled for ordinary continuous kernels (0..9).
*/
static inline int np_cont_largeh_kernel_supported(const int kernel){
  return (kernel >= 0 && kernel <= 9);
}

static inline double np_get_option_double(const char * const name, const double fallback){
  const SEXP sym = Rf_install(name);
  const SEXP val = Rf_GetOption1(sym);

  if(val == R_NilValue)
    return fallback;

  if(TYPEOF(val) == REALSXP && XLENGTH(val) > 0)
    return REAL(val)[0];

  if(TYPEOF(val) == INTSXP && XLENGTH(val) > 0)
    return (double)INTEGER(val)[0];

  if(TYPEOF(val) == LGLSXP && XLENGTH(val) > 0)
    return (double)LOGICAL(val)[0];

  return fallback;
}

/*
  Lean LL/CVLS hot path for fixed bandwidth with one continuous regressor.
  Returns 1 on success and writes *cv_out, else 0 to signal fallback.
*/
static int np_runtime_tol_cache_ready = 0;
static double np_largeh_rel_tol_cache = 1e-3;
static double np_disc_rel_tol_cache = 1e-2;

static inline void np_refresh_runtime_tolerances(void);

static inline double np_cont_largeh_k0(const int kernel){
  switch(kernel){
    case 0: return np_gauss2(0.0);
    case 1: return np_gauss4(0.0);
    case 2: return np_gauss6(0.0);
    case 3: return np_gauss8(0.0);
    case 4: return np_epan2(0.0);
    case 5: return np_epan4(0.0);
    case 6: return np_epan6(0.0);
    case 7: return np_epan8(0.0);
    case 8: return np_rect(0.0);
    case 9: return np_tgauss2(0.0);
    default: return 0.0;
  }
}

static inline double np_cont_largeh_utol(const int kernel, const double rel_tol){
  const double rt = (rel_tol <= 0.0) ? 1e-3 : rel_tol;
  switch(kernel){
    case 0: case 1: case 2: case 3: case 9:
      /* For Gaussian-like kernels, relative deviation near 0 is ~ u^2/2. */
      return sqrt(-2.0*log(1.0-rt));
    case 4: case 5: case 6: case 7:
      /* For Epanechnikov family, relative deviation near 0 is ~ u^2. */
      return sqrt(rt);
    case 8:
      /* Rectangular kernel is exactly constant for |u| < 1. */
      return 1.0 - 32.0*DBL_EPSILON;
    default:
      return 0.0;
  }
}

static inline double np_cont_largeh_rel_tol(void){
  if(!np_runtime_tol_cache_ready)
    np_refresh_runtime_tolerances();
  return np_largeh_rel_tol_cache;
}

static int np_cont_largeh_build_params(const int num_obs_train,
                                       const int num_obs_eval,
                                       const int num_reg_continuous,
                                       const int * const kernel_c,
                                       double * const * const x_train,
                                       double * const * const x_eval,
                                       const double rel_tol,
                                       int **out_ok,
                                       double **out_hmin,
                                       double **out_k0){
  int i, ii;
  int *cont_ok = NULL;
  double *cont_hmin = NULL, *cont_k0 = NULL;

  if((out_ok == NULL) || (out_hmin == NULL) || (out_k0 == NULL))
    return 0;

  *out_ok = NULL;
  *out_hmin = NULL;
  *out_k0 = NULL;

  if(num_reg_continuous <= 0 || kernel_c == NULL || x_train == NULL || x_eval == NULL)
    return 0;

  cont_ok = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
  cont_hmin = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
  cont_k0 = (double *)malloc((size_t)num_reg_continuous*sizeof(double));

  if((cont_ok == NULL) || (cont_hmin == NULL) || (cont_k0 == NULL))
    goto fail;

  for(i = 0; i < num_reg_continuous; i++){
    const int kern = kernel_c[i];
    double xmin = DBL_MAX, xmax = -DBL_MAX;

    cont_ok[i] = 0;
    cont_hmin[i] = DBL_MAX;
    cont_k0[i] = 0.0;

    if(!np_cont_largeh_kernel_supported(kern))
      continue;

    for(ii = 0; ii < num_obs_train; ii++){
      const double v = x_train[i][ii];
      if(!isfinite(v)) continue;
      xmin = MIN(xmin, v);
      xmax = MAX(xmax, v);
    }

    for(ii = 0; ii < num_obs_eval; ii++){
      const double v = x_eval[i][ii];
      if(!isfinite(v)) continue;
      xmin = MIN(xmin, v);
      xmax = MAX(xmax, v);
    }

    if(xmax >= xmin){
      const double utol = np_cont_largeh_utol(kern, rel_tol);
      if(utol > 0.0 && isfinite(utol)){
        cont_ok[i] = 1;
        cont_hmin[i] = (xmax - xmin)/utol;
        cont_k0[i] = np_cont_largeh_k0(kern);
      }
    }
  }

  *out_ok = cont_ok;
  *out_hmin = cont_hmin;
  *out_k0 = cont_k0;
  return 1;

 fail:
  if(cont_ok != NULL) free(cont_ok);
  if(cont_hmin != NULL) free(cont_hmin);
  if(cont_k0 != NULL) free(cont_k0);
  return 0;
}

static int np_cont_largeh_cache_match(const int num_obs_train,
                                      const int num_obs_eval,
                                      const int num_reg_continuous,
                                      const int * const kernel_c,
                                      double * const * const x_train,
                                      double * const * const x_eval,
                                      const double rel_tol){
  int i;
  if(!np_cont_largeh_cache.valid) return 0;
  if(np_cont_largeh_cache.num_obs_train != num_obs_train) return 0;
  if(np_cont_largeh_cache.num_obs_eval != num_obs_eval) return 0;
  if(np_cont_largeh_cache.num_reg_continuous != num_reg_continuous) return 0;
  if(np_cont_largeh_cache.rel_tol != rel_tol) return 0;
  if(num_reg_continuous <= 0) return 0;
  if(np_cont_largeh_cache.kernel_c == NULL || np_cont_largeh_cache.x_train_rows == NULL || np_cont_largeh_cache.x_eval_rows == NULL)
    return 0;
  for(i = 0; i < num_reg_continuous; i++){
    if(np_cont_largeh_cache.kernel_c[i] != kernel_c[i]) return 0;
    if(np_cont_largeh_cache.x_train_rows[i] != x_train[i]) return 0;
    if(np_cont_largeh_cache.x_eval_rows[i] != x_eval[i]) return 0;
  }
  return 1;
}

static int np_cont_largeh_cache_get_or_build(const int num_obs_train,
                                             const int num_obs_eval,
                                             const int num_reg_continuous,
                                             const int * const kernel_c,
                                             double * const * const x_train,
                                             double * const * const x_eval,
                                             const double rel_tol,
                                             int **out_ok,
                                             double **out_hmin,
                                             double **out_k0){
  int i;
  int *cont_ok = NULL, *kernel_copy = NULL;
  double *cont_hmin = NULL, *cont_k0 = NULL;
  const double **x_train_rows = NULL, **x_eval_rows = NULL;

  if((out_ok == NULL) || (out_hmin == NULL) || (out_k0 == NULL))
    return 0;

  *out_ok = NULL;
  *out_hmin = NULL;
  *out_k0 = NULL;

  if(np_cont_largeh_cache_match(num_obs_train, num_obs_eval, num_reg_continuous, kernel_c, x_train, x_eval, rel_tol)){
    *out_ok = np_cont_largeh_cache.cont_ok;
    *out_hmin = np_cont_largeh_cache.cont_hmin;
    *out_k0 = np_cont_largeh_cache.cont_k0;
    return 1;
  }

  if(!np_cont_largeh_build_params(num_obs_train, num_obs_eval, num_reg_continuous, kernel_c, x_train, x_eval, rel_tol, &cont_ok, &cont_hmin, &cont_k0))
    return 0;

  x_train_rows = (const double **)malloc((size_t)num_reg_continuous*sizeof(double *));
  x_eval_rows = (const double **)malloc((size_t)num_reg_continuous*sizeof(double *));
  kernel_copy = (int *)malloc((size_t)num_reg_continuous*sizeof(int));
  if((x_train_rows == NULL) || (x_eval_rows == NULL) || (kernel_copy == NULL))
    goto fail;

  for(i = 0; i < num_reg_continuous; i++){
    x_train_rows[i] = x_train[i];
    x_eval_rows[i] = x_eval[i];
    kernel_copy[i] = kernel_c[i];
  }

  np_cont_largeh_cache_clear();
  np_cont_largeh_cache.valid = 1;
  np_cont_largeh_cache.num_obs_train = num_obs_train;
  np_cont_largeh_cache.num_obs_eval = num_obs_eval;
  np_cont_largeh_cache.num_reg_continuous = num_reg_continuous;
  np_cont_largeh_cache.rel_tol = rel_tol;
  np_cont_largeh_cache.x_train_rows = x_train_rows;
  np_cont_largeh_cache.x_eval_rows = x_eval_rows;
  np_cont_largeh_cache.kernel_c = kernel_copy;
  np_cont_largeh_cache.cont_ok = cont_ok;
  np_cont_largeh_cache.cont_hmin = cont_hmin;
  np_cont_largeh_cache.cont_k0 = cont_k0;

  *out_ok = np_cont_largeh_cache.cont_ok;
  *out_hmin = np_cont_largeh_cache.cont_hmin;
  *out_k0 = np_cont_largeh_cache.cont_k0;
  return 1;

 fail:
  if(x_train_rows != NULL) free((void *)x_train_rows);
  if(x_eval_rows != NULL) free((void *)x_eval_rows);
  if(kernel_copy != NULL) free(kernel_copy);
  if(cont_ok != NULL) free(cont_ok);
  if(cont_hmin != NULL) free(cont_hmin);
  if(cont_k0 != NULL) free(cont_k0);
  return 0;
}

static inline int np_cont_largeh_is_active(const int kernel,
                                           const double * const xt,
                                           const int num_xt,
                                           const double x,
                                           const double h,
                                           const XL * const xl){
  if((h == 0.0) || !isfinite(h) || !np_cont_largeh_kernel_supported(kernel))
    return 0;

  const double utol = np_cont_largeh_utol(kernel, np_cont_largeh_rel_tol());
  if(!(utol > 0.0) || !isfinite(utol))
    return 0;

  const double maxabs = utol * fabs(h);

  if(xl == NULL){
    for(int i = 0; i < num_xt; i++)
      if(fabs(x - xt[i]) > maxabs) return 0;
  } else {
    for(int m = 0; m < xl->n; m++){
      const int istart = xl->istart[m];
      const int nlev = xl->nlev[m];
      for(int i = istart; i < istart+nlev; i++)
        if(fabs(x - xt[i]) > maxabs) return 0;
    }
  }

  return 1;
}

/*
  Discrete "near upper-bound lambda" shortcut:
  if unordered same/different kernel values are nearly identical, replace
  per-observation category checks by a constant multiply.
*/
static inline double np_disc_rel_tol(void){
  if(!np_runtime_tol_cache_ready)
    np_refresh_runtime_tolerances();
  return np_disc_rel_tol_cache;
}

static inline void np_refresh_runtime_tolerances(void){
  const double largeh_dflt = 1e-3;
  const double largeh_optv = np_get_option_double("np.largeh.rel.tol", largeh_dflt);
  if(isfinite(largeh_optv) && largeh_optv > 0.0 && largeh_optv < 0.1) {
    np_largeh_rel_tol_cache = largeh_optv;
  } else {
    np_largeh_rel_tol_cache = largeh_dflt;
    /* fallback for legacy/developer workflows */
    {
      const char *rt_env = getenv("NP_LARGEH_REL_TOL");
      if(rt_env != NULL && rt_env[0] != '\0'){
        const double v = atof(rt_env);
        if(isfinite(v) && v > 0.0 && v < 0.1)
          np_largeh_rel_tol_cache = v;
      }
    }
  }

  {
    const double disc_dflt = 1e-2;
    const double disc_optv = np_get_option_double("np.disc.upper.rel.tol", disc_dflt);
    if(isfinite(disc_optv) && disc_optv > 0.0 && disc_optv < 0.5) {
      np_disc_rel_tol_cache = disc_optv;
    } else {
      np_disc_rel_tol_cache = disc_dflt;
      /* fallback for legacy/developer workflows */
      {
        const char *rt_env = getenv("NP_DISC_UPPER_REL_TOL");
        if(rt_env != NULL && rt_env[0] != '\0'){
          const double v = atof(rt_env);
          if(isfinite(v) && v > 0.0 && v < 0.5)
            np_disc_rel_tol_cache = v;
        }
      }
    }
  }

  np_runtime_tol_cache_ready = 1;
}

static inline int np_disc_unordered_has_upper(const int kernel){
  return (kernel >= 0 && kernel <= 5);
}

static inline double np_disc_unordered_upper(const int kernel, const int ncat){
  const double c = (double)ncat;
  switch(kernel){
    case 0: /* AA */
    case 2: /* econvol AA */
    case 4: /* score AA */
      return (c > 0.0) ? ((c - 1.0)/c) : 0.0;
    case 1: /* Li-Racine */
    case 3: /* econvol Li-Racine */
    case 5: /* score Li-Racine */
      return 1.0;
    default:
      return 1.0;
  }
}

static inline int np_disc_near_upper(const int kernel, const double lambda, const int ncat){
  if(!isfinite(lambda) || !np_disc_unordered_has_upper(kernel))
    return 0;
  const double up = np_disc_unordered_upper(kernel, ncat);
  double tol = np_disc_rel_tol()*fabs(up);
  /* Keep a tiny absolute floor for numerical robustness near machine precision. */
  tol = fmax(tol, 16.0*DBL_EPSILON*fmax(1.0, fabs(up)));
  return fabs(lambda - up) <= tol;
}

static inline int np_disc_near_const_kernel(const double k_same, const double k_diff){
  const double scale = fmax(1.0, fmax(fabs(k_same), fabs(k_diff)));
  return fabs(k_same - k_diff) <= np_disc_rel_tol()*scale;
}

/*
  Gate override behavior for kernel_weighted_sum_np_ctx:
  - `NP_GATE_CTX_INACTIVE`: no override; build gate/caches locally.
  - `NP_GATE_CTX_OVERRIDE`: caller supplies precomputed gate state.
  - `NP_GATE_CTX_DISABLE`: bypass gate shortcuts and run generic path.
*/
typedef struct {
  int active;
  int num_reg_continuous;
  int num_reg_unordered;
  int num_reg_ordered;
  const int *kernel_c;
  const int *kernel_u;
  const int *kernel_o;
  const int *operator;
  const int *cont_ok;
  const double *cont_hmin;
  const double *cont_k0;
  const int *disc_uno_ok;
  const double *disc_uno_const;
  const int *disc_ord_ok;
  const double *disc_ord_const;
  int disc_prof_num_xt;
  int disc_nprof;
  const int *disc_prof_id;
  const int *disc_prof_rep;
} NP_GateOverrideCtx;

#define NP_GATE_CTX_INACTIVE 0
#define NP_GATE_CTX_OVERRIDE 1
#define NP_GATE_CTX_DISABLE 2

static NP_GateOverrideCtx np_gate_override_ctx = {0};

static inline void np_gate_ctx_clear(NP_GateOverrideCtx * const ctx){
  if(ctx != NULL)
    memset(ctx, 0, sizeof(*ctx));
}

static inline void np_gate_override_clear(void){
  memset(&np_gate_override_ctx, 0, sizeof(np_gate_override_ctx));
  if(!np_partial_gate_features_enabled)
    np_gate_override_ctx.active = NP_GATE_CTX_DISABLE;
}

static inline int np_gate_int_array_equal(const int * const a,
                                          const int * const b,
                                          const int n){
  int i;
  if(n <= 0)
    return 1;
  if(a == b)
    return 1;
  if(a == NULL || b == NULL)
    return 0;
  for(i = 0; i < n; i++)
    if(a[i] != b[i])
      return 0;
  return 1;
}

static inline int np_gate_ptr_pair_valid(const void * const a,
                                         const void * const b){
  return ((a == NULL) == (b == NULL));
}

static inline int np_gate_ptr_triplet_valid(const void * const a,
                                            const void * const b,
                                            const void * const c){
  return ((a == NULL) && (b == NULL) && (c == NULL)) ||
         ((a != NULL) && (b != NULL) && (c != NULL));
}

static inline int np_gate_ctx_is_sane(const NP_GateOverrideCtx * const ctx){
  if(ctx == NULL)
    return 1;

  if((ctx->active == NP_GATE_CTX_INACTIVE) || (ctx->active == NP_GATE_CTX_DISABLE))
    return 1;

  if(ctx->active != NP_GATE_CTX_OVERRIDE)
    return 0;

  if((ctx->num_reg_continuous < 0) ||
     (ctx->num_reg_unordered < 0) ||
     (ctx->num_reg_ordered < 0) ||
     (ctx->disc_nprof < 0))
    return 0;

  if(!np_gate_ptr_triplet_valid(ctx->cont_ok, ctx->cont_hmin, ctx->cont_k0))
    return 0;
  if(!np_gate_ptr_pair_valid(ctx->disc_uno_ok, ctx->disc_uno_const))
    return 0;
  if(!np_gate_ptr_pair_valid(ctx->disc_ord_ok, ctx->disc_ord_const))
    return 0;
  if(!np_gate_ptr_pair_valid(ctx->disc_prof_id, ctx->disc_prof_rep))
    return 0;

  return 1;
}

static inline int np_gate_ctx_signature_matches(const NP_GateOverrideCtx * const ctx,
                                                const int num_reg_continuous,
                                                const int num_reg_unordered,
                                                const int num_reg_ordered,
                                                const int * const kernel_c,
                                                const int * const kernel_u,
                                                const int * const kernel_o,
                                                const int * const operator){
  const int total_regs = num_reg_continuous + num_reg_unordered + num_reg_ordered;

  if((ctx == NULL) ||
     (ctx->num_reg_continuous != num_reg_continuous) ||
     (ctx->num_reg_unordered != num_reg_unordered) ||
     (ctx->num_reg_ordered != num_reg_ordered))
    return 0;

  if(!np_gate_int_array_equal(ctx->operator, operator, total_regs))
    return 0;
  if(!np_gate_int_array_equal(ctx->kernel_c, kernel_c, num_reg_continuous))
    return 0;
  if(!np_gate_int_array_equal(ctx->kernel_u, kernel_u, num_reg_unordered))
    return 0;
  if(!np_gate_int_array_equal(ctx->kernel_o, kernel_o, num_reg_ordered))
    return 0;

  return 1;
}

static inline void np_gate_ctx_set(NP_GateOverrideCtx * const ctx,
                                   const int num_reg_continuous,
                                   const int num_reg_unordered,
                                   const int num_reg_ordered,
                                   const int * const kernel_c,
                                   const int * const kernel_u,
                                   const int * const kernel_o,
                                   const int * const operator,
                                   const int * const cont_ok,
                                   const double * const cont_hmin,
                                   const double * const cont_k0,
                                   const int * const disc_uno_ok,
                                   const double * const disc_uno_const,
                                   const int * const disc_ord_ok,
                                   const double * const disc_ord_const){
  if(ctx == NULL) return;
  ctx->active = NP_GATE_CTX_OVERRIDE;
  ctx->num_reg_continuous = num_reg_continuous;
  ctx->num_reg_unordered = num_reg_unordered;
  ctx->num_reg_ordered = num_reg_ordered;
  ctx->kernel_c = kernel_c;
  ctx->kernel_u = kernel_u;
  ctx->kernel_o = kernel_o;
  ctx->operator = operator;
  ctx->cont_ok = cont_ok;
  ctx->cont_hmin = cont_hmin;
  ctx->cont_k0 = cont_k0;
  ctx->disc_uno_ok = disc_uno_ok;
  ctx->disc_uno_const = disc_uno_const;
  ctx->disc_ord_ok = disc_ord_ok;
  ctx->disc_ord_const = disc_ord_const;
}

static inline int np_disc_ordered_has_upper(const int kernel){
  /* Exclude WvR family; support Li-Racine and transformed variants. */
  switch(kernel){
    case 1: case 2: case 3:
    case 6: case 7:
    case 9: case 10: case 11:
    case 13: case 14: case 15:
      return 1;
    default:
      return 0;
  }
}

static inline int np_disc_ordered_near_upper(const int kernel, const double lambda){
  if(!isfinite(lambda) || !np_disc_ordered_has_upper(kernel))
    return 0;
  const double up = 1.0;
  double tol = np_disc_rel_tol()*fabs(up);
  tol = fmax(tol, 16.0*DBL_EPSILON*fmax(1.0, fabs(up)));
  return fabs(lambda - up) <= tol;
}

static inline void np_ckernelv_mul_const(const double c,
                                         const int num_xt,
                                         const int do_xw,
                                         double * const result,
                                         const XL * const xl){
  const int bin_do_xw = do_xw > 0;
  int i;

  if(xl == NULL){
    if(!bin_do_xw){
      for(i = 0; i < num_xt; i++)
        result[i] = c;
    } else {
      for(i = 0; i < num_xt; i++){
        if(result[i] == 0.0) continue;
        result[i] *= c;
      }
    }
  } else {
    if(!bin_do_xw){
      for(int m = 0; m < xl->n; m++){
        const int istart = xl->istart[m];
        const int nlev = xl->nlev[m];
        for(i = istart; i < istart+nlev; i++)
          result[i] = c;
      }
    } else {
      for(int m = 0; m < xl->n; m++){
        const int istart = xl->istart[m];
        const int nlev = xl->nlev[m];
        for(i = istart; i < istart+nlev; i++){
          if(result[i] == 0.0) continue;
          result[i] *= c;
        }
      }
    }
  }
}

// not so simple truncated gaussian convolution kernels
//   In general for our truncated Gaussian kernel the convolution kernel will be a polynomial of the form:
// z < 0: 
// a0*erf*(z/2 + b)*exp(-z^2/4) + a1*z +a2*erf(z/sqrt(2) + b/sqrt(2)) + a3
// z > 0
// -a0*erf*(z/2 - b)*exp(-z^2/4) - a1*z - a2*erf(z/sqrt(2) - b/sqrt(2)) + a3
double np_econvol_rect(const double z){
  return ((fabs(z) < 2.0) ? 0.25 : 0.0);
}

double np_econvol_tgauss2(const double z){
  const double az = fabs(z);
  if(az >= 2*np_tgauss2_b)
    return 0.0;
  else {
    return(-np_tgauss2_a0*erfun(0.5*az - np_tgauss2_b)*exp(-0.25*az*az) - np_tgauss2_a1*az -
           np_tgauss2_a2*erfun(0.7071067810*(az - np_tgauss2_b)) - np_tgauss2_c0);

  }
}

// the simple convolution kernels
double np_econvol_gauss2(const double z){
  return(0.28209479177387814348*exp(-0.25*z*z));
}

double np_econvol_gauss4(const double z){
  const double z2 = z*z;
  return(0.0044077311214668459918*exp(-0.25*z2)*(108.0+z2*(z2-28.0)));
}

double np_econvol_gauss6(const double z){
  const double z2 = z*z;
  return(0.00001721769969*exp(-0.25*z2)*(36240.0+z2*(-19360.0+z2*(2312.0+z2*(-88.0+z2)))));
}

double np_econvol_gauss8(const double z){
  const double z2 = z*z;
  return(0.2989183974E-7*exp(-0.25*z2)*(25018560.0+z2*(-20462400.0+z2*(4202352.0+z2*(-331680.0+z2*(11604.0+z2*(-180.0+z2)))))));
}

// These kernels are in error jracine 27 5 2009
// need to edit both kernel.c and this file
// note we require an additional test for negative (z) as well 

double np_econvol_epan2(const double z){
  const double z2 = z*z;
  return((z2 < 20.0) ? 
         ((z < 0.0) ? 
          (5.579734404642339E-9*(26883*ipow(z,5)-2688300*ipow(z,3)-12022443*z2+48089773)) : 
          (-5.579734404642339E-9*(26883*ipow(z,5)-2688300*ipow(z,3)+12022443*z2-48089773)))
         : 0.0);

}

double np_econvol_epan4(const double z){
  const double z2 = z*z;
  return((z2 < 20.0) ? 
         ((z < 0.0) ? 
          (3.756009615384615e-9*(1456*ipow(z,9)-124800*ipow(z,7)+5491200*ipow(z,5)+15627432*ipow(z,4)-24960000*ipow(z,3)-111624513*z2+148832684))           :
          (-3.756009615384615e-9*(1456*ipow(z,9)-124800*ipow(z,7)+5491200*ipow(z,5)-15627432*ipow(z,4)-24960000*ipow(z,3)+111624513*z2-148832684)))
         : 0.0);
}

double np_econvol_epan6(const double z){
  const double z2 = z*z;
  return((z2 < 20.0) ? 
         ((z < 0.0) ? 
          (9.390024038461537E-11*(2079*ipow(z,13)-206388*ipow(z,11)+8867040*ipow(z,9)-255528000*ipow(z,7)-515705252*ipow(z,6)+1681680000*ipow(z,5)+4922641042*ipow(z,4)-3057600000*ipow(z,3)-13674002896*z2+9015826085)) :
          (-9.390024038461537E-11*(2079*ipow(z,13)-206388*ipow(z,11)+8867040*ipow(z,9)-255528000*ipow(z,7)+515705252*ipow(z,6)+1681680000*ipow(z,5)-4922641042*ipow(z,4)-3057600000*ipow(z,3)+13674002896*z2-9015826085)))
         : 0.0);
}

double np_econvol_epan8(const double z){
  const double z2 = z*z;
  return((z2 < 20.0) ? 
         ((z < 0.0) ? 
          (1.121969784007353E-13*(63063*ipow(z,17)-7351344*ipow(z,15)+373222080*ipow(z,13)-11040382080*ipow(z,11)+241727270400*ipow(z,9)+350679571413*ipow(z,8)-1900039680000*ipow(z,7)-4208154856956*ipow(z,6)+5757696000000*ipow(z,5)+16994471537707*ipow(z,4)-5757696000000*ipow(z,3)-25749199299557*z2+10097725215512)) :
          (-1.121969784007353E-13*(63063*ipow(z,17)-7351344*ipow(z,15)+373222080*ipow(z,13)-11040382080*ipow(z,11)+241727270400*ipow(z,9)-350679571413*ipow(z,8)-1900039680000*ipow(z,7)+4208154856956*ipow(z,6)+5757696000000*ipow(z,5)-16994471537707*ipow(z,4)-5757696000000*ipow(z,3)+25749199299557*z2-10097725215512)))
         : 0.0);
}

double np_econvol_uaa(const int same_cat, const double lambda, const int c){
  if(c < 2)
    return same_cat ? 1.0 : 0.0;

  const double dcm1 = (double)(c-1);
  const double oml = 1.0-lambda;
  return (same_cat)?(oml*oml+lambda*lambda/dcm1):(lambda/dcm1*(2.0*oml+((double)(c-2))*lambda/dcm1));
}

double np_econvol_uli_racine(const int same_cat, const double lambda, const int c){
  return (same_cat)?(1.0 + (double)(c-1)*lambda*lambda):(lambda*(2.0+(double)(c-2)*lambda));
}

double np_econvol_unli_racine(const int same_cat, const double lambda, const int c){
  const double inorm = 1.0/((c-1.0)*lambda + 1.0);
  return ((same_cat)?(1.0 + (double)(c-1)*lambda*lambda):(lambda*(2.0+(double)(c-2)*lambda)))*inorm*inorm;
}

double np_econvol_onli_racine(const double x, const double y, const double lambda, const double cl, const double ch){
  if(lambda == 1.0)
    return 0.0;

  const int cxy = (int)fabs(x-y);
  const double lnorm = (1.0 - lambda)/(1.0 + lambda);
  const double l2 = lambda*lambda;
  return lnorm*lnorm*R_pow_di(lambda, cxy)*((1.0 + l2)/(1.0 - l2) + cxy);

}

double np_econvol_oracine_li_yan(const double x, const double y, const double lambda, const double cl, const double ch){
  double out = 0.0;
  int z;
  const double denx = np_orly_denom_support(x, lambda, NULL, 0, cl, ch);
  const double deny = np_orly_denom_support(y, lambda, NULL, 0, cl, ch);
  const double den = denx*deny;

  if(!(den > 0.0))
    return 0.0;

  for(z = (int)cl; z <= (int)ch; z++){
    const double zz = (double)z;
    out += R_pow_di(lambda, (int)fabs(x-zz)) * R_pow_di(lambda, (int)fabs(y-zz));
  }

  return out/den;
}

double np_econvol_owang_van_ryzin(const double x, const double y, const double lambda, const double cl, const double ch){
  if(lambda == 1.0)
    return 0.0;

  if(x == y) return 0.5*(1.0-lambda)*(1.0-lambda)*(1.0 + 1.0/(1.0-lambda*lambda));

  const int cxy = (int)fabs(x-y);
  const double lnorm = 0.5*(1.0 - lambda);
  const double l2 = lambda*lambda;
  return lnorm*lnorm*R_pow_di(lambda, cxy)*(1.0 + cxy + 2.0/(1.0-l2));
}
// derivative kernels


double np_deriv_tgauss2(const double z){
  return (fabs(z) >= np_tgauss2_b) ? 0.0 : np_tgauss2_alpha*(-z*ONE_OVER_SQRT_TWO_PI*exp(-0.5*z*z));
}


double np_deriv_gauss2(const double z){
  return (-z*ONE_OVER_SQRT_TWO_PI*exp(-0.5*z*z));
}

double np_deriv_gauss4(const double z){
  const double z2 = z*z;
  return (-ONE_OVER_SQRT_TWO_PI*z*(2.5-0.5*z2)*exp(-0.5*z2));
}

double np_deriv_gauss6(const double z){
  const double z2 = z*z;
  return (-0.049867785050179084743*z*exp(-0.5*z2)*(35.0+z2*(-14.0+z2)));
}

double np_deriv_gauss8(const double z){
  const double z2 = z*z;
  return (-ONE_OVER_SQRT_TWO_PI*z*(6.5625+z2*(-3.9375+z2*(0.5625-0.02083333333*z2)))*exp(-0.5*z2));
}

double np_deriv_epan2(const double z){
  return (z*z < 5.0)?(-0.13416407864998738178*z):0.0;
}

double np_deriv_epan4(const double z){
  const double z2 = z*z;
  return (z2 < 5.0)?(z*(2.347871374742824e-1*z2-8.385254921942804e-1)):0.0;
}

double np_deriv_epan6(const double z){
  const double z2 = z*z;
  return (z2 < 5.0)?(z*(-2.567984320334919+z2*(1.848948710641142-2.905490831007508e-1*z2))):0.0;
}

double np_deriv_epan8(const double z){
  const double z2 = z*z;
  return (z2 < 5.0)?(z*(-5.777964720753567+z2*(7.626913431394709+z2*(-2.83285356023232+3.147615066924801e-1*z2)))):0.0;
}

double np_deriv_rect(const double z){
  return (0.0);
}

// cdf kernels

double np_cdf_tgauss2(const double z){
  return (z <= -np_tgauss2_b) ? 0.0 : ((z >= np_tgauss2_b) ? 1.0 : (np_tgauss2_alpha*0.5*erfun(0.7071067810*z)-np_tgauss2_c0*z + 0.5));
}


double np_cdf_gauss2(const double z){
  return (0.5*erfun(0.7071067810*z)+0.5);
}

double np_cdf_gauss4(const double z){
  return (0.5*erfun(0.7071067810*z)+0.1994711401*z*exp(-0.5*z*z)+0.5);
}

double np_cdf_gauss6(const double z){
  const double z2 = z*z;
  return (0.5*erfun(0.7071067810*z)+z*exp(-0.5*z2)*
          (0.3490744952-0.04986778504*z2)+0.5);
}

double np_cdf_gauss8(const double z){
  const double z2 = z*z;
  return (0.5*erfun(0.7071067810*z)+z*exp(-0.5*z2)*
          (0.4737439578+z2*(-0.1329807601+0.008311297511*z2)) + 0.5);
}

double np_cdf_epan2(const double z){
  return (z < -SQRT_5) ? 0.0 : (z > SQRT_5) ? 1.0 : 
    (z*(0.3354101967-0.02236067978*z*z)+0.5);
}

double np_cdf_epan4(const double z){
  const double z2 = z*z;
  return (z < -SQRT_5) ? 0.0 : (z > SQRT_5) ? 1.0 : 
    (0.5+z*(0.6288941188+z2*(-0.1397542486+0.01173935688*z2)));
}

double np_cdf_epan6(const double z){
  const double z2 = z*z;
  return (z < -SQRT_5) ? 0.0 : (z > SQRT_5) ? 1.0 : 
    (0.5+z*(0.9171372566+z2*(-0.4279973864+z2*(0.09244743547-0.006917835307*z2))));
}

double np_cdf_epan8(const double z){
  const double z2 = z*z;
  return (z < -SQRT_5) ? 0.0 : (z > SQRT_5) ? 1.0 : 
    (0.5+z*(1.203742649+z2*(-0.9629941194+z2*(0.3813456714+z2*(-0.06744889424+0.004371687590*z2)))));
}

double np_cdf_rect(const double z){
  return (z < -1.0) ? 0.0 : (z > 1.0) ? 1.0 : (0.5+0.5*z);
}

static inline void np_activate_bounds_x(void){
  int_cker_bound_extern = int_cxker_bound_extern;
  vector_ckerlb_extern = vector_cxkerlb_extern;
  vector_ckerub_extern = vector_cxkerub_extern;
}

static inline void np_activate_bounds_y(void){
  int_cker_bound_extern = int_cyker_bound_extern;
  vector_ckerlb_extern = vector_cykerlb_extern;
  vector_ckerub_extern = vector_cykerub_extern;
}

static inline void np_activate_bounds_xy(void){
  int_cker_bound_extern = int_cxyker_bound_extern;
  vector_ckerlb_extern = vector_cxykerlb_extern;
  vector_ckerub_extern = vector_cxykerub_extern;
}

static inline double np_cker_invnorm(const int kernel,
                                     const double x,
                                     const double h,
                                     const double lb,
                                     const double ub){
  static double (* const cdfk[])(double) = {
    np_cdf_gauss2, np_cdf_gauss4, np_cdf_gauss6, np_cdf_gauss8,
    np_cdf_epan2, np_cdf_epan4, np_cdf_epan6, np_cdf_epan8,
    np_cdf_rect, np_cdf_tgauss2
  };
  static double (* const kbase[])(double) = {
    np_gauss2, np_gauss4, np_gauss6, np_gauss8,
    np_epan2, np_epan4, np_epan6, np_epan8,
    np_rect, np_tgauss2
  };
  int k0 = kernel % 10;
  double zu, zl, den, du, zmid;

  if(!(h > 0.0) || !isfinite(h))
    return 1.0;
  if(!isfinite(lb) && !isfinite(ub))
    return 1.0;

  if(k0 < 0 || k0 > 9)
    k0 = 0;

  zu = isfinite(ub) ? ((ub - x)/h) : R_PosInf;
  zl = isfinite(lb) ? ((lb - x)/h) : R_NegInf;
  du = zu - zl;
  zmid = 0.5*(zu + zl);

  if(isfinite(du) && fabs(du) < 1.0e-5) {
    den = kbase[k0](zmid)*du;
  } else {
    den = cdfk[k0](zu) - cdfk[k0](zl);
  }

  if(!(den > DBL_MIN) && isfinite(du) && (du > 0.0)) {
    den = fmax(kbase[k0](zmid)*du, kbase[k0](0.0)*du);
  }

  den = NZD_POS(den);

  return 1.0/den;
}

static inline double np_cker_base_eval(const int kernel,
                                       const double z){
  static double (* const kbase[])(double) = {
    np_gauss2, np_gauss4, np_gauss6, np_gauss8,
    np_epan2, np_epan4, np_epan6, np_epan8,
    np_rect, np_tgauss2
  };
  int k0 = kernel % 10;

  if(k0 < 0 || k0 > 9)
    k0 = 0;

  return kbase[k0](z);
}

double np_cdf_owang_van_ryzin(const double y, const double x, const double lambda, const double cl, const double ch){
  if(x == y) return 1.0 - 0.5*lambda;
  const int cxy = (int)fabs(x-y);
  const double gee = R_pow_di(lambda, cxy);
  return (x < y) ? 0.5*gee : (1.0 - gee);
}

static inline double np_geom_sum_nonneg_lambda(const int n, const double lambda){
  if(n < 0)
    return 0.0;
  if(lambda == 1.0)
    return (double)(n + 1);
  return (1.0 - R_pow_di(lambda, n + 1))/(1.0 - lambda);
}

double np_cdf_oli_racine(const double y, const double x, const double lambda, const double cl, const double ch){
  const int xh = (x > ch) ? (int)ch : (int)x;
  const int yi = (int)y;
  const int cli = (int)cl;
  const int cxy = abs(xh - yi);

  if(x < y){
    const int nx = (int)x - cli;
    if(x < cl)
      return 0.0;
    return R_pow_di(lambda, cxy) * np_geom_sum_nonneg_lambda(nx, lambda);
  } else{
    const int n1 = yi - cli;
    const int n2 = xh - yi;
    return np_geom_sum_nonneg_lambda(n1, lambda) + np_geom_sum_nonneg_lambda(n2, lambda) - 1.0;
  }
}

double np_cdf_onli_racine(const double y, const double x, const double lambda, const double cl, const double ch){
  const int cxy = (int)fabs(x-y);
  const double gee = R_pow_di(lambda, cxy)/(1.0+lambda);
  return (x < y) ? gee : 1.0 - lambda*gee;
}

double np_cdf_oracine_li_yan(const double y, const double x, const double lambda, const double cl, const double ch){
  double out = 0.0;
  int z;
  const int xh = (x > ch) ? (int)ch : (int)x;
  const double den = np_orly_denom_support(y, lambda, NULL, 0, cl, ch);

  if(x < cl || !(den > 0.0))
    return 0.0;

  for(z = (int)cl; z <= xh; z++)
    out += R_pow_di(lambda, (int)fabs(y - (double)z));

  return out/den;
}

// this is a null kernel, it is a placeholder kernel for testing
double np_onull(const double x, const double y, const double lambda, const double cl, const double ch){
  return(0.0);
}

// adaptive convolution kernels
double np_aconvol_gauss2(const double x, const double y,const double hx,const double hy){
  const double h2 = hx*hx + hy*hy;
  const double z2 = (x-y)*(x-y)/h2;

  return(0.3989422803*hx*hy*exp(-0.5*z2)/sqrt(h2));
}

double np_aconvol_gauss4(const double x, const double y,const double hx,const double hy){
  const double hx2 = hx*hx;
  const double hy2 = hy*hy;
  const double hxy2 = hx2+hy2;
  const double x2 = x*x;
  const double y2 = y*y;
  const double hx3 = hx2*hx;
  const double hy3 = hy2*hy;
  const double hx5 = hx3*hx2;
  const double hy5 = hy3*hy2;
  const double hx7 = hx5*hx2;
  const double hy7 = hy5*hy2;
  const double hx9 = hx7*hx2;
  const double hy9 = hy7*hy2;
  
  return((hx3*hy3*(y2*y2 - 4*x*y*y2 + x2*x2)
                     + (6*hx3*hy3*x2 - 2*hx*hy7 - 6*hx3*hy5 - 12*hx5*hy3 - 2*hx7*hy)*y2
                     + ((4*hx*hy7 + 24*hx3*hy5 + 24*hx5*hy3 + 4*hx7*hy)*x - 4*hx3*hy3*x2*x)*y
                     + ( - 2*hx*hy7 - 12*hx3*hy5 - 12*hx5*hy3 - 2*hx7*hy)*x2
                     + 6*hx*hy9 + 27*hx3*hy7 + 42*hx5*hy5 + 27*hx7*hy3 + 6*hx9*hy)*
         exp(-0.5*(x-y)*(x-y)/hxy2)*ONE_OVER_SQRT_TWO_PI/(sqrt(hy2 + hx2)*4*hxy2*hxy2*hxy2*hxy2));
}

double np_aconvol_gauss6(const double x, const double y,const double hx,const double hy){
  const double hx2 = hx*hx;
  const double hx4 = hx2*hx2;
  const double hx6 = hx4*hx2;
  const double hx8 = hx4*hx4;
  const double hx10 = hx8*hx2;
  const double hx12 = hx10*hx2;
  const double hx14 = hx12*hx2;
  const double hx16 = hx8*hx8;

  const double x2 = x*x;
  const double x3 = x*x2;
  const double x4 = x2*x2;
  const double x5 = x*x4;
  const double x6 = x3*x3;
  const double x7 = x*x6;
  const double x8 = x4*x4;

  
  const double hy2 = hy*hy;
  const double hy4 = hy2*hy2;
  const double hy6 = hy4*hy2;
  const double hy8 = hy4*hy4;
  const double hy10 = hy8*hy2;
  const double hy12 = hy10*hy2;
  const double hy14 = hy12*hy2;
  const double hy16 = hy8*hy8;

  const double y2 = y*y;
  const double y3 = y*y2;
  const double y4 = y2*y2;
  const double y5 = y*y4;
  const double y6 = y3*y3;
  const double y7 = y*y6;
  const double y8 = y4*y4;
  
  const double hxy2 = hx2+hy2;
  const double hxy4 = hxy2*hxy2;
  const double hxy8 = hxy4*hxy4;
    
  return(hx*hy*(hx4*hy4*y8-8*hx4*hy4*x*y7+28*hx4*hy4*x2*y6-4*hx2*hy8*y6
                -40*hx4*hy6*y6-40*hx6*hy4*y6-4*hx8*hy2*y6
                -56*hx4*hy4*x3*y5+24*hx2*hy8*x*y5+240*hx4*hy6*x*y5
                +240*hx6*hy4*x*y5+24*hx8*hy2*x*y5+70*hx4*hy4*x4*y4
                -60*hx2*hy8*x2*y4-600*hx4*hy6*x2*y4
                -600*hx6*hy4*x2*y4-60*hx8*hy2*x2*y4+8*hy12*y4
                +108*hx2*hy10*y4+570*hx4*hy8*y4+940*hx6*hy6*y4
                +570*hx8*hy4*y4+108*hx10*hy2*y4+8*hx12*y4
                -56*hx4*hy4*x5*y3+80*hx2*hy8*x3*y3
                +800*hx4*hy6*x3*y3+800*hx6*hy4*x3*y3
                +80*hx8*hy2*x3*y3-32*hy12*x*y3-432*hx2*hy10*x*y3
                -2280*hx4*hy8*x*y3-3760*hx6*hy6*x*y3
                -2280*hx8*hy4*x*y3-432*hx10*hy2*x*y3-32*hx12*x*y3
                +28*hx4*hy4*x6*y2-60*hx2*hy8*x4*y2
                -600*hx4*hy6*x4*y2-600*hx6*hy4*x4*y2
                -60*hx8*hy2*x4*y2+48*hy12*x2*y2+648*hx2*hy10*x2*y2
                +3420*hx4*hy8*x2*y2+5640*hx6*hy6*x2*y2
                +3420*hx8*hy4*x2*y2+648*hx10*hy2*x2*y2
                +48*hx12*x2*y2-80*hy14*y2-740*hx2*hy12*y2
                -3000*hx4*hy10*y2-5860*hx6*hy8*y2-5860*hx8*hy6*y2
                -3000*hx10*hy4*y2-740*hx12*hy2*y2-80*hx14*y2
                -8*hx4*hy4*x7*y+24*hx2*hy8*x5*y+240*hx4*hy6*x5*y
                +240*hx6*hy4*x5*y+24*hx8*hy2*x5*y-32*hy12*x3*y
                -432*hx2*hy10*x3*y-2280*hx4*hy8*x3*y
                -3760*hx6*hy6*x3*y-2280*hx8*hy4*x3*y
                -432*hx10*hy2*x3*y-32*hx12*x3*y+160*hy14*x*y
                +1480*hx2*hy12*x*y+6000*hx4*hy10*x*y+11720*hx6*hy8*x*y
                +11720*hx8*hy6*x*y+6000*hx10*hy4*x*y+1480*hx12*hy2*x*y
                +160*hx14*x*y+hx4*hy4*x8-4*hx2*hy8*x6-40*hx4*hy6*x6
                -40*hx6*hy4*x6-4*hx8*hy2*x6+8*hy12*x4
                +108*hx2*hy10*x4+570*hx4*hy8*x4+940*hx6*hy6*x4
                +570*hx8*hy4*x4+108*hx10*hy2*x4+8*hx12*x4
                -80*hy14*x2-740*hx2*hy12*x2-3000*hx4*hy10*x2
                -5860*hx6*hy8*x2-5860*hx8*hy6*x2-3000*hx10*hy4*x2
                -740*hx12*hy2*x2-80*hx14*x2+120*hy16+1020*hx2*hy14
                +3825*hx4*hy12+8040*hx6*hy10+10230*hx8*hy8
                +8040*hx10*hy6+3825*hx12*hy4+1020*hx14*hy2+120*hx16)*
         exp(-0.5*(x-y)*(x-y)/hxy2)*ONE_OVER_SQRT_TWO_PI/(sqrt(hxy2)*64*hxy8));
}

double np_aconvol_gauss8(const double x, const double y,const double hx,const double hy){
  const double hx2 = hx*hx;
  const double hx4 = hx2*hx2;
  const double hx6 = hx4*hx2;
  const double hx8 = hx6*hx2;
  const double hx10 = hx8*hx2;
  const double hx12 = hx10*hx2;
  const double hx14 = hx12*hx2;
  const double hx16 = hx14*hx2;
  const double hx18 = hx16*hx2;
  const double hx20 = hx18*hx2;
  const double hx22 = hx20*hx2;
  const double hx24 = hx22*hx2;

  const double hy2 = hy*hy;
  const double hy4 = hy2*hy2;
  const double hy6 = hy4*hy2;
  const double hy8 = hy6*hy2;
  const double hy10 = hy8*hy2;
  const double hy12 = hy10*hy2;
  const double hy14 = hy12*hy2;
  const double hy16 = hy14*hy2;
  const double hy18 = hy16*hy2;
  const double hy20 = hy18*hy2;
  const double hy22 = hy20*hy2;
  const double hy24 = hy22*hy2;

  const double x2 = x*x;
  const double x3 = x2*x;
  const double x4 = x3*x;
  const double x5 = x4*x;
  const double x6 = x5*x;
  const double x7 = x6*x;
  const double x8 = x7*x;
  const double x9 = x8*x;
  const double x10 = x9*x;
  const double x11 = x10*x;
  const double x12 = x11*x;

  const double y2 = y*y;
  const double y3 = y2*y;
  const double y4 = y3*y;
  const double y5 = y4*y;
  const double y6 = y5*y;
  const double y7 = y6*y;
  const double y8 = y7*y;
  const double y9 = y8*y;
  const double y10 = y9*y;
  const double y11 = y10*y;
  const double y12 = y11*y;

  const double hxy2 = hx2+hy2;

  const double s1 = ONE_OVER_SQRT_TWO_PI*hx*hy*exp(-0.5*(x-y)*(x-y)/hxy2)/(9*256*sqrt(hxy2)*ipow(hxy2,12));
  
  const double s2 = (hx6*hy6*y12-12*hx6*hy6*x*y11+66*hx6*hy6*x2*y10-6*hx4*hy10*y10
                     -84*hx6*hy8*y10-84*hx8*hy6*y10-6*hx10*hy4*y10
                     -220*hx6*hy6*x3*y9+60*hx4*hy10*x*y9
                     +840*hx6*hy8*x*y9+840*hx8*hy6*x*y9+60*hx10*hy4*x*y9
                     +495*hx6*hy6*x4*y8-270*hx4*hy10*x2*y8
                     -3780*hx6*hy8*x2*y8-3780*hx8*hy6*x2*y8
                     -270*hx10*hy4*x2*y8+24*hx2*hy14*y8+402*hx4*hy12*y8
                     +2877*hx6*hy10*y8+4998*hx8*hy8*y8+2877*hx10*hy6*y8
                     +402*hx12*hy4*y8+24*hx14*hy2*y8-792*hx6*hy6*x5*y7
                     +720*hx4*hy10*x3*y7+10080*hx6*hy8*x3*y7
                     +10080*hx8*hy6*x3*y7+720*hx10*hy4*x3*y7
                     -192*hx2*hy14*x*y7-3216*hx4*hy12*x*y7
                     -23016*hx6*hy10*x*y7-39984*hx8*hy8*x*y7
                     -23016*hx10*hy6*x*y7-3216*hx12*hy4*x*y7
                     -192*hx14*hy2*x*y7+924*hx6*hy6*x6*y6
                     -1260*hx4*hy10*x4*y6-17640*hx6*hy8*x4*y6
                     -17640*hx8*hy6*x4*y6-1260*hx10*hy4*x4*y6
                     +672*hx2*hy14*x2*y6+11256*hx4*hy12*x2*y6
                     +80556*hx6*hy10*x2*y6+139944*hx8*hy8*x2*y6
                     +80556*hx10*hy6*x2*y6+11256*hx12*hy4*x2*y6
                     +672*hx14*hy2*x2*y6-48*hy18*y6-1104*hx2*hy16*y6
                     -9876*hx4*hy14*y6-49224*hx6*hy12*y6
                     -105588*hx8*hy10*y6-105588*hx10*hy8*y6
                     -49224*hx12*hy6*y6-9876*hx14*hy4*y6
                     -1104*hx16*hy2*y6-48*hx18*y6-792*hx6*hy6*x7*y5
                     +1512*hx4*hy10*x5*y5+21168*hx6*hy8*x5*y5
                     +21168*hx8*hy6*x5*y5+1512*hx10*hy4*x5*y5
                     -1344*hx2*hy14*x3*y5-22512*hx4*hy12*x3*y5
                     -161112*hx6*hy10*x3*y5-279888*hx8*hy8*x3*y5
                     -161112*hx10*hy6*x3*y5-22512*hx12*hy4*x3*y5
                     -1344*hx14*hy2*x3*y5+288*hy18*x*y5
                     +6624*hx2*hy16*x*y5+59256*hx4*hy14*x*y5
                     +295344*hx6*hy12*x*y5+633528*hx8*hy10*x*y5
                     +633528*hx10*hy8*x*y5+295344*hx12*hy6*x*y5
                     +59256*hx14*hy4*x*y5+6624*hx16*hy2*x*y5
                     +288*hx18*x*y5+495*hx6*hy6*x8*y4
                     -1260*hx4*hy10*x6*y4-17640*hx6*hy8*x6*y4
                     -17640*hx8*hy6*x6*y4-1260*hx10*hy4*x6*y4
                     +1680*hx2*hy14*x4*y4+28140*hx4*hy12*x4*y4
                     +201390*hx6*hy10*x4*y4+349860*hx8*hy8*x4*y4
                     +201390*hx10*hy6*x4*y4+28140*hx12*hy4*x4*y4
                     +1680*hx14*hy2*x4*y4-720*hy18*x2*y4
                     -16560*hx2*hy16*x2*y4-148140*hx4*hy14*x2*y4
                     -738360*hx6*hy12*x2*y4-1583820*hx8*hy10*x2*y4
                     -1583820*hx10*hy8*x2*y4-738360*hx12*hy6*x2*y4
                     -148140*hx14*hy4*x2*y4-16560*hx16*hy2*x2*y4);

    const double s3 = (-720*hx18*x2*y4+1008*hy20*y4+15120*hx2*hy18*y4
                       +102060*hx4*hy16*y4+412335*hx6*hy14*y4
                       +947520*hx8*hy12*y4+1246266*hx10*hy10*y4
                       +947520*hx12*hy8*y4+412335*hx14*hy6*y4
                       +102060*hx16*hy4*y4+15120*hx18*hy2*y4+1008*hx20*y4
                       -220*hx6*hy6*x9*y3+720*hx4*hy10*x7*y3
                       +10080*hx6*hy8*x7*y3+10080*hx8*hy6*x7*y3
                       +720*hx10*hy4*x7*y3-1344*hx2*hy14*x5*y3
                       -22512*hx4*hy12*x5*y3-161112*hx6*hy10*x5*y3
                       -279888*hx8*hy8*x5*y3-161112*hx10*hy6*x5*y3
                       -22512*hx12*hy4*x5*y3-1344*hx14*hy2*x5*y3
                       +960*hy18*x3*y3+22080*hx2*hy16*x3*y3
                       +197520*hx4*hy14*x3*y3+984480*hx6*hy12*x3*y3
                       +2111760*hx8*hy10*x3*y3+2111760*hx10*hy8*x3*y3
                       +984480*hx12*hy6*x3*y3+197520*hx14*hy4*x3*y3
                       +22080*hx16*hy2*x3*y3+960*hx18*x3*y3-4032*hy20*x*y3
                       -60480*hx2*hy18*x*y3-408240*hx4*hy16*x*y3
                       -1649340*hx6*hy14*x*y3-3790080*hx8*hy12*x*y3
                       -4985064*hx10*hy10*x*y3-3790080*hx12*hy8*x*y3
                       -1649340*hx14*hy6*x*y3-408240*hx16*hy4*x*y3
                       -60480*hx18*hy2*x*y3-4032*hx20*x*y3
                       +66*hx6*hy6*x10*y2-270*hx4*hy10*x8*y2
                       -3780*hx6*hy8*x8*y2-3780*hx8*hy6*x8*y2
                       -270*hx10*hy4*x8*y2+672*hx2*hy14*x6*y2
                       +11256*hx4*hy12*x6*y2+80556*hx6*hy10*x6*y2
                       +139944*hx8*hy8*x6*y2+80556*hx10*hy6*x6*y2
                       +11256*hx12*hy4*x6*y2+672*hx14*hy2*x6*y2
                       -720*hy18*x4*y2-16560*hx2*hy16*x4*y2
                       -148140*hx4*hy14*x4*y2-738360*hx6*hy12*x4*y2
                       -1583820*hx8*hy10*x4*y2-1583820*hx10*hy8*x4*y2
                       -738360*hx12*hy6*x4*y2-148140*hx14*hy4*x4*y2
                       -16560*hx16*hy2*x4*y2-720*hx18*x4*y2
                       +6048*hy20*x2*y2+90720*hx2*hy18*x2*y2
                       +612360*hx4*hy16*x2*y2+2474010*hx6*hy14*x2*y2
                       +5685120*hx8*hy12*x2*y2+7477596*hx10*hy10*x2*y2
                       +5685120*hx12*hy8*x2*y2+2474010*hx14*hy6*x2*y2
                       +612360*hx16*hy4*x2*y2+90720*hx18*hy2*x2*y2
                       +6048*hx20*x2*y2-5040*hy22*y2-65520*hx2*hy20*y2
                       -391230*hx4*hy18*y2-1420020*hx6*hy16*y2
                       -3311280*hx8*hy14*y2-5038110*hx10*hy12*y2
                       -5038110*hx12*hy10*y2-3311280*hx14*hy8*y2
                       -1420020*hx16*hy6*y2-391230*hx18*hy4*y2
                       -65520*hx20*hy2*y2-5040*hx22*y2-12*hx6*hy6*x11*y
                       +60*hx4*hy10*x9*y+840*hx6*hy8*x9*y+840*hx8*hy6*x9*y
                       +60*hx10*hy4*x9*y-192*hx2*hy14*x7*y);

  const double s4 = (-3216*hx4*hy12*x7*y-23016*hx6*hy10*x7*y
                     -39984*hx8*hy8*x7*y-23016*hx10*hy6*x7*y
                     -3216*hx12*hy4*x7*y-192*hx14*hy2*x7*y+288*hy18*x5*y
                     +6624*hx2*hy16*x5*y+59256*hx4*hy14*x5*y
                     +295344*hx6*hy12*x5*y+633528*hx8*hy10*x5*y
                     +633528*hx10*hy8*x5*y+295344*hx12*hy6*x5*y
                     +59256*hx14*hy4*x5*y+6624*hx16*hy2*x5*y
                     +288*hx18*x5*y-4032*hy20*x3*y-60480*hx2*hy18*x3*y
                     -408240*hx4*hy16*x3*y-1649340*hx6*hy14*x3*y
                     -3790080*hx8*hy12*x3*y-4985064*hx10*hy10*x3*y
                     -3790080*hx12*hy8*x3*y-1649340*hx14*hy6*x3*y
                     -408240*hx16*hy4*x3*y-60480*hx18*hy2*x3*y
                     -4032*hx20*x3*y+10080*hy22*x*y+131040*hx2*hy20*x*y
                     +782460*hx4*hy18*x*y+2840040*hx6*hy16*x*y
                     +6622560*hx8*hy14*x*y+10076220*hx10*hy12*x*y
                     +10076220*hx12*hy10*x*y+6622560*hx14*hy8*x*y
                     +2840040*hx16*hy6*x*y+782460*hx18*hy4*x*y
                     +131040*hx20*hy2*x*y+10080*hx22*x*y+hx6*hy6*x12
                     -6*hx4*hy10*x10-84*hx6*hy8*x10-84*hx8*hy6*x10
                     -6*hx10*hy4*x10+24*hx2*hy14*x8+402*hx4*hy12*x8
                     +2877*hx6*hy10*x8+4998*hx8*hy8*x8+2877*hx10*hy6*x8
                     +402*hx12*hy4*x8+24*hx14*hy2*x8-48*hy18*x6
                     -1104*hx2*hy16*x6-9876*hx4*hy14*x6
                     -49224*hx6*hy12*x6-105588*hx8*hy10*x6
                     -105588*hx10*hy8*x6-49224*hx12*hy6*x6
                     -9876*hx14*hy4*x6-1104*hx16*hy2*x6-48*hx18*x6
                     +1008*hy20*x4+15120*hx2*hy18*x4+102060*hx4*hy16*x4
                     +412335*hx6*hy14*x4+947520*hx8*hy12*x4
                     +1246266*hx10*hy10*x4+947520*hx12*hy8*x4
                     +412335*hx14*hy6*x4+102060*hx16*hy4*x4
                     +15120*hx18*hy2*x4+1008*hx20*x4-5040*hy22*x2
                     -65520*hx2*hy20*x2-391230*hx4*hy18*x2
                     -1420020*hx6*hy16*x2-3311280*hx8*hy14*x2
                     -5038110*hx10*hy12*x2-5038110*hx12*hy10*x2
                     -3311280*hx14*hy8*x2-1420020*hx16*hy6*x2
                     -391230*hx18*hy4*x2-65520*hx20*hy2*x2-5040*hx22*x2
                     +5040*hy24+63000*hx2*hy22+362250*hx4*hy20
                     +1267875*hx6*hy18+2983050*hx8*hy16+4923765*hx10*hy14
                     +5808600*hx12*hy12+4923765*hx14*hy10+2983050*hx16*hy8
                     +1267875*hx18*hy6+362250*hx20*hy4+63000*hx22*hy2
                     +5040*hx24);
  return(s1*(s2+s3+s4));
}

double np_aconvol_epan2_total(const double x, const double y,const double hx,const double hy){
  const double a = 3.0*sqrt(5.0);
  const double hl = MAX(hx,hy);
  const double hs = MIN(hx,hy);
  return((-a*y*y + 2.0*a*x*y - a*x*x + 5.0*a*hl*hl - a*hs*hs)*hs/(100.0*hl*hl));
}

double np_aconvol_epan2_indefinite(const double l, const double x, const double y,const double hx,const double hy){
  const double hxs = hx*hx;
  const double hys = hy*hy;
  const double xs = x*x;
  const double ys = y*y;
  const double a = 3.0/(20000.0*hxs*hys);
    
  return(a*l*(((30.0*xs - 150.0*hxs)*ys + hys*(-150.0*xs + 750.0*hxs)) +
              l*(((150.0*hxs - 30.0*xs)*y + (150.0*hys - 30.0*ys)*x) +
                 l*(10.0*(ys + 4.0*x*y + xs - 5.0*hys - 5.0*hxs) +
                    l*((y + x)*(-15.0) + 6.0*l)))));

}

double np_aconvol_epan2(const double x, const double y,const double hx,const double hy){
  const double a = sqrt(5.0);
  const double dxy = fabs(x-y);

  if(dxy >= a*(hx+hy)){
    return 0;
  } else if(dxy > a*fabs(hx-hy)){
    return (np_aconvol_epan2_indefinite(MIN(x+a*hx,y+a*hy),x,y,hx,hy) - 
            np_aconvol_epan2_indefinite(MAX(x-a*hx,y-a*hy),x,y,hx,hy));
  } else {
    return (np_aconvol_epan2_total(x,y,hx,hy));
  }
}

double np_aconvol_epan4_total(const double x, const double y,const double hx,const double hy){
  const double hl = MAX(hx,hy);
  const double hs = MIN(hx,hy);

  const double x2 = x*x;
  const double x3 = x2*x;
  const double x4 = x2*x2;

  const double y2 = y*y;
  const double y3 = y2*y;
  const double y4 = y2*y2;

  const double hl2 = hl*hl;
  const double hl4 = hl2*hl2;

  const double hs2 = hs*hs;
  const double hs4 = hs2*hs2;

  return(hs*(21*y4-84*x*y3+126*x2*y2-150*hl2*y2-84*x3*y+300*hl2*x*y+21*x4
             -150*hl2*x2+225*hl4-25*hs4)/(32*5*sqrt(5)*hl4));
}

double np_aconvol_epan4_indefinite(const double l,const double x, const double y,const double hx,const double hy){
  const double l2 = l*l;
  const double l3 = l2*l;
  const double l4 = l3*l;
  const double l5 = l4*l;
  const double l6 = l5*l;
  const double l7 = l6*l;
  const double l8 = l7*l;

  const double x2 = x*x;
  const double x3 = x2*x;
  const double x4 = x3*x;

  const double y2 = y*y;
  const double y3 = y2*y;
  const double y4 = y3*y;

  const double hx2 = hx*hx;
  const double hx4 = hx2*hx2;

  const double hy2 = hy*hy;
  const double hy4 = hy2*hy2;

  return(l*(4410*x4*y4-8820*l*x3*y4+8820*l2*x2*y4-31500*hx2*x2*y4
            -4410*l3*x*y4+31500*hx2*l*x*y4+882*l4*y4
            -10500*hx2*l2*y4+47250*hx4*y4-8820*l*x4*y3
            +23520*l2*x3*y3-26460*l3*x2*y3+63000*hx2*l*x2*y3
            +14112*l4*x*y3-84000*hx2*l2*x*y3-2940*l5*y3
            +31500*hx2*l3*y3-94500*hx4*l*y3+8820*l2*x4*y2
            -31500*hy2*x4*y2-26460*l3*x3*y2+63000*hy2*l*x3*y2
            +31752*l4*x2*y2-63000*hy2*l2*x2*y2
            -63000*hx2*l2*x2*y2+225000*hx2*hy2*x2*y2
            -17640*l5*x*y2+31500*hy2*l3*x*y2+94500*hx2*l3*x*y2
            -225000*hx2*hy2*l*x*y2+3780*l6*y2-6300*hy2*l4*y2
            -37800*hx2*l4*y2+75000*hx2*hy2*l2*y2+94500*hx4*l2*y2
            -337500*hx4*hy2*y2-4410*l3*x4*y+31500*hy2*l*x4*y
            +14112*l4*x3*y-84000*hy2*l2*x3*y-17640*l5*x2*y
            +94500*hy2*l3*x2*y+31500*hx2*l3*x2*y
            -225000*hx2*hy2*l*x2*y+10080*l6*x*y-50400*hy2*l4*x*y
            -50400*hx2*l4*x*y+300000*hx2*hy2*l2*x*y-2205*l7*y
            +10500*hy2*l5*y+21000*hx2*l5*y-112500*hx2*hy2*l3*y
            -47250*hx4*l3*y+337500*hx4*hy2*l*y+882*l4*x4
            -10500*hy2*l2*x4+47250*hy4*x4-2940*l5*x3
            +31500*hy2*l3*x3-94500*hy4*l*x3+3780*l6*x2
            -37800*hy2*l4*x2-6300*hx2*l4*x2+94500*hy4*l2*x2
            +75000*hx2*hy2*l2*x2-337500*hx2*hy4*x2-2205*l7*x
            +21000*hy2*l5*x+10500*hx2*l5*x-47250*hy4*l3*x
            -112500*hx2*hy2*l3*x+337500*hx2*hy4*l*x+490*l8
            -4500*hy2*l6-4500*hx2*l6+9450*hy4*l4+45000*hx2*hy2*l4
            +9450*hx4*l4-112500*hx2*hy4*l2-112500*hx4*hy2*l2
            +506250*hx4*hy4)
         /(1280000*hx4*hy4));
}

double np_aconvol_epan4(const double x, const double y,const double hx,const double hy){
  const double a = sqrt(5.0);
  const double dxy = fabs(x-y);

  if(dxy >= a*(hx+hy)){
    return 0;
  } else if(dxy > a*fabs(hx-hy)){
    return (np_aconvol_epan4_indefinite(MIN(x+a*hx,y+a*hy),x,y,hx,hy) - 
            np_aconvol_epan4_indefinite(MAX(x-a*hx,y-a*hy),x,y,hx,hy));
  } else {
    return (np_aconvol_epan4_total(x,y,hx,hy));
  }
}

double np_aconvol_epan6_total(const double x, const double y,const double hx,const double hy){
  const double hl = MAX(hx,hy);
  const double hs = MIN(hx,hy);

  const double x2 = x*x;
  const double x3 = x2*x;
  const double x4 = x2*x2;
  const double x5 = x3*x2;
  const double x6 = x3*x3;

  const double y2 = y*y;
  const double y3 = y2*y;
  const double y4 = y2*y2;
  const double y5 = y3*y2;
  const double y6 = y3*y3;

  const double hl2 = hl*hl;
  const double hl4 = hl2*hl2;
  const double hl6 = hl4*hl2;

  const double hs2 = hs*hs;
  const double hs4 = hs2*hs2;
  const double hs6 = hs4*hs2;

  return(-21*hs
         *(429*y6-2574*x*y5+6435*x2*y4-4095*hl2*y4-8580*x3*y3
           +16380*hl2*x*y3+6435*x4*y2-24570*hl2*x2*y2+11375*hl4*y2
           -2574*x5*y+16380*hl2*x3*y-22750*hl4*x*y+429*x6-4095*hl2*x4
           +11375*hl4*x2-8125*hl6+625*hs6)
         /(3328*sqrt(5)*25*hl6));

}

double np_aconvol_epan6_indefinite(const double u, const double x, const double y,const double hx,const double hy){

  const double x2 = x*x;
  const double x3 = x2*x;
  const double x4 = x2*x2;
  const double x5 = x3*x2;
  const double x6 = x3*x3;

  const double y2 = y*y;
  const double y3 = y2*y;
  const double y4 = y2*y2;
  const double y5 = y3*y2;
  const double y6 = y3*y3;

  const double hx2 = hx*hx;
  const double hx4 = hx2*hx2;
  const double hx6 = hx4*hx2;

  const double hy2 = hy*hy;
  const double hy4 = hy2*hy2;
  const double hy6 = hy4*hy2;

  const double u2 = u*u;
  const double u3 = u2*u;
  const double u4 = u3*u;
  const double u5 = u4*u;
  const double u6 = u5*u;
  const double u7 = u6*u;
  const double u8 = u7*u;
  const double u9 = u8*u;
  const double u10 = u9*u;
  const double u11 = u10*u;
  const double u12 = u11*u;

  return(21*u
         *(1189188*x6*y6-3567564*u*x5*y6+5945940*u2*x4*y6
           -11351340*hx2*x4*y6-5945940*u3*x3*y6
           +22702680*hx2*u*x3*y6+3567564*u4*x2*y6
           -22702680*hx2*u2*x2*y6+31531500*hx4*x2*y6
           -1189188*u5*x*y6+11351340*hx2*u3*x*y6
           -31531500*hx4*u*x*y6+169884*u6*y6-2270268*hx2*u4*y6
           +10510500*hx4*u2*y6-22522500*hx6*y6-3567564*u*x6*y5
           +14270256*u2*x5*y5-26756730*u3*x4*y5
           +34054020*hx2*u*x4*y5+28540512*u4*x3*y5
           -90810720*hx2*u2*x3*y5-17837820*u5*x2*y5
           +102162060*hx2*u3*x2*y5-94594500*hx4*u*x2*y5
           +6115824*u6*x*y5-54486432*hx2*u4*x*y5
           +126126000*hx4*u2*x*y5-891891*u7*y5
           +11351340*hx2*u5*y5-47297250*hx4*u3*y5
           +67567500*hx6*u*y5+5945940*u2*x6*y4
           -11351340*hy2*x6*y4-26756730*u3*x5*y4
           +34054020*hy2*u*x5*y4+53513460*u4*x4*y4
           -56756700*hy2*u2*x4*y4-56756700*hx2*u2*x4*y4
           +108353700*hx2*hy2*x4*y4-59459400*u5*x3*y4
           +56756700*hy2*u3*x3*y4+170270100*hx2*u3*x3*y4
           -216707400*hx2*hy2*u*x3*y4+38223900*u6*x2*y4
           -34054020*hy2*u4*x2*y4-204324120*hx2*u4*x2*y4
           +216707400*hx2*hy2*u2*x2*y4+157657500*hx4*u2*x2*y4
           -300982500*hx4*hy2*x2*y4-13378365*u7*x*y4
           +11351340*hy2*u5*x*y4+113513400*hx2*u5*x*y4
           -108353700*hx2*hy2*u3*x*y4-236486250*hx4*u3*x*y4
           +300982500*hx4*hy2*u*x*y4+1981980*u8*y4
           -1621620*hy2*u6*y4-24324300*hx2*u6*y4
           +21670740*hx2*hy2*u4*y4+94594500*hx4*u4*y4
           -100327500*hx4*hy2*u2*y4-112612500*hx6*u2*y4
           +214987500*hx6*hy2*y4-5945940*u3*x6*y3
           +22702680*hy2*u*x6*y3+28540512*u4*x5*y3
           -90810720*hy2*u2*x5*y3-59459400*u5*x4*y3
           +170270100*hy2*u3*x4*y3+56756700*hx2*u3*x4*y3
           -216707400*hx2*hy2*u*x4*y3+67953600*u6*x3*y3
           -181621440*hy2*u4*x3*y3-181621440*hx2*u4*x3*y3
           +577886400*hx2*hy2*u2*x3*y3-44594550*u7*x2*y3
           +113513400*hy2*u5*x2*y3+227026800*hx2*u5*x2*y3
           -650122200*hx2*hy2*u3*x2*y3-157657500*hx4*u3*x2*y3
           +601965000*hx4*hy2*u*x2*y3+15855840*u8*x*y3
           -38918880*hy2*u6*x*y3-129729600*hx2*u6*x*y3
           +346731840*hx2*hy2*u4*x*y3+252252000*hx4*u4*x*y3
           -802620000*hx4*hy2*u2*x*y3-2378376*u9*y3
           +5675670*hy2*u7*y3+28378350*hx2*u7*y3
           -72235800*hx2*hy2*u5*y3-105105000*hx4*u5*y3
           +300982500*hx4*hy2*u3*y3+112612500*hx6*u3*y3
           -429975000*hx6*hy2*u*y3+3567564*u4*x6*y2
           -22702680*hy2*u2*x6*y2+31531500*hy4*x6*y2
           -17837820*u5*x5*y2+102162060*hy2*u3*x5*y2
           -94594500*hy4*u*x5*y2+38223900*u6*x4*y2
           -204324120*hy2*u4*x4*y2-34054020*hx2*u4*x4*y2
           +157657500*hy4*u2*x4*y2+216707400*hx2*hy2*u2*x4*y2
           -300982500*hx2*hy4*x4*y2-44594550*u7*x3*y2
           +227026800*hy2*u5*x3*y2+113513400*hx2*u5*x3*y2
           -157657500*hy4*u3*x3*y2-650122200*hx2*hy2*u3*x3*y2
           +601965000*hx2*hy4*u*x3*y2+29729700*u8*x2*y2
           -145945800*hy2*u6*x2*y2-145945800*hx2*u6*x2*y2
           +94594500*hy4*u4*x2*y2+780146640*hx2*hy2*u4*x2*y2
           +94594500*hx4*u4*x2*y2-601965000*hx2*hy4*u2*x2*y2
           -601965000*hx4*hy2*u2*x2*y2
           +836062500*hx4*hy4*x2*y2-10702692*u9*x*y2
           +51081030*hy2*u7*x*y2+85135050*hx2*u7*x*y2
           -31531500*hy4*u5*x*y2-433414800*hx2*hy2*u5*x*y2
           -157657500*hx4*u5*x*y2+300982500*hx2*hy4*u3*x*y2
           +902947500*hx4*hy2*u3*x*y2-836062500*hx4*hy4*u*x*y2
           +1621620*u10*y2-7567560*hy2*u8*y2
           -18918900*hx2*u8*y2+4504500*hy4*u6*y2
           +92874600*hx2*hy2*u6*y2+67567500*hx4*u6*y2
           -60196500*hx2*hy4*u4*y2-361179000*hx4*hy2*u4*y2
           -67567500*hx6*u4*y2+278687500*hx4*hy4*u2*y2
           +429975000*hx6*hy2*u2*y2-597187500*hx6*hy4*y2
           -1189188*u5*x6*y+11351340*hy2*u3*x6*y
           -31531500*hy4*u*x6*y+6115824*u6*x5*y
           -54486432*hy2*u4*x5*y+126126000*hy4*u2*x5*y
           -13378365*u7*x4*y+113513400*hy2*u5*x4*y
           +11351340*hx2*u5*x4*y-236486250*hy4*u3*x4*y
           -108353700*hx2*hy2*u3*x4*y+300982500*hx2*hy4*u*x4*y
           +15855840*u8*x3*y-129729600*hy2*u6*x3*y
           -38918880*hx2*u6*x3*y+252252000*hy4*u4*x3*y
           +346731840*hx2*hy2*u4*x3*y
           -802620000*hx2*hy4*u2*x3*y-10702692*u9*x2*y
           +85135050*hy2*u7*x2*y+51081030*hx2*u7*x2*y
           -157657500*hy4*u5*x2*y-433414800*hx2*hy2*u5*x2*y
           -31531500*hx4*u5*x2*y+902947500*hx2*hy4*u3*x2*y
           +300982500*hx4*hy2*u3*x2*y-836062500*hx4*hy4*u*x2*y
           +3891888*u10*x*y-30270240*hy2*u8*x*y
           -30270240*hx2*u8*x*y+54054000*hy4*u6*x*y
           +247665600*hx2*hy2*u6*x*y+54054000*hx4*u6*x*y
           -481572000*hx2*hy4*u4*x*y-481572000*hx4*hy2*u4*x*y
           +1114750000*hx4*hy4*u2*x*y-594594*u11*y
           +4540536*hy2*u9*y+6810804*hx2*u9*y-7882875*hy4*u7*y
           -54176850*hx2*hy2*u7*y-23648625*hx4*u7*y
           +100327500*hx2*hy4*u5*y+200655000*hx4*hy2*u5*y
           +22522500*hx6*u5*y-418031250*hx4*hy4*u3*y
           -214987500*hx6*hy2*u3*y+597187500*hx6*hy4*u*y
           +169884*u6*x6-2270268*hy2*u4*x6+10510500*hy4*u2*x6
           -22522500*hy6*x6-891891*u7*x5+11351340*hy2*u5*x5
           -47297250*hy4*u3*x5+67567500*hy6*u*x5+1981980*u8*x4
           -24324300*hy2*u6*x4-1621620*hx2*u6*x4
           +94594500*hy4*u4*x4+21670740*hx2*hy2*u4*x4
           -112612500*hy6*u2*x4-100327500*hx2*hy4*u2*x4
           +214987500*hx2*hy6*x4-2378376*u9*x3
           +28378350*hy2*u7*x3+5675670*hx2*u7*x3
           -105105000*hy4*u5*x3-72235800*hx2*hy2*u5*x3
           +112612500*hy6*u3*x3+300982500*hx2*hy4*u3*x3
           -429975000*hx2*hy6*u*x3+1621620*u10*x2
           -18918900*hy2*u8*x2-7567560*hx2*u8*x2
           +67567500*hy4*u6*x2+92874600*hx2*hy2*u6*x2
           +4504500*hx4*u6*x2-67567500*hy6*u4*x2
           -361179000*hx2*hy4*u4*x2-60196500*hx4*hy2*u4*x2
           +429975000*hx2*hy6*u2*x2+278687500*hx4*hy4*u2*x2
           -597187500*hx4*hy6*x2-594594*u11*x+6810804*hy2*u9*x
           +4540536*hx2*u9*x-23648625*hy4*u7*x
           -54176850*hx2*hy2*u7*x-7882875*hx4*u7*x
           +22522500*hy6*u5*x+200655000*hx2*hy4*u5*x
           +100327500*hx4*hy2*u5*x-214987500*hx2*hy6*u3*x
           -418031250*hx4*hy4*u3*x+597187500*hx4*hy6*u*x
           +91476*u12-1031940*hy2*u10-1031940*hx2*u10
           +3503500*hy4*u8+12039300*hx2*hy2*u8+3503500*hx4*u8
           -3217500*hy6*u6-42997500*hx2*hy4*u6
           -42997500*hx4*hy2*u6-3217500*hx6*u6
           +42997500*hx2*hy6*u4+167212500*hx4*hy4*u4
           +42997500*hx6*hy2*u4-199062500*hx4*hy6*u2
           -199062500*hx6*hy4*u2+426562500*hx6*hy6)
         /(10649600000*hx6*hy6));
}

double np_aconvol_epan6(const double x, const double y,const double hx,const double hy){
  const double a = sqrt(5.0);
  const double dxy = fabs(x-y);

  if(dxy >= a*(hx+hy)){
    return 0;
  } else if(dxy > a*fabs(hx-hy)){
    return (np_aconvol_epan6_indefinite(MIN(x+a*hx,y+a*hy),x,y,hx,hy) - 
            np_aconvol_epan6_indefinite(MAX(x-a*hx,y-a*hy),x,y,hx,hy));
  } else {
    return (np_aconvol_epan6_total(x,y,hx,hy));
  }
}

double np_aconvol_epan8_total(const double x, const double y,const double hx,const double hy){
  const double hl = MAX(hx,hy);
  const double hs = MIN(hx,hy);

  const double x2 = x*x;
  const double x3 = x2*x;
  const double x4 = x2*x2;
  const double x5 = x3*x2;
  const double x6 = x3*x3;
  const double x7 = x4*x3;
  const double x8 = x4*x4;

  const double y2 = y*y;
  const double y3 = y2*y;
  const double y4 = y2*y2;
  const double y5 = y3*y2;
  const double y6 = y3*y3;
  const double y7 = y4*y3;
  const double y8 = y4*y4;

  const double hl2 = hl*hl;
  const double hl4 = hl2*hl2;
  const double hl6 = hl4*hl2;
  const double hl8 = hl4*hl4;

  const double hs2 = hs*hs;
  const double hs4 = hs2*hs2;
  const double hs8 = hs4*hs4;

  return(63*hs
         *(2431*y8-19448*x*y7+68068*x2*y6-29172*hl2*y6-136136*x3*y5
           +175032*hl2*x*y5+170170*x4*y4-437580*hl2*x2*y4
           +117810*hl4*y4-136136*x5*y3+583440*hl2*x3*y3
           -471240*hl4*x*y3+68068*x6*y2-437580*hl2*x4*y2
           +706860*hl4*x2*y2-178500*hl6*y2-19448*x7*y+175032*hl2*x5*y
           -471240*hl4*x3*y+357000*hl6*x*y+2431*x8-29172*hl2*x6
           +117810*hl4*x4-178500*hl6*x2+74375*hl8-4375*hs8)
         /(69632*sqrt(5)*25*hl8));
}

double np_aconvol_epan8_xlessy(const double x, const double y,const double hx,const double hy){
  const double a = sqrt(5);
  const double y2 = y*y;
  const double y3 = y2*y;
  const double y4 = y2*y2;
  const double y5 = y3*y2;
  const double y6 = y3*y3;
  const double y7 = y4*y3;
  const double y8 = y4*y4;
  const double y9 = y5*y4;
  const double y10 = y6*y4;
  const double y11 = y7*y4;
  const double y12 = y8*y4;
  const double y13 = y9*y4;
  const double y14 = y10*y4;
  const double y15 = y11*y4;
  const double y16 = y12*y4;
  const double y17 = y13*y4;
 
  const double x2 = x*x;
  const double x3 = x2*x;
  const double x4 = x2*x2;
  const double x5 = x3*x2;
  const double x6 = x3*x3;
  const double x7 = x4*x3;
  const double x8 = x4*x4;
  const double x9 = x5*x4;
  const double x10 = x6*x4;
  const double x11 = x7*x4;
  const double x12 = x8*x4;
  const double x13 = x9*x4;
  const double x14 = x10*x4;
  const double x15 = x11*x4;
  const double x16 = x12*x4;
  const double x17 = x13*x4;

  const double hx2 = hx*hx;
  const double hx3 = hx2*hx;
  const double hx4 = hx2*hx2;
  const double hx5 = hx3*hx2;
  const double hx6 = hx3*hx3;
  const double hx8 = hx4*hx4;
  const double hx9 = hx5*hx4;
  const double hx10 = hx6*hx4;
  const double hx12 = hx8*hx4;
  const double hx13 = hx9*hx4;
  const double hx14 = hx10*hx4;
  const double hx16 = hx12*hx4;
  const double hx17 = hx13*hx4;

  const double hy2 = hy*hy;
  const double hy3 = hy2*hy;
  const double hy4 = hy2*hy2;
  const double hy5 = hy3*hy2;
  const double hy6 = hy3*hy3;
  const double hy8 = hy4*hy4;
  const double hy9 = hy5*hy4;
  const double hy10 = hy6*hy4;
  const double hy12 = hy8*hy4;
  const double hy13 = hy9*hy4;
  const double hy14 = hy10*hy4;
  const double hy16 = hy12*hy4;
  const double hy17 = hy13*hy4;


  return(-63*(1001*y17-17017*x*y16+136136*x2*y15-58344*hy2*y15-58344*hx2*y15
              -680680*x3*y14+875160*hy2*x*y14+875160*hx2*x*y14
              +2382380*x4*y13-6126120*hy2*x2*y13-6126120*hx2*x2*y13
              +1649340*hy4*y13+2625480*hx2*hy2*y13+1649340*hx4*y13
              -6194188*x5*y12+26546520*hy2*x3*y12+26546520*hx2*x3*y12
              -21441420*hy4*x*y12-34131240*hx2*hy2*x*y12
              -21441420*hx4*x*y12+12388376*x6*y11-79639560*hy2*x4*y11
              -79639560*hx2*x4*y11+128648520*hy4*x2*y11
              +204787440*hx2*hy2*x2*y11+128648520*hx4*x2*y11
              -32487000*hy6*y11-55135080*hx2*hy4*y11
              -55135080*hx4*hy2*y11-32487000*hx6*y11-19467448*x7*y10
              +175207032*hy2*x5*y10+175207032*hx2*x5*y10
              -471711240*hy4*x3*y10-750887280*hx2*hy2*x3*y10
              -471711240*hx4*x3*y10+357357000*hy6*x*y10
              +606485880*hx2*hy4*x*y10+606485880*hx4*hy2*x*y10
              +357357000*hx6*x*y10+24334310*x8*y9-292011720*hy2*x6*y9
              -292011720*hx2*x6*y9+1179278100*hy4*x4*y9
              +1877218200*hx2*hy2*x4*y9+1179278100*hx4*x4*y9
              -1786785000*hy6*x2*y9-3032429400*hx2*hy4*x2*y9
              -3032429400*hx4*hy2*x2*y9-1786785000*hx6*x2*y9
              +744493750*hy8*y9+765765000*hx2*hy6*y9
              +816423300*hx4*hy4*y9+765765000*hx6*hy2*y9
              +744493750*hx8*y9-24334310*x9*y8+375443640*hy2*x7*y8
              +375443640*hx2*x7*y8-2122700580*hy4*x5*y8
              -3378992760*hx2*hy2*x5*y8-2122700580*hx4*x5*y8
              +5360355000*hy6*x3*y8+9097288200*hx2*hy4*x3*y8
              +9097288200*hx4*hy2*x3*y8+5360355000*hx6*x3*y8
              -6700443750*hy8*x*y8-6891885000*hx2*hy6*x*y8
              -7347809700*hx4*hy4*x*y8-6891885000*hx6*hy2*x*y8
              -6700443750*hx8*x*y8-9957376*a*125*hy9*y8
              -9957376*a*125*hx9*y8+19467448*x10*y7
              -375443640*hy2*x8*y7-375443640*hx2*x8*y7
              +2830267440*hy4*x6*y7+4505323680*hx2*hy2*x6*y7
              +2830267440*hx4*x6*y7-10720710000*hy6*x4*y7
              -18194576400*hx2*hy4*x4*y7-18194576400*hx4*hy2*x4*y7
              -10720710000*hx6*x4*y7+26801775000*hy8*x2*y7
              +27567540000*hx2*hy6*x2*y7+29391238800*hx4*hy4*x2*y7
              +27567540000*hx6*hy2*x2*y7+26801775000*hx8*x2*y7
              +79659008*a*125*hy9*x*y7+79659008*a*125*hx9*x*y7
              +3828825000*hy10*y7-11486475000*hx2*hy8*y7
              -7422030000*hx4*hy6*y7-7422030000*hx6*hy4*y7
              -11486475000*hx8*hy2*y7+3828825000*hx10*y7
              -12388376*x11*y6+292011720*hy2*x9*y6+292011720*hx2*x9*y6
              -2830267440*hy4*x7*y6-4505323680*hx2*hy2*x7*y6
              -2830267440*hx4*x7*y6+15008994000*hy6*x5*y6
              +25472406960*hx2*hy4*x5*y6+25472406960*hx4*hy2*x5*y6
              +15008994000*hx6*x5*y6-62537475000*hy8*x3*y6
              -64324260000*hx2*hy6*x3*y6-68579557200*hx4*hy4*x3*y6
              -64324260000*hx6*hy2*x3*y6-62537475000*hx8*x3*y6
              -278806528*a*125*hy9*x2*y6-278806528*a*125*hx9*x2*y6
              -26801775000*hy10*x*y6+80405325000*hx2*hy8*x*y6
              +51954210000*hx4*hy6*x*y6+51954210000*hx6*hy4*x*y6
              +80405325000*hx8*hy2*x*y6-26801775000*hx10*x*y6
              +119488512*a*125*hx2*hy9*y6+119488512*a*125*hx9*hy2*y6
              +6194188*x12*y5-175207032*hy2*x10*y5
              -175207032*hx2*x10*y5+2122700580*hy4*x8*y5
              +3378992760*hx2*hy2*x8*y5+2122700580*hx4*x8*y5
              -15008994000*hy6*x6*y5-25472406960*hx2*hy4*x6*y5
              -25472406960*hx4*hy2*x6*y5-15008994000*hx6*x6*y5
              +93806212500*hy8*x4*y5+96486390000*hx2*hy6*x4*y5
              +102869335800*hx4*hy4*x4*y5+96486390000*hx6*hy2*x4*y5
              +93806212500*hx8*x4*y5+557613056*a*125*hy9*x3*y5
              +557613056*a*125*hx9*x3*y5+80405325000*hy10*x2*y5
              -241215975000*hx2*hy8*x2*y5-155862630000*hx4*hy6*x2*y5
              -155862630000*hx6*hy4*x2*y5-241215975000*hx8*hy2*x2*y5
              +80405325000*hx10*x2*y5-716931072*a*125*hx2*hy9*x*y5
              -716931072*a*125*hx9*hy2*x*y5-4466962500*hy12*y5
              -34459425000*hx2*hy10*y5+64942762500*hx4*hy8*y5
              +39359250000*hx6*hy6*y5+64942762500*hx8*hy4*y5
              -34459425000*hx10*hy2*y5-4466962500*hx12*y5
              -2382380*x13*y4+79639560*hy2*x11*y4+79639560*hx2*x11*y4
              -1179278100*hy4*x9*y4-1877218200*hx2*hy2*x9*y4
              -1179278100*hx4*x9*y4+10720710000*hy6*x7*y4
              +18194576400*hx2*hy4*x7*y4+18194576400*hx4*hy2*x7*y4
              +10720710000*hx6*x7*y4-93806212500*hy8*x5*y4
              -96486390000*hx2*hy6*x5*y4-102869335800*hx4*hy4*x5*y4
              -96486390000*hx6*hy2*x5*y4-93806212500*hx8*x5*y4
              -139403264*a*625*hy9*x4*y4-139403264*a*625*hx9*x4*y4
              -134008875000*hy10*x3*y4+402026625000*hx2*hy8*x3*y4
              +259771050000*hx4*hy6*x3*y4+259771050000*hx6*hy4*x3*y4
              +402026625000*hx8*hy2*x3*y4-134008875000*hx10*x3*y4
              +358465536*a*625*hx2*hy9*x2*y4
              +358465536*a*625*hx9*hy2*x2*y4+22334812500*hy12*x*y4
              +172297125000*hx2*hy10*x*y4-324713812500*hx4*hy8*x*y4
              -196796250000*hx6*hy6*x*y4-324713812500*hx8*hy4*x*y4
              +172297125000*hx10*hy2*x*y4+22334812500*hx12*x*y4
              -96509952*a*625*hx4*hy9*y4-96509952*a*625*hx9*hy4*y4
              +680680*x14*y3-26546520*hy2*x12*y3-26546520*hx2*x12*y3
              +471711240*hy4*x10*y3+750887280*hx2*hy2*x10*y3
              +471711240*hx4*x10*y3-5360355000*hy6*x8*y3
              -9097288200*hx2*hy4*x8*y3-9097288200*hx4*hy2*x8*y3
              -5360355000*hx6*x8*y3+62537475000*hy8*x6*y3
              +64324260000*hx2*hy6*x6*y3+68579557200*hx4*hy4*x6*y3
              +64324260000*hx6*hy2*x6*y3+62537475000*hx8*x6*y3
              +557613056*a*125*hy9*x5*y3+557613056*a*125*hx9*x5*y3
              +134008875000*hy10*x4*y3-402026625000*hx2*hy8*x4*y3
              -259771050000*hx4*hy6*x4*y3-259771050000*hx6*hy4*x4*y3
              -402026625000*hx8*hy2*x4*y3+134008875000*hx10*x4*y3
              -477954048*a*625*hx2*hy9*x3*y3
              -477954048*a*625*hx9*hy2*x3*y3-44669625000*hy12*x2*y3
              -344594250000*hx2*hy10*x2*y3+649427625000*hx4*hy8*x2*y3
              +393592500000*hx6*hy6*x2*y3+649427625000*hx8*hy4*x2*y3
              -344594250000*hx10*hy2*x2*y3-44669625000*hx12*x2*y3
              +386039808*a*625*hx4*hy9*x*y3
              +386039808*a*625*hx9*hy4*x*y3+6381375000*hy14*y3
              +19144125000*hx2*hy12*y3+92775375000*hx4*hy10*y3
              -163996875000*hx6*hy8*y3-163996875000*hx8*hy6*y3
              +92775375000*hx10*hy4*y3+19144125000*hx12*hy2*y3
              +6381375000*hx14*y3-136136*x15*y2+6126120*hy2*x13*y2
              +6126120*hx2*x13*y2-128648520*hy4*x11*y2
              -204787440*hx2*hy2*x11*y2-128648520*hx4*x11*y2
              +1786785000*hy6*x9*y2+3032429400*hx2*hy4*x9*y2
              +3032429400*hx4*hy2*x9*y2+1786785000*hx6*x9*y2
              -26801775000*hy8*x7*y2-27567540000*hx2*hy6*x7*y2
              -29391238800*hx4*hy4*x7*y2-27567540000*hx6*hy2*x7*y2
              -26801775000*hx8*x7*y2-278806528*a*125*hy9*x6*y2
              -278806528*a*125*hx9*x6*y2-80405325000*hy10*x5*y2
              +241215975000*hx2*hy8*x5*y2+155862630000*hx4*hy6*x5*y2
              +155862630000*hx6*hy4*x5*y2+241215975000*hx8*hy2*x5*y2
              -80405325000*hx10*x5*y2+358465536*a*625*hx2*hy9*x4*y2
              +358465536*a*625*hx9*hy2*x4*y2+44669625000*hy12*x3*y2
              +344594250000*hx2*hy10*x3*y2-649427625000*hx4*hy8*x3*y2
              -393592500000*hx6*hy6*x3*y2-649427625000*hx8*hy4*x3*y2
              +344594250000*hx10*hy2*x3*y2+44669625000*hx12*x3*y2
              -579059712*a*625*hx4*hy9*x2*y2
              -579059712*a*625*hx9*hy4*x2*y2-19144125000*hy14*x*y2
              -57432375000*hx2*hy12*x*y2-278326125000*hx4*hy10*x*y2
              +491990625000*hx6*hy8*x*y2+491990625000*hx8*hy6*x*y2
              -278326125000*hx10*hy4*x*y2-57432375000*hx12*hy2*x*y2
              -19144125000*hx14*x*y2+5849088*a*125*125*hx6*hy9*y2
              +5849088*a*125*125*hx9*hy6*y2+17017*x16*y-875160*hy2*x14*y
              -875160*hx2*x14*y+21441420*hy4*x12*y
              +34131240*hx2*hy2*x12*y+21441420*hx4*x12*y
              -357357000*hy6*x10*y-606485880*hx2*hy4*x10*y
              -606485880*hx4*hy2*x10*y-357357000*hx6*x10*y
              +6700443750*hy8*x8*y+6891885000*hx2*hy6*x8*y
              +7347809700*hx4*hy4*x8*y+6891885000*hx6*hy2*x8*y
              +6700443750*hx8*x8*y+79659008*a*125*hy9*x7*y
              +79659008*a*125*hx9*x7*y+26801775000*hy10*x6*y
              -80405325000*hx2*hy8*x6*y-51954210000*hx4*hy6*x6*y
              -51954210000*hx6*hy4*x6*y-80405325000*hx8*hy2*x6*y
              +26801775000*hx10*x6*y-716931072*a*125*hx2*hy9*x5*y
              -716931072*a*125*hx9*hy2*x5*y-22334812500*hy12*x4*y
              -172297125000*hx2*hy10*x4*y+324713812500*hx4*hy8*x4*y
              +196796250000*hx6*hy6*x4*y+324713812500*hx8*hy4*x4*y
              -172297125000*hx10*hy2*x4*y-22334812500*hx12*x4*y
              +386039808*a*625*hx4*hy9*x3*y
              +386039808*a*625*hx9*hy4*x3*y+19144125000*hy14*x2*y
              +57432375000*hx2*hy12*x2*y+278326125000*hx4*hy10*x2*y
              -491990625000*hx6*hy8*x2*y-491990625000*hx8*hy6*x2*y
              +278326125000*hx10*hy4*x2*y+57432375000*hx12*hy2*x2*y
              +19144125000*hx14*x2*y-11698176*a*125*125*hx6*hy9*x*y
              -11698176*a*125*125*hx9*hy6*x*y-8546484375*hy16*y
              -8204625000*hx2*hy14*y-15462562500*hx4*hy12*y
              -70284375000*hx6*hy10*y+204996093750*hx8*hy8*y
              -70284375000*hx10*hy6*y-15462562500*hx12*hy4*y
              -8204625000*hx14*hy2*y-8546484375*hx16*y-1001*x17
              +58344*hy2*x15+58344*hx2*x15-1649340*hy4*x13
              -2625480*hx2*hy2*x13-1649340*hx4*x13+32487000*hy6*x11
              +55135080*hx2*hy4*x11+55135080*hx4*hy2*x11
              +32487000*hx6*x11-744493750*hy8*x9-765765000*hx2*hy6*x9
              -816423300*hx4*hy4*x9-765765000*hx6*hy2*x9
              -744493750*hx8*x9-9957376*a*125*hy9*x8
              -9957376*a*125*hx9*x8-3828825000*hy10*x7
              +11486475000*hx2*hy8*x7+7422030000*hx4*hy6*x7
              +7422030000*hx6*hy4*x7+11486475000*hx8*hy2*x7
              -3828825000*hx10*x7+119488512*a*125*hx2*hy9*x6
              +119488512*a*125*hx9*hy2*x6+4466962500*hy12*x5
              +34459425000*hx2*hy10*x5-64942762500*hx4*hy8*x5
              -39359250000*hx6*hy6*x5-64942762500*hx8*hy4*x5
              +34459425000*hx10*hy2*x5+4466962500*hx12*x5
              -96509952*a*625*hx4*hy9*x4-96509952*a*625*hx9*hy4*x4
              -6381375000*hy14*x3-19144125000*hx2*hy12*x3
              -92775375000*hx4*hy10*x3+163996875000*hx6*hy8*x3
              +163996875000*hx8*hy6*x3-92775375000*hx10*hy4*x3
              -19144125000*hx12*hy2*x3-6381375000*hx14*x3
              +5849088*a*125*125*hx6*hy9*x2+5849088*a*125*125*hx9*hy6*x2
              +8546484375*hy16*x+8204625000*hx2*hy14*x
              +15462562500*hx4*hy12*x+70284375000*hx6*hy10*x
              -204996093750*hx8*hy8*x+70284375000*hx10*hy6*x
              +15462562500*hx12*hy4*x+8204625000*hx14*hy2*x
              +8546484375*hx16*x+28672*a*125*625*hy17
              -487424*a*125*625*hx8*hy9-487424*a*125*625*hx9*hy8
              +28672*a*125*625*hx17)
         /(8912896000000*hx8*hy8));
}

double np_aconvol_epan8_ylessx(const double x, const double y,const double hx,const double hy){
  const double a = sqrt(5);
  const double y2 = y*y;
  const double y3 = y2*y;
  const double y4 = y2*y2;
  const double y5 = y3*y2;
  const double y6 = y3*y3;
  const double y7 = y4*y3;
  const double y8 = y4*y4;
  const double y9 = y5*y4;
  const double y10 = y6*y4;
  const double y11 = y7*y4;
  const double y12 = y8*y4;
  const double y13 = y9*y4;
  const double y14 = y10*y4;
  const double y15 = y11*y4;
  const double y16 = y12*y4;
  const double y17 = y13*y4;
 
  const double x2 = x*x;
  const double x3 = x2*x;
  const double x4 = x2*x2;
  const double x5 = x3*x2;
  const double x6 = x3*x3;
  const double x7 = x4*x3;
  const double x8 = x4*x4;
  const double x9 = x5*x4;
  const double x10 = x6*x4;
  const double x11 = x7*x4;
  const double x12 = x8*x4;
  const double x13 = x9*x4;
  const double x14 = x10*x4;
  const double x15 = x11*x4;
  const double x16 = x12*x4;
  const double x17 = x13*x4;

  const double hx2 = hx*hx;
  const double hx3 = hx2*hx;
  const double hx4 = hx2*hx2;
  const double hx5 = hx3*hx2;
  const double hx6 = hx3*hx3;
  const double hx8 = hx4*hx4;
  const double hx9 = hx5*hx4;
  const double hx10 = hx6*hx4;
  const double hx12 = hx8*hx4;
  const double hx13 = hx9*hx4;
  const double hx14 = hx10*hx4;
  const double hx16 = hx12*hx4;
  const double hx17 = hx13*hx4;

  const double hy2 = hy*hy;
  const double hy3 = hy2*hy;
  const double hy4 = hy2*hy2;
  const double hy5 = hy3*hy2;
  const double hy6 = hy3*hy3;
  const double hy8 = hy4*hy4;
  const double hy9 = hy5*hy4;
  const double hy10 = hy6*hy4;
  const double hy12 = hy8*hy4;
  const double hy13 = hy9*hy4;
  const double hy14 = hy10*hy4;
  const double hy16 = hy12*hy4;
  const double hy17 = hy13*hy4;

  return(63*(1001*y17-17017*x*y16+136136*x2*y15-58344*hy2*y15-58344*hx2*y15
             -680680*x3*y14+875160*hy2*x*y14+875160*hx2*x*y14
             +2382380*x4*y13-6126120*hy2*x2*y13-6126120*hx2*x2*y13
             +1649340*hy4*y13+2625480*hx2*hy2*y13+1649340*hx4*y13
             -6194188*x5*y12+26546520*hy2*x3*y12+26546520*hx2*x3*y12
             -21441420*hy4*x*y12-34131240*hx2*hy2*x*y12
             -21441420*hx4*x*y12+12388376*x6*y11-79639560*hy2*x4*y11
             -79639560*hx2*x4*y11+128648520*hy4*x2*y11
             +204787440*hx2*hy2*x2*y11+128648520*hx4*x2*y11
             -32487000*hy6*y11-55135080*hx2*hy4*y11
             -55135080*hx4*hy2*y11-32487000*hx6*y11-19467448*x7*y10
             +175207032*hy2*x5*y10+175207032*hx2*x5*y10
             -471711240*hy4*x3*y10-750887280*hx2*hy2*x3*y10
             -471711240*hx4*x3*y10+357357000*hy6*x*y10
             +606485880*hx2*hy4*x*y10+606485880*hx4*hy2*x*y10
             +357357000*hx6*x*y10+24334310*x8*y9-292011720*hy2*x6*y9
             -292011720*hx2*x6*y9+1179278100*hy4*x4*y9
             +1877218200*hx2*hy2*x4*y9+1179278100*hx4*x4*y9
             -1786785000*hy6*x2*y9-3032429400*hx2*hy4*x2*y9
             -3032429400*hx4*hy2*x2*y9-1786785000*hx6*x2*y9
             +744493750*hy8*y9+765765000*hx2*hy6*y9
             +816423300*hx4*hy4*y9+765765000*hx6*hy2*y9
             +744493750*hx8*y9-24334310*x9*y8+375443640*hy2*x7*y8
             +375443640*hx2*x7*y8-2122700580*hy4*x5*y8
             -3378992760*hx2*hy2*x5*y8-2122700580*hx4*x5*y8
             +5360355000*hy6*x3*y8+9097288200*hx2*hy4*x3*y8
             +9097288200*hx4*hy2*x3*y8+5360355000*hx6*x3*y8
             -6700443750*hy8*x*y8-6891885000*hx2*hy6*x*y8
             -7347809700*hx4*hy4*x*y8-6891885000*hx6*hy2*x*y8
             -6700443750*hx8*x*y8+9957376*a*125*hy9*y8
             +9957376*a*125*hx9*y8+19467448*x10*y7
             -375443640*hy2*x8*y7-375443640*hx2*x8*y7
             +2830267440*hy4*x6*y7+4505323680*hx2*hy2*x6*y7
             +2830267440*hx4*x6*y7-10720710000*hy6*x4*y7
             -18194576400*hx2*hy4*x4*y7-18194576400*hx4*hy2*x4*y7
             -10720710000*hx6*x4*y7+26801775000*hy8*x2*y7
             +27567540000*hx2*hy6*x2*y7+29391238800*hx4*hy4*x2*y7
             +27567540000*hx6*hy2*x2*y7+26801775000*hx8*x2*y7
             -79659008*a*125*hy9*x*y7-79659008*a*125*hx9*x*y7
             +3828825000*hy10*y7-11486475000*hx2*hy8*y7
             -7422030000*hx4*hy6*y7-7422030000*hx6*hy4*y7
             -11486475000*hx8*hy2*y7+3828825000*hx10*y7-12388376*x11*y6
             +292011720*hy2*x9*y6+292011720*hx2*x9*y6
             -2830267440*hy4*x7*y6-4505323680*hx2*hy2*x7*y6
             -2830267440*hx4*x7*y6+15008994000*hy6*x5*y6
             +25472406960*hx2*hy4*x5*y6+25472406960*hx4*hy2*x5*y6
             +15008994000*hx6*x5*y6-62537475000*hy8*x3*y6
             -64324260000*hx2*hy6*x3*y6-68579557200*hx4*hy4*x3*y6
             -64324260000*hx6*hy2*x3*y6-62537475000*hx8*x3*y6
             +278806528*a*125*hy9*x2*y6+278806528*a*125*hx9*x2*y6
             -26801775000*hy10*x*y6+80405325000*hx2*hy8*x*y6
             +51954210000*hx4*hy6*x*y6+51954210000*hx6*hy4*x*y6
             +80405325000*hx8*hy2*x*y6-26801775000*hx10*x*y6
             -119488512*a*125*hx2*hy9*y6-119488512*a*125*hx9*hy2*y6
             +6194188*x12*y5-175207032*hy2*x10*y5-175207032*hx2*x10*y5
             +2122700580*hy4*x8*y5+3378992760*hx2*hy2*x8*y5
             +2122700580*hx4*x8*y5-15008994000*hy6*x6*y5
             -25472406960*hx2*hy4*x6*y5-25472406960*hx4*hy2*x6*y5
             -15008994000*hx6*x6*y5+93806212500*hy8*x4*y5
             +96486390000*hx2*hy6*x4*y5+102869335800*hx4*hy4*x4*y5
             +96486390000*hx6*hy2*x4*y5+93806212500*hx8*x4*y5
             -557613056*a*125*hy9*x3*y5-557613056*a*125*hx9*x3*y5
             +80405325000*hy10*x2*y5-241215975000*hx2*hy8*x2*y5
             -155862630000*hx4*hy6*x2*y5-155862630000*hx6*hy4*x2*y5
             -241215975000*hx8*hy2*x2*y5+80405325000*hx10*x2*y5
             +716931072*a*125*hx2*hy9*x*y5
             +716931072*a*125*hx9*hy2*x*y5-4466962500*hy12*y5
             -34459425000*hx2*hy10*y5+64942762500*hx4*hy8*y5
             +39359250000*hx6*hy6*y5+64942762500*hx8*hy4*y5
             -34459425000*hx10*hy2*y5-4466962500*hx12*y5-2382380*x13*y4
             +79639560*hy2*x11*y4+79639560*hx2*x11*y4
             -1179278100*hy4*x9*y4-1877218200*hx2*hy2*x9*y4
             -1179278100*hx4*x9*y4+10720710000*hy6*x7*y4
             +18194576400*hx2*hy4*x7*y4+18194576400*hx4*hy2*x7*y4
             +10720710000*hx6*x7*y4-93806212500*hy8*x5*y4
             -96486390000*hx2*hy6*x5*y4-102869335800*hx4*hy4*x5*y4
             -96486390000*hx6*hy2*x5*y4-93806212500*hx8*x5*y4
             +139403264*a*625*hy9*x4*y4+139403264*a*625*hx9*x4*y4
             -134008875000*hy10*x3*y4+402026625000*hx2*hy8*x3*y4
             +259771050000*hx4*hy6*x3*y4+259771050000*hx6*hy4*x3*y4
             +402026625000*hx8*hy2*x3*y4-134008875000*hx10*x3*y4
             -358465536*a*625*hx2*hy9*x2*y4
             -358465536*a*625*hx9*hy2*x2*y4+22334812500*hy12*x*y4
             +172297125000*hx2*hy10*x*y4-324713812500*hx4*hy8*x*y4
             -196796250000*hx6*hy6*x*y4-324713812500*hx8*hy4*x*y4
             +172297125000*hx10*hy2*x*y4+22334812500*hx12*x*y4
             +96509952*a*625*hx4*hy9*y4+96509952*a*625*hx9*hy4*y4
             +680680*x14*y3-26546520*hy2*x12*y3-26546520*hx2*x12*y3
             +471711240*hy4*x10*y3+750887280*hx2*hy2*x10*y3
             +471711240*hx4*x10*y3-5360355000*hy6*x8*y3
             -9097288200*hx2*hy4*x8*y3-9097288200*hx4*hy2*x8*y3
             -5360355000*hx6*x8*y3+62537475000*hy8*x6*y3
             +64324260000*hx2*hy6*x6*y3+68579557200*hx4*hy4*x6*y3
             +64324260000*hx6*hy2*x6*y3+62537475000*hx8*x6*y3
             -557613056*a*125*hy9*x5*y3-557613056*a*125*hx9*x5*y3
             +134008875000*hy10*x4*y3-402026625000*hx2*hy8*x4*y3
             -259771050000*hx4*hy6*x4*y3-259771050000*hx6*hy4*x4*y3
             -402026625000*hx8*hy2*x4*y3+134008875000*hx10*x4*y3
             +477954048*a*625*hx2*hy9*x3*y3
             +477954048*a*625*hx9*hy2*x3*y3-44669625000*hy12*x2*y3
             -344594250000*hx2*hy10*x2*y3+649427625000*hx4*hy8*x2*y3
             +393592500000*hx6*hy6*x2*y3+649427625000*hx8*hy4*x2*y3
             -344594250000*hx10*hy2*x2*y3-44669625000*hx12*x2*y3
             -386039808*a*625*hx4*hy9*x*y3
             -386039808*a*625*hx9*hy4*x*y3+6381375000*hy14*y3
             +19144125000*hx2*hy12*y3+92775375000*hx4*hy10*y3
             -163996875000*hx6*hy8*y3-163996875000*hx8*hy6*y3
             +92775375000*hx10*hy4*y3+19144125000*hx12*hy2*y3
             +6381375000*hx14*y3-136136*x15*y2+6126120*hy2*x13*y2
             +6126120*hx2*x13*y2-128648520*hy4*x11*y2
             -204787440*hx2*hy2*x11*y2-128648520*hx4*x11*y2
             +1786785000*hy6*x9*y2+3032429400*hx2*hy4*x9*y2
             +3032429400*hx4*hy2*x9*y2+1786785000*hx6*x9*y2
             -26801775000*hy8*x7*y2-27567540000*hx2*hy6*x7*y2
             -29391238800*hx4*hy4*x7*y2-27567540000*hx6*hy2*x7*y2
             -26801775000*hx8*x7*y2+278806528*a*125*hy9*x6*y2
             +278806528*a*125*hx9*x6*y2-80405325000*hy10*x5*y2
             +241215975000*hx2*hy8*x5*y2+155862630000*hx4*hy6*x5*y2
             +155862630000*hx6*hy4*x5*y2+241215975000*hx8*hy2*x5*y2
             -80405325000*hx10*x5*y2-358465536*a*625*hx2*hy9*x4*y2
             -358465536*a*625*hx9*hy2*x4*y2+44669625000*hy12*x3*y2
             +344594250000*hx2*hy10*x3*y2-649427625000*hx4*hy8*x3*y2
             -393592500000*hx6*hy6*x3*y2-649427625000*hx8*hy4*x3*y2
             +344594250000*hx10*hy2*x3*y2+44669625000*hx12*x3*y2
             +579059712*a*625*hx4*hy9*x2*y2
             +579059712*a*625*hx9*hy4*x2*y2-19144125000*hy14*x*y2
             -57432375000*hx2*hy12*x*y2-278326125000*hx4*hy10*x*y2
             +491990625000*hx6*hy8*x*y2+491990625000*hx8*hy6*x*y2
             -278326125000*hx10*hy4*x*y2-57432375000*hx12*hy2*x*y2
             -19144125000*hx14*x*y2-5849088*a*125*125*hx6*hy9*y2
             -5849088*a*125*125*hx9*hy6*y2+17017*x16*y-875160*hy2*x14*y
             -875160*hx2*x14*y+21441420*hy4*x12*y
             +34131240*hx2*hy2*x12*y+21441420*hx4*x12*y
             -357357000*hy6*x10*y-606485880*hx2*hy4*x10*y
             -606485880*hx4*hy2*x10*y-357357000*hx6*x10*y
             +6700443750*hy8*x8*y+6891885000*hx2*hy6*x8*y
             +7347809700*hx4*hy4*x8*y+6891885000*hx6*hy2*x8*y
             +6700443750*hx8*x8*y-79659008*a*125*hy9*x7*y
             -79659008*a*125*hx9*x7*y+26801775000*hy10*x6*y
             -80405325000*hx2*hy8*x6*y-51954210000*hx4*hy6*x6*y
             -51954210000*hx6*hy4*x6*y-80405325000*hx8*hy2*x6*y
             +26801775000*hx10*x6*y+716931072*a*125*hx2*hy9*x5*y
             +716931072*a*125*hx9*hy2*x5*y-22334812500*hy12*x4*y
             -172297125000*hx2*hy10*x4*y+324713812500*hx4*hy8*x4*y
             +196796250000*hx6*hy6*x4*y+324713812500*hx8*hy4*x4*y
             -172297125000*hx10*hy2*x4*y-22334812500*hx12*x4*y
             -386039808*a*625*hx4*hy9*x3*y
             -386039808*a*625*hx9*hy4*x3*y+19144125000*hy14*x2*y
             +57432375000*hx2*hy12*x2*y+278326125000*hx4*hy10*x2*y
             -491990625000*hx6*hy8*x2*y-491990625000*hx8*hy6*x2*y
             +278326125000*hx10*hy4*x2*y+57432375000*hx12*hy2*x2*y
             +19144125000*hx14*x2*y+11698176*a*125*125*hx6*hy9*x*y
             +11698176*a*125*125*hx9*hy6*x*y-8546484375*hy16*y
             -8204625000*hx2*hy14*y-15462562500*hx4*hy12*y
             -70284375000*hx6*hy10*y+204996093750*hx8*hy8*y
             -70284375000*hx10*hy6*y-15462562500*hx12*hy4*y
             -8204625000*hx14*hy2*y-8546484375*hx16*y-1001*x17
             +58344*hy2*x15+58344*hx2*x15-1649340*hy4*x13
             -2625480*hx2*hy2*x13-1649340*hx4*x13+32487000*hy6*x11
             +55135080*hx2*hy4*x11+55135080*hx4*hy2*x11
             +32487000*hx6*x11-744493750*hy8*x9-765765000*hx2*hy6*x9
             -816423300*hx4*hy4*x9-765765000*hx6*hy2*x9
             -744493750*hx8*x9+9957376*a*125*hy9*x8
             +9957376*a*125*hx9*x8-3828825000*hy10*x7
             +11486475000*hx2*hy8*x7+7422030000*hx4*hy6*x7
             +7422030000*hx6*hy4*x7+11486475000*hx8*hy2*x7
             -3828825000*hx10*x7-119488512*a*125*hx2*hy9*x6
             -119488512*a*125*hx9*hy2*x6+4466962500*hy12*x5
             +34459425000*hx2*hy10*x5-64942762500*hx4*hy8*x5
             -39359250000*hx6*hy6*x5-64942762500*hx8*hy4*x5
             +34459425000*hx10*hy2*x5+4466962500*hx12*x5
             +96509952*a*625*hx4*hy9*x4+96509952*a*625*hx9*hy4*x4
             -6381375000*hy14*x3-19144125000*hx2*hy12*x3
             -92775375000*hx4*hy10*x3+163996875000*hx6*hy8*x3
             +163996875000*hx8*hy6*x3-92775375000*hx10*hy4*x3
             -19144125000*hx12*hy2*x3-6381375000*hx14*x3
             -5849088*a*125*125*hx6*hy9*x2-5849088*a*125*125*hx9*hy6*x2
             +8546484375*hy16*x+8204625000*hx2*hy14*x
             +15462562500*hx4*hy12*x+70284375000*hx6*hy10*x
             -204996093750*hx8*hy8*x+70284375000*hx10*hy6*x
             +15462562500*hx12*hy4*x+8204625000*hx14*hy2*x
             +8546484375*hx16*x-28672*a*625*125*hy17
             +487424*a*625*125*hx8*hy9+487424*a*625*125*hx9*hy8
             -28672*a*625*125*hx17)
         /(8912896000000*hx8*hy8));
}

double np_aconvol_epan8(const double x, const double y,const double hx,const double hy){
  const double a = sqrt(5.0);
  const double dxy = fabs(x-y);

  if(dxy >= a*(hx+hy)){
    return 0;
  } else if(dxy > a*fabs(hx-hy)){
    if(x<y)
      return (np_aconvol_epan8_xlessy(x,y,hx,hy));
    else
      return (np_aconvol_epan8_ylessx(x,y,hx,hy));
  } else {
    return (np_aconvol_epan8_total(x,y,hx,hy));
  }
}

double np_aconvol_rect(const double x, const double y,const double hx,const double hy){
  return (fabs(x-y) >= (hx+hy)) ? 0.0 : 0.25/(hx*hy)*(MIN(x+hx,y+hy) - MAX(x-hx,y-hy));
}

double np_aconvol_tgauss2_total(const double x, const double y,const double hx,const double hy){
  const double x2 = x*x;
  const double y2 = y*y;

  const double hx2 = hx*hx;
  const double hx4 = hx2*hx2;

  const double hy2 = hy*hy;
  const double hy4 = hy2*hy2;

  const double a = sqrt(2);
  const double b = sqrt(M_PI);
  const double c = sqrt(hy2+hx2);

  return(exp(-y2/(2*hy2)-x2/(2*hx2)-9)*
         (b*hx*hy
          *exp(hx2*y2/(2*hy4+2*hx2*hy2)
               +x*y/(hy2+hx2)
               +hy2*x2/(2*hx2*hy2+2*hx4)
               +9)
          *erfun((hx*y-hx*x+(hy2+hx2)*np_tgauss2_b)
                 /(a*hy*c))
          -b*hx*hy
          *exp(hx2*y2/(2*hy4+2*hx2*hy2)
               +x*y/(hy2+hx2)
               +hy2*x2/(2*hx2*hy2+2*hx4)
               +9)
          *erfun((hx*y-hx*x+(-hy2-hx2)*np_tgauss2_b)
                 /(a*hy*c))
          -b*hy*c
          *exp(y2/(2*hy2)+x2/(2*hx2)+9/2)
          *erfun((y-x+hx*np_tgauss2_b)/(a*hy))
          +b*hy*c
          *exp(y2/(2*hy2)+x2/(2*hx2)+9/2)
          *erfun((y-x-hx*np_tgauss2_b)/(a*hy))
          -2*b*hx*c
          *erfun(np_tgauss2_b/a)
          *exp(y2/(2*hy2)+x2/(2*hx2)+9/2)
          +a*2*hx*c*np_tgauss2_b
          *exp(y2/(2*hy2)+x2/(2*hx2)))
         /(a*2*M_PI*c*np_tgauss2_alpha*np_tgauss2_alpha));
}

double np_aconvol_tgauss2_indefinite(const double u, const double x, const double y,const double hx,const double hy){
  const double x2 = x*x;
  const double y2 = y*y;

  const double hx2 = hx*hx;
  const double hx4 = hx2*hx2;

  const double hy2 = hy*hy;
  const double hy4 = hy2*hy2;

  const double a = sqrt(2);
  const double b = sqrt(M_PI);
  const double c = sqrt(hy2+hx2);

  return(-exp(-y2/(2*hy2)-x2/(2*hx2)-19/2)*
         (b*hx*hy
          *exp(hx2*y2
            /(2*hy4+2*hx2*hy2)
            +x*y/(hy2+hx2)
            +hy2*x2
            /(2*hx2*hy2+2*hx4)+19/2)
          *erfun(
               (hx2*y+hy2*x
                +(-hy2-hx2)*u)
               /(a*hx*hy
                 *c))
          -b*hy*c
          *exp(y2/(2*hy2)+x2/(2*hx2)
            +5)
          *erfun((y-u)/(a*hy))
          -b*hx*c
          *erfun((x-u)/(a*hx))
          *exp(y2/(2*hy2)+x2/(2*hx2)
            +5)
          -a*c*u
          *exp(y2/(2*hy2)+x2/(2*hx2)
            +1/2))
         /(a*2*M_PI*c*np_tgauss2_alpha*np_tgauss2_alpha));
}

double np_aconvol_tgauss2(const double x, const double y,const double hx,const double hy){
  const double a = np_tgauss2_b;
  const double dxy = fabs(x-y);

  if(dxy >= a*(hx+hy)){
    return 0;
  } else if(dxy > a*fabs(hx-hy)){
    return (np_aconvol_tgauss2_indefinite(MIN(x+a*hx,y+a*hy),x,y,hx,hy) - 
            np_aconvol_tgauss2_indefinite(MAX(x-a*hx,y-a*hy),x,y,hx,hy));
  } else {
    return (np_aconvol_tgauss2_total(x,y,hx,hy));
  }
}
// end kernels

double (* const allck[])(double) = { np_gauss2, np_gauss4, np_gauss6, np_gauss8, 
                                     np_epan2, np_epan4, np_epan6, np_epan8, 
                                     np_rect, np_tgauss2 };
double (* const allok[])(double, double, double, double, double) = { np_owang_van_ryzin, np_oli_racine };
double (* const alluk[])(int, double, int) = { np_uaa, np_unli_racine };

// in cksup we define a scale length for all kernels, outside of which the kernel is 0
// for fixed bandwidths and finite support kernels, only points within += 1 scale length
// contribute at a given point of evaluation
// when constructing the search box, in 1D the size of the box is
// x_eval + [-RIGHT_SUPPORT,-LEFT_SUPPORT], the interval is reversed and negated because you want to know if 
// the evaluation point is in the support of other points, not vice versa :)
// of course in the case of adaptive bandwidths the situation is opposite: x_train + [LEFT_SUPPORT,RIGHT_SUPPORT]

// as it stands, CDF's are only partially accelerated by the tree, but that is fixable
#define SQRT5 2.23606797749979
#define SQRT20 4.47213595499958

double cksup[OP_NCFUN][2] = { {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX}, 
                              {-SQRT5, SQRT5}, {-SQRT5, SQRT5}, {-SQRT5, SQRT5}, {-SQRT5, SQRT5},
                              {-1.0, 1.0}, {-3.0, 3.0},
                              {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX},
                              {-SQRT20, SQRT20},{-SQRT20, SQRT20},{-SQRT20, SQRT20},{-SQRT20, SQRT20},
                              {-2.0, 2.0}, {-6.0, 6.0},
                              {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX},
                              {-SQRT5, SQRT5}, {-SQRT5, SQRT5}, {-SQRT5, SQRT5}, {-SQRT5, SQRT5},
                              {-0.0, 0.0}, // PLEASE DON'T EVER USE THIS
                              {-3.0, 3.0},
                              {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX}, {-DBL_MAX, DBL_MAX}, 
                              {-SQRT5, DBL_MAX}, {-SQRT5, DBL_MAX}, {-SQRT5, DBL_MAX}, {-SQRT5, DBL_MAX}, 
                              {-1.0, DBL_MAX},
                              {-3.0, DBL_MAX} };

/* 
   np_kernelv does weighted products of vectors - this is useful for 
   product kernels, where each kernel in each dimension acts as a weight.
*/

/* xt = training data */
/* xw = x weights */

void np_p_ckernelv(const int KERNEL, 
                   const int P_KERNEL,
                   const int P_IDX,
                   const int P_NIDX,
                   const double * const xt, const int num_xt, 
                   const int do_xw,
                   const double x, const double h, 
                   double * const result,
                   double * const p_result,
                   const XL * const xl,
                   const XL * const p_xl,
                   const int swap_xxt,
                   const int do_perm,
                   const int do_score,
                   double * const scratch_kbuf,
                   const int use_largeh_pre,
                   const double largeh_k0_pre,
                   const double invnorm,
                   const double p_invnorm){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i,j,r,l; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  const double sgn = swap_xxt ? -1.0 : 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double * const pxw = (bin_do_xw ? p_result : &unit_weight);

  double (* const k[])(double) = { np_gauss2, np_gauss4, np_gauss6, np_gauss8, //ordinary kernels
                                   np_epan2, np_epan4, np_epan6, np_epan8, 
                                   np_rect, np_tgauss2, 
                                   np_econvol_gauss2, np_econvol_gauss4, np_econvol_gauss6, np_econvol_gauss8, // convolution kernels
                                   np_econvol_epan2, np_econvol_epan4, np_econvol_epan6, np_econvol_epan8,
                                   np_econvol_rect, np_econvol_tgauss2,
                                   np_deriv_gauss2, np_deriv_gauss4, np_deriv_gauss6, np_deriv_gauss8, // derivative kernels
                                   np_deriv_epan2, np_deriv_epan4, np_deriv_epan6, np_deriv_epan8, 
                                   np_deriv_rect, np_deriv_tgauss2,
                                   np_cdf_gauss2, np_cdf_gauss4, np_cdf_gauss6, np_cdf_gauss8, // cdfative kernels
                                   np_cdf_epan2, np_cdf_epan4, np_cdf_epan6, np_cdf_epan8, 
                                   np_cdf_rect, np_cdf_tgauss2 };

  double *kbuf = scratch_kbuf;
  const int own_kbuf = (kbuf == NULL);
  if(own_kbuf){
    kbuf = (double *)malloc(num_xt*sizeof(double));
    if(kbuf == NULL) error("memory allocation failed");
  }

  /* Base-kernel large-h decision is precomputed by caller once per (j,i). */
  const int use_largeh = use_largeh_pre;
  const int use_largeh_perm = use_largeh &&
    (!do_score) &&
    do_perm &&
    np_cont_largeh_kernel_supported(P_KERNEL) &&
    np_cont_largeh_is_active(P_KERNEL, xt, num_xt, x, h, p_xl);

  if(use_largeh){
    const double kn = largeh_k0_pre*invnorm;
    np_ckernelv_mul_const(kn, num_xt, do_xw, result, xl);

    if(do_perm){
      if(use_largeh_perm){
        const double pkn = np_cont_largeh_k0(P_KERNEL)*p_invnorm;
        np_ckernelv_mul_const(pkn, num_xt, do_xw, p_result + P_IDX*num_xt, p_xl);
      } else {
        if(p_xl == NULL){
          for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
            p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*p_invnorm*k[P_KERNEL]((x-xt[i])*sgn/h)*(do_score ? ((xt[i]-x)*sgn/h) : 1.0);
          }
        } else {
          for (int m = 0; m < p_xl->n; m++){
            const int istart = p_xl->istart[m];
            const int nlev = p_xl->nlev[m];
            for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
              p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*p_invnorm*k[P_KERNEL]((x-xt[i])*sgn/h)*(do_score ? ((xt[i]-x)*sgn/h) : 1.0);
            }
          }
        }
      }
    }

    for(l = 0, r = 0; l < P_NIDX; l++, r += bin_do_xw){
      if((l == P_IDX) && do_perm) continue;
      np_ckernelv_mul_const(kn, num_xt, bin_do_xw, p_result + l*num_xt, xl);
    }

    if(own_kbuf)
      free(kbuf);
    return;
  }

  if(xl == NULL){
    for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
      const double kn = invnorm*k[KERNEL]((x-xt[i])*sgn/h);

      result[i] = xw[j]*kn;
      kbuf[i] = kn;
      
      if(do_perm)
        p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*p_invnorm*k[P_KERNEL]((x-xt[i])*sgn/h)*(do_score ? ((xt[i]-x)*sgn/h) : 1.0);

    }

    for(l = 0, r = 0; l < P_NIDX; l++, r += bin_do_xw){
      if((l == P_IDX) && do_perm) continue;
      for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
        p_result[l*num_xt + i] = pxw[r*num_xt + j]*kbuf[i];
      }
    }

  }
  else{
    for (int m = 0; m < xl->n; m++){
      const int istart = xl->istart[m];
      const int nlev = xl->nlev[m];
      for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
        const double kn = invnorm*k[KERNEL]((x-xt[i])*sgn/h);

        result[i] = xw[j]*kn;
        kbuf[i] = kn;
      }
    }

    if(do_perm){
      for (int m = 0; m < p_xl->n; m++){
        const int istart = p_xl->istart[m];
        const int nlev = p_xl->nlev[m];
        for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
          p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*p_invnorm*k[P_KERNEL]((x-xt[i])*sgn/h)*(do_score ? ((xt[i]-x)*sgn/h) : 1.0);
        }
      }
    }

    for(l = 0, r = 0; l < P_NIDX; l++, r+=bin_do_xw){
      if((l == P_IDX) && do_perm) continue;
      for (int m = 0; m < xl->n; m++){
        const int istart = xl->istart[m];
        const int nlev = xl->nlev[m];
        for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
          p_result[l*num_xt + i] = pxw[r*num_xt + j]*kbuf[i];
        }
      }
    }
  }

  if(own_kbuf)
    free(kbuf);
}

void np_ckernelv(const int KERNEL, 
                 const double * const xt, const int num_xt, 
                 const int do_xw,
                 const double x, const double h, 
                 double * const result,
                 const XL * const xl,
                 const int swap_xxt,
                 const int skip_largeh_check,
                 const double invnorm){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  const double sgn = swap_xxt ? -1.0 : 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  if((!skip_largeh_check) && np_cont_largeh_is_active(KERNEL, xt, num_xt, x, h, xl)){
    np_ckernelv_mul_const(np_cont_largeh_k0(KERNEL)*invnorm, num_xt, do_xw, result, xl);
    return;
  }

  /*
    Hot path: avoid indirect function-pointer calls and avoid branching on
    do_xw/xl inside the innermost kernel loop. This improves CPU efficiency
    without changing the numerical kernel functions.
  */

#define NP_CKERNELV_APPLY(fn)                                                      \
  do {                                                                             \
    if(xl == NULL){                                                                \
      if(!bin_do_xw){                                                              \
        for(i = 0; i < num_xt; i++){                                               \
          result[i] = invnorm*fn((x-xt[i])*sgn/h);                                 \
        }                                                                          \
      } else {                                                                     \
        for(i = 0; i < num_xt; i++){                                               \
          if(xw[i] == 0.0) continue;                                               \
          result[i] = xw[i]*invnorm*fn((x-xt[i])*sgn/h);                           \
        }                                                                          \
      }                                                                            \
    } else {                                                                       \
      if(!bin_do_xw){                                                              \
        for(int m = 0; m < xl->n; m++){                                            \
          const int istart = xl->istart[m];                                        \
          const int nlev = xl->nlev[m];                                            \
          for(i = istart; i < istart+nlev; i++){                                   \
            result[i] = invnorm*fn((x-xt[i])*sgn/h);                               \
          }                                                                        \
        }                                                                          \
      } else {                                                                     \
        for(int m = 0; m < xl->n; m++){                                            \
          const int istart = xl->istart[m];                                        \
          const int nlev = xl->nlev[m];                                            \
          for(i = istart; i < istart+nlev; i++){                                   \
            if(xw[i] == 0.0) continue;                                             \
            result[i] = xw[i]*invnorm*fn((x-xt[i])*sgn/h);                         \
          }                                                                        \
        }                                                                          \
      }                                                                            \
    }                                                                              \
  } while(0)

  switch(KERNEL){
    case 0: NP_CKERNELV_APPLY(np_gauss2); break;
    case 1: NP_CKERNELV_APPLY(np_gauss4); break;
    case 2: NP_CKERNELV_APPLY(np_gauss6); break;
    case 3: NP_CKERNELV_APPLY(np_gauss8); break;

    case 4: NP_CKERNELV_APPLY(np_epan2); break;
    case 5: NP_CKERNELV_APPLY(np_epan4); break;
    case 6: NP_CKERNELV_APPLY(np_epan6); break;
    case 7: NP_CKERNELV_APPLY(np_epan8); break;

    case 8: NP_CKERNELV_APPLY(np_rect); break;
    case 9: NP_CKERNELV_APPLY(np_tgauss2); break;

    case 10: NP_CKERNELV_APPLY(np_econvol_gauss2); break;
    case 11: NP_CKERNELV_APPLY(np_econvol_gauss4); break;
    case 12: NP_CKERNELV_APPLY(np_econvol_gauss6); break;
    case 13: NP_CKERNELV_APPLY(np_econvol_gauss8); break;

    case 14: NP_CKERNELV_APPLY(np_econvol_epan2); break;
    case 15: NP_CKERNELV_APPLY(np_econvol_epan4); break;
    case 16: NP_CKERNELV_APPLY(np_econvol_epan6); break;
    case 17: NP_CKERNELV_APPLY(np_econvol_epan8); break;

    case 18: NP_CKERNELV_APPLY(np_econvol_rect); break;
    case 19: NP_CKERNELV_APPLY(np_econvol_tgauss2); break;

    case 20: NP_CKERNELV_APPLY(np_deriv_gauss2); break;
    case 21: NP_CKERNELV_APPLY(np_deriv_gauss4); break;
    case 22: NP_CKERNELV_APPLY(np_deriv_gauss6); break;
    case 23: NP_CKERNELV_APPLY(np_deriv_gauss8); break;

    case 24: NP_CKERNELV_APPLY(np_deriv_epan2); break;
    case 25: NP_CKERNELV_APPLY(np_deriv_epan4); break;
    case 26: NP_CKERNELV_APPLY(np_deriv_epan6); break;
    case 27: NP_CKERNELV_APPLY(np_deriv_epan8); break;

    case 28: NP_CKERNELV_APPLY(np_deriv_rect); break;
    case 29: NP_CKERNELV_APPLY(np_deriv_tgauss2); break;

    case 30: NP_CKERNELV_APPLY(np_cdf_gauss2); break;
    case 31: NP_CKERNELV_APPLY(np_cdf_gauss4); break;
    case 32: NP_CKERNELV_APPLY(np_cdf_gauss6); break;
    case 33: NP_CKERNELV_APPLY(np_cdf_gauss8); break;

    case 34: NP_CKERNELV_APPLY(np_cdf_epan2); break;
    case 35: NP_CKERNELV_APPLY(np_cdf_epan4); break;
    case 36: NP_CKERNELV_APPLY(np_cdf_epan6); break;
    case 37: NP_CKERNELV_APPLY(np_cdf_epan8); break;

    case 38: NP_CKERNELV_APPLY(np_cdf_rect); break;
    case 39: NP_CKERNELV_APPLY(np_cdf_tgauss2); break;

    default: {
      // Defensive: preserve behavior for unexpected kernel codes.
      double (* const k[])(double) = { np_gauss2, np_gauss4, np_gauss6, np_gauss8,
                                       np_epan2, np_epan4, np_epan6, np_epan8,
                                       np_rect, np_tgauss2,
                                       np_econvol_gauss2, np_econvol_gauss4, np_econvol_gauss6, np_econvol_gauss8,
                                       np_econvol_epan2, np_econvol_epan4, np_econvol_epan6, np_econvol_epan8,
                                       np_econvol_rect, np_econvol_tgauss2,
                                       np_deriv_gauss2, np_deriv_gauss4, np_deriv_gauss6, np_deriv_gauss8,
                                       np_deriv_epan2, np_deriv_epan4, np_deriv_epan6, np_deriv_epan8,
                                       np_deriv_rect, np_deriv_tgauss2,
                                       np_cdf_gauss2, np_cdf_gauss4, np_cdf_gauss6, np_cdf_gauss8,
                                       np_cdf_epan2, np_cdf_epan4, np_cdf_epan6, np_cdf_epan8,
                                       np_cdf_rect, np_cdf_tgauss2 };
      const int kernel = (KERNEL >= 0 && KERNEL < (int)(sizeof(k)/sizeof(k[0]))) ? KERNEL : 0;

      if(xl == NULL){
        if(!bin_do_xw){
          for(i = 0; i < num_xt; i++)
            result[i] = invnorm*k[kernel]((x-xt[i])*sgn/h);
        } else {
          for(i = 0; i < num_xt; i++){
            if(xw[i] == 0.0) continue;
            result[i] = xw[i]*invnorm*k[kernel]((x-xt[i])*sgn/h);
          }
        }
      } else {
        for(int m = 0; m < xl->n; m++){
          const int istart = xl->istart[m];
          const int nlev = xl->nlev[m];
          if(!bin_do_xw){
            for(i = istart; i < istart+nlev; i++)
              result[i] = invnorm*k[kernel]((x-xt[i])*sgn/h);
          } else {
            for(i = istart; i < istart+nlev; i++){
              if(xw[i] == 0.0) continue;
              result[i] = xw[i]*invnorm*k[kernel]((x-xt[i])*sgn/h);
            }
          }
        }
      }
    }
  }

#undef NP_CKERNELV_APPLY

}

void np_convol_ckernelv(const int KERNEL, 
                        const double * const xt, const int num_xt, 
                        const int do_xw,
                        const double x, 
                        double * xt_h,
                        const int xt_h_is_scalar,
                        const double h, 
                        double * const result,
                        const int power){

  int i,j; 
  const int bin_do_xw = do_xw > 0;

  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double (* const k[])(double,double,double,double) = { 
    np_aconvol_gauss2, np_aconvol_gauss4, np_aconvol_gauss6, np_aconvol_gauss8,
    np_aconvol_epan2, np_aconvol_epan4, np_aconvol_epan6, np_aconvol_epan8,
    np_aconvol_rect, np_aconvol_tgauss2
  };

  for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
    double kval;
    const double hy = xt_h_is_scalar ? xt_h[0] : xt_h[i];

    if(xw[j] == 0.0) continue;

    kval = k[KERNEL](x, xt[i], h, hy);

    result[i] = xw[j]*kval/ipow(hy, power);
  }
}


void np_p_ukernelv(const int KERNEL, 
                   const int P_KERNEL,
                   const int P_IDX,
                   const int P_NIDX,
                   const double * const xt, const int num_xt, 
                   const int do_xw,
                   const double x, const double lambda, const int ncat,
                   const double cat,
                   double * const result,
                   double * const p_result,
                   const XL * const xl,
                   const XL * const p_xl,
                   const int swap_xxt,
                   const int do_ocg,
                   double * const scratch_kbuf){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i,j,r,l; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double * const pxw = (bin_do_xw ? p_result : &unit_weight);

  const double ex = do_ocg ? cat : x;

  double (* const k[])(int, double, int) = { np_uaa, np_unli_racine,
                                             np_econvol_uaa, np_econvol_unli_racine,
                                             np_score_uaa, np_score_unli_racine };
  const int nk = (int)(sizeof(k)/sizeof(k[0]));
  const int kernel = (KERNEL >= 0 && KERNEL < nk) ? KERNEL : 0;
  const int p_kernel = (P_KERNEL >= 0 && P_KERNEL < nk) ? P_KERNEL : 0;

  /* Unordered kernels depend only on same/different category; cache both values once. */
  const double kn_same = k[kernel](1, lambda, ncat);
  const double kn_diff = k[kernel](0, lambda, ncat);
  const double pkn_same = k[p_kernel](1, lambda, ncat);
  const double pkn_diff = k[p_kernel](0, lambda, ncat);
  const int use_const_k = np_disc_near_upper(kernel, lambda, ncat) &&
    np_disc_near_const_kernel(kn_same, kn_diff);
  const int use_const_pk = np_disc_near_upper(p_kernel, lambda, ncat) &&
    np_disc_near_const_kernel(pkn_same, pkn_diff);
  const double kn_const = 0.5*(kn_same + kn_diff);
  const double pkn_const = 0.5*(pkn_same + pkn_diff);
  const int p_iscat_const = (swap_xxt && do_ocg);
  const int p_iscat_const_val = (cat == x);

  double *kbuf = scratch_kbuf;
  const int own_kbuf = (kbuf == NULL);
  if(own_kbuf){
    kbuf = (double *)malloc(num_xt*sizeof(double));
    if(kbuf == NULL) error("memory allocation failed");
  }

  if(xl == NULL){
    for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
      const double kn = use_const_k ?
        kn_const :
        (((xt[i] == x) ? kn_same : kn_diff));

      result[i] = xw[j]*kn;
      kbuf[i] = kn;

      const double pkn = use_const_pk ?
        pkn_const :
        (((p_iscat_const ? p_iscat_const_val : (xt[i] == ex)) ? pkn_same : pkn_diff));
      p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*pkn;
    }

    for(l = 0, r = 0; l < P_NIDX; l++, r += bin_do_xw){
      if(l == P_IDX) continue;
      for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
        p_result[l*num_xt + i] = pxw[r*num_xt + j]*kbuf[i];
      }
    }

  }
  else{
    for (int m = 0; m < xl->n; m++){
      const int istart = xl->istart[m];
      const int nlev = xl->nlev[m];
      for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
        const double kn = use_const_k ?
          kn_const :
          (((xt[i] == x) ? kn_same : kn_diff));

        result[i] = xw[j]*kn;
        kbuf[i] = kn;
      }
    }

    for (int m = 0; m < p_xl->n; m++){
      const int istart = p_xl->istart[m];
      const int nlev = p_xl->nlev[m];
      for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
        const double pkn = use_const_pk ?
          pkn_const :
          (((p_iscat_const ? p_iscat_const_val : (xt[i] == ex)) ? pkn_same : pkn_diff));
        p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*pkn;
      }
    }


    for(l = 0, r = 0; l < P_NIDX; l++, r+=bin_do_xw){
      if(l == P_IDX) continue;
      for (int m = 0; m < xl->n; m++){
        const int istart = xl->istart[m];
        const int nlev = xl->nlev[m];
        for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
          p_result[l*num_xt + i] = pxw[r*num_xt + j]*kbuf[i];
        }
      }
    }
  }

  if(own_kbuf)
    free(kbuf);
}

void np_ukernelv(const int KERNEL, 
                 const double * const xt, const int num_xt, 
                 const int do_xw,
                 const double x, const double lambda, const int ncat,
                 double * const result,
                 const XL * const xl,
                 const int skip_upper_gate){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);
  double (* const k[])(int, double, int) = { np_uaa, np_unli_racine,
                                             np_econvol_uaa, np_econvol_unli_racine };
  const int nk = (int)(sizeof(k)/sizeof(k[0]));
  const int kernel = (KERNEL >= 0 && KERNEL < nk) ? KERNEL : 0;
  const double kn_same = k[kernel](1, lambda, ncat);
  const double kn_diff = k[kernel](0, lambda, ncat);
  const int use_const_k = (!skip_upper_gate) &&
    np_disc_near_upper(kernel, lambda, ncat) &&
    np_disc_near_const_kernel(kn_same, kn_diff);
  const double kn_const = 0.5*(kn_same + kn_diff);

  if(xl == NULL){
    if(!bin_do_xw){
      for(i = 0; i < num_xt; i++){
        result[i] = use_const_k ? kn_const : ((xt[i] == x) ? kn_same : kn_diff);
      }
    } else {
      for(i = 0; i < num_xt; i++){
        if(xw[i] == 0.0) continue;
        result[i] = xw[i]*(use_const_k ? kn_const : ((xt[i] == x) ? kn_same : kn_diff));
      }
    }
  } else {
    for(int m = 0; m < xl->n; m++){
      const int istart = xl->istart[m];
      const int nlev = xl->nlev[m];
      if(!bin_do_xw){
        for(i = istart; i < istart+nlev; i++){
          result[i] = use_const_k ? kn_const : ((xt[i] == x) ? kn_same : kn_diff);
        }
      } else {
        for(i = istart; i < istart+nlev; i++){
          if(xw[i] == 0.0) continue;
          result[i] = xw[i]*(use_const_k ? kn_const : ((xt[i] == x) ? kn_same : kn_diff));
        }
      }
    }
  }
}

void np_convol_okernelv(const int KERNEL, 
                        const double * const xt, const int num_xt, 
                        const int do_xw,
                        const double x, const double lambda,
                        int ncat, double * cat,
                        double * const result,
                        const int swap_xxt){

  int i; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  if(!swap_xxt){
    for (i = 0; i < num_xt; i++){
      if(xw[bin_do_xw ? i : 0] == 0.0) continue;
      result[i] = xw[bin_do_xw ? i : 0]*kernel_ordered_convolution(KERNEL, xt[i], x, lambda, ncat, cat);
    }
  } else {
    for (i = 0; i < num_xt; i++){
      if(xw[bin_do_xw ? i : 0] == 0.0) continue;
      result[i] = xw[bin_do_xw ? i : 0]*kernel_ordered_convolution(KERNEL, x, xt[i], lambda, ncat, cat);
    }
  }

}


void np_p_okernelv(const int KERNEL, 
                   const int P_KERNEL,
                   const int P_IDX,
                   const int P_NIDX,
                   const double * const xt, const int num_xt, 
                   const int do_xw,
                   const double x, const double lambda, 
                   const double * cats, const int ncat,
                   double * const result,
                   double * const p_result,
                   const XL * const xl,
                   const XL * const p_xl,
                   const int swap_xxt,
                   const int do_ocg,
                   const int * const ordered_indices,
                   const int swapped_index,
                   double * const scratch_kbuf){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i,j,r,l; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double * const pxw = (bin_do_xw ? p_result : &unit_weight);

  double (* const k[])(double, double, double, double, double) = { 
    np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
    np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
    np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
    np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
  };

  double *kbuf = scratch_kbuf;
  const int own_kbuf = (kbuf == NULL);
  if(own_kbuf){
    kbuf = (double *)malloc(num_xt*sizeof(double));
    if(kbuf == NULL) error("memory allocation failed");
  }

  double s_cat = 0.0;

  if((!swap_xxt) && do_ocg){
    s_cat = cats[abs(swapped_index - 1)];
  }
  
  const double cl = (cats != NULL)? cats[0] : 0.0;
  const double ch = (cats != NULL)? cats[ncat - 1] : 0.0;
  const int max_cxy = (int)fabs(ch-cl);
  const int fast_kernel = (KERNEL >= 0 && KERNEL <= 3 && cats != NULL);
  const int fast_p_kernel = (P_KERNEL >= 0 && P_KERNEL <= 3 && cats != NULL);
  double *lpow = NULL;

  if((fast_kernel || fast_p_kernel) && max_cxy >= 0){
    lpow = (double *)malloc((size_t)(max_cxy+1)*sizeof(double));
    if(lpow == NULL) error("memory allocation failed");
    lpow[0] = 1.0;
    for(int c = 1; c <= max_cxy; c++)
      lpow[c] = lpow[c-1]*lambda;
  }

    if(xl == NULL){
      for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
        const double cat = do_ocg ? (swap_xxt ? cats[abs(ordered_indices[i] - 1)] : s_cat) : 0.0;
        const double c1 = swap_xxt ? x : xt[i];
        const double c2 = swap_xxt ? xt[i] : x;
        const double c3 = do_ocg ? cat : (swap_xxt ? xt[i] : x);

        const double kn = fast_kernel
          ? np_ordered_eval_kernel(KERNEL, c1, c2, lambda, max_cxy, lpow, cats, ncat, cl, ch)
          : k[KERNEL](c1, c2, lambda, cl, ch);

        result[i] = xw[j]*kn;
        kbuf[i] = kn;

        p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*
          (fast_p_kernel
            ? np_ordered_eval_kernel(P_KERNEL, c1, c3, lambda, max_cxy, lpow, cats, ncat, cl, ch)
            : k[P_KERNEL](c1, c3, lambda, cl, ch));
      }

      for(l = 0, r = 0; l < P_NIDX; l++, r += bin_do_xw){
        if(l == P_IDX) continue;
        for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
          p_result[l*num_xt + i] = pxw[r*num_xt + j]*kbuf[i];
        }
      }

    } else {
      for (int m = 0; m < xl->n; m++){
        const int istart = xl->istart[m];
        const int nlev = xl->nlev[m];
        for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
          const double c1 = swap_xxt ? x : xt[i];
          const double c2 = swap_xxt ? xt[i] : x;

          const double kn = fast_kernel
            ? np_ordered_eval_kernel(KERNEL, c1, c2, lambda, max_cxy, lpow, cats, ncat, cl, ch)
            : k[KERNEL](c1, c2, lambda, cl, ch);

          result[i] = xw[j]*kn;
          kbuf[i] = kn;
        }
      }

      for (int m = 0; m < p_xl->n; m++){
        const int istart = p_xl->istart[m];
        const int nlev = p_xl->nlev[m];
        for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
          const double cat = do_ocg ? (swap_xxt ? cats[abs(ordered_indices[i] - 1)] : s_cat) : 0.0;
          const double c1 = swap_xxt ? x : xt[i];
          const double c3 = do_ocg ? cat : (swap_xxt ? xt[i] : x);

          p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*
            (fast_p_kernel
              ? np_ordered_eval_kernel(P_KERNEL, c1, c3, lambda, max_cxy, lpow, cats, ncat, cl, ch)
              : k[P_KERNEL](c1, c3, lambda, cl, ch));
        }
      }


      for(l = 0, r = 0; l < P_NIDX; l++, r+=bin_do_xw){
        if(l == P_IDX) continue;
        for (int m = 0; m < xl->n; m++){
          const int istart = xl->istart[m];
          const int nlev = xl->nlev[m];
          for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
            p_result[l*num_xt + i] = pxw[r*num_xt + j]*kbuf[i];
          }
        }
      }
    }

  if(own_kbuf)
    free(kbuf);
  if(lpow != NULL)
    free(lpow);
}

void np_okernelv(const int KERNEL, 
                 const double * const xt, const int num_xt, 
                 const int do_xw,
                 const double x, const double lambda,
                 const double * cats, const int ncat,
                 double * const result,
                 const XL * const xl,
                 const int swap_xxt){
  
  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  const double cl = (cats != NULL)? cats[0] : 0.0;
  const double ch = (cats != NULL)? cats[ncat - 1] : 0.0;
  const int max_cxy = (int)fabs(ch-cl);
  const int fast_kernel = (KERNEL >= 0 && KERNEL <= 3 && cats != NULL);

  if(fast_kernel && max_cxy >= 0){
    double *lpow = (double *)malloc((size_t)(max_cxy+1)*sizeof(double));
    if(lpow == NULL) error("memory allocation failed");
    lpow[0] = 1.0;
    for(int c = 1; c <= max_cxy; c++)
      lpow[c] = lpow[c-1]*lambda;

    if(!swap_xxt){
      if(xl == NULL){
        if(!bin_do_xw){
          for(i = 0; i < num_xt; i++)
            result[i] = np_ordered_eval_kernel(KERNEL, xt[i], x, lambda, max_cxy, lpow, cats, ncat, cl, ch);
        } else {
          for(i = 0; i < num_xt; i++){
            if(xw[i] == 0.0) continue;
            result[i] = xw[i]*np_ordered_eval_kernel(KERNEL, xt[i], x, lambda, max_cxy, lpow, cats, ncat, cl, ch);
          }
        }
      } else {
        for(int m = 0; m < xl->n; m++){
          const int istart = xl->istart[m];
          const int nlev = xl->nlev[m];
          if(!bin_do_xw){
            for(i = istart; i < istart+nlev; i++)
              result[i] = np_ordered_eval_kernel(KERNEL, xt[i], x, lambda, max_cxy, lpow, cats, ncat, cl, ch);
          } else {
            for(i = istart; i < istart+nlev; i++){
              if(xw[i] == 0.0) continue;
              result[i] = xw[i]*np_ordered_eval_kernel(KERNEL, xt[i], x, lambda, max_cxy, lpow, cats, ncat, cl, ch);
            }
          }
        }
      }
    } else {
      if(xl == NULL){
        if(!bin_do_xw){
          for(i = 0; i < num_xt; i++)
            result[i] = np_ordered_eval_kernel(KERNEL, x, xt[i], lambda, max_cxy, lpow, cats, ncat, cl, ch);
        } else {
          for(i = 0; i < num_xt; i++){
            if(xw[i] == 0.0) continue;
            result[i] = xw[i]*np_ordered_eval_kernel(KERNEL, x, xt[i], lambda, max_cxy, lpow, cats, ncat, cl, ch);
          }
        }
      } else {
        for(int m = 0; m < xl->n; m++){
          const int istart = xl->istart[m];
          const int nlev = xl->nlev[m];
          if(!bin_do_xw){
            for(i = istart; i < istart+nlev; i++)
              result[i] = np_ordered_eval_kernel(KERNEL, x, xt[i], lambda, max_cxy, lpow, cats, ncat, cl, ch);
          } else {
            for(i = istart; i < istart+nlev; i++){
              if(xw[i] == 0.0) continue;
              result[i] = xw[i]*np_ordered_eval_kernel(KERNEL, x, xt[i], lambda, max_cxy, lpow, cats, ncat, cl, ch);
            }
          }
        }
      }
    }

    free(lpow);
    return;
  }


#define NP_OKERNELV_APPLY(fn)                                                      \
  do {                                                                             \
    if(!swap_xxt){                                                                 \
      if(xl == NULL){                                                              \
        if(!bin_do_xw){                                                            \
          for(i = 0; i < num_xt; i++)                                              \
            result[i] = fn(xt[i], x, lambda, cl, ch);                              \
        } else {                                                                   \
          for(i = 0; i < num_xt; i++){                                             \
            if(xw[i] == 0.0) continue;                                             \
            result[i] = xw[i]*fn(xt[i], x, lambda, cl, ch);                        \
          }                                                                        \
        }                                                                          \
      } else {                                                                     \
        for(int m = 0; m < xl->n; m++){                                            \
          const int istart = xl->istart[m];                                        \
          const int nlev = xl->nlev[m];                                            \
          if(!bin_do_xw){                                                          \
            for(i = istart; i < istart+nlev; i++)                                  \
              result[i] = fn(xt[i], x, lambda, cl, ch);                            \
          } else {                                                                 \
            for(i = istart; i < istart+nlev; i++){                                 \
              if(xw[i] == 0.0) continue;                                           \
              result[i] = xw[i]*fn(xt[i], x, lambda, cl, ch);                      \
            }                                                                      \
          }                                                                        \
        }                                                                          \
      }                                                                            \
    } else {                                                                       \
      if(xl == NULL){                                                              \
        if(!bin_do_xw){                                                            \
          for(i = 0; i < num_xt; i++)                                              \
            result[i] = fn(x, xt[i], lambda, cl, ch);                              \
        } else {                                                                   \
          for(i = 0; i < num_xt; i++){                                             \
            if(xw[i] == 0.0) continue;                                             \
            result[i] = xw[i]*fn(x, xt[i], lambda, cl, ch);                        \
          }                                                                        \
        }                                                                          \
      } else {                                                                     \
        for(int m = 0; m < xl->n; m++){                                            \
          const int istart = xl->istart[m];                                        \
          const int nlev = xl->nlev[m];                                            \
          if(!bin_do_xw){                                                          \
            for(i = istart; i < istart+nlev; i++)                                  \
              result[i] = fn(x, xt[i], lambda, cl, ch);                            \
          } else {                                                                 \
            for(i = istart; i < istart+nlev; i++){                                 \
              if(xw[i] == 0.0) continue;                                           \
              result[i] = xw[i]*fn(x, xt[i], lambda, cl, ch);                      \
            }                                                                      \
          }                                                                        \
        }                                                                          \
      }                                                                            \
    }                                                                              \
  } while(0)

  switch(KERNEL){
    case 0: NP_OKERNELV_APPLY(np_owang_van_ryzin); break;
    case 1: NP_OKERNELV_APPLY(np_oli_racine); break;
    case 2: NP_OKERNELV_APPLY(np_onli_racine); break;
    case 3: NP_OKERNELV_APPLY(np_oracine_li_yan); break;
    case 4: NP_OKERNELV_APPLY(np_econvol_owang_van_ryzin); break;
    case 5: NP_OKERNELV_APPLY(np_onull); break;
    case 6: NP_OKERNELV_APPLY(np_econvol_onli_racine); break;
    case 7: NP_OKERNELV_APPLY(np_econvol_oracine_li_yan); break;
    case 8: NP_OKERNELV_APPLY(np_score_owang_van_ryzin); break;
    case 9: NP_OKERNELV_APPLY(np_score_oli_racine); break;
    case 10: NP_OKERNELV_APPLY(np_score_onli_racine); break;
    case 11: NP_OKERNELV_APPLY(np_score_oracine_li_yan); break;
    case 12: NP_OKERNELV_APPLY(np_cdf_owang_van_ryzin); break;
    case 13: NP_OKERNELV_APPLY(np_cdf_oli_racine); break;
    case 14: NP_OKERNELV_APPLY(np_cdf_onli_racine); break;
    case 15: NP_OKERNELV_APPLY(np_cdf_oracine_li_yan); break;
    default: {
      double (* const k[])(double, double, double, double, double) = {
        np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
      };
      const int kernel = (KERNEL >= 0 && KERNEL < (int)(sizeof(k)/sizeof(k[0]))) ? KERNEL : 0;

      if(!swap_xxt){
        if(xl == NULL){
          if(!bin_do_xw){
            for(i = 0; i < num_xt; i++)
              result[i] = k[kernel](xt[i], x, lambda, cl, ch);
          } else {
            for(i = 0; i < num_xt; i++){
              if(xw[i] == 0.0) continue;
              result[i] = xw[i]*k[kernel](xt[i], x, lambda, cl, ch);
            }
          }
        } else {
          for(int m = 0; m < xl->n; m++){
            const int istart = xl->istart[m];
            const int nlev = xl->nlev[m];
            if(!bin_do_xw){
              for(i = istart; i < istart+nlev; i++)
                result[i] = k[kernel](xt[i], x, lambda, cl, ch);
            } else {
              for(i = istart; i < istart+nlev; i++){
                if(xw[i] == 0.0) continue;
                result[i] = xw[i]*k[kernel](xt[i], x, lambda, cl, ch);
              }
            }
          }
        }
      } else {
        if(xl == NULL){
          if(!bin_do_xw){
            for(i = 0; i < num_xt; i++)
              result[i] = k[kernel](x, xt[i], lambda, cl, ch);
          } else {
            for(i = 0; i < num_xt; i++){
              if(xw[i] == 0.0) continue;
              result[i] = xw[i]*k[kernel](x, xt[i], lambda, cl, ch);
            }
          }
        } else {
          for(int m = 0; m < xl->n; m++){
            const int istart = xl->istart[m];
            const int nlev = xl->nlev[m];
            if(!bin_do_xw){
              for(i = istart; i < istart+nlev; i++)
                result[i] = k[kernel](x, xt[i], lambda, cl, ch);
            } else {
              for(i = istart; i < istart+nlev; i++){
                if(xw[i] == 0.0) continue;
                result[i] = xw[i]*k[kernel](x, xt[i], lambda, cl, ch);
              }
            }
          }
        }
      }
    }
  }

#undef NP_OKERNELV_APPLY
}

// W = A
// Y = B
// outer product = AB' (assuming column vector convention)
void np_outer_weighted_sum(double * const * const mat_A, double * const sgn_A, const int ncol_A, 
                           double * const * const mat_B, const int ncol_B,
                           double * const weights, const int num_weights,
                           const int do_leave_one_out, const int which_k,
                           const int kpow,
                           const int parallel_sum, const int which_l,
                           const int symmetric,
                           const int gather_scatter,
                           const int bandwidth_divide, const double dband,
                           double * const result,
                           const XL * const xl){

  int i,j,k, l = parallel_sum?which_l:0;
  const int kstride = (parallel_sum ? (MAX(ncol_A, 1)*MAX(ncol_B, 1)) : 0);
  const int max_A = MAX(ncol_A, 1);
  const int max_B = MAX(ncol_B, 1);
  const int max_AB = max_A*max_B;
  const int have_A = (ncol_A == 0 ? 0 : 1);
  const int have_B =  (ncol_B == 0 ? 0 : 1);
  const int have_sgn = (sgn_A != NULL);
  const int linc = !parallel_sum;

  double unit_weight = 1.0;
  double * const punit_weight = &unit_weight;

  double * const * const pmat_A = (ncol_A == 0 ? 0 : 1)?
    mat_A:&punit_weight;

  double * const p_sgn = (have_sgn == 0 ? 0 : 1)?
    sgn_A:&unit_weight;

  double * const * const pmat_B = (ncol_B == 0 ? 0 : 1)?
    mat_B:&punit_weight;

  const double db = (bandwidth_divide ? dband : unit_weight);
  double temp = DBL_MAX;
  double *wbuf = NULL;
  const int use_wpow = (kpow != 1);
  const int scalar_sum_fast =
    (!have_sgn) &&
    (!symmetric) &&
    (max_A == 1) &&
    (max_B == 1) &&
    (kpow == 1);

  if(use_wpow){
    wbuf = (double *)malloc((size_t)num_weights*sizeof(double));
    if(wbuf == NULL) error("memory allocation failed");
    for(k = 0; k < num_weights; k++){
      wbuf[k] = (weights[k] == 0.0) ? 0.0 : ipow(weights[k]/db, kpow);
    }
  }

  if (do_leave_one_out) {
    temp = weights[which_k];
    weights[which_k] = 0.0;
  }

  if(scalar_sum_fast){
    if(xl == NULL){
      if(!parallel_sum){
        double acc = 0.0;

        for(k = 0; k < num_weights; k++){
          if(weights[k] == 0.0) continue;
          acc += pmat_A[0][k*have_A]*pmat_B[0][k*have_B]*weights[k]/db;
        }

        result[0] += acc;
      } else {
        l = which_l;
        for(k = 0; k < num_weights; k++){
          if(weights[k] == 0.0) continue;
          result[k] += pmat_A[0][l*have_A]*pmat_B[0][l*have_B]*weights[k]/db;
        }
      }
    } else {
      if(!parallel_sum){
        double acc = 0.0;

        for(int m = 0; m < xl->n; m++){
          const int istart = xl->istart[m];
          const int nlev = xl->nlev[m];

          for(k = istart; k < istart+nlev; k++){
            if(weights[k] == 0.0) continue;
            acc += pmat_A[0][k*have_A]*pmat_B[0][k*have_B]*weights[k]/db;
          }
        }

        result[0] += acc;
      } else {
        l = which_l;
        for(int m = 0; m < xl->n; m++){
          const int istart = xl->istart[m];
          const int nlev = xl->nlev[m];

          for(k = istart; k < istart+nlev; k++){
            if(weights[k] == 0.0) continue;
            result[k] += pmat_A[0][l*have_A]*pmat_B[0][l*have_B]*weights[k]/db;
          }
        }
      }
    }

    if(do_leave_one_out)
      weights[which_k] = temp;

    return;
  }
  
  if(xl == NULL){
    if(!gather_scatter){
      if(!symmetric){
        if(kpow == 1){
          for (k = 0; k < num_weights; k++, l+=linc){
            if(weights[k] == 0.0) continue;
            for (j = 0; j < max_A; j++)
              for (i = 0; i < max_B; i++)
                result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*weights[k]/db;
          }
        } else { // kpow != 1
          for (k = 0; k < num_weights; k++, l+=linc){
            if(weights[k] == 0.0) continue;
            for (j = 0; j < max_A; j++)
              for (i = 0; i < max_B; i++)
                result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*wbuf[k];
          }
        }
      } else { // symmetric
        if(kpow == 1){
          for (k = 0; k < num_weights; k++, l+=linc){
            if(weights[k] == 0.0) continue;
            for (j = 0; j < max_A; j++){
              for (i = 0; i <= j; i++){
                result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*weights[k]/db;
              }
            }
          }
        } else { // kpow != 1
          for (k = 0; k < num_weights; k++, l+=linc){
            if(weights[k] == 0.0) continue;
            for (j = 0; j < max_A; j++){
              for (i = 0; i <= j; i++){
                result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*wbuf[k];
              }
            }
          }
        }


        for (j = 0; j < max_A; j++){
          for (i = (max_A-1); i > j; i--){
            result[j*max_B+i] = result[i*max_B+j];
          }
        }
      }
    } else { // gather_scatter
      if(!symmetric){
        if(kpow == 1){
          for (k = 0; k < num_weights; k++, l+=linc){
            if(weights[k] == 0.0) continue;
            for (j = 0; j < max_A; j++)
              for (i = 0; i < max_B; i++)
                result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*weights[k]/db;
          }
        } else { // kpow != 1
          for (k = 0; k < num_weights; k++, l+=linc){
            if(weights[k] == 0.0) continue;
            for (j = 0; j < max_A; j++)
              for (i = 0; i < max_B; i++)
                result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*wbuf[k];
          }
        }
      } else { // symmetric
        if(kpow == 1){
          for (k = 0; k < num_weights; k++){
            if(weights[k] != 0.0){
              for (j = 0; j < max_A; j++){
                for (i = 0; i <= j; i++){
                  const double tp = pmat_A[j][k*have_A]*pmat_B[i][k*have_B]*weights[k]/db;
                  const double sgnp = tp*p_sgn[j*have_sgn]*p_sgn[i*have_sgn];
                  result[k*kstride+j*max_B+i] += tp;
                  result[(k+1)*max_AB+j*max_B+i] += sgnp;
                  //fprintf(stderr,"\nmax ab %d\tsgnp %e\tres %e",max_AB,sgnp,result[(k+1)*max_AB+j*max_B+i]);
                }
              }
            }
            // insert reference weight in top right corner
            result[(k+1)*max_AB+max_B-1] = weights[k]/db;
          }
        } else { // kpow != 1
          for (k = 0; k < num_weights; k++){
            if(weights[k] == 0.0) continue;
            for (j = 0; j < max_A; j++){
              for (i = 0; i <= j; i++){
                result[k*kstride+j*max_B+i] += pmat_A[j][k*have_A]*pmat_B[i][k*have_B]*wbuf[k];
              }
            }
          }
        }

        for (j = 0; j < max_A; j++){
          for (i = (max_A-1); i > j; i--){
            result[j*max_B+i] = result[i*max_B+j];
          }
        }

      }
    }
  } else { // TREE SUPPORT
    if(!gather_scatter){
      if(!symmetric){
        if(kpow == 1){
          for (int m = 0; m < xl->n; m++){
            const int istart = xl->istart[m];
            const int nlev = xl->nlev[m];
            l = linc ? istart : l;
            for (k = istart; k < istart+nlev; k++, l+=linc){
              if(weights[k] == 0.0) continue;
              for (j = 0; j < max_A; j++)
                for (i = 0; i < max_B; i++)
                  result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*weights[k]/db;
            }
          }
        } else { // kpow != 1
          for (int m = 0; m < xl->n; m++){
            const int istart = xl->istart[m];
            const int nlev = xl->nlev[m];
            l = linc ? istart : l;
            for (k = istart; k < istart+nlev; k++, l+=linc){
              if(weights[k] == 0.0) continue;
              for (j = 0; j < max_A; j++)
                for (i = 0; i < max_B; i++)
                  result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*wbuf[k];
            }
          }
        }
      } else { // symmetric
        if(kpow == 1){
          for (int m = 0; m < xl->n; m++){
            const int istart = xl->istart[m];
            const int nlev = xl->nlev[m];
            l = linc ? istart : l;
            for (k = istart; k < istart+nlev; k++, l+=linc){
              if(weights[k] == 0.0) continue;
              for (j = 0; j < max_A; j++){
                for (i = 0; i <= j; i++){
                  result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*weights[k]/db;
                }
              }
            }
          }
        } else { // kpow != 1
          for (int m = 0; m < xl->n; m++){
            const int istart = xl->istart[m];
            const int nlev = xl->nlev[m];
            l = linc ? istart : l;
            for (k = istart; k < istart+nlev; k++, l+=linc){
              if(weights[k] == 0.0) continue;
              for (j = 0; j < max_A; j++){
                for (i = 0; i <= j; i++){
                  result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*wbuf[k];
                }
              }
            }
          }
        }


        for (j = 0; j < max_A; j++){
          for (i = (max_A-1); i > j; i--){
            result[j*max_B+i] = result[i*max_B+j];
          }
        }
      }
    } else { // gather_scatter
      if(!symmetric){
        if(kpow == 1){
          for (int m = 0; m < xl->n; m++){
            const int istart = xl->istart[m];
            const int nlev = xl->nlev[m];
            l = linc ? istart : l;
            for (k = istart; k < istart+nlev; k++, l+=linc){
              if(weights[k] == 0.0) continue;
              for (j = 0; j < max_A; j++)
                for (i = 0; i < max_B; i++)
                  result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*weights[k]/db;
            }
          }
        } else { // kpow != 1
          for (int m = 0; m < xl->n; m++){
            const int istart = xl->istart[m];
            const int nlev = xl->nlev[m];
            l = linc ? istart : l;
            for (k = istart; k < istart+nlev; k++, l+=linc){
              if(weights[k] == 0.0) continue;
              for (j = 0; j < max_A; j++)
                for (i = 0; i < max_B; i++)
                  result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*wbuf[k];
            }
          }
        }
      } else { // symmetric
        if(kpow == 1){
          for (int m = 0; m < xl->n; m++){
            const int istart = xl->istart[m];
            const int nlev = xl->nlev[m];
            l = linc ? istart : l;
            for (k = istart; k < istart+nlev; k++){
              if(weights[k] != 0.0){
                for (j = 0; j < max_A; j++){
                  for (i = 0; i <= j; i++){
                    const double tp = pmat_A[j][k*have_A]*pmat_B[i][k*have_B]*weights[k]/db;
                    const double sgnp = tp*p_sgn[j*have_sgn]*p_sgn[i*have_sgn];
                    result[k*kstride+j*max_B+i] += tp;
                    result[(k+1)*max_AB+j*max_B+i] += sgnp;
                    //fprintf(stderr,"\nmax ab %d\tsgnp %e\tres %e",max_AB,sgnp,result[(k+1)*max_AB+j*max_B+i]);
                  }
                }
              }
              // insert reference weight in top right corner
              result[(k+1)*max_AB+max_B-1] = weights[k]/db;
            }
          }
        } else { // kpow != 1
          for (int m = 0; m < xl->n; m++){
            const int istart = xl->istart[m];
            const int nlev = xl->nlev[m];
            l = linc ? istart : l;
            for (k = istart; k < istart+nlev; k++){
              if(weights[k] == 0.0) continue;
              for (j = 0; j < max_A; j++){
                for (i = 0; i <= j; i++){
                  result[k*kstride+j*max_B+i] += pmat_A[j][k*have_A]*pmat_B[i][k*have_B]*wbuf[k];
                }
              }
            }
          }
        }

        for (j = 0; j < max_A; j++){
          for (i = (max_A-1); i > j; i--){
            result[j*max_B+i] = result[i*max_B+j];
          }
        }

      }
    }
  }

  if (do_leave_one_out)
    weights[which_k] = temp;

  safe_free(wbuf);
}



//Warning: the MPI operations used by this function assume that the
//receiving buffers, i.e. weighted_sum have enough space allocated such that a write
//consisting of nproc*stride will not segfault. It is up to the caller
//to ensure that this is the case.

// this will be fixed by using the mpi*v functions

static double * kernel_weighted_sum_pkw_extern = NULL;
static int kernel_weighted_sum_pkw_nvar_extern = 0;

static int kernel_weighted_sum_np_ctx_ex(
int * KERNEL_reg,
int * KERNEL_unordered_reg,
int * KERNEL_ordered_reg,
const int BANDWIDTH_reg,
const int num_obs_train,
const int num_obs_eval,
const int num_reg_unordered,
const int num_reg_ordered,
const int num_reg_continuous,
const int leave_one_out,
const int leave_one_out_offset,
const int kernel_pow,
const int bandwidth_divide,
const int bandwidth_divide_weights,
const int symmetric,
const int gather_scatter,
const int drop_one_train,
const int drop_which_train,
const int * const operator,
const int permutation_operator,
int do_score,
int do_ocg,
int * bpso,
const int suppress_parallel,
const int ncol_Y,
const int ncol_W,
const int int_TREE,
const int do_partial_tree,
KDT * const kdt,
NL * const inl,
int * const nld,
int * const idx,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double **matrix_Y,
double **matrix_W,
double * sgn,
double *vector_scale_factor,
int bandwidth_provided,
double ** matrix_bw_train,
double ** matrix_bw_eval,
double * lambda_pre,
int *num_categories,
double **matrix_categorical_vals,
int ** matrix_ordered_indices,
double * const weighted_sum,
double * const weighted_permutation_sum,
double * const kw,
const NP_GateOverrideCtx * const gate_override_ctx,
const int keep_kw_owner_local){
  const NP_GateOverrideCtx * const gate_ctx_raw =
    (gate_override_ctx != NULL) ? gate_override_ctx : &np_gate_override_ctx;
  const NP_GateOverrideCtx gate_ctx_empty = {0};
  const NP_GateOverrideCtx * const gate_ctx =
    np_gate_ctx_is_sane(gate_ctx_raw) ? gate_ctx_raw : &gate_ctx_empty;
  const int caller_override_active =
    ((gate_ctx != NULL) && (gate_ctx->active == NP_GATE_CTX_OVERRIDE));
  const int disable_gate_features =
    ((!np_partial_gate_features_enabled) && (!caller_override_active)) ||
    ((gate_ctx != NULL) && (gate_ctx->active == NP_GATE_CTX_DISABLE));
  assert(np_gate_ctx_is_sane(gate_ctx));
  
  /* This function takes a vector Y and returns a kernel weighted
     leave-one-out sum. By default Y should be a vector of ones
     (simply compute the kernel sum). This function will allow users
     to `roll their own' with mixed data leave-one out kernel sums. */

  /* Declarations */

  int i, ii, j, kk, k, l, mstep, js, je, num_obs_eval_alloc, sum_element_length, ip;
  int status = 0;
  int do_psum, swap_xxt;

  int * permutation_kernel = NULL;
  int doscoreocg = do_score || do_ocg;
  int do_perm = permutation_operator != OP_NOOP; 

  
  const int no_bpso = (NULL == bpso);

  int p_nvar;

  if(!np_runtime_tol_cache_ready)
    np_refresh_runtime_tolerances();

  if(no_bpso){
    bpso = (int *)malloc((num_reg_unordered + num_reg_ordered + num_reg_continuous)*sizeof(int));

    for(i = 0; i < num_reg_continuous; i++)
      bpso[i] = do_perm;

    for(i = num_reg_continuous; i < num_reg_continuous+num_reg_unordered+num_reg_ordered; i++)
      bpso[i] = doscoreocg;

    p_nvar = (do_perm ? num_reg_continuous : 0) + (doscoreocg ? num_reg_unordered + num_reg_ordered : 0);

  } else {
    // it's possible that do_perm and do_ocg could be set, but no continuous or categorical
    // variables are enabled via bpso, e.g. conditional density gradients

    for(i = 0, ii = 0, l = 0; i < num_reg_continuous; i++){
      ii |= bpso[i];
      l += bpso[i];
    }

    do_perm &= ii;

    for(i = num_reg_continuous, ii = 0; i < num_reg_continuous + num_reg_unordered + num_reg_ordered; i++){
      ii |= bpso[i];
      l += bpso[i];
    }

    do_ocg &= ii;

    doscoreocg = do_score || do_ocg;
    p_nvar = l;
  }

  int * ps_ukernel = NULL, * ps_okernel = NULL;

  int ps_ok_nli = (num_reg_ordered != 0) && (KERNEL_ordered_reg[0] != 1);

  /* Trees are currently not compatible with all operations */
  int np_ks_tree_use = (int_TREE == NP_TREE_TRUE);
  int any_convolution = 0;
  int is_adaptive = (BANDWIDTH_reg == BW_ADAP_NN);

  int lod = 0;

  const int nws = (weighted_sum == NULL);
  
  assert(!(do_score && do_ocg));
  assert(!(gather_scatter && is_adaptive));

  for(i = 0; (i < (num_reg_unordered + num_reg_ordered + num_reg_continuous)); i++){
    any_convolution |= (operator[i] == OP_CONVOLUTION);
  }

  if(any_convolution && is_adaptive) np_ks_tree_use = 0;

  np_ks_tree_use &= (num_reg_continuous != 0);

  int p_ipow = 0, * bpow = NULL;

  NL nls = {.node = NULL, .n = 0, .nalloc = 0};
  XL xl = {.istart = NULL, .nlev = NULL, .n = 0, .nalloc = 0};

  XL * pxl=  np_ks_tree_use ? &xl : NULL;
  XL * p_pxl =  NULL;

  if((np_ks_tree_use) && (p_nvar > 0)){
    p_pxl = (XL *)malloc(p_nvar*sizeof(XL));
  }

  // root node
  nls.node = (int *)malloc(sizeof(int));
  nls.node[0] = 0;
  nls.nalloc = nls.n = 1;

#ifdef MPI2
  // switch parallelisation strategies based on biggest stride
 
  int stride = MAX((int)ceil((double) num_obs_eval / (double) iNum_Processors),1);
  num_obs_eval_alloc = suppress_parallel ? num_obs_eval : (stride*iNum_Processors);

  int * igatherv = NULL, * idisplsv = NULL;

  igatherv = (int *)malloc(iNum_Processors*sizeof(int));
  if(igatherv == NULL) error("memory allocation failed");
  idisplsv = (int *)malloc(iNum_Processors*sizeof(int));
  if(idisplsv == NULL) error("memory allocation failed");
#else
  num_obs_eval_alloc = num_obs_eval;
#endif

  // we now allow one to use derivative, integral, convolution or standard kernels on a regressor by regressor basis
  int * KERNEL_reg_np = NULL, * KERNEL_unordered_reg_np = NULL, * KERNEL_ordered_reg_np = NULL;

  KERNEL_reg_np = (int *)malloc(sizeof(int)*num_reg_continuous);
  KERNEL_unordered_reg_np = (int *)malloc(sizeof(int)*num_reg_unordered);
  KERNEL_ordered_reg_np = (int *)malloc(sizeof(int)*num_reg_ordered);

  for(l = 0; l < num_reg_continuous; l++)
    KERNEL_reg_np[l] = KERNEL_reg[l] + OP_CFUN_OFFSETS[operator[l]];

  for(l = num_reg_continuous; l < (num_reg_continuous + num_reg_unordered); l++)
    KERNEL_unordered_reg_np[l - num_reg_continuous] = KERNEL_unordered_reg[l - num_reg_continuous] + OP_UFUN_OFFSETS[operator[l]];

  // todo - add (better) support for ordered integral / convolution kernels
  for(l = (num_reg_continuous+num_reg_unordered); l < (num_reg_continuous + num_reg_unordered + num_reg_ordered); l++)
    KERNEL_ordered_reg_np[l - (num_reg_continuous + num_reg_unordered)] = KERNEL_ordered_reg[l - (num_reg_continuous + num_reg_unordered)] + OP_OFUN_OFFSETS[operator[l]];

  const int num_xt = is_adaptive?num_obs_eval:num_obs_train;
  const int progress_total = is_adaptive ? num_obs_train : num_obs_eval;
  const int ws_step = is_adaptive? 0 :
                                 (MAX(ncol_Y, 1) * MAX(ncol_W, 1));

  double *lambda = NULL, **matrix_bandwidth = NULL, **matrix_alt_bandwidth = NULL, **m = NULL;
  double *tprod = NULL, dband, *ws, * p_ws, * tprod_mp = NULL, * p_dband = NULL;
  double *perm_kbuf = NULL, *kw_work = NULL;
  int use_disc_profile_cache = 0, disc_nprof = 0, disc_mark_token = 1;
  int disc_profile_from_override = 0;
  int disc_profile_from_global_cache = 0;
  int *disc_prof_id = NULL, *disc_prof_rep = NULL, *disc_prof_mark = NULL;
  int *disc_prof_list = NULL, *disc_active_idx = NULL;
  uint64_t *disc_prof_hash = NULL;
  double *disc_prof_val = NULL, *disc_ord_cl = NULL, *disc_ord_ch = NULL;
  int *disc_uno_const_ok = NULL, *disc_ord_const_ok = NULL;
  double *disc_uno_const = NULL, *disc_ord_const = NULL;
  int *cont_largeh_ok = NULL;
  double *cont_largeh_hmin = NULL, *cont_largeh_k0 = NULL;
  double cont_largeh_rel_tol = 1e-3;
  int *cont_largeh_active = NULL, *cont_largeh_active_fixed = NULL;
  int *tree_active_dims = NULL;
  int cont_largeh_all_fixed = 0, cont_largeh_fixed_ready = 0;
  int cont_largeh_any_fixed = 0;
  int tree_alllarge_bypass = 0;
  int cont_largeh_from_override = 0;
  int cont_largeh_from_global_cache = 0;
  int disc_uno_from_override = 0, disc_ord_from_override = 0;
  const int lean_reg_cont_loop =
    (!doscoreocg) &&
    (p_nvar == 0) &&
    (!int_cker_bound_extern) &&
    (!any_convolution);

  double * const * const xtc = is_adaptive?
    matrix_X_continuous_eval:matrix_X_continuous_train;
  double * const * const xtu = is_adaptive?
    matrix_X_unordered_eval:matrix_X_unordered_train;
  double * const * const xto = is_adaptive?
    matrix_X_ordered_eval:matrix_X_ordered_train;

  double * const * const xc = is_adaptive?
    matrix_X_continuous_train:matrix_X_continuous_eval;
  double * const * const xu = is_adaptive?
    matrix_X_unordered_train:matrix_X_unordered_eval;
  double * const * const xo = is_adaptive?
    matrix_X_ordered_train:matrix_X_ordered_eval;

#ifdef MPI2
  if((kw != NULL) && (!suppress_parallel)){
    kw_work = (double *)calloc((size_t)num_obs_eval_alloc*(size_t)num_xt, sizeof(double));
    if(kw_work == NULL){
      status = KWSNP_ERR_BADINVOC;
      goto cleanup;
    }
  } else {
    kw_work = kw;
  }
#else
  kw_work = kw;
#endif

  if (num_obs_eval == 0) {
    status = KWSNP_ERR_NOEVAL;
    goto cleanup;
  }

  do_psum = BANDWIDTH_reg == BW_ADAP_NN;
  swap_xxt = BANDWIDTH_reg == BW_ADAP_NN;
  /* Allocate memory for objects */

  mstep = (BANDWIDTH_reg==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_train:1);


  if(bandwidth_provided){
    if(BANDWIDTH_reg == BW_GEN_NN){
      matrix_bandwidth = matrix_bw_eval;
      if(any_convolution)
        matrix_alt_bandwidth = matrix_bw_train;
    }
    else if (is_adaptive){
      if (any_convolution){
        matrix_alt_bandwidth = matrix_bw_eval;
      }
      matrix_bandwidth = matrix_bw_train;
    } else {
      matrix_bandwidth = matrix_bw_train;
      if(any_convolution)
        matrix_alt_bandwidth = matrix_bw_train;
    }
    lambda = lambda_pre;
  } else {
    matrix_bandwidth = alloc_tmatd(mstep, num_reg_continuous);  
    lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);
    if(any_convolution && (BANDWIDTH_reg == BW_FIXED))
      matrix_alt_bandwidth = matrix_bandwidth;
  } 


  tprod = alloc_vecd((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_eval:num_obs_train);

  sum_element_length = MAX(ncol_Y, 1) * 
    MAX(ncol_W, 1);

  /* assert(!(BANDWIDTH_reg == BW_ADAP_NN)); */
  /* Conduct the estimation */

  /* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

  bpow = (int *) malloc(num_reg_continuous*sizeof(int));
  if(bpow == NULL){
    status = KWSNP_ERR_BADINVOC;
    goto cleanup;
  }

  for(i = 0; i < num_reg_continuous; i++){
    bpow[i] = (bandwidth_divide ? 1 : 0) + ((operator[i] == OP_DERIVATIVE) ? 1 : ((operator[i] == OP_INTEGRAL) ? -1 : 0));
  }

  if(p_nvar > 0){
    if(do_perm){
      permutation_kernel = (int *) malloc(num_reg_continuous*sizeof(int));
      for(i = 0; i < num_reg_continuous; i++){
        permutation_kernel[i] = KERNEL_reg[i] + OP_CFUN_OFFSETS[permutation_operator];
      }
    }
    if(do_score){
      ps_ukernel = (int *) malloc(num_reg_unordered*sizeof(int));
      ps_okernel = (int *) malloc(num_reg_ordered*sizeof(int));

      for(i = 0; i < num_reg_unordered; i++){
        ps_ukernel[i] = KERNEL_unordered_reg[i] + OP_UFUN_OFFSETS[permutation_operator];
      }

      for(i = 0; i < num_reg_ordered; i++){
        ps_okernel[i] = KERNEL_ordered_reg[i] + OP_OFUN_OFFSETS[permutation_operator];
      }
    } else if(do_ocg) {
      ps_ukernel = KERNEL_unordered_reg;
      ps_okernel = KERNEL_ordered_reg;
    }

    tprod_mp = (double *)malloc(((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_eval:num_obs_train)*p_nvar*sizeof(double));
    if(tprod_mp == NULL){
      status = KWSNP_ERR_BADINVOC;
      goto cleanup;
    }

    p_dband = (double *)malloc(p_nvar*sizeof(double));
    if(p_dband == NULL){
      status = KWSNP_ERR_BADINVOC;
      goto cleanup;
    }

    p_ipow = (bandwidth_divide ? 1 : 0) + ((permutation_operator == OP_DERIVATIVE) ? 1 : ((permutation_operator == OP_INTEGRAL) ? -1 : 0));

    perm_kbuf = (double *)malloc((size_t)num_xt*sizeof(double));
    if(perm_kbuf == NULL){
      status = KWSNP_ERR_BADINVOC;
      goto cleanup;
    }

  }

  if(!bandwidth_provided){
    if(kernel_bandwidth_mean((num_reg_continuous != 0) ? KERNEL_reg[0]: 0,
                             BANDWIDTH_reg,
                             num_obs_train,
                             num_obs_eval,
                             0,
                             0,
                             0,
                             num_reg_continuous,
                             num_reg_unordered,
                             num_reg_ordered,
                             suppress_parallel, // suppress_parallel if requested
                             vector_scale_factor,
                             matrix_X_continuous_train,				 /* Not used */
                             matrix_X_continuous_eval,				 /* Not used */
                             matrix_X_continuous_train,
                             matrix_X_continuous_eval,
                             matrix_bandwidth,						 /* Not used */
                             matrix_bandwidth,
                             lambda)==1){

      status = KWSNP_ERR_BADBW;
      goto cleanup;
    }
  }

  if(!bandwidth_provided && (is_adaptive && any_convolution)){ // need additional bandwidths
    matrix_alt_bandwidth = alloc_tmatd(num_obs_eval, num_reg_continuous);  

    // Adaptive convolution requires an auxiliary BW_GEN_NN bandwidth matrix.
    if(kernel_bandwidth_mean((num_reg_continuous != 0) ? KERNEL_reg[0]: 0,
                             BW_GEN_NN, // this is not an error!
                             num_obs_train,
                             num_obs_eval,
                             0,
                             0,
                             0,
                             num_reg_continuous,
                             num_reg_unordered,
                             num_reg_ordered,
                             suppress_parallel, // suppress_parallel if requested
                             vector_scale_factor,
                             NULL,				 /* Not used */
                             NULL,				 /* Not used */
                             matrix_X_continuous_train,
                             matrix_X_continuous_eval,
                             NULL,						 /* Not used */
                             matrix_alt_bandwidth,
                             lambda)==1){

      status = KWSNP_ERR_BADBW;
      goto cleanup;
    }

  }

  if((!disable_gate_features) && (num_reg_continuous > 0)){
    if((gate_ctx->active == NP_GATE_CTX_OVERRIDE) &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->cont_ok != NULL &&
       gate_ctx->cont_hmin != NULL &&
       gate_ctx->cont_k0 != NULL){
      cont_largeh_ok = (int *)gate_ctx->cont_ok;
      cont_largeh_hmin = (double *)gate_ctx->cont_hmin;
      cont_largeh_k0 = (double *)gate_ctx->cont_k0;
      cont_largeh_from_override = 1;
    } else if(gate_ctx->active &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->cont_ok != NULL &&
       gate_ctx->cont_hmin != NULL &&
       gate_ctx->cont_k0 != NULL){
      cont_largeh_ok = (int *)gate_ctx->cont_ok;
      cont_largeh_hmin = (double *)gate_ctx->cont_hmin;
      cont_largeh_k0 = (double *)gate_ctx->cont_k0;
      cont_largeh_from_override = 1;
    }

    if(!cont_largeh_from_override){
      cont_largeh_rel_tol = np_cont_largeh_rel_tol();
      if(np_cont_largeh_cache_get_or_build(num_obs_train,
                                           num_obs_eval,
                                           num_reg_continuous,
                                           KERNEL_reg_np,
                                           matrix_X_continuous_train,
                                           matrix_X_continuous_eval,
                                           cont_largeh_rel_tol,
                                           &cont_largeh_ok,
                                           &cont_largeh_hmin,
                                           &cont_largeh_k0)){
        cont_largeh_from_global_cache = 1;
      } else {
        if(!np_cont_largeh_build_params(num_obs_train,
                                        num_obs_eval,
                                        num_reg_continuous,
                                        KERNEL_reg_np,
                                        matrix_X_continuous_train,
                                        matrix_X_continuous_eval,
                                        cont_largeh_rel_tol,
                                        &cont_largeh_ok,
                                        &cont_largeh_hmin,
                                        &cont_largeh_k0)){
          cont_largeh_ok = NULL;
          cont_largeh_hmin = NULL;
          cont_largeh_k0 = NULL;
        }
      }
    }
  }

  if((!disable_gate_features) && (num_reg_continuous > 0)){
    cont_largeh_active = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
    if(BANDWIDTH_reg == BW_FIXED){
      cont_largeh_active_fixed = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
      if(cont_largeh_active != NULL && cont_largeh_active_fixed != NULL){
        cont_largeh_all_fixed = 1;
        cont_largeh_any_fixed = 0;
        for(i = 0; i < num_reg_continuous; i++){
          const int is_largeh = (cont_largeh_ok != NULL) &&
            cont_largeh_ok[i] &&
            (operator[i] != OP_CONVOLUTION) &&
            isfinite(matrix_bandwidth[i][0]) &&
            (fabs(matrix_bandwidth[i][0]) >= cont_largeh_hmin[i]);

          cont_largeh_active_fixed[i] = is_largeh;
          cont_largeh_active[i] = is_largeh;
          cont_largeh_all_fixed &= is_largeh;
          cont_largeh_any_fixed |= is_largeh;
        }
        cont_largeh_fixed_ready = 1;
      }
    }
  }

  if(np_ks_tree_use &&
     cont_largeh_fixed_ready &&
     cont_largeh_all_fixed &&
     (p_nvar == 0) &&
     (!do_partial_tree)){
    /* All continuous kernels collapse to constants for every eval point.
       Bypass per-j tree list construction and use dense accumulation. */
    tree_alllarge_bypass = 1;
    pxl = NULL;
  }

  if(np_ks_tree_use && (!tree_alllarge_bypass) && (num_reg_continuous > 0)){
    tree_active_dims = (int *)malloc((size_t)num_reg_continuous*sizeof(int));
  }

  if((!disable_gate_features) && (num_reg_unordered > 0)){
    if((gate_ctx->active == NP_GATE_CTX_OVERRIDE) &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->disc_uno_ok != NULL &&
       gate_ctx->disc_uno_const != NULL){
      disc_uno_const_ok = (int *)gate_ctx->disc_uno_ok;
      disc_uno_const = (double *)gate_ctx->disc_uno_const;
      disc_uno_from_override = 1;
    } else if(gate_ctx->active &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->disc_uno_ok != NULL &&
       gate_ctx->disc_uno_const != NULL){
      disc_uno_const_ok = (int *)gate_ctx->disc_uno_ok;
      disc_uno_const = (double *)gate_ctx->disc_uno_const;
      disc_uno_from_override = 1;
    }

    if(!disc_uno_from_override){
    disc_uno_const_ok = (int *)calloc((size_t)num_reg_unordered, sizeof(int));
    disc_uno_const = (double *)malloc((size_t)num_reg_unordered*sizeof(double));

    if(disc_uno_const_ok != NULL && disc_uno_const != NULL){
      double (* const ukf[])(int, double, int) = {
        np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
        np_score_uaa, np_score_unli_racine
      };
      const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));

      for(i = 0; i < num_reg_unordered; i++){
        const int ku = KERNEL_unordered_reg_np[i];
        const int ncat = (num_categories != NULL) ? num_categories[i] : 0;
        const double lam = lambda[i];

        disc_uno_const_ok[i] = 0;
        disc_uno_const[i] = 0.0;

        if(ku < 0 || ku >= nuk) continue;
        if(!np_disc_near_upper(ku, lam, ncat)) continue;

        const double kn_same = ukf[ku](1, lam, ncat);
        const double kn_diff = ukf[ku](0, lam, ncat);
        if(np_disc_near_const_kernel(kn_same, kn_diff)){
          disc_uno_const_ok[i] = 1;
          disc_uno_const[i] = 0.5*(kn_same + kn_diff);
        }
      }
    } else {
      if(disc_uno_const_ok != NULL){ free(disc_uno_const_ok); disc_uno_const_ok = NULL; }
      if(disc_uno_const != NULL){ free(disc_uno_const); disc_uno_const = NULL; }
    }
    }
  }

  if((!disable_gate_features) && (num_reg_ordered > 0)){
    if((gate_ctx->active == NP_GATE_CTX_OVERRIDE) &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->disc_ord_ok != NULL &&
       gate_ctx->disc_ord_const != NULL){
      disc_ord_const_ok = (int *)gate_ctx->disc_ord_ok;
      disc_ord_const = (double *)gate_ctx->disc_ord_const;
      disc_ord_from_override = 1;
    } else if(gate_ctx->active &&
       np_gate_ctx_signature_matches(gate_ctx,
                                     num_reg_continuous,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     KERNEL_reg_np,
                                     KERNEL_unordered_reg_np,
                                     KERNEL_ordered_reg_np,
                                     operator) &&
       gate_ctx->disc_ord_ok != NULL &&
       gate_ctx->disc_ord_const != NULL){
      disc_ord_const_ok = (int *)gate_ctx->disc_ord_ok;
      disc_ord_const = (double *)gate_ctx->disc_ord_const;
      disc_ord_from_override = 1;
    }

    if(!disc_ord_from_override){
    disc_ord_const_ok = (int *)calloc((size_t)num_reg_ordered, sizeof(int));
    disc_ord_const = (double *)malloc((size_t)num_reg_ordered*sizeof(double));

    if(disc_ord_const_ok != NULL && disc_ord_const != NULL){
      double (* const okf[])(double, double, double, double, double) = {
        np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
      };
      const int nok = (int)(sizeof(okf)/sizeof(okf[0]));

      for(i = 0; i < num_reg_ordered; i++){
        const int oi = i + num_reg_unordered;
        const int ko = KERNEL_ordered_reg_np[i];
        const int ncat = (num_categories != NULL) ? num_categories[oi] : 0;
        const double lam = lambda[oi];

        disc_ord_const_ok[i] = 0;
        disc_ord_const[i] = 0.0;

        if(ko < 0 || ko >= nok) continue;
        if(ncat <= 0 || matrix_categorical_vals == NULL) continue;
        if(!np_disc_ordered_near_upper(ko, lam)) continue;

        const double cl = matrix_categorical_vals[oi][0];
        const double ch = matrix_categorical_vals[oi][ncat - 1];
        const double k0 = okf[ko](cl, cl, lam, cl, ch);
        const double k1 = okf[ko](cl, ch, lam, cl, ch);
        if(np_disc_near_const_kernel(k0, k1)){
          disc_ord_const_ok[i] = 1;
          disc_ord_const[i] = 0.5*(k0 + k1);
        }
      }
    } else {
      if(disc_ord_const_ok != NULL){ free(disc_ord_const_ok); disc_ord_const_ok = NULL; }
      if(disc_ord_const != NULL){ free(disc_ord_const); disc_ord_const = NULL; }
    }
    }
  }
  
  if (leave_one_out && drop_one_train) {
    REprintf("\n error, leave one out estimator and drop-one estimator can't be enabled simultaneously");
    status = KWSNP_ERR_BADINVOC;
    goto cleanup;

  }

  if ((num_obs_train < (num_obs_eval + leave_one_out_offset)) && leave_one_out){
    
    REprintf("\nnumber of training points must be >= number of evaluation points to use the 'leave one out' estimator");
    status = KWSNP_ERR_BADINVOC;
    goto cleanup;
  }

  if(!gather_scatter && !nws)
    for(i = 0; i < num_obs_eval_alloc*sum_element_length; i++){
      weighted_sum[i] = 0.0;
    }

  if(!gather_scatter){
    for(i = 0; i < num_obs_eval_alloc*sum_element_length*p_nvar; i++){
      weighted_permutation_sum[i] = 0.0;
    }
  }

  if (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN){
#ifdef MPI2
    if((!suppress_parallel) && (!gather_scatter)){
      js = stride * my_rank;
      je = MIN(num_obs_eval - 1, js + stride - 1);

      ws = nws ? NULL : weighted_sum + js*sum_element_length;

      if(weighted_permutation_sum != NULL){
        p_ws = weighted_permutation_sum +js*sum_element_length;
      } else {
        p_ws = NULL;
      }

      // hopefully this now handles the case where there are more processors than observations
      for(i = 0, ii = 0; ii*stride < num_obs_eval; ii++){
        igatherv[ii] = stride*sum_element_length;
        idisplsv[ii] = i;
        i += stride*sum_element_length;
      }
      
      if(i < num_obs_eval*sum_element_length){
        const int de1 = (num_obs_eval - (ii-1)*stride)*sum_element_length;
        igatherv[ii] = de1;
        idisplsv[ii++] = i;
        i += de1;
      }

      for(; ii < iNum_Processors; ii++){
        igatherv[ii] = 0;
        idisplsv[ii] = i;
      }


    } else {
      js = 0;
      je = num_obs_eval - 1;
      ws = weighted_sum;
      p_ws = weighted_permutation_sum;

      for(ii = 0; ii < iNum_Processors; ii++){
        igatherv[ii] = -1;
        idisplsv[ii] = -1;
      }
    }
#else
    js = 0;
    je = num_obs_eval - 1;
    ws = weighted_sum;
    p_ws = weighted_permutation_sum;
#endif
    

  } else {
#ifdef MPI2
    if((!suppress_parallel) && (!gather_scatter)){
      js = stride * my_rank;
      je = MIN(num_obs_train - 1, js + stride - 1);
      ws = weighted_sum;
      p_ws = weighted_permutation_sum;

      for(ii = 0; ii < iNum_Processors; ii++){
        igatherv[ii] = -1;
        idisplsv[ii] = -1;
      }

    } else {
      js = 0;
      je = num_obs_train - 1;
      ws = weighted_sum;
      p_ws = weighted_permutation_sum;

      for(ii = 0; ii < iNum_Processors; ii++){
        igatherv[ii] = -1;
        idisplsv[ii] = -1;
      }

    }
#else
    js = 0;
    je = num_obs_train - 1;
    ws = weighted_sum;
    p_ws = weighted_permutation_sum;
#endif
  }

  if((!disable_gate_features) && (!doscoreocg) && (p_nvar == 0) &&
     ((num_reg_unordered + num_reg_ordered) > 1) &&
     ((num_reg_ordered > 0) || ((num_reg_unordered + num_reg_ordered) > 2)) &&
     (num_xt >= 128)){
    if(gate_ctx->disc_prof_id != NULL &&
       gate_ctx->disc_prof_rep != NULL &&
       gate_ctx->disc_nprof > 0 &&
       gate_ctx->disc_prof_num_xt == num_xt){
      disc_prof_id = (int *)gate_ctx->disc_prof_id;
      disc_prof_rep = (int *)gate_ctx->disc_prof_rep;
      disc_nprof = gate_ctx->disc_nprof;
      disc_profile_from_override = 1;
      disc_profile_from_global_cache = 0;
    } else {
      disc_profile_from_override = 0;
      if(np_disc_profile_cache_get_or_build(num_xt,
                                            num_reg_unordered,
                                            num_reg_ordered,
                                            xtu,
                                            xto,
                                            &disc_prof_id,
                                            &disc_prof_rep,
                                            &disc_nprof)){
        disc_profile_from_global_cache = 1;
      } else {
        disc_prof_id = NULL;
        disc_prof_rep = NULL;
        disc_nprof = 0;
        disc_profile_from_global_cache = 0;
      }
    }

    if((disc_nprof > 0) && (4*disc_nprof <= 3*num_xt)){
      disc_prof_mark = (int *)calloc((size_t)disc_nprof, sizeof(int));
      disc_prof_list = (int *)malloc((size_t)disc_nprof*sizeof(int));
      disc_active_idx = (int *)malloc((size_t)num_xt*sizeof(int));
      disc_prof_val = (double *)malloc((size_t)disc_nprof*sizeof(double));

      if(num_reg_ordered > 0){
        disc_ord_cl = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
        disc_ord_ch = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
      }

      if(disc_prof_mark != NULL && disc_prof_list != NULL &&
         disc_active_idx != NULL && disc_prof_val != NULL &&
         ((num_reg_ordered == 0) || (disc_ord_cl != NULL && disc_ord_ch != NULL))){
        for(i = 0; i < num_reg_ordered; i++){
          const int oi = i + num_reg_unordered;
          disc_ord_cl[i] = (matrix_categorical_vals != NULL) ? matrix_categorical_vals[oi][0] : 0.0;
          disc_ord_ch[i] = (matrix_categorical_vals != NULL) ? matrix_categorical_vals[oi][num_categories[oi] - 1] : 0.0;
        }
        use_disc_profile_cache = 1;
      }
    }

    if(!use_disc_profile_cache){
      if(disc_prof_mark != NULL){ free(disc_prof_mark); disc_prof_mark = NULL; }
      if(disc_prof_list != NULL){ free(disc_prof_list); disc_prof_list = NULL; }
      if(disc_active_idx != NULL){ free(disc_active_idx); disc_active_idx = NULL; }
      if(disc_prof_val != NULL){ free(disc_prof_val); disc_prof_val = NULL; }
      if(disc_ord_cl != NULL){ free(disc_ord_cl); disc_ord_cl = NULL; }
      if(disc_ord_ch != NULL){ free(disc_ord_ch); disc_ord_ch = NULL; }
      if((disc_prof_id != NULL) && (!disc_profile_from_override) && (!disc_profile_from_global_cache)){ free(disc_prof_id); disc_prof_id = NULL; }
      if((disc_prof_rep != NULL) && (!disc_profile_from_override) && (!disc_profile_from_global_cache)){ free(disc_prof_rep); disc_prof_rep = NULL; }
      disc_nprof = 0;
    }
    if(use_disc_profile_cache && ((disc_prof_id == NULL) || (disc_prof_rep == NULL)))
      error("discrete profile cache setup failed");
  }

  const int leave_or_drop = leave_one_out || (drop_one_train && (BANDWIDTH_reg != BW_ADAP_NN));
  if(drop_one_train && (BANDWIDTH_reg != BW_ADAP_NN)) lod = drop_which_train;

  m = matrix_bandwidth;

    /* do sums */
  for(j=js; j <= je; j++, ws = (ws==NULL ? NULL : ws+ws_step), p_ws=(p_ws == NULL ? NULL : p_ws+ws_step)){
    R_CheckUserInterrupt();

    dband = 1.0;
    const int jbw = (BANDWIDTH_reg != BW_FIXED) ? j:0;
    int tprod_has_vals = 0;
    int deferred_const_active = 0;
    double deferred_const = 1.0;
    int all_cont_largeh = (num_reg_continuous > 0);
    int any_cont_largeh = 0;
    int tree_active_n = num_reg_continuous;
    const int tree_use_active_dims =
      np_ks_tree_use &&
      (!tree_alllarge_bypass) &&
      (p_nvar == 0) &&
      (!do_partial_tree) &&
      (cont_largeh_active != NULL) &&
      (tree_active_dims != NULL);

    if(cont_largeh_fixed_ready){
      all_cont_largeh = cont_largeh_all_fixed;
      any_cont_largeh = cont_largeh_any_fixed;
    } else if(cont_largeh_active != NULL){
      for(i = 0; i < num_reg_continuous; i++){
        const int is_largeh = (cont_largeh_ok != NULL) &&
          cont_largeh_ok[i] &&
          (operator[i] != OP_CONVOLUTION) &&
          isfinite(m[i][jbw]) &&
          (fabs(m[i][jbw]) >= cont_largeh_hmin[i]);

        cont_largeh_active[i] = is_largeh;
        all_cont_largeh &= is_largeh;
        any_cont_largeh |= is_largeh;
      }
    } else {
      all_cont_largeh = 0;
      any_cont_largeh = 0;
    }

    if(tree_use_active_dims){
      tree_active_n = 0;
      for(i = 0; i < num_reg_continuous; i++){
        if(!cont_largeh_active[i]){
          tree_active_dims[tree_active_n++] = i;
        }
      }
    }

    for (ii = 0; ii < p_nvar; ii++)
      p_dband[ii] = 1.0;

    // if we are consistently dropping one obs from training data, and we are doing a parallel sum, then that means
    // we need to skip here

    if(leave_one_out) lod = j + leave_one_out_offset;

    // do a hail mary, then generate the interaction list
    // anything but a fixed bandwidth is not yet supported
    // that includes convolutions 

    if(np_ks_tree_use && (!tree_alllarge_bypass)){

      if(kw != NULL)
        for(i = 0; i < num_xt; i++)
          tprod[i] = 0.0;

      double bb[kdt->ndim*2];

      // reset the interaction node list
      xl.n = 0;
      if(all_cont_largeh || (tree_use_active_dims && (tree_active_n == 0))){
        /* Global large-h: all continuous kernels are effectively K(0), so support is full. */
        merge_end_xl(pxl, &kdt->kdn[0]);
      } else {
        if(!do_partial_tree){
          if(tree_use_active_dims && (tree_active_n < num_reg_continuous)){
            for(kk = 0; kk < tree_active_n; kk++){
              const int id = tree_active_dims[kk];
              const double sf = m[id][jbw];
              if(!is_adaptive){
                bb[2*id] = -cksup[KERNEL_reg_np[id]][1];
                bb[2*id+1] = -cksup[KERNEL_reg_np[id]][0];
              } else {
                bb[2*id] = cksup[KERNEL_reg_np[id]][0];
                bb[2*id+1] = cksup[KERNEL_reg_np[id]][1];
              }
              bb[2*id] = (fabs(bb[2*id]) == DBL_MAX) ? bb[2*id] : (xc[id][j] + bb[2*id]*sf);
              bb[2*id+1] = (fabs(bb[2*id+1]) == DBL_MAX) ? bb[2*id+1] : (xc[id][j] + bb[2*id+1]*sf);
            }

            boxSearchNLPartial(kdt, &nls, bb, NULL, pxl, tree_active_dims, tree_active_n);
          } else {
            for(i = 0; i < num_reg_continuous; i++){
              const double sf = m[i][jbw];
              if(!is_adaptive){
                bb[2*i] = -cksup[KERNEL_reg_np[i]][1];
                bb[2*i+1] = -cksup[KERNEL_reg_np[i]][0];
              }else{
                bb[2*i] = cksup[KERNEL_reg_np[i]][0];
                bb[2*i+1] = cksup[KERNEL_reg_np[i]][1];
              }
              bb[2*i] = (fabs(bb[2*i]) == DBL_MAX) ? bb[2*i] : (xc[i][j] + bb[2*i]*sf);
              bb[2*i+1] = (fabs(bb[2*i+1]) == DBL_MAX) ? bb[2*i+1] : (xc[i][j] + bb[2*i+1]*sf);
            }

            boxSearchNL(kdt, &nls, bb, NULL, pxl);
          }
        } else {
          for(i = 0; i < num_reg_continuous; i++){
            const double sf = m[i][jbw];
            if(!is_adaptive){
              bb[2*nld[i]] = -cksup[KERNEL_reg_np[i]][1];
              bb[2*nld[i]+1] = -cksup[KERNEL_reg_np[i]][0];
            }else{
              bb[2*nld[i]] = cksup[KERNEL_reg_np[i]][0];
              bb[2*nld[i]+1] = cksup[KERNEL_reg_np[i]][1];
            }

            bb[2*nld[i]] = (fabs(bb[2*nld[i]]) == DBL_MAX) ? bb[2*nld[i]] : (xc[i][j] + bb[2*nld[i]]*sf);
            bb[2*nld[i]+1] = (fabs(bb[2*nld[i]+1]) == DBL_MAX) ? bb[2*nld[i]+1] : (xc[i][j] + bb[2*nld[i]+1]*sf);
          }
          if(idx == NULL)
            boxSearchNLPartial(kdt, inl, bb, NULL, pxl, nld, num_reg_continuous);
          else {
            boxSearchNLPartialIdx(kdt, inl, bb, NULL, pxl, nld, num_reg_continuous, idx);
          }
        }
      }

      if(do_perm && (permutation_operator == OP_INTEGRAL)){
        for(ii = 0, k = 0; ii < num_reg_continuous; ii++){
          if(bpso[ii]){
            // reset the interaction node list
            p_pxl[k].n = 0;

            if(all_cont_largeh){
              p_pxl[k] = pxl[0];
            } else if(!do_partial_tree){
              for(i = 0; i < num_reg_continuous; i++){
                const int knp = (i == ii) ? permutation_kernel[i] : KERNEL_reg_np[i];
                if(!is_adaptive){
                  bb[2*i] = -cksup[knp][1];
                  bb[2*i+1] = -cksup[knp][0];
                }else{
                  bb[2*i] = cksup[knp][0];
                  bb[2*i+1] = cksup[knp][1];
                }
                const double sf = m[i][jbw];
                bb[2*i] = (fabs(bb[2*i]) == DBL_MAX) ? bb[2*i] : (xc[i][j] + bb[2*i]*sf);
                bb[2*i+1] = (fabs(bb[2*i+1]) == DBL_MAX) ? bb[2*i+1] : (xc[i][j] + bb[2*i+1]*sf);
              }

              boxSearchNL(kdt, &nls, bb, NULL, p_pxl + k);
            } else {
              for(i = 0; i < num_reg_continuous; i++){
                const int knp = (i == ii) ? permutation_kernel[i] : KERNEL_reg_np[i];
                if(!is_adaptive){
                  bb[2*nld[i]] = -cksup[knp][1];
                  bb[2*nld[i]+1] = -cksup[knp][0];
                }else{
                  bb[2*nld[i]] = cksup[knp][0];
                  bb[2*nld[i]+1] = cksup[knp][1];
                }
                const double sf = m[i][jbw];
                bb[2*nld[i]] = (fabs(bb[2*nld[i]]) == DBL_MAX) ? bb[2*nld[i]] : (xc[i][j] + bb[2*nld[i]]*sf);
                bb[2*nld[i]+1] = (fabs(bb[2*nld[i]+1]) == DBL_MAX) ? bb[2*nld[i]+1] : (xc[i][j] + bb[2*nld[i]+1]*sf);
              }
              if(idx == NULL)
                boxSearchNLPartial(kdt, inl, bb, NULL, p_pxl + k, nld, num_reg_continuous);
              else{
                boxSearchNLPartialIdx(kdt, inl, bb, NULL, p_pxl + k, nld, num_reg_continuous, idx);
              }
            }

            k++;
          }
        }
        
        for(ii = k; ii < p_nvar; ii++){
          p_pxl[ii] = pxl[0];
        }
                  
      } else {
        for(ii = 0; ii < p_nvar; ii++){
          p_pxl[ii] = pxl[0];
        }
      } 

      /* No continuous support from the tree: all kernel products are zero. */
      if(xl.n == 0){
        int any_p_support = 0;
        for(ii = 0; ii < p_nvar; ii++){
          any_p_support |= (p_pxl[ii].n > 0);
        }

        if(!any_p_support){
          if(kw != NULL){
            for(i = 0; i < num_xt; i++)
              kw[j*num_xt + i] = 0.0;
          }
          continue;
        }
      }

    }
    /* continuous first */

    /* for the first iteration, no weights */
    /* for the rest, the accumulated products are the weights */
    if(lean_reg_cont_loop){
      for(i = 0, l = 0, ip = 0, k = 0; i < num_reg_continuous; i++, l++, ip += do_perm){
        const int use_largeh = any_cont_largeh && (cont_largeh_active != NULL) ? cont_largeh_active[i] : 0;

        if(use_largeh){
          deferred_const *= cont_largeh_k0[i];
          deferred_const_active = 1;
        } else {
          np_ckernelv(KERNEL_reg_np[i], xtc[i], num_xt, tprod_has_vals, xc[i][j], m[i][jbw], tprod, pxl, swap_xxt, 1, 1.0);
          tprod_has_vals = 1;
        }

        dband *= ipow(m[i][jbw], bpow[i]);
        k += bpso[l];
      }
    } else {
      for(i = 0, l = 0, ip = 0, k = 0; i < num_reg_continuous; i++, l++, ip += do_perm){
        const int kbase_i = KERNEL_reg_np[i] % 10;
        const int p_kbase_i = (do_perm ? permutation_kernel[i] : KERNEL_reg_np[i]) % 10;
        const int use_bounds_i = int_cker_bound_extern &&
          vector_ckerlb_extern != NULL &&
          vector_ckerub_extern != NULL &&
          (kbase_i >= 0 && kbase_i <= 9) &&
          ((isfinite(vector_ckerlb_extern[i]) && (fabs(vector_ckerlb_extern[i]) < 0.5*DBL_MAX)) ||
           (isfinite(vector_ckerub_extern[i]) && (fabs(vector_ckerub_extern[i]) < 0.5*DBL_MAX)));
        const int use_p_bounds_i = int_cker_bound_extern &&
          vector_ckerlb_extern != NULL &&
          vector_ckerub_extern != NULL &&
          (p_kbase_i >= 0 && p_kbase_i <= 9) &&
          ((isfinite(vector_ckerlb_extern[i]) && (fabs(vector_ckerlb_extern[i]) < 0.5*DBL_MAX)) ||
           (isfinite(vector_ckerub_extern[i]) && (fabs(vector_ckerub_extern[i]) < 0.5*DBL_MAX)));
        const double invnorm = use_bounds_i ?
          np_cker_invnorm(KERNEL_reg_np[i], xc[i][j], m[i][jbw],
                          vector_ckerlb_extern[i], vector_ckerub_extern[i]) : 1.0;
        const double p_invnorm = use_p_bounds_i ?
          np_cker_invnorm((do_perm ? permutation_kernel[i] : KERNEL_reg_np[i]),
                          xc[i][j], m[i][jbw],
                          vector_ckerlb_extern[i], vector_ckerub_extern[i]) : 1.0;

        if(operator[l] != OP_CONVOLUTION){
          if(p_nvar == 0){
            const int use_largeh = any_cont_largeh && (cont_largeh_active != NULL) ? cont_largeh_active[i] : 0;

            if(use_largeh){
              deferred_const *= cont_largeh_k0[i]*invnorm;
              deferred_const_active = 1;
            } else {
              np_ckernelv(KERNEL_reg_np[i], xtc[i], num_xt, tprod_has_vals, xc[i][j], m[i][jbw], tprod, pxl, swap_xxt, 1, invnorm);
              tprod_has_vals = 1;
            }
          } else {
            const int use_largeh = any_cont_largeh && (cont_largeh_active != NULL) ? cont_largeh_active[i] : 0;

            np_p_ckernelv(KERNEL_reg_np[i], (do_perm ? permutation_kernel[i] : KERNEL_reg_np[i]), k, p_nvar, xtc[i], num_xt, l, xc[i][j], m[i][jbw], tprod, tprod_mp, pxl, (p_pxl==NULL?NULL : p_pxl+k), swap_xxt, bpso[l], do_score, perm_kbuf, use_largeh, (use_largeh ? cont_largeh_k0[i] : 0.0), invnorm, p_invnorm);
          }
        }
        else if(p_nvar == 0){
          np_convol_ckernelv(KERNEL_reg[i], xtc[i], num_xt, tprod_has_vals, xc[i][j],
                             matrix_alt_bandwidth[i], (BANDWIDTH_reg == BW_FIXED),
                             m[i][jbw],
                             tprod, bpow[i]);
          tprod_has_vals = 1;
        } else
        {
          np_convol_ckernelv(KERNEL_reg[i], xtc[i], num_xt, l, xc[i][j],
                             matrix_alt_bandwidth[i], (BANDWIDTH_reg == BW_FIXED),
                             m[i][jbw],
                             tprod, bpow[i]);
        }
        dband *= ipow(m[i][jbw], bpow[i]);

        if(do_perm){
          for(ii = 0, kk = 0; ii < num_reg_continuous; ii++){
            if(bpso[ii]){
              if (i != ii){
                p_dband[kk] *= ipow(m[i][jbw], bpow[i]);              
              } else {
                p_dband[kk] *= ipow(m[i][jbw], p_ipow);
                if(((BANDWIDTH_reg == BW_FIXED) && (int_LARGE_SF == 0) && do_score)){
                  p_dband[kk] *= vector_scale_factor[ii]/(m[i][jbw]);
                }
              }
              kk++;
            }
          }
        }
        k += bpso[l];
      }
    }


    for(ii = k; ii < p_nvar; ii++){
      p_dband[ii] = dband;
    }

    if(use_disc_profile_cache){
      int nactive = 0, nplist = 0;
      double disc_profile_const = 1.0;
      int disc_profile_const_active = 0;
      const int profile_base_empty = !tprod_has_vals;

      if(disc_uno_const_ok != NULL){
        for(kk = 0; kk < num_reg_unordered; kk++){
          if(!disc_uno_const_ok[kk]) continue;
          disc_profile_const *= disc_uno_const[kk];
          disc_profile_const_active = 1;
        }
      }

      if(disc_ord_const_ok != NULL){
        for(kk = 0; kk < num_reg_ordered; kk++){
          const int opidx = num_reg_continuous + num_reg_unordered + kk;
          if(!disc_ord_const_ok[kk]) continue;
          if(!(ps_ok_nli || (operator[opidx] != OP_CONVOLUTION))) continue;
          disc_profile_const *= disc_ord_const[kk];
          disc_profile_const_active = 1;
        }
      }

      if(profile_base_empty){
        if(pxl == NULL){
          for(i = 0; i < num_xt; i++)
            disc_active_idx[nactive++] = i;
        } else {
          for(ii = 0; ii < pxl->n; ii++){
            const int istart = pxl->istart[ii];
            const int iend = istart + pxl->nlev[ii];
            for(i = istart; i < iend; i++)
              disc_active_idx[nactive++] = i;
          }
        }
      } else if(pxl == NULL){
        for(i = 0; i < num_xt; i++){
          if(tprod[i] == 0.0) continue;
          disc_active_idx[nactive++] = i;
        }
      } else {
        for(ii = 0; ii < pxl->n; ii++){
          const int istart = pxl->istart[ii];
          const int iend = istart + pxl->nlev[ii];
          for(i = istart; i < iend; i++){
            if(tprod[i] == 0.0) continue;
            disc_active_idx[nactive++] = i;
          }
        }
      }

      if(disc_mark_token == INT_MAX){
        memset(disc_prof_mark, 0, (size_t)disc_nprof*sizeof(int));
        disc_mark_token = 1;
      } else {
        disc_mark_token++;
      }

      for(i = 0; i < nactive; i++){
        const int pid = disc_prof_id[disc_active_idx[i]];
        if(disc_prof_mark[pid] == disc_mark_token) continue;
        disc_prof_mark[pid] = disc_mark_token;
        disc_prof_list[nplist++] = pid;
      }

      double (* const ukf[])(int, double, int) = {
        np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
        np_score_uaa, np_score_unli_racine
      };
      double (* const okf[])(double, double, double, double, double) = {
        np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
      };
      const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
      const int nok = (int)(sizeof(okf)/sizeof(okf[0]));

      for(i = 0; i < nplist; i++){
        const int pid = disc_prof_list[i];
        const int ridx = disc_prof_rep[pid];
        double dprod = 1.0;

        for(kk = 0; kk < num_reg_unordered; kk++){
          if((disc_uno_const_ok != NULL) && disc_uno_const_ok[kk]) continue;
          const int ku = (KERNEL_unordered_reg_np[kk] >= 0 && KERNEL_unordered_reg_np[kk] < nuk)
            ? KERNEL_unordered_reg_np[kk] : 0;
          dprod *= ukf[ku]((xtu[kk][ridx] == xu[kk][j]), lambda[kk], num_categories[kk]);
          if(dprod == 0.0) break;
        }

        if(dprod != 0.0){
          for(kk = 0; kk < num_reg_ordered; kk++){
            const int opidx = num_reg_continuous + num_reg_unordered + kk;
            if((disc_ord_const_ok != NULL) && disc_ord_const_ok[kk] &&
               (ps_ok_nli || (operator[opidx] != OP_CONVOLUTION)))
              continue;
            const int ko = (KERNEL_ordered_reg_np[kk] >= 0 && KERNEL_ordered_reg_np[kk] < nok)
              ? KERNEL_ordered_reg_np[kk] : 0;
            const double c1 = swap_xxt ? xo[kk][j] : xto[kk][ridx];
            const double c2 = swap_xxt ? xto[kk][ridx] : xo[kk][j];
            dprod *= okf[ko](c1, c2, lambda[num_reg_unordered + kk], disc_ord_cl[kk], disc_ord_ch[kk]);
            if(dprod == 0.0) break;
          }
        }

        disc_prof_val[pid] = dprod;
      }

      if(profile_base_empty){
        for(i = 0; i < nactive; i++){
          const int idx_i = disc_active_idx[i];
          tprod[idx_i] = disc_prof_val[disc_prof_id[idx_i]];
        }
        tprod_has_vals = (nactive > 0);
      } else {
        for(i = 0; i < nactive; i++){
          const int idx_i = disc_active_idx[i];
          tprod[idx_i] *= disc_prof_val[disc_prof_id[idx_i]];
        }
      }

      if(disc_profile_const_active){
        deferred_const *= disc_profile_const;
        deferred_const_active = 1;
      }

      for(ii = 0; ii < (num_reg_unordered + num_reg_ordered); ii++){
        k += bpso[l + ii];
      }
      l += (num_reg_unordered + num_reg_ordered);
      ip += doscoreocg*(num_reg_unordered + num_reg_ordered);
    } else {
      /* unordered second */
      for(i=0; i < num_reg_unordered; i++, l++, ip += doscoreocg){
        if(doscoreocg){
          np_p_ukernelv(KERNEL_unordered_reg_np[i], ps_ukernel[i], k, p_nvar, xtu[i], num_xt, l, xu[i][j], 
                        lambda[i], num_categories[i], matrix_categorical_vals[i][0], tprod, tprod_mp, pxl, p_pxl + k, swap_xxt, (bpso[l] ? do_ocg : 0), perm_kbuf);
        } else {
          if((p_nvar == 0) && (disc_uno_const_ok != NULL) && disc_uno_const_ok[i]){
            deferred_const *= disc_uno_const[i];
            deferred_const_active = 1;
          } else if(disc_uno_const_ok != NULL && disc_uno_const_ok[i]){
            np_ckernelv_mul_const(disc_uno_const[i], num_xt, l, tprod, pxl);
          } else {
            np_ukernelv(KERNEL_unordered_reg_np[i], xtu[i], num_xt, tprod_has_vals, xu[i][j], 
                        lambda[i], num_categories[i], tprod, pxl,
                        (disable_gate_features ? 1 : (disc_uno_const_ok != NULL)));
            tprod_has_vals = 1;
          }
        }
        k += bpso[l];
      }

      /* ordered third */
      for(i=0; i < num_reg_ordered; i++, l++, ip += doscoreocg){
        if(!doscoreocg){
          if((p_nvar == 0) &&
             (disc_ord_const_ok != NULL) && disc_ord_const_ok[i] &&
             (ps_ok_nli || (operator[l] != OP_CONVOLUTION))){
            deferred_const *= disc_ord_const[i];
            deferred_const_active = 1;
          } else if((disc_ord_const_ok != NULL) && disc_ord_const_ok[i] &&
                    (ps_ok_nli || (operator[l] != OP_CONVOLUTION))){
            np_ckernelv_mul_const(disc_ord_const[i], num_xt, l, tprod, pxl);
          } else if(ps_ok_nli || (operator[l] != OP_CONVOLUTION)){
            np_okernelv(KERNEL_ordered_reg_np[i], xto[i], num_xt, tprod_has_vals,
                        xo[i][j], lambda[num_reg_unordered+i], 
                        (matrix_categorical_vals != NULL) ? matrix_categorical_vals[i+num_reg_unordered] : NULL, 
                        (num_categories != NULL) ? num_categories[i+num_reg_unordered] : 0,
                        tprod, pxl, swap_xxt);      
            tprod_has_vals = 1;
          } else {
            np_convol_okernelv(KERNEL_ordered_reg[i], xto[i], num_xt, tprod_has_vals,
                               xo[i][j], lambda[num_reg_unordered+i], 
                               num_categories[i+num_reg_unordered],
                               matrix_categorical_vals[i+num_reg_unordered],
                               tprod, swap_xxt);
            tprod_has_vals = 1;
          }
        } else {
          np_p_okernelv(KERNEL_ordered_reg_np[i], ps_okernel[i], k, p_nvar, xto[i], num_xt, l,
                        xo[i][j], lambda[num_reg_unordered+i], 
                        (matrix_categorical_vals != NULL) ? matrix_categorical_vals[i+num_reg_unordered] : NULL, 
                        (num_categories != NULL) ? num_categories[i+num_reg_unordered] : 0,
                        tprod, tprod_mp, pxl, p_pxl + k, swap_xxt, (bpso[l] ? do_ocg : 0),
                        matrix_ordered_indices[i], (swap_xxt ? 0 : matrix_ordered_indices[i][j]),
                        perm_kbuf);
        }
        k += bpso[l];
      }
    }
    if((p_nvar == 0) && deferred_const_active){
      np_ckernelv_mul_const(deferred_const, num_xt, tprod_has_vals, tprod, pxl);
      tprod_has_vals = 1;
    }

    /* expand matrix outer product, multiply by kernel weights, etc, do sum */

    if (!(drop_one_train && do_psum && (j == drop_which_train))){
      if(!nws){
        np_outer_weighted_sum(matrix_W, sgn, ncol_W, 
                              matrix_Y, ncol_Y,
                              tprod, num_xt,
                              leave_or_drop, lod,
                              kernel_pow,
                              do_psum, j,
                              symmetric,
                              gather_scatter,
                              1, dband,
                              ws, pxl);
      }


        for(ii = 0; ii < p_nvar; ii++){
          np_outer_weighted_sum(matrix_W, sgn, ncol_W, 
                                matrix_Y, ncol_Y,
                                tprod_mp+ii*num_xt, num_xt,
                                leave_or_drop, lod,
                                kernel_pow,
                                do_psum, j,
                                symmetric,
                                gather_scatter,
                                1, p_dband[ii],
                                p_ws + ii*num_obs_eval*sum_element_length, (p_pxl == NULL) ? NULL : (p_pxl+ii));
        }

    }

    if(kw_work != NULL){ 
      // if using adaptive bandwidths, kw is returned transposed
      if(bandwidth_divide_weights)
        for(i = 0; i < num_xt; i++)
          kw_work[j*num_xt + i] = tprod[i]/dband;
      else
        for(i = 0; i < num_xt; i++)
          kw_work[j*num_xt + i] = tprod[i];
    }

    if((kernel_weighted_sum_pkw_extern != NULL) && (kernel_weighted_sum_pkw_nvar_extern > 0)){
      for(ii = 0; ii < kernel_weighted_sum_pkw_nvar_extern; ii++){
        if(bandwidth_divide_weights)
          for(i = 0; i < num_xt; i++)
            kernel_weighted_sum_pkw_extern[ii*num_obs_eval*num_xt + j*num_xt + i] =
              tprod_mp[ii*num_xt + i]/p_dband[ii];
        else
          for(i = 0; i < num_xt; i++)
            kernel_weighted_sum_pkw_extern[ii*num_obs_eval*num_xt + j*num_xt + i] =
              tprod_mp[ii*num_xt + i];
      }
    }

    np_progress_fit_loop_step(j + 1, progress_total);
    
  }

  if((!gather_scatter) && (!suppress_parallel)){ 
    // gather_scatter is only used for the local-linear cv
    // note: ll cv + adaptive_nn does not work in parallel
#ifdef MPI2
    if(!nws){
      if (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN){
        MPI_Allgather(MPI_IN_PLACE, stride * sum_element_length, MPI_DOUBLE, weighted_sum, stride * sum_element_length, MPI_DOUBLE, comm[1]);
      } else if(BANDWIDTH_reg == BW_ADAP_NN){
        MPI_Allreduce(MPI_IN_PLACE, weighted_sum, num_obs_eval*sum_element_length, MPI_DOUBLE, MPI_SUM, comm[1]);
      }
    }

    if((kw_work != NULL) && (!keep_kw_owner_local)){
      MPI_Allgather(MPI_IN_PLACE, stride * num_xt, MPI_DOUBLE, kw_work, stride * num_xt, MPI_DOUBLE, comm[1]);
    }

    if(p_nvar > 0){
      if (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN){
        for(ii = 0; ii < p_nvar; ii++){
          MPI_Allgatherv(MPI_IN_PLACE, igatherv[my_rank], MPI_DOUBLE, weighted_permutation_sum + ii*num_obs_eval*sum_element_length, igatherv, idisplsv, MPI_DOUBLE, comm[1]);
        }
      } else if(BANDWIDTH_reg == BW_ADAP_NN){
        MPI_Allreduce(MPI_IN_PLACE, weighted_permutation_sum, p_nvar*num_obs_eval*sum_element_length, MPI_DOUBLE, MPI_SUM, comm[1]);
      }
    }
#endif
  }

#ifdef MPI2
  if((kw != NULL) && (kw_work != kw)){
    memcpy(kw, kw_work, (size_t)num_obs_eval*(size_t)num_xt*sizeof(double));
  }
#endif

cleanup:
#ifdef MPI2
  free(igatherv);
  free(idisplsv);
#endif

  free(KERNEL_reg_np);
  free(KERNEL_unordered_reg_np);
  free(KERNEL_ordered_reg_np);
  
  if(!bandwidth_provided){
    free(lambda);
    free_tmat(matrix_bandwidth);
  }

  if(!bandwidth_provided && ((BANDWIDTH_reg == BW_ADAP_NN) && any_convolution))
    free_tmat(matrix_alt_bandwidth);

  free(tprod);
  free(bpow);
  
  clean_xl(pxl);
  clean_nl(&nls);

  if(p_nvar > 0){
    free(tprod_mp);
    free(perm_kbuf);

    if(do_perm)
      free(permutation_kernel);

    if(do_score){
      free(ps_ukernel);
      free(ps_okernel);
    }

    if(np_ks_tree_use){
      if(permutation_operator == OP_INTEGRAL){
        for(l = 0, k = 0; l < num_reg_continuous; l++){
          if(bpso[l]){
            clean_xl(p_pxl+k);
            k++;
          }
        }
      }
      free(p_pxl);
    }

    free(p_dband);
  }

  if((disc_prof_id != NULL) && (!disc_profile_from_override) && (!disc_profile_from_global_cache)) free(disc_prof_id);
  if((disc_prof_rep != NULL) && (!disc_profile_from_override) && (!disc_profile_from_global_cache)) free(disc_prof_rep);
  if(disc_prof_hash != NULL) free(disc_prof_hash);
  if(disc_prof_mark != NULL) free(disc_prof_mark);
  if(disc_prof_list != NULL) free(disc_prof_list);
  if(disc_active_idx != NULL) free(disc_active_idx);
  if(disc_prof_val != NULL) free(disc_prof_val);
  if(disc_ord_cl != NULL) free(disc_ord_cl);
  if(disc_ord_ch != NULL) free(disc_ord_ch);
  if((disc_uno_const_ok != NULL) && (!disc_uno_from_override)) free(disc_uno_const_ok);
  if((disc_uno_const != NULL) && (!disc_uno_from_override)) free(disc_uno_const);
  if((disc_ord_const_ok != NULL) && (!disc_ord_from_override)) free(disc_ord_const_ok);
  if((disc_ord_const != NULL) && (!disc_ord_from_override)) free(disc_ord_const);
  if((cont_largeh_ok != NULL) && (!cont_largeh_from_override) && (!cont_largeh_from_global_cache)) free(cont_largeh_ok);
  if((cont_largeh_hmin != NULL) && (!cont_largeh_from_override) && (!cont_largeh_from_global_cache)) free(cont_largeh_hmin);
  if((cont_largeh_k0 != NULL) && (!cont_largeh_from_override) && (!cont_largeh_from_global_cache)) free(cont_largeh_k0);
  if(cont_largeh_active != NULL) free(cont_largeh_active);
  if(cont_largeh_active_fixed != NULL) free(cont_largeh_active_fixed);
  if(tree_active_dims != NULL) free(tree_active_dims);

#ifdef MPI2
  if((kw_work != NULL) && (kw_work != kw)) free(kw_work);
#endif

  if(no_bpso)
    free(bpso);


  return(status);
}

int kernel_weighted_sum_np_ctx(
int * KERNEL_reg,
int * KERNEL_unordered_reg,
int * KERNEL_ordered_reg,
const int BANDWIDTH_reg,
const int num_obs_train,
const int num_obs_eval,
const int num_reg_unordered,
const int num_reg_ordered,
const int num_reg_continuous,
const int leave_one_out,
const int leave_one_out_offset,
const int kernel_pow,
const int bandwidth_divide,
const int bandwidth_divide_weights,
const int symmetric,
const int gather_scatter,
const int drop_one_train,
const int drop_which_train,
const int * const operator,
const int permutation_operator,
int do_score,
int do_ocg,
int * bpso,
const int suppress_parallel,
const int ncol_Y,
const int ncol_W,
const int int_TREE,
const int do_partial_tree,
KDT * const kdt,
NL * const inl,
int * const nld,
int * const idx,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double **matrix_Y,
double **matrix_W,
double * sgn,
double *vector_scale_factor,
int bandwidth_provided,
double ** matrix_bw_train,
double ** matrix_bw_eval,
double * lambda_pre,
int *num_categories,
double **matrix_categorical_vals,
int ** matrix_ordered_indices,
double * const weighted_sum,
double * const weighted_permutation_sum,
double * const kw,
const NP_GateOverrideCtx * const gate_override_ctx){
  return kernel_weighted_sum_np_ctx_ex(
    KERNEL_reg,
    KERNEL_unordered_reg,
    KERNEL_ordered_reg,
    BANDWIDTH_reg,
    num_obs_train,
    num_obs_eval,
    num_reg_unordered,
    num_reg_ordered,
    num_reg_continuous,
    leave_one_out,
    leave_one_out_offset,
    kernel_pow,
    bandwidth_divide,
    bandwidth_divide_weights,
    symmetric,
    gather_scatter,
    drop_one_train,
    drop_which_train,
    operator,
    permutation_operator,
    do_score,
    do_ocg,
    bpso,
    suppress_parallel,
    ncol_Y,
    ncol_W,
    int_TREE,
    do_partial_tree,
    kdt,
    inl,
    nld,
    idx,
    matrix_X_unordered_train,
    matrix_X_ordered_train,
    matrix_X_continuous_train,
    matrix_X_unordered_eval,
    matrix_X_ordered_eval,
    matrix_X_continuous_eval,
    matrix_Y,
    matrix_W,
    sgn,
    vector_scale_factor,
    bandwidth_provided,
    matrix_bw_train,
    matrix_bw_eval,
    lambda_pre,
    num_categories,
    matrix_categorical_vals,
    matrix_ordered_indices,
    weighted_sum,
    weighted_permutation_sum,
    kw,
    gate_override_ctx,
    0);
}

int kernel_weighted_sum_np(
int * KERNEL_reg,
int * KERNEL_unordered_reg,
int * KERNEL_ordered_reg,
const int BANDWIDTH_reg,
const int num_obs_train,
const int num_obs_eval,
const int num_reg_unordered,
const int num_reg_ordered,
const int num_reg_continuous,
const int leave_one_out,
const int leave_one_out_offset,
const int kernel_pow,
const int bandwidth_divide,
const int bandwidth_divide_weights,
const int symmetric,
const int gather_scatter,
const int drop_one_train,
const int drop_which_train,
const int * const operator,
const int permutation_operator,
int do_score,
int do_ocg,
int * bpso,
const int suppress_parallel,
const int ncol_Y,
const int ncol_W,
const int int_TREE,
const int do_partial_tree,
KDT * const kdt,
NL * const inl,
int * const nld,
int * const idx,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double **matrix_Y,
double **matrix_W,
double * sgn,
double *vector_scale_factor,
int bandwidth_provided,
double ** matrix_bw_train,
double ** matrix_bw_eval,
double * lambda_pre,
int *num_categories,
double **matrix_categorical_vals,
int ** matrix_ordered_indices,
double * const weighted_sum,
double * const weighted_permutation_sum,
double * const kw,
double * const pkw){
  double * old_pkw = kernel_weighted_sum_pkw_extern;
  int old_pkw_nvar = kernel_weighted_sum_pkw_nvar_extern;
  int status = 0;

  kernel_weighted_sum_pkw_extern = pkw;
  kernel_weighted_sum_pkw_nvar_extern = (pkw == NULL) ? 0 : (((permutation_operator != OP_NOOP) ? num_reg_continuous : 0) + ((do_score || do_ocg) ? num_reg_unordered + num_reg_ordered : 0));

  status = kernel_weighted_sum_np_ctx(KERNEL_reg,
                                    KERNEL_unordered_reg,
                                    KERNEL_ordered_reg,
                                    BANDWIDTH_reg,
                                    num_obs_train,
                                    num_obs_eval,
                                    num_reg_unordered,
                                    num_reg_ordered,
                                    num_reg_continuous,
                                    leave_one_out,
                                    leave_one_out_offset,
                                    kernel_pow,
                                    bandwidth_divide,
                                    bandwidth_divide_weights,
                                    symmetric,
                                    gather_scatter,
                                    drop_one_train,
                                    drop_which_train,
                                    operator,
                                    permutation_operator,
                                    do_score,
                                    do_ocg,
                                    bpso,
                                    suppress_parallel,
                                    ncol_Y,
                                    ncol_W,
                                    int_TREE,
                                    do_partial_tree,
                                    kdt,
                                    inl,
                                    nld,
                                    idx,
                                    matrix_X_unordered_train,
                                    matrix_X_ordered_train,
                                    matrix_X_continuous_train,
                                    matrix_X_unordered_eval,
                                    matrix_X_ordered_eval,
                                    matrix_X_continuous_eval,
                                    matrix_Y,
                                    matrix_W,
                                    sgn,
                                    vector_scale_factor,
                                    bandwidth_provided,
                                    matrix_bw_train,
                                    matrix_bw_eval,
                                    lambda_pre,
                                    num_categories,
                                    matrix_categorical_vals,
                                    matrix_ordered_indices,
                                    weighted_sum,
                                    weighted_permutation_sum,
                                    kw,
                                    NULL);
  kernel_weighted_sum_pkw_extern = old_pkw;
  kernel_weighted_sum_pkw_nvar_extern = old_pkw_nvar;
  return status;
}

int np_kernel_estimate_con_density_categorical_convolution_cv(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int num_obs,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_Y_unordered,
double **matrix_Y_ordered,
double **matrix_Y_continuous,
double **matrix_X_unordered,
double **matrix_X_ordered,
double **matrix_X_continuous,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *cv){

  int i = 0, j = 0, k = 0;
  int l, m, n, o, ib, ob, ms, ns, os;

  const int CBW_MAXBLKLEN = 64;
  int blklen = MIN(CBW_MAXBLKLEN, num_obs);

  int blj, bli, blk;

  double * sum_ker_convol, * sum_ker_marginal, * sum_ker;

  double *lambda;
  double **matrix_bandwidth_var, **matrix_bandwidth_reg;

  double * blk_xi, * blk_xj, * blk_yij, ts;

  double (* const yck)(double) = allck[KERNEL_den];
  double (* const xck)(double) = allck[KERNEL_reg];
  double (* const yok)(double, double, double, double, double) = allok[KERNEL_ordered_den];
  double (* const xok)(double, double, double, double, double) = allok[KERNEL_ordered_reg];
  double (* const yuk)(int, double, int) = alluk[KERNEL_unordered_den];
  double (* const xuk)(int, double, int) = alluk[KERNEL_unordered_reg];

  // load balancing / allocation
  // very simple : set k,j,i
#ifdef MPI2

  double * sum_ker_convolf, * sum_ker_marginalf, *sum_kerf;

  const uint64_t Np = (uint64_t)iNum_Processors;
  const uint64_t P = (uint64_t)my_rank;

  const uint64_t Nb = num_obs/blklen + (num_obs%blklen != 0);
  const uint64_t Q = (Nb*Nb*Nb)/Np;
  const uint64_t Nr = (Nb*Nb*Nb)%Np;

  const uint64_t Wi = P*Q + ((P > Nr)?Nr:P);
  const uint64_t Wf = Wi + Q + (P < Nr);
  const uint64_t uki = Wi/(Nb*Nb);
  const uint64_t uji = (Wi%(Nb*Nb))/Nb;
  const uint64_t uii = Wi%Nb;

  uint64_t W = Wi;

  k = (int)(uki*(uint64_t)blklen);
  j = (int)(uji*(uint64_t)blklen);
  i = (int)(uii*(uint64_t)blklen);

  sum_kerf = (double *)malloc(num_obs*sizeof(double));
  if(!(sum_kerf != NULL))
    error("!(sum_kerf != NULL)");

  sum_ker_convolf = (double *)malloc(num_obs*sizeof(double));
  if(!(sum_ker_convolf != NULL))
    error("!(sum_ker_convolf != NULL)");

  sum_ker_marginalf = (double *)malloc(num_obs*sizeof(double));
  if(!(sum_ker_marginalf != NULL))
    error("!(sum_ker_marginalf != NULL)");
#endif

  lambda = alloc_vecd(num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered);
  matrix_bandwidth_var = alloc_matd(num_obs,num_var_continuous);
  matrix_bandwidth_reg = alloc_matd(num_obs,num_reg_continuous);

  if(kernel_bandwidth_mean(KERNEL_den,
                           BANDWIDTH_den,
                           num_obs,
                           num_obs,
                           num_var_continuous,
                           num_var_unordered,
                           num_var_ordered,
                           num_reg_continuous,
                           num_reg_unordered,
                           num_reg_ordered,
                           0, // do not suppress_parallel
                           vector_scale_factor,
                           matrix_Y_continuous,
                           matrix_Y_continuous,
                           matrix_X_continuous,
                           matrix_X_continuous,
                           matrix_bandwidth_var,
                           matrix_bandwidth_reg,
                           lambda)==1){
    free(lambda);
    free_mat(matrix_bandwidth_var,num_var_continuous);
    free_mat(matrix_bandwidth_reg,num_reg_continuous);
    return(1);
  }


  blk_xi = (double *)malloc(sizeof(double)*blklen*blklen);
  if(!(blk_xi != NULL)) 
    error("!(blk_xi != NULL)");

  blk_xj = (double *)malloc(sizeof(double)*blklen*blklen);
  if(!(blk_xj != NULL))
    error("!(blk_xj != NULL)");

  blk_yij = (double *)malloc(sizeof(double)*blklen*blklen);
  if(!(blk_yij != NULL)) 
    error("!(blk_yij != NULL)");

  sum_ker = (double *)malloc(num_obs*sizeof(double));
  if(!(sum_ker != NULL))
    error("!(sum_ker != NULL)");

  sum_ker_convol = (double *)malloc(num_obs*sizeof(double));
  if(!(sum_ker_convol != NULL))
    error("!(sum_ker_convol != NULL)");

  sum_ker_marginal = (double *)malloc(num_obs*sizeof(double));
  if(!(sum_ker_marginal != NULL))
    error("!(sum_ker_marginal != NULL)");

  if(!(BANDWIDTH_den == BW_FIXED))
    error("!(BANDWIDTH_den == BW_FIXED)");

  for(int ii = 0; ii < num_obs; ii++){
    sum_ker[ii] = 0.0;
    sum_ker_convol[ii] = 0.0;
    sum_ker_marginal[ii] = 0.0;
  }

  *cv = 0.0;

  // top level loop corresponds to a block of X_l's 
  // (leave one out evaluation points)
#ifdef MPI2
  for(; (k < num_obs) && (W < Wf); k+=blklen)
#else
  for(; k < num_obs; k+=blklen)
#endif
    {
    blk = MIN(num_obs-blklen, k);

    // we generate blocks of X_j, X_i, and Y_ji
    // outermost loop determines the X_j's
    // innermost loop we generate new X_i and Y_ji blocks

    // one improvement would be a flip-flop algorithm where the roles
    // of X_j and X_i are swapped and only a new Y_ji need be generated
    R_CheckUserInterrupt();
#ifdef MPI2
    for(; (j < num_obs) && (W < Wf); j+=blklen)
#else
    for(; j < num_obs; j+=blklen)
#endif
      {
      blj = MIN(num_obs-blklen, j);

      // first: fill in blk_xj array with kernel evals
      // k(xl-xj)
      for(n=0, ib=0; n < blklen; n++){
        for(m=0; m < blklen; m++,ib++){
          blk_xj[ib] = 1.0;


          for(l = 0; l < num_reg_continuous; l++){
            blk_xj[ib] *= xck((matrix_X_continuous[l][blk+n]-matrix_X_continuous[l][blj+m])/matrix_bandwidth_reg[l][0])/matrix_bandwidth_reg[l][0];
          }

          for(l = 0; l < num_reg_unordered; l++){
            blk_xj[ib] *= xuk(matrix_X_unordered[l][blk+n]==matrix_X_unordered[l][blj+m],
                              lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
          }

          for(l = 0; l < num_reg_ordered; l++){
            const int lcat = l+num_var_unordered+num_var_ordered+num_reg_unordered;
            const double cl = matrix_categorical_vals[lcat][0];
            const double ch = matrix_categorical_vals[lcat][num_categories[lcat]-1];

            blk_xj[ib] *= xok(matrix_X_ordered[l][blk+n],matrix_X_ordered[l][blj+m],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered],cl,ch);
          }

          // k(yj-yi)
          ts = 1.0;

          for(l = 0; l < num_var_continuous; l++){
            ts *= yck((matrix_Y_continuous[l][blk+n]-matrix_Y_continuous[l][blj+m])/matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];
          }


          for(l = 0; l < num_var_unordered; l++){
            ts *= yuk(matrix_Y_unordered[l][blk+n]==matrix_Y_unordered[l][blj+m],lambda[l],num_categories[l]);
          }

          for(l = 0; l < num_var_ordered; l++){
            const int lcat = l+num_var_unordered;
            const double cl = matrix_categorical_vals[lcat][0];
            const double ch = matrix_categorical_vals[lcat][num_categories[lcat]-1];

            ts *= yok(matrix_Y_ordered[l][blk+n],matrix_Y_ordered[l][blj+m],lambda[l+num_var_unordered],cl,ch);
          }

          // accumulate marginals and kernel sums along the way
          if((blk+n >= k) && (blj+m >= j) && (blk+n != blj+m)){
            sum_ker[blk+n] += blk_xj[ib]*ts;
            sum_ker_marginal[blk+n] += blk_xj[ib];
          }

        }
      }

#ifdef MPI2
      for(; (i < num_obs) && (W < Wf); i+=blklen, W++)
#else
      for(; i < num_obs; i+=blklen)
#endif      
        {
        bli = MIN(num_obs-blklen, i);
        
        // second: fill in k(xl-xi) and k(yj-yi) arrays
        for(n=0, ib=0; n < blklen; n++){
          for(m=0; m < blklen; m++,ib++){
            // k(xl-xi) 
            blk_xi[ib] = 1.0;

            for(l = 0; l < num_reg_continuous; l++){
              blk_xi[ib] *= xck((matrix_X_continuous[l][blk+n]-matrix_X_continuous[l][bli+m])/matrix_bandwidth_reg[l][0])/matrix_bandwidth_reg[l][0];
            }

            for(l = 0; l < num_reg_unordered; l++){
              blk_xi[ib] *= xuk(matrix_X_unordered[l][blk+n]==matrix_X_unordered[l][bli+m],
                                lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
            }

            for(l = 0; l < num_reg_ordered; l++){
              const int lcat = l+num_var_unordered+num_var_ordered+num_reg_unordered;
              const double cl = matrix_categorical_vals[lcat][0];
              const double ch = matrix_categorical_vals[lcat][num_categories[lcat]-1];

              blk_xi[ib] *= xok(matrix_X_ordered[l][blk+n],matrix_X_ordered[l][bli+m],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered],cl,ch);
            }

            // k(2)(yj-yi)
            blk_yij[ib] = 1.0;

            for(l = 0; l < num_var_continuous; l++){
              blk_yij[ib] *= kernel_convol(KERNEL_den,BANDWIDTH_den,
                                           (matrix_Y_continuous[l][blj+n]-matrix_Y_continuous[l][bli+m])/matrix_bandwidth_var[l][0],matrix_bandwidth_var[l][0],
                                           matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];
            }

            for(l = 0; l < num_var_unordered; l++){
              blk_yij[ib] *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_Y_unordered[l][blj+n],
                                                          matrix_Y_unordered[l][bli+m],lambda[l], num_categories[l], matrix_categorical_vals[l]);
            }

            for(l = 0; l < num_var_ordered; l++){
              blk_yij[ib] *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_Y_ordered[l][blj+n],matrix_Y_ordered[l][bli+m], lambda[l+num_var_unordered], 
                                                        num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
            }
          }
        }

        // adjustments for re-centering on num_obs % 64 != 0 data
        ms = i-bli;
        ns = j-blj;
        os = k-blk;

        // third: do partial convolution
        // here there be dragons

        for(o=os, ob=os*blklen; o < blklen; o++, ob+=blklen){ //controls l-indexing
          for(n=ns, ib=ns*blklen; n < blklen; n++, ib+=blklen){ //controls j-indexing
            if(blj+n == blk+o)
              continue;

            ts = 0.0;
            for(m=ms; m < blklen; m++){ //controls i-indexing
              if(bli+m == blk+o)
                continue;
              ts += blk_xi[ob+m]*blk_yij[ib+m];
            }
            sum_ker_convol[blk+o] += blk_xj[ob+n]*ts;
          }
        }
      }
      i = 0;
    }
    j = 0;
  }

#ifdef MPI2
  MPI_Allreduce(sum_ker, sum_kerf, num_obs, MPI_DOUBLE, MPI_SUM, comm[1]);
  MPI_Allreduce(sum_ker_convol, sum_ker_convolf, num_obs, MPI_DOUBLE, MPI_SUM, comm[1]);
  MPI_Allreduce(sum_ker_marginal, sum_ker_marginalf, num_obs, MPI_DOUBLE, MPI_SUM, comm[1]);

  for(j = 0; j < num_obs; j++){
    /*    if(sum_ker_marginalf[j] <= 0.0){
      *cv = DBL_MAX;
      break;
      } jracine 16/05/10, test for zero respecting sign */
    sum_ker_marginalf[j] =  NZD_POS(sum_ker_marginalf[j]);
    *cv += (sum_ker_convolf[j]/sum_ker_marginalf[j]-2.0*sum_kerf[j])/sum_ker_marginalf[j];
  }
#else
  for(j = 0; j < num_obs; j++){
    /*    if(sum_ker_marginal[j] <= 0.0){
      *cv = DBL_MAX;
      break;
      } jracine 16/05/10 */
    sum_ker_marginal[j] =  NZD_POS(sum_ker_marginal[j]);
    *cv += (sum_ker_convol[j]/sum_ker_marginal[j]-2.0*sum_ker[j])/sum_ker_marginal[j];
  }

#endif
  if (*cv != DBL_MAX)
    *cv /= (double) num_obs;
    
  free(lambda);
  free_mat(matrix_bandwidth_var,num_var_continuous);
  free_mat(matrix_bandwidth_reg,num_reg_continuous);

  free(blk_xi);
  free(blk_xj);
  free(blk_yij);

  free(sum_ker);
  free(sum_ker_convol);
  free(sum_ker_marginal);

#ifdef MPI2
  free(sum_kerf);
  free(sum_ker_convolf);
  free(sum_ker_marginalf);
#endif

  return(0);
}

static int np_glp_max_degree(const int ncon, const int *deg){
  int j, dmax = 0;
  if((ncon <= 0) || (deg == NULL)) return 0;
  for(j = 0; j < ncon; j++)
    if(deg[j] > dmax) dmax = deg[j];
  return dmax;
}

static int np_glp_store_term(const int ncon,
                             int **terms_ptr,
                             int *nterms,
                             int *cap,
                             const int *cur){
  int j, k;
  if(*nterms >= *cap){
    int newcap = (*cap <= 0) ? 32 : 2*(*cap);
    int *tmp = (int *)realloc(*terms_ptr, (size_t)newcap*(size_t)ncon*sizeof(int));
    if(tmp == NULL) return 0;
    *terms_ptr = tmp;
    *cap = newcap;
  }
  k = (*nterms)*ncon;
  for(j = 0; j < ncon; j++)
    (*terms_ptr)[k + j] = cur[j];
  (*nterms)++;
  return 1;
}

static int np_glp_enum_terms_rec(const int idx,
                                 const int ncon,
                                 const int basis_mode,
                                 const int dmax,
                                 const int *deg,
                                 int *cur,
                                 int **terms_ptr,
                                 int *nterms,
                                 int *cap){
  int k;
  if(idx == ncon){
    int s = 0;
    int nz = 0;
    for(k = 0; k < ncon; k++)
      s += cur[k];
    for(k = 0; k < ncon; k++)
      if(cur[k] > 0) nz++;

    if(basis_mode == 2){ /* tensor */
      if(s <= 0) return 1;
    } else if(basis_mode == 0){ /* additive */
      if((s <= 0) || (nz != 1)) return 1;
    } else { /* glp */
      if((s <= 0) || (s > dmax)) return 1;
      for(k = 0; k < ncon; k++){
        const int dk = deg[k];
        if((dk > 0) && (dk < dmax) && (s > dk) && (cur[k] == dk))
          return 1;
      }
    }

    if(!np_glp_store_term(ncon, terms_ptr, nterms, cap, cur))
      return 0;
    return 1;
  }

  for(k = 0; k <= deg[idx]; k++){
    cur[idx] = k;
    if(np_glp_enum_terms_rec(idx + 1, ncon, basis_mode, dmax, deg, cur, terms_ptr, nterms, cap) == 0)
      return 0;
  }
  return 1;
}

static int np_glp_build_terms(const int ncon,
                              const int *deg,
                              const int basis_mode,
                              int **terms_out,
                              int *nterms_out){
  int dmax, j;
  double tcount_bound = 1.0;
  int *terms = NULL;
  int *cur = NULL;
  int *zero = NULL;
  int nterms = 0, cap = 0;

  *terms_out = NULL;
  *nterms_out = 0;

  if(ncon <= 0){
    *nterms_out = 1;
    return 1;
  }

  for(j = 0; j < ncon; j++){
    if((deg[j] < 0) || (deg[j] > 12))
      return 0;
    tcount_bound *= (double)(deg[j] + 1);
    if(tcount_bound > 100000.0)
      return 0;
  }

  dmax = np_glp_max_degree(ncon, deg);
  if(dmax > 12)
    return 0;

  zero = (int *)calloc((size_t)ncon, sizeof(int));
  if(zero == NULL)
    return 0;
  if(!np_glp_store_term(ncon, &terms, &nterms, &cap, zero)){
    free(zero);
    return 0;
  }
  free(zero);

  cur = (int *)calloc((size_t)ncon, sizeof(int));
  if(cur == NULL){
    free(terms);
    return 0;
  }

  if(dmax > 0){
    if(np_glp_enum_terms_rec(0, ncon, basis_mode, dmax, deg, cur, &terms, &nterms, &cap) == 0){
      free(cur);
      free(terms);
      return 0;
    }
  }

  free(cur);
  *terms_out = terms;
  *nterms_out = nterms;
  return 1;
}

static void np_glp_fill_basis_raw_train(const int ncon,
                                        const int *terms,
                                        const int nterms,
                                        double **matrix_X_continuous_train,
                                        const int num_obs_train,
                                        double **basis){
  int t, i, j;
  for(t = 0; t < nterms; t++){
    const int *et = terms + t*ncon;
    for(i = 0; i < num_obs_train; i++){
      double b = 1.0;
      for(j = 0; j < ncon; j++){
        const int p = et[j];
        if(p > 0)
          b *= ipow(matrix_X_continuous_train[j][i], p);
      }
      basis[t][i] = b;
    }
  }
}

static void np_glp_fill_basis_eval_raw(const int ncon,
                                       const int *terms,
                                       const int nterms,
                                       double **matrix_X_continuous_eval,
                                       const int eval_index,
                                       double *eval_basis){
  int t, j;
  for(t = 0; t < nterms; t++){
    const int *et = terms + t*ncon;
    double b = 1.0;
    for(j = 0; j < ncon; j++){
      const int p = et[j];
      if(p > 0)
        b *= ipow(matrix_X_continuous_eval[j][eval_index], p);
    }
    eval_basis[t] = b;
  }
}

static inline double np_glp_falling_factorial(const int p, const int r){
  int k;
  double c = 1.0;
  if(r <= 0) return 1.0;
  if(p < r) return 0.0;
  for(k = 0; k < r; k++)
    c *= (double)(p - k);
  return c;
}

static void np_glp_fill_basis_eval_deriv_raw(const int which_var,
                                             const int deriv_order,
                                             const int ncon,
                                             const int *terms,
                                             const int nterms,
                                             double **matrix_X_continuous_eval,
                                             const int eval_index,
                                             double *eval_deriv){
  int t, j;
  for(t = 0; t < nterms; t++){
    const int *et = terms + t*ncon;
    double b = 1.0;
    for(j = 0; j < ncon; j++){
      const int p = et[j];
      const double x = matrix_X_continuous_eval[j][eval_index];
      if(j == which_var){
        if(p < deriv_order){
          b = 0.0;
          break;
        }
        b *= np_glp_falling_factorial(p, deriv_order)*ipow(x, p - deriv_order);
      } else if(p > 0){
        b *= ipow(x, p);
      }
    }
    eval_deriv[t] = b;
  }
}

typedef struct {
  int degree;
  int use_basis;
  int max_deriv;
  double xmin;
  double xmax;
  gsl_bspline_workspace *bw;
  gsl_bspline_deriv_workspace *dw;
  gsl_vector *B;
  gsl_matrix *dB;
} NPGLPBasisCtx;

static void np_glp_basis_ctx_free(NPGLPBasisCtx *ctx);

static int np_glp_basis_ctx_init(NPGLPBasisCtx *ctx, const int degree, const double xmin, const double xmax){
  const size_t k = (size_t)(degree + 1);
  const size_t nbreak = 2;
  const size_t ncoeff = (size_t)(degree + 1);
  const size_t nderiv = (size_t)(degree + 1);
  ctx->degree = degree;
  ctx->use_basis = (degree > 0) && (xmax > xmin);
  ctx->max_deriv = degree;
  ctx->xmin = xmin;
  ctx->xmax = xmax;
  ctx->bw = NULL;
  ctx->dw = NULL;
  ctx->B = NULL;
  ctx->dB = NULL;
  if(!ctx->use_basis) return 1;
  ctx->bw = gsl_bspline_alloc(k, nbreak);
  if(ctx->bw == NULL) return 0;
  if(gsl_bspline_knots_uniform(xmin, xmax, ctx->bw) != GSL_SUCCESS){
    np_glp_basis_ctx_free(ctx);
    return 0;
  }
  ctx->B = gsl_vector_alloc(ncoeff);
  ctx->dw = gsl_bspline_deriv_alloc(k);
  ctx->dB = gsl_matrix_alloc(ncoeff, nderiv);
  if((ctx->B == NULL) || (ctx->dw == NULL) || (ctx->dB == NULL)){
    np_glp_basis_ctx_free(ctx);
    return 0;
  }
  return 1;
}

static void np_glp_basis_ctx_free(NPGLPBasisCtx *ctx){
  if(ctx->dB != NULL) gsl_matrix_free(ctx->dB);
  if(ctx->dw != NULL) gsl_bspline_deriv_free(ctx->dw);
  if(ctx->B != NULL) gsl_vector_free(ctx->B);
  if(ctx->bw != NULL) gsl_bspline_free(ctx->bw);
  ctx->dB = NULL;
  ctx->dw = NULL;
  ctx->B = NULL;
  ctx->bw = NULL;
}

static inline double np_glp_basis_factor_value(const NPGLPBasisCtx *ctx, const int idx){
  if(idx == 0) return 1.0;
  if((ctx->degree <= 0) || (idx < 1) || (idx > ctx->degree) || !ctx->use_basis) return 0.0;
  return gsl_vector_get(ctx->B, (size_t)idx);
}

static void np_glp_fill_basis_train(const int ncon,
                                    const int *terms,
                                    const int nterms,
                                    double **matrix_X_continuous_train,
                                    const int num_obs_train,
                                    NPGLPBasisCtx *ctx,
                                    double **basis){
  int t, i, j;
  for(t = 0; t < nterms; t++){
    const int *et = terms + t*ncon;
    for(i = 0; i < num_obs_train; i++){
      double b = 1.0;
      for(j = 0; j < ncon; j++){
        if(ctx[j].use_basis)
          gsl_bspline_eval(matrix_X_continuous_train[j][i], ctx[j].B, ctx[j].bw);
        b *= np_glp_basis_factor_value(&ctx[j], et[j]);
      }
      basis[t][i] = b;
    }
  }
}

static void np_glp_fill_basis_eval(const int ncon,
                                   const int *terms,
                                   const int nterms,
                                   double **matrix_X_continuous_eval,
                                   const int eval_index,
                                   NPGLPBasisCtx *ctx,
                                   double *eval_basis){
  int t, j;
  for(j = 0; j < ncon; j++)
    if(ctx[j].use_basis)
      gsl_bspline_eval(matrix_X_continuous_eval[j][eval_index], ctx[j].B, ctx[j].bw);
  for(t = 0; t < nterms; t++){
    const int *et = terms + t*ncon;
    double b = 1.0;
    for(j = 0; j < ncon; j++)
      b *= np_glp_basis_factor_value(&ctx[j], et[j]);
    eval_basis[t] = b;
  }
}

static void np_glp_fill_basis_eval_deriv(const int which_var,
                                         const int deriv_order,
                                         const int ncon,
                                         const int *terms,
                                         const int nterms,
                                         double **matrix_X_continuous_eval,
                                         const int eval_index,
                                         NPGLPBasisCtx *ctx,
                                         double *eval_deriv){
  int t, j;
  for(j = 0; j < ncon; j++)
    if(ctx[j].use_basis)
      gsl_bspline_eval(matrix_X_continuous_eval[j][eval_index], ctx[j].B, ctx[j].bw);
  if((deriv_order > 0) && ctx[which_var].use_basis){
    if(deriv_order <= ctx[which_var].max_deriv){
      gsl_bspline_deriv_eval(matrix_X_continuous_eval[which_var][eval_index],
                             deriv_order,
                             ctx[which_var].dB,
                             ctx[which_var].bw,
                             ctx[which_var].dw);
    }
  }
  for(t = 0; t < nterms; t++){
    const int *et = terms + t*ncon;
    double b = 1.0;
    for(j = 0; j < ncon; j++){
      const int idx = et[j];
      if(j == which_var){
        if(deriv_order <= 0){
          b *= np_glp_basis_factor_value(&ctx[j], idx);
        } else if(idx == 0){
          b = 0.0;
          break;
        } else if((ctx[j].degree <= 0) || (idx < 1) || (idx > ctx[j].degree) || !ctx[j].use_basis || (deriv_order > ctx[j].max_deriv)){
          b = 0.0;
          break;
        } else {
          b *= gsl_matrix_get(ctx[j].dB, (size_t)idx, (size_t)deriv_order);
        }
      } else {
        b *= np_glp_basis_factor_value(&ctx[j], idx);
      }
    }
    eval_deriv[t] = b;
  }
}

typedef struct {
  int ready;
  int use_bernstein;
  int basis_mode;
  int num_obs;
  int ncon;
  int nterms;
  int *terms;
  double **basis;
  NPGLPBasisCtx *basis_ctx;
  double **matrix_X_continuous_train_ptr;
} NPGLPCVCache;

static NPGLPCVCache np_glp_cv_cache = {0, 0, 1, 0, 0, 0, NULL, NULL, NULL, NULL};

static void np_glp_cv_cache_clear(void){
  int l;
  if(np_glp_cv_cache.basis != NULL)
    free_mat(np_glp_cv_cache.basis, np_glp_cv_cache.nterms);
  if(np_glp_cv_cache.basis_ctx != NULL){
    for(l = 0; l < np_glp_cv_cache.ncon; l++)
      np_glp_basis_ctx_free(&np_glp_cv_cache.basis_ctx[l]);
    free(np_glp_cv_cache.basis_ctx);
  }
  free(np_glp_cv_cache.terms);
  np_glp_cv_cache.ready = 0;
  np_glp_cv_cache.use_bernstein = 0;
  np_glp_cv_cache.basis_mode = 1;
  np_glp_cv_cache.num_obs = 0;
  np_glp_cv_cache.ncon = 0;
  np_glp_cv_cache.nterms = 0;
  np_glp_cv_cache.terms = NULL;
  np_glp_cv_cache.basis = NULL;
  np_glp_cv_cache.basis_ctx = NULL;
  np_glp_cv_cache.matrix_X_continuous_train_ptr = NULL;
}

static int np_glp_cv_cache_prepare(const int int_ll,
                                   const int num_obs,
                                   const int ncon,
                                   double **matrix_X_continuous_train){
  int l, i;
  int *terms = NULL;
  int nterms = 0;
  int implicit_degree[MAX(1, ncon)];
  const int *degree_vec = vector_glp_degree_extern;
  const int use_bernstein =
    ((int_ll == LL_LP) && (int_ll_extern == LL_LL)) ? 0 : (int_glp_bernstein_extern != 0);
  const int basis_mode =
    ((int_ll == LL_LP) && (int_ll_extern == LL_LL)) ? 1 : int_glp_basis_extern;
  double **basis = NULL;
  NPGLPBasisCtx *basis_ctx = NULL;

  np_glp_cv_cache_clear();

  if(int_ll != LL_LP) return 1;
  if((int_ll_extern == LL_LL) && (ncon > 0)){
    for(l = 0; l < ncon; l++)
      implicit_degree[l] = 1;
    degree_vec = implicit_degree;
  }
  if((degree_vec == NULL) || (ncon <= 0) || (num_obs <= 0))
    return 0;
  if(!np_glp_build_terms(ncon, degree_vec, basis_mode, &terms, &nterms))
    return 0;
  if(nterms <= 0){
    free(terms);
    return 0;
  }

  basis = alloc_matd(num_obs, nterms);
  if(basis == NULL){
    free(terms);
    return 0;
  }

  if(use_bernstein){
    basis_ctx = (NPGLPBasisCtx *)calloc((size_t)ncon, sizeof(NPGLPBasisCtx));
    if(basis_ctx == NULL){
      free_mat(basis, nterms);
      free(terms);
      return 0;
    }
    for(l = 0; l < ncon; l++){
      double xmin = matrix_X_continuous_train[l][0];
      double xmax = matrix_X_continuous_train[l][0];
      for(i = 1; i < num_obs; i++){
        const double xi = matrix_X_continuous_train[l][i];
        if(xi < xmin) xmin = xi;
        if(xi > xmax) xmax = xi;
      }
      if(!np_glp_basis_ctx_init(&basis_ctx[l], degree_vec[l], xmin, xmax)){
        for(i = 0; i <= l; i++) np_glp_basis_ctx_free(&basis_ctx[i]);
        free(basis_ctx);
        free_mat(basis, nterms);
        free(terms);
        return 0;
      }
    }
    np_glp_fill_basis_train(ncon, terms, nterms, matrix_X_continuous_train, num_obs, basis_ctx, basis);
  } else {
    np_glp_fill_basis_raw_train(ncon, terms, nterms, matrix_X_continuous_train, num_obs, basis);
  }

  np_glp_cv_cache.ready = 1;
  np_glp_cv_cache.use_bernstein = use_bernstein;
  np_glp_cv_cache.basis_mode = basis_mode;
  np_glp_cv_cache.num_obs = num_obs;
  np_glp_cv_cache.ncon = ncon;
  np_glp_cv_cache.nterms = nterms;
  np_glp_cv_cache.terms = terms;
  np_glp_cv_cache.basis = basis;
  np_glp_cv_cache.basis_ctx = basis_ctx;
  np_glp_cv_cache.matrix_X_continuous_train_ptr = matrix_X_continuous_train;
  return 1;
}

int np_glp_cv_prepare_extern(const int int_ll,
                             const int num_obs,
                             const int ncon,
                             double **matrix_X_continuous_train){
  return np_glp_cv_cache_prepare(int_ll, num_obs, ncon, matrix_X_continuous_train);
}

void np_glp_cv_clear_extern(void){
  np_glp_cv_cache_clear();
}

typedef struct {
  int ready;
  int ncon;
  int nuno;
  int nord;
  int bw_rows;
  int total_reg;
  int kernel_reg;
  int kernel_unordered_reg;
  int kernel_ordered_reg;
  int *operator;
  int *kernel_c;
  int *kernel_u;
  int *kernel_o;
  double *lambda;
  double **matrix_bandwidth;
} NPRegCVCachedCore;

static NPRegCVCachedCore np_reg_cv_core_cache =
  {0, 0, 0, 0, 0, 0, -1, -1, -1, NULL, NULL, NULL, NULL, NULL, NULL};

static void np_reg_cv_core_cache_clear(void){
  if(np_reg_cv_core_cache.operator != NULL) free(np_reg_cv_core_cache.operator);
  if(np_reg_cv_core_cache.kernel_c != NULL) free(np_reg_cv_core_cache.kernel_c);
  if(np_reg_cv_core_cache.kernel_u != NULL) free(np_reg_cv_core_cache.kernel_u);
  if(np_reg_cv_core_cache.kernel_o != NULL) free(np_reg_cv_core_cache.kernel_o);
  if(np_reg_cv_core_cache.lambda != NULL) free(np_reg_cv_core_cache.lambda);
  if(np_reg_cv_core_cache.matrix_bandwidth != NULL) free_tmat(np_reg_cv_core_cache.matrix_bandwidth);
  np_reg_cv_core_cache.ready = 0;
  np_reg_cv_core_cache.ncon = 0;
  np_reg_cv_core_cache.nuno = 0;
  np_reg_cv_core_cache.nord = 0;
  np_reg_cv_core_cache.bw_rows = 0;
  np_reg_cv_core_cache.total_reg = 0;
  np_reg_cv_core_cache.kernel_reg = -1;
  np_reg_cv_core_cache.kernel_unordered_reg = -1;
  np_reg_cv_core_cache.kernel_ordered_reg = -1;
  np_reg_cv_core_cache.operator = NULL;
  np_reg_cv_core_cache.kernel_c = NULL;
  np_reg_cv_core_cache.kernel_u = NULL;
  np_reg_cv_core_cache.kernel_o = NULL;
  np_reg_cv_core_cache.lambda = NULL;
  np_reg_cv_core_cache.matrix_bandwidth = NULL;
}

static int np_reg_cv_core_cache_prepare(const int KERNEL_reg,
                                        const int KERNEL_unordered_reg,
                                        const int KERNEL_ordered_reg,
                                        const int BANDWIDTH_reg,
                                        const int num_obs,
                                        const int num_reg_continuous,
                                        const int num_reg_unordered,
                                        const int num_reg_ordered){
  const int total_reg = num_reg_continuous + num_reg_unordered + num_reg_ordered;
  const int bw_rows = (BANDWIDTH_reg == BW_FIXED) ? 1 : num_obs;
  int i;
  int need_realloc = 0;

  if(!np_reg_cv_core_cache.ready) {
    need_realloc = 1;
  } else if(np_reg_cv_core_cache.ncon != num_reg_continuous ||
            np_reg_cv_core_cache.nuno != num_reg_unordered ||
            np_reg_cv_core_cache.nord != num_reg_ordered ||
            np_reg_cv_core_cache.bw_rows != bw_rows) {
    need_realloc = 1;
  }

  if(need_realloc){
    np_reg_cv_core_cache_clear();
    np_reg_cv_core_cache.operator = (int *)malloc(sizeof(int)*MAX(1,total_reg));
    np_reg_cv_core_cache.kernel_c = (int *)malloc(sizeof(int)*MAX(1,num_reg_continuous));
    np_reg_cv_core_cache.kernel_u = (int *)malloc(sizeof(int)*MAX(1,num_reg_unordered));
    np_reg_cv_core_cache.kernel_o = (int *)malloc(sizeof(int)*MAX(1,num_reg_ordered));
    np_reg_cv_core_cache.lambda = alloc_vecd(MAX(1, num_reg_unordered + num_reg_ordered));
    np_reg_cv_core_cache.matrix_bandwidth = alloc_tmatd(bw_rows, num_reg_continuous);
    if(np_reg_cv_core_cache.operator == NULL ||
       np_reg_cv_core_cache.kernel_c == NULL ||
       np_reg_cv_core_cache.kernel_u == NULL ||
       np_reg_cv_core_cache.kernel_o == NULL ||
       np_reg_cv_core_cache.lambda == NULL ||
       np_reg_cv_core_cache.matrix_bandwidth == NULL){
      np_reg_cv_core_cache_clear();
      return 0;
    }
    np_reg_cv_core_cache.ready = 1;
    np_reg_cv_core_cache.ncon = num_reg_continuous;
    np_reg_cv_core_cache.nuno = num_reg_unordered;
    np_reg_cv_core_cache.nord = num_reg_ordered;
    np_reg_cv_core_cache.bw_rows = bw_rows;
    np_reg_cv_core_cache.total_reg = total_reg;
  }

  np_reg_cv_core_cache.kernel_reg = KERNEL_reg;
  np_reg_cv_core_cache.kernel_unordered_reg = KERNEL_unordered_reg;
  np_reg_cv_core_cache.kernel_ordered_reg = KERNEL_ordered_reg;

  for(i = 0; i < total_reg; i++)
    np_reg_cv_core_cache.operator[i] = OP_NORMAL;
  for(i = 0; i < num_reg_continuous; i++)
    np_reg_cv_core_cache.kernel_c[i] = KERNEL_reg;
  for(i = 0; i < num_reg_unordered; i++)
    np_reg_cv_core_cache.kernel_u[i] = KERNEL_unordered_reg;
  for(i = 0; i < num_reg_ordered; i++)
    np_reg_cv_core_cache.kernel_o[i] = KERNEL_ordered_reg;

  return 1;
}

void np_reg_cv_core_clear_extern(void){
  np_reg_cv_core_cache_clear();
}

static inline int np_fastcv_disc_unordered_all_large(const int num_reg_unordered,
                                                     const int * const kernel_u,
                                                     const double * const lambda,
                                                     const int * const num_categories){
  int i;
  double (* const ukf[])(int, double, int) = {
    np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
    np_score_uaa, np_score_unli_racine
  };
  const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));

  if(num_reg_unordered <= 0)
    return 1;
  if((kernel_u == NULL) || (lambda == NULL))
    return 0;

  for(i = 0; i < num_reg_unordered; i++){
    const int ku = kernel_u[i];
    const int ncat = (num_categories != NULL) ? num_categories[i] : 0;
    const double lam = lambda[i];
    double ks, kd;

    if(ku < 0 || ku >= nuk)
      return 0;
    if(!np_disc_near_upper(ku, lam, ncat))
      return 0;

    ks = ukf[ku](1, lam, ncat);
    kd = ukf[ku](0, lam, ncat);
    if(!np_disc_near_const_kernel(ks, kd))
      return 0;
  }

  return 1;
}

static inline int np_fastcv_disc_ordered_all_large(const int num_reg_unordered,
                                                   const int num_reg_ordered,
                                                   const int * const kernel_o,
                                                   const double * const lambda,
                                                   const int * const num_categories,
                                                   double **matrix_categorical_vals){
  int i;
  double (* const okf[])(double, double, double, double, double) = {
    np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
    np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
    np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
    np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
  };
  const int nok = (int)(sizeof(okf)/sizeof(okf[0]));

  if(num_reg_ordered <= 0)
    return 1;
  if((kernel_o == NULL) || (lambda == NULL) || (matrix_categorical_vals == NULL))
    return 0;

  for(i = 0; i < num_reg_ordered; i++){
    const int oi = i + num_reg_unordered;
    const int ko = kernel_o[i];
    const int ncat = (num_categories != NULL) ? num_categories[oi] : 0;
    const double lam = lambda[oi];
    double cl, ch, k0, k1;

    if(ko < 0 || ko >= nok || ncat <= 0)
      return 0;
    if(!np_disc_ordered_near_upper(ko, lam))
      return 0;

    cl = matrix_categorical_vals[oi][0];
    ch = matrix_categorical_vals[oi][ncat - 1];
    k0 = okf[ko](cl, cl, lam, cl, ch);
    k1 = okf[ko](cl, ch, lam, cl, ch);
    if(!np_disc_near_const_kernel(k0, k1))
      return 0;
  }

  return 1;
}

static int np_reg_cv_all_large_gate(const int BANDWIDTH_reg,
                                    const int num_obs,
                                    const int num_reg_continuous,
                                    const int num_reg_unordered,
                                    const int num_reg_ordered,
                                    const int * const kernel_c,
                                    const int * const kernel_u,
                                    const int * const kernel_o,
                                    double **matrix_X_continuous,
                                    double **matrix_bandwidth,
                                    double *lambda,
                                    int *num_categories,
                                    double **matrix_categorical_vals,
                                    int **ov_cont_ok,
                                    double **ov_cont_hmin,
                                    double **ov_cont_k0,
                                    int *ov_cont_from_cache){
  int i;
  int all_large_gate = (BANDWIDTH_reg == BW_FIXED);
  if((ov_cont_ok == NULL) || (ov_cont_hmin == NULL) || (ov_cont_k0 == NULL) ||
     (ov_cont_from_cache == NULL))
    return 0;

  if(!all_large_gate)
    return 0;

  *ov_cont_ok = NULL;
  *ov_cont_hmin = NULL;
  *ov_cont_k0 = NULL;
  *ov_cont_from_cache = 0;

  if(num_reg_continuous > 0){
    const double rel_tol = np_cont_largeh_rel_tol();
    if(np_cont_largeh_cache_get_or_build(num_obs,
                                         num_obs,
                                         num_reg_continuous,
                                         kernel_c,
                                         matrix_X_continuous,
                                         matrix_X_continuous,
                                         rel_tol,
                                         ov_cont_ok,
                                         ov_cont_hmin,
                                         ov_cont_k0)){
      *ov_cont_from_cache = 1;
    } else if(!np_cont_largeh_build_params(num_obs,
                                           num_obs,
                                           num_reg_continuous,
                                           kernel_c,
                                           matrix_X_continuous,
                                           matrix_X_continuous,
                                           rel_tol,
                                           ov_cont_ok,
                                           ov_cont_hmin,
                                           ov_cont_k0)){
      return 0;
    }

    for(i = 0; i < num_reg_continuous; i++){
      const double h = matrix_bandwidth[i][0];
      if((*ov_cont_ok == NULL) || (!(*ov_cont_ok)[i]) || (!isfinite(h)) ||
         (fabs(h) < (*ov_cont_hmin)[i]))
        return 0;
    }
  }

  if(!np_fastcv_disc_unordered_all_large(num_reg_unordered,
                                         kernel_u,
                                         lambda,
                                         num_categories))
    return 0;

  if(!np_fastcv_disc_ordered_all_large(num_reg_unordered,
                                       num_reg_ordered,
                                       kernel_o,
                                       lambda,
                                       num_categories,
                                       matrix_categorical_vals))
    return 0;

  return 1;
}

/* Canonical selector for CVLS route between full symmetric drop-one and reduced branch. */
static inline int np_reg_cv_use_symmetric_dropone_path(const int bwm,
                                                       const int ks_tree_use,
                                                       const int BANDWIDTH_reg){
  return (bwm == RBWM_CVLS) || ks_tree_use || (BANDWIDTH_reg == BW_ADAP_NN);
}

/* Canonical selector for density CV tree bypass in all-large/adaptive regimes. */
static inline int np_den_cv_use_tree_bypass_path(const int gate_x_all_large_fixed,
                                                 const int int_TREE_XY,
                                                 const int BANDWIDTH_den){
  return gate_x_all_large_fixed || !int_TREE_XY || (BANDWIDTH_den == BW_ADAP_NN);
}

typedef struct {
  double cv;
  double traceH;
  int ok;
} NPRegCvLpResult;

static inline int np_reg_cv_use_canonical_glp_fixed_kernel(const int int_ll,
                                                           const int bwm,
                                                           const int BANDWIDTH_reg,
                                                           const int num_reg_continuous,
                                                           const int ks_tree_use,
                                                           const int use_bernstein){
  if((bwm != RBWM_CVLS) ||
     (BANDWIDTH_reg != BW_FIXED) ||
     ks_tree_use ||
     (num_reg_continuous <= 0))
    return 0;

  if(int_ll == LL_LL)
    return 1;

  if((int_ll == LL_LP) &&
     (!use_bernstein) &&
     (int_glp_basis_extern == 1))
    return 1;

  return 0;
}

static inline int np_reg_cv_use_canonical_ll_degree1_lp_objective(const int int_ll,
                                                                  const int bwm,
                                                                  const int BANDWIDTH_reg,
                                                                  const int num_reg_continuous,
                                                                  const int ks_tree_use){
  return (int_ll == LL_LL) &&
    (bwm == RBWM_CVLS) &&
    (BANDWIDTH_reg == BW_FIXED) &&
    (num_reg_continuous > 0) &&
    (!ks_tree_use);
}

static inline int np_reg_use_canonical_glp_degree1_estimation(const int int_ll,
                                                              const int BANDWIDTH_reg,
                                                              const int num_reg_continuous){
  int i;

  if((int_ll != LL_LP) ||
     (BANDWIDTH_reg != BW_GEN_NN) ||
     (num_reg_continuous <= 0) ||
     (vector_glp_degree_extern == NULL) ||
     (int_glp_basis_extern != 1) ||
     (int_glp_bernstein_extern != 0))
    return 0;

  for(i = 0; i < num_reg_continuous; i++){
    if(vector_glp_degree_extern[i] != 1)
      return 0;
  }

  return 1;
}

static inline double np_glp_binom_coeff(const int n, const int k){
  int i, kk;
  double c = 1.0;
  if((k < 0) || (k > n))
    return 0.0;
  kk = MIN(k, n-k);
  for(i = 1; i <= kk; i++)
    c *= ((double)(n - kk + i))/((double)i);
  return c;
}

static void np_glp_fill_shift_raw_from_center(const int ncon,
                                              const int *terms,
                                              const int nterms,
                                              const double *xj,
                                              MATRIX SHIFT){
  int r, c, j;
  for(r = 0; r < nterms; r++){
    const int *alpha = terms + r*ncon;
    for(c = 0; c < nterms; c++){
      const int *beta = terms + c*ncon;
      double coef = 1.0;
      int ok = 1;
      for(j = 0; j < ncon; j++){
        const int aj = alpha[j];
        const int bj = beta[j];
        if(bj > aj){
          ok = 0;
          break;
        }
        coef *= np_glp_binom_coeff(aj, bj)*ipow(xj[j], aj - bj);
      }
      SHIFT[r][c] = ok ? coef : 0.0;
    }
  }
}

static int np_glp_center_raw_moments_at_eval(const int ncon,
                                             const int *terms,
                                             const int nterms,
                                             const double *xj,
                                             const double *raw_s,
                                             const double *raw_t,
                                             MATRIX SHIFT,
                                             MATRIX SHIFTINV,
                                             MATRIX TMP,
                                             MATRIX KWM,
                                             MATRIX XTKY){
  int a, b, c;

  np_glp_fill_shift_raw_from_center(ncon, terms, nterms, xj, SHIFT);
  if(mat_inv(SHIFT, SHIFTINV) == NULL)
    return 0;

  for(a = 0; a < nterms; a++){
    double ta = 0.0;
    for(c = 0; c < nterms; c++)
      ta += SHIFTINV[a][c]*raw_t[c];
    XTKY[a][0] = ta;
  }

  for(a = 0; a < nterms; a++){
    for(b = 0; b < nterms; b++){
      double sab = 0.0;
      for(c = 0; c < nterms; c++)
        sab += SHIFTINV[a][c]*raw_s[c*nterms+b];
      TMP[a][b] = sab;
    }
  }

  for(a = 0; a < nterms; a++){
    for(b = 0; b < nterms; b++){
      double sab = 0.0;
      for(c = 0; c < nterms; c++)
        sab += TMP[a][c]*SHIFTINV[b][c];
      KWM[a][b] = sab;
    }
  }

  return 1;
}

static NPRegCvLpResult np_regression_cv_glp_rawbasis_fixed(
    const int int_ll,
    const int num_obs,
    const int num_reg_unordered,
    const int num_reg_ordered,
    const int num_reg_continuous,
    double **matrix_X_unordered,
    double **matrix_X_ordered,
    double **matrix_X_continuous,
    double *vector_Y,
    double *vector_scale_factor,
    int *num_categories,
    int *kernel_c,
    int *kernel_u,
    int *kernel_o,
    int *operator,
    double *lambda,
    double **matrix_bandwidth,
    const int *glp_terms_in,
    const int glp_nterms_in,
    double **glp_basis_in){
  NPRegCvLpResult result = {DBL_MAX, 0.0, 0};
  const double epsilon = 1.0/(double)MAX(1, num_obs);
  int i, j, a, b, l, sf_flag = 0;
  int local_fail = 0;
  int nterms = glp_nterms_in;
  const int *terms = glp_terms_in;
  double **basis = glp_basis_in;
  int *terms_local = NULL;
  double **basis_local = NULL;
  double *ones = NULL;
  double *moments = NULL, *rhs = NULL, *kw = NULL, *xj = NULL;
  double *moments_local = NULL, *rhs_local = NULL;
  double *vsf = NULL;
  double **train_u = NULL, **train_o = NULL, **train_c = NULL;
  MATRIX eval_u = NULL, eval_o = NULL, eval_c = NULL;
  MATRIX matrix_bandwidth_eval = NULL;
  MATRIX KWM = NULL, XTKY = NULL, DELTA = NULL, SHIFT = NULL, SHIFTINV = NULL, TMP = NULL;
  double mean_dummy = 0.0;
#ifdef MPI2
  const int use_mpi_transport = (iNum_Processors > 1);
#else
  const int use_mpi_transport = 0;
#endif

#define NP_GLP_CV_FAIL() do { local_fail = 1; goto glp_cv_collective_gate; } while(0)

  if((num_obs <= 0) || (num_reg_continuous <= 0))
    return result;

  if(int_ll == LL_LL){
    nterms = num_reg_continuous + 1;
    terms_local = (int *)calloc((size_t)nterms*(size_t)num_reg_continuous, sizeof(int));
    basis_local = (double **)malloc((size_t)nterms*sizeof(double *));
    ones = alloc_vecd(num_obs);
    if((terms_local == NULL) || (basis_local == NULL) || (ones == NULL))
      NP_GLP_CV_FAIL();
    for(i = 0; i < num_obs; i++)
      ones[i] = 1.0;
    basis_local[0] = ones;
    for(l = 0; l < num_reg_continuous; l++){
      terms_local[(l+1)*num_reg_continuous + l] = 1;
      basis_local[l+1] = matrix_X_continuous[l];
    }
    terms = terms_local;
    basis = basis_local;
  }

  if((terms == NULL) || (basis == NULL) || (nterms <= 0))
    NP_GLP_CV_FAIL();

  if((sf_flag = (int_LARGE_SF == 0))){
    int_LARGE_SF = 1;
    vsf = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
    if(vsf == NULL)
      NP_GLP_CV_FAIL();
    for(l = 0; l < num_reg_continuous; l++)
      vsf[l] = matrix_bandwidth[l][0];
  } else {
    vsf = vector_scale_factor;
  }

  moments = (double *)calloc((size_t)num_obs, (size_t)nterms*(size_t)nterms*sizeof(double));
  rhs = (double *)calloc((size_t)num_obs, (size_t)nterms*sizeof(double));
  kw = alloc_vecd(MAX(1, num_obs));
  xj = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
  train_u = (double **)malloc((size_t)MAX(1, num_reg_unordered)*sizeof(double *));
  train_o = (double **)malloc((size_t)MAX(1, num_reg_ordered)*sizeof(double *));
  train_c = (double **)malloc((size_t)MAX(1, num_reg_continuous)*sizeof(double *));
  eval_u = mat_creat(num_reg_unordered, 1, UNDEFINED);
  eval_o = mat_creat(num_reg_ordered, 1, UNDEFINED);
  eval_c = mat_creat(num_reg_continuous, 1, UNDEFINED);
  matrix_bandwidth_eval = alloc_tmatd(1, num_reg_continuous);
  KWM = mat_creat(nterms, nterms, UNDEFINED);
  XTKY = mat_creat(nterms, 1, UNDEFINED);
  DELTA = mat_creat(nterms, 1, UNDEFINED);
  SHIFT = mat_creat(nterms, nterms, UNDEFINED);
  SHIFTINV = mat_creat(nterms, nterms, UNDEFINED);
  TMP = mat_creat(nterms, nterms, UNDEFINED);

  if((moments == NULL) || (rhs == NULL) || (kw == NULL) || (xj == NULL) ||
     (train_u == NULL) || (train_o == NULL) || (train_c == NULL) ||
     (eval_u == NULL) || (eval_o == NULL) || (eval_c == NULL) ||
     (matrix_bandwidth_eval == NULL) || (KWM == NULL) || (XTKY == NULL) ||
     (DELTA == NULL) || (SHIFT == NULL) || (SHIFTINV == NULL) || (TMP == NULL))
    NP_GLP_CV_FAIL();

#ifdef MPI2
  if(use_mpi_transport){
    moments_local = (double *)calloc((size_t)num_obs, (size_t)nterms*(size_t)nterms*sizeof(double));
    rhs_local = (double *)calloc((size_t)num_obs, (size_t)nterms*sizeof(double));
    if((moments_local == NULL) || (rhs_local == NULL))
      NP_GLP_CV_FAIL();
  }
#endif

  for(j = 0; j < (use_mpi_transport ? num_obs : (num_obs - 1)); j++){
    const double yj = vector_Y[j];
    const int nsub = num_obs - j - 1;
#ifdef MPI2
    double * const moments_acc = use_mpi_transport ? moments_local : moments;
    double * const rhs_acc = use_mpi_transport ? rhs_local : rhs;

    if(use_mpi_transport && ((j % iNum_Processors) != my_rank))
      continue;
#else
    double * const moments_acc = moments;
    double * const rhs_acc = rhs;
#endif

    for(l = 0; l < num_reg_unordered; l++){
      eval_u[l][0] = matrix_X_unordered[l][j];
      train_u[l] = use_mpi_transport ? matrix_X_unordered[l] : (matrix_X_unordered[l] + j + 1);
    }
    for(l = 0; l < num_reg_ordered; l++){
      eval_o[l][0] = matrix_X_ordered[l][j];
      train_o[l] = use_mpi_transport ? matrix_X_ordered[l] : (matrix_X_ordered[l] + j + 1);
    }
    for(l = 0; l < num_reg_continuous; l++){
      eval_c[l][0] = matrix_X_continuous[l][j];
      matrix_bandwidth_eval[l][0] = matrix_bandwidth[l][0];
      train_c[l] = use_mpi_transport ? matrix_X_continuous[l] : (matrix_X_continuous[l] + j + 1);
    }

    if(kernel_weighted_sum_np_ctx(kernel_c,
                                  kernel_u,
                                  kernel_o,
                                  BW_FIXED,
                                  use_mpi_transport ? num_obs : nsub,
                                  1,
                                  num_reg_unordered,
                                  num_reg_ordered,
                                  num_reg_continuous,
                                  0,
                                  0,
                                  1,
                                  1,
                                  1,
                                  0,
                                  0,
                                  use_mpi_transport ? 1 : 0,
                                  j,
                                  operator,
                                  OP_NOOP,
                                  0,
                                  0,
                                  NULL,
                                  1,
                                  0,
                                  0,
                                  NP_TREE_FALSE,
                                  0,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  train_u,
                                  train_o,
                                  train_c,
                                  eval_u,
                                  eval_o,
                                  eval_c,
                                  NULL,
                                  NULL,
                                  NULL,
                                  vsf,
                                  1,
                                  matrix_bandwidth,
                                  matrix_bandwidth_eval,
                                  lambda,
                                  num_categories,
                                  matrix_categorical_vals_extern,
                                  NULL,
                                  &mean_dummy,
                                  NULL,
                                  kw,
                                  NULL) != 0)
      NP_GLP_CV_FAIL();

    if(use_mpi_transport){
      double * const sj = moments_acc + (size_t)j*(size_t)nterms*(size_t)nterms;
      double * const tj = rhs_acc + (size_t)j*(size_t)nterms;

      for(i = 0; i < num_obs; i++){
        const double w = kw[i];
        const double yi = vector_Y[i];

        if((i == j) || (w == 0.0))
          continue;

        for(a = 0; a < nterms; a++){
          const double bia = basis[a][i];
          tj[a] += w*bia*yi;
          for(b = 0; b < nterms; b++)
            sj[a*nterms+b] += w*bia*basis[b][i];
        }
      }
    } else {
      for(i = 0; i < nsub; i++){
        const int ii = j + 1 + i;
        const double w = kw[i];
        const double yi = vector_Y[ii];
        double * const sj = moments_acc + (size_t)j*(size_t)nterms*(size_t)nterms;
        double * const si = moments_acc + (size_t)ii*(size_t)nterms*(size_t)nterms;
        double * const tj = rhs_acc + (size_t)j*(size_t)nterms;
        double * const ti = rhs_acc + (size_t)ii*(size_t)nterms;

        if(w == 0.0)
          continue;

        for(a = 0; a < nterms; a++){
          const double bia = basis[a][ii];
          const double bja = basis[a][j];
          tj[a] += w*bia*yi;
          ti[a] += w*bja*yj;
          for(b = 0; b < nterms; b++){
            sj[a*nterms+b] += w*bia*basis[b][ii];
            si[a*nterms+b] += w*bja*basis[b][j];
          }
        }
      }
    }
  }

glp_cv_collective_gate:
#ifdef MPI2
  if(use_mpi_transport && np_mpi_rank_failure_injected("NP_RMPI_INJECT_GLP_CV_FAIL_RANK"))
    local_fail = 1;
#endif

#ifdef MPI2
  if(use_mpi_transport){
    int any_fail = 0;

    MPI_Allreduce(&local_fail, &any_fail, 1, MPI_INT, MPI_MAX, comm[1]);
    if(any_fail)
      goto cleanup_glp_cv;

    MPI_Allreduce(moments_local, moments, num_obs*nterms*nterms, MPI_DOUBLE, MPI_SUM, comm[1]);
    MPI_Allreduce(rhs_local, rhs, num_obs*nterms, MPI_DOUBLE, MPI_SUM, comm[1]);
  }
#endif

  if(local_fail)
    goto cleanup_glp_cv;

  result.cv = 0.0;
  result.traceH = 0.0;

  for(j = 0; j < num_obs; j++){
    const double * const sj = moments + (size_t)j*(size_t)nterms*(size_t)nterms;
    const double * const tj = rhs + (size_t)j*(size_t)nterms;
    double nepsilon = 0.0;
    double fit = 0.0;

    for(l = 0; l < num_reg_continuous; l++)
      xj[l] = matrix_X_continuous[l][j];

    if(!np_glp_center_raw_moments_at_eval(num_reg_continuous,
                                          terms,
                                          nterms,
                                          xj,
                                          sj,
                                          tj,
                                          SHIFT,
                                          SHIFTINV,
                                          TMP,
                                          KWM,
                                          XTKY))
      goto cleanup_glp_cv;

    while(mat_solve(KWM, XTKY, DELTA) == NULL){
      for(a = 0; a < nterms; a++)
        KWM[a][a] += epsilon;
      nepsilon += epsilon;
      if(nepsilon > 128.0*epsilon)
        goto cleanup_glp_cv;
    }

    XTKY[0][0] += nepsilon*XTKY[0][0]/NZD_POS(KWM[0][0]);
    if(nepsilon > 0.0){
      if(mat_solve(KWM, XTKY, DELTA) == NULL)
        goto cleanup_glp_cv;
    }

    fit = DELTA[0][0];

    {
      const double dy = vector_Y[j] - fit;
      result.cv += dy*dy;
    }
  }

  result.ok = 1;

cleanup_glp_cv:
  if(train_u != NULL) free(train_u);
  if(train_o != NULL) free(train_o);
  if(train_c != NULL) free(train_c);
  if(eval_u != NULL) mat_free(eval_u);
  if(eval_o != NULL) mat_free(eval_o);
  if(eval_c != NULL) mat_free(eval_c);
  if(matrix_bandwidth_eval != NULL) free_tmat(matrix_bandwidth_eval);
  if(KWM != NULL) mat_free(KWM);
  if(XTKY != NULL) mat_free(XTKY);
  if(DELTA != NULL) mat_free(DELTA);
  if(SHIFT != NULL) mat_free(SHIFT);
  if(SHIFTINV != NULL) mat_free(SHIFTINV);
  if(TMP != NULL) mat_free(TMP);
  if(moments != NULL) free(moments);
  if(moments_local != NULL) free(moments_local);
  if(rhs != NULL) free(rhs);
  if(rhs_local != NULL) free(rhs_local);
  if(xj != NULL) free(xj);
  if(kw != NULL) free(kw);
  if(terms_local != NULL) free(terms_local);
  if(basis_local != NULL) free(basis_local);
  if(ones != NULL) free(ones);
  if(sf_flag){
    int_LARGE_SF = 0;
    free(vsf);
  }

  if(!result.ok){
    result.cv = DBL_MAX;
    result.traceH = 0.0;
  }

#undef NP_GLP_CV_FAIL

  return result;
}

// Regression CV objective for local polynomial regression:
// lc (degree 0), ll (degree 1), and lp (general degree vector).
// The LL/LP branches solve weighted normal equations with ridge fallback if singular.

double np_kernel_estimate_regression_categorical_ls_aic(
int int_ll,
int bwm,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_reg,
int num_obs,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_X_unordered,
double **matrix_X_ordered,
double **matrix_X_continuous,
double *vector_Y,
double *vector_scale_factor,
int *num_categories){

  np_gate_override_clear();

  // note that mean has 2*num_obs allocated for npksum
  int i, j, l, sf_flag = 0, num_obs_eval_alloc, tsf;
  double cv = 0.0;
  double * lambda = NULL, * vsf = NULL;
  double ** matrix_bandwidth = NULL;

  double aicc = 0.0;
  double traceH = 0.0;
int * operator = NULL;
int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;
  int *ov_cont_ok = NULL;
  double *ov_cont_hmin = NULL, *ov_cont_k0 = NULL;
  int ov_cont_from_cache = 0;

  const int leave_one_out = (bwm == RBWM_CVLS)?1:0;

  if(!np_reg_cv_core_cache_prepare(KERNEL_reg,
                                   KERNEL_unordered_reg,
                                   KERNEL_ordered_reg,
                                   BANDWIDTH_reg,
                                   num_obs,
                                   num_reg_continuous,
                                   num_reg_unordered,
                                   num_reg_ordered))
    return DBL_MAX;
  /* Canonical CVAIC parity: LP(degree=1) should follow the LL-equivalent
     objective branch for bandwidth search. */
  int int_ll_cv = int_ll;
  if((int_ll == LL_LP) &&
     (bwm == RBWM_CVAIC) &&
     (vector_glp_degree_extern != NULL) &&
     (num_reg_continuous > 0)){
    int all_deg_one = 1;
    for(i = 0; i < num_reg_continuous; i++){
      if(vector_glp_degree_extern[i] != 1){
        all_deg_one = 0;
        break;
      }
    }
    if(all_deg_one)
      int_ll_cv = LL_LL;
  }

  operator = np_reg_cv_core_cache.operator;
  kernel_c = np_reg_cv_core_cache.kernel_c;
  kernel_u = np_reg_cv_core_cache.kernel_u;
  kernel_o = np_reg_cv_core_cache.kernel_o;
  lambda = np_reg_cv_core_cache.lambda;
  matrix_bandwidth = np_reg_cv_core_cache.matrix_bandwidth;

#ifdef MPI2
    int stride = MAX((int)ceil((double) num_obs / (double) iNum_Processors),1);
    num_obs_eval_alloc = stride*iNum_Processors;
#else
    num_obs_eval_alloc = num_obs;
#endif

    int ks_tree_use = (int_TREE_X == NP_TREE_TRUE) && (!((BANDWIDTH_reg == BW_ADAP_NN) && (int_ll_cv == LL_LL)));
    if(np_reg_cv_use_canonical_ll_degree1_lp_objective(int_ll,
                                                       bwm,
                                                       BANDWIDTH_reg,
                                                       num_reg_continuous,
                                                       ks_tree_use))
      int_ll_cv = LL_LP;

  if(kernel_bandwidth_mean(KERNEL_reg,
                           BANDWIDTH_reg,
                           num_obs,
                           num_obs,
                           0,
                           0,
                           0,
                           num_reg_continuous,
                           num_reg_unordered,
                           num_reg_ordered,
                           0, // do not suppress_parallel
                           vector_scale_factor,
                           NULL,			 // Not used 
                           NULL,			 // Not used 
                           matrix_X_continuous,
                           matrix_X_continuous,
                           NULL,					 // Not used 
                           matrix_bandwidth,
                           lambda)==1){
    
    return(DBL_MAX);
  }
  if(int_ll_cv == LL_LP){
    const int use_bernstein = (int_glp_bernstein_extern != 0);
    const int *glp_terms = NULL;
    int glp_nterms = 0;
    double **basis = NULL;
    MATRIX XTKY = NULL, DELTA = NULL, KWM = NULL;
    MATRIX TCON = NULL, TUNO = NULL, TORD = NULL;
    double ** matrix_bandwidth_eval = NULL;
    double ** XTKX = NULL;
    int glp_ok = 1;

    if((vector_glp_degree_extern == NULL) || (num_reg_continuous <= 0)){
      cv = DBL_MAX;
      goto finish_cv_path;
    }

    if(!np_glp_cv_cache.ready ||
       (np_glp_cv_cache.use_bernstein != use_bernstein) ||
       (np_glp_cv_cache.basis_mode != int_glp_basis_extern) ||
       (np_glp_cv_cache.num_obs != num_obs) ||
       (np_glp_cv_cache.ncon != num_reg_continuous) ||
       (np_glp_cv_cache.matrix_X_continuous_train_ptr != matrix_X_continuous)){
      if(!np_glp_cv_cache_prepare(int_ll_cv, num_obs, num_reg_continuous, matrix_X_continuous)){
        cv = DBL_MAX;
        goto finish_cv_path;
      }
    }
    if(!np_glp_cv_cache.ready){
      cv = DBL_MAX;
      goto finish_cv_path;
    }

    glp_terms = np_glp_cv_cache.terms;
    glp_nterms = np_glp_cv_cache.nterms;
    basis = np_glp_cv_cache.basis;
    if((glp_terms == NULL) || (basis == NULL) || (glp_nterms <= 0)){
      cv = DBL_MAX;
      goto finish_cv_path;
    }

    {
      const int use_canonical_glp_kernel =
        np_reg_cv_use_canonical_glp_fixed_kernel(LL_LP,
                                                 bwm,
                                                 BANDWIDTH_reg,
                                                 num_reg_continuous,
                                                 ks_tree_use,
                                                 use_bernstein);
      const int all_large_gate = np_reg_cv_all_large_gate(BANDWIDTH_reg,
                                                          num_obs,
                                                          num_reg_continuous,
                                                          num_reg_unordered,
                                                          num_reg_ordered,
                                                          kernel_c,
                                                          kernel_u,
                                                          kernel_o,
                                                          matrix_X_continuous,
                                                          matrix_bandwidth,
                                                          lambda,
                                                          num_categories,
                                                          matrix_categorical_vals_extern,
                                                          &ov_cont_ok,
                                                          &ov_cont_hmin,
                                                          &ov_cont_k0,
                                                          &ov_cont_from_cache);
      if(use_canonical_glp_kernel && !all_large_gate){
        NPRegCvLpResult glp_result =
          np_regression_cv_glp_rawbasis_fixed(LL_LP,
                                              num_obs,
                                              num_reg_unordered,
                                              num_reg_ordered,
                                              num_reg_continuous,
                                              matrix_X_unordered,
                                              matrix_X_ordered,
                                              matrix_X_continuous,
                                              vector_Y,
                                              vector_scale_factor,
                                              num_categories,
                                              kernel_c,
                                              kernel_u,
                                              kernel_o,
                                              operator,
                                              lambda,
                                              matrix_bandwidth,
                                              glp_terms,
                                              glp_nterms,
                                              basis);
        cv = glp_result.cv;
        traceH = glp_result.traceH;
        goto finish_cv_path;
      }

      if(all_large_gate){
        const int k = glp_nterms;
        MATRIX XtX = mat_creat(k, k, UNDEFINED);
        MATRIX XtXINV = mat_creat(k, k, UNDEFINED);
        MATRIX XtY = mat_creat(k, 1, UNDEFINED);
        MATRIX BETA = mat_creat(k, 1, UNDEFINED);
        int fast_ok = (XtX != NULL) && (XtXINV != NULL) && (XtY != NULL) && (BETA != NULL);

        if(fast_ok){
          const double ridge_eps = 1.0/(double)MAX(1, num_obs);
          int ridge_it = 0;

          for(i = 0; i < k; i++){
            XtY[i][0] = 0.0;
            BETA[i][0] = 0.0;
            for(j = 0; j < k; j++)
              XtX[i][j] = 0.0;
          }

          for(i = 0; i < num_obs; i++){
            const double yi = vector_Y[i];
            for(int a = 0; a < k; a++){
              const double za = basis[a][i];
              XtY[a][0] += za*yi;
              for(int b = a; b < k; b++){
                const double zb = basis[b][i];
                XtX[a][b] += za*zb;
                if(b != a) XtX[b][a] += za*zb;
              }
            }
          }

          while(mat_inv(XtX, XtXINV) == NULL){
            for(i = 0; i < k; i++) XtX[i][i] += ridge_eps;
            ridge_it++;
            if(ridge_it > 64){
              fast_ok = 0;
              break;
            }
          }

          if(fast_ok){
            for(i = 0; i < k; i++){
              double s = 0.0;
              for(j = 0; j < k; j++) s += XtXINV[i][j]*XtY[j][0];
              BETA[i][0] = s;
            }

            cv = 0.0;
            traceH = 0.0;
            for(i = 0; i < num_obs; i++){
              double yhat = 0.0;
              double hii = 0.0;
              for(j = 0; j < k; j++){
                const double zj = basis[j][i];
                yhat += BETA[j][0]*zj;
              }
              for(j = 0; j < k; j++){
                const double zj = basis[j][i];
                for(int b = 0; b < k; b++)
                  hii += zj*XtXINV[j][b]*basis[b][i];
              }
              {
                const double err = vector_Y[i] - yhat;
                if(bwm == RBWM_CVLS){
                  const double den = NZD_POS(1.0 - hii);
                  const double err_loo = err/den;
                  cv += err_loo*err_loo;
                } else {
                  cv += err*err;
                  traceH += hii;
                }
              }
            }
          }
        }

        if(XtX != NULL) mat_free(XtX);
        if(XtXINV != NULL) mat_free(XtXINV);
        if(XtY != NULL) mat_free(XtY);
        if(BETA != NULL) mat_free(BETA);

        if(fast_ok){
          np_fastcv_alllarge_hits++;
          goto finish_cv_path;
        }
      }
    }

    const int nrc1 = glp_nterms;
    const int nrc2 = nrc1 + 1;
    const int nrcc22 = nrc2*nrc2;
    double * PKWM[nrc1], * PXTKY[nrc1], * PXTKX[nrc2];

    double * PXC[MAX(1,num_reg_continuous)];
    double * PXU[MAX(1,num_reg_unordered)];
    double * PXO[MAX(1,num_reg_ordered)];

    PXC[0] = NULL;
    PXU[0] = NULL;
    PXO[0] = NULL;

    for(l = 0; l < num_reg_continuous; l++)
      PXC[l] = matrix_X_continuous[l];
    for(l = 0; l < num_reg_unordered; l++)
      PXU[l] = matrix_X_unordered[l];
    for(l = 0; l < num_reg_ordered; l++)
      PXO[l] = matrix_X_ordered[l];

    if((sf_flag = (int_LARGE_SF == 0))){
      int_LARGE_SF = 1;
      vsf = (double *)malloc(num_reg_continuous*sizeof(double));
      for(i = 0; i < num_reg_continuous; i++)
        vsf[i] = matrix_bandwidth[i][0];
    } else {
      vsf = vector_scale_factor;
    }

    XTKY = mat_creat(nrc1, 1, UNDEFINED);
    DELTA = mat_creat(nrc1, 1, UNDEFINED);
    KWM = mat_creat(nrc1, nrc1, UNDEFINED);
    TCON = mat_creat(num_reg_continuous, 1, UNDEFINED);
    TUNO = mat_creat(num_reg_unordered, 1, UNDEFINED);
    TORD = mat_creat(num_reg_ordered, 1, UNDEFINED);
    matrix_bandwidth_eval = alloc_tmatd(1, num_reg_continuous);
    XTKX = (double **)malloc((size_t)nrc2*sizeof(double *));

    const size_t kwm_len = (size_t)nrcc22*(size_t)num_obs_eval_alloc;
    double * kwm = (double *)malloc(kwm_len*sizeof(double));
    double * sgn = (double *)malloc((size_t)nrc2*sizeof(double));
    double * evalv = (double *)malloc((size_t)nrc1*sizeof(double));

    glp_ok = (XTKY != NULL) && (DELTA != NULL) && (KWM != NULL) &&
      (TCON != NULL) && (TUNO != NULL) && (TORD != NULL) &&
      (matrix_bandwidth_eval != NULL) && (XTKX != NULL) &&
      (kwm != NULL) && (sgn != NULL) && (evalv != NULL);

    if(!glp_ok){
      cv = DBL_MAX;
    } else {
      for(size_t ii = 0; ii < kwm_len; ii++)
        kwm[ii] = 0.0;

      sgn[0] = 1.0;
      for(i = 1; i < nrc2; i++) sgn[i] = 1.0;

      XTKX[0] = vector_Y;
      for(i = 0; i < nrc1; i++)
        XTKX[i+1] = basis[i];

      for(i = 0; i < nrc2; i++)
        PXTKX[i] = XTKX[i];

      for(i = 0; i < nrc1; i++){
        PKWM[i] = KWM[i];
        PXTKY[i] = XTKY[i];
        KWM[i] = &kwm[(i+1)*nrc2 + 1];
        XTKY[i] = &kwm[i+1];
      }

      if(bwm == RBWM_CVAIC){
        tsf = int_LARGE_SF;
        int_LARGE_SF = 1;
        kernel_weighted_sum_np_ctx(kernel_c,
                                   kernel_u,
                                   kernel_o,
                                   BANDWIDTH_reg,
                                   1,
                                   1,
                                   num_reg_unordered,
                                   num_reg_ordered,
                                   num_reg_continuous,
                                   0,
                                   0,
                                   1,
                                   0,
                                   0,
                                   0,
                                   0,
                                   0,
                                   0,
                                   operator,
                                   OP_NOOP,
                                   0,
                                   0,
                                   NULL,
                                   1,
                                   0,
                                   0,
                                   NP_TREE_FALSE,
                                   0,
                                   NULL, NULL, NULL, NULL,
                                   matrix_X_unordered,
                                   matrix_X_ordered,
                                   matrix_X_continuous,
                                   matrix_X_unordered,
                                   matrix_X_ordered,
                                   matrix_X_continuous,
                                   NULL,
                                   NULL,
                                   NULL,
                                   vector_scale_factor,
                                   1,
                                   matrix_bandwidth,
                                   matrix_bandwidth,
                                   lambda,
                                   num_categories,
                                   NULL,
                                   NULL,
                                   &aicc,
                                   NULL,
                                   NULL,
                                   NULL);
        int_LARGE_SF = tsf;
      }

      const double epsilon = 1.0/(double)MAX(1, num_obs);
      for(j = 0; j < num_obs; j++){
        double nepsilon = 0.0;
        double pnh = 1.0;

        for(i = 0; i < nrc1; i++){
          KWM[i] = &kwm[j*nrcc22+(i+1)*nrc2+1];
          XTKY[i] = &kwm[j*nrcc22+i+1];
        }

#ifdef MPI2
        if(np_reg_cv_use_symmetric_dropone_path(bwm, ks_tree_use, BANDWIDTH_reg)){
          if((j % iNum_Processors) == 0){
            if((j+my_rank) < (num_obs)){
              for(l = 0; l < num_reg_continuous; l++){
                TCON[l][0] = matrix_X_continuous[l][j+my_rank];
                if(BANDWIDTH_reg == BW_GEN_NN)
                  matrix_bandwidth_eval[l][0] = matrix_bandwidth[l][j+my_rank];
              }
              for(l = 0; l < num_reg_unordered; l++)
                TUNO[l][0] = matrix_X_unordered[l][j+my_rank];
              for(l = 0; l < num_reg_ordered; l++)
                TORD[l][0] = matrix_X_ordered[l][j+my_rank];

              kernel_weighted_sum_np_ctx(kernel_c,
                                         kernel_u,
                                         kernel_o,
                                         BANDWIDTH_reg,
                                         num_obs,
                                         1,
                                         num_reg_unordered,
                                         num_reg_ordered,
                                         num_reg_continuous,
                                         0,
                                         0,
                                         1,
                                         (BANDWIDTH_reg == BW_ADAP_NN)?1:0,
                                         0,
                                         1,
                                         0,
                                         1,
                                         j+my_rank,
                                         operator,
                                         OP_NOOP,
                                         0,
                                         0,
                                         NULL,
                                         1,
                                         nrc2,
                                         nrc2,
                                         (BANDWIDTH_reg == BW_ADAP_NN) ? NP_TREE_FALSE : int_TREE_X,
                                         0,
                                         (BANDWIDTH_reg == BW_ADAP_NN) ? NULL : kdt_extern_X,
                                         NULL, NULL, NULL,
                                         PXU,
                                         PXO,
                                         PXC,
                                         TUNO,
                                         TORD,
                                         TCON,
                                         XTKX,
                                         XTKX,
                                         NULL,
                                         vsf,
                                         1,
                                         matrix_bandwidth,
                                         matrix_bandwidth_eval,
                                         lambda,
                                         num_categories,
                                         NULL,
                                         NULL,
                                         kwm+(j+my_rank)*nrcc22,
                                         NULL,
                                         NULL,
                                         NULL);
            }
            MPI_Allgather(MPI_IN_PLACE, nrcc22, MPI_DOUBLE, kwm+j*nrcc22, nrcc22, MPI_DOUBLE, comm[1]);
          }
        } else {
          if((j % iNum_Processors) == 0){
            if((j+my_rank) < (num_obs-1)){
              for(l = 0; l < nrc2; l++)
                XTKX[l] = PXTKX[l] + j + my_rank + 1;
              for(l = 0; l < num_reg_continuous; l++)
                PXC[l] = matrix_X_continuous[l] + j + my_rank + 1;
              for(l = 0; l < num_reg_unordered; l++)
                PXU[l] = matrix_X_unordered[l] + j + my_rank + 1;
              for(l = 0; l < num_reg_ordered; l++)
                PXO[l] = matrix_X_ordered[l] + j + my_rank + 1;

              for(l = 0; l < num_reg_continuous; l++){
                TCON[l][0] = matrix_X_continuous[l][j+my_rank];
                if(BANDWIDTH_reg == BW_GEN_NN)
                  matrix_bandwidth_eval[l][0] = matrix_bandwidth[l][j+my_rank];
              }
              for(l = 0; l < num_reg_unordered; l++)
                TUNO[l][0] = matrix_X_unordered[l][j+my_rank];
              for(l = 0; l < num_reg_ordered; l++)
                TORD[l][0] = matrix_X_ordered[l][j+my_rank];

              kernel_weighted_sum_np_ctx(kernel_c,
                                         kernel_u,
                                         kernel_o,
                                         BANDWIDTH_reg,
                                         num_obs-j-my_rank-1,
                                         1,
                                         num_reg_unordered,
                                         num_reg_ordered,
                                         num_reg_continuous,
                                         0,
                                         0,
                                         1,
                                         (BANDWIDTH_reg == BW_ADAP_NN)?1:0,
                                         0,
                                         1,
                                         1,
                                         0,
                                         0,
                                         operator,
                                         OP_NOOP,
                                         0,
                                         0,
                                         NULL,
                                         1,
                                         nrc2,
                                         nrc2,
                                         0,
                                         0,
                                         NULL, NULL, NULL, NULL,
                                         PXU,
                                         PXO,
                                         PXC,
                                         TUNO,
                                         TORD,
                                         TCON,
                                         XTKX,
                                         XTKX,
                                         sgn,
                                         vsf,
                                         1,
                                         matrix_bandwidth,
                                         matrix_bandwidth_eval,
                                         lambda,
                                         num_categories,
                                         NULL,
                                         NULL,
                                         kwm+(j+my_rank)*nrcc22,
                                         NULL,
                                         NULL,
                                         NULL);

              for(int jj = j+my_rank+1; jj < num_obs; jj++){
                const double RW = kwm[jj*nrcc22+nrc1]*(XTKX[0][-1]-XTKX[0][jj-j-my_rank-1]);
                for(int ii = 1; ii < nrc2; ii++)
                  kwm[jj*nrcc22+ii*nrc2] += RW*XTKX[ii][jj-j-my_rank-1]*sgn[ii];
              }
            }

            const int nrem = num_obs % iNum_Processors;
            const int nred = ((j+iNum_Processors) > num_obs) ? nrem : iNum_Processors;
            MPI_Allreduce(MPI_IN_PLACE, kwm+j*nrcc22, nred*nrcc22, MPI_DOUBLE, MPI_SUM, comm[1]);
          }

          double * const tpk = kwm+j*nrcc22;
          for (int jj = 0; jj < nrc2; jj++){
            for (int ii = nrc1; ii > jj; ii--)
              tpk[jj*nrc2+ii] = tpk[ii*nrc2+jj];
          }
        }
#else
        if(np_reg_cv_use_symmetric_dropone_path(bwm, ks_tree_use, BANDWIDTH_reg)){
          for(l = 0; l < num_reg_continuous; l++){
            TCON[l][0] = matrix_X_continuous[l][j];
            if(BANDWIDTH_reg == BW_GEN_NN)
              matrix_bandwidth_eval[l][0] = matrix_bandwidth[l][j];
          }
          for(l = 0; l < num_reg_unordered; l++)
            TUNO[l][0] = matrix_X_unordered[l][j];
          for(l = 0; l < num_reg_ordered; l++)
            TORD[l][0] = matrix_X_ordered[l][j];

          kernel_weighted_sum_np_ctx(kernel_c,
                                     kernel_u,
                                     kernel_o,
                                     BANDWIDTH_reg,
                                     num_obs,
                                     1,
                                     num_reg_unordered,
                                     num_reg_ordered,
                                     num_reg_continuous,
                                     0,
                                     0,
                                     1,
                                     (BANDWIDTH_reg == BW_ADAP_NN)?1:0,
                                     0,
                                     1,
                                     0,
                                     1,
                                     j,
                                     operator,
                                     OP_NOOP,
                                     0,
                                     0,
                                     NULL,
                                     0,
                                     nrc2,
                                     nrc2,
                                     (BANDWIDTH_reg == BW_ADAP_NN) ? NP_TREE_FALSE : int_TREE_X,
                                     0,
                                     (BANDWIDTH_reg == BW_ADAP_NN) ? NULL : kdt_extern_X,
                                     NULL, NULL, NULL,
                                     PXU,
                                     PXO,
                                     PXC,
                                     TUNO,
                                     TORD,
                                     TCON,
                                     XTKX,
                                     XTKX,
                                     NULL,
                                     vsf,
                                     1,
                                     matrix_bandwidth,
                                     matrix_bandwidth_eval,
                                     lambda,
                                     num_categories,
                                     NULL,
                                     NULL,
                                     kwm+j*nrcc22,
                                     NULL,
                                     NULL,
                                     NULL);
        } else {
          if(j < (num_obs-1)){
            for(l = 0; l < nrc2; l++)
              XTKX[l]++;
            for(l = 0; l < num_reg_continuous; l++)
              PXC[l]++;
            for(l = 0; l < num_reg_unordered; l++)
              PXU[l]++;
            for(l = 0; l < num_reg_ordered; l++)
              PXO[l]++;

            for(l = 0; l < num_reg_continuous; l++){
              TCON[l][0] = matrix_X_continuous[l][j];
              if(BANDWIDTH_reg == BW_GEN_NN)
                matrix_bandwidth_eval[l][0] = matrix_bandwidth[l][j];
            }
            for(l = 0; l < num_reg_unordered; l++)
              TUNO[l][0] = matrix_X_unordered[l][j];
            for(l = 0; l < num_reg_ordered; l++)
              TORD[l][0] = matrix_X_ordered[l][j];

            kernel_weighted_sum_np_ctx(kernel_c,
                                       kernel_u,
                                       kernel_o,
                                       BANDWIDTH_reg,
                                       num_obs-j-1,
                                       1,
                                       num_reg_unordered,
                                       num_reg_ordered,
                                       num_reg_continuous,
                                       0,
                                       0,
                                       1,
                                       (BANDWIDTH_reg == BW_ADAP_NN)?1:0,
                                       0,
                                       1,
                                       1,
                                       0,
                                       0,
                                       operator,
                                       OP_NOOP,
                                       0,
                                       0,
                                       NULL,
                                       0,
                                       nrc2,
                                       nrc2,
                                       0,
                                       0,
                                       NULL, NULL, NULL, NULL,
                                       PXU,
                                       PXO,
                                       PXC,
                                       TUNO,
                                       TORD,
                                       TCON,
                                       XTKX,
                                       XTKX,
                                       sgn,
                                       vsf,
                                       1,
                                       matrix_bandwidth,
                                       matrix_bandwidth_eval,
                                       lambda,
                                       num_categories,
                                       NULL,
                                       NULL,
                                       kwm+j*nrcc22,
                                       NULL,
                                       NULL,
                                       NULL);

            for(int jj = j+1; jj < num_obs; jj++){
              const double RW = kwm[jj*nrcc22+nrc1]*(XTKX[0][-1]-XTKX[0][jj-j-1]);
              for(int ii = 1; ii < nrc2; ii++)
                kwm[jj*nrcc22+ii*nrc2] += RW*XTKX[ii][jj-j-1]*sgn[ii];
            }
          } else {
            double * const tpk = kwm+j*nrcc22;
            for (int jj = 0; jj < nrc2; jj++){
              for (int ii = nrc1; ii > jj; ii--)
                tpk[jj*nrc2+ii] = tpk[ii*nrc2+jj];
            }
          }
        }
#endif

        if((BANDWIDTH_reg == BW_ADAP_NN)&&(bwm == RBWM_CVAIC)){
          for(l = 0; l < num_reg_continuous; l++)
            pnh /= matrix_bandwidth[l][j];
        }

        if(bwm == RBWM_CVAIC){
          KWM[0][0] += pnh*aicc;
          XTKY[0][0] += pnh*aicc*vector_Y[j];
        }

        while(mat_solve(KWM, XTKY, DELTA) == NULL){
          for(i = 0; i < nrc1; i++)
            KWM[i][i] += epsilon;
          nepsilon += epsilon;
        }

        XTKY[0][0] += nepsilon*XTKY[0][0]/NZD_POS(KWM[0][0]);
        if(nepsilon > 0.0){
          if(mat_solve(KWM, XTKY, DELTA) == NULL){
            glp_ok = 0;
            break;
          }
        }

        for(i = 0; i < nrc1; i++)
          evalv[i] = basis[i][j];
        {
          double mhat = 0.0;
          for(i = 0; i < nrc1; i++)
            mhat += evalv[i]*DELTA[i][0];
          const double dy = vector_Y[j]-mhat;
          const double d2 = dy*dy;
          cv += d2;
        }

        if(bwm == RBWM_CVAIC){
          for(i = 0; i < nrc1; i++)
            XTKY[i][0] = evalv[i];
          if(mat_solve(KWM, XTKY, DELTA) == NULL){
            glp_ok = 0;
            break;
          }
          {
            double hii = 0.0;
            for(i = 0; i < nrc1; i++)
              hii += evalv[i]*DELTA[i][0];
            traceH += hii*pnh*aicc;
          }
        }
      }
    }

    if((XTKY != NULL) && (KWM != NULL)){
      for(i = 0; i < nrc1; i++){
        KWM[i] = PKWM[i];
        XTKY[i] = PXTKY[i];
      }
    }

    if(XTKY != NULL) mat_free(XTKY);
    if(DELTA != NULL) mat_free(DELTA);
    if(KWM != NULL) mat_free(KWM);
    if(TCON != NULL) mat_free(TCON);
    if(TUNO != NULL) mat_free(TUNO);
    if(TORD != NULL) mat_free(TORD);
    if(matrix_bandwidth_eval != NULL) free_tmat(matrix_bandwidth_eval);
    if(XTKX != NULL) free(XTKX);
    if(kwm != NULL) free(kwm);
    if(sgn != NULL) free(sgn);
    if(evalv != NULL) free(evalv);

    if(sf_flag){
      int_LARGE_SF = 0;
      free(vsf);
      vsf = NULL;
    }

    if(!glp_ok){
      cv = DBL_MAX;
      goto finish_cv_path;
    }

    goto finish_cv_path;
  }

  /*
    All-large gate shortcut:
    when every active kernel component is effectively constant, the estimator
    collapses to a global least-squares fit.

    - LC: intercept-only global mean model
    - LL: global linear model on continuous regressors (matches current LL path)

    Then:
    - CVLS: exact LOOCV via e_i/(1-h_ii)
    - CVAIC: SSE term from in-sample residuals, trace(H)=sum(h_ii)
  */
  {
    const int all_large_gate = np_reg_cv_all_large_gate(BANDWIDTH_reg,
                                                        num_obs,
                                                        num_reg_continuous,
                                                        num_reg_unordered,
                                                        num_reg_ordered,
                                                        kernel_c,
                                                        kernel_u,
                                                        kernel_o,
                                                        matrix_X_continuous,
                                                        matrix_bandwidth,
                                                        lambda,
                                                        num_categories,
                                                        matrix_categorical_vals_extern,
                                                        &ov_cont_ok,
                                                        &ov_cont_hmin,
                                                        &ov_cont_k0,
                                                        &ov_cont_from_cache);

    if(all_large_gate){
      const int k = (int_ll_cv == LL_LC) ? 1 : (num_reg_continuous + 1);
      MATRIX XtX = mat_creat(k, k, UNDEFINED);
      MATRIX XtXINV = mat_creat(k, k, UNDEFINED);
      MATRIX XtY = mat_creat(k, 1, UNDEFINED);
      MATRIX BETA = mat_creat(k, 1, UNDEFINED);
      int fast_ok = (XtX != NULL) && (XtXINV != NULL) && (XtY != NULL) && (BETA != NULL);

      if(fast_ok){
        const double ridge_eps = 1.0/(double)MAX(1, num_obs);
        int ridge_it = 0;

        for(i = 0; i < k; i++){
          XtY[i][0] = 0.0;
          BETA[i][0] = 0.0;
          for(j = 0; j < k; j++)
            XtX[i][j] = 0.0;
        }

        for(i = 0; i < num_obs; i++){
          XtX[0][0] += 1.0;
          XtY[0][0] += vector_Y[i];
          if(k > 1){
            for(j = 0; j < num_reg_continuous; j++){
              const double xj = matrix_X_continuous[j][i];
              const int cj = j + 1;
              XtX[0][cj] += xj;
              XtX[cj][0] += xj;
              XtY[cj][0] += xj*vector_Y[i];
            }
            for(int a = 0; a < num_reg_continuous; a++){
              const double xa = matrix_X_continuous[a][i];
              const int ca = a + 1;
              for(int b = a; b < num_reg_continuous; b++){
                const double xb = matrix_X_continuous[b][i];
                const int cb = b + 1;
                XtX[ca][cb] += xa*xb;
                if(cb != ca) XtX[cb][ca] += xa*xb;
              }
            }
          }
        }

        while(mat_inv(XtX, XtXINV) == NULL){
          for(i = 0; i < k; i++)
            XtX[i][i] += ridge_eps;
          ridge_it++;
          if(ridge_it > 64){
            fast_ok = 0;
            break;
          }
        }

        if(fast_ok){
          for(i = 0; i < k; i++){
            double s = 0.0;
            for(j = 0; j < k; j++)
              s += XtXINV[i][j]*XtY[j][0];
            BETA[i][0] = s;
          }

          cv = 0.0;
          traceH = 0.0;
          for(i = 0; i < num_obs; i++){
            double yhat = BETA[0][0];
            double hii = XtXINV[0][0];
            if(k > 1){
              for(j = 0; j < num_reg_continuous; j++){
                const double xj = matrix_X_continuous[j][i];
                yhat += BETA[j+1][0]*xj;
              }
              for(j = 0; j < num_reg_continuous; j++){
                const double xj = matrix_X_continuous[j][i];
                hii += 2.0*xj*XtXINV[0][j+1];
              }
              for(int a = 0; a < num_reg_continuous; a++){
                const double xa = matrix_X_continuous[a][i];
                for(int b = 0; b < num_reg_continuous; b++){
                  const double xb = matrix_X_continuous[b][i];
                  hii += xa*XtXINV[a+1][b+1]*xb;
                }
              }
            }

            const double err = vector_Y[i] - yhat;
            if(bwm == RBWM_CVLS){
              const double den = NZD_POS(1.0 - hii);
              const double err_loo = err/den;
              cv += err_loo*err_loo;
            } else {
              cv += err*err;
              traceH += hii;
            }
          }
        }
      }

      if(XtX != NULL) mat_free(XtX);
      if(XtXINV != NULL) mat_free(XtXINV);
      if(XtY != NULL) mat_free(XtY);
      if(BETA != NULL) mat_free(BETA);

      if(fast_ok){
        np_fastcv_alllarge_hits++;
        goto finish_cv_path;
      }
    }
  }

  if(np_reg_cv_use_canonical_glp_fixed_kernel(int_ll_cv,
                                              bwm,
                                              BANDWIDTH_reg,
                                              num_reg_continuous,
                                              ks_tree_use,
                                              0)){
    NPRegCvLpResult glp_result =
      np_regression_cv_glp_rawbasis_fixed(int_ll_cv,
                                          num_obs,
                                          num_reg_unordered,
                                          num_reg_ordered,
                                          num_reg_continuous,
                                          matrix_X_unordered,
                                          matrix_X_ordered,
                                          matrix_X_continuous,
                                          vector_Y,
                                          vector_scale_factor,
                                          num_categories,
                                          kernel_c,
                                          kernel_u,
                                          kernel_o,
                                          operator,
                                          lambda,
                                          matrix_bandwidth,
                                          NULL,
                                          0,
                                          NULL);
    cv = glp_result.cv;
    traceH = glp_result.traceH;
    goto finish_cv_path;
  }

  if(bwm == RBWM_CVAIC){
    // compute normalisation constant

    // workaround for bwscaling = TRUE
    // really just want to get the full product kernel evaluated at zero
    tsf = int_LARGE_SF;

    int_LARGE_SF = 1;

    kernel_weighted_sum_np_ctx(kernel_c,
                           kernel_u,
                           kernel_o,
                           BANDWIDTH_reg,
                           1,
                           1,
                           num_reg_unordered,
                           num_reg_ordered,
                           num_reg_continuous,
                           0, // do not leave out 
                           0,
                           1, // kernel_pow = 1
                           0, // bandwidth_divide = FALSE 
                           0, 
                           0, // not symmetric
                           0, // do not gather-scatter
                           0, // do not drop train
                           0, // do not drop train
                           operator, // all regressors use the normal kernels (not cdf or derivative ones) 
                           OP_NOOP, // no permutations
                           0, // no score
                           0, // no ocg
                           NULL,
                           1, // explicitly suppress parallel
                           0, // no Y
                           0, // no weights
                           NP_TREE_FALSE, // disable the tree because we are fiddling with the data
                           0,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           matrix_X_unordered, // TRAIN
                           matrix_X_ordered,
                           matrix_X_continuous,
                           matrix_X_unordered, // EVAL
                           matrix_X_ordered,
                           matrix_X_continuous,
                           NULL,
                           NULL,
                           NULL,
                           vector_scale_factor,
                           1,
                           matrix_bandwidth,
                           matrix_bandwidth,
                           lambda,
                           num_categories,
                           NULL,
                           NULL,
                           &aicc,
                           NULL, // no permutations
                           NULL, // do not return kernel weights
                           NULL);
    int_LARGE_SF = tsf;

    //fprintf(stderr,"\n%e\n",aicc);
  }


  // Conduct the estimation 

  if(int_ll_cv == LL_LC) { // local constant
    // Nadaraya-Watson
    // Generate bandwidth vector given scale factors, nearest neighbors, or lambda 

    double * lc_Y[2];
    double * mean = (double *)malloc(2*num_obs_eval_alloc*sizeof(double));

    lc_Y[0] = vector_Y;
      
    lc_Y[1] = (double *)malloc(num_obs*sizeof(double));
    for(int ii = 0; ii < num_obs; ii++)
      lc_Y[1][ii] = 1.0;

    kernel_weighted_sum_np_ctx(kernel_c,
                           kernel_u,
                           kernel_o,
                           BANDWIDTH_reg,
                           num_obs,
                           num_obs,
                           num_reg_unordered,
                           num_reg_ordered,
                           num_reg_continuous,
                           leave_one_out, 
                           0,
                           1, // kernel_pow = 1
                           (BANDWIDTH_reg == BW_ADAP_NN)?1:0, // bandwidth_divide = FALSE when not adaptive
                           0, 
                           0, // not symmetric
                           0, // do not gather-scatter
                           0, // do not drop train
                           0, // do not drop train
                           operator, // no special operators being used
                           OP_NOOP, // no permutations
                           0, // no score
                           0, // no ocg
                           NULL,
                           0, // don't explicitly suppress parallel
                           2, // 2 cols in Y
                           0, // 0 cols in W
                           int_TREE_X,
                           0,
                           kdt_extern_X,
                           NULL,
                           NULL,
                           NULL,
                           matrix_X_unordered, // TRAIN
                           matrix_X_ordered,
                           matrix_X_continuous,
                           matrix_X_unordered, // EVAL
                           matrix_X_ordered,
                           matrix_X_continuous,
                           lc_Y,
                           NULL,
                           NULL,
                           vector_scale_factor,
                           1,
                           matrix_bandwidth,
                           matrix_bandwidth,
                           lambda,
                           num_categories,
                           NULL,
                           NULL,
                           mean,
                           NULL, // no permutations
                           NULL, // do not return kernel weights
                           NULL);

    // every even entry in mean is sum(y*kij)
    // every odd is sum(kij)

    for(int ii = 0; ii < num_obs; ii++){
      const int ii2 = 2*ii;
      const double sk = copysign(DBL_MIN, mean[ii2+1]) + mean[ii2+1];
      const double dy = vector_Y[ii]-mean[ii2]/sk;
      cv += dy*dy;
      if(bwm == RBWM_CVAIC){
        if(BANDWIDTH_reg != BW_ADAP_NN){
          traceH += aicc/sk;
        }else{
          double pnh = 1.0;
          for(int jj = 0; jj < num_reg_continuous; jj++)
            pnh /= matrix_bandwidth[jj][ii];
          traceH += pnh*aicc/sk;
        }
        
      }
      //fprintf(stderr,"mj: %e\n",mean[ii2]/(MAX(DBL_MIN, mean[ii2+1])));
    }

    free(lc_Y[1]);
    free(mean);
  } else { // Local Linear 

    // because we manipulate the training data scale factors can be wrong

    if((sf_flag = (int_LARGE_SF == 0))){ 
      int_LARGE_SF = 1;
      vsf = (double *)malloc(num_reg_continuous*sizeof(double));
      for(int ii = 0; ii < num_reg_continuous; ii++)
        vsf[ii] = matrix_bandwidth[ii][0];
    } else {
      vsf = vector_scale_factor;
    }

    MATRIX XTKX = mat_creat( num_reg_continuous + 2, num_obs, UNDEFINED );
    MATRIX XTKXINV = mat_creat( num_reg_continuous + 1, num_reg_continuous + 1, UNDEFINED );
    MATRIX XTKY = mat_creat( num_reg_continuous + 1, 1, UNDEFINED );
    MATRIX DELTA = mat_creat( num_reg_continuous + 1, 1, UNDEFINED );

    MATRIX KWM = mat_creat( num_reg_continuous + 1, num_reg_continuous + 1, UNDEFINED );
    // Generate bandwidth vector given scale factors, nearest neighbors, or lambda 
    
    MATRIX TCON = mat_creat(num_reg_continuous, 1, UNDEFINED);
    MATRIX TUNO = mat_creat(num_reg_unordered, 1, UNDEFINED);
    MATRIX TORD = mat_creat(num_reg_ordered, 1, UNDEFINED);

    const int nrc2 = (num_reg_continuous+2);
    const int nrc1 = (num_reg_continuous+1);
    const int nrcc22 = nrc2*nrc2;

    double ** matrix_bandwidth_eval = NULL;

    double * PKWM[nrc1], * PXTKY[nrc1], * PXTKX[nrc2];

    double * PXC[MAX(1,num_reg_continuous)]; 
    double * PXU[MAX(1,num_reg_unordered)];
    double * PXO[MAX(1,num_reg_ordered)];

    PXC[0] = NULL;
    PXU[0] = NULL;
    PXO[0] = NULL;

    for(l = 0; l < num_reg_continuous; l++)
      PXC[l] = matrix_X_continuous[l];

    for(l = 0; l < num_reg_unordered; l++)
      PXU[l] = matrix_X_unordered[l];

    for(l = 0; l < num_reg_ordered; l++)
      PXO[l] = matrix_X_ordered[l];

    const size_t kwm_len = (size_t)nrcc22*(size_t)num_obs_eval_alloc;
    double * kwm = (double *)malloc(kwm_len*sizeof(double));

    for(size_t ii = 0; ii < kwm_len; ii++)
      kwm[ii] = 0.0;

    double * sgn = (double *)malloc((nrc2)*sizeof(double));

    sgn[0] = sgn[1] = 1.0;
    
    for(int ii = 0; ii < (num_reg_continuous); ii++)
      sgn[ii+2] = -1.0;
    
    for(int ii = 0; ii < (nrc2); ii++)
      PXTKX[ii] = XTKX[ii];
    
    for(int ii = 0; ii < (nrc1); ii++){
      PKWM[ii] = KWM[ii];
      PXTKY[ii] = XTKY[ii];

      KWM[ii] = &kwm[(ii+1)*(nrc2)+1];
      XTKY[ii] = &kwm[ii+1];

    }

    matrix_bandwidth_eval = alloc_tmatd(1,num_reg_continuous);


    const double epsilon = 1.0/num_obs;
    double nepsilon;

    //    matrix_bandwidth = alloc_matd(num_obs,num_reg_continuous);

    // populate the xtkx matrix first 
    
    for(i = 0; i < num_obs; i++){
      XTKX[0][i] = vector_Y[i];
      XTKX[1][i] = 1.0;
    }


    for(j = 0; j < num_obs; j++){ // main loop
      nepsilon = 0.0;

      for(l = 0; l < (nrc1); l++){
        KWM[l] = &kwm[j*nrcc22+(l+1)*(nrc2)+1];
        XTKY[l] = &kwm[j*nrcc22+l+1];
      }

#ifdef MPI2
      if(np_reg_cv_use_symmetric_dropone_path(bwm, ks_tree_use, BANDWIDTH_reg)){
        if((j % iNum_Processors) == 0){
          if((j+my_rank) < (num_obs)){
            for(l = 0; l < num_reg_continuous; l++){
          
              for(i = 0; i < num_obs; i++){
                XTKX[l+2][i] = matrix_X_continuous[l][i]-matrix_X_continuous[l][j+my_rank];
              }
              TCON[l][0] = matrix_X_continuous[l][j+my_rank]; // temporary storage

              if(BANDWIDTH_reg == BW_GEN_NN)
                matrix_bandwidth_eval[l][0] = matrix_bandwidth[l][j+my_rank]; // temporary storage
            }


            for(l = 0; l < num_reg_unordered; l++)
              TUNO[l][0] = matrix_X_unordered[l][j+my_rank];

            for(l = 0; l < num_reg_ordered; l++)
              TORD[l][0] = matrix_X_ordered[l][j+my_rank];

            kernel_weighted_sum_np_ctx(kernel_c,
                                   kernel_u,
                                   kernel_o,
                                   BANDWIDTH_reg,
                                   num_obs,
                                   1,
                                   num_reg_unordered,
                                   num_reg_ordered,
                                   num_reg_continuous,
                                   0, 
                                   0,
                                   1, // kernel_pow = 1
                                   (BANDWIDTH_reg == BW_ADAP_NN)?1:0, // bandwidth_divide = FALSE when not adaptive
                                   0, 
                                   1, // symmetric
                                   0, // NO gather-scatter sum
                                   1, // drop train
                                   j+my_rank, // drop this training datum
                                   operator, // no convolution
                                   OP_NOOP, // no permutations
                                   0, // no score
                                   0, // no ocg
                                   NULL,
                                   1, // explicitly suppress parallel
                                   nrc2, // nrc2 cols in Y
                                   nrc2, // nrc2 cols in W
                                   (BANDWIDTH_reg == BW_ADAP_NN) ? NP_TREE_FALSE : int_TREE_X,
                                   0,
                                   (BANDWIDTH_reg == BW_ADAP_NN) ? NULL : kdt_extern_X,
                                   NULL,
                                   NULL,
                                   NULL,
                                   PXU, // TRAIN
                                   PXO, 
                                   PXC,
                                   TUNO, // EVAL
                                   TORD,
                                   TCON,
                                   XTKX,
                                   XTKX,
                                   NULL,
                                   vsf,
                                   1,
                                   matrix_bandwidth,
                                   matrix_bandwidth_eval,
                                   lambda,
                                   num_categories,
                                   NULL,
                                   NULL,
                                   kwm+(j+my_rank)*nrcc22,  // weighted sum
                                   NULL, // no permutations
                                   NULL, // do not return kernel weights
                                   NULL);

          }
          // synchro step
          MPI_Allgather(MPI_IN_PLACE, nrcc22, MPI_DOUBLE, kwm+j*nrcc22, nrcc22, MPI_DOUBLE, comm[1]);
        }
      } else {
        if((j % iNum_Processors) == 0){

          // some guys have to sit out the last calculation, but
          // they still sync up afterwards
          if((j+my_rank) < (num_obs-1)){

            for(l = 0; l < (nrc2); l++){
              XTKX[l] = PXTKX[l] + j + my_rank + 1;
            }

            for(l = 0; l < num_reg_continuous; l++)
              PXC[l] = matrix_X_continuous[l] + j + my_rank + 1;

            for(l = 0; l < num_reg_unordered; l++)
              PXU[l] = matrix_X_unordered[l] + j + my_rank + 1;

            for(l = 0; l < num_reg_ordered; l++)
              PXO[l] = matrix_X_ordered[l] + j + my_rank + 1;
        
            for(l = 0; l < num_reg_continuous; l++){
        
              for(i = 0; i < (num_obs-j-1-my_rank); i++){
                XTKX[l+2][i] = matrix_X_continuous[l][i+j+1+my_rank]-matrix_X_continuous[l][j+my_rank];
              }
              TCON[l][0] = matrix_X_continuous[l][j+my_rank]; // temporary storage

              if(BANDWIDTH_reg == BW_GEN_NN)
                matrix_bandwidth_eval[l][0] = matrix_bandwidth[l][j+my_rank]; // temporary storage
            }

            for(l = 0; l < num_reg_unordered; l++)
              TUNO[l][0] = matrix_X_unordered[l][j+my_rank];

            for(l = 0; l < num_reg_ordered; l++)
              TORD[l][0] = matrix_X_ordered[l][j+my_rank];

            kernel_weighted_sum_np_ctx(kernel_c,
                                   kernel_u,
                                   kernel_o,
                                   BANDWIDTH_reg,
                                   num_obs-j-my_rank-1,
                                   1,
                                   num_reg_unordered,
                                   num_reg_ordered,
                                   num_reg_continuous,
                                   0, // we leave one out via the weight matrix
                                   0,
                                   1, // kernel_pow = 1
                                   (BANDWIDTH_reg == BW_ADAP_NN)?1:0, // bandwidth_divide = FALSE when not adaptive
                                   0, 
                                   1, // symmetric
                                   1, // gather-scatter sum
                                   0, // do not drop train
                                   0, // do not drop train
                                   operator, // no convolution
                                   OP_NOOP, // no permutations
                                   0, // no score
                                   0, // no ocg
                                   NULL,
                                   1, // explicitly suppress parallel
                                   nrc2, // cols in Y
                                   nrc2, // cols in W
                                   0, // no tree?
                                   0,
                                   NULL, NULL, NULL, NULL,
                                   PXU, // TRAIN
                                   PXO, 
                                   PXC,
                                   TUNO, // EVAL
                                   TORD,
                                   TCON,
                                   XTKX,
                                   XTKX,
                                   sgn,
                                   vsf,
                                   1,
                                   matrix_bandwidth,
                                   matrix_bandwidth_eval,
                                   lambda,
                                   num_categories,
                                   NULL,
                                   NULL,
                                   kwm+(j+my_rank)*nrcc22,  // weighted sum
                                   NULL, // no permutations
                                   NULL, // do not return kernel weights
                                   NULL);

            // need to use reference weight to fix weight sum
            for(int jj = j+my_rank+1; jj < num_obs; jj++){
              const double RW = kwm[jj*nrcc22+nrc1]*(XTKX[0][-1]-XTKX[0][jj-j-my_rank-1]);
              for(int ii = 1; ii < nrc2; ii++){
                kwm[jj*nrcc22+ii*nrc2] += RW*XTKX[ii][jj-j-my_rank-1]*sgn[ii];
              }
            }
          }
          // reduce all work arrays
          const int nrem = num_obs % iNum_Processors;
          const int nred = ((j+iNum_Processors) > num_obs) ? nrem : iNum_Processors;

          MPI_Allreduce(MPI_IN_PLACE, kwm+j*nrcc22, nred*nrcc22, MPI_DOUBLE, MPI_SUM, comm[1]);
        }

        // due to a quirk of the algorithm in parallel, always need to re-symmetrise array
        
        double * const tpk = kwm+j*nrcc22;
        for (int jj = 0; jj < (nrc2); jj++){
          for (int ii = (nrc1); ii > jj; ii--){
            tpk[jj*(nrc2)+ii] = tpk[ii*(nrc2)+jj];
          }
        }
      }

#else
      if(np_reg_cv_use_symmetric_dropone_path(bwm, ks_tree_use, BANDWIDTH_reg)){

        for(l = 0; l < num_reg_continuous; l++){
          
          for(i = 0; i < num_obs; i++){
            XTKX[l+2][i] = matrix_X_continuous[l][i]-matrix_X_continuous[l][j];
          }
          TCON[l][0] = matrix_X_continuous[l][j]; // temporary storage

          if(BANDWIDTH_reg == BW_GEN_NN)
            matrix_bandwidth_eval[l][0] = matrix_bandwidth[l][j]; // temporary storage
        }


        for(l = 0; l < num_reg_unordered; l++)
          TUNO[l][0] = matrix_X_unordered[l][j];

        for(l = 0; l < num_reg_ordered; l++)
          TORD[l][0] = matrix_X_ordered[l][j];

        kernel_weighted_sum_np_ctx(kernel_c,
                               kernel_u,
                               kernel_o,
                               BANDWIDTH_reg,
                               num_obs,
                               1,
                               num_reg_unordered,
                               num_reg_ordered,
                               num_reg_continuous,
                               0, // we leave one out via the weight matrix
                               0, 
                               1, // kernel_pow = 1
                               (BANDWIDTH_reg == BW_ADAP_NN)?1:0, // bandwidth_divide = FALSE when not adaptive
                               0, 
                               1, // symmetric
                               0, // gather-scatter sum
                               1, // do not drop train
                               j, // do not drop train
                               operator, // no convolution
                               OP_NOOP, // no permutations
                               0, // no score
                               0, // no ocg
                               NULL,
                               0, // don't explicity suppress parallel
                               nrc2,
                               nrc2,
                               (BANDWIDTH_reg == BW_ADAP_NN) ? NP_TREE_FALSE : int_TREE_X,
                               0,
                               (BANDWIDTH_reg == BW_ADAP_NN) ? NULL : kdt_extern_X,
                               NULL, NULL, NULL,
                               PXU, // TRAIN
                               PXO, 
                               PXC,
                               TUNO, // EVAL
                               TORD,
                               TCON,
                               XTKX,
                               XTKX,
                               NULL,
                               vsf,
                               1,
                               matrix_bandwidth,
                               matrix_bandwidth_eval,
                               lambda,
                               num_categories,
                               NULL,
                               NULL,
                               kwm+j*nrcc22,  // weighted sum
                               NULL, // no permutations
                               NULL, // do not return kernel weights
                               NULL);

      } else {
        if(j < (num_obs-1)){

          for(l = 0; l < (nrc2); l++){
            XTKX[l]++;
          }

          for(l = 0; l < num_reg_continuous; l++)
            PXC[l]++;

          for(l = 0; l < num_reg_unordered; l++)
            PXU[l]++;

          for(l = 0; l < num_reg_ordered; l++)
            PXO[l]++;

          for(l = 0; l < num_reg_continuous; l++){
          
            for(i = 0; i < (num_obs-j-1); i++){
              XTKX[l+2][i] = matrix_X_continuous[l][i+j+1]-matrix_X_continuous[l][j];
            }
            TCON[l][0] = matrix_X_continuous[l][j]; // temporary storage

            if(BANDWIDTH_reg == BW_GEN_NN)
              matrix_bandwidth_eval[l][0] = matrix_bandwidth[l][j]; // temporary storage
          }

      
          for(l = 0; l < num_reg_unordered; l++)
            TUNO[l][0] = matrix_X_unordered[l][j];

          for(l = 0; l < num_reg_ordered; l++)
            TORD[l][0] = matrix_X_ordered[l][j];
      
          kernel_weighted_sum_np_ctx(kernel_c,
                                 kernel_u,
                                 kernel_o,
                                 BANDWIDTH_reg,
                                 num_obs-j-1,
                                 1,
                                 num_reg_unordered,
                                 num_reg_ordered,
                                 num_reg_continuous,
                                 0, // we leave one out via the weight matrix
                                 0,
                                 1, // kernel_pow = 1
                                 (BANDWIDTH_reg == BW_ADAP_NN)?1:0, // bandwidth_divide = FALSE when not adaptive
                                 0, 
                                 1, // symmetric
                                 1, // gather-scatter sum
                                 0, // do not drop train
                                 0, // do not drop train
                                 operator, // no convolution
                                 OP_NOOP, // no permutations
                                 0, // no score
                                 0, // no ocg
                                 NULL,
                                 0, // don't explicity suppress parallel
                                 nrc2,
                                 nrc2,
                                 0, // no trees
                                 0,
                                 NULL, NULL, NULL, NULL,
                                 PXU, // TRAIN
                                 PXO, 
                                 PXC,
                                 TUNO, // EVAL
                                 TORD,
                                 TCON,
                                 XTKX,
                                 XTKX,
                                 sgn,
                                 vsf,
                                 1,
                                 matrix_bandwidth,
                                 matrix_bandwidth_eval,
                                 lambda,
                                 num_categories,
                                 NULL,
                                 NULL,
                                 kwm+j*nrcc22, // weighted sum
                                 NULL, // no permutations
                                 NULL,  // no kernel weights
                                 NULL);

          // need to use reference weight to fix weight sum
          for(int jj = j+1; jj < num_obs; jj++){
            const double RW = kwm[jj*nrcc22+nrc1]*(XTKX[0][-1]-XTKX[0][jj-j-1]);
            for(int ii = 1; ii < nrc2; ii++){
              kwm[jj*nrcc22+ii*nrc2] += RW*XTKX[ii][jj-j-1]*sgn[ii];
            }
          }
        } else { // because we skip the last call to npksum, we need to copy L to U for last observation
          double * const tpk = kwm+j*nrcc22;
          for (int jj = 0; jj < (nrc2); jj++){
            for (int ii = (nrc1); ii > jj; ii--){
              tpk[jj*(nrc2)+ii] = tpk[ii*(nrc2)+jj];
            }
          }
        }
      }
#endif
      double pnh = 1.0;

      if((BANDWIDTH_reg == BW_ADAP_NN)&&(bwm == RBWM_CVAIC)){
        for(int jj = 0; jj < num_reg_continuous; jj++)
          pnh /= matrix_bandwidth[jj][j];
      }

      // need to manipulate KWM pointers and XTKY - done
      if(bwm == RBWM_CVAIC){
        KWM[0][0] += pnh*aicc;
        XTKY[0][0] += pnh*aicc*vector_Y[j];
      }

      while(mat_solve(KWM, XTKY, DELTA) == NULL){ // singular = ridge about
        for(int ii = 0; ii < (nrc1); ii++)
          KWM[ii][ii] += epsilon;
        nepsilon += epsilon;
      }

      if(bwm == RBWM_CVAIC){
        int ok00 = 0;
        const double inv00 = mat_inv00(KWM, &ok00);
        if(!ok00)
          error("mat_inv00 failed after ridge adjustment");
        traceH += inv00*pnh*aicc;
      }

      XTKY[0][0] += nepsilon*XTKY[0][0]/NZD(KWM[0][0]);
      if(nepsilon > 0.0){
        if(mat_solve(KWM, XTKY, DELTA) == NULL)
          error("mat_solve failed after ridge adjustment");
      }
      const double dy = vector_Y[j]-DELTA[0][0];
      cv += dy*dy; 
    }

    for(int ii = 0; ii < (nrc1); ii++){
      KWM[ii] = PKWM[ii];
      XTKY[ii] = PXTKY[ii];
    }

    for(int ii = 0; ii < (nrc2); ii++)
      XTKX[ii] = PXTKX[ii];

    if(sf_flag){
      int_LARGE_SF = 0;
      free(vsf);
    }

    
    mat_free(XTKX);
    mat_free(XTKXINV);
    mat_free(XTKY);
    mat_free(DELTA);
    mat_free(KWM);

    mat_free(TCON);
    mat_free(TUNO);
    mat_free(TORD);

    free(kwm);
    free(sgn);
    free_tmat(matrix_bandwidth_eval);
  }

finish_cv_path:
  if((ov_cont_ok != NULL) && (!ov_cont_from_cache)) free(ov_cont_ok);
  if((ov_cont_hmin != NULL) && (!ov_cont_from_cache)) free(ov_cont_hmin);
  if((ov_cont_k0 != NULL) && (!ov_cont_from_cache)) free(ov_cont_k0);

	/* Negative penalties are treated as infinite: Hurvich et al pg 277 */

  cv /= (double)num_obs;
  if(bwm == RBWM_CVAIC){
    if((1.0+traceH/((double)num_obs))/(1.0-(traceH+2.0)/((double)num_obs)) < 0) {
      cv = DBL_MAX;
    } else {
      cv = log(cv) + (1.0+traceH/((double)num_obs))/(1.0-(traceH+2.0)/((double)num_obs));
    }
  }

  //Rprintf("cv: %3.15g ",cv);
  //for(int ii = 0; ii < num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern; ii++)
  //  Rprintf("%3.15g ", vector_scale_factor[ii]);
  //Rprintf("\n");

  return(cv);
}

double np_kernel_estimate_distribution_ls_cv( 
int KERNEL_den,
int KERNEL_den_unordered,
int KERNEL_den_ordered,
int BANDWIDTH_den,
int num_obs_train,
int num_obs_eval,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
int cdfontrain,
double memfac,
double ** matrix_X_unordered_train,
double ** matrix_X_ordered_train,
double ** matrix_X_continuous_train,
double ** matrix_X_unordered_eval,
double ** matrix_X_ordered_eval,
double ** matrix_X_continuous_eval,
double * vsf,
int * num_categories,
double ** matrix_categorical_vals,
double * cv){
  NP_GateOverrideCtx gate_ctx_local;
  np_gate_ctx_clear(&gate_ctx_local);
  const int bwmdim = (BANDWIDTH_den==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_den==BW_ADAP_NN)?num_obs_train:1);

  int indy;

  int64_t i,j,l,iwx;

  int * operator = NULL;
  int gate_override_active = 0;
  int all_large_gate = 0;
  int *ov_cont_ok = NULL;
  double *ov_cont_hmin = NULL, *ov_cont_k0 = NULL;
  double **matrix_bandwidth = NULL;
  double *lambda = NULL;

  double **matrix_wX_unordered_eval=NULL;
  double **matrix_wX_ordered_eval=NULL;
  double **matrix_wX_continuous_eval=NULL;

  int64_t N, num_obs_eval_alloc, num_obs_train_alloc, num_obs_wx_alloc;
  int64_t wx, nwx;

  size_t Nm = MIN((size_t)ceil(memfac*300000.0), (size_t)SIZE_MAX/10);

#ifdef MPI2
  int64_t stride_t = MAX((int64_t)ceil((double) num_obs_train / (double) iNum_Processors),1);
  int64_t stride_e = MAX((int64_t)ceil((double) num_obs_eval / (double) iNum_Processors),1);
  
  num_obs_train_alloc = stride_t*iNum_Processors;
  num_obs_eval_alloc = stride_e*iNum_Processors;
#else
  num_obs_train_alloc = num_obs_train;
  num_obs_eval_alloc = num_obs_eval;
#endif

  // blocking algo calculations
  N = num_obs_eval_alloc*(num_obs_train_alloc+1);
  
  const int64_t sa = num_obs_eval_alloc*num_obs_train_alloc*sizeof(double);

  if((N > Nm) || (sa > (((int64_t)1<<31)-1))){
    const int64_t wx0 = Nm/(1+num_obs_train_alloc);
    wx = (wx0 > num_obs_eval_alloc) ? num_obs_eval_alloc : wx0;
    nwx = num_obs_eval_alloc/wx + (((num_obs_eval_alloc % wx) > 0) ? 1 : 0);
  } else {
    wx = num_obs_eval_alloc;
    nwx = 1;
  }

#ifdef MPI2
  int64_t stride_wx = MAX((int64_t)ceil((double)wx / (double) iNum_Processors),1);

  num_obs_wx_alloc = stride_wx*iNum_Processors;
#else
  num_obs_wx_alloc = wx;
#endif

#ifdef MPI2
#endif

  // allocate some pointers
  matrix_wX_continuous_eval = (double **)np_jksum_malloc_array_or_die((size_t)num_reg_continuous, sizeof(double *), "np_kernel_estimate_density_categorical_leave_one_out_cv matrix_wX_continuous_eval");
  matrix_wX_unordered_eval = (double **)np_jksum_malloc_array_or_die((size_t)num_reg_unordered, sizeof(double *), "np_kernel_estimate_density_categorical_leave_one_out_cv matrix_wX_unordered_eval");
  matrix_wX_ordered_eval = (double **)np_jksum_malloc_array_or_die((size_t)num_reg_ordered, sizeof(double *), "np_kernel_estimate_density_categorical_leave_one_out_cv matrix_wX_ordered_eval");
 
  double * mean = (double *)np_jksum_malloc_array_or_die((size_t)num_obs_wx_alloc, sizeof(double), "np_kernel_estimate_density_categorical_leave_one_out_cv mean");

  double ofac = num_obs_train - 1.0;

  operator = (int *)np_jksum_malloc_array_or_die((size_t)(num_reg_continuous+num_reg_unordered+num_reg_ordered), sizeof(int), "np_kernel_estimate_density_categorical_leave_one_out_cv operator");

  for(i = 0; i < (num_reg_continuous+num_reg_unordered+num_reg_ordered); i++)
    operator[i] = OP_INTEGRAL;

  int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;

  kernel_c = (int *)np_jksum_malloc_array_or_die((size_t)num_reg_continuous, sizeof(int), "np_kernel_estimate_density_categorical_leave_one_out_cv kernel_c");

  for(i = 0; i < num_reg_continuous; i++)
    kernel_c[i] = KERNEL_den;

  kernel_u = (int *)np_jksum_malloc_array_or_die((size_t)num_reg_unordered, sizeof(int), "np_kernel_estimate_density_categorical_leave_one_out_cv kernel_u");

  for(i = 0; i < num_reg_unordered; i++)
    kernel_u[i] = KERNEL_den_unordered;

  kernel_o = (int *)np_jksum_malloc_array_or_die((size_t)num_reg_ordered, sizeof(int), "np_kernel_estimate_density_categorical_leave_one_out_cv kernel_o");

  for(i = 0; i < num_reg_ordered; i++)
    kernel_o[i] = KERNEL_den_ordered;

  matrix_bandwidth = alloc_matd(bwmdim,num_reg_continuous);
  lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

  if(kernel_bandwidth_mean(KERNEL_den,
                           BANDWIDTH_den,
                           num_obs_train,
                           num_obs_eval,
                           0,0,0,
                           num_reg_continuous,
                           num_reg_unordered,
                           num_reg_ordered,
                           0,
                           vsf,
                           NULL,NULL,
                           matrix_X_continuous_train,
                           matrix_X_continuous_eval,
                           NULL,
                           matrix_bandwidth,
                           lambda)==1){
    error("\n** Error: invalid bandwidth.");
  }

  if(num_reg_continuous > 0){
    int ok_all = 1;
    ov_cont_ok = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
    ov_cont_hmin = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
    ov_cont_k0 = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
    ok_all = (ov_cont_ok != NULL) && (ov_cont_hmin != NULL) && (ov_cont_k0 != NULL);

    if(ok_all){
      const double rel_tol = np_cont_largeh_rel_tol();
      for(i = 0; i < num_reg_continuous; i++){
        const int kern = kernel_c[i];
        double xmin = DBL_MAX, xmax = -DBL_MAX;

        ov_cont_ok[i] = 0;
        ov_cont_hmin[i] = DBL_MAX;
        ov_cont_k0[i] = 0.0;
        if(!np_cont_largeh_kernel_supported(kern)) continue;

        for(j = 0; j < num_obs_train; j++){
          const double v = matrix_X_continuous_train[i][j];
          if(!isfinite(v)) continue;
          xmin = MIN(xmin, v);
          xmax = MAX(xmax, v);
        }

        for(j = 0; j < num_obs_eval; j++){
          const double v = matrix_X_continuous_eval[i][j];
          if(!isfinite(v)) continue;
          xmin = MIN(xmin, v);
          xmax = MAX(xmax, v);
        }

        if(xmax >= xmin){
          const double utol = np_cont_largeh_utol(kern, rel_tol);
          if(utol > 0.0 && isfinite(utol)){
            ov_cont_ok[i] = 1;
            ov_cont_hmin[i] = (xmax - xmin)/utol;
            ov_cont_k0[i] = np_cont_largeh_k0(kern);
          }
        }
      }

      np_gate_ctx_set(&gate_ctx_local,
                      num_reg_continuous,
                      0,
                      0,
                      kernel_c,
                      kernel_u,
                      kernel_o,
                      operator,
                      ov_cont_ok,
                      ov_cont_hmin,
                      ov_cont_k0,
                      NULL,
                      NULL,
                      NULL,
                      NULL);
      gate_override_active = 1;
    }
  }

  all_large_gate = (BANDWIDTH_den == BW_FIXED) && gate_override_active;
  if(all_large_gate){
    for(i = 0; i < num_reg_continuous; i++){
      const double bw = matrix_bandwidth[i][0];
      if((ov_cont_ok == NULL) || (!ov_cont_ok[i]) || (!isfinite(bw)) ||
         (bw <= 0.0) || (bw < ov_cont_hmin[i])){
        all_large_gate = 0;
        break;
      }
    }
  }
  if(all_large_gate)
    np_fastcv_alllarge_hits++;

  
  *cv = 0;

  double * kwx = (double *)np_jksum_malloc_array3_or_die((size_t)num_obs_train_alloc, (size_t)num_obs_wx_alloc, sizeof(double), "np_kernel_estimate_density_categorical_leave_one_out_cv kwx");

  for(iwx = 0; iwx < nwx; iwx++){
    const int64_t wxo = iwx*wx;
    const int64_t dwx = (iwx != (nwx - 1)) ? wx : num_obs_eval - (nwx - 1)*wx;

    for(l = 0; l < num_reg_continuous; l++)
      matrix_wX_continuous_eval[l] = matrix_X_continuous_eval[l] + wxo;

    for(l = 0; l < num_reg_unordered; l++)
      matrix_wX_unordered_eval[l] = matrix_X_unordered_eval[l] + wxo;

    for(l = 0; l < num_reg_ordered; l++)
      matrix_wX_ordered_eval[l] = matrix_X_ordered_eval[l] + wxo;


    kernel_weighted_sum_np_ctx_ex(kernel_c,
                              kernel_u,
                              kernel_o,
                              BANDWIDTH_den,
                              num_obs_train,
                              dwx,
                              num_reg_unordered,
                              num_reg_ordered,
                              num_reg_continuous,
                              0,
                              0,
                              1,
                              1,
                              0, 
                              0,
                              0,
                              0,
                              0,
                              operator,
                              OP_NOOP, // no permutations
                              0, // no score
                              0, // no ocg
                              NULL,
                              0, // don't explicity suppress parallel
                              0,
                              0,
                              int_TREE_X,
                              0,
                              kdt_extern_X, 
                              NULL, NULL, NULL,
                              matrix_X_unordered_train,
                              matrix_X_ordered_train,
                              matrix_X_continuous_train,
                              matrix_wX_unordered_eval,
                              matrix_wX_ordered_eval,
                              matrix_wX_continuous_eval,
                              NULL,
                              NULL,
                              NULL,
                              vsf,
                              1,matrix_bandwidth,matrix_bandwidth,lambda,
                              num_categories,
                              matrix_categorical_vals,
                              NULL,
                              mean,
                              NULL, // no permutations
                              kwx,
                              &gate_ctx_local,
                              1);
    
#ifdef MPI2
    {
      const int64_t stride_local_eval = MAX((int64_t)ceil((double) dwx / (double) iNum_Processors), 1);
      const int64_t js_local_eval = stride_local_eval * my_rank;
      const int64_t je_local_eval = MIN(dwx - 1, js_local_eval + stride_local_eval - 1);

      for(j = js_local_eval; j <= je_local_eval; j++){
        const int64_t j_global = wxo + j;

        for(i = 0; i < num_obs_train; i++){
          if(cdfontrain && (j_global == i)) continue;
          indy = 1;
          for(l = 0; (l < num_reg_ordered) && (indy != 0); l++){
            indy *= (matrix_X_ordered_train[l][i] <= matrix_X_ordered_eval[l][j_global]);
          }
          for(l = 0; (l < num_reg_continuous) && (indy != 0); l++){
            indy *= (matrix_X_continuous_train[l][i] <= matrix_X_continuous_eval[l][j_global]);
          }
          if(BANDWIDTH_den != BW_ADAP_NN){
            const double tvd = (indy - mean[j]/ofac + kwx[j*num_obs_train + i]/ofac);
            *cv += tvd*tvd;
          } else {
            const double tvd = (indy - mean[j]/ofac + kwx[i*dwx + j]/ofac);
            *cv += tvd*tvd;
          }
        }
      }
    }
#else
    for(i = 0; i < num_obs_train; i++){
      for(j = wxo; j < (wxo + dwx); j++){
        const int64_t jo = j - wxo;
        if(cdfontrain && (j == i)) continue;
        indy = 1;
        for(l = 0; (l < num_reg_ordered) && (indy != 0); l++){
          indy *= (matrix_X_ordered_train[l][i] <= matrix_X_ordered_eval[l][j]);
        }
        for(l = 0; (l < num_reg_continuous) && (indy != 0); l++){
          indy *= (matrix_X_continuous_train[l][i] <= matrix_X_continuous_eval[l][j]);
        }
        if(BANDWIDTH_den != BW_ADAP_NN){
          const double tvd = (indy - mean[jo]/ofac + kwx[jo*num_obs_train + i]/ofac);
          *cv += tvd*tvd;
        } else {
          const double tvd = (indy - mean[jo]/ofac + kwx[i*dwx + jo]/ofac);
          *cv += tvd*tvd;
        }
      }
    }
#endif
  }
#ifdef MPI2
  MPI_Allreduce(MPI_IN_PLACE, cv, 1, MPI_DOUBLE, MPI_SUM, comm[1]);
#endif

  *cv /= (double) num_obs_train*num_obs_eval;

  free(kwx);

  free(operator);
  free(kernel_c);
  free(kernel_u);
  free(kernel_o);
  free(mean);
  free(lambda);
  free_mat(matrix_bandwidth, num_reg_continuous);

  free(matrix_wX_continuous_eval);
  free(matrix_wX_unordered_eval);
  free(matrix_wX_ordered_eval);
  if(gate_override_active)
    np_gate_ctx_clear(&gate_ctx_local);
  if(ov_cont_ok != NULL) free(ov_cont_ok);
  if(ov_cont_hmin != NULL) free(ov_cont_hmin);
  if(ov_cont_k0 != NULL) free(ov_cont_k0);

  return(0);
}

int np_conditional_distribution_cvls_lp_stream(double *vector_scale_factor,
                                               double *cv);

int np_kernel_estimate_con_distribution_categorical_leave_one_out_ls_cv(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int64_t num_obs_train,
int64_t num_obs_eval,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
int cdfontrain,
double memfac,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_XY_unordered_train, 
double **matrix_XY_ordered_train, 
double **matrix_XY_continuous_train, 
double **matrix_Y_unordered_eval,
double **matrix_Y_ordered_eval,
double **matrix_Y_continuous_eval,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *cv){
  NP_GateOverrideCtx gate_x_ctx, gate_y_ctx;
  np_gate_ctx_clear(&gate_x_ctx);
  np_gate_ctx_clear(&gate_y_ctx);

  if(((BANDWIDTH_den == BW_FIXED) || (BANDWIDTH_den == BW_GEN_NN) || (BANDWIDTH_den == BW_ADAP_NN)) &&
     (int_ll_extern == LL_LP))
    return np_conditional_distribution_cvls_lp_stream(vector_scale_factor, cv);

  int indy;
  int64_t i,j,l,iwx,iwy;

  const int num_reg_tot = num_reg_continuous+num_reg_unordered+num_reg_ordered;
  const int num_var_tot = num_var_continuous+num_var_unordered+num_var_ordered;
  const int num_all_var = num_reg_tot + num_var_tot;

  const int num_all_cvar = num_reg_continuous + num_var_continuous;
  const int num_all_uvar = num_reg_unordered + num_var_unordered;

  size_t Nm = MIN((size_t)ceil(memfac*300000.0), (size_t)SIZE_MAX/10);

  int64_t N, num_obs_eval_alloc, num_obs_train_alloc, num_obs_wx_alloc, num_obs_wy_alloc;
  int64_t wx, wy, nwx, nwy;

  int * x_operator = NULL, * y_operator = NULL, * xy_operator = NULL;
  int gate_x_active = 0, gate_y_active = 0;
  int gate_x_all_large_fixed = 0;
  int *x_cont_ok = NULL, *x_disc_uno_ok = NULL, *x_disc_ord_ok = NULL;
  int *y_cont_ok = NULL, *y_disc_uno_ok = NULL, *y_disc_ord_ok = NULL;
  double *x_cont_hmin = NULL, *x_cont_k0 = NULL, *x_disc_uno_const = NULL, *x_disc_ord_const = NULL;
  double x_all_large_fixed_const = 1.0;
  double *y_cont_hmin = NULL, *y_cont_k0 = NULL, *y_disc_uno_const = NULL, *y_disc_ord_const = NULL;

  double vsfx[num_reg_tot];
  double vsfy[num_var_tot];
  double vsfxy[num_var_tot+num_reg_tot];
  double lambdax[MAX(1,num_reg_unordered+num_reg_ordered)];
  double lambday[MAX(1,num_var_unordered+num_var_ordered)];
  double xyj;

  double **matrix_wY_unordered_train;
  double **matrix_wY_ordered_train;
  double **matrix_wY_continuous_train;
  double **matrix_wX_unordered_train;
  double **matrix_wX_ordered_train;
  double **matrix_wX_continuous_train;
  double **matrix_wY_unordered_eval;
  double **matrix_wY_ordered_eval;
  double **matrix_wY_continuous_eval;

  double ** matrix_bandwidth_y, ** matrix_bandwidth_x;


  const int nbwmy = (BANDWIDTH_den == BW_FIXED) ? 1 : ((BANDWIDTH_den == BW_GEN_NN) ? num_obs_eval : num_obs_train);
  const int nbwmx = (BANDWIDTH_den == BW_FIXED) ? 1 : num_obs_train;

  int64_t js, je;

#ifdef MPI2
  int64_t stride_t = MAX((int64_t)ceil((double) num_obs_train / (double) iNum_Processors),1);
  int64_t stride_e = MAX((int64_t)ceil((double) num_obs_eval / (double) iNum_Processors),1);

  num_obs_train_alloc = stride_t*iNum_Processors;
  num_obs_eval_alloc = stride_e*iNum_Processors;

#else
  num_obs_train_alloc = num_obs_train;
  num_obs_eval_alloc = num_obs_eval;

  js = 0;
  je = num_obs_eval;
#endif


  // blocking algo calculations
  N = num_obs_train_alloc*num_obs_train_alloc + num_obs_eval_alloc*num_obs_train_alloc + num_obs_train_alloc + (num_obs_train_alloc + num_obs_eval_alloc)*num_obs_train_alloc;
  
  const int64_t sa = num_obs_train_alloc*num_obs_train_alloc*sizeof(double);
  const int64_t sb = num_obs_eval_alloc*num_obs_train_alloc*sizeof(double);

  if((N > Nm) || (sa > (((int64_t)1<<31)-1)) || (sb > (((int64_t)1<<31)-1))){
    const int64_t wy0 = (Nm - 2*num_obs_train_alloc - 1)/(2*num_obs_train_alloc);
    wy = ((wy0 > num_obs_eval_alloc) || (wy0 <= 0)) ? num_obs_eval_alloc : wy0;
    wx = (wy0 > 0) ? (Nm - 2*num_obs_train_alloc*wy)/(1 + 2*num_obs_train_alloc) : 1;
    nwx = num_obs_train_alloc/wx + (((num_obs_train_alloc % wx) > 0) ? 1 : 0);
    nwy = num_obs_eval_alloc/wy + (((num_obs_eval_alloc % wy) > 0) ? 1 : 0);
  } else {
    wx = num_obs_train_alloc;
    wy = num_obs_eval_alloc;
    nwx = 1;
    nwy = 1;
  }

#ifdef MPI2
  int64_t stride_wx = MAX((int64_t)ceil((double)wx / (double) iNum_Processors),1);
  int64_t stride_wy = MAX((int64_t)ceil((double)wy / (double) iNum_Processors),1);

  num_obs_wx_alloc = stride_wx*iNum_Processors;
  num_obs_wy_alloc = stride_wy*iNum_Processors;

  js = stride_wy*my_rank;
  je = MIN(wy, js + stride_wy);
#else
  num_obs_wx_alloc = wx;
  num_obs_wy_alloc = wy;
#endif

  int * kernel_cx = NULL, * kernel_ux = NULL, * kernel_ox = NULL;
  int * kernel_cy = NULL, * kernel_uy = NULL, * kernel_oy = NULL;

  // first the x-kernels
  kernel_cx = (int *)malloc(sizeof(int)*num_reg_continuous);

  for(i = 0; i < num_reg_continuous; i++)
    kernel_cx[i] = KERNEL_reg;

  kernel_ux = (int *)malloc(sizeof(int)*num_reg_unordered);

  for(i = 0; i < num_reg_unordered; i++)
    kernel_ux[i] = KERNEL_unordered_reg;

  kernel_ox = (int *)malloc(sizeof(int)*num_reg_ordered);

  for(i = 0; i < num_reg_ordered; i++)
    kernel_ox[i] = KERNEL_ordered_reg;

  // then the y-kernels
  kernel_cy = (int *)malloc(sizeof(int)*num_var_continuous);

  for(i = 0; i < num_var_continuous; i++)
    kernel_cy[i] = KERNEL_den;

  kernel_uy = (int *)malloc(sizeof(int)*num_var_unordered);

  for(i = 0; i < num_var_unordered; i++)
    kernel_uy[i] = KERNEL_unordered_den;

  kernel_oy = (int *)malloc(sizeof(int)*num_var_ordered);

  for(i = 0; i < num_var_ordered; i++)
    kernel_oy[i] = KERNEL_ordered_den;

  // allocate some pointers
  matrix_wX_continuous_train = (double **)malloc(num_reg_continuous*sizeof(double *));
  matrix_wX_unordered_train = (double **)malloc(num_reg_unordered*sizeof(double *));
  matrix_wX_ordered_train = (double **)malloc(num_reg_ordered*sizeof(double *));

  matrix_wY_continuous_train = (double **)malloc(num_var_continuous*sizeof(double *));
  matrix_wY_unordered_train = (double **)malloc(num_var_unordered*sizeof(double *));
  matrix_wY_ordered_train = (double **)malloc(num_var_ordered*sizeof(double *));

  matrix_wY_continuous_eval = (double **)malloc(num_var_continuous*sizeof(double *));
  matrix_wY_unordered_eval = (double **)malloc(num_var_unordered*sizeof(double *));
  matrix_wY_ordered_eval = (double **)malloc(num_var_ordered*sizeof(double *));


  double * mean = (double *)malloc(num_obs_wx_alloc*sizeof(double));

  
  if(mean == NULL)
    error("failed to allocate mean");

  np_splitxy_vsf_mcv_nc(num_var_unordered, num_var_ordered, num_var_continuous,
                        num_reg_unordered, num_reg_ordered, num_reg_continuous,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        vsfx,
                        vsfy,
                        vsfxy,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);
  

  x_operator = (int *)malloc(sizeof(int)*(num_reg_continuous+num_reg_unordered+num_reg_ordered));

  if(x_operator == NULL)
    error("failed to allocate x_operator");

  for(i = 0; i < (num_reg_continuous+num_reg_unordered+num_reg_ordered); i++)
    x_operator[i] = OP_NORMAL;

  y_operator = (int *)malloc(sizeof(int)*(num_var_continuous+num_var_unordered+num_var_ordered));

  if(y_operator == NULL)
    error("failed to allocate y_operator");

  for(i = 0; i < (num_var_continuous+num_var_unordered+num_var_ordered); i++)
    y_operator[i] = OP_INTEGRAL;

  xy_operator = (int *)malloc(sizeof(int)*num_all_var);

  if(xy_operator == NULL)
    error("failed to allocate xy_operator");

  for(i = 0; i < num_reg_continuous; i++)
    xy_operator[i] = OP_NORMAL;

  for(i = num_reg_continuous; i < num_all_cvar; i++)
    xy_operator[i] = OP_INTEGRAL;

  // no cdf for unordered
  for(i = num_all_cvar; i < (num_all_cvar+num_all_uvar+num_reg_ordered); i++)
    xy_operator[i] = OP_NORMAL;

  for(i = num_all_cvar+num_all_uvar+num_reg_ordered; i < num_all_var; i++)
    xy_operator[i] = OP_INTEGRAL;


  matrix_bandwidth_x = alloc_matd(nbwmx, num_reg_continuous);
  matrix_bandwidth_y = alloc_matd(nbwmy, num_var_continuous);

  kernel_bandwidth_mean(KERNEL_den,
                        BANDWIDTH_den,
                        num_obs_train,
                        nbwmx,
                        0,
                        0,
                        0,
                        num_reg_continuous,
                        num_reg_unordered,
                        num_reg_ordered,
                        0, // do not suppress_parallel
                        vsfx,
                        NULL,
                        NULL,
                        matrix_X_continuous_train,
                        matrix_X_continuous_train,
                        NULL,					 // Not used 
                        matrix_bandwidth_x,
                        lambdax);

  kernel_bandwidth_mean(KERNEL_reg,
                        BANDWIDTH_den,
                        num_obs_train,
                        nbwmy,
                        0,
                        0,
                        0,
                        num_var_continuous,
                        num_var_unordered,
                        num_var_ordered,
                        0, // do not suppress_parallel
                        vsfy,
                        NULL,
                        NULL,
                        matrix_Y_continuous_train,
                        matrix_Y_continuous_eval,
                        NULL,					 // Not used 
                        matrix_bandwidth_y,
                        lambday);

  if(num_reg_continuous > 0 || num_reg_unordered > 0 || num_reg_ordered > 0){
    int ok_all = 1;

    if(num_reg_continuous > 0){
      x_cont_ok = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
      x_cont_hmin = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
      x_cont_k0 = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
      ok_all = (x_cont_ok != NULL) && (x_cont_hmin != NULL) && (x_cont_k0 != NULL);
      if(ok_all){
        const double rel_tol = np_cont_largeh_rel_tol();
        for(i = 0; i < num_reg_continuous; i++){
          const int kern = kernel_cx[i];
          double xmin = DBL_MAX, xmax = -DBL_MAX;
          x_cont_ok[i] = 0; x_cont_hmin[i] = DBL_MAX; x_cont_k0[i] = 0.0;
          if(!np_cont_largeh_kernel_supported(kern)) continue;
          for(j = 0; j < num_obs_train; j++){
            const double v = matrix_X_continuous_train[i][j];
            if(!isfinite(v)) continue;
            xmin = MIN(xmin, v); xmax = MAX(xmax, v);
          }
          if(xmax >= xmin){
            const double utol = np_cont_largeh_utol(kern, rel_tol);
            if(utol > 0.0 && isfinite(utol)){
              x_cont_ok[i] = 1;
              x_cont_hmin[i] = (xmax - xmin)/utol;
              x_cont_k0[i] = np_cont_largeh_k0(kern);
            }
          }
        }
      }
    }

    if(ok_all && num_reg_unordered > 0){
      x_disc_uno_ok = (int *)calloc((size_t)num_reg_unordered, sizeof(int));
      x_disc_uno_const = (double *)malloc((size_t)num_reg_unordered*sizeof(double));
      ok_all = (x_disc_uno_ok != NULL) && (x_disc_uno_const != NULL);
      if(ok_all){
        double (* const ukf[])(int, double, int) = {
          np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
          np_score_uaa, np_score_unli_racine
        };
        const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
        for(i = 0; i < num_reg_unordered; i++){
          const int ku = kernel_ux[i];
          const int ncat = (num_categories_extern_X != NULL) ? num_categories_extern_X[i] : 0;
          const double lam = lambdax[i];
          x_disc_uno_ok[i] = 0; x_disc_uno_const[i] = 0.0;
          if(ku < 0 || ku >= nuk) continue;
          if(!np_disc_near_upper(ku, lam, ncat)) continue;
          {
            const double ks = ukf[ku](1, lam, ncat);
            const double kd = ukf[ku](0, lam, ncat);
            if(np_disc_near_const_kernel(ks, kd)){
              x_disc_uno_ok[i] = 1;
              x_disc_uno_const[i] = 0.5*(ks + kd);
            }
          }
        }
      }
    }

    if(ok_all && num_reg_ordered > 0){
      x_disc_ord_ok = (int *)calloc((size_t)num_reg_ordered, sizeof(int));
      x_disc_ord_const = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
      ok_all = (x_disc_ord_ok != NULL) && (x_disc_ord_const != NULL);
      if(ok_all){
        double (* const okf[])(double, double, double, double, double) = {
          np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
        };
        const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
        for(i = 0; i < num_reg_ordered; i++){
          const int oi = i + num_reg_unordered;
          const int ko = kernel_ox[i];
          const int ncat = (num_categories_extern_X != NULL) ? num_categories_extern_X[oi] : 0;
          const double lam = lambdax[oi];
          x_disc_ord_ok[i] = 0; x_disc_ord_const[i] = 0.0;
          if(ko < 0 || ko >= nok) continue;
          if(ncat <= 0 || matrix_categorical_vals_extern_X == NULL) continue;
          if(!np_disc_ordered_near_upper(ko, lam)) continue;
          {
            const double cl = matrix_categorical_vals_extern_X[oi][0];
            const double ch = matrix_categorical_vals_extern_X[oi][ncat - 1];
            const double k0 = okf[ko](cl, cl, lam, cl, ch);
            const double k1 = okf[ko](cl, ch, lam, cl, ch);
            if(np_disc_near_const_kernel(k0, k1)){
              x_disc_ord_ok[i] = 1;
              x_disc_ord_const[i] = 0.5*(k0 + k1);
            }
          }
        }
      }
    }

    if(ok_all){
      gate_x_active = 1;
    }
  }

  if((BANDWIDTH_den == BW_FIXED) && (num_reg_tot > 0)){
    int ok_all_large = 1;
    double kconst = 1.0;

    for(l = 0; l < num_reg_continuous; l++){
      const double bw = matrix_bandwidth_x[l][0];
      if((x_cont_ok == NULL) || (!x_cont_ok[l]) || (!isfinite(bw)) ||
         (bw <= 0.0) || (bw < x_cont_hmin[l])){
        ok_all_large = 0;
        break;
      }
      kconst *= x_cont_k0[l]/bw;
    }

    if(ok_all_large){
      double (* const ukf[])(int, double, int) = {
        np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
        np_score_uaa, np_score_unli_racine
      };
      const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
      for(l = 0; l < num_reg_unordered; l++){
        const int ku = kernel_ux[l];
        const int ncat = (num_categories_extern_X != NULL) ? num_categories_extern_X[l] : 0;
        const double lam = lambdax[l];
        if((ku < 0) || (ku >= nuk) || (!isfinite(lam)) || (!np_disc_near_upper(ku, lam, ncat))){
          ok_all_large = 0;
          break;
        }
        {
          const double ks = ukf[ku](1, lam, ncat);
          const double kd = ukf[ku](0, lam, ncat);
          if(!np_disc_near_const_kernel(ks, kd)){
            ok_all_large = 0;
            break;
          }
          kconst *= 0.5*(ks + kd);
        }
      }
    }

    if(ok_all_large){
      double (* const okf[])(double, double, double, double, double) = {
        np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
      };
      const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
      for(l = 0; l < num_reg_ordered; l++){
        const int oi = num_reg_unordered + l;
        const int ko = kernel_ox[l];
        const int ncat = (num_categories_extern_X != NULL) ? num_categories_extern_X[oi] : 0;
        const double lam = lambdax[oi];
        if((ko < 0) || (ko >= nok) || (!isfinite(lam)) || (!np_disc_ordered_near_upper(ko, lam)) ||
           (ncat <= 0) || (matrix_categorical_vals_extern_X == NULL)){
          ok_all_large = 0;
          break;
        }
        {
          const double cl = matrix_categorical_vals_extern_X[oi][0];
          const double ch = matrix_categorical_vals_extern_X[oi][ncat - 1];
          const double k0 = okf[ko](cl, cl, lam, cl, ch);
          const double k1 = okf[ko](cl, ch, lam, cl, ch);
          if(!np_disc_near_const_kernel(k0, k1)){
            ok_all_large = 0;
            break;
          }
          kconst *= 0.5*(k0 + k1);
        }
      }
    }

    if(ok_all_large && isfinite(kconst) && (kconst > 0.0)){
      gate_x_all_large_fixed = 1;
      x_all_large_fixed_const = kconst;
    }
  }
  if(gate_x_all_large_fixed)
    np_fastcv_alllarge_hits++;

  if(num_var_continuous > 0 || num_var_unordered > 0 || num_var_ordered > 0){
    int ok_all = 1;

    if(num_var_continuous > 0){
      y_cont_ok = (int *)calloc((size_t)num_var_continuous, sizeof(int));
      y_cont_hmin = (double *)malloc((size_t)num_var_continuous*sizeof(double));
      y_cont_k0 = (double *)malloc((size_t)num_var_continuous*sizeof(double));
      ok_all = (y_cont_ok != NULL) && (y_cont_hmin != NULL) && (y_cont_k0 != NULL);
      if(ok_all){
        const double rel_tol = np_cont_largeh_rel_tol();
        for(i = 0; i < num_var_continuous; i++){
          const int kern = kernel_cy[i];
          double xmin = DBL_MAX, xmax = -DBL_MAX;
          y_cont_ok[i] = 0; y_cont_hmin[i] = DBL_MAX; y_cont_k0[i] = 0.0;
          if(!np_cont_largeh_kernel_supported(kern)) continue;
          for(j = 0; j < num_obs_train; j++){
            const double v = matrix_Y_continuous_train[i][j];
            if(!isfinite(v)) continue;
            xmin = MIN(xmin, v); xmax = MAX(xmax, v);
          }
          for(j = 0; j < num_obs_eval; j++){
            const double v = matrix_Y_continuous_eval[i][j];
            if(!isfinite(v)) continue;
            xmin = MIN(xmin, v); xmax = MAX(xmax, v);
          }
          if(xmax >= xmin){
            const double utol = np_cont_largeh_utol(kern, rel_tol);
            if(utol > 0.0 && isfinite(utol)){
              y_cont_ok[i] = 1;
              y_cont_hmin[i] = (xmax - xmin)/utol;
              y_cont_k0[i] = np_cont_largeh_k0(kern);
            }
          }
        }
      }
    }

    if(ok_all && num_var_unordered > 0){
      y_disc_uno_ok = (int *)calloc((size_t)num_var_unordered, sizeof(int));
      y_disc_uno_const = (double *)malloc((size_t)num_var_unordered*sizeof(double));
      ok_all = (y_disc_uno_ok != NULL) && (y_disc_uno_const != NULL);
      if(ok_all){
        double (* const ukf[])(int, double, int) = {
          np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
          np_score_uaa, np_score_unli_racine
        };
        const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
        for(i = 0; i < num_var_unordered; i++){
          const int ku = kernel_uy[i];
          const int ncat = (num_categories_extern_Y != NULL) ? num_categories_extern_Y[i] : 0;
          const double lam = lambday[i];
          y_disc_uno_ok[i] = 0; y_disc_uno_const[i] = 0.0;
          if(ku < 0 || ku >= nuk) continue;
          if(!np_disc_near_upper(ku, lam, ncat)) continue;
          {
            const double ks = ukf[ku](1, lam, ncat);
            const double kd = ukf[ku](0, lam, ncat);
            if(np_disc_near_const_kernel(ks, kd)){
              y_disc_uno_ok[i] = 1;
              y_disc_uno_const[i] = 0.5*(ks + kd);
            }
          }
        }
      }
    }

    if(ok_all && num_var_ordered > 0){
      y_disc_ord_ok = (int *)calloc((size_t)num_var_ordered, sizeof(int));
      y_disc_ord_const = (double *)malloc((size_t)num_var_ordered*sizeof(double));
      ok_all = (y_disc_ord_ok != NULL) && (y_disc_ord_const != NULL);
      if(ok_all){
        double (* const okf[])(double, double, double, double, double) = {
          np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
        };
        const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
        for(i = 0; i < num_var_ordered; i++){
          const int oi = i + num_var_unordered;
          const int ko = kernel_oy[i];
          const int ncat = (num_categories_extern_Y != NULL) ? num_categories_extern_Y[oi] : 0;
          const double lam = lambday[oi];
          y_disc_ord_ok[i] = 0; y_disc_ord_const[i] = 0.0;
          if(ko < 0 || ko >= nok) continue;
          if(ncat <= 0 || matrix_categorical_vals_extern_Y == NULL) continue;
          if(!np_disc_ordered_near_upper(ko, lam)) continue;
          {
            const double cl = matrix_categorical_vals_extern_Y[oi][0];
            const double ch = matrix_categorical_vals_extern_Y[oi][ncat - 1];
            const double k0 = okf[ko](cl, cl, lam, cl, ch);
            const double k1 = okf[ko](cl, ch, lam, cl, ch);
            if(np_disc_near_const_kernel(k0, k1)){
              y_disc_ord_ok[i] = 1;
              y_disc_ord_const[i] = 0.5*(k0 + k1);
            }
          }
        }
      }
    }

    if(ok_all){
      gate_y_active = 1;
    }
  }

  *cv = 0;

  if(!int_TREE_XY || gate_x_all_large_fixed){
    double * kwx = NULL;

    double * kwy = (double *)malloc(num_obs_train_alloc*num_obs_wy_alloc*sizeof(double));

    if(kwy == NULL)
      error("failed to allocate kwy, try reducing num_obs_eval, tried to allocate: %" PRIi64 "bytes\n", num_obs_train_alloc*num_obs_wy_alloc*sizeof(double));

    double *kwy_row_sum = NULL;

    if(!gate_x_all_large_fixed){
      kwx = (double *)malloc(num_obs_wx_alloc*num_obs_train_alloc*sizeof(double));
      if(kwx == NULL)
        error("failed to allocate kwx, tried to allocate: %" PRIi64 "bytes\n", num_obs_wx_alloc*num_obs_train_alloc*sizeof(double));
    } else {
      kwy_row_sum = (double *)malloc(num_obs_wy_alloc*sizeof(double));
      if(kwy_row_sum == NULL)
        error("failed to allocate kwy_row_sum");
    }
    
    for(iwx = 0; iwx < nwx; iwx++){
      const int64_t wxo = iwx*wx;
      const int64_t dwx = (iwx != (nwx - 1)) ? wx : num_obs_train - (nwx - 1)*wx;

      for(l = 0; l < num_reg_continuous; l++)
        matrix_wX_continuous_train[l] = matrix_X_continuous_train[l] + wxo;

      for(l = 0; l < num_reg_unordered; l++)
        matrix_wX_unordered_train[l] = matrix_X_unordered_train[l] + wxo;

      for(l = 0; l < num_reg_ordered; l++)
        matrix_wX_ordered_train[l] = matrix_X_ordered_train[l] + wxo;


      if(!gate_x_all_large_fixed){
        if(gate_x_active){
          np_gate_ctx_set(&gate_x_ctx,
                          num_reg_continuous,
                          num_reg_unordered,
                          num_reg_ordered,
                          kernel_cx,
                          kernel_ux,
                          kernel_ox,
                          x_operator,
                          x_cont_ok,
                          x_cont_hmin,
                          x_cont_k0,
                          x_disc_uno_ok,
                          x_disc_uno_const,
                          x_disc_ord_ok,
                          x_disc_ord_const);
        } else {
          np_gate_ctx_clear(&gate_x_ctx);
        }
        np_activate_bounds_x();
        kernel_weighted_sum_np_ctx(kernel_cx,
                               kernel_ux,
                               kernel_ox,
                               BANDWIDTH_den,
                               num_obs_train,
                               dwx,
                               num_reg_unordered,
                               num_reg_ordered,
                               num_reg_continuous,
                               0, // (do not) compute the leave-one-out marginals
                               0, // '' offset
                               1, // kernel power
                               0, // bandwidth_divide
                               0, // '' weights
                               0, // symmetric
                               0, // gather_scatter sum
                               0, // drop train
                               0, // drop which train
                               x_operator,
                               OP_NOOP, // no permutations
                               0, // no score
                               0, // no ocg
                               NULL, // explicit bpso
                               0, // don't explicity suppress parallel
                               0, // ncol y 
                               0, // ncol w
                               int_TREE_X, // do tree
                               0, // do partial tree 
                               kdt_extern_X, // which tree
                               NULL, NULL, NULL, // partial tree data
                               matrix_X_unordered_train, 
                               matrix_X_ordered_train,
                               matrix_X_continuous_train,
                               matrix_wX_unordered_train,
                               matrix_wX_ordered_train,
                               matrix_wX_continuous_train,
                               NULL, // matrix y
                               NULL, // matrix w
                               NULL, // sgn
                               vsfx,
                               1,matrix_bandwidth_x,matrix_bandwidth_x,lambdax,
                               num_categories_extern_X,
                               matrix_categorical_vals_extern_X,
                               NULL, // moo
                               mean,
                               NULL, // no permutations
                               kwx,
                               &gate_x_ctx);
      }

      for(iwy = 0; iwy < nwy; iwy++){

        const int64_t wyo = iwy*wy;
        const int64_t dwy = (iwy != (nwy - 1)) ? wy : num_obs_eval - (nwy - 1)*wy;

        for(l = 0; l < num_var_continuous; l++)
          matrix_wY_continuous_eval[l] = matrix_Y_continuous_eval[l] + wyo;

        for(l = 0; l < num_var_unordered; l++)
          matrix_wY_unordered_eval[l] = matrix_Y_unordered_eval[l] + wyo;

        for(l = 0; l < num_var_ordered; l++)
          matrix_wY_ordered_eval[l] = matrix_Y_ordered_eval[l] + wyo;

        // compute y weights first
        if(gate_y_active){
          np_gate_ctx_set(&gate_y_ctx,
                          num_var_continuous,
                          num_var_unordered,
                          num_var_ordered,
                          kernel_cy,
                          kernel_uy,
                          kernel_oy,
                          y_operator,
                          y_cont_ok,
                          y_cont_hmin,
                          y_cont_k0,
                          y_disc_uno_ok,
                          y_disc_uno_const,
                          y_disc_ord_ok,
                          y_disc_ord_const);
        } else {
          np_gate_ctx_clear(&gate_y_ctx);
        }
        np_activate_bounds_y();
        kernel_weighted_sum_np_ctx(kernel_cy,
                               kernel_uy,
                               kernel_oy,
                               BANDWIDTH_den,
                               num_obs_train,
                               dwy,
                               num_var_unordered,
                               num_var_ordered,
                               num_var_continuous,
                               0,
                               0,
                               1,
                               0,
                               0, 
                               0,
                               0,
                               0,
                               0,
                               y_operator,
                               OP_NOOP, // no permutations
                               0, // no score
                               0, // no ocg
                               NULL,
                               0, // don't explicity suppress parallel
                               0,
                               0,
                               int_TREE_Y,
                               0,
                               kdt_extern_Y,
                               NULL, NULL, NULL,
                               matrix_Y_unordered_train,
                               matrix_Y_ordered_train,
                               matrix_Y_continuous_train,
                               matrix_wY_unordered_eval,
                               matrix_wY_ordered_eval,
                               matrix_wY_continuous_eval,
                               NULL,
                               NULL,
                               NULL,
                               vsfy,
                               1,matrix_bandwidth_y,matrix_bandwidth_y,lambday,
                               num_categories_extern_Y,
                               matrix_categorical_vals_extern_Y,
                               NULL,
                               NULL,
                               NULL, // no permutations
                               kwy,
                               &gate_y_ctx);

        const int64_t je_dwy = MIN(je,dwy);

        if(gate_x_all_large_fixed){
          const double inv_nmo = 1.0/(((double)(num_obs_train - 1)) + DBL_MIN);
          for(j = (wyo + js); j < (wyo + je_dwy); j++){
            const int64_t jo = j - wyo;
            double rows = 0.0;
            for(l = 0; l < num_obs_train; l++)
              rows += kwy[jo*num_obs_train+l];
            kwy_row_sum[jo] = rows;
          }

          for(i = wxo; i < (wxo + dwx); i++){
            for(j = (wyo + js); j < (wyo + je_dwy); j++){
              if(cdfontrain && (j == i)) continue;
              const int64_t jo = j - wyo;
              indy = 1;
              for(l = 0; l < num_var_ordered; l++){
                indy *= (matrix_Y_ordered_train[l][i] <= matrix_Y_ordered_eval[l][j]);
              }
              for(l = 0; l < num_var_continuous; l++){
                indy *= (matrix_Y_continuous_train[l][i] <= matrix_Y_continuous_eval[l][j]);
              }
              xyj = kwy_row_sum[jo] - kwy[jo*num_obs_train+i];
              {
                const double tvd = (indy - xyj*inv_nmo);
                *cv += tvd*tvd;
              }
            }
          }
        } else {
          for(i = wxo; i < (wxo + dwx); i++){
            const int64_t io = i - wxo;

            for(j = (wyo + js); j < (wyo + je_dwy); j++){
              if(cdfontrain && (j == i)) continue;
              const int64_t jo = j - wyo;
              indy = 1;
              for(l = 0; l < num_var_ordered; l++){
                indy *= (matrix_Y_ordered_train[l][i] <= matrix_Y_ordered_eval[l][j]);
              }
              for(l = 0; l < num_var_continuous; l++){
                indy *= (matrix_Y_continuous_train[l][i] <= matrix_Y_continuous_eval[l][j]);
              }
              xyj = 0.0;

              if(BANDWIDTH_den != BW_ADAP_NN){
                // leave-one-out joint density

                for(l = 0; l < num_obs_train; l++)
                  xyj += kwy[jo*num_obs_train+l]*kwx[io*num_obs_train+l];
                xyj -= kwy[jo*num_obs_train+i]*kwx[io*num_obs_train+i];

                const double tvd = (indy - xyj/(mean[io] - kwx[io*num_obs_train+i] + DBL_MIN));
                *cv += tvd*tvd;
              } else {
                // leave-one-out joint density

                for(l = 0; l < num_obs_train; l++)
                  xyj += kwy[l*dwy+jo]*kwx[l*dwx+io];
                xyj -= kwy[i*dwy+jo]*kwx[i*dwx+io];

                const double tvd = (indy - xyj/(mean[io] - kwx[i*dwx + io] + DBL_MIN));
                *cv += tvd*tvd;
              }
            }
          }
        }
      }
    }

#ifdef MPI2
    MPI_Allreduce(MPI_IN_PLACE, cv, 1, MPI_DOUBLE, MPI_SUM, comm[1]);
#endif
    *cv /= (double) num_obs_train*num_obs_eval;

    free(kwx);
    free(kwy);
    free(kwy_row_sum);
  } else {
    NL nls = {.node = NULL, .n = 0, .nalloc = 0};
    NL nlps = {.node = NULL, .n = 0, .nalloc = 0};

    XL xl = {.istart = NULL, .nlev = NULL, .n = 0, .nalloc = 0};

    int icx[MAX(num_reg_continuous,1)], icy[MAX(num_var_continuous,1)];

    double bb[MAX(1,2*num_all_cvar)];

    int KERNEL_XY[MAX(1,num_all_cvar)], m;

    for(l = 0; l < num_reg_continuous; l++)
      KERNEL_XY[l] = KERNEL_reg + OP_CFUN_OFFSETS[xy_operator[l]];

    for(l = num_reg_continuous; l < num_all_cvar; l++)
      KERNEL_XY[l] = KERNEL_den + OP_CFUN_OFFSETS[xy_operator[l]];

    double * kwx = (double *)malloc(num_obs_wx_alloc*num_obs_train_alloc*sizeof(double));

    if(kwx == NULL)
      error("failed to allocate kwx, tried to allocate: %" PRIi64 "bytes\n", num_obs_wx_alloc*num_obs_train_alloc*sizeof(double));

    double * kwy = (double *)malloc(num_obs_train_alloc*num_obs_wy_alloc*sizeof(double));

    if(kwy == NULL)
      error("failed to allocate kwy, try reducing num_obs_eval, tried to allocate: %" PRIi64 "bytes\n", num_obs_train_alloc*num_obs_wy_alloc*sizeof(double));

    nls.node = (int *)malloc(sizeof(int));
    nls.nalloc = 1;

    nls.node[0] = 0;
    nls.n = 1;

    for(l = 0; l < num_reg_continuous; l++)
      icx[l] = l;

    for(l = num_reg_continuous; l < num_all_cvar; l++)
      icy[l-num_reg_continuous] = l;

    for(iwx = 0; iwx < nwx; iwx++){
      const int64_t wxo = iwx*wx;
      const int64_t dwx = (iwx != (nwx - 1)) ? wx : num_obs_train - (nwx - 1)*wx;

      for(l = 0; l < num_reg_continuous; l++)
        matrix_wX_continuous_train[l] = matrix_XY_continuous_train[l] + wxo;

      for(l = 0; l < num_reg_unordered; l++)
        matrix_wX_unordered_train[l] = matrix_XY_unordered_train[l] + wxo;

      for(l = 0; l < num_reg_ordered; l++)
        matrix_wX_ordered_train[l] = matrix_XY_ordered_train[l] + wxo;


      if(gate_x_all_large_fixed){
        const double mconst = ((double)num_obs_train)*x_all_large_fixed_const;
        for(i = 0; i < dwx; i++)
          mean[i] = mconst;
        for(i = 0; i < (dwx*num_obs_train); i++)
          kwx[i] = x_all_large_fixed_const;
      } else {
        if(gate_x_active){
          np_gate_ctx_set(&gate_x_ctx,
                          num_reg_continuous,
                          num_reg_unordered,
                          num_reg_ordered,
                          kernel_cx,
                          kernel_ux,
                          kernel_ox,
                          x_operator,
                          x_cont_ok,
                          x_cont_hmin,
                          x_cont_k0,
                          x_disc_uno_ok,
                          x_disc_uno_const,
                          x_disc_ord_ok,
                          x_disc_ord_const);
        } else {
          np_gate_ctx_clear(&gate_x_ctx);
        }
        np_activate_bounds_x();
        kernel_weighted_sum_np_ctx(kernel_cx,
                               kernel_ux,
                               kernel_ox,
                               BANDWIDTH_den,
                               num_obs_train,
                               dwx,
                               num_reg_unordered,
                               num_reg_ordered,
                               num_reg_continuous,
                               0, // compute the leave-one-out marginals
                               0,
                               1,
                               0,
                               0, 
                               0,
                               0,
                               0,
                               0,
                               x_operator,
                               OP_NOOP, // no permutations
                               0, // no score
                               0, // no ocg
                               NULL,
                               0, // don't explicity suppress parallel
                               0,
                               0,
                               int_TREE_XY,
                               1,
                               kdt_extern_XY,
                               &nls, icx, NULL,
                               matrix_XY_unordered_train,
                               matrix_XY_ordered_train,
                               matrix_XY_continuous_train,
                               matrix_wX_unordered_train,
                               matrix_wX_ordered_train,
                               matrix_wX_continuous_train,
                               NULL,
                               NULL,
                               NULL,
                               vsfx,
                               1,matrix_bandwidth_x,matrix_bandwidth_x,lambdax,
                               num_categories_extern_X,
                               matrix_categorical_vals_extern_X,
                               NULL,
                               mean,
                               NULL, // no permutations
                               kwx,
                               &gate_x_ctx);
      }

      for(iwy = 0; iwy < nwy; iwy++){
        const int64_t wyo = iwy*wy;
        const int64_t dwy = (iwy != (nwy - 1)) ? wy : num_obs_eval - (nwy - 1)*wy;

        for(l = 0; l < num_var_continuous; l++)
          matrix_wY_continuous_eval[l] = matrix_Y_continuous_eval[l] + wyo;

        for(l = 0; l < num_var_unordered; l++)
          matrix_wY_unordered_eval[l] = matrix_Y_unordered_eval[l] + wyo;

        for(l = 0; l < num_var_ordered; l++)
          matrix_wY_ordered_eval[l] = matrix_Y_ordered_eval[l] + wyo;

        // compute y weights first
        if(gate_y_active){
          np_gate_ctx_set(&gate_y_ctx,
                          num_var_continuous,
                          num_var_unordered,
                          num_var_ordered,
                          kernel_cy,
                          kernel_uy,
                          kernel_oy,
                          y_operator,
                          y_cont_ok,
                          y_cont_hmin,
                          y_cont_k0,
                          y_disc_uno_ok,
                          y_disc_uno_const,
                          y_disc_ord_ok,
                          y_disc_ord_const);
        } else {
          np_gate_ctx_clear(&gate_y_ctx);
        }
        np_activate_bounds_y();
        kernel_weighted_sum_np_ctx(kernel_cy,
                               kernel_uy,
                               kernel_oy,
                               BANDWIDTH_den,
                               num_obs_train,
                               dwy,
                               num_var_unordered,
                               num_var_ordered,
                               num_var_continuous,
                               0,
                               0,
                               1,
                               0,
                               0, 
                               0,
                               0,
                               0,
                               0,
                               y_operator,
                               OP_NOOP, // no permutations
                               0, // no score
                               0, // no ocg
                               NULL,
                               0, // don't explicity suppress parallel
                               0,
                               0,
                               int_TREE_XY,
                               1,
                               kdt_extern_XY,
                               &nls, icy, NULL,
                               matrix_XY_unordered_train + num_reg_unordered,
                               matrix_XY_ordered_train + num_reg_ordered,
                               matrix_XY_continuous_train + num_reg_continuous,
                               matrix_wY_unordered_eval,
                               matrix_wY_ordered_eval,
                               matrix_wY_continuous_eval,
                               NULL,
                               NULL,
                               NULL,
                               vsfy,
                               1,matrix_bandwidth_y,matrix_bandwidth_y,lambday,
                               num_categories_extern_Y,
                               matrix_categorical_vals_extern_Y,
                               NULL,
                               NULL,
                               NULL, // no permutations
                               kwy,
                               &gate_y_ctx);

        const int64_t je_dwy = MIN(je,dwy);

        for(i = wxo; i < (wxo + dwx); i++){
          const int64_t io = i - wxo;

          const int ixbw = (BANDWIDTH_den == BW_FIXED) ? 0 : i;

          for(l = 0; l < num_reg_continuous; l++){
            bb[2*l] = -cksup[KERNEL_XY[l]][1];
            bb[2*l+1] = -cksup[KERNEL_XY[l]][0];

            bb[2*l] = (fabs(bb[2*l]) == DBL_MAX) ? bb[2*l] : (matrix_XY_continuous_train[l][i] + bb[2*l]*matrix_bandwidth_x[l][ixbw]);
            bb[2*l+1] = (fabs(bb[2*l+1]) == DBL_MAX) ? bb[2*l+1] : (matrix_XY_continuous_train[l][i] + bb[2*l+1]*matrix_bandwidth_x[l][ixbw]);
          }

          const double mi = mean[io] - kwx[io*num_obs_train + i];

          // search for the point (x_i,y_j) in the xy tree
          // reset the interaction node list
          nlps.n = 0;

          boxSearchNLPartial(kdt_extern_XY, &nls, bb, &nlps, NULL, icx, num_reg_continuous);

          for(j = (wyo + js); j < (wyo + je_dwy); j++){
            if(cdfontrain && (j == i)) continue;

            const int64_t jo = j - wyo;

            const int iybw = (BANDWIDTH_den == BW_FIXED) ? 0 : j;

            indy = 1;
            for(l = 0; l < num_var_ordered; l++){
              indy *= (matrix_XY_ordered_train[l+num_reg_ordered][i] <= matrix_Y_ordered_eval[l][j]);
            }
            for(l = 0; l < num_var_continuous; l++){
              indy *= (matrix_XY_continuous_train[l+num_reg_continuous][i] <= matrix_Y_continuous_eval[l][j]);
            }

            // search for the point (x_i,y_j) in the xy tree
            // reset the interaction node list
            xl.n = 0;

            for(l = num_reg_continuous; l < num_all_cvar; l++){
              bb[2*l] = -cksup[KERNEL_XY[l]][1];
              bb[2*l+1] = -cksup[KERNEL_XY[l]][0];

              bb[2*l] = (fabs(bb[2*l]) == DBL_MAX) ? bb[2*l] : (matrix_Y_continuous_eval[l-num_reg_continuous][j] + bb[2*l]*matrix_bandwidth_y[l-num_reg_continuous][iybw]);
              bb[2*l+1] = (fabs(bb[2*l+1]) == DBL_MAX) ? bb[2*l+1] : (matrix_Y_continuous_eval[l-num_reg_continuous][j] + bb[2*l+1]*matrix_bandwidth_y[l-num_reg_continuous][iybw]);
            }

            boxSearchNLPartial(kdt_extern_XY, &nlps, bb, NULL, &xl, icy, num_var_continuous);

            xyj = 0.0;

            for (m = 0; m < xl.n; m++){
              for (l = xl.istart[m]; l < (xl.istart[m] + xl.nlev[m]); l++){
                xyj += kwy[jo*num_obs_train+l]*kwx[io*num_obs_train+l];
              }
            }
            xyj -= kwy[jo*num_obs_train+i]*kwx[io*num_obs_train+i];
            const double tvd = (indy - xyj/(mi + DBL_MIN));
            *cv += tvd*tvd;
          }
        }
      }
    }

#ifdef MPI2
    MPI_Allreduce(MPI_IN_PLACE, cv, 1, MPI_DOUBLE, MPI_SUM, comm[1]);
#endif
    *cv /= (double) num_obs_train*num_obs_eval;

    free(kwx);
    free(kwy);

    clean_nl(&nls);
    clean_nl(&nlps);

    clean_xl(&xl);
  }

  free(x_operator);
  free(y_operator);
  free(xy_operator);

  free(x_cont_ok);
  free(x_disc_uno_ok);
  free(x_disc_ord_ok);
  free(x_cont_hmin);
  free(x_cont_k0);
  free(x_disc_uno_const);
  free(x_disc_ord_const);

  free(y_cont_ok);
  free(y_disc_uno_ok);
  free(y_disc_ord_ok);
  free(y_cont_hmin);
  free(y_cont_k0);
  free(y_disc_uno_const);
  free(y_disc_ord_const);

  free(kernel_cx);
  free(kernel_cy);

  free(kernel_ux);
  free(kernel_uy);

  free(kernel_ox);
  free(kernel_oy);

  free(mean);

  free(matrix_wX_continuous_train);
  free(matrix_wX_unordered_train);
  free(matrix_wX_ordered_train);

  free(matrix_wY_continuous_train);
  free(matrix_wY_unordered_train);
  free(matrix_wY_ordered_train);

  free(matrix_wY_continuous_eval);
  free(matrix_wY_unordered_eval);
  free(matrix_wY_ordered_eval);

  free_mat(matrix_bandwidth_x, num_reg_continuous);
  free_mat(matrix_bandwidth_y, num_var_continuous);

  np_gate_override_clear();

  return(0);

}

int np_conditional_density_cvls_lp_stream(double *vector_scale_factor,
                                          double *cv);

int np_kernel_estimate_con_density_categorical_leave_one_out_ls_cv(
int KERNEL_var,
int KERNEL_unordered_var,
int KERNEL_ordered_var,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int64_t num_obs_train,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double memfac,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_XY_unordered_train, 
double **matrix_XY_ordered_train, 
double **matrix_XY_continuous_train, 
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *cv){
  NP_GateOverrideCtx gate_x_ctx, gate_y_ctx, gate_xy_ctx;
  np_gate_ctx_clear(&gate_x_ctx);
  np_gate_ctx_clear(&gate_y_ctx);
  np_gate_ctx_clear(&gate_xy_ctx);

  if((BANDWIDTH_den == BW_FIXED) || (BANDWIDTH_den == BW_GEN_NN) || (BANDWIDTH_den == BW_ADAP_NN))
    return np_conditional_density_cvls_lp_stream(vector_scale_factor, cv);

  int64_t i,j,k,l;

  int64_t iwi, iwj, iwk;

  const int num_reg_tot = num_reg_continuous+num_reg_unordered+num_reg_ordered;
  const int num_var_tot = num_var_continuous+num_var_unordered+num_var_ordered;
  const int num_all_var = num_reg_tot + num_var_tot;

  const int num_all_cvar = num_reg_continuous + num_var_continuous;
  const int num_all_uvar = num_reg_unordered + num_var_unordered;
  const int num_all_ovar = num_reg_ordered + num_var_ordered;

  size_t Nm = MIN((size_t)ceil(memfac*300000.0), (size_t)SIZE_MAX/10);

  int64_t N, num_obs_train_alloc, num_obs_wi_alloc, num_obs_wj_alloc, num_obs_wk_alloc;

  int64_t wi, wk, wj, nwi, nwk, nwj;

  int * x_operator = NULL, * y_operator = NULL, * xy_operator = NULL;
  int gate_x_active = 0, gate_y_active = 0, gate_xy_active = 0;
  int gate_x_all_large_fixed = 0;
  int *x_cont_ok = NULL, *x_disc_uno_ok = NULL, *x_disc_ord_ok = NULL;
  int *y_cont_ok = NULL, *y_disc_uno_ok = NULL, *y_disc_ord_ok = NULL;
  int *xy_cont_ok = NULL, *xy_disc_uno_ok = NULL, *xy_disc_ord_ok = NULL;
  double *x_cont_hmin = NULL, *x_cont_k0 = NULL, *x_disc_uno_const = NULL, *x_disc_ord_const = NULL;
  double x_all_large_fixed_const = 1.0;
  double *y_cont_hmin = NULL, *y_cont_k0 = NULL, *y_disc_uno_const = NULL, *y_disc_ord_const = NULL;
  double *xy_cont_hmin = NULL, *xy_cont_k0 = NULL, *xy_disc_uno_const = NULL, *xy_disc_ord_const = NULL;
  int bandwidth_provided = BANDWIDTH_den != BW_FIXED;

  double vsfx[num_reg_tot];
  double vsfy[num_var_tot];
  double vsfxy[num_var_tot+num_reg_tot];

  // we need the various bandwidth matrices to make sure that the correct bandwidths are used
  // when adaptive and/or generalized nearest neighbor bandwidths are selected
  // in combination with the blocking algorithm 

  double * lambdax = NULL;
  double * lambday = NULL;
  double * lambdaxy = NULL;

  double **matrix_bandwidth_x = NULL;
  double **matrix_bandwidth_xy = NULL;
  double **matrix_bandwidth_y = NULL;

  double **matrix_bandwidth_xi = NULL;
  double **matrix_bandwidth_xj = NULL;
  double **matrix_bandwidth_xk = NULL;

  double **matrix_bandwidth_yj = NULL;
  double **matrix_bandwidth_yk = NULL;

  double **matrix_Yk_unordered_train;
  double **matrix_Yk_ordered_train;
  double **matrix_Yk_continuous_train;

  double **matrix_Yj_unordered_train;
  double **matrix_Yj_ordered_train;
  double **matrix_Yj_continuous_train;

  double **matrix_Xi_unordered_train;
  double **matrix_Xi_ordered_train;
  double **matrix_Xi_continuous_train;

  double **matrix_Xj_unordered_train;
  double **matrix_Xj_ordered_train;
  double **matrix_Xj_continuous_train;

  double **matrix_Xk_unordered_train;
  double **matrix_Xk_ordered_train;
  double **matrix_Xk_continuous_train;

  double * pkx_ij = NULL, * pkx_ik = NULL, * pky_jk = NULL;

  double tcvk, tcvj;

  int64_t is, ie;
  int64_t is_i2n, ie_i2n;

  int m_ij, m_ik;

  NL nls = {.node = NULL, .n = 0, .nalloc = 0};

  XL xl_xij = {.istart = NULL, .nlev = NULL, .n = 0, .nalloc = 0};
  XL xl_xik = {.istart = NULL, .nlev = NULL, .n = 0, .nalloc = 0};

  nls.node = (int *)malloc(sizeof(int));
  nls.nalloc = 1;

  nls.node[0] = 0;
  nls.n = 1;
  
  int xyd[MAX(1,num_all_cvar)];
  int idxj[2], idxk[2];

  for(i = 0; i < num_all_cvar; i++)
    xyd[i] = i;

#ifdef MPI2
  int64_t stride_t = MAX((int64_t)ceil((double) num_obs_train / (double) iNum_Processors),1);

  num_obs_train_alloc = stride_t*iNum_Processors;

  is_i2n = stride_t*my_rank;
  ie_i2n = MIN(num_obs_train, is_i2n + stride_t);
#else
  num_obs_train_alloc = num_obs_train;

  is = is_i2n = 0;
  ie = ie_i2n = num_obs_train;
#endif

  // since we don't have an algorithm yet, it's hard to compute N exactly :)
  // blocking algo calculations
  // at the moment all blocks are of the same size
  // 3 weight matrices + one mean vector 
  N = (3*num_obs_train_alloc + 1)*num_obs_train_alloc;
  
  if(N > Nm){
    const int64_t s0 = (int64_t)ceil(sqrt((double)(Nm-num_obs_train_alloc)/3.0));
    wi = wk = wj = MIN(num_obs_train_alloc, s0);


    nwi = nwk = nwj = num_obs_train_alloc/wi + (((num_obs_train_alloc % wi) > 0) ? 1 : 0);
  } else {
    wi = wk = wj = num_obs_train_alloc;
    nwi = nwk = nwj = 1;
  }

#ifdef MPI2
  int64_t stride_wi = MAX((int64_t)ceil((double)wi / (double) iNum_Processors),1);

  num_obs_wi_alloc = num_obs_wj_alloc = num_obs_wk_alloc = stride_wi*iNum_Processors;

  is = stride_wi*my_rank;
  ie = MIN(num_obs_train, is + stride_wi);
#else
  num_obs_wi_alloc = num_obs_wj_alloc = num_obs_wk_alloc = wi;
#endif

  int * kernel_cx = NULL, * kernel_ux = NULL, * kernel_ox = NULL;
  int * kernel_cy = NULL, * kernel_uy = NULL, * kernel_oy = NULL;
  int * kernel_cxy = NULL, * kernel_uxy = NULL, * kernel_oxy = NULL;

  // first the x-kernels
  kernel_cx = (int *)malloc(sizeof(int)*num_reg_continuous);

  for(i = 0; i < num_reg_continuous; i++)
    kernel_cx[i] = KERNEL_reg;

  kernel_ux = (int *)malloc(sizeof(int)*num_reg_unordered);

  for(i = 0; i < num_reg_unordered; i++)
    kernel_ux[i] = KERNEL_unordered_reg;

  kernel_ox = (int *)malloc(sizeof(int)*num_reg_ordered);

  for(i = 0; i < num_reg_ordered; i++)
    kernel_ox[i] = KERNEL_ordered_reg;

  // then the y-kernels
  kernel_cy = (int *)malloc(sizeof(int)*num_var_continuous);

  for(i = 0; i < num_var_continuous; i++)
    kernel_cy[i] = KERNEL_var;

  kernel_uy = (int *)malloc(sizeof(int)*num_var_unordered);

  for(i = 0; i < num_var_unordered; i++)
    kernel_uy[i] = KERNEL_unordered_var;

  kernel_oy = (int *)malloc(sizeof(int)*num_var_ordered);

  for(i = 0; i < num_var_ordered; i++)
    kernel_oy[i] = KERNEL_ordered_var;

  // finally the xy-kernels
  kernel_cxy = (int *)malloc(sizeof(int)*num_all_cvar);

  for(i = 0; i < num_reg_continuous; i++)
    kernel_cxy[i] = KERNEL_reg;

  for(i = num_reg_continuous; i < num_all_cvar; i++)
    kernel_cxy[i] = KERNEL_var;

  kernel_uxy = (int *)malloc(sizeof(int)*num_all_uvar);

  for(i = 0; i < num_reg_unordered; i++)
    kernel_uxy[i] = KERNEL_unordered_reg;

  for(i = num_reg_unordered; i < num_all_uvar; i++)
    kernel_uxy[i] = KERNEL_unordered_var;

  kernel_oxy = (int *)malloc(sizeof(int)*num_all_ovar);

  for(i = 0; i < num_reg_ordered; i++)
    kernel_oxy[i] = KERNEL_ordered_reg;

  for(i = num_reg_ordered; i < num_all_ovar; i++)
    kernel_oxy[i] = KERNEL_ordered_var;

  // allocate some pointers
  matrix_Xi_continuous_train = (double **)malloc(num_reg_continuous*sizeof(double *));
  matrix_Xi_unordered_train = (double **)malloc(num_reg_unordered*sizeof(double *));
  matrix_Xi_ordered_train = (double **)malloc(num_reg_ordered*sizeof(double *));

  matrix_Xj_continuous_train = (double **)malloc(num_reg_continuous*sizeof(double *));
  matrix_Xj_unordered_train = (double **)malloc(num_reg_unordered*sizeof(double *));
  matrix_Xj_ordered_train = (double **)malloc(num_reg_ordered*sizeof(double *));

  matrix_Xk_continuous_train = (double **)malloc(num_reg_continuous*sizeof(double *));
  matrix_Xk_unordered_train = (double **)malloc(num_reg_unordered*sizeof(double *));
  matrix_Xk_ordered_train = (double **)malloc(num_reg_ordered*sizeof(double *));

  matrix_Yj_continuous_train = (double **)malloc(num_var_continuous*sizeof(double *));
  matrix_Yj_unordered_train = (double **)malloc(num_var_unordered*sizeof(double *));
  matrix_Yj_ordered_train = (double **)malloc(num_var_ordered*sizeof(double *));

  matrix_Yk_continuous_train = (double **)malloc(num_var_continuous*sizeof(double *));
  matrix_Yk_unordered_train = (double **)malloc(num_var_unordered*sizeof(double *));
  matrix_Yk_ordered_train = (double **)malloc(num_var_ordered*sizeof(double *));

  double * mean = (double *)malloc(num_obs_train_alloc*sizeof(double));
  double * jmean = (double *)malloc(num_obs_train_alloc*sizeof(double));

  if(mean == NULL)
    error("failed to allocate mean");

  np_splitxy_vsf_mcv_nc(num_var_unordered, num_var_ordered, num_var_continuous,
                        num_reg_unordered, num_reg_ordered, num_reg_continuous,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        vsfx,
                        vsfy,
                        vsfxy,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  x_operator = (int *)malloc(sizeof(int)*(num_reg_continuous+num_reg_unordered+num_reg_ordered));

  if(x_operator == NULL)
    error("failed to allocate x_operator");

  for(i = 0; i < (num_reg_continuous+num_reg_unordered+num_reg_ordered); i++)
    x_operator[i] = OP_NORMAL;

  y_operator = (int *)malloc(sizeof(int)*(num_var_continuous+num_var_unordered+num_var_ordered));

  if(y_operator == NULL)
    error("failed to allocate y_operator");

  for(i = 0; i < (num_var_continuous+num_var_unordered+num_var_ordered); i++)
    y_operator[i] = OP_CONVOLUTION;

  xy_operator = (int *)malloc(sizeof(int)*num_all_var);

  if(xy_operator == NULL)
    error("failed to allocate xy_operator");

  for(i = 0; i < num_reg_continuous; i++)
    xy_operator[i] = OP_NORMAL;

  for(i = num_reg_continuous; i < num_all_cvar; i++)
    xy_operator[i] = OP_NORMAL;

  // no cdf for unordered
  for(i = num_all_cvar; i < (num_all_cvar+num_all_uvar+num_reg_ordered); i++)
    xy_operator[i] = OP_NORMAL;

  for(i = num_all_cvar+num_all_uvar+num_reg_ordered; i < num_all_var; i++)
    xy_operator[i] = OP_NORMAL;

  // special bandwidths
  if(BANDWIDTH_den != BW_FIXED){
    matrix_bandwidth_x = alloc_tmatd(num_obs_train, num_reg_continuous);
    lambdax = alloc_vecd(num_reg_unordered+num_reg_ordered);

    kernel_bandwidth_mean(KERNEL_reg,
                          BANDWIDTH_den,
                          num_obs_train,
                          num_obs_train,
                          0,
                          0,
                          0,
                          num_reg_continuous,
                          num_reg_unordered,
                          num_reg_ordered,
                          0, // do not suppress_parallel
                          vsfx,
                          NULL,
                          NULL,
                          matrix_X_continuous_train,
                          matrix_X_continuous_train,
                          NULL,					 // Not used 
                          matrix_bandwidth_x,
                          lambdax);



    matrix_bandwidth_y = alloc_tmatd(num_obs_train, num_var_continuous);
    lambday = alloc_vecd(num_var_unordered+num_var_ordered);

    kernel_bandwidth_mean(KERNEL_var,
                          BANDWIDTH_den,
                          num_obs_train,
                          num_obs_train,
                          0,
                          0,
                          0,
                          num_var_continuous,
                          num_var_unordered,
                          num_var_ordered,
                          0, // do not suppress_parallel
                          vsfy,
                          NULL,
                          NULL,
                          matrix_Y_continuous_train,
                          matrix_Y_continuous_train,
                          NULL,					 // Not used 
                          matrix_bandwidth_y,
                          lambday);

    matrix_bandwidth_xy = alloc_tmatd(num_obs_train, num_all_cvar);
    lambdaxy = alloc_vecd(num_all_ovar+num_all_uvar);

    kernel_bandwidth_mean(KERNEL_reg,
                          BANDWIDTH_den,
                          num_obs_train,
                          num_obs_train,
                          0,
                          0,
                          0,
                          num_all_cvar,
                          num_all_uvar,
                          num_all_ovar,
                          0, // do not suppress_parallel
                          vsfxy,
                          NULL,
                          NULL,
                          matrix_XY_continuous_train,
                          matrix_XY_continuous_train,
                          NULL,					 // Not used 
                          matrix_bandwidth_xy,
                          lambdaxy);

    matrix_bandwidth_yj = (double **)malloc(sizeof(double *)*num_var_continuous);
    matrix_bandwidth_yk = (double **)malloc(sizeof(double *)*num_var_continuous);

    if(BANDWIDTH_den == BW_ADAP_NN){
      matrix_bandwidth_xi = (double **)malloc(sizeof(double *)*num_reg_continuous);
    }else{
      matrix_bandwidth_xj = (double **)malloc(sizeof(double *)*num_reg_continuous);
      matrix_bandwidth_xk = (double **)malloc(sizeof(double *)*num_reg_continuous);
    }

  } else if (BANDWIDTH_den == BW_FIXED) {
    matrix_bandwidth_x = (double **)malloc(sizeof(double*));
    matrix_bandwidth_x[0] = vsfx;
  }

  if(num_reg_continuous > 0 || num_reg_unordered > 0 || num_reg_ordered > 0){
    int ok_all = 1;

    if(num_reg_continuous > 0){
      x_cont_ok = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
      x_cont_hmin = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
      x_cont_k0 = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
      ok_all = (x_cont_ok != NULL) && (x_cont_hmin != NULL) && (x_cont_k0 != NULL);
      if(ok_all){
        const double rel_tol = np_cont_largeh_rel_tol();
        for(i = 0; i < num_reg_continuous; i++){
          const int kern = kernel_cx[i];
          double xmin = DBL_MAX, xmax = -DBL_MAX;
          x_cont_ok[i] = 0; x_cont_hmin[i] = DBL_MAX; x_cont_k0[i] = 0.0;
          if(!np_cont_largeh_kernel_supported(kern)) continue;
          for(j = 0; j < num_obs_train; j++){
            const double v = matrix_X_continuous_train[i][j];
            if(!isfinite(v)) continue;
            xmin = MIN(xmin, v); xmax = MAX(xmax, v);
          }
          if(xmax >= xmin){
            const double utol = np_cont_largeh_utol(kern, rel_tol);
            if(utol > 0.0 && isfinite(utol)){
              x_cont_ok[i] = 1;
              x_cont_hmin[i] = (xmax - xmin)/utol;
              x_cont_k0[i] = np_cont_largeh_k0(kern);
            }
          }
        }
      }
    }

    if(ok_all && (lambdax != NULL) && num_reg_unordered > 0){
      x_disc_uno_ok = (int *)calloc((size_t)num_reg_unordered, sizeof(int));
      x_disc_uno_const = (double *)malloc((size_t)num_reg_unordered*sizeof(double));
      ok_all = (x_disc_uno_ok != NULL) && (x_disc_uno_const != NULL);
      if(ok_all){
        double (* const ukf[])(int, double, int) = {
          np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
          np_score_uaa, np_score_unli_racine
        };
        const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
        for(i = 0; i < num_reg_unordered; i++){
          const int ku = kernel_ux[i];
          const int ncat = (num_categories_extern_X != NULL) ? num_categories_extern_X[i] : 0;
          const double lam = lambdax[i];
          x_disc_uno_ok[i] = 0; x_disc_uno_const[i] = 0.0;
          if(ku < 0 || ku >= nuk) continue;
          if(!np_disc_near_upper(ku, lam, ncat)) continue;
          {
            const double ks = ukf[ku](1, lam, ncat);
            const double kd = ukf[ku](0, lam, ncat);
            if(np_disc_near_const_kernel(ks, kd)){
              x_disc_uno_ok[i] = 1;
              x_disc_uno_const[i] = 0.5*(ks + kd);
            }
          }
        }
      }
    }

    if(ok_all && (lambdax != NULL) && num_reg_ordered > 0){
      x_disc_ord_ok = (int *)calloc((size_t)num_reg_ordered, sizeof(int));
      x_disc_ord_const = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
      ok_all = (x_disc_ord_ok != NULL) && (x_disc_ord_const != NULL);
      if(ok_all){
        double (* const okf[])(double, double, double, double, double) = {
          np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
        };
        const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
        for(i = 0; i < num_reg_ordered; i++){
          const int oi = i + num_reg_unordered;
          const int ko = kernel_ox[i];
          const int ncat = (num_categories_extern_X != NULL) ? num_categories_extern_X[oi] : 0;
          const double lam = lambdax[oi];
          x_disc_ord_ok[i] = 0; x_disc_ord_const[i] = 0.0;
          if(ko < 0 || ko >= nok) continue;
          if(ncat <= 0 || matrix_categorical_vals_extern_X == NULL) continue;
          if(!np_disc_ordered_near_upper(ko, lam)) continue;
          {
            const double cl = matrix_categorical_vals_extern_X[oi][0];
            const double ch = matrix_categorical_vals_extern_X[oi][ncat - 1];
            const double k0 = okf[ko](cl, cl, lam, cl, ch);
            const double k1 = okf[ko](cl, ch, lam, cl, ch);
            if(np_disc_near_const_kernel(k0, k1)){
              x_disc_ord_ok[i] = 1;
              x_disc_ord_const[i] = 0.5*(k0 + k1);
            }
          }
        }
      }
    }

    if(ok_all){
      gate_x_active = 1;
    }
  }

  if((BANDWIDTH_den == BW_FIXED) && (num_reg_tot > 0)){
    int ok_all_large = 1;
    double kconst = 1.0;

    for(l = 0; l < num_reg_continuous; l++){
      const double bw = vsfx[l];
      if((x_cont_ok == NULL) || (!x_cont_ok[l]) || (!isfinite(bw)) ||
         (bw <= 0.0) || (bw < x_cont_hmin[l])){
        ok_all_large = 0;
        break;
      }
      kconst *= x_cont_k0[l]/bw;
    }

    if(ok_all_large){
      double (* const ukf[])(int, double, int) = {
        np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
        np_score_uaa, np_score_unli_racine
      };
      const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
      for(l = 0; l < num_reg_unordered; l++){
        const int ku = kernel_ux[l];
        const int oi = num_reg_continuous + l;
        const int ncat = (num_categories_extern_X != NULL) ? num_categories_extern_X[l] : 0;
        const double lam = vsfx[oi];
        if((ku < 0) || (ku >= nuk) || (!isfinite(lam)) || (!np_disc_near_upper(ku, lam, ncat))){
          ok_all_large = 0;
          break;
        }
        {
          const double ks = ukf[ku](1, lam, ncat);
          const double kd = ukf[ku](0, lam, ncat);
          if(!np_disc_near_const_kernel(ks, kd)){
            ok_all_large = 0;
            break;
          }
          kconst *= 0.5*(ks + kd);
        }
      }
    }

    if(ok_all_large){
      double (* const okf[])(double, double, double, double, double) = {
        np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
      };
      const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
      for(l = 0; l < num_reg_ordered; l++){
        const int ko = kernel_ox[l];
        const int oi = num_reg_unordered + l;
        const int vi = num_reg_continuous + oi;
        const int ncat = (num_categories_extern_X != NULL) ? num_categories_extern_X[oi] : 0;
        const double lam = vsfx[vi];
        if((ko < 0) || (ko >= nok) || (!isfinite(lam)) || (!np_disc_ordered_near_upper(ko, lam)) ||
           (ncat <= 0) || (matrix_categorical_vals_extern_X == NULL)){
          ok_all_large = 0;
          break;
        }
        {
          const double cl = matrix_categorical_vals_extern_X[oi][0];
          const double ch = matrix_categorical_vals_extern_X[oi][ncat - 1];
          const double k0 = okf[ko](cl, cl, lam, cl, ch);
          const double k1 = okf[ko](cl, ch, lam, cl, ch);
          if(!np_disc_near_const_kernel(k0, k1)){
            ok_all_large = 0;
            break;
          }
          kconst *= 0.5*(k0 + k1);
        }
      }
    }

    if(ok_all_large && isfinite(kconst) && (kconst > 0.0)){
      gate_x_all_large_fixed = 1;
      x_all_large_fixed_const = kconst;
    }
  }
  if(gate_x_all_large_fixed)
    np_fastcv_alllarge_hits++;

  if(num_var_continuous > 0 || num_var_unordered > 0 || num_var_ordered > 0){
    int ok_all = 1;

    if(num_var_continuous > 0){
      y_cont_ok = (int *)calloc((size_t)num_var_continuous, sizeof(int));
      y_cont_hmin = (double *)malloc((size_t)num_var_continuous*sizeof(double));
      y_cont_k0 = (double *)malloc((size_t)num_var_continuous*sizeof(double));
      ok_all = (y_cont_ok != NULL) && (y_cont_hmin != NULL) && (y_cont_k0 != NULL);
      if(ok_all){
        const double rel_tol = np_cont_largeh_rel_tol();
        for(i = 0; i < num_var_continuous; i++){
          const int kern = kernel_cy[i];
          double xmin = DBL_MAX, xmax = -DBL_MAX;
          y_cont_ok[i] = 0; y_cont_hmin[i] = DBL_MAX; y_cont_k0[i] = 0.0;
          if(!np_cont_largeh_kernel_supported(kern)) continue;
          for(j = 0; j < num_obs_train; j++){
            const double v = matrix_Y_continuous_train[i][j];
            if(!isfinite(v)) continue;
            xmin = MIN(xmin, v); xmax = MAX(xmax, v);
          }
          if(xmax >= xmin){
            const double utol = np_cont_largeh_utol(kern, rel_tol);
            if(utol > 0.0 && isfinite(utol)){
              y_cont_ok[i] = 1;
              y_cont_hmin[i] = (xmax - xmin)/utol;
              y_cont_k0[i] = np_cont_largeh_k0(kern);
            }
          }
        }
      }
    }

    if(ok_all && (lambday != NULL) && num_var_unordered > 0){
      y_disc_uno_ok = (int *)calloc((size_t)num_var_unordered, sizeof(int));
      y_disc_uno_const = (double *)malloc((size_t)num_var_unordered*sizeof(double));
      ok_all = (y_disc_uno_ok != NULL) && (y_disc_uno_const != NULL);
      if(ok_all){
        double (* const ukf[])(int, double, int) = {
          np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
          np_score_uaa, np_score_unli_racine
        };
        const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
        for(i = 0; i < num_var_unordered; i++){
          const int ku = kernel_uy[i];
          const int ncat = (num_categories_extern_Y != NULL) ? num_categories_extern_Y[i] : 0;
          const double lam = lambday[i];
          y_disc_uno_ok[i] = 0; y_disc_uno_const[i] = 0.0;
          if(ku < 0 || ku >= nuk) continue;
          if(!np_disc_near_upper(ku, lam, ncat)) continue;
          {
            const double ks = ukf[ku](1, lam, ncat);
            const double kd = ukf[ku](0, lam, ncat);
            if(np_disc_near_const_kernel(ks, kd)){
              y_disc_uno_ok[i] = 1;
              y_disc_uno_const[i] = 0.5*(ks + kd);
            }
          }
        }
      }
    }

    if(ok_all && (lambday != NULL) && num_var_ordered > 0){
      y_disc_ord_ok = (int *)calloc((size_t)num_var_ordered, sizeof(int));
      y_disc_ord_const = (double *)malloc((size_t)num_var_ordered*sizeof(double));
      ok_all = (y_disc_ord_ok != NULL) && (y_disc_ord_const != NULL);
      if(ok_all){
        double (* const okf[])(double, double, double, double, double) = {
          np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
        };
        const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
        for(i = 0; i < num_var_ordered; i++){
          const int oi = i + num_var_unordered;
          const int ko = kernel_oy[i];
          const int ncat = (num_categories_extern_Y != NULL) ? num_categories_extern_Y[oi] : 0;
          const double lam = lambday[oi];
          y_disc_ord_ok[i] = 0; y_disc_ord_const[i] = 0.0;
          if(ko < 0 || ko >= nok) continue;
          if(ncat <= 0 || matrix_categorical_vals_extern_Y == NULL) continue;
          if(!np_disc_ordered_near_upper(ko, lam)) continue;
          {
            const double cl = matrix_categorical_vals_extern_Y[oi][0];
            const double ch = matrix_categorical_vals_extern_Y[oi][ncat - 1];
            const double k0 = okf[ko](cl, cl, lam, cl, ch);
            const double k1 = okf[ko](cl, ch, lam, cl, ch);
            if(np_disc_near_const_kernel(k0, k1)){
              y_disc_ord_ok[i] = 1;
              y_disc_ord_const[i] = 0.5*(k0 + k1);
            }
          }
        }
      }
    }

    if(ok_all){
      gate_y_active = 1;
    }
  }

  if(num_all_cvar > 0 || num_all_uvar > 0 || num_all_ovar > 0){
    int ok_all = 1;

    if(num_all_cvar > 0){
      xy_cont_ok = (int *)calloc((size_t)num_all_cvar, sizeof(int));
      xy_cont_hmin = (double *)malloc((size_t)num_all_cvar*sizeof(double));
      xy_cont_k0 = (double *)malloc((size_t)num_all_cvar*sizeof(double));
      ok_all = (xy_cont_ok != NULL) && (xy_cont_hmin != NULL) && (xy_cont_k0 != NULL);
      if(ok_all){
        const double rel_tol = np_cont_largeh_rel_tol();
        for(i = 0; i < num_all_cvar; i++){
          const int kern = kernel_cxy[i];
          double xmin = DBL_MAX, xmax = -DBL_MAX;
          xy_cont_ok[i] = 0; xy_cont_hmin[i] = DBL_MAX; xy_cont_k0[i] = 0.0;
          if(!np_cont_largeh_kernel_supported(kern)) continue;
          for(j = 0; j < num_obs_train; j++){
            const double v = matrix_XY_continuous_train[i][j];
            if(!isfinite(v)) continue;
            xmin = MIN(xmin, v); xmax = MAX(xmax, v);
          }
          if(xmax >= xmin){
            const double utol = np_cont_largeh_utol(kern, rel_tol);
            if(utol > 0.0 && isfinite(utol)){
              xy_cont_ok[i] = 1;
              xy_cont_hmin[i] = (xmax - xmin)/utol;
              xy_cont_k0[i] = np_cont_largeh_k0(kern);
            }
          }
        }
      }
    }

    if(ok_all && (lambdaxy != NULL) && num_all_uvar > 0){
      xy_disc_uno_ok = (int *)calloc((size_t)num_all_uvar, sizeof(int));
      xy_disc_uno_const = (double *)malloc((size_t)num_all_uvar*sizeof(double));
      ok_all = (xy_disc_uno_ok != NULL) && (xy_disc_uno_const != NULL);
      if(ok_all){
        double (* const ukf[])(int, double, int) = {
          np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
          np_score_uaa, np_score_unli_racine
        };
        const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
        for(i = 0; i < num_all_uvar; i++){
          const int ku = kernel_uxy[i];
          const int ncat = (num_categories_extern_XY != NULL) ? num_categories_extern_XY[i] : 0;
          const double lam = lambdaxy[i];
          xy_disc_uno_ok[i] = 0; xy_disc_uno_const[i] = 0.0;
          if(ku < 0 || ku >= nuk) continue;
          if(!np_disc_near_upper(ku, lam, ncat)) continue;
          {
            const double ks = ukf[ku](1, lam, ncat);
            const double kd = ukf[ku](0, lam, ncat);
            if(np_disc_near_const_kernel(ks, kd)){
              xy_disc_uno_ok[i] = 1;
              xy_disc_uno_const[i] = 0.5*(ks + kd);
            }
          }
        }
      }
    }

    if(ok_all && (lambdaxy != NULL) && num_all_ovar > 0){
      xy_disc_ord_ok = (int *)calloc((size_t)num_all_ovar, sizeof(int));
      xy_disc_ord_const = (double *)malloc((size_t)num_all_ovar*sizeof(double));
      ok_all = (xy_disc_ord_ok != NULL) && (xy_disc_ord_const != NULL);
      if(ok_all){
        double (* const okf[])(double, double, double, double, double) = {
          np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
        };
        const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
        for(i = 0; i < num_all_ovar; i++){
          const int oi = i + num_all_uvar;
          const int ko = kernel_oxy[i];
          const int ncat = (num_categories_extern_XY != NULL) ? num_categories_extern_XY[oi] : 0;
          const double lam = lambdaxy[oi];
          xy_disc_ord_ok[i] = 0; xy_disc_ord_const[i] = 0.0;
          if(ko < 0 || ko >= nok) continue;
          if(ncat <= 0 || matrix_categorical_vals_extern_XY == NULL) continue;
          if(!np_disc_ordered_near_upper(ko, lam)) continue;
          {
            const double cl = matrix_categorical_vals_extern_XY[oi][0];
            const double ch = matrix_categorical_vals_extern_XY[oi][ncat - 1];
            const double k0 = okf[ko](cl, cl, lam, cl, ch);
            const double k1 = okf[ko](cl, ch, lam, cl, ch);
            if(np_disc_near_const_kernel(k0, k1)){
              xy_disc_ord_ok[i] = 1;
              xy_disc_ord_const[i] = 0.5*(k0 + k1);
            }
          }
        }
      }
    }

    if(ok_all){
      gate_xy_active = 1;
    }
  }

  // extra kernel bookkeeping for trees
  int KERNEL_XY[MAX(1,num_all_cvar)];
  double bb[MAX(1,2*num_all_cvar)];

  for(l = 0; l < num_reg_continuous; l++)
    KERNEL_XY[l] = KERNEL_reg + OP_CFUN_OFFSETS[x_operator[l]];

  for(l = num_reg_continuous; l < num_all_cvar; l++)
    KERNEL_XY[l] = KERNEL_var + OP_CFUN_OFFSETS[y_operator[l-num_reg_continuous]];
  
  *cv = 0;

  int fast_joint_y_only = 0;
  int y_joint_operator[MAX(1,num_var_tot)];
  // joint density
  if(gate_x_all_large_fixed && (BANDWIDTH_den == BW_FIXED)){
    fast_joint_y_only = 1;
    for(i = 0; i < (num_var_continuous+num_var_unordered+num_var_ordered); i++)
      y_joint_operator[i] = OP_NORMAL;

    np_gate_ctx_clear(&gate_y_ctx);
    np_activate_bounds_y();
    kernel_weighted_sum_np_ctx(kernel_cy,
                           kernel_uy,
                           kernel_oy,
                           BANDWIDTH_den,
                           num_obs_train,
                           num_obs_train,
                           num_var_unordered,
                           num_var_ordered,
                           num_var_continuous,
                           1, // compute the leave-one-out marginals
                           0,
                           1, // kpow = 1
                           1, // bw divide
                           0,
                           0,
                           0,
                           0,
                           0,
                           y_joint_operator,
                           OP_NOOP, // no permutations
                           0, // no score
                           0, // no ocg
                           NULL,
                           0, // don't explicity suppress parallel
                           0,
                           0,
                           0, // do not use tree to keep index order direct
                           0,
                           NULL,
                           NULL, NULL, NULL,
                           matrix_Y_unordered_train,
                           matrix_Y_ordered_train,
                           matrix_Y_continuous_train,
                           matrix_Y_unordered_train,
                           matrix_Y_ordered_train,
                           matrix_Y_continuous_train,
                           NULL,
                           NULL,
                           NULL,
                           vsfy,
                           bandwidth_provided,
                           matrix_bandwidth_y,
                           matrix_bandwidth_y,
                           lambday,
                           num_categories_extern_Y,
                           matrix_categorical_vals_extern_Y,
                           NULL,
                           jmean,
                           NULL, // no permutations
                           NULL,
                           &gate_y_ctx);
  } else {
    if(gate_xy_active){
      np_gate_ctx_set(&gate_xy_ctx,
                      num_all_cvar,
                      num_all_uvar,
                      num_all_ovar,
                      kernel_cxy,
                      kernel_uxy,
                      kernel_oxy,
                      xy_operator,
                      xy_cont_ok,
                      xy_cont_hmin,
                      xy_cont_k0,
                      xy_disc_uno_ok,
                      xy_disc_uno_const,
                      xy_disc_ord_ok,
                      xy_disc_ord_const);
    } else {
      np_gate_ctx_clear(&gate_xy_ctx);
    }
    np_activate_bounds_xy();
    kernel_weighted_sum_np_ctx(kernel_cxy,
                           kernel_uxy,
                           kernel_oxy,
                           BANDWIDTH_den,
                           num_obs_train,
                           num_obs_train,
                           num_all_uvar,
                           num_all_ovar,
                           num_all_cvar,
                           1, // compute the leave-one-out marginals
                           0,
                           1, // kpow = 1
                           1, // bw divide
                           0, 
                           0,
                           0,
                           0,
                           0,
                           xy_operator,
                           OP_NOOP, // no permutations
                           0, // no score
                           0, // no ocg
                           NULL,
                           0, // don't explicity suppress parallel
                           0,
                           0,
                           int_TREE_XY,
                           0,
                           kdt_extern_XY,
                           NULL, NULL, NULL,
                           matrix_XY_unordered_train,
                           matrix_XY_ordered_train,
                           matrix_XY_continuous_train,
                           matrix_XY_unordered_train,
                           matrix_XY_ordered_train,
                           matrix_XY_continuous_train,
                           NULL,
                           NULL,
                           NULL,
                           vsfxy,
                           bandwidth_provided,
                           matrix_bandwidth_xy,
                           matrix_bandwidth_xy,
                           lambdaxy,
                           num_categories_extern_XY,
                           matrix_categorical_vals_extern_XY,
                           NULL,
                           jmean,
                           NULL, // no permutations
                           NULL,
                           &gate_xy_ctx);
  }

  // X density
  if(gate_x_all_large_fixed){
    const double cmean = ((double)(num_obs_train - 1))*x_all_large_fixed_const;
    for(i = 0; i < num_obs_train_alloc; i++)
      mean[i] = cmean;
  } else {
    if(gate_x_active){
      np_gate_ctx_set(&gate_x_ctx,
                      num_reg_continuous,
                      num_reg_unordered,
                      num_reg_ordered,
                      kernel_cx,
                      kernel_ux,
                      kernel_ox,
                      x_operator,
                      x_cont_ok,
                      x_cont_hmin,
                      x_cont_k0,
                      x_disc_uno_ok,
                      x_disc_uno_const,
                      x_disc_ord_ok,
                      x_disc_ord_const);
    } else {
      np_gate_ctx_clear(&gate_x_ctx);
    }
    np_activate_bounds_x();
    kernel_weighted_sum_np_ctx(kernel_cx,
                           kernel_ux,
                           kernel_ox,
                           BANDWIDTH_den,
                           num_obs_train,
                           num_obs_train,
                           num_reg_unordered,
                           num_reg_ordered,
                           num_reg_continuous,
                           1, // compute the leave-one-out marginals
                           0,
                           1, // kpow
                           1, // bw divide
                           0, 
                           0,
                           0,
                           0,
                           0,
                           x_operator,
                           OP_NOOP, // no permutations
                           0, // no score
                           0, // no ocg
                           NULL,
                           0, // don't explicity suppress parallel
                           0,
                           0,
                           int_TREE_X,
                           0,
                           kdt_extern_X,
                           NULL, NULL, NULL,
                           matrix_X_unordered_train,
                           matrix_X_ordered_train,
                           matrix_X_continuous_train,
                           matrix_X_unordered_train,
                           matrix_X_ordered_train,
                           matrix_X_continuous_train,
                           NULL,
                           NULL,
                           NULL,
                           vsfx,
                           bandwidth_provided,
                           matrix_bandwidth_x,
                           matrix_bandwidth_x,
                           lambdax,
                           num_categories_extern_X,
                           matrix_categorical_vals_extern_X,
                           NULL,
                           mean,
                           NULL, // no permutations
                           NULL,
                           &gate_x_ctx);
  }

  if(fast_joint_y_only){
    const double inv_nmo = 1.0/(((double)(num_obs_train - 1)) + DBL_MIN);
    for(i = is_i2n; i < ie_i2n; i++)
      *cv -= 2.0*jmean[i]*inv_nmo;
  } else if((!int_TREE_XY) && (!int_TREE_X)){
    for(i = is_i2n; i < ie_i2n; i++)
      *cv -= 2.0*jmean[i]/(mean[i] + DBL_MIN);
  } else {
    for(i = is_i2n; i < ie_i2n; i++)
      *cv -= 2.0*jmean[ipt_lookup_extern_XY[ipt_extern_X[i]]]/(mean[i] + DBL_MIN);
  }

  free(jmean);

  double *ky_jk = NULL;
  ky_jk = (double *)malloc(num_obs_wj_alloc*num_obs_wk_alloc*sizeof(double));
  pky_jk = ky_jk;

  if(ky_jk == NULL)
    error("failed to allocate ky_jk, tried to allocate: %" PRIi64 "bytes\n", num_obs_wj_alloc*num_obs_wk_alloc*sizeof(double));

  if(gate_x_all_large_fixed && (BANDWIDTH_den == BW_FIXED)){
    double *row_sum = (double *)calloc((size_t)num_obs_train_alloc, sizeof(double));
    double *col_sum = (double *)calloc((size_t)num_obs_train_alloc, sizeof(double));
    double *diag_sum = (double *)calloc((size_t)num_obs_train_alloc, sizeof(double));
    double total_sum = 0.0;

    if((row_sum == NULL) || (col_sum == NULL) || (diag_sum == NULL))
      error("failed to allocate fast-path aggregate buffers");

    if(gate_y_active){
      np_gate_ctx_set(&gate_y_ctx,
                      num_var_continuous,
                      num_var_unordered,
                      num_var_ordered,
                      kernel_cy,
                      kernel_uy,
                      kernel_oy,
                      y_operator,
                      y_cont_ok,
                      y_cont_hmin,
                      y_cont_k0,
                      y_disc_uno_ok,
                      y_disc_uno_const,
                      y_disc_ord_ok,
                      y_disc_ord_const);
    } else {
      np_gate_ctx_clear(&gate_y_ctx);
    }

    for(iwj = 0; iwj < nwj; iwj++){
      const int64_t wjo = iwj*wj;
      const int64_t dwj = (iwj != (nwj - 1)) ? wj : num_obs_train - (nwj - 1)*wj;

      idxj[0] = wjo;
      idxj[1] = wjo + dwj -1;

      for(l = 0; l < num_var_continuous; l++)
        matrix_Yj_continuous_train[l] = matrix_XY_continuous_train[l+num_reg_continuous] + wjo;

      for(l = 0; l < num_var_unordered; l++)
        matrix_Yj_unordered_train[l] = matrix_XY_unordered_train[l+num_reg_unordered] + wjo;

      for(l = 0; l < num_var_ordered; l++)
        matrix_Yj_ordered_train[l] = matrix_XY_ordered_train[l+num_reg_ordered] + wjo;

      for(iwk = 0; iwk < nwk; iwk++){
        const int64_t wko = iwk*wk;
        const int64_t dwk = (iwk != (nwk - 1)) ? wk : num_obs_train - (nwk - 1)*wk;

        idxk[0] = wko;
        idxk[1] = wko + dwk -1;

        for(l = 0; l < num_var_continuous; l++)
          matrix_Yk_continuous_train[l] = matrix_XY_continuous_train[l+num_reg_continuous] + wko;

        for(l = 0; l < num_var_unordered; l++)
          matrix_Yk_unordered_train[l] = matrix_XY_unordered_train[l+num_reg_unordered] + wko;

        for(l = 0; l < num_var_ordered; l++)
          matrix_Yk_ordered_train[l] = matrix_XY_ordered_train[l+num_reg_ordered] + wko;

        np_activate_bounds_y();
        kernel_weighted_sum_np_ctx(kernel_cy,
                               kernel_uy,
                               kernel_oy,
                               BANDWIDTH_den,
                               dwk,
                               dwj,
                               num_var_unordered,
                               num_var_ordered,
                               num_var_continuous,
                               0,
                               0,
                               1,
                               1,
                               1,
                               0,
                               0,
                               0,
                               0,
                               y_operator,
                               OP_NOOP,
                               0,
                               0,
                               NULL,
                               0,
                               0,
                               0,
                               int_TREE_XY,
                               1,
                               kdt_extern_XY,
                               &nls,
                               xyd + num_reg_continuous,
                               idxk,
                               matrix_Yk_unordered_train,
                               matrix_Yk_ordered_train,
                               matrix_Yk_continuous_train,
                               matrix_Yj_unordered_train,
                               matrix_Yj_ordered_train,
                               matrix_Yj_continuous_train,
                               NULL,
                               NULL,
                               NULL,
                               vsfy,
                               bandwidth_provided,
                               matrix_bandwidth_yj,
                               matrix_bandwidth_yk,
                               lambday,
                               num_categories_extern_Y,
                               matrix_categorical_vals_extern_Y,
                               NULL,
                               NULL,
                               NULL,
                               ky_jk,
                               &gate_y_ctx);

        for(j = 0; j < dwj; j++){
          const int64_t jg = wjo + j;
          const int64_t joff = j*dwk;
          for(k = 0; k < dwk; k++){
            const int64_t kg = wko + k;
            const double v = ky_jk[joff + k];
            total_sum += v;
            row_sum[jg] += v;
            col_sum[kg] += v;
            if(jg == kg) diag_sum[jg] += v;
          }
        }
      }
    }

    {
      const double inv_m2 = 1.0/(((double)(num_obs_train - 1))*((double)(num_obs_train - 1)) + DBL_MIN);
      for(i = is_i2n; i < ie_i2n; i++){
        *cv += (total_sum - row_sum[i] - col_sum[i] + diag_sum[i])*inv_m2;
      }
    }

    free(row_sum);
    free(col_sum);
    free(diag_sum);
  } else {
    double * kx_ij = (double *)malloc(num_obs_wi_alloc*num_obs_wj_alloc*sizeof(double));
    pkx_ij = kx_ij;

    if(kx_ij == NULL)
      error("failed to allocate kx_ij, tried to allocate: %" PRIi64 "bytes\n", num_obs_wi_alloc*num_obs_wj_alloc*sizeof(double));

    double * kx_ik = (double *)malloc(num_obs_wi_alloc*num_obs_wk_alloc*sizeof(double));
    pkx_ik = kx_ik;

    if(kx_ik == NULL)
      error("failed to allocate kx_ik, tried to allocate: %" PRIi64 "bytes\n", num_obs_wi_alloc*num_obs_wk_alloc*sizeof(double));

    if(gate_x_active){
      np_gate_ctx_set(&gate_x_ctx,
                      num_reg_continuous,
                      num_reg_unordered,
                      num_reg_ordered,
                      kernel_cx,
                      kernel_ux,
                      kernel_ox,
                      x_operator,
                      x_cont_ok,
                      x_cont_hmin,
                      x_cont_k0,
                      x_disc_uno_ok,
                      x_disc_uno_const,
                      x_disc_ord_ok,
                      x_disc_ord_const);
    } else {
      np_gate_ctx_clear(&gate_x_ctx);
    }

    if(gate_y_active){
      np_gate_ctx_set(&gate_y_ctx,
                      num_var_continuous,
                      num_var_unordered,
                      num_var_ordered,
                      kernel_cy,
                      kernel_uy,
                      kernel_oy,
                      y_operator,
                      y_cont_ok,
                      y_cont_hmin,
                      y_cont_k0,
                      y_disc_uno_ok,
                      y_disc_uno_const,
                      y_disc_ord_ok,
                      y_disc_ord_const);
    } else {
      np_gate_ctx_clear(&gate_y_ctx);
    }

    for(iwi = 0; iwi < nwi; iwi++){
    const int64_t wio = iwi*wi;
    const int64_t dwi = (iwi != (nwi - 1)) ? wi : num_obs_train - (nwi - 1)*wi;

    for(l = 0; l < num_reg_continuous; l++)
      matrix_Xi_continuous_train[l] = matrix_XY_continuous_train[l] + wio;

    for(l = 0; l < num_reg_unordered; l++)
      matrix_Xi_unordered_train[l] = matrix_XY_unordered_train[l] + wio;

    for(l = 0; l < num_reg_ordered; l++)
      matrix_Xi_ordered_train[l] = matrix_XY_ordered_train[l] + wio;

    // offset the appropriate bandwidths
    if(BANDWIDTH_den == BW_ADAP_NN){
      for(l = 0; l < num_reg_continuous; l++)
        matrix_bandwidth_xi[l] = matrix_bandwidth_x[l] + wio;
    }

    for(iwj = 0; iwj < nwj; iwj++){
      const int64_t wjo = iwj*wj;
      const int64_t dwj = (iwj != (nwj - 1)) ? wj : num_obs_train - (nwj - 1)*wj;

      idxj[0] = wjo;
      idxj[1] = wjo + dwj -1;

      for(l = 0; l < num_reg_continuous; l++)
        matrix_Xj_continuous_train[l] = matrix_XY_continuous_train[l] + wjo;

      for(l = 0; l < num_reg_unordered; l++)
        matrix_Xj_unordered_train[l] = matrix_XY_unordered_train[l] + wjo;

      for(l = 0; l < num_reg_ordered; l++)
        matrix_Xj_ordered_train[l] = matrix_XY_ordered_train[l] + wjo;


      for(l = 0; l < num_var_continuous; l++)
        matrix_Yj_continuous_train[l] = matrix_XY_continuous_train[l+num_reg_continuous] + wjo;

      for(l = 0; l < num_var_unordered; l++)
        matrix_Yj_unordered_train[l] = matrix_XY_unordered_train[l+num_reg_unordered] + wjo;

      for(l = 0; l < num_var_ordered; l++)
        matrix_Yj_ordered_train[l] = matrix_XY_ordered_train[l+num_reg_ordered] + wjo;

        // offset the appropriate bandwidths
      if(BANDWIDTH_den != BW_FIXED){
        for(l = 0; l < num_var_continuous; l++)
          matrix_bandwidth_yj[l] = matrix_bandwidth_y[l] + wjo;

        if(BANDWIDTH_den == BW_GEN_NN)
          for(l = 0; l < num_reg_continuous; l++)
            matrix_bandwidth_xj[l] = matrix_bandwidth_x[l] + wjo;
      }

      // compute block kx_ij
      // i is eval, j is train

      if(gate_x_all_large_fixed){
        for(i = 0; i < (dwi*dwj); i++)
          kx_ij[i] = x_all_large_fixed_const;
      } else {
        np_activate_bounds_x();
        kernel_weighted_sum_np_ctx(kernel_cx,
                               kernel_ux,
                               kernel_ox,
                               BANDWIDTH_den,
                               dwj,
                               dwi,
                               num_reg_unordered,
                               num_reg_ordered,
                               num_reg_continuous,
                               0, // (do not) compute the leave-one-out marginals
                               0,
                               1, // kpow
                               1, // bw divide
                               1, // divide weights 
                               0,
                               0,
                               0,
                               0,
                               x_operator,
                               OP_NOOP, // no permutations
                               0, // no score
                               0, // no ocg
                               NULL,
                               0, // don't explicity suppress parallel
                               0,
                               0,
                               int_TREE_XY,
                               1,
                               kdt_extern_XY,
                               &nls, 
                               xyd,
                               idxj,
                               matrix_Xj_unordered_train,
                               matrix_Xj_ordered_train,
                               matrix_Xj_continuous_train,
                               matrix_Xi_unordered_train,
                               matrix_Xi_ordered_train,
                               matrix_Xi_continuous_train,
                               NULL,
                               NULL,
                               NULL,
                               vsfx,
                               bandwidth_provided,
                               matrix_bandwidth_xi,
                               matrix_bandwidth_xj,
                               lambdax,
                               num_categories_extern_X,
                               matrix_categorical_vals_extern_X,
                               NULL,
                               NULL, // no mean
                               NULL, // no permutations
                               kx_ij,
                               &gate_x_ctx);
      }

      for(iwk = 0; iwk < nwk; iwk++){
        const int64_t wko = iwk*wk;
        const int64_t dwk = (iwk != (nwk - 1)) ? wk : num_obs_train - (nwk - 1)*wk;

        idxk[0] = wko;
        idxk[1] = wko + dwk -1;

        kx_ik = pkx_ik;

        for(l = 0; l < num_reg_continuous; l++)
          matrix_Xk_continuous_train[l] = matrix_XY_continuous_train[l] + wko;

        for(l = 0; l < num_reg_unordered; l++)
          matrix_Xk_unordered_train[l] = matrix_XY_unordered_train[l] + wko;

        for(l = 0; l < num_reg_ordered; l++)
          matrix_Xk_ordered_train[l] = matrix_XY_ordered_train[l] + wko;


        for(l = 0; l < num_var_continuous; l++)
          matrix_Yk_continuous_train[l] = matrix_XY_continuous_train[l+num_reg_continuous] + wko;

        for(l = 0; l < num_var_unordered; l++)
          matrix_Yk_unordered_train[l] = matrix_XY_unordered_train[l+num_reg_unordered] + wko;

        for(l = 0; l < num_var_ordered; l++)
          matrix_Yk_ordered_train[l] = matrix_XY_ordered_train[l+num_reg_ordered] + wko;

        // offset the appropriate bandwidths
        if(BANDWIDTH_den != BW_FIXED){
          for(l = 0; l < num_var_continuous; l++)
            matrix_bandwidth_yk[l] = matrix_bandwidth_y[l] + wko;

          if(BANDWIDTH_den == BW_GEN_NN)
            for(l = 0; l < num_reg_continuous; l++)
              matrix_bandwidth_xk[l] = matrix_bandwidth_x[l] + wko;
        }
        // compute block kx_ik

        if (iwk != iwj) {
          if(gate_x_all_large_fixed){
            for(i = 0; i < (dwi*dwk); i++)
              kx_ik[i] = x_all_large_fixed_const;
          } else {
            np_activate_bounds_x();
            kernel_weighted_sum_np_ctx(kernel_cx,
                                   kernel_ux,
                                   kernel_ox,
                                   BANDWIDTH_den,
                                   dwk,
                                   dwi,
                                   num_reg_unordered,
                                   num_reg_ordered,
                                   num_reg_continuous,
                                   0, // (do not) compute the leave-one-out marginals
                                   0,
                                   1, // kpow
                                   1, // bw divide
                                   1, // divide weights
                                   0,
                                   0,
                                   0,
                                   0,
                                   x_operator,
                                   OP_NOOP, // no permutations
                                   0, // no score
                                   0, // no ocg
                                   NULL,
                                   0, // don't explicity suppress parallel
                                   0,
                                   0,
                                   int_TREE_XY,
                                   1,
                                   kdt_extern_XY,
                                   &nls, 
                                   xyd,
                                   idxk,
                                   matrix_Xk_unordered_train,
                                   matrix_Xk_ordered_train,
                                   matrix_Xk_continuous_train,
                                   matrix_Xi_unordered_train,
                                   matrix_Xi_ordered_train,
                                   matrix_Xi_continuous_train,
                                   NULL,
                                   NULL,
                                   NULL,
                                   vsfx,
                                   bandwidth_provided,
                                   matrix_bandwidth_xi,
                                   matrix_bandwidth_xk,
                                   lambdax,
                                   num_categories_extern_X,
                                   matrix_categorical_vals_extern_X,
                                   NULL,
                                   NULL, // no mean
                                   NULL, // no permutations
                                   kx_ik,
                                   &gate_x_ctx);
          }

        } else {
          kx_ik = kx_ij;
        }
        // compute block ky_jk

        np_activate_bounds_y();
        kernel_weighted_sum_np_ctx(kernel_cy,
                               kernel_uy,
                               kernel_oy,
                               BANDWIDTH_den,
                               dwk,
                               dwj,
                               num_var_unordered,
                               num_var_ordered,
                               num_var_continuous,
                               0, // (do not) compute the leave-one-out marginals
                               0,
                               1,
                               1, // divide bw for the convolution kernels 
                               1, // divide weights 
                               0,
                               0,
                               0,
                               0,
                               y_operator,
                               OP_NOOP, // no permutations
                               0, // no score
                               0, // no ocg
                               NULL,
                               0, // don't explicity suppress parallel
                               0,
                               0,
                               int_TREE_XY,
                               1,
                               kdt_extern_XY,
                               &nls, 
                               xyd + num_reg_continuous,
                               idxk,
                               matrix_Yk_unordered_train,
                               matrix_Yk_ordered_train,
                               matrix_Yk_continuous_train,
                               matrix_Yj_unordered_train,
                               matrix_Yj_ordered_train,
                               matrix_Yj_continuous_train,
                               NULL,
                               NULL,
                               NULL,
                               vsfy,
                               bandwidth_provided,
                               matrix_bandwidth_yj,
                               matrix_bandwidth_yk,
                               lambday,
                               num_categories_extern_Y,
                               matrix_categorical_vals_extern_Y,
                               NULL,
                               NULL, // no mean
                               NULL, // no permutations
                               ky_jk,
                               &gate_y_ctx);

        if(np_den_cv_use_tree_bypass_path(gate_x_all_large_fixed, int_TREE_XY, BANDWIDTH_den)){
          const int64_t ie_dwi = MIN(ie,dwi);
          if(BANDWIDTH_den != BW_ADAP_NN){
            for(i = is+wio; i < (wio+ie_dwi); i++){
              const int64_t j_end = wjo + dwj;
              const int64_t j_stop_before_i = MIN(j_end, i);
              const int64_t j_start_after_i = MAX(wjo, i + 1);
              const int64_t k_end = wko + dwk;
              const int64_t k_stop_before_i = MIN(k_end, i);
              const int64_t k_start_after_i = MAX(wko, i + 1);
              tcvj = 0.0;

              for(j = wjo; j < j_stop_before_i; j++){
                const double tkxij = kx_ij[(i-wio)*dwj + j-wjo];
                if (tkxij == 0.0) continue;
                tcvk = 0.0;
                for(k = wko; k < k_stop_before_i; k++)
                  tcvk += kx_ik[(i-wio)*dwk + k-wko]*ky_jk[(j-wjo)*dwk + k-wko];
                for(k = k_start_after_i; k < k_end; k++)
                  tcvk += kx_ik[(i-wio)*dwk + k-wko]*ky_jk[(j-wjo)*dwk + k-wko];
                tcvj += tkxij*tcvk;
              }

              for(j = j_start_after_i; j < j_end; j++){
                const double tkxij = kx_ij[(i-wio)*dwj + j-wjo];
                if (tkxij == 0.0) continue;
                tcvk = 0.0;
                for(k = wko; k < k_stop_before_i; k++)
                  tcvk += kx_ik[(i-wio)*dwk + k-wko]*ky_jk[(j-wjo)*dwk + k-wko];
                for(k = k_start_after_i; k < k_end; k++)
                  tcvk += kx_ik[(i-wio)*dwk + k-wko]*ky_jk[(j-wjo)*dwk + k-wko];
                tcvj += tkxij*tcvk;
              }
              *cv += tcvj/(mean[i]*mean[i] + DBL_MIN);
            }
          } else {
            for(i = is+wio; i < (wio+ie_dwi); i++){
              const int64_t j_end = wjo + dwj;
              const int64_t j_stop_before_i = MIN(j_end, i);
              const int64_t j_start_after_i = MAX(wjo, i + 1);
              const int64_t k_end = wko + dwk;
              const int64_t k_stop_before_i = MIN(k_end, i);
              const int64_t k_start_after_i = MAX(wko, i + 1);
              tcvj = 0.0;

              for(j = wjo; j < j_stop_before_i; j++){
                const double tkxij = kx_ij[(j-wjo)*dwi + i-wio];
                if (tkxij == 0.0) continue;
                tcvk = 0.0;
                for(k = wko; k < k_stop_before_i; k++)
                  tcvk += kx_ik[(k-wko)*dwi + i-wio]*ky_jk[(j-wjo)*dwk + k-wko];
                for(k = k_start_after_i; k < k_end; k++)
                  tcvk += kx_ik[(k-wko)*dwi + i-wio]*ky_jk[(j-wjo)*dwk + k-wko];
                tcvj += tkxij*tcvk;
              }

              for(j = j_start_after_i; j < j_end; j++){
                const double tkxij = kx_ij[(j-wjo)*dwi + i-wio];
                if (tkxij == 0.0) continue;
                tcvk = 0.0;
                for(k = wko; k < k_stop_before_i; k++)
                  tcvk += kx_ik[(k-wko)*dwi + i-wio]*ky_jk[(j-wjo)*dwk + k-wko];
                for(k = k_start_after_i; k < k_end; k++)
                  tcvk += kx_ik[(k-wko)*dwi + i-wio]*ky_jk[(j-wjo)*dwk + k-wko];
                tcvj += tkxij*tcvk;
              }
              *cv += tcvj/(mean[i]*mean[i] + DBL_MIN);
              //              Rprintf("i, cv: %d, %3.15g \n",i, *cv);
            }
          }
        } else {
          const int64_t ie_dwi = MIN(ie,dwi);

          for(i = is+wio; i < (wio+ie_dwi); i++){

            // reset interaction lists
            xl_xij.n = 0;
            xl_xik.n = 0;

            for(l = 0; l < num_reg_continuous; l++){
              const double tbw = (BANDWIDTH_den == BW_FIXED) ? vsfx[l] : matrix_bandwidth_x[l][i];

              bb[2*l] = -cksup[KERNEL_XY[l]][1];
              bb[2*l+1] = -cksup[KERNEL_XY[l]][0];

              bb[2*l] = (fabs(bb[2*l]) == DBL_MAX) ? bb[2*l] : (matrix_Xi_continuous_train[l][i-wio] + bb[2*l]*tbw);
              bb[2*l+1] = (fabs(bb[2*l+1]) == DBL_MAX) ? bb[2*l+1] : (matrix_Xi_continuous_train[l][i-wio] + bb[2*l+1]*tbw);
            }

            boxSearchNLPartialIdx(kdt_extern_XY, &nls, bb, NULL, &xl_xij, xyd, num_reg_continuous, idxj);

            if((idxj[0] != idxk[0]) || (idxj[1] != idxk[1])){
              boxSearchNLPartialIdx(kdt_extern_XY, &nls, bb, NULL, &xl_xik, xyd, num_reg_continuous, idxk);
            } else {
              mirror_xl(&xl_xij,&xl_xik);
            }

            tcvj = 0.0;
            for (m_ij = 0; m_ij < xl_xij.n; m_ij++){
              // j is offset by wjo!
              const int64_t jstart = xl_xij.istart[m_ij];
              const int64_t jnlev = xl_xij.nlev[m_ij];

              for(j = jstart; j < (jstart + jnlev); j++){              
                tcvk = 0.0;
                if((j+wjo) == i){
                  continue;
                }

                const double tkxij = kx_ij[(i-wio)*dwj + j];
                if (tkxij != 0.0) {
                  for (m_ik = 0; m_ik < xl_xik.n; m_ik++){
                    const int64_t kstart = xl_xik.istart[m_ik];
                    const int64_t knlev = xl_xik.nlev[m_ik];

                    for(k = kstart; k < (kstart + knlev); k++){
                      if((k+wko) == i) continue;
                      tcvk += kx_ik[(i-wio)*dwk + k]*ky_jk[j*dwk + k];
                    }
                  }
                  tcvj += tkxij*tcvk;
                }
              }
            }
            const int tixy = ipt_lookup_extern_X[ipt_extern_XY[i]];
            *cv += tcvj/(mean[tixy]*mean[tixy] + DBL_MIN);
          }

        }
      }
    }
  }
  }


#ifdef MPI2
  MPI_Allreduce(MPI_IN_PLACE, cv, 1, MPI_DOUBLE, MPI_SUM, comm[1]);
#endif

  *cv /= (double)num_obs_train;

  free(mean);

  free(pkx_ij);
  free(pkx_ik);
  free(pky_jk);

  free(x_operator);
  free(y_operator);
  free(xy_operator);

  free(x_cont_ok);
  free(x_disc_uno_ok);
  free(x_disc_ord_ok);
  free(x_cont_hmin);
  free(x_cont_k0);
  free(x_disc_uno_const);
  free(x_disc_ord_const);

  free(y_cont_ok);
  free(y_disc_uno_ok);
  free(y_disc_ord_ok);
  free(y_cont_hmin);
  free(y_cont_k0);
  free(y_disc_uno_const);
  free(y_disc_ord_const);

  free(xy_cont_ok);
  free(xy_disc_uno_ok);
  free(xy_disc_ord_ok);
  free(xy_cont_hmin);
  free(xy_cont_k0);
  free(xy_disc_uno_const);
  free(xy_disc_ord_const);

  free(kernel_cx);
  free(kernel_cy);
  free(kernel_cxy);

  free(kernel_ox);
  free(kernel_oy);
  free(kernel_oxy);

  free(kernel_ux);
  free(kernel_uy);
  free(kernel_uxy);

  if(BANDWIDTH_den != BW_FIXED){
    if(BANDWIDTH_den == BW_ADAP_NN){
      free(matrix_bandwidth_xi);
    }else{
      free(matrix_bandwidth_xj);      
      free(matrix_bandwidth_xk);      
    }

      free(matrix_bandwidth_yj);
      free(matrix_bandwidth_yk);

      free_tmat(matrix_bandwidth_x);
      free_tmat(matrix_bandwidth_xy);
      free_tmat(matrix_bandwidth_y);
      free(lambdax);
      free(lambday);
      free(lambdaxy);
  }else{
    free(matrix_bandwidth_x);
  }

  free(matrix_Xi_continuous_train);
  free(matrix_Xi_unordered_train);
  free(matrix_Xi_ordered_train);

  free(matrix_Xj_continuous_train);
  free(matrix_Xj_unordered_train);
  free(matrix_Xj_ordered_train);

  free(matrix_Xk_continuous_train);
  free(matrix_Xk_unordered_train);
  free(matrix_Xk_ordered_train);

  free(matrix_Yj_continuous_train);
  free(matrix_Yj_unordered_train);
  free(matrix_Yj_ordered_train);

  free(matrix_Yk_continuous_train);
  free(matrix_Yk_unordered_train);
  free(matrix_Yk_ordered_train);

  clean_nl(&nls);

  clean_xl(&xl_xij);
  clean_xl(&xl_xik);

  np_gate_ctx_clear(&gate_x_ctx);
  np_gate_ctx_clear(&gate_y_ctx);
  np_gate_ctx_clear(&gate_xy_ctx);

  return(0);
}

// estimation functions

int kernel_estimate_regression_categorical_tree_np(
int int_ll,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_reg,
int num_obs_train,
int num_obs_eval,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_Y,
double *vector_Y_eval,
double *vector_scale_factor,
int *num_categories,
double ** matrix_categorical_vals,
double * const restrict mean,
double **gradient,
double * const restrict mean_stderr,
double **gradient_stderr,
double *R_squared,
double *MSE,
double *MAE,
double *MAPE,
double *CORR,
double *SIGN){

  // note that mean has 2*num_obs allocated for npksum
  int i, j, l, sf_flag = 0;

	double INT_KERNEL_P;					 /* Integral of K(z)^2 */
	double K_INT_KERNEL_P;				 /*  K^p */
	/* Integral of K(z-0.5)*K(z+0.5) */
	double INT_KERNEL_PM_HALF = 0.0;
	double DIFF_KER_PPM = 0.0;		 /* Difference between int K(z)^p and int K(z-.5)K(z+.5) */
  double hprod;

  double * lambda = NULL, * vsf = NULL;
  double ** matrix_bandwidth = NULL;
  double ** matrix_bandwidth_deriv = NULL;
  int * operator = NULL;
  NP_GateOverrideCtx gate_ctx_local;
  int gate_override_active = 0;
  int *ov_cont_ok = NULL, *ov_disc_uno_ok = NULL, *ov_disc_ord_ok = NULL;
  double *ov_cont_hmin = NULL, *ov_cont_k0 = NULL;
  double *ov_disc_uno_const = NULL, *ov_disc_ord_const = NULL;
  int ov_cont_from_cache = 0;
  int estimation_shortcut_done = 0;

  operator = (int *)malloc(sizeof(int)*(num_reg_continuous+num_reg_unordered+num_reg_ordered));

  for(i = 0; i < (num_reg_continuous+num_reg_unordered+num_reg_ordered); i++)
    operator[i] = OP_NORMAL;

#ifdef MPI2
  int stride_e = MAX((int)ceil((double) num_obs_eval / (double) iNum_Processors),1);
  int num_obs_eval_alloc = stride_e*iNum_Processors;

#else
  int num_obs_eval_alloc = num_obs_eval;
#endif

  const int do_grad = (gradient != NULL); 
  const int do_gerr = (gradient_stderr != NULL);
  const int int_ll_est =
    np_reg_use_canonical_glp_degree1_estimation(int_ll,
                                                BANDWIDTH_reg,
                                                num_reg_continuous) ? LL_LL : int_ll;
  np_gate_ctx_clear(&gate_ctx_local);
  const NP_GateOverrideCtx * const est_gate_ctx_ptr = &gate_ctx_local;

  struct th_table * otabs = NULL;
  struct th_entry * ret = NULL;
  int ** matrix_ordered_indices = NULL;

  const int bwmdim = (BANDWIDTH_reg==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_train:1);
  const int fit_progress_total =
    (BANDWIDTH_reg == BW_ADAP_NN) ? num_obs_train : num_obs_eval;
  const int fit_progress_active = np_progress_fit_is_active();

  int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;

  kernel_c = (int *)malloc(sizeof(int)*num_reg_continuous);

  for(i = 0; i < num_reg_continuous; i++)
    kernel_c[i] = KERNEL_reg;

  kernel_u = (int *)malloc(sizeof(int)*num_reg_unordered);

  for(i = 0; i < num_reg_unordered; i++)
    kernel_u[i] = KERNEL_unordered_reg;

  kernel_o = (int *)malloc(sizeof(int)*num_reg_ordered);

  for(i = 0; i < num_reg_ordered; i++)
    kernel_o[i] = KERNEL_ordered_reg;


  // assert(BANDWIDTH_reg == BW_FIXED);

  // Allocate memory for objects 

  lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);
  matrix_bandwidth = alloc_tmatd(bwmdim,num_reg_continuous);

  matrix_bandwidth_deriv = alloc_matd(bwmdim,num_reg_continuous);

  if(kernel_bandwidth(KERNEL_reg,
                      BANDWIDTH_reg,
                      num_obs_train,
                      num_obs_eval,
                      0,
                      0,
                      0,
                      num_reg_continuous,
                      num_reg_unordered,
                      num_reg_ordered,
                      vector_scale_factor,
                      NULL,			 // Not used 
                      NULL,			 // Not used 
                      matrix_X_continuous_train,
                      matrix_X_continuous_eval,
                      NULL,					 // Not used 
	                      matrix_bandwidth,
	                      lambda,
	                      matrix_bandwidth_deriv)==1){
	    error("\n** Error: invalid bandwidth.");
	  }

  for(l = 0, hprod = 1.0; l < num_reg_continuous; l++)
    hprod *= matrix_bandwidth[l][0];


  if(num_reg_continuous != 0) {
    initialize_kernel_regression_asymptotic_constants(KERNEL_reg,
                                                      num_reg_continuous,
                                                      &INT_KERNEL_P,
                                                      &K_INT_KERNEL_P,
                                                      &INT_KERNEL_PM_HALF,
                                                      &DIFF_KER_PPM);
  } else {
    INT_KERNEL_P = 1.0;
    K_INT_KERNEL_P = 1.0;
  }

  const double gfac = sqrt(DIFF_KER_PPM/K_INT_KERNEL_P);

  // compute hash stuff here if necessary

  if(do_grad && (num_reg_ordered > 0)){
    otabs = (struct th_table *)malloc(num_reg_ordered*sizeof(struct th_table));
    matrix_ordered_indices = (int **)malloc(num_reg_ordered*sizeof(int *));
    int * tc = (int *)malloc(num_reg_ordered*num_obs_eval*sizeof(int));
    for(l = 0; l < num_reg_ordered; l++)
      matrix_ordered_indices[l] = tc + l*num_obs_eval;

    for(l = 0; l < num_reg_ordered; l++){
      if(thcreate_r((size_t)ceil(1.6*num_categories[l+num_reg_unordered]), otabs + l) == TH_ERROR)
        error("hash table creation failed");

      for(i = 0; i < num_categories[l+num_reg_unordered]; i++){
        struct th_entry centry;
        centry.key.dkey = matrix_categorical_vals[l+num_reg_unordered][i];
        centry.data = i;

        if(thsearch_r(&centry, TH_ENTER, &ret, otabs+l) == TH_FAILURE)
          error("insertion into hash table failed");
      }

      // now do lookups
      struct th_entry te;
      te.key.dkey = 0.0;
      te.data = -1;
      ret = &te;

      for(i = 0; i < num_obs_eval; i++){
        if(ret->key.dkey != matrix_X_ordered_eval[l][i]){
          te.key.dkey = matrix_X_ordered_eval[l][i];
          if(thsearch_r(&te, TH_SEARCH, &ret, otabs+l) == TH_FAILURE)
            error("hash table lookup failed (which should be impossible)");
        } 

        matrix_ordered_indices[l][i] = ret->data;
      }
    }
  }

  if((num_reg_continuous + num_reg_unordered + num_reg_ordered) > 0){
    int ok_all = 1;

    if(num_reg_continuous > 0){
      const double rel_tol = np_cont_largeh_rel_tol();
      if(np_cont_largeh_cache_get_or_build(num_obs_train,
                                           num_obs_eval,
                                           num_reg_continuous,
                                           kernel_c,
                                           matrix_X_continuous_train,
                                           matrix_X_continuous_eval,
                                           rel_tol,
                                           &ov_cont_ok,
                                           &ov_cont_hmin,
                                           &ov_cont_k0)){
        ov_cont_from_cache = 1;
      } else {
        ov_cont_ok = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
        ov_cont_hmin = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
        ov_cont_k0 = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
        ok_all = (ov_cont_ok != NULL) && (ov_cont_hmin != NULL) && (ov_cont_k0 != NULL);
        if(ok_all){
          for(i = 0; i < num_reg_continuous; i++){
            const int kern = kernel_c[i];
            double xmin = DBL_MAX, xmax = -DBL_MAX;
            ov_cont_ok[i] = 0; ov_cont_hmin[i] = DBL_MAX; ov_cont_k0[i] = 0.0;
            if(!np_cont_largeh_kernel_supported(kern)) continue;
            for(j = 0; j < num_obs_train; j++){
              const double v = matrix_X_continuous_train[i][j];
              if(!isfinite(v)) continue;
              xmin = MIN(xmin, v); xmax = MAX(xmax, v);
            }
            if(xmax >= xmin){
              const double utol = np_cont_largeh_utol(kern, rel_tol);
              if(utol > 0.0 && isfinite(utol)){
                ov_cont_ok[i] = 1;
                ov_cont_hmin[i] = (xmax - xmin)/utol;
                ov_cont_k0[i] = np_cont_largeh_k0(kern);
              }
            }
          }
        }
      }
    }

    if(ok_all && num_reg_unordered > 0){
      ov_disc_uno_ok = (int *)calloc((size_t)num_reg_unordered, sizeof(int));
      ov_disc_uno_const = (double *)malloc((size_t)num_reg_unordered*sizeof(double));
      ok_all = (ov_disc_uno_ok != NULL) && (ov_disc_uno_const != NULL);
      if(ok_all){
        double (* const ukf[])(int, double, int) = {
          np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
          np_score_uaa, np_score_unli_racine
        };
        const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
        for(i = 0; i < num_reg_unordered; i++){
          const int ku = kernel_u[i];
          const int ncat = (num_categories != NULL) ? num_categories[i] : 0;
          const double lam = lambda[i];
          ov_disc_uno_ok[i] = 0; ov_disc_uno_const[i] = 0.0;
          if(ku < 0 || ku >= nuk) continue;
          if(!np_disc_near_upper(ku, lam, ncat)) continue;
          {
            const double ks = ukf[ku](1, lam, ncat);
            const double kd = ukf[ku](0, lam, ncat);
            if(np_disc_near_const_kernel(ks, kd)){
              ov_disc_uno_ok[i] = 1;
              ov_disc_uno_const[i] = 0.5*(ks + kd);
            }
          }
        }
      }
    }

    if(ok_all && num_reg_ordered > 0){
      ov_disc_ord_ok = (int *)calloc((size_t)num_reg_ordered, sizeof(int));
      ov_disc_ord_const = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
      ok_all = (ov_disc_ord_ok != NULL) && (ov_disc_ord_const != NULL);
      if(ok_all){
        double (* const okf[])(double, double, double, double, double) = {
          np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
          np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
          np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
          np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
        };
        const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
        for(i = 0; i < num_reg_ordered; i++){
          const int oi = i + num_reg_unordered;
          const int ko = kernel_o[i];
          const int ncat = (num_categories != NULL) ? num_categories[oi] : 0;
          const double lam = lambda[oi];
          ov_disc_ord_ok[i] = 0; ov_disc_ord_const[i] = 0.0;
          if(ko < 0 || ko >= nok) continue;
          if(ncat <= 0 || matrix_categorical_vals == NULL) continue;
          if(!np_disc_ordered_near_upper(ko, lam)) continue;
          {
            const double cl = matrix_categorical_vals[oi][0];
            const double ch = matrix_categorical_vals[oi][ncat - 1];
            const double k0 = okf[ko](cl, cl, lam, cl, ch);
            const double k1 = okf[ko](cl, ch, lam, cl, ch);
            if(np_disc_near_const_kernel(k0, k1)){
              ov_disc_ord_ok[i] = 1;
              ov_disc_ord_const[i] = 0.5*(k0 + k1);
            }
          }
        }
      }
    }

    if(ok_all){
      np_gate_ctx_set(&gate_ctx_local,
                      num_reg_continuous,
                      num_reg_unordered,
                      num_reg_ordered,
                      kernel_c,
                      kernel_u,
                      kernel_o,
                      operator,
                      ov_cont_ok,
                      ov_cont_hmin,
                      ov_cont_k0,
                      ov_disc_uno_ok,
                      ov_disc_uno_const,
                      ov_disc_ord_ok,
                      ov_disc_ord_const);
      gate_override_active = 1;
    }
  }

  {
    int all_large_gate = (BANDWIDTH_reg == BW_FIXED) && gate_override_active;
    if(all_large_gate){
      for(i = 0; i < num_reg_continuous; i++){
        const double h = matrix_bandwidth[i][0];
        if((ov_cont_ok == NULL) || (!ov_cont_ok[i]) || (!isfinite(h)) ||
           (fabs(h) < ov_cont_hmin[i])){
          all_large_gate = 0;
          break;
        }
      }
    }
    if(all_large_gate){
      for(i = 0; i < num_reg_unordered; i++){
        if((ov_disc_uno_ok == NULL) || (!ov_disc_uno_ok[i])){
          all_large_gate = 0;
          break;
        }
      }
    }
    if(all_large_gate){
      for(i = 0; i < num_reg_ordered; i++){
        if((ov_disc_ord_ok == NULL) || (!ov_disc_ord_ok[i])){
          all_large_gate = 0;
          break;
        }
      }
    }

    if(all_large_gate &&
       ((int_ll_est == LL_LC) || (int_ll_est == LL_LL) || (int_ll_est == LL_LP))){
      double kconst = 1.0;
      int kconst_ok = 1;
      const double ridge_eps = 1.0/(double)MAX(1, num_obs_train);
      double sigma2hat = 0.0;
      const double ymean = meand(num_obs_train, vector_Y);

      for(i = 0; i < num_obs_train; i++){
        const double dy = vector_Y[i] - ymean;
        sigma2hat += dy*dy;
      }
      sigma2hat /= (double)MAX(1, num_obs_train);

      for(i = 0; i < num_reg_continuous; i++){
        const double h = matrix_bandwidth[i][0];
        if(!isfinite(h) || (h == 0.0)){
          kconst_ok = 0;
          break;
        }
        kconst *= ov_cont_k0[i]/h;
      }
      if(kconst_ok){
        for(i = 0; i < num_reg_unordered; i++)
          kconst *= ov_disc_uno_const[i];
      }
      if(kconst_ok){
        for(i = 0; i < num_reg_ordered; i++)
          kconst *= ov_disc_ord_const[i];
      }
      if(!isfinite(kconst) || (kconst <= 0.0))
        kconst_ok = 0;

      if(kconst_ok && (int_ll_est == LL_LC) && (!do_grad)){
        const double sk = ((double)num_obs_train)*kconst;
        const double sefac = (sk*hprod > 0.0) ?
          sqrt(MAX(0.0, sigma2hat) * K_INT_KERNEL_P / (sk*hprod)) : 0.0;

        if (fit_progress_active) {
          for(i = 0; i < num_obs_eval; i++){
            mean[i] = ymean;
            mean_stderr[i] = sefac;
            np_progress_fit_loop_step(i + 1, fit_progress_total);
          }
        } else {
          for(i = 0; i < num_obs_eval; i++){
            mean[i] = ymean;
            mean_stderr[i] = sefac;
          }
        }

        estimation_shortcut_done = 1;
      } else if(kconst_ok && int_ll_est == LL_LL){
        const int k = num_reg_continuous + 1;
        MATRIX XtX = mat_creat(k, k, UNDEFINED);
        MATRIX XtXINV = mat_creat(k, k, UNDEFINED);
        MATRIX XtY = mat_creat(k, 1, UNDEFINED);
        MATRIX BETA = mat_creat(k, 1, UNDEFINED);
        int fast_ok = (XtX != NULL) && (XtXINV != NULL) && (XtY != NULL) && (BETA != NULL);

        if(fast_ok){
          int ridge_it = 0;
          const double sk = ((double)num_obs_train)*kconst;
          const double sefac = (sk*hprod > 0.0) ? sqrt(MAX(0.0, sigma2hat) * K_INT_KERNEL_P / (sk*hprod)) : 0.0;

          for(i = 0; i < k; i++){
            XtY[i][0] = 0.0;
            BETA[i][0] = 0.0;
            for(j = 0; j < k; j++)
              XtX[i][j] = 0.0;
          }

          for(i = 0; i < num_obs_train; i++){
            const double yi = vector_Y[i];
            XtX[0][0] += 1.0;
            XtY[0][0] += yi;
            for(j = 0; j < num_reg_continuous; j++){
              const double xj = matrix_X_continuous_train[j][i];
              const int cj = j + 1;
              XtX[0][cj] += xj;
              XtX[cj][0] += xj;
              XtY[cj][0] += xj*yi;
            }
            for(int a = 0; a < num_reg_continuous; a++){
              const double xa = matrix_X_continuous_train[a][i];
              const int ca = a + 1;
              for(int b = a; b < num_reg_continuous; b++){
                const double xb = matrix_X_continuous_train[b][i];
                const int cb = b + 1;
                XtX[ca][cb] += xa*xb;
                if(cb != ca) XtX[cb][ca] += xa*xb;
              }
            }
          }

          while(mat_inv(XtX, XtXINV) == NULL){
            for(i = 0; i < k; i++)
              XtX[i][i] += ridge_eps;
            ridge_it++;
            if(ridge_it > 64){
              fast_ok = 0;
              break;
            }
          }

          if(fast_ok){
            for(i = 0; i < k; i++){
              double s = 0.0;
              for(j = 0; j < k; j++)
                s += XtXINV[i][j]*XtY[j][0];
              BETA[i][0] = s;
            }

            for(i = 0; i < num_obs_eval; i++){
              double yhat = BETA[0][0];
              for(j = 0; j < num_reg_continuous; j++)
                yhat += BETA[j+1][0]*matrix_X_continuous_eval[j][i];
              mean[i] = yhat;
              mean_stderr[i] = sefac;
              if (fit_progress_active)
                np_progress_fit_loop_step(i + 1, fit_progress_total);
            }

            if(do_grad){
              const int nvars = num_reg_continuous + num_reg_unordered + num_reg_ordered;
              for(j = 0; j < num_reg_continuous; j++){
                const double bj = BETA[j+1][0];
                for(i = 0; i < num_obs_eval; i++){
                  gradient[j][i] = bj;
                  if(do_gerr){
                    gradient_stderr[j][i] = gfac*mean_stderr[i]/
                      ((BANDWIDTH_reg == BW_ADAP_NN) ? 1.0 :
                       ((BANDWIDTH_reg == BW_GEN_NN) ? matrix_bandwidth[j][i] : matrix_bandwidth[j][0]));
                  }
                }
              }

              for(j = num_reg_continuous; j < nvars; j++){
                for(i = 0; i < num_obs_eval; i++){
                  gradient[j][i] = 0.0;
                  if(do_gerr) gradient_stderr[j][i] = 0.0;
                }
              }
            }
            estimation_shortcut_done = 1;
          }
        }

        if(XtX != NULL) mat_free(XtX);
        if(XtXINV != NULL) mat_free(XtXINV);
        if(XtY != NULL) mat_free(XtY);
        if(BETA != NULL) mat_free(BETA);
      } else if(kconst_ok && int_ll_est == LL_LP &&
                (vector_glp_degree_extern != NULL) && (num_reg_continuous > 0)){
        const int use_bernstein = (int_glp_bernstein_extern != 0);
        int *glp_terms = NULL;
        int glp_nterms = 0;
        double **basis = NULL;
        NPGLPBasisCtx *basis_ctx = NULL;
        double *eval_basis = NULL, *eval_deriv = NULL;
        MATRIX XtX = NULL, XtXINV = NULL, XtY = NULL, BETA = NULL;
        int fast_ok = np_glp_build_terms(num_reg_continuous,
                                         vector_glp_degree_extern,
                                         int_glp_basis_extern,
                                         &glp_terms,
                                         &glp_nterms);
        if(fast_ok && (glp_nterms > 0)){
          basis = alloc_matd(num_obs_train, glp_nterms);
          XtX = mat_creat(glp_nterms, glp_nterms, UNDEFINED);
          XtXINV = mat_creat(glp_nterms, glp_nterms, UNDEFINED);
          XtY = mat_creat(glp_nterms, 1, UNDEFINED);
          BETA = mat_creat(glp_nterms, 1, UNDEFINED);
          eval_basis = (double *)malloc((size_t)glp_nterms*sizeof(double));
          if(do_grad)
            eval_deriv = (double *)malloc((size_t)glp_nterms*sizeof(double));
          if(use_bernstein)
            basis_ctx = (NPGLPBasisCtx *)calloc((size_t)num_reg_continuous, sizeof(NPGLPBasisCtx));
          fast_ok = (basis != NULL) && (XtX != NULL) && (XtXINV != NULL) &&
            (XtY != NULL) && (BETA != NULL) && (eval_basis != NULL) &&
            ((!do_grad) || (eval_deriv != NULL)) &&
            (!use_bernstein || (basis_ctx != NULL));
        } else {
          fast_ok = 0;
        }

        if(fast_ok){
          for(l = 0; l < num_reg_continuous; l++){
            if(use_bernstein){
              double xmin = matrix_X_continuous_train[l][0];
              double xmax = matrix_X_continuous_train[l][0];
              for(i = 1; i < num_obs_train; i++){
                const double xi = matrix_X_continuous_train[l][i];
                if(xi < xmin) xmin = xi;
                if(xi > xmax) xmax = xi;
              }
              if(!np_glp_basis_ctx_init(&basis_ctx[l], vector_glp_degree_extern[l], xmin, xmax)){
                fast_ok = 0;
                break;
              }
            }
          }
        }

        if(fast_ok){
          if(use_bernstein){
            np_glp_fill_basis_train(num_reg_continuous,
                                    glp_terms,
                                    glp_nterms,
                                    matrix_X_continuous_train,
                                    num_obs_train,
                                    basis_ctx,
                                    basis);
          } else {
            np_glp_fill_basis_raw_train(num_reg_continuous,
                                        glp_terms,
                                        glp_nterms,
                                        matrix_X_continuous_train,
                                        num_obs_train,
                                        basis);
          }
        }

        if(fast_ok){
          int ridge_it = 0;
          const double sk = ((double)num_obs_train)*kconst;
          const double se_default = (sk*hprod > 0.0) ?
            sqrt(MAX(0.0, sigma2hat) * K_INT_KERNEL_P / (sk*hprod)) : 0.0;

          for(i = 0; i < glp_nterms; i++){
            XtY[i][0] = 0.0;
            BETA[i][0] = 0.0;
            for(j = 0; j < glp_nterms; j++)
              XtX[i][j] = 0.0;
          }

          for(i = 0; i < num_obs_train; i++){
            const double yi = vector_Y[i];
            for(int a = 0; a < glp_nterms; a++){
              const double za = basis[a][i];
              XtY[a][0] += za*yi;
              for(int b = a; b < glp_nterms; b++){
                const double zb = basis[b][i];
                XtX[a][b] += za*zb;
                if(b != a) XtX[b][a] += za*zb;
              }
            }
          }

          while(mat_inv(XtX, XtXINV) == NULL){
            for(i = 0; i < glp_nterms; i++)
              XtX[i][i] += ridge_eps;
            ridge_it++;
            if(ridge_it > 64){
              fast_ok = 0;
              break;
            }
          }

          if(fast_ok){
            for(i = 0; i < glp_nterms; i++){
              double s = 0.0;
              for(j = 0; j < glp_nterms; j++)
                s += XtXINV[i][j]*XtY[j][0];
              BETA[i][0] = s;
            }

            for(i = 0; i < num_obs_eval; i++){
              double yhat = 0.0;
              double q = 0.0;

              if(use_bernstein){
                np_glp_fill_basis_eval(num_reg_continuous,
                                       glp_terms,
                                       glp_nterms,
                                       matrix_X_continuous_eval,
                                       i,
                                       basis_ctx,
                                       eval_basis);
              } else {
                np_glp_fill_basis_eval_raw(num_reg_continuous,
                                           glp_terms,
                                           glp_nterms,
                                           matrix_X_continuous_eval,
                                           i,
                                           eval_basis);
              }

              for(j = 0; j < glp_nterms; j++)
                yhat += eval_basis[j]*BETA[j][0];
              for(j = 0; j < glp_nterms; j++){
                const double zj = eval_basis[j];
                for(int b = 0; b < glp_nterms; b++)
                  q += zj*XtXINV[j][b]*eval_basis[b];
              }

              mean[i] = yhat;
              {
                const double mv = sigma2hat*q;
                mean_stderr[i] = (mv > 0.0 && isfinite(mv)) ? sqrt(mv) : se_default;
              }

              if(do_grad){
                const int nvars = num_reg_continuous + num_reg_unordered + num_reg_ordered;
                for(l = 0; l < num_reg_continuous; l++){
                  const int grad_order =
                    (vector_glp_gradient_order_extern != NULL) ?
                    MAX(1, vector_glp_gradient_order_extern[l]) : 1;
                  double qg = 0.0;
                  double dg = 0.0;

                  if(use_bernstein){
                    np_glp_fill_basis_eval_deriv(l,
                                                 grad_order,
                                                 num_reg_continuous,
                                                 glp_terms,
                                                 glp_nterms,
                                                 matrix_X_continuous_eval,
                                                 i,
                                                 basis_ctx,
                                                 eval_deriv);
                  } else {
                    np_glp_fill_basis_eval_deriv_raw(l,
                                                     grad_order,
                                                     num_reg_continuous,
                                                     glp_terms,
                                                     glp_nterms,
                                                     matrix_X_continuous_eval,
                                                     i,
                                                     eval_deriv);
                  }

                  for(j = 0; j < glp_nterms; j++)
                    dg += eval_deriv[j]*BETA[j][0];
                  gradient[l][i] = dg;

                  if(do_gerr){
                    for(j = 0; j < glp_nterms; j++){
                      const double dj = eval_deriv[j];
                      for(int b = 0; b < glp_nterms; b++)
                        qg += dj*XtXINV[j][b]*eval_deriv[b];
                    }
                    {
                      const double gv = sigma2hat*qg;
                      gradient_stderr[l][i] = (gv > 0.0 && isfinite(gv)) ? sqrt(gv) : 0.0;
                    }
                  }
                }

                for(l = num_reg_continuous; l < nvars; l++){
                  gradient[l][i] = 0.0;
                  if(do_gerr) gradient_stderr[l][i] = 0.0;
                }
              }
              if (fit_progress_active)
                np_progress_fit_loop_step(i + 1, fit_progress_total);
            }
            estimation_shortcut_done = 1;
          }
        }

        if(use_bernstein && (basis_ctx != NULL)){
          for(l = 0; l < num_reg_continuous; l++) np_glp_basis_ctx_free(&basis_ctx[l]);
          free(basis_ctx);
        }
        if(eval_basis != NULL) free(eval_basis);
        if(eval_deriv != NULL) free(eval_deriv);
        if(BETA != NULL) mat_free(BETA);
        if(XtY != NULL) mat_free(XtY);
        if(XtXINV != NULL) mat_free(XtXINV);
        if(XtX != NULL) mat_free(XtX);
        if(basis != NULL) free_mat(basis, glp_nterms);
        if(glp_terms != NULL) free(glp_terms);
      }
    }
  }

  if(estimation_shortcut_done)
    goto finish_regression_estimation;

  // Conduct the estimation 

  if(int_ll_est == LL_LC) { // local constant
    // Nadaraya-Watson
    // Generate bandwidth vector given scale factors, nearest neighbors, or lambda 

#define NCOL_Y 3

    double * lc_Y[NCOL_Y] = {NULL,NULL,NULL};
    double * meany = (double *)malloc(NCOL_Y*num_obs_eval_alloc*sizeof(double));
    double * permy = NULL;
    int pop = OP_NOOP;
    int p_nvar = do_grad ? (num_reg_continuous + num_reg_unordered + num_reg_ordered) : 0;
  
    if(do_grad){
      permy = (double *)malloc(NCOL_Y*num_obs_eval_alloc*p_nvar*sizeof(double));
      if(permy == NULL)
        error("\n** Error: memory allocation failed.");
      pop = OP_DERIVATIVE;
    }
    
    if(meany == NULL)
      error("\n** Error: memory allocation failed.");

    lc_Y[0] = vector_Y;
      
    lc_Y[1] = (double *)malloc(num_obs_train*sizeof(double));
    lc_Y[2] = (double *)malloc(num_obs_train*sizeof(double));

    if((lc_Y[1] == NULL) || (lc_Y[2] == NULL))
      error("\n** Error: memory allocation failed.");

    for(int ii = 0; ii < num_obs_train; ii++){
      lc_Y[1][ii] = 1.0;
      lc_Y[2][ii] = vector_Y[ii]*vector_Y[ii];
    }
    
    kernel_weighted_sum_np_ctx(kernel_c,
                           kernel_u,
                           kernel_o,
                           BANDWIDTH_reg,
                           num_obs_train,
                           num_obs_eval,
                           num_reg_unordered,
                           num_reg_ordered,
                           num_reg_continuous,
                           0, // no leave one out 
                           0,
                           1, // kernel_pow = 1
                           1, // bandwidth_divide = TRUE, always
                           0, 
                           0, // not symmetric
                           0, // do not gather-scatter
                           0, // do not drop train
                           0, // do not drop train
                           operator, // no special operators being used
                           pop, // permutations used for gradients
                           0, // no score
                           do_grad, // ocg if grad 
                           NULL,
                           0, // don't explicity suppress parallel
                           NCOL_Y,
                           0,
                           int_TREE_X,
                           0,
                           kdt_extern_X,
                           NULL, NULL, NULL,
                           matrix_X_unordered_train, // TRAIN
                           matrix_X_ordered_train,
                           matrix_X_continuous_train,
                           matrix_X_unordered_eval, // EVAL
                           matrix_X_ordered_eval,
                           matrix_X_continuous_eval,
                           lc_Y,
                           NULL, // no W matrix
                           NULL, // no sgn 
                           vector_scale_factor,
                           1,
                           matrix_bandwidth,
                           matrix_bandwidth,
                           lambda,
                           num_categories,
                           matrix_categorical_vals,
                           matrix_ordered_indices, 
                           meany,
                           permy, // permutations used for gradients
                           NULL, // do not return kernel weights
                           est_gate_ctx_ptr);

    for(int ii = 0; ii < num_obs_eval; ii++){
      const int ii3 = NCOL_Y*ii;
      const double sk = copysign(DBL_MIN, meany[ii3+1]) + meany[ii3+1];
      mean[ii] = meany[ii3]/sk;

      mean_stderr[ii] = meany[ii3+2]/sk - mean[ii]*mean[ii];

      mean_stderr[ii] = (mean_stderr[ii] <= 0.0) ? 0.0 : sqrt(mean_stderr[ii] * K_INT_KERNEL_P / (sk*hprod));
    }
   
    if(do_grad){
      for(l = 0; l < num_reg_continuous; l++){
        for(i = 0; i < num_obs_eval; i++){
          // ordered y, 1, y*y
          const int ii3 = NCOL_Y*i;
          const int li3 = l*num_obs_eval*NCOL_Y + ii3;
          const double sk = copysign(DBL_MIN, meany[ii3+1]) + meany[ii3+1];
          gradient[l][i] = (permy[li3] - mean[i]*permy[li3+1])/sk;
          
          if(do_gerr){
            gradient_stderr[l][i] = gfac*mean_stderr[i]/((BANDWIDTH_reg == BW_ADAP_NN) ? 1.0 : ((BANDWIDTH_reg == BW_GEN_NN) ? matrix_bandwidth[l][i]:matrix_bandwidth[l][0]));
          }
        }
      }

      for(l = num_reg_continuous; l < (num_reg_continuous + num_reg_unordered); l++){
        for(i = 0; i < num_obs_eval; i++){
          const int ii3 = NCOL_Y*i;
          const int li3 = l*num_obs_eval*NCOL_Y + ii3;
          const double sk = copysign(DBL_MIN, permy[li3+1]) + permy[li3+1];
          const double s1 = permy[li3]/sk;

          gradient[l][i] = mean[i] - s1;
          
          if(do_gerr && (num_reg_continuous > 0)){
            const double se = permy[li3+2]/sk - s1*s1;
            const double senn = (se <= 0.0) ? 0.0 : se;
            gradient_stderr[l][i] = sqrt(mean_stderr[i]*mean_stderr[i] + senn*K_INT_KERNEL_P/(sk*hprod));
          } else {
            gradient_stderr[l][i] = 0.0;
          }
        }
      }

      for(l = num_reg_continuous + num_reg_unordered; l < p_nvar; l++){
        for(i = 0; i < num_obs_eval; i++){
          const int ii3 = NCOL_Y*i;
          const int li3 = l*num_obs_eval*NCOL_Y + ii3;
          const double sk = copysign(DBL_MIN, permy[li3+1]) + permy[li3+1];
          const double s1 = permy[li3]/sk;
          
          gradient[l][i] = (mean[i] - s1)*((matrix_ordered_indices[l - num_reg_continuous - num_reg_unordered][i] != 0) ? 1.0 : -1.0);

          if(do_gerr && (num_reg_continuous > 0)){
            const double se = permy[li3+2]/sk - s1*s1;
            const double senn = (se <= 0.0) ? 0.0 : se;
            gradient_stderr[l][i] = sqrt(mean_stderr[i]*mean_stderr[i] + senn*K_INT_KERNEL_P/(sk*hprod));
          } else {
            gradient_stderr[l][i] = 0.0;
          }
        }
      }
    }


    free(meany);

    if(do_grad && (p_nvar > 0)){
      free(permy);
    }

    free(lc_Y[1]);
    free(lc_Y[2]);
#undef NCOL_Y
  } else if(int_ll_est == LL_LP) { // local polynomial (regtype = "lp")
    int *glp_terms = NULL;
    int glp_nterms = 0;
    const int use_bernstein = (int_glp_bernstein_extern != 0);
    MATRIX KWM = NULL, XTKY = NULL, DELTA = NULL;
    MATRIX KWM2 = NULL, KWM_INV = NULL, IDEN = NULL;
    double **basis = NULL;
    double **TCON = NULL, **TUNO = NULL, **TORD = NULL;
    double **Ycols = NULL, **Wcols = NULL;
    double *y2 = NULL, *out = NULL, *out2 = NULL;
    NPGLPBasisCtx *basis_ctx = NULL;
    double *eval_basis = NULL, *eval_deriv = NULL;
    double *tmp_v = NULL, *tmp_w = NULL;
    const double epsilon = 1.0/(double)MAX(1, num_obs_train);

    if((vector_glp_degree_extern == NULL) || (num_reg_continuous <= 0))
      error("glp degree vector unavailable");

    if(!np_glp_build_terms(num_reg_continuous, vector_glp_degree_extern, int_glp_basis_extern, &glp_terms, &glp_nterms))
      error("failed to build glp basis terms");
    if(glp_nterms <= 0)
      error("invalid glp basis dimension");

    KWM = mat_creat(glp_nterms, glp_nterms, UNDEFINED);
    XTKY = mat_creat(glp_nterms, 1, UNDEFINED);
    DELTA = mat_creat(glp_nterms, 1, UNDEFINED);
    KWM2 = mat_creat(glp_nterms, glp_nterms, UNDEFINED);
    KWM_INV = mat_creat(glp_nterms, glp_nterms, UNDEFINED);
    IDEN = mat_creat(glp_nterms, glp_nterms, UNDEFINED);
    basis = alloc_matd(num_obs_train, glp_nterms);
    TCON = alloc_matd(1, num_reg_continuous);
    TUNO = alloc_matd(1, num_reg_unordered);
    TORD = alloc_matd(1, num_reg_ordered);
    Ycols = (double **)malloc((size_t)(glp_nterms + 2)*sizeof(double *));
    Wcols = (double **)malloc((size_t)glp_nterms*sizeof(double *));
    y2 = (double *)malloc((size_t)num_obs_train*sizeof(double));
    out = (double *)malloc((size_t)(glp_nterms + 2)*(size_t)glp_nterms*sizeof(double));
    out2 = (double *)malloc((size_t)glp_nterms*(size_t)glp_nterms*sizeof(double));
    tmp_v = (double *)malloc((size_t)glp_nterms*sizeof(double));
    tmp_w = (double *)malloc((size_t)glp_nterms*sizeof(double));
    eval_basis = (double *)malloc((size_t)glp_nterms*sizeof(double));
    eval_deriv = (double *)malloc((size_t)glp_nterms*sizeof(double));
    if(use_bernstein)
      basis_ctx = (NPGLPBasisCtx *)calloc((size_t)num_reg_continuous, sizeof(NPGLPBasisCtx));

    if(!((KWM != NULL) && (XTKY != NULL) && (DELTA != NULL) &&
      (KWM2 != NULL) && (KWM_INV != NULL) && (IDEN != NULL) &&
      (basis != NULL) &&
      ((num_reg_continuous == 0) || (TCON != NULL)) &&
      ((num_reg_unordered == 0) || (TUNO != NULL)) &&
      ((num_reg_ordered == 0) || (TORD != NULL)) &&
      (Ycols != NULL) && (Wcols != NULL) && (y2 != NULL) && (out != NULL) &&
      (out2 != NULL) && (tmp_v != NULL) && (tmp_w != NULL) &&
      (eval_basis != NULL) && (eval_deriv != NULL) &&
      (!use_bernstein || (basis_ctx != NULL))))
      error("memory allocation failed in glp path");

    for(i = 0; i < num_obs_train; i++)
      y2[i] = vector_Y[i]*vector_Y[i];

    if(use_bernstein){
      for(l = 0; l < num_reg_continuous; l++){
        double xmin = matrix_X_continuous_train[l][0];
        double xmax = matrix_X_continuous_train[l][0];
        for(i = 1; i < num_obs_train; i++){
          const double xi = matrix_X_continuous_train[l][i];
          if(xi < xmin) xmin = xi;
          if(xi > xmax) xmax = xi;
        }
        if(!np_glp_basis_ctx_init(&basis_ctx[l], vector_glp_degree_extern[l], xmin, xmax))
          error("failed to initialize glp Bernstein basis");
      }

      np_glp_fill_basis_train(num_reg_continuous,
                              glp_terms,
                              glp_nterms,
                              matrix_X_continuous_train,
                              num_obs_train,
                              basis_ctx,
                              basis);
    } else {
      np_glp_fill_basis_raw_train(num_reg_continuous,
                                  glp_terms,
                                  glp_nterms,
                                  matrix_X_continuous_train,
                                  num_obs_train,
                                  basis);
    }

    for(j = 0; j < num_obs_eval; j++){
      double nepsilon = 0.0;
      double sk, ey, ey2, sigma2hat;
      int have_vcov = 0;

      for(l = 0; l < num_reg_continuous; l++)
        TCON[l][0] = matrix_X_continuous_eval[l][j];
      for(l = 0; l < num_reg_unordered; l++)
        TUNO[l][0] = matrix_X_unordered_eval[l][j];
      for(l = 0; l < num_reg_ordered; l++)
        TORD[l][0] = matrix_X_ordered_eval[l][j];

      Ycols[0] = y2;
      Ycols[1] = vector_Y;
      for(l = 0; l < glp_nterms; l++){
        Ycols[l + 2] = basis[l];
        Wcols[l] = basis[l];
      }

      kernel_weighted_sum_np_ctx(kernel_c,
                             kernel_u,
                             kernel_o,
                             BANDWIDTH_reg,
                             num_obs_train,
                             1,
                             num_reg_unordered,
                             num_reg_ordered,
                             num_reg_continuous,
                             0,
                             0,
                             1,
                             1,
                             0,
                             0,
                             0,
                             0,
                             0,
                             operator,
                             OP_NOOP,
                             0,
                             0,
                             NULL,
                             1,
                             glp_nterms + 2,
                             glp_nterms,
                             (BANDWIDTH_reg == BW_ADAP_NN) ? NP_TREE_FALSE : int_TREE_X,
                             0,
                             (BANDWIDTH_reg == BW_ADAP_NN) ? NULL : kdt_extern_X,
                             NULL, NULL, NULL,
                             matrix_X_unordered_train,
                             matrix_X_ordered_train,
                             matrix_X_continuous_train,
                             TUNO,
                             TORD,
                             TCON,
                             Ycols,
                             Wcols,
                             NULL,
                             vector_scale_factor,
                             1,
                             matrix_bandwidth,
                             matrix_bandwidth,
                             lambda,
                             num_categories,
                             matrix_categorical_vals,
                             NULL,
                             out,
                             NULL,
                             NULL,
                             est_gate_ctx_ptr);

      for(i = 0; i < glp_nterms; i++){
        /* np_outer_weighted_sum lays out result as [W x Y], row-major by W term. */
        const int base = i*(glp_nterms + 2);
        XTKY[i][0] = out[base + 1]; /* Y column 1 is y */
        for(l = 0; l < glp_nterms; l++)
          KWM[i][l] = out[base + (l + 2)]; /* Y columns 2.. are basis terms */
      }

      while(mat_solve(KWM, XTKY, DELTA) == NULL){
        for(i = 0; i < glp_nterms; i++)
          KWM[i][i] += epsilon;
        nepsilon += epsilon;
      }

      XTKY[0][0] += nepsilon*XTKY[0][0]/NZD_POS(KWM[0][0]);
      if(nepsilon > 0.0){
        if(mat_solve(KWM, XTKY, DELTA) == NULL)
          error("mat_solve failed in glp path");
      }

      if(use_bernstein)
        np_glp_fill_basis_eval(num_reg_continuous,
                               glp_terms,
                               glp_nterms,
                               matrix_X_continuous_eval,
                               j,
                               basis_ctx,
                               eval_basis);
      else
        np_glp_fill_basis_eval_raw(num_reg_continuous,
                                   glp_terms,
                                   glp_nterms,
                                   matrix_X_continuous_eval,
                                   j,
                                   eval_basis);
      mean[j] = 0.0;
      for(i = 0; i < glp_nterms; i++)
        mean[j] += eval_basis[i]*DELTA[i][0];
      /* Row 0 corresponds to the constant basis term W0. */
      sk = copysign(DBL_MIN, out[2]) + out[2]; /* sum K * W0 * W0 */
      ey = out[1]/sk;  /* sum K * W0 * y / sk */
      ey2 = out[0]/sk; /* sum K * W0 * y^2 / sk */
      sigma2hat = ey2 - ey*ey;
      mean_stderr[j] = (sigma2hat <= 0.0) ? 0.0 : sqrt(sigma2hat * K_INT_KERNEL_P / (sk*hprod));
      sigma2hat = (sigma2hat <= 0.0) ? 0.0 : sigma2hat;

      for(l = 0; l < glp_nterms; l++){
        Ycols[l] = basis[l];
        Wcols[l] = basis[l];
      }

      kernel_weighted_sum_np_ctx(kernel_c,
                             kernel_u,
                             kernel_o,
                             BANDWIDTH_reg,
                             num_obs_train,
                             1,
                             num_reg_unordered,
                             num_reg_ordered,
                             num_reg_continuous,
                             0,
                             0,
                             2,
                             1,
                             0,
                             0,
                             0,
                             0,
                             0,
                             operator,
                             OP_NOOP,
                             0,
                             0,
                             NULL,
                             1,
                             glp_nterms,
                             glp_nterms,
                             (BANDWIDTH_reg == BW_ADAP_NN) ? NP_TREE_FALSE : int_TREE_X,
                             0,
                             (BANDWIDTH_reg == BW_ADAP_NN) ? NULL : kdt_extern_X,
                             NULL, NULL, NULL,
                             matrix_X_unordered_train,
                             matrix_X_ordered_train,
                             matrix_X_continuous_train,
                             TUNO,
                             TORD,
                             TCON,
                             Ycols,
                             Wcols,
                             NULL,
                             vector_scale_factor,
                             1,
                             matrix_bandwidth,
                             matrix_bandwidth,
                             lambda,
                             num_categories,
                             matrix_categorical_vals,
                             NULL,
                             out2,
                             NULL,
                             NULL,
                             est_gate_ctx_ptr);

      for(i = 0; i < glp_nterms; i++){
        const int base = i*glp_nterms;
        for(l = 0; l < glp_nterms; l++){
          KWM2[i][l] = out2[base + l];
          IDEN[i][l] = (i == l) ? 1.0 : 0.0;
        }
      }

      if(mat_solve(KWM, IDEN, KWM_INV) != NULL)
        have_vcov = 1;

      if(have_vcov){
        double q = 0.0;
        for(i = 0; i < glp_nterms; i++){
          double vv = 0.0;
          for(int ii = 0; ii < glp_nterms; ii++)
            vv += KWM_INV[i][ii]*eval_basis[ii];
          tmp_v[i] = vv;
        }
        for(i = 0; i < glp_nterms; i++){
          double ww = 0.0;
          for(int ii = 0; ii < glp_nterms; ii++)
            ww += KWM2[i][ii]*tmp_v[ii];
          tmp_w[i] = ww;
        }
        for(i = 0; i < glp_nterms; i++)
          q += tmp_v[i]*tmp_w[i];
        {
          const double mv = sigma2hat*q;
          if((mv > 0.0) && isfinite(mv))
            mean_stderr[j] = sqrt(mv);
        }
      }

      if(do_grad){
        for(l = 0; l < num_reg_continuous; l++){
          const int grad_order = (vector_glp_gradient_order_extern != NULL) ? MAX(1, vector_glp_gradient_order_extern[l]) : 1;
          if(use_bernstein)
            np_glp_fill_basis_eval_deriv(l,
                                         grad_order,
                                         num_reg_continuous,
                                         glp_terms,
                                         glp_nterms,
                                         matrix_X_continuous_eval,
                                         j,
                                         basis_ctx,
                                         eval_deriv);
          else
            np_glp_fill_basis_eval_deriv_raw(l,
                                             grad_order,
                                             num_reg_continuous,
                                             glp_terms,
                                             glp_nterms,
                                             matrix_X_continuous_eval,
                                             j,
                                             eval_deriv);
          gradient[l][j] = 0.0;
          for(i = 0; i < glp_nterms; i++)
            gradient[l][j] += eval_deriv[i]*DELTA[i][0];
          if(do_gerr){
            if(have_vcov){
              double q = 0.0;
              for(i = 0; i < glp_nterms; i++){
                double vv = 0.0;
                for(int ii = 0; ii < glp_nterms; ii++)
                  vv += KWM_INV[i][ii]*eval_deriv[ii];
                tmp_v[i] = vv;
              }
              for(i = 0; i < glp_nterms; i++){
                double ww = 0.0;
                for(int ii = 0; ii < glp_nterms; ii++)
                  ww += KWM2[i][ii]*tmp_v[ii];
                tmp_w[i] = ww;
              }
              for(i = 0; i < glp_nterms; i++)
                q += tmp_v[i]*tmp_w[i];
              {
                const double gv = sigma2hat*q;
                gradient_stderr[l][j] = (gv > 0.0 && isfinite(gv)) ? sqrt(gv) : 0.0;
              }
            } else {
              gradient_stderr[l][j] = 0.0;
            }
          }
        }
        for(l = num_reg_continuous; l < (num_reg_continuous + num_reg_unordered + num_reg_ordered); l++){
          gradient[l][j] = 0.0;
          if(do_gerr) gradient_stderr[l][j] = 0.0;
        }
      }

      if (fit_progress_active)
        np_progress_fit_loop_step(j + 1, fit_progress_total);
    }

    mat_free(KWM);
    mat_free(XTKY);
    mat_free(DELTA);
    mat_free(KWM2);
    mat_free(KWM_INV);
    mat_free(IDEN);
    free_mat(basis, glp_nterms);
    if((TCON != NULL) && (num_reg_continuous > 0)) free_mat(TCON, num_reg_continuous);
    if((TUNO != NULL) && (num_reg_unordered > 0)) free_mat(TUNO, num_reg_unordered);
    if((TORD != NULL) && (num_reg_ordered > 0)) free_mat(TORD, num_reg_ordered);
    free(Ycols);
    free(Wcols);
    free(y2);
    free(out);
    free(out2);
    free(tmp_v);
    free(tmp_w);
    if(use_bernstein){
      for(l = 0; l < num_reg_continuous; l++) np_glp_basis_ctx_free(&basis_ctx[l]);
      free(basis_ctx);
    }
    free(eval_basis);
    free(eval_deriv);
    free(glp_terms);

  } else { // Local Linear 

    // because we manipulate the training data scale factors can be wrong

    if((sf_flag = (int_LARGE_SF == 0)) && (BANDWIDTH_reg == BW_FIXED)){ 
      int_LARGE_SF = 1;
      vsf = (double *)malloc(num_reg_continuous*sizeof(double));
      for(int ii = 0; ii < num_reg_continuous; ii++)
        vsf[ii] = matrix_bandwidth[ii][0];
    } else {
      vsf = vector_scale_factor;
    }

    MATRIX XTKX = mat_creat( num_reg_continuous + 3, num_obs_train, UNDEFINED );
    MATRIX XTKY = mat_creat( num_reg_continuous + 1, 1, UNDEFINED );
    MATRIX DELTA = mat_creat( num_reg_continuous + 1, 1, UNDEFINED );

    MATRIX KWM = mat_creat( num_reg_continuous + 1, num_reg_continuous + 1, UNDEFINED );
    // Generate bandwidth vector given scale factors, nearest neighbors, or lambda 
    
    MATRIX TCON = mat_creat(num_reg_continuous, 1, UNDEFINED);
    MATRIX TUNO = mat_creat(num_reg_unordered, 1, UNDEFINED);
    MATRIX TORD = mat_creat(num_reg_ordered, 1, UNDEFINED);

    const int nrc3 = (num_reg_continuous+3);
    const int nrc1 = (num_reg_continuous+1);
    const int nrcc33 = nrc3*nrc3;

    double ** matrix_bandwidth_eval = NULL;

    double * PKWM[nrc1], * PXTKY[nrc1], * PXTKX[nrc3];

    double * PXC[MAX(1,num_reg_continuous)]; 
    double * PXU[MAX(1,num_reg_unordered)];
    double * PXO[MAX(1,num_reg_ordered)];

    PXC[0] = NULL;
    PXU[0] = NULL;
    PXO[0] = NULL;

    for(l = 0; l < num_reg_continuous; l++)
      PXC[l] = matrix_X_continuous_train[l];

    for(l = 0; l < num_reg_unordered; l++)
      PXU[l] = matrix_X_unordered_train[l];

    for(l = 0; l < num_reg_ordered; l++)
      PXO[l] = matrix_X_ordered_train[l];

    const size_t kwm_len = (size_t)nrcc33*(size_t)num_obs_eval_alloc;
    double * kwm = (double *)malloc(kwm_len*sizeof(double));

    for(size_t ii = 0; ii < kwm_len; ii++)
      kwm[ii] = 0.0;

    // with local linear, we already have the gradients of the continuous components
    // so we only need to worry about unordered + ordered comps

    double * permy = NULL;
    int ** moo = NULL;

    int p_nvar = do_grad ? (num_reg_unordered + num_reg_ordered) : 0;
    int do_ocg = do_grad && (p_nvar > 0);

    if(do_ocg){
      const size_t permy_len = kwm_len*(size_t)p_nvar;
      permy = (double *)malloc(permy_len*sizeof(double));
      if(permy == NULL)
        error("\n** Error: memory allocation failed.");
    }

    double * sgn = (double *)malloc((nrc3)*sizeof(double));

    sgn[0] = sgn[1] = sgn[2] = 1.0;
    
    for(int ii = 0; ii < (num_reg_continuous); ii++)
      sgn[ii+3] = -1.0;
    
    for(int ii = 0; ii < (nrc3); ii++)
      PXTKX[ii] = XTKX[ii];
    
    for(int ii = 0; ii < (nrc1); ii++){
      PKWM[ii] = KWM[ii];
      PXTKY[ii] = XTKY[ii];

      KWM[ii] = &kwm[(ii+2)*(nrc3)+2];
      XTKY[ii] = &kwm[ii+nrc3+2];
    }

    matrix_bandwidth_eval = alloc_tmatd(1,num_reg_continuous);

    const double epsilon = 1.0/num_obs_train;
    double nepsilon;

    // populate the xtkx matrix first 
    
    for(i = 0; i < num_obs_train; i++){
      const double vyi = vector_Y[i];
      XTKX[0][i] = vyi*vyi;
      XTKX[1][i] = vyi;
      XTKX[2][i] = 1.0;
    }

    if(do_ocg){
      moo = (int **)malloc(num_reg_ordered*sizeof(int *));
      for(l = 0; l < num_reg_ordered; l++){
        moo[l] = matrix_ordered_indices[l];
      }
    }
    for(j = 0; j < num_obs_eval; j++){ // main loop
      nepsilon = 0.0;

      for(l = 0; l < (nrc1); l++){
        KWM[l] = &kwm[j*nrcc33+(l+2)*(nrc3)+2];
        XTKY[l] = &kwm[j*nrcc33+l+nrc3+2];
      }

#ifdef MPI2
      
      if((j % iNum_Processors) == 0){
        if((j+my_rank) < (num_obs_eval)){
          for(l = 0; l < num_reg_continuous; l++){
          
            for(i = 0; i < num_obs_train; i++){
              XTKX[l+3][i] = matrix_X_continuous_train[l][i]-matrix_X_continuous_eval[l][j+my_rank];
            }
            TCON[l][0] = matrix_X_continuous_eval[l][j+my_rank]; // temporary storage

            if(BANDWIDTH_reg == BW_GEN_NN)
              matrix_bandwidth_eval[l][0] = matrix_bandwidth[l][j+my_rank]; // temporary storage

          }


          for(l = 0; l < num_reg_unordered; l++)
            TUNO[l][0] = matrix_X_unordered_eval[l][j+my_rank];

          for(l = 0; l < num_reg_ordered; l++)
            TORD[l][0] = matrix_X_ordered_eval[l][j+my_rank];

          kernel_weighted_sum_np_ctx(kernel_c,
                                 kernel_u,
                                 kernel_o,
                                 BANDWIDTH_reg,
                                 num_obs_train,
                                 1,
                                 num_reg_unordered,
                                 num_reg_ordered,
                                 num_reg_continuous,
                                 0, 
                                 0,
                                 1, // kernel_pow = 1
                                 1, // bandwidth_divide = FALSE when not adaptive
                                 0, 
                                 1, // symmetric
                                 0, // NO gather-scatter sum
                                 0, // do not drop train
                                 0, // do not drop train
                                 operator, // no convolution
                                 OP_NOOP, // no permutations
                                 0, // no score
                                 do_ocg, // ocg
                                 NULL,
                                 1, // explicity suppress parallel
                                 nrc3,
                                 nrc3,
                                 (BANDWIDTH_reg == BW_ADAP_NN) ? NP_TREE_FALSE : int_TREE_X,
                                 0,
                                 (BANDWIDTH_reg == BW_ADAP_NN) ? NULL : kdt_extern_X,
                                 NULL, NULL, NULL,
                                 PXU, // TRAIN
                                 PXO, 
                                 PXC,
                                 TUNO, // EVAL
                                 TORD,
                                 TCON,
                                 XTKX,
                                 XTKX,
                                 NULL,
                                 vsf,
                                 1,
                                 matrix_bandwidth,
                                 matrix_bandwidth_eval,
                                 lambda,
                                 num_categories,
                                 matrix_categorical_vals,
                                 moo,
                                 kwm+(j+my_rank)*nrcc33,  // weighted sum
                                 do_ocg ? (permy+(j+my_rank)*nrcc33*p_nvar) : NULL, // ocg
                                 NULL, // do not return kernel weights
                                 est_gate_ctx_ptr);

        }
        // synchro step
        if(!np_mpi_local_regression_active()){
          MPI_Allgather(MPI_IN_PLACE, nrcc33, MPI_DOUBLE, kwm+j*nrcc33, nrcc33, MPI_DOUBLE, comm[1]);          
          if(do_ocg){
            MPI_Allgather(MPI_IN_PLACE, nrcc33*p_nvar, MPI_DOUBLE, permy+j*nrcc33*p_nvar, nrcc33*p_nvar, MPI_DOUBLE, comm[1]);
          }
        }
      }
      
#else


      for(l = 0; l < num_reg_continuous; l++){
          
        for(i = 0; i < num_obs_train; i++){
          XTKX[l+3][i] = matrix_X_continuous_train[l][i]-matrix_X_continuous_eval[l][j];
        }
        TCON[l][0] = matrix_X_continuous_eval[l][j]; // temporary storage

        if(BANDWIDTH_reg == BW_GEN_NN)
          matrix_bandwidth_eval[l][0] = matrix_bandwidth[l][j]; // temporary storage

      }


      for(l = 0; l < num_reg_unordered; l++)
        TUNO[l][0] = matrix_X_unordered_eval[l][j];

      for(l = 0; l < num_reg_ordered; l++)
        TORD[l][0] = matrix_X_ordered_eval[l][j];

      kernel_weighted_sum_np_ctx(kernel_c,
                             kernel_u,
                             kernel_o,
                             BANDWIDTH_reg,
                             num_obs_train,
                             1,
                             num_reg_unordered,
                             num_reg_ordered,
                             num_reg_continuous,
                             0, // we leave one out via the weight matrix
                             0,
                             1, // kernel_pow = 1
                             1, // bandwidth_divide = FALSE when not adaptive
                             0, 
                             1, // symmetric
                             0, // gather-scatter sum
                             0, // do not drop train
                             0, // do not drop train
                             operator, // no convolution
                             OP_NOOP, // no permutations
                             0, // no score
                             do_ocg, // no ocg
                             NULL,
                             1, //  explicity suppress parallel
                             nrc3,
                             nrc3,
                             (BANDWIDTH_reg == BW_ADAP_NN) ? NP_TREE_FALSE : int_TREE_X,
                             0,
                             (BANDWIDTH_reg == BW_ADAP_NN) ? NULL : kdt_extern_X,
                             NULL, NULL, NULL,
                             PXU, // TRAIN
                             PXO, 
                             PXC,
                             TUNO, // EVAL
                             TORD,
                             TCON,
                             XTKX,
                             XTKX,
                             NULL,
                             vsf,
                             1,
                             matrix_bandwidth,
                             matrix_bandwidth_eval,
                             lambda,
                             num_categories,
                             matrix_categorical_vals,
                             moo,
                             kwm+j*nrcc33,  // weighted sum
                             do_ocg ? (permy+j*nrcc33*p_nvar) : NULL, // no permutations
                             NULL, // do not return kernel weights
                             est_gate_ctx_ptr);

#endif
      
      if(do_ocg){
        for(l = 0; l < num_reg_ordered; l++){
          moo[l]++;
        }
      }

      while(mat_solve(KWM, XTKY, DELTA) == NULL){ // singular = ridge about
        for(int ii = 0; ii < (nrc1); ii++)
          KWM[ii][ii] += epsilon;
        nepsilon += epsilon;
      }

      XTKY[0][0] += nepsilon*XTKY[0][0]/NZD_POS(KWM[0][0]);
      if(nepsilon > 0.0){
        if(mat_solve(KWM, XTKY, DELTA) == NULL)
          error("mat_solve failed after ridge adjustment");
      }
      mean[j] = DELTA[0][0];

      const double sk = copysign(DBL_MIN, (kwm+j*nrcc33)[2*nrc3+2]) + (kwm+j*nrcc33)[2*nrc3+2];
      const double ey = (kwm+j*nrcc33)[nrc3+2]/sk;
      const double ey2 = (kwm+j*nrcc33)[nrc3+1]/sk;

      {
        const double v = (ey2 - ey*ey)*K_INT_KERNEL_P / (sk*hprod);
        mean_stderr[j] = (v <= 0.0) ? 0.0 : sqrt(v);
      }

      if(do_grad){
        for(int ii = 0; ii < num_reg_continuous; ii++){
          gradient[ii][j] = DELTA[ii+1][0];
          if(do_gerr)
            gradient_stderr[ii][j] = gfac*mean_stderr[j]/((BANDWIDTH_reg == BW_ADAP_NN) ? 1.0 : ((BANDWIDTH_reg == BW_GEN_NN) ? matrix_bandwidth[ii][j]:matrix_bandwidth[ii][0]));
        }
        
        // we need to do new matrix inversions here for the unordered + ordered data
        // we can safely taint KWM , DELTA, and XTKY here
        for(l = num_reg_continuous; l < (num_reg_continuous + num_reg_unordered); l++){
          const int dl = l - num_reg_continuous;
          const int ojp = j*nrcc33*p_nvar + dl*nrcc33;

          for(int ii = 0; ii < nrc1; ii++){
            KWM[ii] = &permy[ojp +(ii+2)*(nrc3)+2];
            XTKY[ii] = &permy[ojp + ii + nrc3 + 2];
          }

          nepsilon = 0.0;
          while(mat_solve(KWM, XTKY, DELTA) == NULL){ // singular = ridge about
            for(int ii = 0; ii < (nrc1); ii++)
              KWM[ii][ii] += epsilon;
            nepsilon += epsilon;
          }

          XTKY[0][0] += nepsilon*XTKY[0][0]/NZD_POS(KWM[0][0]);
          if(nepsilon > 0.0){
            if(mat_solve(KWM, XTKY, DELTA) == NULL)
              error("mat_solve failed after ridge adjustment");
          }

          gradient[l][j] = mean[j] - DELTA[0][0];
          
          if(do_gerr){
            const double skg = copysign(DBL_MIN, (permy+ojp)[2*nrc3+2]) + (permy+ojp)[2*nrc3+2];
            const double eyg = (permy+ojp)[nrc3+2]/skg;
            const double ey2g = (permy+ojp)[nrc3+1]/skg;

            const double se2 = (ey2g - eyg*eyg)*K_INT_KERNEL_P / (skg*hprod);

            gradient_stderr[l][j] = sqrt(mean_stderr[j]*mean_stderr[j] + se2);
          } else {
            gradient_stderr[l][j] = 0.0;
          }
        }

        // we need to do new matrix inversions here for the unordered + ordered data
        // we can safely taint KWM , DELTA, and XTKY here
        for(l = num_reg_continuous + num_reg_unordered; l < (num_reg_continuous + num_reg_unordered + num_reg_ordered); l++){

          const int dl = l - num_reg_continuous;
          const int ojp = j*nrcc33*p_nvar + dl*nrcc33;

          for(int ii = 0; ii < nrc1; ii++){
            KWM[ii] = &permy[ojp +(ii+2)*(nrc3)+2];
            XTKY[ii] = &permy[ojp + ii + nrc3 + 2];
          }

          nepsilon = 0.0;
          while(mat_solve(KWM, XTKY, DELTA) == NULL){ // singular = ridge about
            for(int ii = 0; ii < (nrc1); ii++)
              KWM[ii][ii] += epsilon;
            nepsilon += epsilon;
          }

          XTKY[0][0] += nepsilon*XTKY[0][0]/NZD_POS(KWM[0][0]);
          if(nepsilon > 0.0){
            if(mat_solve(KWM, XTKY, DELTA) == NULL)
              error("mat_solve failed after ridge adjustment");
          }

          gradient[l][j] = (mean[j] - DELTA[0][0])*((matrix_ordered_indices[l - num_reg_continuous - num_reg_unordered][j] != 0) ? 1.0 : -1.0);
          
          if(do_gerr){

            const double skg = copysign(DBL_MIN, (permy+ojp)[2*nrc3+2]) + (permy+ojp)[2*nrc3+2];
            const double eyg = (permy+ojp)[nrc3+2]/skg;
            const double ey2g = (permy+ojp)[nrc3+1]/skg;

            const double se2 = (ey2g - eyg*eyg)*K_INT_KERNEL_P / (skg*hprod);

            gradient_stderr[l][j] = sqrt(mean_stderr[j]*mean_stderr[j] + se2);
          } else {
            gradient_stderr[l][j] = 0.0;
          }
        }

      }

      if (fit_progress_active)
        np_progress_fit_loop_step(j + 1, fit_progress_total);
    }
    
    for(int ii = 0; ii < (nrc1); ii++){
      KWM[ii] = PKWM[ii];
      XTKY[ii] = PXTKY[ii];
    }

    for(int ii = 0; ii < (nrc3); ii++)
      XTKX[ii] = PXTKX[ii];

    if(sf_flag){
      int_LARGE_SF = 0;
      free(vsf);
    }

    
    mat_free(XTKX);
    mat_free(XTKY);
    mat_free(DELTA);
    mat_free(KWM);

    mat_free(TCON);
    mat_free(TUNO);
    mat_free(TORD);

    free(kwm);
    free(sgn);

    if(do_ocg){
      free(permy);
      free(moo);
    }

    free_tmat(matrix_bandwidth_eval);
  }

finish_regression_estimation:
  // clean up hash stuff
  if(do_grad && (num_reg_ordered > 0)){
    for(l = 0; l < num_reg_ordered; l++)
      thdestroy_r(otabs+l);
    free(otabs);
    free(matrix_ordered_indices[0]);
    free(matrix_ordered_indices);
  }

  if(gate_override_active) np_gate_ctx_clear(&gate_ctx_local);
  if((ov_cont_ok != NULL) && (!ov_cont_from_cache)) free(ov_cont_ok);
  if((ov_cont_hmin != NULL) && (!ov_cont_from_cache)) free(ov_cont_hmin);
  if((ov_cont_k0 != NULL) && (!ov_cont_from_cache)) free(ov_cont_k0);
  if(ov_disc_uno_ok != NULL) free(ov_disc_uno_ok);
  if(ov_disc_uno_const != NULL) free(ov_disc_uno_const);
  if(ov_disc_ord_ok != NULL) free(ov_disc_ord_ok);
  if(ov_disc_ord_const != NULL) free(ov_disc_ord_const);

  free(operator);
  free(kernel_c);
  free(kernel_u);
  free(kernel_o);
  free(lambda);
  free_tmat(matrix_bandwidth);
  free_mat(matrix_bandwidth_deriv,num_reg_continuous);

	if(vector_Y_eval != NULL)
	{
		*R_squared = fGoodness_of_Fit(num_obs_eval, vector_Y_eval, mean);
		*MSE = fMSE(num_obs_eval, vector_Y_eval, mean);
		*MAE = fMAE(num_obs_eval, vector_Y_eval, mean);
		*MAPE = fMAPE(num_obs_eval, vector_Y_eval, mean);
		*CORR = fCORR(num_obs_eval, vector_Y_eval, mean);
		*SIGN = fSIGN(num_obs_eval, vector_Y_eval, mean);
	}
	else
	{
		*R_squared = 0.0;
		*MSE = 0.0;
		*MAE = 0.0;
		*MAPE = 0.0;
		*CORR = 0.0;
		*SIGN = 0.0;
	}

  return(0);
}

static int np_conditional_indicator_row_core(const int train_idx,
                                             const int eval_idx,
                                             const int cdfontrain,
                                             double **matrix_Y_ordered_train,
                                             double **matrix_Y_continuous_train,
                                             double **matrix_Y_ordered_eval,
                                             double **matrix_Y_continuous_eval,
                                             const int num_var_ordered,
                                             const int num_var_continuous){
  int l;
  int indy = 1;

  if(cdfontrain && (train_idx == eval_idx))
    return 0;

  for(l = 0; (l < num_var_ordered) && (indy != 0); l++)
    indy *= (matrix_Y_ordered_train[l][train_idx] <= matrix_Y_ordered_eval[l][eval_idx]);

  for(l = 0; (l < num_var_continuous) && (indy != 0); l++)
    indy *= (matrix_Y_continuous_train[l][train_idx] <= matrix_Y_continuous_eval[l][eval_idx]);

  return indy;
}

static int np_shadow_conditional_indicator_row(const int train_idx,
                                               const int eval_idx,
                                               const int cdfontrain,
                                               double **matrix_Y_ordered_train,
                                               double **matrix_Y_continuous_train,
                                               double **matrix_Y_ordered_eval,
                                               double **matrix_Y_continuous_eval,
                                               const int num_var_ordered,
                                               const int num_var_continuous){
  return np_conditional_indicator_row_core(train_idx,
                                           eval_idx,
                                           cdfontrain,
                                           matrix_Y_ordered_train,
                                           matrix_Y_continuous_train,
                                           matrix_Y_ordered_eval,
                                           matrix_Y_continuous_eval,
                                           num_var_ordered,
                                           num_var_continuous);
}

static int np_shadow_conditional_kernel_row_core(const int *kernel_c,
                                                 const int *kernel_u,
                                                 const int *kernel_o,
                                                 const int *operator,
                                                 const int BANDWIDTH_den,
                                                 const int num_train,
                                                 const int num_unordered,
                                                 const int num_ordered,
                                                 const int num_continuous,
                                                 double **matrix_train_unordered,
                                                 double **matrix_train_ordered,
                                                 double **matrix_train_continuous,
                                                 double **matrix_eval_unordered_one,
                                                 double **matrix_eval_ordered_one,
                                                 double **matrix_eval_continuous_one,
                                                 double *vector_scale_factor,
                                                 int bandwidth_provided,
                                                 double **matrix_bandwidth_train,
                                                 double **matrix_bandwidth_eval_one,
                                                 double *lambda,
                                                 int *num_categories,
                                                 double **matrix_categorical_vals,
                                                 const int int_tree,
                                                 KDT *kdt,
                                                 const int bandwidth_divide_weights,
                                                 double *kw_out,
                                                 double *mean_out){
  return kernel_weighted_sum_np_ctx((int *)kernel_c,
                                    (int *)kernel_u,
                                    (int *)kernel_o,
                                    BANDWIDTH_den,
                                    num_train,
                                    1,
                                    num_unordered,
                                    num_ordered,
                                    num_continuous,
                                    0,
                                    0,
                                    1,
                                    1,
                                    bandwidth_divide_weights,
                                    0,
                                    0,
                                    0,
                                    0,
                                    operator,
                                    OP_NOOP,
                                    0,
                                    0,
                                    NULL,
                                    1,
                                    0,
                                    0,
                                    int_tree,
                                    0,
                                    kdt,
                                    NULL,
                                    NULL,
                                    NULL,
                                    matrix_train_unordered,
                                    matrix_train_ordered,
                                    matrix_train_continuous,
                                    matrix_eval_unordered_one,
                                    matrix_eval_ordered_one,
                                    matrix_eval_continuous_one,
                                    NULL,
                                    NULL,
                                    NULL,
                                    vector_scale_factor,
                                    bandwidth_provided,
                                    matrix_bandwidth_train,
                                    matrix_bandwidth_eval_one,
                                    lambda,
                                    num_categories,
                                    matrix_categorical_vals,
                                    NULL,
                                    mean_out,
                                    NULL,
                                    kw_out,
                                    NULL);
}

static int np_shadow_conditional_kernel_row(const int *kernel_c,
                                            const int *kernel_u,
                                            const int *kernel_o,
                                            const int *operator,
                                            const int BANDWIDTH_den,
                                            const int num_train,
                                            const int num_unordered,
                                            const int num_ordered,
                                            const int num_continuous,
                                            double **matrix_train_unordered,
                                            double **matrix_train_ordered,
                                            double **matrix_train_continuous,
                                            double **matrix_eval_unordered_one,
                                            double **matrix_eval_ordered_one,
                                            double **matrix_eval_continuous_one,
                                            double *vector_scale_factor,
                                            int bandwidth_provided,
                                            double **matrix_bandwidth_train,
                                            double **matrix_bandwidth_eval_one,
                                            double *lambda,
                                            int *num_categories,
                                            double **matrix_categorical_vals,
                                            const int int_tree,
                                            KDT *kdt,
                                            double *kw_out,
                                            double *mean_out){
  return np_shadow_conditional_kernel_row_core(kernel_c,
                                               kernel_u,
                                               kernel_o,
                                               operator,
                                               BANDWIDTH_den,
                                               num_train,
                                               num_unordered,
                                               num_ordered,
                                               num_continuous,
                                               matrix_train_unordered,
                                               matrix_train_ordered,
                                               matrix_train_continuous,
                                               matrix_eval_unordered_one,
                                               matrix_eval_ordered_one,
                                               matrix_eval_continuous_one,
                                               vector_scale_factor,
                                               bandwidth_provided,
                                               matrix_bandwidth_train,
                                               matrix_bandwidth_eval_one,
                                               lambda,
                                               num_categories,
                                               matrix_categorical_vals,
                                               int_tree,
                                               kdt,
                                               1,
                                               kw_out,
                                               mean_out);
}

static int np_shadow_conditional_kernel_row_raw(const int *kernel_c,
                                                const int *kernel_u,
                                                const int *kernel_o,
                                                const int *operator,
                                                const int BANDWIDTH_den,
                                                const int num_train,
                                                const int num_unordered,
                                                const int num_ordered,
                                                const int num_continuous,
                                                double **matrix_train_unordered,
                                                double **matrix_train_ordered,
                                                double **matrix_train_continuous,
                                                double **matrix_eval_unordered_one,
                                                double **matrix_eval_ordered_one,
                                                double **matrix_eval_continuous_one,
                                                double *vector_scale_factor,
                                                int bandwidth_provided,
                                                double **matrix_bandwidth_train,
                                                double **matrix_bandwidth_eval_one,
                                                double *lambda,
                                                int *num_categories,
                                                double **matrix_categorical_vals,
                                                const int int_tree,
                                                KDT *kdt,
                                                double *kw_out,
                                                double *mean_out){
  return np_shadow_conditional_kernel_row_core(kernel_c,
                                               kernel_u,
                                               kernel_o,
                                               operator,
                                               BANDWIDTH_den,
                                               num_train,
                                               num_unordered,
                                               num_ordered,
                                               num_continuous,
                                               matrix_train_unordered,
                                               matrix_train_ordered,
                                               matrix_train_continuous,
                                               matrix_eval_unordered_one,
                                               matrix_eval_ordered_one,
                                               matrix_eval_continuous_one,
                                               vector_scale_factor,
                                               bandwidth_provided,
                                               matrix_bandwidth_train,
                                               matrix_bandwidth_eval_one,
                                               lambda,
                                               num_categories,
                                               matrix_categorical_vals,
                                               int_tree,
                                               kdt,
                                               0,
                                               kw_out,
                                               mean_out);
}

static int np_glp_qr_drop_row_bkcde(double **basis,
                                    int num_train,
                                    int k,
                                    const double *kw,
                                    int eval_pos,
                                    double *row_out){
  const double tol = 1.0e-7;
  double *xqr = NULL, *qraux = NULL, *work = NULL, *y = NULL, *qy = NULL;
  int *pivot = NULL;
  int ldx = num_train, n = num_train, p = k, rank = 0, ny = 1;
  int i, j;
  int status = 1;

  if((basis == NULL) || (kw == NULL) || (row_out == NULL) || (num_train <= 0) || (k <= 0))
    return 1;
  if((eval_pos < 0) || (eval_pos >= num_train))
    return 1;

  xqr = (double *)malloc((size_t)num_train*(size_t)k*sizeof(double));
  qraux = (double *)malloc((size_t)k*sizeof(double));
  work = (double *)malloc((size_t)(2*k)*sizeof(double));
  y = (double *)calloc((size_t)num_train, sizeof(double));
  qy = (double *)malloc((size_t)num_train*sizeof(double));
  pivot = (int *)malloc((size_t)k*sizeof(int));
  if((xqr == NULL) || (qraux == NULL) || (work == NULL) ||
     (y == NULL) || (qy == NULL) || (pivot == NULL))
    goto cleanup_glp_qr_drop;

  for(j = 0; j < k; j++){
    pivot[j] = j + 1;
    for(i = 0; i < num_train; i++){
      const double w = kw[i];
      xqr[i + j*num_train] = ((w > 0.0) ? sqrt(w) : 0.0) * basis[j][i];
    }
  }

  F77_NAME(dqrdc2)(xqr, &ldx, &n, &p, (double *)&tol, &rank, qraux, pivot, work);
  if((rank < 0) || (rank > k))
    goto cleanup_glp_qr_drop;
  if(rank < k)
    goto cleanup_glp_qr_drop;

  for(i = 0; i < rank; i++){
    double s = basis[i][eval_pos];
    for(j = 0; j < i; j++)
      s -= xqr[j + i*num_train]*y[j];
    if(fabs(xqr[i + i*num_train]) <= DBL_MIN)
      goto cleanup_glp_qr_drop;
    y[i] = s/xqr[i + i*num_train];
  }

  F77_NAME(dqrqy)(xqr, &n, &p, qraux, y, &ny, qy);

  for(i = 0; i < num_train; i++){
    const double w = kw[i];
    row_out[i] = ((w > 0.0) ? sqrt(w) : 0.0) * qy[i];
  }

  status = 0;

cleanup_glp_qr_drop:
  if(xqr != NULL) free(xqr);
  if(qraux != NULL) free(qraux);
  if(work != NULL) free(work);
  if(y != NULL) free(y);
  if(qy != NULL) free(qy);
  if(pivot != NULL) free(pivot);
  return status;
}

static int np_mat_bad_rcond_sym(MATRIX A, double min_rcond){
  const int n = MatRow(A);
  char jobz = 'N';
  char uplo = 'U';
  int i, j;
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  double max_eval = 0.0, min_eval = DBL_MAX;
  double *Ac = NULL, *eval = NULL, *work = NULL;
  int bad = 1;

  if((A == NULL) || (MatCol(A) != n) || (n <= 0))
    return 1;

  Ac = (double *)malloc((size_t)n*(size_t)n*sizeof(double));
  eval = (double *)malloc((size_t)n*sizeof(double));
  if((Ac == NULL) || (eval == NULL))
    goto cleanup_bad_rcond;

  for(j = 0; j < n; j++)
    for(i = 0; i < n; i++)
      Ac[i + j*n] = A[i][j];

  F77_CALL(dsyev)(&jobz, &uplo, &n, Ac, &n, eval, &work_query, &lwork, &info FCONE FCONE);
  if(info != 0)
    goto cleanup_bad_rcond;

  lwork = MAX(1, (int)work_query);
  work = (double *)malloc((size_t)lwork*sizeof(double));
  if(work == NULL)
    goto cleanup_bad_rcond;

  F77_CALL(dsyev)(&jobz, &uplo, &n, Ac, &n, eval, work, &lwork, &info FCONE FCONE);
  if(info != 0)
    goto cleanup_bad_rcond;

  for(i = 0; i < n; i++){
    if(!R_FINITE(eval[i]))
      goto cleanup_bad_rcond;
    max_eval = MAX(max_eval, fabs(eval[i]));
    min_eval = MIN(min_eval, fabs(eval[i]));
  }

  if(!(max_eval > 0.0))
    goto cleanup_bad_rcond;

  bad = ((min_eval / max_eval) < min_rcond);

cleanup_bad_rcond:
  if(Ac != NULL) free(Ac);
  if(eval != NULL) free(eval);
  if(work != NULL) free(work);
  return bad;
}

typedef struct {
  int int_cker_bound;
  double *ckerlb;
  double *ckerub;
} NPConditionalBoundState;

static void np_conditional_push_bounds(int new_bound_flag,
                                       double *new_lb,
                                       double *new_ub,
                                       NPConditionalBoundState *saved);
static void np_conditional_pop_bounds(const NPConditionalBoundState *saved);

static int np_shadow_conditional_build_x_weights_core(double *vector_scale_factor,
                                                      int drop_eval_self,
                                                      double *weights_out){
  const int num_train = num_obs_train_extern;
  const int num_reg_tot = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
  const int ll_mode = (int_ll_extern == LL_LP) ? LL_LP : LL_LC;
  const int bw_rows = (BANDWIDTH_den_extern == BW_FIXED) ? 1 : num_train;
  int *kernel_cx = NULL, *kernel_ux = NULL, *kernel_ox = NULL, *x_operator = NULL;
  double *vsfx = NULL, *lambdax = NULL, *kw = NULL, *mean_row = NULL;
  double **matrix_bandwidth_x = NULL, **matrix_bandwidth_eval_one = NULL;
  double **eval_xuno_one = NULL, **eval_xord_one = NULL, **eval_xcon_one = NULL;
  MATRIX KWM = NULL, RHS = NULL, SOL = NULL;
  NPConditionalBoundState bounds_state;
  int i, j, l;
  int status = 1;

  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;

  if(num_train <= 0)
    return 1;

  memset(weights_out, 0, (size_t)num_train*(size_t)num_train*sizeof(double));

  if(num_reg_tot <= 0){
    const double denom = drop_eval_self ? ((double)(num_train - 1)) : ((double)num_train);
    const double w = (denom > 0.0) ? 1.0/denom : 0.0;
    for(i = 0; i < num_train; i++)
      for(j = 0; j < num_train; j++)
        if((!drop_eval_self) || (i != j))
          weights_out[(size_t)i*(size_t)num_train + (size_t)j] = w;
    return 0;
  }

  vsfx = alloc_vecd(MAX(1, num_reg_tot));
  lambdax = alloc_vecd(MAX(1, num_reg_unordered_extern + num_reg_ordered_extern));
  kw = alloc_vecd(MAX(1, num_train));
  mean_row = alloc_vecd(MAX(1, num_train));
  matrix_bandwidth_x = alloc_tmatd(bw_rows, num_reg_continuous_extern);
  matrix_bandwidth_eval_one = alloc_tmatd(1, num_reg_continuous_extern);
  if(num_reg_unordered_extern > 0) eval_xuno_one = alloc_matd(1, num_reg_unordered_extern);
  if(num_reg_ordered_extern > 0) eval_xord_one = alloc_matd(1, num_reg_ordered_extern);
  if(num_reg_continuous_extern > 0) eval_xcon_one = alloc_matd(1, num_reg_continuous_extern);

  kernel_cx = (int *)calloc((size_t)MAX(1, num_reg_continuous_extern), sizeof(int));
  kernel_ux = (int *)calloc((size_t)MAX(1, num_reg_unordered_extern), sizeof(int));
  kernel_ox = (int *)calloc((size_t)MAX(1, num_reg_ordered_extern), sizeof(int));
  x_operator = (int *)calloc((size_t)MAX(1, num_reg_tot), sizeof(int));

  if((vsfx == NULL) || (lambdax == NULL) || (kw == NULL) || (mean_row == NULL) ||
     ((num_reg_continuous_extern > 0) && (matrix_bandwidth_x == NULL)) ||
     ((num_reg_continuous_extern > 0) && (matrix_bandwidth_eval_one == NULL)) ||
     ((num_reg_unordered_extern > 0) && (eval_xuno_one == NULL)) ||
     ((num_reg_ordered_extern > 0) && (eval_xord_one == NULL)) ||
     ((num_reg_continuous_extern > 0) && (eval_xcon_one == NULL)) ||
     (kernel_cx == NULL) || (kernel_ux == NULL) || (kernel_ox == NULL) || (x_operator == NULL))
    goto cleanup_xweights;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        vsfx,
                        NULL,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  for(i = 0; i < num_reg_continuous_extern; i++) kernel_cx[i] = KERNEL_reg_extern;
  for(i = 0; i < num_reg_unordered_extern; i++) kernel_ux[i] = KERNEL_reg_unordered_extern;
  for(i = 0; i < num_reg_ordered_extern; i++) kernel_ox[i] = KERNEL_reg_ordered_extern;
  for(i = 0; i < num_reg_tot; i++) x_operator[i] = OP_NORMAL;

  if(kernel_bandwidth_mean(KERNEL_reg_extern,
                           BANDWIDTH_den_extern,
                           num_train,
                           num_train,
                           0,
                           0,
                           0,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           vsfx,
                           NULL,
                           NULL,
                           matrix_X_continuous_train_extern,
                           matrix_X_continuous_train_extern,
                           NULL,
                           matrix_bandwidth_x,
                           lambdax) == 1)
    goto cleanup_xweights;

  if(ll_mode == LL_LP){
    if((vector_glp_degree_extern == NULL) || (num_reg_continuous_extern <= 0))
      goto cleanup_xweights;
    if(!np_glp_cv_prepare_extern(LL_LP,
                                 num_train,
                                 num_reg_continuous_extern,
                                 matrix_X_continuous_train_extern))
      goto cleanup_xweights;
  }

  for(i = 0; i < num_train; i++){
    const int orig_i = (int_TREE_X == NP_TREE_TRUE) ? ipt_extern_X[i] : i;
    double row_sum = 0.0;

    for(l = 0; l < num_reg_unordered_extern; l++)
      eval_xuno_one[l][0] = matrix_X_unordered_train_extern[l][i];
    for(l = 0; l < num_reg_ordered_extern; l++)
      eval_xord_one[l][0] = matrix_X_ordered_train_extern[l][i];
    for(l = 0; l < num_reg_continuous_extern; l++){
      eval_xcon_one[l][0] = matrix_X_continuous_train_extern[l][i];
      matrix_bandwidth_eval_one[l][0] =
        (BANDWIDTH_den_extern == BW_GEN_NN) ? matrix_bandwidth_x[l][i] : matrix_bandwidth_x[l][0];
    }

    np_conditional_push_bounds(int_cxker_bound_extern,
                               vector_cxkerlb_extern,
                               vector_cxkerub_extern,
                               &bounds_state);
    if(np_shadow_conditional_kernel_row_raw(kernel_cx,
                                            kernel_ux,
                                            kernel_ox,
                                            x_operator,
                                            BANDWIDTH_den_extern,
                                            num_train,
                                            num_reg_unordered_extern,
                                            num_reg_ordered_extern,
                                            num_reg_continuous_extern,
                                            matrix_X_unordered_train_extern,
                                            matrix_X_ordered_train_extern,
                                            matrix_X_continuous_train_extern,
                                            eval_xuno_one,
                                            eval_xord_one,
                                            eval_xcon_one,
                                            vsfx,
                                            1,
                                            matrix_bandwidth_x,
                                            matrix_bandwidth_eval_one,
                                            lambdax,
                                            num_categories_extern_X,
                                            matrix_categorical_vals_extern_X,
                                            int_TREE_X,
                                            kdt_extern_X,
                                            kw,
                                            mean_row) != 0){
      np_conditional_pop_bounds(&bounds_state);
      goto cleanup_xweights;
    }
    np_conditional_pop_bounds(&bounds_state);

    if(drop_eval_self)
      kw[i] = 0.0;

    if(ll_mode == LL_LC){
      for(j = 0; j < num_train; j++)
        row_sum += kw[j];
      if(!(row_sum > DBL_MIN))
        goto cleanup_xweights;
      for(j = 0; j < num_train; j++){
        const int orig_j = (int_TREE_X == NP_TREE_TRUE) ? ipt_extern_X[j] : j;
        weights_out[(size_t)orig_i*(size_t)num_train + (size_t)orig_j] = kw[j]/row_sum;
      }
    } else {
      const int k = np_glp_cv_cache.nterms;

      if((k <= 0) || (np_glp_cv_cache.basis == NULL))
        goto cleanup_xweights;

      KWM = mat_creat(k, k, UNDEFINED);
      RHS = mat_creat(k, 1, UNDEFINED);
      SOL = mat_creat(k, 1, UNDEFINED);
      if((KWM == NULL) || (RHS == NULL) || (SOL == NULL))
        goto cleanup_xweights;

      if(drop_eval_self){
        if(np_glp_qr_drop_row_bkcde(np_glp_cv_cache.basis,
                                    num_train,
                                    k,
                                    kw,
                                    i,
                                    mean_row) != 0)
          goto cleanup_xweights;

        for(j = 0; j < num_train; j++){
          const int orig_j = (int_TREE_X == NP_TREE_TRUE) ? ipt_extern_X[j] : j;
          weights_out[(size_t)orig_i*(size_t)num_train + (size_t)orig_j] = mean_row[j];
        }
      } else {
        for(l = 0; l < k; l++){
          RHS[l][0] = np_glp_cv_cache.basis[l][i];
          for(j = 0; j < k; j++)
            KWM[l][j] = 0.0;
        }

        for(j = 0; j < num_train; j++){
          const double wj = kw[j];
          if(wj == 0.0)
            continue;
          for(int a = 0; a < k; a++){
            const double za = np_glp_cv_cache.basis[a][j];
            for(int b = a; b < k; b++){
              const double zb = np_glp_cv_cache.basis[b][j];
              KWM[a][b] += wj*za*zb;
              if(b != a) KWM[b][a] += wj*za*zb;
            }
          }
        }

        if(np_mat_bad_rcond_sym(KWM, 1.0e-10))
          goto cleanup_xweights;
        if(mat_solve(KWM, RHS, SOL) == NULL)
          goto cleanup_xweights;

        for(j = 0; j < num_train; j++){
          double zju = 0.0;
          const int orig_j = (int_TREE_X == NP_TREE_TRUE) ? ipt_extern_X[j] : j;
          for(l = 0; l < k; l++)
            zju += np_glp_cv_cache.basis[l][j]*SOL[l][0];
          weights_out[(size_t)orig_i*(size_t)num_train + (size_t)orig_j] = kw[j]*zju;
        }
      }

      mat_free(KWM);
      mat_free(RHS);
      mat_free(SOL);
      KWM = RHS = SOL = NULL;
    }
  }

  status = 0;

cleanup_xweights:
  if(KWM != NULL) mat_free(KWM);
  if(RHS != NULL) mat_free(RHS);
  if(SOL != NULL) mat_free(SOL);
  if(vsfx != NULL) free(vsfx);
  if(lambdax != NULL) free(lambdax);
  if(kw != NULL) free(kw);
  if(mean_row != NULL) free(mean_row);
  if(matrix_bandwidth_x != NULL) free_tmat(matrix_bandwidth_x);
  if(matrix_bandwidth_eval_one != NULL) free_tmat(matrix_bandwidth_eval_one);
  if(eval_xuno_one != NULL) free_mat(eval_xuno_one, num_reg_unordered_extern);
  if(eval_xord_one != NULL) free_mat(eval_xord_one, num_reg_ordered_extern);
  if(eval_xcon_one != NULL) free_mat(eval_xcon_one, num_reg_continuous_extern);
  if(kernel_cx != NULL) free(kernel_cx);
  if(kernel_ux != NULL) free(kernel_ux);
  if(kernel_ox != NULL) free(kernel_ox);
  if(x_operator != NULL) free(x_operator);
  return status;
}

static int np_shadow_conditional_build_x_weights(double *vector_scale_factor,
                                                 double *weights_out){
  return np_shadow_conditional_build_x_weights_core(vector_scale_factor, 1, weights_out);
}

static int np_shadow_conditional_build_x_weights_full(double *vector_scale_factor,
                                                      double *weights_out){
  return np_shadow_conditional_build_x_weights_core(vector_scale_factor, 0, weights_out);
}

typedef struct {
  int ready;
  int ll_mode;
  int num_train;
  int num_reg_tot;
  int *kernel_cx;
  int *kernel_ux;
  int *kernel_ox;
  int *x_operator;
  double *vsfx;
  double *lambdax;
  double *kw;
  double *mean_row;
  double **matrix_bandwidth_x;
  double **matrix_bandwidth_eval_one;
  double **eval_xuno_one;
  double **eval_xord_one;
  double **eval_xcon_one;
} NPConditionalXRowCtx;

typedef struct {
  int ready;
  int num_train;
  int num_var_tot;
  int *kernel_cy;
  int *kernel_uy;
  int *kernel_oy;
  int *operator_y;
  double *vsfy;
  double *lambday;
  double *kw;
  double **matrix_bandwidth_y;
  double **matrix_bandwidth_eval_one;
  double **eval_yuno_one;
  double **eval_yord_one;
  double **eval_ycon_one;
} NPConditionalYRowCtx;

static void np_conditional_xrow_ctx_clear(NPConditionalXRowCtx *ctx){
  if(ctx == NULL)
    return;
  if(ctx->vsfx != NULL) free(ctx->vsfx);
  if(ctx->lambdax != NULL) free(ctx->lambdax);
  if(ctx->kw != NULL) free(ctx->kw);
  if(ctx->mean_row != NULL) free(ctx->mean_row);
  if(ctx->matrix_bandwidth_x != NULL) free_tmat(ctx->matrix_bandwidth_x);
  if(ctx->matrix_bandwidth_eval_one != NULL) free_tmat(ctx->matrix_bandwidth_eval_one);
  if(ctx->eval_xuno_one != NULL) free_mat(ctx->eval_xuno_one, num_reg_unordered_extern);
  if(ctx->eval_xord_one != NULL) free_mat(ctx->eval_xord_one, num_reg_ordered_extern);
  if(ctx->eval_xcon_one != NULL) free_mat(ctx->eval_xcon_one, num_reg_continuous_extern);
  if(ctx->kernel_cx != NULL) free(ctx->kernel_cx);
  if(ctx->kernel_ux != NULL) free(ctx->kernel_ux);
  if(ctx->kernel_ox != NULL) free(ctx->kernel_ox);
  if(ctx->x_operator != NULL) free(ctx->x_operator);
  memset(ctx, 0, sizeof(*ctx));
}

static int np_conditional_xrow_ctx_prepare(double *vector_scale_factor,
                                           NPConditionalXRowCtx *ctx){
  const int num_train = num_obs_train_extern;
  const int num_reg_tot = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
  const int ll_mode = (int_ll_extern == LL_LP) ? LL_LP : LL_LC;
  const int bw_rows = (BANDWIDTH_den_extern == BW_FIXED) ? 1 : num_train;
  int i;

  if((ctx == NULL) || (vector_scale_factor == NULL))
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;
  if(num_train <= 0)
    return 1;

  memset(ctx, 0, sizeof(*ctx));
  ctx->ll_mode = ll_mode;
  ctx->num_train = num_train;
  ctx->num_reg_tot = num_reg_tot;

  if(num_reg_tot <= 0){
    ctx->ready = 1;
    return 0;
  }

  ctx->vsfx = alloc_vecd(MAX(1, num_reg_tot));
  ctx->lambdax = alloc_vecd(MAX(1, num_reg_unordered_extern + num_reg_ordered_extern));
  ctx->kw = alloc_vecd(MAX(1, num_train));
  ctx->mean_row = alloc_vecd(MAX(1, num_train));
  ctx->matrix_bandwidth_x = alloc_tmatd(bw_rows, num_reg_continuous_extern);
  ctx->matrix_bandwidth_eval_one = alloc_tmatd(1, num_reg_continuous_extern);
  if(num_reg_unordered_extern > 0) ctx->eval_xuno_one = alloc_matd(1, num_reg_unordered_extern);
  if(num_reg_ordered_extern > 0) ctx->eval_xord_one = alloc_matd(1, num_reg_ordered_extern);
  if(num_reg_continuous_extern > 0) ctx->eval_xcon_one = alloc_matd(1, num_reg_continuous_extern);

  ctx->kernel_cx = (int *)calloc((size_t)MAX(1, num_reg_continuous_extern), sizeof(int));
  ctx->kernel_ux = (int *)calloc((size_t)MAX(1, num_reg_unordered_extern), sizeof(int));
  ctx->kernel_ox = (int *)calloc((size_t)MAX(1, num_reg_ordered_extern), sizeof(int));
  ctx->x_operator = (int *)calloc((size_t)MAX(1, num_reg_tot), sizeof(int));

  if((ctx->vsfx == NULL) || (ctx->lambdax == NULL) || (ctx->kw == NULL) || (ctx->mean_row == NULL) ||
     ((num_reg_continuous_extern > 0) && (ctx->matrix_bandwidth_x == NULL)) ||
     ((num_reg_continuous_extern > 0) && (ctx->matrix_bandwidth_eval_one == NULL)) ||
     ((num_reg_unordered_extern > 0) && (ctx->eval_xuno_one == NULL)) ||
     ((num_reg_ordered_extern > 0) && (ctx->eval_xord_one == NULL)) ||
     ((num_reg_continuous_extern > 0) && (ctx->eval_xcon_one == NULL)) ||
     (ctx->kernel_cx == NULL) || (ctx->kernel_ux == NULL) || (ctx->kernel_ox == NULL) || (ctx->x_operator == NULL))
    goto fail_xrow_ctx_prepare;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        ctx->vsfx,
                        NULL,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  for(i = 0; i < num_reg_continuous_extern; i++) ctx->kernel_cx[i] = KERNEL_reg_extern;
  for(i = 0; i < num_reg_unordered_extern; i++) ctx->kernel_ux[i] = KERNEL_reg_unordered_extern;
  for(i = 0; i < num_reg_ordered_extern; i++) ctx->kernel_ox[i] = KERNEL_reg_ordered_extern;
  for(i = 0; i < num_reg_tot; i++) ctx->x_operator[i] = OP_NORMAL;

  if(kernel_bandwidth_mean(KERNEL_reg_extern,
                           BANDWIDTH_den_extern,
                           num_train,
                           num_train,
                           0,
                           0,
                           0,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           1,
                           ctx->vsfx,
                           NULL,
                           NULL,
                           matrix_X_continuous_train_extern,
                           matrix_X_continuous_train_extern,
                           NULL,
                           ctx->matrix_bandwidth_x,
                           ctx->lambdax) == 1)
    goto fail_xrow_ctx_prepare;

  if(ll_mode == LL_LP){
    const int use_bernstein = (int_glp_bernstein_extern != 0);

    if((vector_glp_degree_extern == NULL) || (num_reg_continuous_extern <= 0))
      goto fail_xrow_ctx_prepare;
    if(!np_glp_cv_cache.ready ||
       (np_glp_cv_cache.use_bernstein != use_bernstein) ||
       (np_glp_cv_cache.basis_mode != int_glp_basis_extern) ||
       (np_glp_cv_cache.num_obs != num_train) ||
       (np_glp_cv_cache.ncon != num_reg_continuous_extern) ||
       (np_glp_cv_cache.matrix_X_continuous_train_ptr != matrix_X_continuous_train_extern)){
      if(!np_glp_cv_cache_prepare(LL_LP,
                                  num_train,
                                  num_reg_continuous_extern,
                                  matrix_X_continuous_train_extern))
        goto fail_xrow_ctx_prepare;
    }
  }

  ctx->ready = 1;
  return 0;

fail_xrow_ctx_prepare:
  np_conditional_xrow_ctx_clear(ctx);
  return 1;
}

static int np_conditional_xrow_from_ctx(NPConditionalXRowCtx *ctx,
                                        int eval_idx,
                                        double *row_out){
  const int num_train = num_obs_train_extern;
  int eval_pos = eval_idx;
  MATRIX KWM = NULL, RHS = NULL, SOL = NULL;
  NPConditionalBoundState bounds_state;
  int j, l;
  int status = 1;

  if((ctx == NULL) || (!ctx->ready) || (row_out == NULL))
    return 1;
  if((eval_idx < 0) || (eval_idx >= num_train))
    return 1;

  memset(row_out, 0, (size_t)num_train*sizeof(double));

  if((int_TREE_X == NP_TREE_TRUE) && (ipt_lookup_extern_X != NULL))
    eval_pos = ipt_lookup_extern_X[eval_idx];

  if(ctx->num_reg_tot <= 0){
    const double w = (num_train > 1) ? 1.0/((double)(num_train - 1)) : 0.0;
    for(j = 0; j < num_train; j++)
      if(j != eval_idx)
        row_out[j] = w;
    return 0;
  }

  for(l = 0; l < num_reg_unordered_extern; l++)
    ctx->eval_xuno_one[l][0] = matrix_X_unordered_train_extern[l][eval_pos];
  for(l = 0; l < num_reg_ordered_extern; l++)
    ctx->eval_xord_one[l][0] = matrix_X_ordered_train_extern[l][eval_pos];
  for(l = 0; l < num_reg_continuous_extern; l++){
    ctx->eval_xcon_one[l][0] = matrix_X_continuous_train_extern[l][eval_pos];
    ctx->matrix_bandwidth_eval_one[l][0] =
      (BANDWIDTH_den_extern == BW_GEN_NN) ? ctx->matrix_bandwidth_x[l][eval_pos] : ctx->matrix_bandwidth_x[l][0];
  }

  np_conditional_push_bounds(int_cxker_bound_extern,
                             vector_cxkerlb_extern,
                             vector_cxkerub_extern,
                             &bounds_state);
  if(np_shadow_conditional_kernel_row_raw(ctx->kernel_cx,
                                          ctx->kernel_ux,
                                          ctx->kernel_ox,
                                          ctx->x_operator,
                                          BANDWIDTH_den_extern,
                                          num_train,
                                          num_reg_unordered_extern,
                                          num_reg_ordered_extern,
                                          num_reg_continuous_extern,
                                          matrix_X_unordered_train_extern,
                                          matrix_X_ordered_train_extern,
                                          matrix_X_continuous_train_extern,
                                          ctx->eval_xuno_one,
                                          ctx->eval_xord_one,
                                          ctx->eval_xcon_one,
                                          ctx->vsfx,
                                          1,
                                          ctx->matrix_bandwidth_x,
                                          ctx->matrix_bandwidth_eval_one,
                                          ctx->lambdax,
                                          num_categories_extern_X,
                                          matrix_categorical_vals_extern_X,
                                          int_TREE_X,
                                          kdt_extern_X,
                                          ctx->kw,
                                          ctx->mean_row) != 0){
    np_conditional_pop_bounds(&bounds_state);
    goto cleanup_xrow_from_ctx;
  }
  np_conditional_pop_bounds(&bounds_state);

  ctx->kw[eval_pos] = 0.0;

  if(ctx->ll_mode == LL_LC){
    double row_sum = 0.0;
    for(j = 0; j < num_train; j++)
      row_sum += ctx->kw[j];
    if(!(row_sum > DBL_MIN))
      goto cleanup_xrow_from_ctx;
    for(j = 0; j < num_train; j++){
      const int orig_j = (int_TREE_X == NP_TREE_TRUE) ? ipt_extern_X[j] : j;
      row_out[orig_j] = ctx->kw[j]/row_sum;
    }
  } else {
    const int k = np_glp_cv_cache.nterms;

    if((k <= 0) || (np_glp_cv_cache.basis == NULL))
      goto cleanup_xrow_from_ctx;

    KWM = mat_creat(k, k, UNDEFINED);
    RHS = mat_creat(k, 1, UNDEFINED);
    SOL = mat_creat(k, 1, UNDEFINED);
    if((KWM == NULL) || (RHS == NULL) || (SOL == NULL))
      goto cleanup_xrow_from_ctx;

    if(np_glp_qr_drop_row_bkcde(np_glp_cv_cache.basis,
                                num_train,
                                k,
                                ctx->kw,
                                eval_pos,
                                ctx->mean_row) != 0)
      goto cleanup_xrow_from_ctx;

    for(j = 0; j < num_train; j++){
      const int orig_j = (int_TREE_X == NP_TREE_TRUE) ? ipt_extern_X[j] : j;
      row_out[orig_j] = ctx->mean_row[j];
    }
  }

  status = 0;

cleanup_xrow_from_ctx:
  if(KWM != NULL) mat_free(KWM);
  if(RHS != NULL) mat_free(RHS);
  if(SOL != NULL) mat_free(SOL);
  return status;
}

int np_regression_lp_apply_matrix(double *vector_scale_factor,
                                  double **rhs_cols,
                                  int n_rhs,
                                  double *fitted_out){
  const int num_train = num_obs_train_extern;
  const int num_eval = num_obs_eval_extern;
  const int num_reg_tot = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
  const int bw_rows = (BANDWIDTH_den_extern == BW_FIXED) ? 1 : num_eval;
  const int use_bernstein = (int_glp_bernstein_extern != 0);
  const double epsilon = 1.0/(double)MAX(1, num_train);
  int *kernel_cx = NULL, *kernel_ux = NULL, *kernel_ox = NULL, *x_operator = NULL;
  double *vsfx = NULL, *lambdax = NULL;
  double **matrix_bandwidth_x = NULL;
  double **TCON = NULL, **TUNO = NULL, **TORD = NULL;
  MATRIX KWM = NULL, XTKY = NULL, DELTA = NULL;
  double **Ycols = NULL, **Wcols = NULL;
  double *out = NULL, *eval_basis = NULL;
  int i, j, l;
  int status = 1;
  int target_deriv = -1;
  int target_order = 0;

  if((vector_scale_factor == NULL) || (rhs_cols == NULL) || (fitted_out == NULL))
    return 1;
  if(int_ll_extern != LL_LP)
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN))
    return 1;
  if((num_train <= 0) || (num_eval <= 0) || (num_reg_continuous_extern <= 0) || (n_rhs <= 0))
    return 1;

  if(vector_glp_gradient_order_extern != NULL){
    for(i = 0; i < num_reg_continuous_extern; i++){
      const int ord = vector_glp_gradient_order_extern[i];
      if(ord > 0){
        if(target_deriv >= 0)
          return 1;
        target_deriv = i;
        target_order = ord;
      }
    }
  }

  memset(fitted_out, 0, (size_t)num_eval*(size_t)n_rhs*sizeof(double));

  vsfx = alloc_vecd(MAX(1, num_reg_tot));
  lambdax = alloc_vecd(MAX(1, num_reg_unordered_extern + num_reg_ordered_extern));
  matrix_bandwidth_x = alloc_tmatd(bw_rows, num_reg_continuous_extern);
  TCON = alloc_matd(1, num_reg_continuous_extern);
  TUNO = alloc_matd(1, num_reg_unordered_extern);
  TORD = alloc_matd(1, num_reg_ordered_extern);
  kernel_cx = (int *)calloc((size_t)MAX(1, num_reg_continuous_extern), sizeof(int));
  kernel_ux = (int *)calloc((size_t)MAX(1, num_reg_unordered_extern), sizeof(int));
  kernel_ox = (int *)calloc((size_t)MAX(1, num_reg_ordered_extern), sizeof(int));
  x_operator = (int *)calloc((size_t)MAX(1, num_reg_tot), sizeof(int));

  if((vsfx == NULL) || (lambdax == NULL) ||
     ((num_reg_continuous_extern > 0) && (matrix_bandwidth_x == NULL)) ||
     ((num_reg_continuous_extern > 0) && (TCON == NULL)) ||
     ((num_reg_unordered_extern > 0) && (TUNO == NULL)) ||
     ((num_reg_ordered_extern > 0) && (TORD == NULL)) ||
     (kernel_cx == NULL) || (kernel_ux == NULL) || (kernel_ox == NULL) || (x_operator == NULL))
    goto cleanup_lp_apply;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        vsfx,
                        NULL,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  for(i = 0; i < num_reg_continuous_extern; i++) kernel_cx[i] = KERNEL_reg_extern;
  for(i = 0; i < num_reg_unordered_extern; i++) kernel_ux[i] = KERNEL_reg_unordered_extern;
  for(i = 0; i < num_reg_ordered_extern; i++) kernel_ox[i] = KERNEL_reg_ordered_extern;
  for(i = 0; i < num_reg_tot; i++) x_operator[i] = OP_NORMAL;

  if(kernel_bandwidth_mean(KERNEL_reg_extern,
                           BANDWIDTH_den_extern,
                           num_train,
                           num_eval,
                           0,
                           0,
                           0,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           1,
                           vsfx,
                           NULL,
                           NULL,
                           matrix_X_continuous_train_extern,
                           matrix_X_continuous_eval_extern,
                           NULL,
                           matrix_bandwidth_x,
                           lambdax) == 1)
    goto cleanup_lp_apply;

  if(!np_glp_cv_cache.ready ||
     (np_glp_cv_cache.use_bernstein != use_bernstein) ||
     (np_glp_cv_cache.basis_mode != int_glp_basis_extern) ||
     (np_glp_cv_cache.num_obs != num_train) ||
     (np_glp_cv_cache.ncon != num_reg_continuous_extern) ||
     (np_glp_cv_cache.matrix_X_continuous_train_ptr != matrix_X_continuous_train_extern)){
    if(!np_glp_cv_cache_prepare(LL_LP,
                                num_train,
                                num_reg_continuous_extern,
                                matrix_X_continuous_train_extern))
      goto cleanup_lp_apply;
  }

  if((np_glp_cv_cache.nterms <= 0) || (np_glp_cv_cache.basis == NULL) || (np_glp_cv_cache.terms == NULL))
    goto cleanup_lp_apply;

  KWM = mat_creat(np_glp_cv_cache.nterms, np_glp_cv_cache.nterms, UNDEFINED);
  XTKY = mat_creat(np_glp_cv_cache.nterms, n_rhs, UNDEFINED);
  DELTA = mat_creat(np_glp_cv_cache.nterms, n_rhs, UNDEFINED);
  Ycols = (double **)malloc((size_t)(n_rhs + np_glp_cv_cache.nterms)*sizeof(double *));
  Wcols = (double **)malloc((size_t)np_glp_cv_cache.nterms*sizeof(double *));
  out = (double *)malloc((size_t)(n_rhs + np_glp_cv_cache.nterms)*(size_t)np_glp_cv_cache.nterms*sizeof(double));
  eval_basis = (double *)malloc((size_t)np_glp_cv_cache.nterms*sizeof(double));

  if((KWM == NULL) || (XTKY == NULL) || (DELTA == NULL) ||
     (Ycols == NULL) || (Wcols == NULL) || (out == NULL) || (eval_basis == NULL))
    goto cleanup_lp_apply;

  for(i = 0; i < n_rhs; i++)
    Ycols[i] = rhs_cols[i];
  for(l = 0; l < np_glp_cv_cache.nterms; l++){
    Ycols[n_rhs + l] = np_glp_cv_cache.basis[l];
    Wcols[l] = np_glp_cv_cache.basis[l];
  }

  for(j = 0; j < num_eval; j++){
    double nepsilon = 0.0;

    for(l = 0; l < num_reg_continuous_extern; l++)
      TCON[l][0] = matrix_X_continuous_eval_extern[l][j];
    for(l = 0; l < num_reg_unordered_extern; l++)
      TUNO[l][0] = matrix_X_unordered_eval_extern[l][j];
    for(l = 0; l < num_reg_ordered_extern; l++)
      TORD[l][0] = matrix_X_ordered_eval_extern[l][j];

    if(kernel_weighted_sum_np_ctx(kernel_cx,
                                  kernel_ux,
                                  kernel_ox,
                                  BANDWIDTH_den_extern,
                                  num_train,
                                  1,
                                  num_reg_unordered_extern,
                                  num_reg_ordered_extern,
                                  num_reg_continuous_extern,
                                  0,
                                  0,
                                  1,
                                  1,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  x_operator,
                                  OP_NOOP,
                                  0,
                                  0,
                                  NULL,
                                  1,
                                  n_rhs + np_glp_cv_cache.nterms,
                                  np_glp_cv_cache.nterms,
                                  int_TREE_X,
                                  0,
                                  kdt_extern_X,
                                  NULL, NULL, NULL,
                                  matrix_X_unordered_train_extern,
                                  matrix_X_ordered_train_extern,
                                  matrix_X_continuous_train_extern,
                                  TUNO,
                                  TORD,
                                  TCON,
                                  Ycols,
                                  Wcols,
                                  NULL,
                                  vector_scale_factor,
                                  1,
                                  matrix_bandwidth_x,
                                  matrix_bandwidth_x,
                                  lambdax,
                                  num_categories_extern_X,
                                  matrix_categorical_vals_extern_X,
                                  NULL,
                                  out,
                                  NULL,
                                  NULL,
                                  NULL) != 0)
      goto cleanup_lp_apply;

    for(i = 0; i < np_glp_cv_cache.nterms; i++){
      const int base = i*(n_rhs + np_glp_cv_cache.nterms);
      for(int rhs_idx = 0; rhs_idx < n_rhs; rhs_idx++)
        XTKY[i][rhs_idx] = out[base + rhs_idx];
      for(l = 0; l < np_glp_cv_cache.nterms; l++)
        KWM[i][l] = out[base + n_rhs + l];
    }

    while(mat_solve(KWM, XTKY, DELTA) == NULL){
      for(i = 0; i < np_glp_cv_cache.nterms; i++)
        KWM[i][i] += epsilon;
      nepsilon += epsilon;
    }

    if(nepsilon > 0.0){
      for(int rhs_idx = 0; rhs_idx < n_rhs; rhs_idx++)
        XTKY[0][rhs_idx] += nepsilon*XTKY[0][rhs_idx]/NZD_POS(KWM[0][0]);
      if(mat_solve(KWM, XTKY, DELTA) == NULL)
        goto cleanup_lp_apply;
    }

    if(target_deriv >= 0){
      if(use_bernstein)
        np_glp_fill_basis_eval_deriv(target_deriv,
                                     target_order,
                                     num_reg_continuous_extern,
                                     np_glp_cv_cache.terms,
                                     np_glp_cv_cache.nterms,
                                     matrix_X_continuous_eval_extern,
                                     j,
                                     np_glp_cv_cache.basis_ctx,
                                     eval_basis);
      else
        np_glp_fill_basis_eval_deriv_raw(target_deriv,
                                         target_order,
                                         num_reg_continuous_extern,
                                         np_glp_cv_cache.terms,
                                         np_glp_cv_cache.nterms,
                                         matrix_X_continuous_eval_extern,
                                         j,
                                         eval_basis);
    } else if(use_bernstein) {
      np_glp_fill_basis_eval(num_reg_continuous_extern,
                             np_glp_cv_cache.terms,
                             np_glp_cv_cache.nterms,
                             matrix_X_continuous_eval_extern,
                             j,
                             np_glp_cv_cache.basis_ctx,
                             eval_basis);
    } else {
      np_glp_fill_basis_eval_raw(num_reg_continuous_extern,
                                 np_glp_cv_cache.terms,
                                 np_glp_cv_cache.nterms,
                                 matrix_X_continuous_eval_extern,
                                 j,
                                 eval_basis);
    }

    for(int rhs_idx = 0; rhs_idx < n_rhs; rhs_idx++){
      double fit = 0.0;
      for(i = 0; i < np_glp_cv_cache.nterms; i++)
        fit += eval_basis[i]*DELTA[i][rhs_idx];
      fitted_out[(size_t)j + (size_t)num_eval*(size_t)rhs_idx] = fit;
    }
  }

  status = 0;

cleanup_lp_apply:
  if(KWM != NULL) mat_free(KWM);
  if(XTKY != NULL) mat_free(XTKY);
  if(DELTA != NULL) mat_free(DELTA);
  if(matrix_bandwidth_x != NULL) free_tmat(matrix_bandwidth_x);
  if((TCON != NULL) && (num_reg_continuous_extern > 0)) free_mat(TCON, num_reg_continuous_extern);
  if((TUNO != NULL) && (num_reg_unordered_extern > 0)) free_mat(TUNO, num_reg_unordered_extern);
  if((TORD != NULL) && (num_reg_ordered_extern > 0)) free_mat(TORD, num_reg_ordered_extern);
  if(kernel_cx != NULL) free(kernel_cx);
  if(kernel_ux != NULL) free(kernel_ux);
  if(kernel_ox != NULL) free(kernel_ox);
  if(x_operator != NULL) free(x_operator);
  if(vsfx != NULL) free(vsfx);
  if(lambdax != NULL) free(lambdax);
  if(Ycols != NULL) free(Ycols);
  if(Wcols != NULL) free(Wcols);
  if(out != NULL) free(out);
  if(eval_basis != NULL) free(eval_basis);
  return status;
}

static void np_conditional_yrow_ctx_clear(NPConditionalYRowCtx *ctx){
  if(ctx == NULL)
    return;
  if(ctx->vsfy != NULL) free(ctx->vsfy);
  if(ctx->lambday != NULL) free(ctx->lambday);
  if(ctx->kw != NULL) free(ctx->kw);
  if(ctx->matrix_bandwidth_y != NULL) free_tmat(ctx->matrix_bandwidth_y);
  if(ctx->matrix_bandwidth_eval_one != NULL) free_tmat(ctx->matrix_bandwidth_eval_one);
  if(ctx->eval_yuno_one != NULL) free_mat(ctx->eval_yuno_one, num_var_unordered_extern);
  if(ctx->eval_yord_one != NULL) free_mat(ctx->eval_yord_one, num_var_ordered_extern);
  if(ctx->eval_ycon_one != NULL) free_mat(ctx->eval_ycon_one, num_var_continuous_extern);
  if(ctx->kernel_cy != NULL) free(ctx->kernel_cy);
  if(ctx->kernel_uy != NULL) free(ctx->kernel_uy);
  if(ctx->kernel_oy != NULL) free(ctx->kernel_oy);
  if(ctx->operator_y != NULL) free(ctx->operator_y);
  memset(ctx, 0, sizeof(*ctx));
}

static int np_conditional_yrow_ctx_prepare(double *vector_scale_factor,
                                           int operator_code,
                                           NPConditionalYRowCtx *ctx){
  const int num_train = num_obs_train_extern;
  const int num_var_tot = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  const int bw_rows = (BANDWIDTH_den_extern == BW_FIXED) ? 1 : num_train;
  int i;

  if((ctx == NULL) || (vector_scale_factor == NULL))
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN))
    return 1;
  if(num_train <= 0)
    return 1;

  memset(ctx, 0, sizeof(*ctx));
  ctx->num_train = num_train;
  ctx->num_var_tot = num_var_tot;

  if(num_var_tot <= 0){
    ctx->ready = 1;
    return 0;
  }

  ctx->vsfy = alloc_vecd(MAX(1, num_var_tot));
  ctx->lambday = alloc_vecd(MAX(1, num_var_unordered_extern + num_var_ordered_extern));
  ctx->kw = alloc_vecd(MAX(1, num_train));
  ctx->matrix_bandwidth_y = alloc_tmatd(bw_rows, num_var_continuous_extern);
  ctx->matrix_bandwidth_eval_one = alloc_tmatd(1, num_var_continuous_extern);
  if(num_var_unordered_extern > 0) ctx->eval_yuno_one = alloc_matd(1, num_var_unordered_extern);
  if(num_var_ordered_extern > 0) ctx->eval_yord_one = alloc_matd(1, num_var_ordered_extern);
  if(num_var_continuous_extern > 0) ctx->eval_ycon_one = alloc_matd(1, num_var_continuous_extern);

  ctx->kernel_cy = (int *)calloc((size_t)MAX(1, num_var_continuous_extern), sizeof(int));
  ctx->kernel_uy = (int *)calloc((size_t)MAX(1, num_var_unordered_extern), sizeof(int));
  ctx->kernel_oy = (int *)calloc((size_t)MAX(1, num_var_ordered_extern), sizeof(int));
  ctx->operator_y = (int *)calloc((size_t)MAX(1, num_var_tot), sizeof(int));

  if((ctx->vsfy == NULL) || (ctx->lambday == NULL) || (ctx->kw == NULL) ||
     ((num_var_continuous_extern > 0) && (ctx->matrix_bandwidth_y == NULL)) ||
     ((num_var_continuous_extern > 0) && (ctx->matrix_bandwidth_eval_one == NULL)) ||
     ((num_var_unordered_extern > 0) && (ctx->eval_yuno_one == NULL)) ||
     ((num_var_ordered_extern > 0) && (ctx->eval_yord_one == NULL)) ||
     ((num_var_continuous_extern > 0) && (ctx->eval_ycon_one == NULL)) ||
     (ctx->kernel_cy == NULL) || (ctx->kernel_uy == NULL) || (ctx->kernel_oy == NULL) || (ctx->operator_y == NULL))
    goto fail_yrow_ctx_prepare;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        NULL,
                        ctx->vsfy,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  for(i = 0; i < num_var_continuous_extern; i++) ctx->kernel_cy[i] = KERNEL_den_extern;
  for(i = 0; i < num_var_unordered_extern; i++) ctx->kernel_uy[i] = KERNEL_den_unordered_extern;
  for(i = 0; i < num_var_ordered_extern; i++) ctx->kernel_oy[i] = KERNEL_den_ordered_extern;
  for(i = 0; i < num_var_tot; i++) ctx->operator_y[i] = operator_code;

  if(kernel_bandwidth_mean(KERNEL_den_extern,
                           BANDWIDTH_den_extern,
                           num_train,
                           num_train,
                           0,
                           0,
                           0,
                           num_var_continuous_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           1,
                           ctx->vsfy,
                           NULL,
                           NULL,
                           matrix_Y_continuous_train_extern,
                           matrix_Y_continuous_train_extern,
                           NULL,
                           ctx->matrix_bandwidth_y,
                           ctx->lambday) == 1)
    goto fail_yrow_ctx_prepare;

  ctx->ready = 1;
  return 0;

fail_yrow_ctx_prepare:
  np_conditional_yrow_ctx_clear(ctx);
  return 1;
}

static int np_conditional_y_scalar_fixed_row_direct(NPConditionalYRowCtx *ctx,
                                                    double eval_y,
                                                    double *row_out);

static void np_conditional_push_bounds(int new_bound_flag,
                                       double *new_lb,
                                       double *new_ub,
                                       NPConditionalBoundState *saved){
  if(saved == NULL)
    return;
  saved->int_cker_bound = int_cker_bound_extern;
  saved->ckerlb = vector_ckerlb_extern;
  saved->ckerub = vector_ckerub_extern;
  int_cker_bound_extern = new_bound_flag;
  vector_ckerlb_extern = new_lb;
  vector_ckerub_extern = new_ub;
}

static void np_conditional_pop_bounds(const NPConditionalBoundState *saved){
  if(saved == NULL)
    return;
  int_cker_bound_extern = saved->int_cker_bound;
  vector_ckerlb_extern = saved->ckerlb;
  vector_ckerub_extern = saved->ckerub;
}

static int np_conditional_yrow_from_ctx(NPConditionalYRowCtx *ctx,
                                        int eval_idx,
                                        double *row_out){
  const int num_train = num_obs_train_extern;
  int eval_pos = eval_idx;
  NPConditionalBoundState bounds_state;
  int j, l;

  if((ctx == NULL) || (!ctx->ready) || (row_out == NULL))
    return 1;
  if((eval_idx < 0) || (eval_idx >= num_train))
    return 1;

  memset(row_out, 0, (size_t)num_train*sizeof(double));

  if((int_TREE_Y == NP_TREE_TRUE) && (ipt_lookup_extern_Y != NULL))
    eval_pos = ipt_lookup_extern_Y[eval_idx];

  if(ctx->num_var_tot <= 0){
    for(j = 0; j < num_train; j++)
      row_out[j] = 1.0;
    return 0;
  }

  if((BANDWIDTH_den_extern == BW_FIXED) &&
     (int_TREE_Y != NP_TREE_TRUE) &&
     (ctx->num_var_tot == 1) &&
     (num_var_continuous_extern == 1) &&
     (num_var_unordered_extern == 0) &&
     (num_var_ordered_extern == 0))
    return np_conditional_y_scalar_fixed_row_direct(ctx,
                                                    matrix_Y_continuous_train_extern[0][eval_idx],
                                                    row_out);

  for(l = 0; l < num_var_unordered_extern; l++)
    ctx->eval_yuno_one[l][0] = matrix_Y_unordered_train_extern[l][eval_pos];
  for(l = 0; l < num_var_ordered_extern; l++)
    ctx->eval_yord_one[l][0] = matrix_Y_ordered_train_extern[l][eval_pos];
  for(l = 0; l < num_var_continuous_extern; l++){
    ctx->eval_ycon_one[l][0] = matrix_Y_continuous_train_extern[l][eval_pos];
    ctx->matrix_bandwidth_eval_one[l][0] =
      (BANDWIDTH_den_extern == BW_GEN_NN) ? ctx->matrix_bandwidth_y[l][eval_pos] : ctx->matrix_bandwidth_y[l][0];
  }

  np_conditional_push_bounds(int_cyker_bound_extern,
                             vector_cykerlb_extern,
                             vector_cykerub_extern,
                             &bounds_state);
  if(np_shadow_conditional_kernel_row(ctx->kernel_cy,
                                      ctx->kernel_uy,
                                      ctx->kernel_oy,
                                      ctx->operator_y,
                                      BANDWIDTH_den_extern,
                                      num_train,
                                      num_var_unordered_extern,
                                      num_var_ordered_extern,
                                      num_var_continuous_extern,
                                      matrix_Y_unordered_train_extern,
                                      matrix_Y_ordered_train_extern,
                                      matrix_Y_continuous_train_extern,
                                      ctx->eval_yuno_one,
                                      ctx->eval_yord_one,
                                      ctx->eval_ycon_one,
                                      ctx->vsfy,
                                      1,
                                      ctx->matrix_bandwidth_y,
                                      ctx->matrix_bandwidth_eval_one,
                                      ctx->lambday,
                                      num_categories_extern_Y,
                                      matrix_categorical_vals_extern_Y,
                                      int_TREE_Y,
                                      kdt_extern_Y,
                                      ctx->kw,
                                      NULL) != 0){
    np_conditional_pop_bounds(&bounds_state);
    return 1;
  }
  np_conditional_pop_bounds(&bounds_state);

  for(j = 0; j < num_train; j++){
    const int orig_j = (int_TREE_Y == NP_TREE_TRUE) ? ipt_extern_Y[j] : j;
    row_out[orig_j] = ctx->kw[j];
  }

  return 0;
}

static int np_conditional_y_eval_row_stream_op_core(double *vector_scale_factor,
                                                    int eval_idx,
                                                    int operator_code,
                                                    double **matrix_Y_unordered_eval,
                                                    double **matrix_Y_ordered_eval,
                                                    double **matrix_Y_continuous_eval,
                                                    int num_eval,
                                                    int map_train_tree_index,
                                                    double *row_out);

static int np_conditional_y_scalar_fixed_row_direct(NPConditionalYRowCtx *ctx,
                                                    double eval_y,
                                                    double *row_out){
  const int num_train = num_obs_train_extern;
  const double *train_y = (matrix_Y_continuous_train_extern != NULL) ?
    matrix_Y_continuous_train_extern[0] : NULL;
  const double lb = (vector_cykerlb_extern != NULL) ? vector_cykerlb_extern[0] : R_NegInf;
  const double ub = (vector_cykerub_extern != NULL) ? vector_cykerub_extern[0] : R_PosInf;
  double h, invnorm;
  int all_largeh = 0;
  double max_largeh_abs = 0.0, largeh_k0 = 0.0;
  int j;

  if((ctx == NULL) || (!ctx->ready) || (row_out == NULL) || (train_y == NULL))
    return 1;
  if(BANDWIDTH_den_extern != BW_FIXED)
    return 1;
  if(int_TREE_Y == NP_TREE_TRUE)
    return 1;
  if((num_var_continuous_extern != 1) ||
     (num_var_unordered_extern != 0) ||
     (num_var_ordered_extern != 0))
    return 1;
  if((ctx->matrix_bandwidth_y == NULL) || (ctx->matrix_bandwidth_y[0] == NULL))
    return 1;

  h = ctx->matrix_bandwidth_y[0][0];
  if(!(h > 0.0) || !R_FINITE(h) || !R_FINITE(eval_y))
    return 1;

  invnorm = np_cker_invnorm(KERNEL_den_extern, eval_y, h, lb, ub);

  if(np_cont_largeh_kernel_supported(KERNEL_den_extern)){
    const double utol = np_cont_largeh_utol(KERNEL_den_extern, np_cont_largeh_rel_tol());
    if((utol > 0.0) && R_FINITE(utol)){
      all_largeh = 1;
      max_largeh_abs = utol*fabs(h);
      largeh_k0 = np_cont_largeh_k0(KERNEL_den_extern);
    }
  }

  for(j = 0; j < num_train; j++){
    const double diff = eval_y - train_y[j];
    if(all_largeh && (fabs(diff) > max_largeh_abs))
      all_largeh = 0;
    row_out[j] = invnorm*np_cker_base_eval(KERNEL_den_extern, diff/h)/h;
  }

  if(all_largeh){
    const double kval = invnorm*largeh_k0/h;
    for(j = 0; j < num_train; j++)
      row_out[j] = kval;
  }

  return 0;
}

static int np_conditional_y_scalar_eval_from_ctx(double *vector_scale_factor,
                                                 NPConditionalYRowCtx *ctx,
                                                 double eval_y,
                                                 double *row_out){
  const int num_train = num_obs_train_extern;
  NPConditionalBoundState bounds_state;
  int j;

  if((ctx == NULL) || (!ctx->ready) || (row_out == NULL) || (vector_scale_factor == NULL))
    return 1;
  if(num_var_unordered_extern != 0)
    return 1;
  if(num_var_ordered_extern != 0)
    return 1;
  if(num_var_continuous_extern != 1)
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;

  if(np_conditional_y_scalar_fixed_row_direct(ctx, eval_y, row_out) == 0)
    return 0;

  ctx->eval_ycon_one[0][0] = eval_y;

  if(BANDWIDTH_den_extern != BW_FIXED){
    return np_conditional_y_eval_row_stream_op_core(vector_scale_factor,
                                                    0,
                                                    OP_NORMAL,
                                                    NULL,
                                                    NULL,
                                                    ctx->eval_ycon_one,
                                                    1,
                                                    0,
                                                    row_out);
  }

  memset(row_out, 0, (size_t)num_train*sizeof(double));

  ctx->matrix_bandwidth_eval_one[0][0] = ctx->matrix_bandwidth_y[0][0];

  np_conditional_push_bounds(int_cyker_bound_extern,
                             vector_cykerlb_extern,
                             vector_cykerub_extern,
                             &bounds_state);
  if(np_shadow_conditional_kernel_row(ctx->kernel_cy,
                                      ctx->kernel_uy,
                                      ctx->kernel_oy,
                                      ctx->operator_y,
                                      BANDWIDTH_den_extern,
                                      num_train,
                                      num_var_unordered_extern,
                                      num_var_ordered_extern,
                                      num_var_continuous_extern,
                                      matrix_Y_unordered_train_extern,
                                      matrix_Y_ordered_train_extern,
                                      matrix_Y_continuous_train_extern,
                                      ctx->eval_yuno_one,
                                      ctx->eval_yord_one,
                                      ctx->eval_ycon_one,
                                      ctx->vsfy,
                                      1,
                                      ctx->matrix_bandwidth_y,
                                      ctx->matrix_bandwidth_eval_one,
                                      ctx->lambday,
                                      num_categories_extern_Y,
                                      matrix_categorical_vals_extern_Y,
                                      int_TREE_Y,
                                      kdt_extern_Y,
                                      ctx->kw,
                                      NULL) != 0){
    np_conditional_pop_bounds(&bounds_state);
    return 1;
  }
  np_conditional_pop_bounds(&bounds_state);

  for(j = 0; j < num_train; j++){
    const int orig_j = (int_TREE_Y == NP_TREE_TRUE) ? ipt_extern_Y[j] : j;
    row_out[orig_j] = ctx->kw[j];
  }

  return 0;
}

static int np_conditional_x_weight_row_stream_core_impl(double *vector_scale_factor,
                                                        int eval_idx,
                                                        int drop_eval_self,
                                                        double *row_out){
  const int num_train = num_obs_train_extern;
  const int num_reg_tot = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
  const int ll_mode = (int_ll_extern == LL_LP) ? LL_LP : LL_LC;
  const int bw_rows = (BANDWIDTH_den_extern == BW_FIXED) ? 1 : num_train;
  int *kernel_cx = NULL, *kernel_ux = NULL, *kernel_ox = NULL, *x_operator = NULL;
  double *vsfx = NULL, *lambdax = NULL, *kw = NULL, *mean_row = NULL;
  double **matrix_bandwidth_x = NULL, **matrix_bandwidth_eval_one = NULL;
  double **eval_xuno_one = NULL, **eval_xord_one = NULL, **eval_xcon_one = NULL;
  MATRIX KWM = NULL, RHS = NULL, SOL = NULL;
  NPConditionalBoundState bounds_state;
  int eval_pos = eval_idx;
  int i, j, l;
  int status = 1;

  if((row_out == NULL) || (vector_scale_factor == NULL))
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;
  if((eval_idx < 0) || (eval_idx >= num_train))
    return 1;

  memset(row_out, 0, (size_t)num_train*sizeof(double));

  if((int_TREE_X == NP_TREE_TRUE) && (ipt_lookup_extern_X != NULL))
    eval_pos = ipt_lookup_extern_X[eval_idx];

  if(num_reg_tot <= 0){
    const double denom = drop_eval_self ? ((double)(num_train - 1)) : ((double)num_train);
    const double w = (denom > 0.0) ? 1.0/denom : 0.0;
    for(j = 0; j < num_train; j++)
      if((!drop_eval_self) || (j != eval_idx))
        row_out[j] = w;
    return 0;
  }

  vsfx = alloc_vecd(MAX(1, num_reg_tot));
  lambdax = alloc_vecd(MAX(1, num_reg_unordered_extern + num_reg_ordered_extern));
  kw = alloc_vecd(MAX(1, num_train));
  mean_row = alloc_vecd(MAX(1, num_train));
  matrix_bandwidth_x = alloc_tmatd(bw_rows, num_reg_continuous_extern);
  matrix_bandwidth_eval_one = alloc_tmatd(1, num_reg_continuous_extern);
  if(num_reg_unordered_extern > 0) eval_xuno_one = alloc_matd(1, num_reg_unordered_extern);
  if(num_reg_ordered_extern > 0) eval_xord_one = alloc_matd(1, num_reg_ordered_extern);
  if(num_reg_continuous_extern > 0) eval_xcon_one = alloc_matd(1, num_reg_continuous_extern);

  kernel_cx = (int *)calloc((size_t)MAX(1, num_reg_continuous_extern), sizeof(int));
  kernel_ux = (int *)calloc((size_t)MAX(1, num_reg_unordered_extern), sizeof(int));
  kernel_ox = (int *)calloc((size_t)MAX(1, num_reg_ordered_extern), sizeof(int));
  x_operator = (int *)calloc((size_t)MAX(1, num_reg_tot), sizeof(int));

  if((vsfx == NULL) || (lambdax == NULL) || (kw == NULL) || (mean_row == NULL) ||
     ((num_reg_continuous_extern > 0) && (matrix_bandwidth_x == NULL)) ||
     ((num_reg_continuous_extern > 0) && (matrix_bandwidth_eval_one == NULL)) ||
     ((num_reg_unordered_extern > 0) && (eval_xuno_one == NULL)) ||
     ((num_reg_ordered_extern > 0) && (eval_xord_one == NULL)) ||
     ((num_reg_continuous_extern > 0) && (eval_xcon_one == NULL)) ||
     (kernel_cx == NULL) || (kernel_ux == NULL) || (kernel_ox == NULL) || (x_operator == NULL))
    goto cleanup_xweight_row;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        vsfx,
                        NULL,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  for(i = 0; i < num_reg_continuous_extern; i++) kernel_cx[i] = KERNEL_reg_extern;
  for(i = 0; i < num_reg_unordered_extern; i++) kernel_ux[i] = KERNEL_reg_unordered_extern;
  for(i = 0; i < num_reg_ordered_extern; i++) kernel_ox[i] = KERNEL_reg_ordered_extern;
  for(i = 0; i < num_reg_tot; i++) x_operator[i] = OP_NORMAL;

  if(kernel_bandwidth_mean(KERNEL_reg_extern,
                           BANDWIDTH_den_extern,
                           num_train,
                           num_train,
                           0,
                           0,
                           0,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           vsfx,
                           NULL,
                           NULL,
                           matrix_X_continuous_train_extern,
                           matrix_X_continuous_train_extern,
                           NULL,
                           matrix_bandwidth_x,
                           lambdax) == 1)
    goto cleanup_xweight_row;

  if(ll_mode == LL_LP){
    if((vector_glp_degree_extern == NULL) || (num_reg_continuous_extern <= 0))
      goto cleanup_xweight_row;
    if(!np_glp_cv_prepare_extern(LL_LP,
                                 num_train,
                                 num_reg_continuous_extern,
                                 matrix_X_continuous_train_extern))
      goto cleanup_xweight_row;
  }

  for(l = 0; l < num_reg_unordered_extern; l++)
    eval_xuno_one[l][0] = matrix_X_unordered_train_extern[l][eval_pos];
  for(l = 0; l < num_reg_ordered_extern; l++)
    eval_xord_one[l][0] = matrix_X_ordered_train_extern[l][eval_pos];
  for(l = 0; l < num_reg_continuous_extern; l++){
    eval_xcon_one[l][0] = matrix_X_continuous_train_extern[l][eval_pos];
    matrix_bandwidth_eval_one[l][0] =
      (BANDWIDTH_den_extern == BW_GEN_NN) ? matrix_bandwidth_x[l][eval_pos] : matrix_bandwidth_x[l][0];
  }

  np_conditional_push_bounds(int_cxker_bound_extern,
                             vector_cxkerlb_extern,
                             vector_cxkerub_extern,
                             &bounds_state);
  if(np_shadow_conditional_kernel_row_raw(kernel_cx,
                                          kernel_ux,
                                          kernel_ox,
                                          x_operator,
                                          BANDWIDTH_den_extern,
                                          num_train,
                                          num_reg_unordered_extern,
                                          num_reg_ordered_extern,
                                          num_reg_continuous_extern,
                                          matrix_X_unordered_train_extern,
                                          matrix_X_ordered_train_extern,
                                          matrix_X_continuous_train_extern,
                                          eval_xuno_one,
                                          eval_xord_one,
                                          eval_xcon_one,
                                          vsfx,
                                          1,
                                          matrix_bandwidth_x,
                                          matrix_bandwidth_eval_one,
                                          lambdax,
                                          num_categories_extern_X,
                                          matrix_categorical_vals_extern_X,
                                          int_TREE_X,
                                          kdt_extern_X,
                                          kw,
                                          mean_row) != 0){
    np_conditional_pop_bounds(&bounds_state);
    goto cleanup_xweight_row;
  }
  np_conditional_pop_bounds(&bounds_state);

  if(drop_eval_self)
    kw[eval_pos] = 0.0;

  if(ll_mode == LL_LC){
    double row_sum = 0.0;
    for(j = 0; j < num_train; j++)
      row_sum += kw[j];
    if(!(row_sum > DBL_MIN))
      goto cleanup_xweight_row;
    for(j = 0; j < num_train; j++){
      const int orig_j = (int_TREE_X == NP_TREE_TRUE) ? ipt_extern_X[j] : j;
      row_out[orig_j] = kw[j]/row_sum;
    }
  } else {
    const int k = np_glp_cv_cache.nterms;

    if((k <= 0) || (np_glp_cv_cache.basis == NULL))
      goto cleanup_xweight_row;

    KWM = mat_creat(k, k, UNDEFINED);
    RHS = mat_creat(k, 1, UNDEFINED);
    SOL = mat_creat(k, 1, UNDEFINED);
    if((KWM == NULL) || (RHS == NULL) || (SOL == NULL))
      goto cleanup_xweight_row;

    if(drop_eval_self){
      if(np_glp_qr_drop_row_bkcde(np_glp_cv_cache.basis,
                                  num_train,
                                  k,
                                  kw,
                                  eval_pos,
                                  mean_row) != 0)
        goto cleanup_xweight_row;

      for(j = 0; j < num_train; j++){
        const int orig_j = (int_TREE_X == NP_TREE_TRUE) ? ipt_extern_X[j] : j;
        row_out[orig_j] = mean_row[j];
      }
    } else {
      for(l = 0; l < k; l++){
        RHS[l][0] = np_glp_cv_cache.basis[l][eval_pos];
        for(j = 0; j < k; j++)
          KWM[l][j] = 0.0;
      }

      for(j = 0; j < num_train; j++){
        const double wj = kw[j];
        if(wj == 0.0)
          continue;
        for(int a = 0; a < k; a++){
          const double za = np_glp_cv_cache.basis[a][j];
          for(int b = a; b < k; b++){
            const double zb = np_glp_cv_cache.basis[b][j];
            KWM[a][b] += wj*za*zb;
            if(b != a) KWM[b][a] += wj*za*zb;
          }
        }
      }

      if(np_mat_bad_rcond_sym(KWM, 1.0e-10))
        goto cleanup_xweight_row;
      if(mat_solve(KWM, RHS, SOL) == NULL)
        goto cleanup_xweight_row;

      for(j = 0; j < num_train; j++){
        double zju = 0.0;
        const int orig_j = (int_TREE_X == NP_TREE_TRUE) ? ipt_extern_X[j] : j;
        for(l = 0; l < k; l++)
          zju += np_glp_cv_cache.basis[l][j]*SOL[l][0];
        row_out[orig_j] = kw[j]*zju;
      }
    }
  }

  status = 0;

cleanup_xweight_row:
  if(KWM != NULL) mat_free(KWM);
  if(RHS != NULL) mat_free(RHS);
  if(SOL != NULL) mat_free(SOL);
  if(vsfx != NULL) free(vsfx);
  if(lambdax != NULL) free(lambdax);
  if(kw != NULL) free(kw);
  if(mean_row != NULL) free(mean_row);
  if(matrix_bandwidth_x != NULL) free_tmat(matrix_bandwidth_x);
  if(matrix_bandwidth_eval_one != NULL) free_tmat(matrix_bandwidth_eval_one);
  if(eval_xuno_one != NULL) free_mat(eval_xuno_one, num_reg_unordered_extern);
  if(eval_xord_one != NULL) free_mat(eval_xord_one, num_reg_ordered_extern);
  if(eval_xcon_one != NULL) free_mat(eval_xcon_one, num_reg_continuous_extern);
  if(kernel_cx != NULL) free(kernel_cx);
  if(kernel_ux != NULL) free(kernel_ux);
  if(kernel_ox != NULL) free(kernel_ox);
  if(x_operator != NULL) free(x_operator);
  np_glp_cv_clear_extern();
  return status;
}

static int np_conditional_x_weight_row_stream_core(double *vector_scale_factor,
                                                   int eval_idx,
                                                   double *row_out){
  return np_conditional_x_weight_row_stream_core_impl(vector_scale_factor,
                                                      eval_idx,
                                                      1,
                                                      row_out);
}

static int np_conditional_x_weight_row_full_stream_core(double *vector_scale_factor,
                                                        int eval_idx,
                                                        double *row_out){
  return np_conditional_x_weight_row_stream_core_impl(vector_scale_factor,
                                                      eval_idx,
                                                      0,
                                                      row_out);
}

int np_shadow_proof_conditional_x_weight_row_stream(double *vector_scale_factor,
                                                    int eval_idx,
                                                    double *row_out){
  return np_conditional_x_weight_row_stream_core(vector_scale_factor, eval_idx, row_out);
}

static int np_conditional_y_eval_row_stream_op_core(double *vector_scale_factor,
                                                    int eval_idx,
                                                    int operator_code,
                                                    double **matrix_Y_unordered_eval,
                                                    double **matrix_Y_ordered_eval,
                                                    double **matrix_Y_continuous_eval,
                                                    int num_eval,
                                                    int map_train_tree_index,
                                                    double *row_out){
  const int num_train = num_obs_train_extern;
  const int num_var_tot = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  const int bw_rows =
    (BANDWIDTH_den_extern == BW_FIXED) ? 1 :
    ((BANDWIDTH_den_extern == BW_GEN_NN) ? num_eval : num_train);
  int *kernel_cy = NULL, *kernel_uy = NULL, *kernel_oy = NULL, *operator_y = NULL;
  double *vsfy = NULL, *lambday = NULL, *kw = NULL;
  double **matrix_bandwidth_y = NULL, **matrix_bandwidth_eval_one = NULL;
  double **eval_yuno_one = NULL, **eval_yord_one = NULL, **eval_ycon_one = NULL;
  NPConditionalBoundState bounds_state;
  int eval_pos = eval_idx;
  int i, j, l;
  int status = 1;

  if((row_out == NULL) || (vector_scale_factor == NULL))
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;
  if((eval_idx < 0) || (eval_idx >= num_eval))
    return 1;

  memset(row_out, 0, (size_t)num_train*sizeof(double));

  if(map_train_tree_index && (int_TREE_Y == NP_TREE_TRUE) && (ipt_lookup_extern_Y != NULL))
    eval_pos = ipt_lookup_extern_Y[eval_idx];

  if(num_var_tot <= 0){
    for(j = 0; j < num_train; j++)
      row_out[j] = 1.0;
    return 0;
  }

  vsfy = alloc_vecd(MAX(1, num_var_tot));
  lambday = alloc_vecd(MAX(1, num_var_unordered_extern + num_var_ordered_extern));
  kw = alloc_vecd(MAX(1, num_train));
  matrix_bandwidth_y = alloc_tmatd(bw_rows, num_var_continuous_extern);
  matrix_bandwidth_eval_one = alloc_tmatd(1, num_var_continuous_extern);
  if(num_var_unordered_extern > 0) eval_yuno_one = alloc_matd(1, num_var_unordered_extern);
  if(num_var_ordered_extern > 0) eval_yord_one = alloc_matd(1, num_var_ordered_extern);
  if(num_var_continuous_extern > 0) eval_ycon_one = alloc_matd(1, num_var_continuous_extern);

  kernel_cy = (int *)calloc((size_t)MAX(1, num_var_continuous_extern), sizeof(int));
  kernel_uy = (int *)calloc((size_t)MAX(1, num_var_unordered_extern), sizeof(int));
  kernel_oy = (int *)calloc((size_t)MAX(1, num_var_ordered_extern), sizeof(int));
  operator_y = (int *)calloc((size_t)MAX(1, num_var_tot), sizeof(int));

  if((vsfy == NULL) || (lambday == NULL) || (kw == NULL) ||
     ((num_var_continuous_extern > 0) && (matrix_bandwidth_y == NULL)) ||
     ((num_var_continuous_extern > 0) && (matrix_bandwidth_eval_one == NULL)) ||
     ((num_var_unordered_extern > 0) && (eval_yuno_one == NULL)) ||
     ((num_var_ordered_extern > 0) && (eval_yord_one == NULL)) ||
     ((num_var_continuous_extern > 0) && (eval_ycon_one == NULL)) ||
     (kernel_cy == NULL) || (kernel_uy == NULL) || (kernel_oy == NULL) || (operator_y == NULL))
    goto cleanup_yweight_row;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        NULL,
                        vsfy,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  for(i = 0; i < num_var_continuous_extern; i++) kernel_cy[i] = KERNEL_den_extern;
  for(i = 0; i < num_var_unordered_extern; i++) kernel_uy[i] = KERNEL_den_unordered_extern;
  for(i = 0; i < num_var_ordered_extern; i++) kernel_oy[i] = KERNEL_den_ordered_extern;
  for(i = 0; i < num_var_tot; i++) operator_y[i] = operator_code;

  if(kernel_bandwidth_mean(KERNEL_den_extern,
                           BANDWIDTH_den_extern,
                           num_train,
                           num_eval,
                           0,
                           0,
                           0,
                           num_var_continuous_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           0,
                           vsfy,
                           NULL,
                           NULL,
                           matrix_Y_continuous_train_extern,
                           matrix_Y_continuous_eval,
                           NULL,
                           matrix_bandwidth_y,
                           lambday) == 1)
    goto cleanup_yweight_row;

  for(l = 0; l < num_var_unordered_extern; l++)
    eval_yuno_one[l][0] = matrix_Y_unordered_eval[l][eval_pos];
  for(l = 0; l < num_var_ordered_extern; l++)
    eval_yord_one[l][0] = matrix_Y_ordered_eval[l][eval_pos];
  for(l = 0; l < num_var_continuous_extern; l++){
    eval_ycon_one[l][0] = matrix_Y_continuous_eval[l][eval_pos];
    matrix_bandwidth_eval_one[l][0] =
      (BANDWIDTH_den_extern == BW_GEN_NN) ? matrix_bandwidth_y[l][eval_pos] : matrix_bandwidth_y[l][0];
  }

  np_conditional_push_bounds(int_cyker_bound_extern,
                             vector_cykerlb_extern,
                             vector_cykerub_extern,
                             &bounds_state);
  if(np_shadow_conditional_kernel_row(kernel_cy,
                                      kernel_uy,
                                      kernel_oy,
                                      operator_y,
                                      BANDWIDTH_den_extern,
                                      num_train,
                                      num_var_unordered_extern,
                                      num_var_ordered_extern,
                                      num_var_continuous_extern,
                                      matrix_Y_unordered_train_extern,
                                      matrix_Y_ordered_train_extern,
                                      matrix_Y_continuous_train_extern,
                                      eval_yuno_one,
                                      eval_yord_one,
                                      eval_ycon_one,
                                      vsfy,
                                      1,
                                      matrix_bandwidth_y,
                                      matrix_bandwidth_eval_one,
                                      lambday,
                                      num_categories_extern_Y,
                                      matrix_categorical_vals_extern_Y,
                                      int_TREE_Y,
                                      kdt_extern_Y,
                                      kw,
                                      NULL) != 0){
    np_conditional_pop_bounds(&bounds_state);
    goto cleanup_yweight_row;
  }
  np_conditional_pop_bounds(&bounds_state);

  for(j = 0; j < num_train; j++){
    const int orig_j = (int_TREE_Y == NP_TREE_TRUE) ? ipt_extern_Y[j] : j;
    row_out[orig_j] = kw[j];
  }

  status = 0;

cleanup_yweight_row:
  if(vsfy != NULL) free(vsfy);
  if(lambday != NULL) free(lambday);
  if(kw != NULL) free(kw);
  if(matrix_bandwidth_y != NULL) free_tmat(matrix_bandwidth_y);
  if(matrix_bandwidth_eval_one != NULL) free_tmat(matrix_bandwidth_eval_one);
  if(eval_yuno_one != NULL) free_mat(eval_yuno_one, num_var_unordered_extern);
  if(eval_yord_one != NULL) free_mat(eval_yord_one, num_var_ordered_extern);
  if(eval_ycon_one != NULL) free_mat(eval_ycon_one, num_var_continuous_extern);
  if(kernel_cy != NULL) free(kernel_cy);
  if(kernel_uy != NULL) free(kernel_uy);
  if(kernel_oy != NULL) free(kernel_oy);
  if(operator_y != NULL) free(operator_y);
  return status;
}

static int np_conditional_y_row_stream_op_core(double *vector_scale_factor,
                                               int eval_idx,
                                               int operator_code,
                                               double *row_out){
  return np_conditional_y_eval_row_stream_op_core(vector_scale_factor,
                                                  eval_idx,
                                                  operator_code,
                                                  matrix_Y_unordered_train_extern,
                                                  matrix_Y_ordered_train_extern,
                                                  matrix_Y_continuous_train_extern,
                                                  num_obs_train_extern,
                                                  1,
                                                  row_out);
}

static int np_conditional_y_row_stream_core(double *vector_scale_factor,
                                            int eval_idx,
                                            double *row_out){
  return np_conditional_y_row_stream_op_core(vector_scale_factor,
                                             eval_idx,
                                             OP_NORMAL,
                                             row_out);
}

static int np_conditional_lp_cvls_block_size(void){
  return 64;
}

static int np_conditional_x_weight_block_stream_core_impl(double *vector_scale_factor,
                                                          int eval_start,
                                                          int block_rows,
                                                          int drop_eval_self,
                                                          double **rows_out){
  const int num_train = num_obs_train_extern;
  const int num_reg_tot = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
  const int ll_mode = (int_ll_extern == LL_LP) ? LL_LP : LL_LC;
  const int bw_rows = (BANDWIDTH_den_extern == BW_FIXED) ? 1 : block_rows;
  int *kernel_cx = NULL, *kernel_ux = NULL, *kernel_ox = NULL, *x_operator = NULL;
  double *vsfx = NULL, *lambdax = NULL, *kw = NULL, *mean_row = NULL;
  double **matrix_bandwidth_x = NULL, **matrix_bandwidth_eval_one = NULL;
  double **eval_xuno_one = NULL, **eval_xord_one = NULL, **eval_xcon_one = NULL;
  double **matrix_X_continuous_eval_block = NULL;
  MATRIX KWM = NULL, RHS = NULL, SOL = NULL;
  NPConditionalBoundState bounds_state;
  int i, j, l;
  int status = 1;

  if((rows_out == NULL) || (vector_scale_factor == NULL))
    return 1;
  if(int_TREE_X == NP_TREE_TRUE)
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN))
    return 1;
  if((eval_start < 0) || (block_rows <= 0) || ((eval_start + block_rows) > num_train))
    return 1;

  if(num_reg_tot <= 0){
    const double denom = drop_eval_self ? ((double)(num_train - 1)) : ((double)num_train);
    const double w = (denom > 0.0) ? 1.0/denom : 0.0;
    for(i = 0; i < block_rows; i++){
      const int eval_pos = eval_start + i;
      memset(rows_out[i], 0, (size_t)num_train*sizeof(double));
      for(j = 0; j < num_train; j++)
        if((!drop_eval_self) || (j != eval_pos))
          rows_out[i][j] = w;
    }
    return 0;
  }

  vsfx = alloc_vecd(MAX(1, num_reg_tot));
  lambdax = alloc_vecd(MAX(1, num_reg_unordered_extern + num_reg_ordered_extern));
  kw = alloc_vecd(MAX(1, num_train));
  mean_row = alloc_vecd(MAX(1, num_train));
  matrix_bandwidth_x = alloc_tmatd(bw_rows, num_reg_continuous_extern);
  matrix_bandwidth_eval_one = alloc_tmatd(1, num_reg_continuous_extern);
  if(num_reg_unordered_extern > 0) eval_xuno_one = alloc_matd(1, num_reg_unordered_extern);
  if(num_reg_ordered_extern > 0) eval_xord_one = alloc_matd(1, num_reg_ordered_extern);
  if(num_reg_continuous_extern > 0) eval_xcon_one = alloc_matd(1, num_reg_continuous_extern);
  if(num_reg_continuous_extern > 0)
    matrix_X_continuous_eval_block = (double **)calloc((size_t)num_reg_continuous_extern, sizeof(double *));

  kernel_cx = (int *)calloc((size_t)MAX(1, num_reg_continuous_extern), sizeof(int));
  kernel_ux = (int *)calloc((size_t)MAX(1, num_reg_unordered_extern), sizeof(int));
  kernel_ox = (int *)calloc((size_t)MAX(1, num_reg_ordered_extern), sizeof(int));
  x_operator = (int *)calloc((size_t)MAX(1, num_reg_tot), sizeof(int));

  if((vsfx == NULL) || (lambdax == NULL) || (kw == NULL) || (mean_row == NULL) ||
     ((num_reg_continuous_extern > 0) && (matrix_bandwidth_x == NULL)) ||
     ((num_reg_continuous_extern > 0) && (matrix_bandwidth_eval_one == NULL)) ||
     ((num_reg_unordered_extern > 0) && (eval_xuno_one == NULL)) ||
     ((num_reg_ordered_extern > 0) && (eval_xord_one == NULL)) ||
     ((num_reg_continuous_extern > 0) && (eval_xcon_one == NULL)) ||
     ((num_reg_continuous_extern > 0) && (matrix_X_continuous_eval_block == NULL)) ||
     (kernel_cx == NULL) || (kernel_ux == NULL) || (kernel_ox == NULL) || (x_operator == NULL))
    goto cleanup_xweight_block;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        vsfx,
                        NULL,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  for(i = 0; i < num_reg_continuous_extern; i++){
    kernel_cx[i] = KERNEL_reg_extern;
    matrix_X_continuous_eval_block[i] = matrix_X_continuous_train_extern[i] + eval_start;
  }
  for(i = 0; i < num_reg_unordered_extern; i++) kernel_ux[i] = KERNEL_reg_unordered_extern;
  for(i = 0; i < num_reg_ordered_extern; i++) kernel_ox[i] = KERNEL_reg_ordered_extern;
  for(i = 0; i < num_reg_tot; i++) x_operator[i] = OP_NORMAL;

  if(kernel_bandwidth_mean(KERNEL_reg_extern,
                           BANDWIDTH_den_extern,
                           num_train,
                           block_rows,
                           0,
                           0,
                           0,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           vsfx,
                           NULL,
                           NULL,
                           matrix_X_continuous_train_extern,
                           matrix_X_continuous_eval_block,
                           NULL,
                           matrix_bandwidth_x,
                           lambdax) == 1)
    goto cleanup_xweight_block;

  if(ll_mode == LL_LP){
    const int use_bernstein = (int_glp_bernstein_extern != 0);

    if((vector_glp_degree_extern == NULL) || (num_reg_continuous_extern <= 0))
      goto cleanup_xweight_block;
    if(!np_glp_cv_cache.ready ||
       (np_glp_cv_cache.use_bernstein != use_bernstein) ||
       (np_glp_cv_cache.basis_mode != int_glp_basis_extern) ||
       (np_glp_cv_cache.num_obs != num_train) ||
       (np_glp_cv_cache.ncon != num_reg_continuous_extern) ||
       (np_glp_cv_cache.matrix_X_continuous_train_ptr != matrix_X_continuous_train_extern)){
      if(!np_glp_cv_cache_prepare(LL_LP,
                                  num_train,
                                  num_reg_continuous_extern,
                                  matrix_X_continuous_train_extern))
        goto cleanup_xweight_block;
    }
    if((np_glp_cv_cache.nterms <= 0) || (np_glp_cv_cache.basis == NULL))
      goto cleanup_xweight_block;

    KWM = mat_creat(np_glp_cv_cache.nterms, np_glp_cv_cache.nterms, UNDEFINED);
    RHS = mat_creat(np_glp_cv_cache.nterms, 1, UNDEFINED);
    SOL = mat_creat(np_glp_cv_cache.nterms, 1, UNDEFINED);
    if((KWM == NULL) || (RHS == NULL) || (SOL == NULL))
      goto cleanup_xweight_block;
  }

  for(i = 0; i < block_rows; i++){
    const int eval_pos = eval_start + i;

    memset(rows_out[i], 0, (size_t)num_train*sizeof(double));
    for(l = 0; l < num_reg_unordered_extern; l++)
      eval_xuno_one[l][0] = matrix_X_unordered_train_extern[l][eval_pos];
    for(l = 0; l < num_reg_ordered_extern; l++)
      eval_xord_one[l][0] = matrix_X_ordered_train_extern[l][eval_pos];
    for(l = 0; l < num_reg_continuous_extern; l++){
      eval_xcon_one[l][0] = matrix_X_continuous_train_extern[l][eval_pos];
      matrix_bandwidth_eval_one[l][0] =
        (BANDWIDTH_den_extern == BW_GEN_NN) ? matrix_bandwidth_x[l][i] : matrix_bandwidth_x[l][0];
    }

    np_conditional_push_bounds(int_cxker_bound_extern,
                               vector_cxkerlb_extern,
                               vector_cxkerub_extern,
                               &bounds_state);
    if(np_shadow_conditional_kernel_row_raw(kernel_cx,
                                            kernel_ux,
                                            kernel_ox,
                                            x_operator,
                                            BANDWIDTH_den_extern,
                                            num_train,
                                            num_reg_unordered_extern,
                                            num_reg_ordered_extern,
                                            num_reg_continuous_extern,
                                            matrix_X_unordered_train_extern,
                                            matrix_X_ordered_train_extern,
                                            matrix_X_continuous_train_extern,
                                            eval_xuno_one,
                                            eval_xord_one,
                                            eval_xcon_one,
                                            vsfx,
                                            1,
                                            matrix_bandwidth_x,
                                            matrix_bandwidth_eval_one,
                                            lambdax,
                                            num_categories_extern_X,
                                            matrix_categorical_vals_extern_X,
                                            int_TREE_X,
                                            kdt_extern_X,
                                            kw,
                                            mean_row) != 0){
      np_conditional_pop_bounds(&bounds_state);
      goto cleanup_xweight_block;
    }
    np_conditional_pop_bounds(&bounds_state);

    if(drop_eval_self)
      kw[eval_pos] = 0.0;

    if(ll_mode == LL_LC){
      double row_sum = 0.0;
      for(j = 0; j < num_train; j++)
        row_sum += kw[j];
      if(!(row_sum > DBL_MIN))
        goto cleanup_xweight_block;
      for(j = 0; j < num_train; j++)
        rows_out[i][j] = kw[j]/row_sum;
    } else {
      const int k = np_glp_cv_cache.nterms;

      if(drop_eval_self){
        if(np_glp_qr_drop_row_bkcde(np_glp_cv_cache.basis,
                                    num_train,
                                    k,
                                    kw,
                                    eval_pos,
                                    rows_out[i]) != 0)
          goto cleanup_xweight_block;
      } else {
        for(l = 0; l < k; l++){
          RHS[l][0] = np_glp_cv_cache.basis[l][eval_pos];
          for(j = 0; j < k; j++)
            KWM[l][j] = 0.0;
        }

        for(j = 0; j < num_train; j++){
          const double wj = kw[j];
          if(wj == 0.0)
            continue;
          for(int a = 0; a < k; a++){
            const double za = np_glp_cv_cache.basis[a][j];
            for(int b = a; b < k; b++){
              const double zb = np_glp_cv_cache.basis[b][j];
              KWM[a][b] += wj*za*zb;
              if(b != a) KWM[b][a] += wj*za*zb;
            }
          }
        }

        if(np_mat_bad_rcond_sym(KWM, 1.0e-10))
          goto cleanup_xweight_block;
        if(mat_solve(KWM, RHS, SOL) == NULL)
          goto cleanup_xweight_block;

        for(j = 0; j < num_train; j++){
          double zju = 0.0;
          for(l = 0; l < k; l++)
            zju += np_glp_cv_cache.basis[l][j]*SOL[l][0];
          rows_out[i][j] = kw[j]*zju;
        }
      }
    }
  }

  status = 0;

cleanup_xweight_block:
  if(KWM != NULL) mat_free(KWM);
  if(RHS != NULL) mat_free(RHS);
  if(SOL != NULL) mat_free(SOL);
  if(vsfx != NULL) free(vsfx);
  if(lambdax != NULL) free(lambdax);
  if(kw != NULL) free(kw);
  if(mean_row != NULL) free(mean_row);
  if(matrix_bandwidth_x != NULL) free_tmat(matrix_bandwidth_x);
  if(matrix_bandwidth_eval_one != NULL) free_tmat(matrix_bandwidth_eval_one);
  if(eval_xuno_one != NULL) free_mat(eval_xuno_one, num_reg_unordered_extern);
  if(eval_xord_one != NULL) free_mat(eval_xord_one, num_reg_ordered_extern);
  if(eval_xcon_one != NULL) free_mat(eval_xcon_one, num_reg_continuous_extern);
  if(matrix_X_continuous_eval_block != NULL) free(matrix_X_continuous_eval_block);
  if(kernel_cx != NULL) free(kernel_cx);
  if(kernel_ux != NULL) free(kernel_ux);
  if(kernel_ox != NULL) free(kernel_ox);
  if(x_operator != NULL) free(x_operator);
  return status;
}

static int np_conditional_x_weight_block_stream_core(double *vector_scale_factor,
                                                     int eval_start,
                                                     int block_rows,
                                                     double **rows_out){
  return np_conditional_x_weight_block_stream_core_impl(vector_scale_factor,
                                                        eval_start,
                                                        block_rows,
                                                        1,
                                                        rows_out);
}

static int np_conditional_x_weight_block_full_stream_core(double *vector_scale_factor,
                                                          int eval_start,
                                                          int block_rows,
                                                          double **rows_out){
  return np_conditional_x_weight_block_stream_core_impl(vector_scale_factor,
                                                        eval_start,
                                                        block_rows,
                                                        0,
                                                        rows_out);
}

static int np_conditional_y_block_stream_op_core(double *vector_scale_factor,
                                                 int eval_start,
                                                 int block_rows,
                                                 int operator_code,
                                                 double **rows_out){
  const int num_train = num_obs_train_extern;
  const int num_var_tot = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  const int bw_rows = (BANDWIDTH_den_extern == BW_FIXED) ? 1 : block_rows;
  int *kernel_cy = NULL, *kernel_uy = NULL, *kernel_oy = NULL, *operator_y = NULL;
  double *vsfy = NULL, *lambday = NULL, *kw = NULL;
  double **matrix_bandwidth_y = NULL, **matrix_bandwidth_eval_one = NULL;
  double **eval_yuno_one = NULL, **eval_yord_one = NULL, **eval_ycon_one = NULL;
  double **matrix_Y_continuous_eval_block = NULL;
  NPConditionalBoundState bounds_state;
  int i, j, l;
  int status = 1;

  if((rows_out == NULL) || (vector_scale_factor == NULL))
    return 1;
  if(int_TREE_Y == NP_TREE_TRUE)
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN))
    return 1;
  if((eval_start < 0) || (block_rows <= 0) || ((eval_start + block_rows) > num_train))
    return 1;

  if(num_var_tot <= 0){
    for(i = 0; i < block_rows; i++)
      for(j = 0; j < num_train; j++)
        rows_out[i][j] = 1.0;
    return 0;
  }

  vsfy = alloc_vecd(MAX(1, num_var_tot));
  lambday = alloc_vecd(MAX(1, num_var_unordered_extern + num_var_ordered_extern));
  kw = alloc_vecd(MAX(1, num_train));
  matrix_bandwidth_y = alloc_tmatd(bw_rows, num_var_continuous_extern);
  matrix_bandwidth_eval_one = alloc_tmatd(1, num_var_continuous_extern);
  if(num_var_unordered_extern > 0) eval_yuno_one = alloc_matd(1, num_var_unordered_extern);
  if(num_var_ordered_extern > 0) eval_yord_one = alloc_matd(1, num_var_ordered_extern);
  if(num_var_continuous_extern > 0) eval_ycon_one = alloc_matd(1, num_var_continuous_extern);
  if(num_var_continuous_extern > 0)
    matrix_Y_continuous_eval_block = (double **)calloc((size_t)num_var_continuous_extern, sizeof(double *));

  kernel_cy = (int *)calloc((size_t)MAX(1, num_var_continuous_extern), sizeof(int));
  kernel_uy = (int *)calloc((size_t)MAX(1, num_var_unordered_extern), sizeof(int));
  kernel_oy = (int *)calloc((size_t)MAX(1, num_var_ordered_extern), sizeof(int));
  operator_y = (int *)calloc((size_t)MAX(1, num_var_tot), sizeof(int));

  if((vsfy == NULL) || (lambday == NULL) || (kw == NULL) ||
     ((num_var_continuous_extern > 0) && (matrix_bandwidth_y == NULL)) ||
     ((num_var_continuous_extern > 0) && (matrix_bandwidth_eval_one == NULL)) ||
     ((num_var_unordered_extern > 0) && (eval_yuno_one == NULL)) ||
     ((num_var_ordered_extern > 0) && (eval_yord_one == NULL)) ||
     ((num_var_continuous_extern > 0) && (eval_ycon_one == NULL)) ||
     ((num_var_continuous_extern > 0) && (matrix_Y_continuous_eval_block == NULL)) ||
     (kernel_cy == NULL) || (kernel_uy == NULL) || (kernel_oy == NULL) || (operator_y == NULL))
    goto cleanup_yweight_block;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        NULL,
                        vsfy,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  for(i = 0; i < num_var_continuous_extern; i++){
    kernel_cy[i] = KERNEL_den_extern;
    matrix_Y_continuous_eval_block[i] = matrix_Y_continuous_train_extern[i] + eval_start;
  }
  for(i = 0; i < num_var_unordered_extern; i++) kernel_uy[i] = KERNEL_den_unordered_extern;
  for(i = 0; i < num_var_ordered_extern; i++) kernel_oy[i] = KERNEL_den_ordered_extern;
  for(i = 0; i < num_var_tot; i++) operator_y[i] = operator_code;

  if(kernel_bandwidth_mean(KERNEL_den_extern,
                           BANDWIDTH_den_extern,
                           num_train,
                           block_rows,
                           0,
                           0,
                           0,
                           num_var_continuous_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           0,
                           vsfy,
                           NULL,
                           NULL,
                           matrix_Y_continuous_train_extern,
                           matrix_Y_continuous_eval_block,
                           NULL,
                           matrix_bandwidth_y,
                           lambday) == 1)
    goto cleanup_yweight_block;

  for(i = 0; i < block_rows; i++){
    const int eval_pos = eval_start + i;

    memset(rows_out[i], 0, (size_t)num_train*sizeof(double));
    for(l = 0; l < num_var_unordered_extern; l++)
      eval_yuno_one[l][0] = matrix_Y_unordered_train_extern[l][eval_pos];
    for(l = 0; l < num_var_ordered_extern; l++)
      eval_yord_one[l][0] = matrix_Y_ordered_train_extern[l][eval_pos];
    for(l = 0; l < num_var_continuous_extern; l++){
      eval_ycon_one[l][0] = matrix_Y_continuous_train_extern[l][eval_pos];
      matrix_bandwidth_eval_one[l][0] =
        (BANDWIDTH_den_extern == BW_GEN_NN) ? matrix_bandwidth_y[l][i] : matrix_bandwidth_y[l][0];
    }

    np_conditional_push_bounds(int_cyker_bound_extern,
                               vector_cykerlb_extern,
                               vector_cykerub_extern,
                               &bounds_state);
    if(np_shadow_conditional_kernel_row(kernel_cy,
                                        kernel_uy,
                                        kernel_oy,
                                        operator_y,
                                        BANDWIDTH_den_extern,
                                        num_train,
                                        num_var_unordered_extern,
                                        num_var_ordered_extern,
                                        num_var_continuous_extern,
                                        matrix_Y_unordered_train_extern,
                                        matrix_Y_ordered_train_extern,
                                        matrix_Y_continuous_train_extern,
                                        eval_yuno_one,
                                        eval_yord_one,
                                        eval_ycon_one,
                                        vsfy,
                                        1,
                                        matrix_bandwidth_y,
                                        matrix_bandwidth_eval_one,
                                        lambday,
                                        num_categories_extern_Y,
                                        matrix_categorical_vals_extern_Y,
                                        int_TREE_Y,
                                        kdt_extern_Y,
                                        kw,
                                        NULL) != 0){
      np_conditional_pop_bounds(&bounds_state);
      goto cleanup_yweight_block;
    }
    np_conditional_pop_bounds(&bounds_state);

    for(j = 0; j < num_train; j++)
      rows_out[i][j] = kw[j];
  }

  status = 0;

cleanup_yweight_block:
  if(vsfy != NULL) free(vsfy);
  if(lambday != NULL) free(lambday);
  if(kw != NULL) free(kw);
  if(matrix_bandwidth_y != NULL) free_tmat(matrix_bandwidth_y);
  if(matrix_bandwidth_eval_one != NULL) free_tmat(matrix_bandwidth_eval_one);
  if(eval_yuno_one != NULL) free_mat(eval_yuno_one, num_var_unordered_extern);
  if(eval_yord_one != NULL) free_mat(eval_yord_one, num_var_ordered_extern);
  if(eval_ycon_one != NULL) free_mat(eval_ycon_one, num_var_continuous_extern);
  if(matrix_Y_continuous_eval_block != NULL) free(matrix_Y_continuous_eval_block);
  if(kernel_cy != NULL) free(kernel_cy);
  if(kernel_uy != NULL) free(kernel_uy);
  if(kernel_oy != NULL) free(kernel_oy);
  if(operator_y != NULL) free(operator_y);
  return status;
}

static int np_conditional_y_eval_block_stream_op_core(double *vector_scale_factor,
                                                      int eval_start,
                                                      int block_rows,
                                                      int operator_code,
                                                      double **matrix_Y_unordered_eval,
                                                      double **matrix_Y_ordered_eval,
                                                      double **matrix_Y_continuous_eval,
                                                      int num_eval,
                                                      double **rows_out){
  const int num_train = num_obs_train_extern;
  const int num_var_tot = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  const int bw_rows = (BANDWIDTH_den_extern == BW_FIXED) ? 1 : block_rows;
  int *kernel_cy = NULL, *kernel_uy = NULL, *kernel_oy = NULL, *operator_y = NULL;
  double *vsfy = NULL, *lambday = NULL, *kw = NULL;
  double **matrix_bandwidth_y = NULL, **matrix_bandwidth_eval_one = NULL;
  double **eval_yuno_one = NULL, **eval_yord_one = NULL, **eval_ycon_one = NULL;
  double **matrix_Y_continuous_eval_block = NULL;
  NPConditionalBoundState bounds_state;
  int i, j, l;
  int status = 1;

  if((rows_out == NULL) || (vector_scale_factor == NULL))
    return 1;
  if(int_TREE_Y == NP_TREE_TRUE)
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN))
    return 1;
  if((eval_start < 0) || (block_rows <= 0) || ((eval_start + block_rows) > num_eval))
    return 1;

  if(num_var_tot <= 0){
    for(i = 0; i < block_rows; i++)
      for(j = 0; j < num_train; j++)
        rows_out[i][j] = 1.0;
    return 0;
  }

  vsfy = alloc_vecd(MAX(1, num_var_tot));
  lambday = alloc_vecd(MAX(1, num_var_unordered_extern + num_var_ordered_extern));
  kw = alloc_vecd(MAX(1, num_train));
  matrix_bandwidth_y = alloc_tmatd(bw_rows, num_var_continuous_extern);
  matrix_bandwidth_eval_one = alloc_tmatd(1, num_var_continuous_extern);
  if(num_var_unordered_extern > 0) eval_yuno_one = alloc_matd(1, num_var_unordered_extern);
  if(num_var_ordered_extern > 0) eval_yord_one = alloc_matd(1, num_var_ordered_extern);
  if(num_var_continuous_extern > 0) eval_ycon_one = alloc_matd(1, num_var_continuous_extern);
  if(num_var_continuous_extern > 0)
    matrix_Y_continuous_eval_block = (double **)calloc((size_t)num_var_continuous_extern, sizeof(double *));

  kernel_cy = (int *)calloc((size_t)MAX(1, num_var_continuous_extern), sizeof(int));
  kernel_uy = (int *)calloc((size_t)MAX(1, num_var_unordered_extern), sizeof(int));
  kernel_oy = (int *)calloc((size_t)MAX(1, num_var_ordered_extern), sizeof(int));
  operator_y = (int *)calloc((size_t)MAX(1, num_var_tot), sizeof(int));

  if((vsfy == NULL) || (lambday == NULL) || (kw == NULL) ||
     ((num_var_continuous_extern > 0) && (matrix_bandwidth_y == NULL)) ||
     ((num_var_continuous_extern > 0) && (matrix_bandwidth_eval_one == NULL)) ||
     ((num_var_unordered_extern > 0) && (eval_yuno_one == NULL)) ||
     ((num_var_ordered_extern > 0) && (eval_yord_one == NULL)) ||
     ((num_var_continuous_extern > 0) && (eval_ycon_one == NULL)) ||
     ((num_var_continuous_extern > 0) && (matrix_Y_continuous_eval_block == NULL)) ||
     (kernel_cy == NULL) || (kernel_uy == NULL) || (kernel_oy == NULL) || (operator_y == NULL))
    goto cleanup_yweight_eval_block;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        NULL,
                        vsfy,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  for(i = 0; i < num_var_continuous_extern; i++){
    kernel_cy[i] = KERNEL_den_extern;
    matrix_Y_continuous_eval_block[i] = matrix_Y_continuous_eval[i] + eval_start;
  }
  for(i = 0; i < num_var_unordered_extern; i++) kernel_uy[i] = KERNEL_den_unordered_extern;
  for(i = 0; i < num_var_ordered_extern; i++) kernel_oy[i] = KERNEL_den_ordered_extern;
  for(i = 0; i < num_var_tot; i++) operator_y[i] = operator_code;

  if(kernel_bandwidth_mean(KERNEL_den_extern,
                           BANDWIDTH_den_extern,
                           num_train,
                           block_rows,
                           0,
                           0,
                           0,
                           num_var_continuous_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           0,
                           vsfy,
                           NULL,
                           NULL,
                           matrix_Y_continuous_train_extern,
                           matrix_Y_continuous_eval_block,
                           NULL,
                           matrix_bandwidth_y,
                           lambday) == 1)
    goto cleanup_yweight_eval_block;

  for(i = 0; i < block_rows; i++){
    const int eval_pos = eval_start + i;

    memset(rows_out[i], 0, (size_t)num_train*sizeof(double));
    for(l = 0; l < num_var_unordered_extern; l++)
      eval_yuno_one[l][0] = matrix_Y_unordered_eval[l][eval_pos];
    for(l = 0; l < num_var_ordered_extern; l++)
      eval_yord_one[l][0] = matrix_Y_ordered_eval[l][eval_pos];
    for(l = 0; l < num_var_continuous_extern; l++){
      eval_ycon_one[l][0] = matrix_Y_continuous_eval[l][eval_pos];
      matrix_bandwidth_eval_one[l][0] =
        (BANDWIDTH_den_extern == BW_GEN_NN) ? matrix_bandwidth_y[l][i] : matrix_bandwidth_y[l][0];
    }

    np_conditional_push_bounds(int_cyker_bound_extern,
                               vector_cykerlb_extern,
                               vector_cykerub_extern,
                               &bounds_state);
    if(np_shadow_conditional_kernel_row(kernel_cy,
                                        kernel_uy,
                                        kernel_oy,
                                        operator_y,
                                        BANDWIDTH_den_extern,
                                        num_train,
                                        num_var_unordered_extern,
                                        num_var_ordered_extern,
                                        num_var_continuous_extern,
                                        matrix_Y_unordered_train_extern,
                                        matrix_Y_ordered_train_extern,
                                        matrix_Y_continuous_train_extern,
                                        eval_yuno_one,
                                        eval_yord_one,
                                        eval_ycon_one,
                                        vsfy,
                                        1,
                                        matrix_bandwidth_y,
                                        matrix_bandwidth_eval_one,
                                        lambday,
                                        num_categories_extern_Y,
                                        matrix_categorical_vals_extern_Y,
                                        int_TREE_Y,
                                        kdt_extern_Y,
                                        kw,
                                        NULL) != 0){
      np_conditional_pop_bounds(&bounds_state);
      goto cleanup_yweight_eval_block;
    }
    np_conditional_pop_bounds(&bounds_state);

    for(j = 0; j < num_train; j++)
      rows_out[i][j] = kw[j];
  }

  status = 0;

cleanup_yweight_eval_block:
  if(vsfy != NULL) free(vsfy);
  if(lambday != NULL) free(lambday);
  if(kw != NULL) free(kw);
  if(matrix_bandwidth_y != NULL) free_tmat(matrix_bandwidth_y);
  if(matrix_bandwidth_eval_one != NULL) free_tmat(matrix_bandwidth_eval_one);
  if(eval_yuno_one != NULL) free_mat(eval_yuno_one, num_var_unordered_extern);
  if(eval_yord_one != NULL) free_mat(eval_yord_one, num_var_ordered_extern);
  if(eval_ycon_one != NULL) free_mat(eval_ycon_one, num_var_continuous_extern);
  if(matrix_Y_continuous_eval_block != NULL) free(matrix_Y_continuous_eval_block);
  if(kernel_cy != NULL) free(kernel_cy);
  if(kernel_uy != NULL) free(kernel_uy);
  if(kernel_oy != NULL) free(kernel_oy);
  if(operator_y != NULL) free(operator_y);
  return status;
}

#define NP_BOUNDED_CVLS_I1_GRID_POINTS 81
#define NP_BOUNDED_CVLS_I1_GRID_POINTS_2D 31
#define NP_BOUNDED_CVLS_I1_MODE_BOOK 1
#define NP_BOUNDED_CVLS_I1_MODE_FULL 0
#define NP_BOUNDED_CVLS_DEFAULT_EXTEND_FACTOR 1.0
#define NP_BOUNDED_CVLS_GRID_UNIFORM 0
#define NP_BOUNDED_CVLS_GRID_HYBRID 1
#define NP_BOUNDED_CVLS_GRID_SAMPLE 2

static int np_bounded_cvls_continuous_support_ok(const int ncon,
                                                 const double *lb,
                                                 const double *ub){
  int d;

  if((ncon < 1) || (ncon > 2))
    return 0;
  if((lb == NULL) || (ub == NULL))
    return 0;

  for(d = 0; d < ncon; d++){
    if((!R_FINITE(lb[d])) || (!R_FINITE(ub[d])) || (!(ub[d] > lb[d])))
      return 0;
  }

  return 1;
}

static int np_bounded_cvls_bound_is_effectively_infinite(const double x){
  const double sentinel = 0.5*DBL_MAX;
  const double route_sentinel = 0.25*DBL_MAX;
  const double ax = fabs(x);
  return (!R_FINITE(x)) ||
    (ax >= sentinel) ||
    (fabs(ax - route_sentinel) <= route_sentinel*(8.0*DBL_EPSILON));
}

static void np_bounded_cvls_conditional_effective_integration_bounds(const int ncon,
                                                                     double **train_continuous,
                                                                     const int num_obs,
                                                                     const double *support_lb,
                                                                     const double *support_ub,
                                                                     double *quad_lb,
                                                                     double *quad_ub){
  int d, i;
  double extend_factor = double_bounded_cvls_quadrature_extend_factor_extern;

  if((!R_FINITE(extend_factor)) || (extend_factor <= 0.0))
    extend_factor = NP_BOUNDED_CVLS_DEFAULT_EXTEND_FACTOR;

  for(d = 0; d < ncon; d++){
    const double * const train_y = train_continuous[d];
    double ymin = train_y[0];
    double ymax = train_y[0];
    double span;
    double ext;

    for(i = 1; i < num_obs; i++){
      if(train_y[i] < ymin)
        ymin = train_y[i];
      if(train_y[i] > ymax)
        ymax = train_y[i];
    }

    /* Keep finite support edges intact and extend only infinite sides
       numerically. This avoids Inf/DBL_MAX quadrature-grid dilution while
       preserving the bounded-kernel evaluation path. */
    span = ymax - ymin;
    if((span > 0.0) && R_FINITE(span))
      ext = extend_factor*span;
    else {
      const double scale = MAX(1.0, MAX(fabs(ymin), fabs(ymax)));
      ext = extend_factor*scale;
    }

    quad_lb[d] = np_bounded_cvls_bound_is_effectively_infinite(support_lb[d]) ?
      (ymin - ext) : support_lb[d];
    quad_ub[d] = np_bounded_cvls_bound_is_effectively_infinite(support_ub[d]) ?
      (ymax + ext) : support_ub[d];
  }
}

static int np_bounded_cvls_grid_points(const int ncon){
  return (ncon == 1) ? NP_BOUNDED_CVLS_I1_GRID_POINTS :
    NP_BOUNDED_CVLS_I1_GRID_POINTS_2D;
}

static int np_bounded_cvls_conditional_grid_points(const int ncon){
  if(int_bounded_cvls_quadrature_points_extern >= 2)
    return int_bounded_cvls_quadrature_points_extern;

  return np_bounded_cvls_grid_points(ncon);
}

static int np_conditional_density_cvls_bounded_scalar_route_ok(void){
  if(int_cyker_bound_extern == 0)
    return 0;
  if(num_var_continuous_extern != 1)
    return 0;
  if(num_var_unordered_extern != 0)
    return 0;
  if(num_var_ordered_extern != 0)
    return 0;
  return np_bounded_cvls_continuous_support_ok(num_var_continuous_extern,
                                               vector_cykerlb_extern,
                                               vector_cykerub_extern);
}

static void np_fill_trapezoid_rule(double a,
                                   double b,
                                   int n,
                                   double *grid,
                                   double *weights){
  const double step = (b - a)/(double)(n - 1);
  int i;

  for(i = 0; i < n; i++){
    grid[i] = a + ((double)i)*step;
    weights[i] = step;
  }
  weights[0] *= 0.5;
  weights[n - 1] *= 0.5;
}

static int np_compare_double_asc(const void *lhs, const void *rhs){
  const double a = *((const double *)lhs);
  const double b = *((const double *)rhs);
  return (a > b) - (a < b);
}

static int np_fill_trapezoid_rule_nonuniform(const double *grid,
                                             int n,
                                             double *weights){
  int i;

  if((grid == NULL) || (weights == NULL) || (n < 2))
    return 1;
  if(!(grid[1] > grid[0]))
    return 1;

  weights[0] = 0.5*(grid[1] - grid[0]);
  for(i = 1; i < (n - 1); i++){
    if(!(grid[i] > grid[i - 1]) || !(grid[i + 1] > grid[i]))
      return 1;
    weights[i] = 0.5*(grid[i + 1] - grid[i - 1]);
  }
  weights[n - 1] = 0.5*(grid[n - 1] - grid[n - 2]);

  return 0;
}

static int np_bounded_cvls_unique_sorted_inplace(double *x, const int n){
  int i, unique_count;

  if((x == NULL) || (n <= 0))
    return 0;

  qsort(x, (size_t)n, sizeof(double), np_compare_double_asc);
  unique_count = 0;
  for(i = 0; i < n; i++){
    if((unique_count == 0) || (x[i] > x[unique_count - 1]))
      x[unique_count++] = x[i];
  }

  return unique_count;
}

static int np_bounded_cvls_collect_sample_nodes_1d(const double *train_y,
                                                   const int num_obs,
                                                   const double lb,
                                                   const double ub,
                                                   double *work){
  int i, count;

  if((train_y == NULL) || (work == NULL) || (num_obs <= 0) || !(ub > lb))
    return 0;

  count = 0;
  for(i = 0; i < num_obs; i++){
    const double y = train_y[i];
    if(R_FINITE(y) && (y >= lb) && (y <= ub))
      work[count++] = y;
  }

  return np_bounded_cvls_unique_sorted_inplace(work, count);
}

static int np_bounded_cvls_ranked_sample_nodes_1d(const double *unique_y,
                                                  const int n_unique,
                                                  const int q_target,
                                                  double *grid_out,
                                                  int *q_out){
  int i, q_use, unique_count;

  if((unique_y == NULL) || (grid_out == NULL) || (q_out == NULL) ||
     (n_unique < 2) || (q_target < 2))
    return 1;

  q_use = MIN(q_target, n_unique);
  for(i = 0; i < q_use; i++){
    const double pos =
      ((double)i)*((double)(n_unique - 1))/((double)(q_use - 1));
    int idx = (int)floor(pos + 0.5);
    idx = MAX(0, MIN(n_unique - 1, idx));
    grid_out[i] = unique_y[idx];
  }

  unique_count = np_bounded_cvls_unique_sorted_inplace(grid_out, q_use);
  if(unique_count < 2)
    return 1;

  *q_out = unique_count;
  return 0;
}

static int np_bounded_cvls_fill_uniform_nodes_1d(const double lb,
                                                 const double ub,
                                                 const int q,
                                                 double *grid_out){
  const double step = (ub - lb)/(double)(q - 1);
  int i;

  if((grid_out == NULL) || (q < 2) || !(ub > lb) || !R_FINITE(step))
    return 1;

  for(i = 0; i < q; i++)
    grid_out[i] = lb + ((double)i)*step;

  return 0;
}

static int np_bounded_cvls_node_exists_1d(const double *grid,
                                          const int n,
                                          const double value){
  int i;

  for(i = 0; i < n; i++){
    if(grid[i] == value)
      return 1;
  }

  return 0;
}

static void np_bounded_cvls_append_unique_node_1d(double *grid,
                                                  int *count,
                                                  const int q_target,
                                                  const double value){
  if((grid == NULL) || (count == NULL) || (*count >= q_target) ||
     !R_FINITE(value))
    return;

  if(!np_bounded_cvls_node_exists_1d(grid, *count, value))
    grid[(*count)++] = value;
}

static void np_bounded_cvls_fill_remaining_uniform_nodes_1d(const double lb,
                                                            const double ub,
                                                            const int q_target,
                                                            double *grid,
                                                            int *count){
  int pass;

  if((grid == NULL) || (count == NULL) || (q_target < 2) || !(ub > lb))
    return;

  for(pass = 1; (*count < q_target) && (pass <= 8); pass++){
    const int q_try = MAX(q_target, q_target*(pass + 1));
    const double step = (ub - lb)/(double)(q_try - 1);
    int i;

    if(!R_FINITE(step) || !(step > 0.0))
      return;

    for(i = 0; (i < q_try) && (*count < q_target); i++){
      const double value = lb + ((double)i)*step;
      np_bounded_cvls_append_unique_node_1d(grid, count, q_target, value);
    }
  }
}

static int np_bounded_cvls_largest_remainder_index(const double rem0,
                                                   const double rem1,
                                                   const double rem2){
  if((rem1 >= rem0) && (rem1 >= rem2))
    return 1;
  if((rem2 >= rem0) && (rem2 >= rem1))
    return 2;
  return 0;
}

static void np_bounded_cvls_split_hybrid_counts(const int q_target,
                                                int *q_uniform,
                                                int *q_sample,
                                                int *q_gl){
  double r0 = double_bounded_cvls_quadrature_ratios_extern[0];
  double r1 = double_bounded_cvls_quadrature_ratios_extern[1];
  double r2 = double_bounded_cvls_quadrature_ratios_extern[2];
  double raw0, raw1, raw2, rem0, rem1, rem2;
  int c0, c1, c2, total;

  if((!R_FINITE(r0)) || (!R_FINITE(r1)) || (!R_FINITE(r2)) ||
     (r0 < 0.0) || (r1 < 0.0) || (r2 < 0.0) ||
     (fabs(r0 + r1 + r2 - 1.0) > 1.0e-8)){
    r0 = 0.20;
    r1 = 0.55;
    r2 = 0.25;
  }

  raw0 = ((double)q_target)*r0;
  raw1 = ((double)q_target)*r1;
  raw2 = ((double)q_target)*r2;

  c0 = (int)floor(raw0);
  c1 = (int)floor(raw1);
  c2 = (int)floor(raw2);
  rem0 = raw0 - (double)c0;
  rem1 = raw1 - (double)c1;
  rem2 = raw2 - (double)c2;

  total = c0 + c1 + c2;
  while(total < q_target){
    const int idx = np_bounded_cvls_largest_remainder_index(rem0, rem1, rem2);
    if(idx == 1){
      c1++;
      rem1 = -1.0;
    } else if(idx == 2){
      c2++;
      rem2 = -1.0;
    } else {
      c0++;
      rem0 = -1.0;
    }
    total++;
  }

  while(total > q_target){
    if((c0 >= c1) && (c0 >= c2) && (c0 > 0))
      c0--;
    else if((c1 >= c0) && (c1 >= c2) && (c1 > 0))
      c1--;
    else if(c2 > 0)
      c2--;
    total--;
  }

  if((c2 > 0) && (c2 < 2)){
    if(c1 > 0)
      c1--;
    else if(c0 > 0)
      c0--;
    c2++;
  }
  if((c2 % 2) != 0){
    c2--;
    if(c1 >= c0)
      c1++;
    else
      c0++;
  }

  *q_uniform = c0;
  *q_sample = c1;
  *q_gl = c2;
}

static int np_bounded_cvls_build_conditional_grid_1d(const double *train_y,
                                                     const int num_obs,
                                                     const double lb,
                                                     const double ub,
                                                     const int q_target,
                                                     double *grid_out,
                                                     double *weights_out,
                                                     int *q_out){
  int grid_mode = int_bounded_cvls_quadrature_grid_extern;
  double *sample_work = NULL;
  double *sample_nodes = NULL;
  double *gl_breaks = NULL;
  int n_unique = 0;
  int status = 1;

  if((grid_out == NULL) || (weights_out == NULL) || (q_out == NULL) ||
     (q_target < 2) || !(ub > lb))
    return 1;

  if((grid_mode < NP_BOUNDED_CVLS_GRID_UNIFORM) ||
     (grid_mode > NP_BOUNDED_CVLS_GRID_SAMPLE))
    grid_mode = NP_BOUNDED_CVLS_GRID_HYBRID;

  if(grid_mode == NP_BOUNDED_CVLS_GRID_UNIFORM){
    np_fill_trapezoid_rule(lb, ub, q_target, grid_out, weights_out);
    *q_out = q_target;
    return 0;
  }

  sample_work = alloc_vecd(num_obs);
  if(sample_work == NULL)
    goto cleanup_build_conditional_grid_1d;

  n_unique = np_bounded_cvls_collect_sample_nodes_1d(train_y,
                                                     num_obs,
                                                     lb,
                                                     ub,
                                                     sample_work);

  if(grid_mode == NP_BOUNDED_CVLS_GRID_SAMPLE){
    if(np_bounded_cvls_ranked_sample_nodes_1d(sample_work,
                                             n_unique,
                                             q_target,
                                             grid_out,
                                             q_out) != 0)
      goto cleanup_build_conditional_grid_1d;
    if(np_fill_trapezoid_rule_nonuniform(grid_out, *q_out, weights_out) != 0)
      goto cleanup_build_conditional_grid_1d;
    status = 0;
    goto cleanup_build_conditional_grid_1d;
  }

  sample_nodes = alloc_vecd(q_target);
  if(sample_nodes == NULL)
    goto cleanup_build_conditional_grid_1d;

  {
    const double inv_sqrt3 = 0.57735026918962576451;
    int q_uniform_target = 0;
    int q_sample_target = 0;
    int q_gl_target = 0;
    int count = 0;
    int q_sample = 0;
    int q_breaks = 0;
    int i;

    np_bounded_cvls_split_hybrid_counts(q_target,
                                        &q_uniform_target,
                                        &q_sample_target,
                                        &q_gl_target);

    if((q_uniform_target >= 2) &&
       (np_bounded_cvls_fill_uniform_nodes_1d(lb, ub, q_uniform_target, grid_out) == 0)){
      count = q_uniform_target;
    }

    if((q_sample_target >= 2) &&
       (np_bounded_cvls_ranked_sample_nodes_1d(sample_work,
                                              n_unique,
                                              q_sample_target,
                                              sample_nodes,
                                              &q_sample) == 0)){
      for(i = 0; i < q_sample; i++)
        np_bounded_cvls_append_unique_node_1d(grid_out, &count, q_target, sample_nodes[i]);
    }

    if(q_gl_target >= 2){
      const int q_intervals = q_gl_target/2;
      gl_breaks = alloc_vecd(q_intervals + 1);
      if((gl_breaks != NULL) &&
         (np_bounded_cvls_ranked_sample_nodes_1d(sample_work,
                                                n_unique,
                                                q_intervals + 1,
                                                gl_breaks,
                                                &q_breaks) == 0)){
        for(i = 0; i < q_breaks - 1; i++){
          const double a = gl_breaks[i];
          const double b = gl_breaks[i + 1];

          if(b > a){
            const double mid = 0.5*(a + b);
            const double half = 0.5*(b - a);

            np_bounded_cvls_append_unique_node_1d(grid_out, &count, q_target,
                                                  mid - inv_sqrt3*half);
            np_bounded_cvls_append_unique_node_1d(grid_out, &count, q_target,
                                                  mid + inv_sqrt3*half);
          }
        }
      }
    }

    count = np_bounded_cvls_unique_sorted_inplace(grid_out, count);
    np_bounded_cvls_fill_remaining_uniform_nodes_1d(lb, ub, q_target, grid_out, &count);
    count = np_bounded_cvls_unique_sorted_inplace(grid_out, count);

    if(count == q_target){
      *q_out = count;
      if(np_fill_trapezoid_rule_nonuniform(grid_out, *q_out, weights_out) != 0)
        goto cleanup_build_conditional_grid_1d;
      status = 0;
      goto cleanup_build_conditional_grid_1d;
    }
  }

  {
    const int q_uniform_start = MAX(2, (q_target + 1)/2);
    int q_uniform;

    for(q_uniform = q_uniform_start; q_uniform <= q_target; q_uniform++){
      const int q_sample_target = q_target - q_uniform;
      int count = 0;
      int q_sample = 0;
      int i;

      if(np_bounded_cvls_fill_uniform_nodes_1d(lb, ub, q_uniform, grid_out) != 0)
        goto cleanup_build_conditional_grid_1d;
      count = q_uniform;

      if((q_sample_target >= 2) &&
         (np_bounded_cvls_ranked_sample_nodes_1d(sample_work,
                                                n_unique,
                                                q_sample_target,
                                                sample_nodes,
                                                &q_sample) == 0)){
        for(i = 0; i < q_sample; i++)
          grid_out[count++] = sample_nodes[i];
      }

      count = np_bounded_cvls_unique_sorted_inplace(grid_out, count);
      if((count >= 2) && ((count >= q_target) || (q_uniform == q_target))){
        *q_out = count;
        if(np_fill_trapezoid_rule_nonuniform(grid_out, *q_out, weights_out) != 0)
          goto cleanup_build_conditional_grid_1d;
        status = 0;
        goto cleanup_build_conditional_grid_1d;
      }
    }
  }

cleanup_build_conditional_grid_1d:
  if(sample_work != NULL) free(sample_work);
  if(sample_nodes != NULL) free(sample_nodes);
  if(gl_breaks != NULL) free(gl_breaks);
  return status;
}

typedef struct {
  int ready;
  int required;
  int q;
  int num_obs;
  double quad_lb[2];
  double quad_ub[2];
  double *base_grid;
  double *base_weights;
} NPBoundedCVLSConditionalQuadContext;

static NPBoundedCVLSConditionalQuadContext np_bounded_cvls_conditional_quad_ctx = {
  0, 0, 0, 0, {0.0, 0.0}, {0.0, 0.0}, NULL, NULL
};

void np_bounded_cvls_conditional_quad_context_clear_extern(void){
  if(np_bounded_cvls_conditional_quad_ctx.base_grid != NULL)
    free(np_bounded_cvls_conditional_quad_ctx.base_grid);
  if(np_bounded_cvls_conditional_quad_ctx.base_weights != NULL)
    free(np_bounded_cvls_conditional_quad_ctx.base_weights);

  np_bounded_cvls_conditional_quad_ctx.ready = 0;
  np_bounded_cvls_conditional_quad_ctx.required = 0;
  np_bounded_cvls_conditional_quad_ctx.q = 0;
  np_bounded_cvls_conditional_quad_ctx.num_obs = 0;
  np_bounded_cvls_conditional_quad_ctx.quad_lb[0] = 0.0;
  np_bounded_cvls_conditional_quad_ctx.quad_lb[1] = 0.0;
  np_bounded_cvls_conditional_quad_ctx.quad_ub[0] = 0.0;
  np_bounded_cvls_conditional_quad_ctx.quad_ub[1] = 0.0;
  np_bounded_cvls_conditional_quad_ctx.base_grid = NULL;
  np_bounded_cvls_conditional_quad_ctx.base_weights = NULL;
}

int np_bounded_cvls_conditional_quad_context_prepare_extern(void){
  const int num_obs = num_obs_train_extern;
  const int q = np_bounded_cvls_conditional_grid_points(1);
  int q_actual = 0;

  np_bounded_cvls_conditional_quad_context_clear_extern();

  if(!np_conditional_density_cvls_bounded_scalar_route_ok())
    return 0;
  if((matrix_Y_continuous_train_extern == NULL) ||
     (matrix_Y_continuous_train_extern[0] == NULL) ||
     (vector_cykerlb_extern == NULL) ||
     (vector_cykerub_extern == NULL) ||
     (num_obs <= 0) ||
     (q < 2))
    return 1;

  np_bounded_cvls_conditional_quad_ctx.base_grid = alloc_vecd(q);
  np_bounded_cvls_conditional_quad_ctx.base_weights = alloc_vecd(q);
  if((np_bounded_cvls_conditional_quad_ctx.base_grid == NULL) ||
     (np_bounded_cvls_conditional_quad_ctx.base_weights == NULL)){
    np_bounded_cvls_conditional_quad_context_clear_extern();
    return 1;
  }

  np_bounded_cvls_conditional_effective_integration_bounds(
    1,
    matrix_Y_continuous_train_extern,
    num_obs,
    vector_cykerlb_extern,
    vector_cykerub_extern,
    np_bounded_cvls_conditional_quad_ctx.quad_lb,
    np_bounded_cvls_conditional_quad_ctx.quad_ub
  );

  if(np_bounded_cvls_build_conditional_grid_1d(
       matrix_Y_continuous_train_extern[0],
       num_obs,
       np_bounded_cvls_conditional_quad_ctx.quad_lb[0],
       np_bounded_cvls_conditional_quad_ctx.quad_ub[0],
       q,
       np_bounded_cvls_conditional_quad_ctx.base_grid,
       np_bounded_cvls_conditional_quad_ctx.base_weights,
       &q_actual) != 0){
    np_bounded_cvls_conditional_quad_context_clear_extern();
    return 1;
  }

  np_bounded_cvls_conditional_quad_ctx.ready = 1;
  np_bounded_cvls_conditional_quad_ctx.required = 1;
  np_bounded_cvls_conditional_quad_ctx.q = q_actual;
  np_bounded_cvls_conditional_quad_ctx.num_obs = num_obs;

  return 0;
}

static int np_density_cvls_bounded_scalar_route_ok(const int num_reg_unordered,
                                                   const int num_reg_ordered,
                                                   const int num_reg_continuous){
  if(int_cker_bound_extern == 0)
    return 0;
  if(num_reg_continuous != 1)
    return 0;
  if(num_reg_unordered != 0)
    return 0;
  if(num_reg_ordered != 0)
    return 0;
  return np_bounded_cvls_continuous_support_ok(num_reg_continuous,
                                               vector_ckerlb_extern,
                                               vector_ckerub_extern);
}

static int np_conditional_density_cvls_bounded_general_route_ok(void){
  if(int_cyker_bound_extern == 0)
    return 0;
  return np_bounded_cvls_continuous_support_ok(num_var_continuous_extern,
                                               vector_cykerlb_extern,
                                               vector_cykerub_extern);
}

static int np_conditional_density_cvls_bounded_i1_eval_on_grid(double *vector_scale_factor,
                                                               NPConditionalXRowCtx *xctx,
                                                               NPConditionalYRowCtx *yctx,
                                                               const double *grid,
                                                               const double *weights,
                                                               const int q,
                                                               const int i1_mode,
                                                               double **xblock,
                                                               double **xblock_full,
                                                               double *yrow,
	                                                               double **ygridblock,
	                                                               double *fit_block,
	                                                               double *lin_block,
	                                                               const int block_size,
	                                                               double *cv,
	                                                               double *i1_mean){
	  const int num_obs = num_obs_train_extern;
	  double local_cv = 0.0, local_i1_mean = 0.0;
	  int i0, m;
	  int local_fail = 0;
#ifdef MPI2
	  const int use_parallel_blocks = (iNum_Processors > 1) && !np_mpi_local_regression_active();
#else
	  const int use_parallel_blocks = 0;
#endif

	  if((vector_scale_factor == NULL) || (xctx == NULL) || (yctx == NULL) ||
	     (grid == NULL) || (weights == NULL) || (xblock == NULL) || (yrow == NULL) ||
	     (ygridblock == NULL) || (fit_block == NULL) || (lin_block == NULL) ||
	     (cv == NULL) || (i1_mean == NULL))
	    return 1;
  if((i1_mode == NP_BOUNDED_CVLS_I1_MODE_FULL) && (xblock_full == NULL))
    return 1;
  if((q < 2) || (block_size <= 0))
    return 1;

	  for(m = 0; m < q; m++){
	    if(np_conditional_y_scalar_eval_from_ctx(vector_scale_factor,
	                                             yctx,
	                                             grid[m],
	                                             ygridblock[m]) != 0)
	      local_fail = 1;
	  }

	  for(i0 = 0; (i0 < num_obs) && !local_fail; i0 += block_size){
	    const int ib = MIN(block_size, num_obs - i0);
	    const int block_id = i0 / block_size;
	    double ** const xuse =
	      (i1_mode == NP_BOUNDED_CVLS_I1_MODE_BOOK) ? xblock : xblock_full;
	    int b;

	    if(use_parallel_blocks && ((block_id % iNum_Processors) != my_rank))
	      continue;

	    for(b = 0; b < ib; b++){
	      const int i = i0 + b;

	      if(np_conditional_xrow_from_ctx(xctx, i, xblock[b]) != 0){
	        local_fail = 1;
	        break;
	      }
	      if((i1_mode == NP_BOUNDED_CVLS_I1_MODE_FULL) &&
	         (np_conditional_x_weight_row_full_stream_core(vector_scale_factor, i, xblock_full[b]) != 0)){
	        local_fail = 1;
	        break;
	      }
	      if(np_conditional_yrow_from_ctx(yctx, i, yrow) != 0){
	        local_fail = 1;
	        break;
	      }

	      lin_block[b] = np_blas_ddot_int(num_obs, xblock[b], yrow);
	    }
	    if(local_fail)
	      break;

	    np_blas_dgemm_tn_int(q, ib, num_obs, ygridblock[0], xuse[0], fit_block);

	    for(b = 0; b < ib; b++){
	      double quad = 0.0;
	      const int offset = b*q;

	      for(m = 0; m < q; m++){
	        const double fit = fit_block[offset + m];
	        quad += weights[m]*fit*fit;
	      }
	      local_i1_mean += quad;
	      local_cv += quad - 2.0*lin_block[b];
	    }
	  }

#ifdef MPI2
	  if(use_parallel_blocks){
	    int any_fail = 0;
	    MPI_Allreduce(&local_fail, &any_fail, 1, MPI_INT, MPI_MAX, comm[1]);
	    if(any_fail)
	      return 1;
	    MPI_Allreduce(&local_cv, cv, 1, MPI_DOUBLE, MPI_SUM, comm[1]);
	    MPI_Allreduce(&local_i1_mean, i1_mean, 1, MPI_DOUBLE, MPI_SUM, comm[1]);
	  } else
#endif
	  {
	    if(local_fail)
	      return 1;
	    *cv = local_cv;
	    *i1_mean = local_i1_mean;
	  }

	  *i1_mean /= (double)num_obs;
	  *cv /= (double)num_obs;

	  return 0;
}

static int np_density_cvls_bounded_general_route_ok(const int num_reg_continuous){
  if(int_cker_bound_extern == 0)
    return 0;
  return np_bounded_cvls_continuous_support_ok(num_reg_continuous,
                                               vector_ckerlb_extern,
                                               vector_ckerub_extern);
}

static int np_bounded_cvls_eval_count(const int ncon,
                                      const int nuno,
                                      const int nord,
                                      const int q,
                                      const int *num_categories,
                                      size_t *total_out){
  size_t total = 1;
  int l;

  if(total_out == NULL)
    return 1;
  if((ncon < 1) || (ncon > 2) || (q < 2))
    return 1;

  for(l = 0; l < ncon; l++){
    if(total > (SIZE_MAX/(size_t)q))
      return 1;
    total *= (size_t)q;
  }

  for(l = 0; l < (nuno + nord); l++){
    const int ncat = num_categories[l];
    if(ncat <= 0)
      return 1;
    if(total > (SIZE_MAX/(size_t)ncat))
      return 1;
    total *= (size_t)ncat;
  }

  *total_out = total;
  return 0;
}

static void np_bounded_cvls_fill_eval_block(const size_t eval_start,
                                            const int block_rows,
                                            const int ncon,
                                            const int nuno,
                                            const int nord,
                                            const int q,
                                            double **cont_grid,
                                            double **cont_weight,
                                            int *num_categories,
                                            double **matrix_categorical_vals,
                                            double **eval_uno,
                                            double **eval_ord,
                                            double **eval_con,
                                            double *eval_weight){
  const int ndisc = nuno + nord;
  size_t disc_total = 1;
  int l, b;

  for(l = 0; l < ndisc; l++)
    disc_total *= (size_t)num_categories[l];

  for(b = 0; b < block_rows; b++){
    size_t idx = eval_start + (size_t)b;
    size_t disc_idx = (ndisc > 0) ? (idx % disc_total) : 0;
    size_t cont_idx = (ndisc > 0) ? (idx / disc_total) : idx;
    double w = 1.0;

    for(l = 0; l < ncon; l++){
      const int g = (int)(cont_idx % (size_t)q);
      cont_idx /= (size_t)q;
      eval_con[l][b] = cont_grid[l][g];
      w *= cont_weight[l][g];
    }

    for(l = 0; l < nuno; l++){
      const int lev = (int)(disc_idx % (size_t)num_categories[l]);
      disc_idx /= (size_t)num_categories[l];
      eval_uno[l][b] = matrix_categorical_vals[l][lev];
    }

    for(l = 0; l < nord; l++){
      const int cat_idx = l + nuno;
      const int lev = (int)(disc_idx % (size_t)num_categories[cat_idx]);
      disc_idx /= (size_t)num_categories[cat_idx];
      eval_ord[l][b] = matrix_categorical_vals[cat_idx][lev];
    }

    eval_weight[b] = w;
  }
}

static int np_density_cvls_bounded_i1_quadrature(const int KERNEL_den,
                                                 const int BANDWIDTH_den,
                                                 const int num_obs,
                                                 const int num_reg_unordered,
                                                 const int num_reg_ordered,
                                                 const int num_reg_continuous,
                                                 double **matrix_X_continuous,
                                                 double *vector_scale_factor,
                                                 double *cv1){
  const int q = NP_BOUNDED_CVLS_I1_GRID_POINTS;
  const double lb = vector_ckerlb_extern[0];
  const double ub = vector_ckerub_extern[0];
  const double * const train_x = matrix_X_continuous[0];
  double *grid = NULL, *weights = NULL;
  double **matrix_eval_continuous = NULL;
  double **matrix_bandwidth = NULL;
  double *lambda = NULL;
  int m;
  int status = 1;

  if((cv1 == NULL) || (vector_scale_factor == NULL) || (num_obs <= 0))
    return 1;
  if((BANDWIDTH_den != BW_FIXED) &&
     (BANDWIDTH_den != BW_GEN_NN) &&
     (BANDWIDTH_den != BW_ADAP_NN))
    return 1;
  if(!np_density_cvls_bounded_scalar_route_ok(num_reg_unordered,
                                              num_reg_ordered,
                                              num_reg_continuous))
    return 1;

  grid = alloc_vecd(MAX(2, q));
  weights = alloc_vecd(MAX(2, q));
  matrix_eval_continuous = alloc_matd(q, 1);
  matrix_bandwidth = alloc_matd((BANDWIDTH_den == BW_ADAP_NN) ? num_obs :
                                ((BANDWIDTH_den == BW_GEN_NN) ? q : 1),
                                1);
  if((num_reg_unordered + num_reg_ordered) > 0)
    lambda = alloc_vecd(num_reg_unordered + num_reg_ordered);
  if((grid == NULL) || (weights == NULL) || (matrix_eval_continuous == NULL) ||
     (matrix_bandwidth == NULL))
    goto cleanup_density_bounded_quad;

  np_fill_trapezoid_rule(lb, ub, q, grid, weights);
  for(m = 0; m < q; m++)
    matrix_eval_continuous[0][m] = grid[m];

  if(kernel_bandwidth_mean(KERNEL_den,
                           BANDWIDTH_den,
                           num_obs,
                           q,
                           0,
                           0,
                           0,
                           1,
                           0,
                           0,
                           0,
                           vector_scale_factor,
                           NULL,
                           NULL,
                           matrix_X_continuous,
                           matrix_eval_continuous,
                           NULL,
                           matrix_bandwidth,
                           lambda) == 1)
    goto cleanup_density_bounded_quad;

  *cv1 = 0.0;
  for(m = 0; m < q; m++){
    double fit = 0.0;

    if((BANDWIDTH_den == BW_FIXED) || (BANDWIDTH_den == BW_GEN_NN)){
      const double h = matrix_bandwidth[0][(BANDWIDTH_den == BW_FIXED) ? 0 : m];
      const double invnorm = np_cker_invnorm(KERNEL_den, grid[m], h, lb, ub);

      if(!(h > 0.0) || !isfinite(h))
        goto cleanup_density_bounded_quad;
      for(int i = 0; i < num_obs; i++)
        fit += invnorm*np_cker_base_eval(KERNEL_den, (grid[m] - train_x[i])/h)/h;
      fit /= (double)num_obs;
    } else {
      for(int i = 0; i < num_obs; i++){
        const double h = matrix_bandwidth[0][i];
        const double invnorm = np_cker_invnorm(KERNEL_den, grid[m], h, lb, ub);

        if(!(h > 0.0) || !isfinite(h))
          goto cleanup_density_bounded_quad;
        fit += invnorm*np_cker_base_eval(KERNEL_den, (grid[m] - train_x[i])/h) /
          ((double)num_obs*h);
      }
    }

    *cv1 += weights[m]*fit*fit;
  }

  status = 0;

cleanup_density_bounded_quad:
  if(grid != NULL) free(grid);
  if(weights != NULL) free(weights);
  if(matrix_eval_continuous != NULL) free_mat(matrix_eval_continuous, 1);
  if(matrix_bandwidth != NULL) free_mat(matrix_bandwidth, 1);
  if(lambda != NULL) free(lambda);
  return status;
}

static int np_density_cvls_bounded_i1_quadrature_general(const int KERNEL_den,
                                                         const int KERNEL_unordered_den,
                                                         const int KERNEL_ordered_den,
                                                         const int BANDWIDTH_den,
                                                         const int num_obs,
                                                         const int num_reg_unordered,
                                                         const int num_reg_ordered,
                                                         const int num_reg_continuous,
                                                         double **matrix_X_unordered,
                                                         double **matrix_X_ordered,
                                                         double **matrix_X_continuous,
                                                         double *vector_scale_factor,
                                                         int *num_categories,
                                                         double **matrix_categorical_vals,
                                                         double *cv1){
  const int ncon = num_reg_continuous;
  const int nuno = num_reg_unordered;
  const int nord = num_reg_ordered;
  const int q = np_bounded_cvls_grid_points(ncon);
  const int block_size = 64;
  size_t total_eval = 0;
  double **cont_grid = NULL, **cont_weight = NULL;
  double **eval_xuno = NULL, **eval_xord = NULL, **eval_xcon = NULL;
  double **matrix_bandwidth = NULL;
  double *lambda = NULL, *eval_weight = NULL;
  int bw_rows, d;
  int status = 1;

  if((cv1 == NULL) || (vector_scale_factor == NULL) || (num_obs <= 0))
    return 1;
  if((BANDWIDTH_den != BW_FIXED) &&
     (BANDWIDTH_den != BW_GEN_NN) &&
     (BANDWIDTH_den != BW_ADAP_NN))
    return 1;
  if(!np_density_cvls_bounded_general_route_ok(num_reg_continuous))
    return 1;
  if(np_bounded_cvls_eval_count(ncon,
                                nuno,
                                nord,
                                q,
                                num_categories,
                                &total_eval) != 0)
    return 1;

  bw_rows =
    (BANDWIDTH_den == BW_FIXED) ? 1 :
    ((BANDWIDTH_den == BW_GEN_NN) ? block_size : num_obs);

  cont_grid = alloc_matd(q, ncon);
  cont_weight = alloc_matd(q, ncon);
  if(nuno > 0) eval_xuno = alloc_matd(block_size, nuno);
  if(nord > 0) eval_xord = alloc_matd(block_size, nord);
  eval_xcon = alloc_matd(block_size, ncon);
  matrix_bandwidth = alloc_matd(bw_rows, ncon);
  lambda = alloc_vecd(nuno + nord);
  eval_weight = alloc_vecd(block_size);

  if((cont_grid == NULL) || (cont_weight == NULL) || (eval_xcon == NULL) ||
     (matrix_bandwidth == NULL) || (lambda == NULL) || (eval_weight == NULL) ||
     ((nuno > 0) && (eval_xuno == NULL)) ||
     ((nord > 0) && (eval_xord == NULL)))
    goto cleanup_density_bounded_quad_general;

  for(d = 0; d < ncon; d++)
    np_fill_trapezoid_rule(vector_ckerlb_extern[d],
                           vector_ckerub_extern[d],
                           q,
                           cont_grid[d],
                           cont_weight[d]);

  *cv1 = 0.0;
  for(size_t eval_start = 0; eval_start < total_eval; ){
    const int eb = (int)MIN((size_t)block_size, total_eval - eval_start);
    int b, i, l;

    np_bounded_cvls_fill_eval_block(eval_start,
                                    eb,
                                    ncon,
                                    nuno,
                                    nord,
                                    q,
                                    cont_grid,
                                    cont_weight,
                                    num_categories,
                                    matrix_categorical_vals,
                                    eval_xuno,
                                    eval_xord,
                                    eval_xcon,
                                    eval_weight);

    if(kernel_bandwidth_mean(KERNEL_den,
                             BANDWIDTH_den,
                             num_obs,
                             eb,
                             0,
                             0,
                             0,
                             ncon,
                             nuno,
                             nord,
                             0,
                             vector_scale_factor,
                             NULL,
                             NULL,
                             matrix_X_continuous,
                             eval_xcon,
                             NULL,
                             matrix_bandwidth,
                             lambda) == 1)
      goto cleanup_density_bounded_quad_general;

    for(b = 0; b < eb; b++){
      double fit = 0.0;

      for(i = 0; i < num_obs; i++){
        double prod = 1.0;

        for(d = 0; d < ncon; d++){
          const int bw_idx =
            (BANDWIDTH_den == BW_FIXED) ? 0 :
            ((BANDWIDTH_den == BW_GEN_NN) ? b : i);
          const double h = matrix_bandwidth[d][bw_idx];
          const double eval_x = eval_xcon[d][b];
          const double invnorm =
            np_cker_invnorm(KERNEL_den,
                            eval_x,
                            h,
                            vector_ckerlb_extern[d],
                            vector_ckerub_extern[d]);

          if(!(h > 0.0) || !isfinite(h))
            goto cleanup_density_bounded_quad_general;
          prod *= invnorm*np_cker_base_eval(KERNEL_den,
                                            (eval_x - matrix_X_continuous[d][i])/h)/h;
        }

        for(l = 0; l < nuno; l++)
          prod *= kernel_unordered(KERNEL_unordered_den,
                                   eval_xuno[l][b],
                                   matrix_X_unordered[l][i],
                                   lambda[l],
                                   num_categories[l]);

        for(l = 0; l < nord; l++)
          prod *= kernel_ordered(KERNEL_ordered_den,
                                 eval_xord[l][b],
                                 matrix_X_ordered[l][i],
                                 lambda[l + nuno]);

        fit += prod;
      }

      fit /= (double)num_obs;
      *cv1 += eval_weight[b]*fit*fit;
    }

    eval_start += (size_t)eb;
  }

  status = 0;

cleanup_density_bounded_quad_general:
  if(cont_grid != NULL) free_mat(cont_grid, ncon);
  if(cont_weight != NULL) free_mat(cont_weight, ncon);
  if(eval_xuno != NULL) free_mat(eval_xuno, nuno);
  if(eval_xord != NULL) free_mat(eval_xord, nord);
  if(eval_xcon != NULL) free_mat(eval_xcon, ncon);
  if(matrix_bandwidth != NULL) free_mat(matrix_bandwidth, ncon);
  if(lambda != NULL) free(lambda);
  if(eval_weight != NULL) free(eval_weight);
  return status;
}

static int np_conditional_density_cvls_bounded_i1_quadrature_row_stream(double *vector_scale_factor,
                                                                         double *cv,
                                                                         int i1_mode){
  const int num_obs = num_obs_train_extern;
  const NPBoundedCVLSConditionalQuadContext *quad_ctx =
    &np_bounded_cvls_conditional_quad_ctx;
  const int block_size = MAX(1, MIN(np_conditional_lp_cvls_block_size(), 64));
  int q = 0;
  int use_quad_context = 0;
  double quad_lb[2] = {0.0, 0.0};
  double quad_ub[2] = {0.0, 0.0};
  NPConditionalXRowCtx xctx = {0};
  NPConditionalYRowCtx yctx = {0};
  double *yrow = NULL, *fit_block = NULL, *lin_block = NULL;
  double **xblock = NULL, **xblock_full = NULL;
  const double *base_grid = NULL, *base_weights = NULL;
  double *base_grid_local = NULL, *base_weights_local = NULL;
  double **ygridblock = NULL;
  double i1_mean = 0.0;
  int q_actual = 0;
  int status = 1;

  if((cv == NULL) || (vector_scale_factor == NULL) || (num_obs <= 0))
    return 1;
  if(!np_conditional_density_cvls_bounded_scalar_route_ok())
    return 1;
  if((i1_mode != NP_BOUNDED_CVLS_I1_MODE_BOOK) &&
     (i1_mode != NP_BOUNDED_CVLS_I1_MODE_FULL))
    return 1;
  if((quad_ctx->ready != 0) &&
     (quad_ctx->num_obs == num_obs) &&
     (quad_ctx->q >= 2) &&
     (quad_ctx->base_grid != NULL) &&
     (quad_ctx->base_weights != NULL)){
    use_quad_context = 1;
    q = quad_ctx->q;
    base_grid = quad_ctx->base_grid;
    base_weights = quad_ctx->base_weights;
    quad_lb[0] = quad_ctx->quad_lb[0];
    quad_ub[0] = quad_ctx->quad_ub[0];
  } else {
    if(quad_ctx->required != 0){
      np_bwm_set_deferred_error("bounded npcdens cv.ls quadrature context is missing");
      return 1;
    }
    q = np_bounded_cvls_conditional_grid_points(1);
  }
  if(q < 2)
    return 1;
  if(q > (INT_MAX / block_size))
    return 1;

  xblock = alloc_tmatd(num_obs, block_size);
  if(i1_mode == NP_BOUNDED_CVLS_I1_MODE_FULL)
    xblock_full = alloc_tmatd(num_obs, block_size);
  yrow = alloc_vecd(MAX(1, num_obs));
  fit_block = alloc_vecd(q*block_size);
  lin_block = alloc_vecd(block_size);
  if(!use_quad_context){
    base_grid_local = alloc_vecd(q);
    base_weights_local = alloc_vecd(q);
    base_grid = base_grid_local;
    base_weights = base_weights_local;
  }
  ygridblock = alloc_tmatd(num_obs, q);
  if((xblock == NULL) || (yrow == NULL) || (fit_block == NULL) || (lin_block == NULL) ||
     (base_grid == NULL) || (base_weights == NULL) ||
     (ygridblock == NULL) ||
     ((i1_mode == NP_BOUNDED_CVLS_I1_MODE_FULL) && (xblock_full == NULL)))
    goto cleanup_bounded_cvls_quad;

  if(np_conditional_xrow_ctx_prepare(vector_scale_factor, &xctx) != 0)
    goto cleanup_bounded_cvls_quad;
  if(np_conditional_yrow_ctx_prepare(vector_scale_factor, OP_NORMAL, &yctx) != 0)
    goto cleanup_bounded_cvls_quad;

  if(!use_quad_context){
    np_bounded_cvls_conditional_effective_integration_bounds(
      1,
      matrix_Y_continuous_train_extern,
      num_obs,
      vector_cykerlb_extern,
      vector_cykerub_extern,
      quad_lb,
      quad_ub
    );

    /* Evaluate the bounded I1 term by integrating the pointwise estimator itself,
       bypassing the suspect bounded-convolution shortcut on this narrow surface. */
    if(np_bounded_cvls_build_conditional_grid_1d(
         matrix_Y_continuous_train_extern[0],
         num_obs,
         quad_lb[0],
         quad_ub[0],
         q,
         base_grid_local,
         base_weights_local,
         &q_actual) != 0)
      goto cleanup_bounded_cvls_quad;
    q = q_actual;
  }
  if(np_conditional_density_cvls_bounded_i1_eval_on_grid(vector_scale_factor,
                                                         &xctx,
                                                         &yctx,
                                                         base_grid,
                                                         base_weights,
                                                         q,
                                                         i1_mode,
                                                         xblock,
                                                         xblock_full,
                                                         yrow,
                                                         ygridblock,
                                                         fit_block,
                                                         lin_block,
                                                         block_size,
                                                         cv,
                                                         &i1_mean) != 0)
    goto cleanup_bounded_cvls_quad;
  status = 0;

cleanup_bounded_cvls_quad:
  np_conditional_xrow_ctx_clear(&xctx);
  np_conditional_yrow_ctx_clear(&yctx);
  np_glp_cv_clear_extern();
  if(xblock != NULL) free_tmat(xblock);
  if(xblock_full != NULL) free_tmat(xblock_full);
  if(yrow != NULL) free(yrow);
  if(fit_block != NULL) free(fit_block);
  if(lin_block != NULL) free(lin_block);
  if(base_grid_local != NULL) free(base_grid_local);
  if(base_weights_local != NULL) free(base_weights_local);
  if(ygridblock != NULL) free_tmat(ygridblock);
  return status;
}

static int np_conditional_y_eval_any_block_stream_core(double *vector_scale_factor,
                                                       int block_rows,
                                                       double **matrix_Y_unordered_eval,
                                                       double **matrix_Y_ordered_eval,
                                                       double **matrix_Y_continuous_eval,
                                                       double **rows_out){
  int b;

  if((rows_out == NULL) || (vector_scale_factor == NULL) || (block_rows <= 0))
    return 1;

  if((int_TREE_Y != NP_TREE_TRUE) && (BANDWIDTH_den_extern != BW_ADAP_NN)){
    return np_conditional_y_eval_block_stream_op_core(vector_scale_factor,
                                                      0,
                                                      block_rows,
                                                      OP_NORMAL,
                                                      matrix_Y_unordered_eval,
                                                      matrix_Y_ordered_eval,
                                                      matrix_Y_continuous_eval,
                                                      block_rows,
                                                      rows_out);
  }

  for(b = 0; b < block_rows; b++){
    if(np_conditional_y_eval_row_stream_op_core(vector_scale_factor,
                                                b,
                                                OP_NORMAL,
                                                matrix_Y_unordered_eval,
                                                matrix_Y_ordered_eval,
                                                matrix_Y_continuous_eval,
                                                block_rows,
                                                0,
                                                rows_out[b]) != 0)
      return 1;
  }

  return 0;
}

static int np_conditional_density_cvls_bounded_i1_quadrature_general_row_stream(double *vector_scale_factor,
                                                                                 double *cv,
                                                                                 int i1_mode){
  const int num_obs = num_obs_train_extern;
  const int ncon = num_var_continuous_extern;
  const int nuno = num_var_unordered_extern;
  const int nord = num_var_ordered_extern;
  int q = np_bounded_cvls_conditional_grid_points(ncon);
  const int block_size = MAX(1, MIN(np_conditional_lp_cvls_block_size(), 64));
  size_t total_eval = 0;
  NPConditionalYRowCtx yctx = {0};
  double *yrow = NULL, *eval_weight = NULL, *fit_block = NULL, *lin_block = NULL, *quad_block = NULL;
  double **xblock = NULL, **xblock_full = NULL;
  double **cont_grid = NULL, **cont_weight = NULL;
  double **eval_yuno = NULL, **eval_yord = NULL, **eval_ycon = NULL, **yevalblock = NULL;
  double quad_lb[2] = {0.0, 0.0};
  double quad_ub[2] = {0.0, 0.0};
  double local_cv = 0.0;
  int i0, d;
  int status = 1;
  int local_fail = 0;
#ifdef MPI2
  const int use_parallel_blocks = (iNum_Processors > 1) && !np_mpi_local_regression_active();
#else
  const int use_parallel_blocks = 0;
#endif

  if((cv == NULL) || (vector_scale_factor == NULL) || (num_obs <= 0))
    return 1;
  if(!np_conditional_density_cvls_bounded_general_route_ok())
    return 1;
  if((i1_mode != NP_BOUNDED_CVLS_I1_MODE_BOOK) &&
     (i1_mode != NP_BOUNDED_CVLS_I1_MODE_FULL))
    return 1;
  if(q < 2)
    return 1;

  xblock = alloc_tmatd(num_obs, block_size);
  if(i1_mode == NP_BOUNDED_CVLS_I1_MODE_FULL)
    xblock_full = alloc_tmatd(num_obs, block_size);
  yrow = alloc_vecd(MAX(1, num_obs));
  eval_weight = alloc_vecd(block_size);
  fit_block = alloc_vecd(block_size*block_size);
  lin_block = alloc_vecd(block_size);
  quad_block = alloc_vecd(block_size);
  cont_grid = alloc_matd(q, ncon);
  cont_weight = alloc_matd(q, ncon);
  if(nuno > 0) eval_yuno = alloc_matd(block_size, nuno);
  if(nord > 0) eval_yord = alloc_matd(block_size, nord);
  eval_ycon = alloc_matd(block_size, ncon);
  yevalblock = alloc_tmatd(num_obs, block_size);

  if((xblock == NULL) || (yrow == NULL) || (eval_weight == NULL) ||
     (fit_block == NULL) || (lin_block == NULL) || (quad_block == NULL) ||
     (cont_grid == NULL) || (cont_weight == NULL) || (eval_ycon == NULL) ||
     (yevalblock == NULL) ||
     ((nuno > 0) && (eval_yuno == NULL)) ||
     ((nord > 0) && (eval_yord == NULL)) ||
     ((i1_mode == NP_BOUNDED_CVLS_I1_MODE_FULL) && (xblock_full == NULL)))
    goto cleanup_bounded_cvls_quad_general;

  np_bounded_cvls_conditional_effective_integration_bounds(
    ncon,
    matrix_Y_continuous_train_extern,
    num_obs,
    vector_cykerlb_extern,
    vector_cykerub_extern,
    quad_lb,
    quad_ub
  );

  if((ncon == 1) &&
     (int_bounded_cvls_quadrature_grid_extern != NP_BOUNDED_CVLS_GRID_UNIFORM)){
    int q_actual = 0;
    if(np_bounded_cvls_build_conditional_grid_1d(
         matrix_Y_continuous_train_extern[0],
         num_obs,
         quad_lb[0],
         quad_ub[0],
         q,
         cont_grid[0],
         cont_weight[0],
         &q_actual) != 0)
      goto cleanup_bounded_cvls_quad_general;
    q = q_actual;
  } else {
    for(d = 0; d < ncon; d++)
      np_fill_trapezoid_rule(quad_lb[d],
                             quad_ub[d],
                             q,
                             cont_grid[d],
                             cont_weight[d]);
  }

  if(np_bounded_cvls_eval_count(ncon,
                                nuno,
                                nord,
                                q,
                                num_categories_extern_Y,
                                &total_eval) != 0)
    goto cleanup_bounded_cvls_quad_general;

  if(np_conditional_yrow_ctx_prepare(vector_scale_factor, OP_NORMAL, &yctx) != 0)
    goto cleanup_bounded_cvls_quad_general;

  for(i0 = 0; (i0 < num_obs) && !local_fail; i0 += block_size){
    const int ib = MIN(block_size, num_obs - i0);
    const int block_id = i0 / block_size;
    double ** const xuse =
      (i1_mode == NP_BOUNDED_CVLS_I1_MODE_BOOK) ? xblock : xblock_full;
    size_t eval_start = 0;
    int b;

    if(use_parallel_blocks && ((block_id % iNum_Processors) != my_rank))
      continue;

    for(b = 0; b < ib; b++){
      const int i = i0 + b;

      if(np_conditional_x_weight_row_stream_core(vector_scale_factor, i, xblock[b]) != 0){
        local_fail = 1;
        break;
      }
      if((i1_mode == NP_BOUNDED_CVLS_I1_MODE_FULL) &&
         (np_conditional_x_weight_row_full_stream_core(vector_scale_factor, i, xblock_full[b]) != 0)){
        local_fail = 1;
        break;
      }
      if(np_conditional_yrow_from_ctx(&yctx, i, yrow) != 0){
        local_fail = 1;
        break;
      }

      lin_block[b] = np_blas_ddot_int(num_obs, xblock[b], yrow);
      quad_block[b] = 0.0;
    }
    if(local_fail)
      break;

    while((eval_start < total_eval) && !local_fail){
      const int eb = (int)MIN((size_t)block_size, total_eval - eval_start);

      np_bounded_cvls_fill_eval_block(eval_start,
                                      eb,
                                      ncon,
                                      nuno,
                                      nord,
                                      q,
                                      cont_grid,
                                      cont_weight,
                                      num_categories_extern_Y,
                                      matrix_categorical_vals_extern_Y,
                                      eval_yuno,
                                      eval_yord,
                                      eval_ycon,
                                      eval_weight);

      if(np_conditional_y_eval_any_block_stream_core(vector_scale_factor,
                                                     eb,
                                                     eval_yuno,
                                                     eval_yord,
                                                     eval_ycon,
                                                     yevalblock) != 0){
        local_fail = 1;
        break;
      }

      np_blas_dgemm_tn_int(eb, ib, num_obs, yevalblock[0], xuse[0], fit_block);
      for(b = 0; b < ib; b++){
        const int offset = b*eb;
        int e;

        for(e = 0; e < eb; e++){
          const double fit = fit_block[offset + e];
          quad_block[b] += eval_weight[e]*fit*fit;
        }
      }

      eval_start += (size_t)eb;
    }
    if(local_fail)
      break;

    for(b = 0; b < ib; b++)
      local_cv += quad_block[b] - 2.0*lin_block[b];
  }

#ifdef MPI2
  if(use_parallel_blocks){
    int any_fail = 0;
    MPI_Allreduce(&local_fail, &any_fail, 1, MPI_INT, MPI_MAX, comm[1]);
    if(any_fail)
      goto cleanup_bounded_cvls_quad_general;
    MPI_Allreduce(&local_cv, cv, 1, MPI_DOUBLE, MPI_SUM, comm[1]);
  } else
#endif
  {
    if(local_fail)
      goto cleanup_bounded_cvls_quad_general;
    *cv = local_cv;
  }

  *cv /= (double)num_obs;
  status = 0;

cleanup_bounded_cvls_quad_general:
  np_conditional_yrow_ctx_clear(&yctx);
  np_glp_cv_clear_extern();
  if(xblock != NULL) free_tmat(xblock);
  if(xblock_full != NULL) free_tmat(xblock_full);
  if(yrow != NULL) free(yrow);
  if(eval_weight != NULL) free(eval_weight);
  if(fit_block != NULL) free(fit_block);
  if(lin_block != NULL) free(lin_block);
  if(quad_block != NULL) free(quad_block);
  if(cont_grid != NULL) free_mat(cont_grid, ncon);
  if(cont_weight != NULL) free_mat(cont_weight, ncon);
  if(eval_yuno != NULL) free_mat(eval_yuno, nuno);
  if(eval_yord != NULL) free_mat(eval_yord, nord);
  if(eval_ycon != NULL) free_mat(eval_ycon, ncon);
  if(yevalblock != NULL) free_tmat(yevalblock);
  return status;
}

typedef struct {
  int ready;
  int num_train;
  int nterms;
  double **basis;
  MATRIX XtXINV;
  double *hdiag;
} NPConditionalLpAllLargeCtx;

static void np_conditional_lp_all_large_ctx_clear(NPConditionalLpAllLargeCtx *ctx){
  if(ctx == NULL)
    return;
  if(ctx->XtXINV != NULL) mat_free(ctx->XtXINV);
  if(ctx->hdiag != NULL) free(ctx->hdiag);
  memset(ctx, 0, sizeof(*ctx));
}

static int np_conditional_lp_all_large_ctx_prepare(double *vector_scale_factor,
                                                   NPConditionalLpAllLargeCtx *ctx){
  const int num_train = num_obs_train_extern;
  const int num_reg_tot = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
  const int use_bernstein = (int_glp_bernstein_extern != 0);
  int *kernel_cx = NULL, *kernel_ux = NULL, *kernel_ox = NULL;
  int *ov_cont_ok = NULL;
  double *vsfx = NULL, *lambdax = NULL, *ov_cont_hmin = NULL, *ov_cont_k0 = NULL;
  double **matrix_bandwidth_x = NULL;
  MATRIX XtX = NULL;
  int ov_cont_from_cache = 0;
  int i, a, b;
  int status = 1;

  if((ctx == NULL) || (vector_scale_factor == NULL))
    return 1;

  np_conditional_lp_all_large_ctx_clear(ctx);

  if(int_ll_extern != LL_LP)
    return 1;
  if(BANDWIDTH_den_extern != BW_FIXED)
    return 1;
  if((int_TREE_X == NP_TREE_TRUE) || (int_TREE_Y == NP_TREE_TRUE))
    return 1;
  if((num_train <= 0) || (num_reg_continuous_extern <= 0))
    return 1;
  if(vector_glp_degree_extern == NULL)
    return 1;

  vsfx = alloc_vecd(MAX(1, num_reg_tot));
  lambdax = alloc_vecd(MAX(1, num_reg_unordered_extern + num_reg_ordered_extern));
  matrix_bandwidth_x = alloc_tmatd(1, num_reg_continuous_extern);
  kernel_cx = (int *)calloc((size_t)MAX(1, num_reg_continuous_extern), sizeof(int));
  kernel_ux = (int *)calloc((size_t)MAX(1, num_reg_unordered_extern), sizeof(int));
  kernel_ox = (int *)calloc((size_t)MAX(1, num_reg_ordered_extern), sizeof(int));

  if((vsfx == NULL) || (lambdax == NULL) ||
     ((num_reg_continuous_extern > 0) && (matrix_bandwidth_x == NULL)) ||
     (kernel_cx == NULL) || (kernel_ux == NULL) || (kernel_ox == NULL))
    goto cleanup_all_large_prepare;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        vsfx,
                        NULL,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  for(i = 0; i < num_reg_continuous_extern; i++) kernel_cx[i] = KERNEL_reg_extern;
  for(i = 0; i < num_reg_unordered_extern; i++) kernel_ux[i] = KERNEL_reg_unordered_extern;
  for(i = 0; i < num_reg_ordered_extern; i++) kernel_ox[i] = KERNEL_reg_ordered_extern;

  if(kernel_bandwidth_mean(KERNEL_reg_extern,
                           BANDWIDTH_den_extern,
                           num_train,
                           num_train,
                           0,
                           0,
                           0,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           vsfx,
                           NULL,
                           NULL,
                           matrix_X_continuous_train_extern,
                           matrix_X_continuous_train_extern,
                           NULL,
                           matrix_bandwidth_x,
                           lambdax) == 1)
    goto cleanup_all_large_prepare;

  if(!np_reg_cv_all_large_gate(BANDWIDTH_den_extern,
                               num_train,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               kernel_cx,
                               kernel_ux,
                               kernel_ox,
                               matrix_X_continuous_train_extern,
                               matrix_bandwidth_x,
                               lambdax,
                               num_categories_extern_X,
                               matrix_categorical_vals_extern_X,
                               &ov_cont_ok,
                               &ov_cont_hmin,
                               &ov_cont_k0,
                               &ov_cont_from_cache))
    goto cleanup_all_large_prepare;

  if(!np_glp_cv_cache.ready ||
     (np_glp_cv_cache.use_bernstein != use_bernstein) ||
     (np_glp_cv_cache.basis_mode != int_glp_basis_extern) ||
     (np_glp_cv_cache.num_obs != num_train) ||
     (np_glp_cv_cache.ncon != num_reg_continuous_extern) ||
     (np_glp_cv_cache.matrix_X_continuous_train_ptr != matrix_X_continuous_train_extern)){
    if(!np_glp_cv_cache_prepare(LL_LP,
                                num_train,
                                num_reg_continuous_extern,
                                matrix_X_continuous_train_extern))
      goto cleanup_all_large_prepare;
  }

  if((np_glp_cv_cache.nterms <= 0) || (np_glp_cv_cache.basis == NULL))
    goto cleanup_all_large_prepare;

  XtX = mat_creat(np_glp_cv_cache.nterms, np_glp_cv_cache.nterms, UNDEFINED);
  ctx->XtXINV = mat_creat(np_glp_cv_cache.nterms, np_glp_cv_cache.nterms, UNDEFINED);
  ctx->hdiag = alloc_vecd(MAX(1, num_train));
  if((XtX == NULL) || (ctx->XtXINV == NULL) || (ctx->hdiag == NULL))
    goto cleanup_all_large_prepare;

  for(a = 0; a < np_glp_cv_cache.nterms; a++){
    for(b = 0; b < np_glp_cv_cache.nterms; b++)
      XtX[a][b] = 0.0;
  }

  for(i = 0; i < num_train; i++){
    for(a = 0; a < np_glp_cv_cache.nterms; a++){
      const double za = np_glp_cv_cache.basis[a][i];
      for(b = a; b < np_glp_cv_cache.nterms; b++){
        const double zb = np_glp_cv_cache.basis[b][i];
        XtX[a][b] += za*zb;
        if(b != a)
          XtX[b][a] += za*zb;
      }
    }
  }

  if(np_mat_bad_rcond_sym(XtX, 1.0e-10))
    goto cleanup_all_large_prepare;
  if(mat_inv(XtX, ctx->XtXINV) == NULL)
    goto cleanup_all_large_prepare;

  for(i = 0; i < num_train; i++){
    double hii = 0.0;
    for(a = 0; a < np_glp_cv_cache.nterms; a++){
      double s = 0.0;
      for(b = 0; b < np_glp_cv_cache.nterms; b++)
        s += ctx->XtXINV[a][b]*np_glp_cv_cache.basis[b][i];
      hii += np_glp_cv_cache.basis[a][i]*s;
    }
    if((!isfinite(hii)) || (hii >= 1.0))
      goto cleanup_all_large_prepare;
    ctx->hdiag[i] = hii;
  }

  ctx->ready = 1;
  ctx->num_train = num_train;
  ctx->nterms = np_glp_cv_cache.nterms;
  ctx->basis = np_glp_cv_cache.basis;
  status = 0;

cleanup_all_large_prepare:
  if(XtX != NULL) mat_free(XtX);
  if(vsfx != NULL) free(vsfx);
  if(lambdax != NULL) free(lambdax);
  if(matrix_bandwidth_x != NULL) free_tmat(matrix_bandwidth_x);
  if(kernel_cx != NULL) free(kernel_cx);
  if(kernel_ux != NULL) free(kernel_ux);
  if(kernel_ox != NULL) free(kernel_ox);
  if((ov_cont_ok != NULL) && (!ov_cont_from_cache)) free(ov_cont_ok);
  if((ov_cont_hmin != NULL) && (!ov_cont_from_cache)) free(ov_cont_hmin);
  if((ov_cont_k0 != NULL) && (!ov_cont_from_cache)) free(ov_cont_k0);

  if(status != 0)
    np_conditional_lp_all_large_ctx_clear(ctx);

  return status;
}

static double np_conditional_lp_all_large_row_fit(const NPConditionalLpAllLargeCtx *ctx,
                                                  const double *rhs_row,
                                                  const int eval_pos,
                                                  double *cross_terms,
                                                  double *beta,
                                                  const int leave_one_out){
  double fit = 0.0;
  int a, b;

  for(a = 0; a < ctx->nterms; a++)
    cross_terms[a] = np_blas_ddot_int(ctx->num_train, rhs_row, ctx->basis[a]);

  for(a = 0; a < ctx->nterms; a++){
    double s = 0.0;
    for(b = 0; b < ctx->nterms; b++)
      s += ctx->XtXINV[a][b]*cross_terms[b];
    beta[a] = s;
    fit += ctx->basis[a][eval_pos]*s;
  }

  if(leave_one_out)
    fit = (fit - ctx->hdiag[eval_pos]*rhs_row[eval_pos]) / NZD_POS(1.0 - ctx->hdiag[eval_pos]);

  return fit;
}

static int np_conditional_lp_all_large_build_conv_quad(double *vector_scale_factor,
                                                       const NPConditionalLpAllLargeCtx *ctx,
                                                       MATRIX quad_mat){
  MATRIX moment = NULL, temp = NULL;
  double *yconv = NULL, *cross_terms = NULL;
  int j, a, b, c;
  int status = 1;

  if((vector_scale_factor == NULL) || (ctx == NULL) || (!ctx->ready) || (quad_mat == NULL))
    return 1;

  moment = mat_creat(ctx->nterms, ctx->nterms, UNDEFINED);
  temp = mat_creat(ctx->nterms, ctx->nterms, UNDEFINED);
  yconv = alloc_vecd(MAX(1, ctx->num_train));
  cross_terms = alloc_vecd(MAX(1, ctx->nterms));
  if((moment == NULL) || (temp == NULL) || (yconv == NULL) || (cross_terms == NULL))
    goto cleanup_all_large_quad;

  for(a = 0; a < ctx->nterms; a++)
    for(b = 0; b < ctx->nterms; b++)
      moment[a][b] = 0.0;

  for(j = 0; j < ctx->num_train; j++){
    if(np_conditional_y_row_stream_op_core(vector_scale_factor, j, OP_CONVOLUTION, yconv) != 0)
      goto cleanup_all_large_quad;

    for(b = 0; b < ctx->nterms; b++)
      cross_terms[b] = np_blas_ddot_int(ctx->num_train, yconv, ctx->basis[b]);

    for(a = 0; a < ctx->nterms; a++){
      const double za = ctx->basis[a][j];
      for(b = 0; b < ctx->nterms; b++)
        moment[a][b] += za*cross_terms[b];
    }
  }

  for(a = 0; a < ctx->nterms; a++){
    for(b = 0; b < ctx->nterms; b++){
      double s = 0.0;
      for(c = 0; c < ctx->nterms; c++)
        s += ctx->XtXINV[a][c]*moment[c][b];
      temp[a][b] = s;
    }
  }

  for(a = 0; a < ctx->nterms; a++){
    for(b = 0; b < ctx->nterms; b++){
      double s = 0.0;
      for(c = 0; c < ctx->nterms; c++)
        s += temp[a][c]*ctx->XtXINV[c][b];
      quad_mat[a][b] = s;
    }
  }

  status = 0;

cleanup_all_large_quad:
  if(moment != NULL) mat_free(moment);
  if(temp != NULL) mat_free(temp);
  if(yconv != NULL) free(yconv);
  if(cross_terms != NULL) free(cross_terms);
  return status;
}

static double np_conditional_lp_all_large_quad_eval(const NPConditionalLpAllLargeCtx *ctx,
                                                    MATRIX quad_mat,
                                                    const int eval_pos){
  double quad = 0.0;
  int a, b;

  for(a = 0; a < ctx->nterms; a++){
    const double za = ctx->basis[a][eval_pos];
    for(b = 0; b < ctx->nterms; b++)
      quad += za*quad_mat[a][b]*ctx->basis[b][eval_pos];
  }

  return quad;
}

static int np_conditional_density_cvml_lp_all_large_stream(double *vector_scale_factor,
                                                           double *cv){
  NPConditionalLpAllLargeCtx ctx = {0};
  double *yrow = NULL, *cross_terms = NULL, *beta = NULL;
  int i;
  int status = 1;

  if((cv == NULL) || (vector_scale_factor == NULL))
    return 1;
  if(np_conditional_lp_all_large_ctx_prepare(vector_scale_factor, &ctx) != 0)
    goto cleanup_cvml_all_large;

  yrow = alloc_vecd(MAX(1, ctx.num_train));
  cross_terms = alloc_vecd(MAX(1, ctx.nterms));
  beta = alloc_vecd(MAX(1, ctx.nterms));
  if((yrow == NULL) || (cross_terms == NULL) || (beta == NULL))
    goto cleanup_cvml_all_large;

  *cv = 0.0;
  for(i = 0; i < ctx.num_train; i++){
    double fit;

    if(np_conditional_y_row_stream_op_core(vector_scale_factor, i, OP_NORMAL, yrow) != 0)
      goto cleanup_cvml_all_large;

    fit = np_conditional_lp_all_large_row_fit(&ctx, yrow, i, cross_terms, beta, 1);

    if(fit > DBL_MIN){
      *cv -= log(fit);
    } else if(fit < -DBL_MIN){
      *cv += log(-fit) - 2.0*log(DBL_MIN);
    } else {
      *cv -= log(DBL_MIN);
    }
  }

  status = 0;

cleanup_cvml_all_large:
  if(yrow != NULL) free(yrow);
  if(cross_terms != NULL) free(cross_terms);
  if(beta != NULL) free(beta);
  np_conditional_lp_all_large_ctx_clear(&ctx);
  return status;
}

static int np_conditional_density_cvls_lp_all_large_stream(double *vector_scale_factor,
                                                           double *cv){
  NPConditionalLpAllLargeCtx ctx = {0};
  MATRIX quad_mat = NULL;
  double *yrow = NULL, *cross_terms = NULL, *beta = NULL;
  int i;
  int status = 1;

  if((cv == NULL) || (vector_scale_factor == NULL))
    return 1;
  if(np_conditional_lp_all_large_ctx_prepare(vector_scale_factor, &ctx) != 0)
    goto cleanup_cvls_all_large;

  quad_mat = mat_creat(ctx.nterms, ctx.nterms, UNDEFINED);
  yrow = alloc_vecd(MAX(1, ctx.num_train));
  cross_terms = alloc_vecd(MAX(1, ctx.nterms));
  beta = alloc_vecd(MAX(1, ctx.nterms));
  if((quad_mat == NULL) || (yrow == NULL) || (cross_terms == NULL) || (beta == NULL))
    goto cleanup_cvls_all_large;

  if(np_conditional_lp_all_large_build_conv_quad(vector_scale_factor, &ctx, quad_mat) != 0)
    goto cleanup_cvls_all_large;

  *cv = 0.0;
  for(i = 0; i < ctx.num_train; i++){
    double lin;
    double quad;

    if(np_conditional_y_row_stream_op_core(vector_scale_factor, i, OP_NORMAL, yrow) != 0)
      goto cleanup_cvls_all_large;

    lin = np_conditional_lp_all_large_row_fit(&ctx, yrow, i, cross_terms, beta, 1);
    quad = np_conditional_lp_all_large_quad_eval(&ctx, quad_mat, i);
    *cv += quad - 2.0*lin;
  }

  *cv /= (double)ctx.num_train;
  status = 0;

cleanup_cvls_all_large:
  if(quad_mat != NULL) mat_free(quad_mat);
  if(yrow != NULL) free(yrow);
  if(cross_terms != NULL) free(cross_terms);
  if(beta != NULL) free(beta);
  np_conditional_lp_all_large_ctx_clear(&ctx);
  return status;
}

static int np_conditional_distribution_cvls_lp_all_large_stream(double *vector_scale_factor,
                                                                double *cv){
  NPConditionalLpAllLargeCtx ctx = {0};
  const int num_eval = num_obs_eval_extern;
  double *yint = NULL, *cross_terms = NULL, *beta = NULL;
  int i, j;
  int status = 1;

  if((cv == NULL) || (vector_scale_factor == NULL))
    return 1;
  if(np_conditional_lp_all_large_ctx_prepare(vector_scale_factor, &ctx) != 0)
    goto cleanup_cdist_all_large;

  yint = alloc_vecd(MAX(1, ctx.num_train));
  cross_terms = alloc_vecd(MAX(1, ctx.nterms));
  beta = alloc_vecd(MAX(1, ctx.nterms));
  if((yint == NULL) || (cross_terms == NULL) || (beta == NULL))
    goto cleanup_cdist_all_large;

  *cv = 0.0;
  for(j = 0; j < num_eval; j++){
    if(np_conditional_y_eval_row_stream_op_core(vector_scale_factor,
                                                j,
                                                OP_INTEGRAL,
                                                matrix_Y_unordered_eval_extern,
                                                matrix_Y_ordered_eval_extern,
                                                matrix_Y_continuous_eval_extern,
                                                num_eval,
                                                cdfontrain_extern && (num_eval == ctx.num_train),
                                                yint) != 0)
      goto cleanup_cdist_all_large;

    (void)np_conditional_lp_all_large_row_fit(&ctx, yint, 0, cross_terms, beta, 0);

    for(i = 0; i < ctx.num_train; i++){
      double fit = 0.0;
      const int indy = np_conditional_indicator_row_core(i,
                                                         j,
                                                         cdfontrain_extern,
                                                         matrix_Y_ordered_train_extern,
                                                         matrix_Y_continuous_train_extern,
                                                         matrix_Y_ordered_eval_extern,
                                                         matrix_Y_continuous_eval_extern,
                                                         num_var_ordered_extern,
                                                         num_var_continuous_extern);

      if(cdfontrain_extern && (i == j))
        continue;

      for(int a = 0; a < ctx.nterms; a++)
        fit += ctx.basis[a][i]*beta[a];
      fit = (fit - ctx.hdiag[i]*yint[i]) / NZD_POS(1.0 - ctx.hdiag[i]);

      {
        const double tvd = ((double)indy) - fit;
        *cv += tvd*tvd;
      }
    }
  }

  *cv /= ((double)ctx.num_train*(double)MAX(1, num_eval));
  status = 0;

cleanup_cdist_all_large:
  if(yint != NULL) free(yint);
  if(cross_terms != NULL) free(cross_terms);
  if(beta != NULL) free(beta);
  np_conditional_lp_all_large_ctx_clear(&ctx);
  return status;
}

int np_conditional_density_cvml_lp_stream(double *vector_scale_factor,
                                          double *cv){
  const int num_obs = num_obs_train_extern;
  const int use_cached_gnn = (BANDWIDTH_den_extern == BW_GEN_NN);
  NPConditionalXRowCtx xctx = {0};
  NPConditionalYRowCtx yctx = {0};
  double *xrow = NULL, *yrow = NULL;
  double local_cv = 0.0;
  int i, j;
  int status = 1;
  int local_fail = 0;
#ifdef MPI2
  const int use_parallel_rows = (iNum_Processors > 1) && !np_mpi_local_regression_active();
#else
  const int use_parallel_rows = 0;
#endif

  if((cv == NULL) || (vector_scale_factor == NULL) || (num_obs <= 0))
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;

  if((BANDWIDTH_den_extern == BW_FIXED) &&
     (np_conditional_density_cvml_lp_all_large_stream(vector_scale_factor, cv) == 0)){
    np_glp_cv_clear_extern();
    np_fastcv_alllarge_hits++;
    return 0;
  }

  xrow = alloc_vecd(MAX(1, num_obs));
  yrow = alloc_vecd(MAX(1, num_obs));
  if((xrow == NULL) || (yrow == NULL))
    goto cleanup_cvml_lp_stream;

  if(use_cached_gnn){
    if(np_conditional_xrow_ctx_prepare(vector_scale_factor, &xctx) != 0)
      goto cleanup_cvml_lp_stream;
    if(np_conditional_yrow_ctx_prepare(vector_scale_factor, OP_NORMAL, &yctx) != 0)
      goto cleanup_cvml_lp_stream;
  }

  for(i = 0; (i < num_obs) && !local_fail; i++){
    double fit = 0.0;

    if(use_parallel_rows && ((i % iNum_Processors) != my_rank))
      continue;

    if(use_cached_gnn){
      if(np_conditional_xrow_from_ctx(&xctx, i, xrow) != 0){
        local_fail = 1;
        break;
      }
      if(np_conditional_yrow_from_ctx(&yctx, i, yrow) != 0){
        local_fail = 1;
        break;
      }
    } else {
      if(np_conditional_x_weight_row_stream_core(vector_scale_factor, i, xrow) != 0){
        local_fail = 1;
        break;
      }
      if(np_conditional_y_row_stream_core(vector_scale_factor, i, yrow) != 0){
        local_fail = 1;
        break;
      }
    }

    for(j = 0; j < num_obs; j++)
      fit += xrow[j]*yrow[j];

    /* Match bkcde's default smooth penalty for negative LP delete-one densities. */
    if(fit > DBL_MIN){
      local_cv -= log(fit);
    } else if(fit < -DBL_MIN){
      local_cv += log(-fit) - 2.0*log(DBL_MIN);
    } else {
      local_cv -= log(DBL_MIN);
    }
  }

#ifdef MPI2
  if(use_parallel_rows){
    int any_fail = 0;
    MPI_Allreduce(&local_fail, &any_fail, 1, MPI_INT, MPI_MAX, comm[1]);
    if(any_fail)
      goto cleanup_cvml_lp_stream;
    MPI_Allreduce(&local_cv, cv, 1, MPI_DOUBLE, MPI_SUM, comm[1]);
  } else
#endif
  {
    if(local_fail)
      goto cleanup_cvml_lp_stream;
    *cv = local_cv;
  }

  status = 0;

cleanup_cvml_lp_stream:
  np_conditional_xrow_ctx_clear(&xctx);
  np_conditional_yrow_ctx_clear(&yctx);
  np_glp_cv_clear_extern();
  if(xrow != NULL) free(xrow);
  if(yrow != NULL) free(yrow);
  return status;
}

static int np_conditional_density_cvls_lp_row_stream(double *vector_scale_factor,
                                                     double *cv){
  const int num_obs = num_obs_train_extern;
  double *xrow = NULL, *xrow_full = NULL, *yrow = NULL, *yconv = NULL;
  int i, j, k;
  int status = 1;

  if((cv == NULL) || (vector_scale_factor == NULL) || (num_obs <= 0))
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;

  xrow = alloc_vecd(MAX(1, num_obs));
  xrow_full = alloc_vecd(MAX(1, num_obs));
  yrow = alloc_vecd(MAX(1, num_obs));
  yconv = alloc_vecd(MAX(1, num_obs));
  if((xrow == NULL) || (xrow_full == NULL) || (yrow == NULL) || (yconv == NULL))
    goto cleanup_cvls_lp_stream;

  *cv = 0.0;
  for(i = 0; i < num_obs; i++){
    double lin = 0.0;
    double quad = 0.0;

    if(np_conditional_x_weight_row_stream_core(vector_scale_factor, i, xrow) != 0)
      goto cleanup_cvls_lp_stream;
    if(np_conditional_x_weight_row_full_stream_core(vector_scale_factor, i, xrow_full) != 0)
      goto cleanup_cvls_lp_stream;
    if(np_conditional_y_row_stream_core(vector_scale_factor, i, yrow) != 0)
      goto cleanup_cvls_lp_stream;

    for(j = 0; j < num_obs; j++)
      lin += xrow[j]*yrow[j];

    for(j = 0; j < num_obs; j++){
      double inner = 0.0;
      if(xrow_full[j] == 0.0)
        continue;
      if(np_conditional_y_row_stream_op_core(vector_scale_factor, j, OP_CONVOLUTION, yconv) != 0)
        goto cleanup_cvls_lp_stream;
      for(k = 0; k < num_obs; k++)
        inner += xrow_full[k]*yconv[k];
      quad += xrow_full[j]*inner;
    }

    *cv += quad - 2.0*lin;
  }

  *cv /= (double)num_obs;
  status = 0;

cleanup_cvls_lp_stream:
  if(xrow != NULL) free(xrow);
  if(xrow_full != NULL) free(xrow_full);
  if(yrow != NULL) free(yrow);
  if(yconv != NULL) free(yconv);
  return status;
}

int np_conditional_density_cvls_lp_stream(double *vector_scale_factor,
                                          double *cv){
  const int num_obs = num_obs_train_extern;
  const int block_size = MIN(np_conditional_lp_cvls_block_size(), MAX(1, num_obs));
  const int nblocks = (num_obs + block_size - 1)/block_size;
  double **xblock = NULL, **xblock_full = NULL, **yblock = NULL, **yconvblock = NULL;
  double *block_terms = NULL;
  int i0, j0, ii, jj;
  int status = 1;
  int local_fail = 0;
#ifdef MPI2
  const int use_parallel_blocks = (iNum_Processors > 1);
#else
  const int use_parallel_blocks = 0;
#endif

  if((cv == NULL) || (vector_scale_factor == NULL) || (num_obs <= 0))
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;

  if(np_conditional_density_cvls_bounded_scalar_route_ok())
    return np_conditional_density_cvls_bounded_i1_quadrature_row_stream(vector_scale_factor,
                                                                        cv,
                                                                        NP_BOUNDED_CVLS_I1_MODE_BOOK);
  if(np_conditional_density_cvls_bounded_general_route_ok())
    return np_conditional_density_cvls_bounded_i1_quadrature_general_row_stream(vector_scale_factor,
                                                                                cv,
                                                                                NP_BOUNDED_CVLS_I1_MODE_BOOK);
  if(int_cyker_bound_extern != 0){
    np_bwm_set_deferred_error("bounded npcdens cv.ls currently supports up to two continuous response variables");
    return 1;
  }

  if((BANDWIDTH_den_extern == BW_FIXED) &&
     (np_conditional_density_cvls_lp_all_large_stream(vector_scale_factor, cv) == 0)){
    np_glp_cv_clear_extern();
    np_fastcv_alllarge_hits++;
    return 0;
  }

  if((BANDWIDTH_den_extern == BW_ADAP_NN) ||
     (int_TREE_X == NP_TREE_TRUE) ||
     (int_TREE_Y == NP_TREE_TRUE))
    return np_conditional_density_cvls_lp_row_stream(vector_scale_factor, cv);

  xblock = alloc_matd(num_obs, block_size);
  xblock_full = alloc_matd(num_obs, block_size);
  yblock = alloc_matd(num_obs, block_size);
  yconvblock = alloc_matd(num_obs, block_size);
  if(use_parallel_blocks){
    block_terms = (double *)calloc((size_t)nblocks, sizeof(double));
  }
  if((xblock == NULL) || (xblock_full == NULL) || (yblock == NULL) || (yconvblock == NULL) ||
     (use_parallel_blocks && (block_terms == NULL)))
    local_fail = 1;

  *cv = 0.0;
  for(i0 = 0; i0 < num_obs && !local_fail; i0 += block_size){
    const int ib = MIN(block_size, num_obs - i0);
    const int block_id = i0 / block_size;
    double lin = 0.0;
    double quad = 0.0;

    if(use_parallel_blocks && ((block_id % iNum_Processors) != my_rank))
      continue;

    if(np_conditional_x_weight_block_stream_core(vector_scale_factor, i0, ib, xblock) != 0){
      local_fail = 1;
      break;
    }

    if(np_conditional_x_weight_block_full_stream_core(vector_scale_factor, i0, ib, xblock_full) != 0){
      local_fail = 1;
      break;
    }

    if(np_conditional_y_block_stream_op_core(vector_scale_factor, i0, ib, OP_NORMAL, yblock) != 0){
      local_fail = 1;
      break;
    }

    for(ii = 0; ii < ib; ii++)
      lin += np_blas_ddot_int(num_obs, xblock[ii], yblock[ii]);

    for(j0 = 0; j0 < num_obs; j0 += block_size){
      const int jb = MIN(block_size, num_obs - j0);

      if(np_conditional_y_block_stream_op_core(vector_scale_factor, j0, jb, OP_CONVOLUTION, yconvblock) != 0){
        local_fail = 1;
        break;
      }

      for(ii = 0; ii < ib; ii++){
        double * const ai = xblock_full[ii];
        for(jj = 0; jj < jb; jj++){
          const double aij = ai[j0 + jj];
          double inner;
          if(aij == 0.0)
            continue;
          inner = np_blas_ddot_int(num_obs, ai, yconvblock[jj]);
          quad += aij*inner;
        }
      }
    }

    if(use_parallel_blocks){
      block_terms[block_id] = quad - 2.0*lin;
    } else {
      *cv += quad - 2.0*lin;
    }
  }

#ifdef MPI2
  if(use_parallel_blocks){
    int any_fail = 0;

    if(np_mpi_rank_failure_injected("NP_RMPI_INJECT_CDEN_CVLS_FAIL_RANK"))
      local_fail = 1;

    MPI_Allreduce(&local_fail, &any_fail, 1, MPI_INT, MPI_MAX, comm[1]);
    if(any_fail)
      goto cleanup_cvls_lp_block;

    MPI_Allreduce(MPI_IN_PLACE, block_terms, nblocks, MPI_DOUBLE, MPI_SUM, comm[1]);
    *cv = 0.0;
    for(ii = 0; ii < nblocks; ii++)
      *cv += block_terms[ii];
  }
#endif

  if(local_fail)
    goto cleanup_cvls_lp_block;

  *cv /= (double)num_obs;
  status = 0;

cleanup_cvls_lp_block:
  if(xblock != NULL) free_mat(xblock, block_size);
  if(xblock_full != NULL) free_mat(xblock_full, block_size);
  if(yblock != NULL) free_mat(yblock, block_size);
  if(yconvblock != NULL) free_mat(yconvblock, block_size);
  if(block_terms != NULL) free(block_terms);
  np_glp_cv_clear_extern();
  return status;
}

static int np_conditional_distribution_cvls_lp_row_stream(double *vector_scale_factor,
                                                          double *cv){
  const int num_train = num_obs_train_extern;
  const int num_eval = num_obs_eval_extern;
  double *xrow = NULL, *yint = NULL;
  int i, j, l;
  int status = 1;

  if((cv == NULL) || (vector_scale_factor == NULL) || (num_train <= 0) || (num_eval <= 0))
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;

  if((BANDWIDTH_den_extern == BW_FIXED) &&
     (np_conditional_distribution_cvls_lp_all_large_stream(vector_scale_factor, cv) == 0)){
    np_glp_cv_clear_extern();
    np_fastcv_alllarge_hits++;
    return 0;
  }

  xrow = alloc_vecd(MAX(1, num_train));
  yint = alloc_vecd(MAX(1, num_train));
  if((xrow == NULL) || (yint == NULL))
    goto cleanup_cdist_lp_stream;

  *cv = 0.0;
  for(i = 0; i < num_train; i++){
    if(np_conditional_x_weight_row_stream_core(vector_scale_factor, i, xrow) != 0)
      goto cleanup_cdist_lp_stream;

    for(j = 0; j < num_eval; j++){
      double fit = 0.0;
      const int indy = np_conditional_indicator_row_core(i,
                                                         j,
                                                         cdfontrain_extern,
                                                         matrix_Y_ordered_train_extern,
                                                         matrix_Y_continuous_train_extern,
                                                         matrix_Y_ordered_eval_extern,
                                                         matrix_Y_continuous_eval_extern,
                                                         num_var_ordered_extern,
                                                         num_var_continuous_extern);

      if(np_conditional_y_eval_row_stream_op_core(vector_scale_factor,
                                                  j,
                                                  OP_INTEGRAL,
                                                  matrix_Y_unordered_eval_extern,
                                                  matrix_Y_ordered_eval_extern,
                                                  matrix_Y_continuous_eval_extern,
                                                  num_eval,
                                                  cdfontrain_extern && (num_eval == num_train),
                                                  yint) != 0)
        goto cleanup_cdist_lp_stream;

      if(cdfontrain_extern && (i == j))
        continue;

      for(l = 0; l < num_train; l++)
        fit += xrow[l]*yint[l];

      {
        const double tvd = ((double)indy) - fit;
        *cv += tvd*tvd;
      }
    }
  }

  *cv /= ((double)num_train*(double)MAX(1, num_eval));
  status = 0;

cleanup_cdist_lp_stream:
  if(xrow != NULL) free(xrow);
  if(yint != NULL) free(yint);
  return status;
}

int np_conditional_distribution_cvls_lp_stream(double *vector_scale_factor,
                                               double *cv){
  const int num_train = num_obs_train_extern;
  const int num_eval = num_obs_eval_extern;
  const int block_size = MIN(np_conditional_lp_cvls_block_size(), MAX(1, num_train));
  const int nblocks = (num_train + block_size - 1)/block_size;
  double **xblock = NULL, **yintblock = NULL;
  double *block_terms = NULL;
  int i0, j0, ii, jj;
  int status = 1;
  int local_fail = 0;
#ifdef MPI2
  const int use_parallel_blocks = (iNum_Processors > 1) && !np_mpi_local_regression_active();
#else
  const int use_parallel_blocks = 0;
#endif

  if((cv == NULL) || (vector_scale_factor == NULL) || (num_train <= 0) || (num_eval <= 0))
    return 1;
  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;

  if((BANDWIDTH_den_extern == BW_FIXED) &&
     (np_conditional_distribution_cvls_lp_all_large_stream(vector_scale_factor, cv) == 0)){
    np_glp_cv_clear_extern();
    np_fastcv_alllarge_hits++;
    return 0;
  }

  if((BANDWIDTH_den_extern == BW_ADAP_NN) ||
     (int_TREE_X == NP_TREE_TRUE) ||
     (int_TREE_Y == NP_TREE_TRUE))
    return np_conditional_distribution_cvls_lp_row_stream(vector_scale_factor, cv);

  xblock = alloc_matd(num_train, block_size);
  yintblock = alloc_matd(num_train, block_size);
  if(use_parallel_blocks)
    block_terms = (double *)calloc((size_t)nblocks, sizeof(double));
  if((xblock == NULL) || (yintblock == NULL) ||
     (use_parallel_blocks && (block_terms == NULL)))
    goto cleanup_cdist_lp_block;

  *cv = 0.0;
  for(i0 = 0; i0 < num_train && !local_fail; i0 += block_size){
    const int ib = MIN(block_size, num_train - i0);
    const int block_id = i0 / block_size;
    double block_sum = 0.0;

    if(use_parallel_blocks && ((block_id % iNum_Processors) != my_rank))
      continue;

    if(np_conditional_x_weight_block_stream_core(vector_scale_factor, i0, ib, xblock) != 0){
      local_fail = 1;
      break;
    }

    for(j0 = 0; j0 < num_eval; j0 += block_size){
      const int jb = MIN(block_size, num_eval - j0);

      if(np_conditional_y_eval_block_stream_op_core(vector_scale_factor,
                                                    j0,
                                                    jb,
                                                    OP_INTEGRAL,
                                                    matrix_Y_unordered_eval_extern,
                                                    matrix_Y_ordered_eval_extern,
                                                    matrix_Y_continuous_eval_extern,
                                                    num_eval,
                                                    yintblock) != 0){
        local_fail = 1;
        break;
      }

      for(ii = 0; ii < ib; ii++){
        const int train_i = i0 + ii;
        double * const ai = xblock[ii];

        for(jj = 0; jj < jb; jj++){
          const int eval_j = j0 + jj;
          double fit;
          int indy;

          if(cdfontrain_extern && (train_i == eval_j))
            continue;

          fit = np_blas_ddot_int(num_train, ai, yintblock[jj]);
          indy = np_conditional_indicator_row_core(train_i,
                                                   eval_j,
                                                   cdfontrain_extern,
                                                   matrix_Y_ordered_train_extern,
                                                   matrix_Y_continuous_train_extern,
                                                   matrix_Y_ordered_eval_extern,
                                                   matrix_Y_continuous_eval_extern,
                                                   num_var_ordered_extern,
                                                   num_var_continuous_extern);
          {
            const double tvd = ((double)indy) - fit;
            block_sum += tvd*tvd;
          }
        }
      }
    }

    if(use_parallel_blocks){
      block_terms[block_id] = block_sum;
    } else {
      *cv += block_sum;
    }
  }

#ifdef MPI2
  if(use_parallel_blocks){
    int any_fail = 0;

    MPI_Allreduce(&local_fail, &any_fail, 1, MPI_INT, MPI_MAX, comm[1]);
    if(any_fail)
      goto cleanup_cdist_lp_block;

    MPI_Allreduce(MPI_IN_PLACE, block_terms, nblocks, MPI_DOUBLE, MPI_SUM, comm[1]);
    *cv = 0.0;
    for(ii = 0; ii < nblocks; ii++)
      *cv += block_terms[ii];
  }
#endif

  if(local_fail)
    goto cleanup_cdist_lp_block;

  *cv /= ((double)num_train*(double)MAX(1, num_eval));
  status = 0;

cleanup_cdist_lp_block:
  if(xblock != NULL) free_mat(xblock, block_size);
  if(yintblock != NULL) free_mat(yintblock, block_size);
  if(block_terms != NULL) free(block_terms);
  np_glp_cv_clear_extern();
  return status;
}

static int np_shadow_conditional_build_y_matrix(const int *operator_y,
                                                double *vector_scale_factor,
                                                double **matrix_Y_unordered_eval,
                                                double **matrix_Y_ordered_eval,
                                                double **matrix_Y_continuous_eval,
                                                const int num_eval,
                                                double *out_matrix){
  const int num_var_tot = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  const int bw_rows =
    (BANDWIDTH_den_extern == BW_FIXED) ? 1 :
    ((BANDWIDTH_den_extern == BW_GEN_NN) ? num_eval : num_obs_train_extern);
  int *kernel_cy = NULL, *kernel_uy = NULL, *kernel_oy = NULL;
  double *vsfy = NULL, *lambday = NULL, *kw = NULL;
  double **matrix_bandwidth_y = NULL, **matrix_bandwidth_eval_one = NULL;
  double **eval_yuno_one = NULL, **eval_yord_one = NULL, **eval_ycon_one = NULL;
  int i, l, status = 1;

  if((BANDWIDTH_den_extern != BW_FIXED) &&
     (BANDWIDTH_den_extern != BW_GEN_NN) &&
     (BANDWIDTH_den_extern != BW_ADAP_NN))
    return 1;

  vsfy = alloc_vecd(MAX(1, num_var_tot));
  lambday = alloc_vecd(MAX(1, num_var_unordered_extern + num_var_ordered_extern));
  kw = alloc_vecd(MAX(1, num_obs_train_extern));
  matrix_bandwidth_y = alloc_tmatd(bw_rows, num_var_continuous_extern);
  matrix_bandwidth_eval_one = alloc_tmatd(1, num_var_continuous_extern);
  if(num_var_unordered_extern > 0) eval_yuno_one = alloc_matd(1, num_var_unordered_extern);
  if(num_var_ordered_extern > 0) eval_yord_one = alloc_matd(1, num_var_ordered_extern);
  if(num_var_continuous_extern > 0) eval_ycon_one = alloc_matd(1, num_var_continuous_extern);
  kernel_cy = (int *)calloc((size_t)MAX(1, num_var_continuous_extern), sizeof(int));
  kernel_uy = (int *)calloc((size_t)MAX(1, num_var_unordered_extern), sizeof(int));
  kernel_oy = (int *)calloc((size_t)MAX(1, num_var_ordered_extern), sizeof(int));

  if((vsfy == NULL) || (lambday == NULL) || (kw == NULL) ||
     ((num_var_continuous_extern > 0) && (matrix_bandwidth_y == NULL)) ||
     ((num_var_continuous_extern > 0) && (matrix_bandwidth_eval_one == NULL)) ||
     ((num_var_unordered_extern > 0) && (eval_yuno_one == NULL)) ||
     ((num_var_ordered_extern > 0) && (eval_yord_one == NULL)) ||
     ((num_var_continuous_extern > 0) && (eval_ycon_one == NULL)) ||
     (kernel_cy == NULL) || (kernel_uy == NULL) || (kernel_oy == NULL))
    goto cleanup_ymat;

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern,
                        num_var_ordered_extern,
                        num_var_continuous_extern,
                        num_reg_unordered_extern,
                        num_reg_ordered_extern,
                        num_reg_continuous_extern,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        NULL,
                        vsfy,
                        NULL,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  for(i = 0; i < num_var_continuous_extern; i++) kernel_cy[i] = KERNEL_den_extern;
  for(i = 0; i < num_var_unordered_extern; i++) kernel_uy[i] = KERNEL_den_unordered_extern;
  for(i = 0; i < num_var_ordered_extern; i++) kernel_oy[i] = KERNEL_den_ordered_extern;

  if(kernel_bandwidth_mean(KERNEL_den_extern,
                           BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_eval,
                           0,
                           0,
                           0,
                           num_var_continuous_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           0,
                           vsfy,
                           NULL,
                           NULL,
                           matrix_Y_continuous_train_extern,
                           matrix_Y_continuous_eval,
                           NULL,
                           matrix_bandwidth_y,
                           lambday) == 1)
    goto cleanup_ymat;

  for(i = 0; i < num_eval; i++){
    const int orig_i = (int_TREE_Y == NP_TREE_TRUE) ? ipt_extern_Y[i] : i;
    for(l = 0; l < num_var_unordered_extern; l++)
      eval_yuno_one[l][0] = matrix_Y_unordered_eval[l][i];
    for(l = 0; l < num_var_ordered_extern; l++)
      eval_yord_one[l][0] = matrix_Y_ordered_eval[l][i];
    for(l = 0; l < num_var_continuous_extern; l++){
      eval_ycon_one[l][0] = matrix_Y_continuous_eval[l][i];
      matrix_bandwidth_eval_one[l][0] =
        (BANDWIDTH_den_extern == BW_GEN_NN) ? matrix_bandwidth_y[l][i] : matrix_bandwidth_y[l][0];
    }

    if(np_shadow_conditional_kernel_row(kernel_cy,
                                        kernel_uy,
                                        kernel_oy,
                                        operator_y,
                                        BANDWIDTH_den_extern,
                                        num_obs_train_extern,
                                        num_var_unordered_extern,
                                        num_var_ordered_extern,
                                        num_var_continuous_extern,
                                        matrix_Y_unordered_train_extern,
                                        matrix_Y_ordered_train_extern,
                                        matrix_Y_continuous_train_extern,
                                        eval_yuno_one,
                                        eval_yord_one,
                                        eval_ycon_one,
                                        vsfy,
                                        1,
                                        matrix_bandwidth_y,
                                        matrix_bandwidth_eval_one,
                                        lambday,
                                        num_categories_extern_Y,
                                        matrix_categorical_vals_extern_Y,
                                        int_TREE_Y,
                                        kdt_extern_Y,
                                        kw,
                                        NULL) != 0)
      goto cleanup_ymat;

    for(l = 0; l < num_obs_train_extern; l++){
      const int orig_l = (int_TREE_Y == NP_TREE_TRUE) ? ipt_extern_Y[l] : l;
      out_matrix[(size_t)orig_i*(size_t)num_obs_train_extern + (size_t)orig_l] = kw[l];
    }
  }

  status = 0;

cleanup_ymat:
  if(vsfy != NULL) free(vsfy);
  if(lambday != NULL) free(lambday);
  if(kw != NULL) free(kw);
  if(matrix_bandwidth_y != NULL) free_tmat(matrix_bandwidth_y);
  if(matrix_bandwidth_eval_one != NULL) free_tmat(matrix_bandwidth_eval_one);
  if(eval_yuno_one != NULL) free_mat(eval_yuno_one, num_var_unordered_extern);
  if(eval_yord_one != NULL) free_mat(eval_yord_one, num_var_ordered_extern);
  if(eval_ycon_one != NULL) free_mat(eval_ycon_one, num_var_continuous_extern);
  if(kernel_cy != NULL) free(kernel_cy);
  if(kernel_uy != NULL) free(kernel_uy);
  if(kernel_oy != NULL) free(kernel_oy);
  return status;
}

int np_shadow_cv_con_density_ml(double *vector_scale_factor, double *cv){
  double *weights = NULL, *yrow = NULL;
  int *operator_y = NULL;
  int i, j, status = 1;

  if((cv == NULL) || (vector_scale_factor == NULL))
    return 1;

  weights = (double *)malloc((size_t)num_obs_train_extern*(size_t)num_obs_train_extern*sizeof(double));
  yrow = (double *)malloc((size_t)num_obs_train_extern*(size_t)num_obs_train_extern*sizeof(double));
  operator_y = (int *)calloc((size_t)MAX(1, num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern), sizeof(int));
  if((weights == NULL) || (yrow == NULL) || (operator_y == NULL))
    goto cleanup_cvml_shadow;

  for(i = 0; i < (num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern); i++)
    operator_y[i] = OP_NORMAL;
  if(np_shadow_conditional_build_y_matrix(operator_y,
                                         vector_scale_factor,
                                         matrix_Y_unordered_train_extern,
                                         matrix_Y_ordered_train_extern,
                                         matrix_Y_continuous_train_extern,
                                         num_obs_train_extern,
                                         yrow) != 0)
    goto cleanup_cvml_shadow;

  *cv = 0.0;
  if(np_shadow_conditional_build_x_weights(vector_scale_factor, weights) != 0)
    goto cleanup_cvml_shadow;

  for(i = 0; i < num_obs_train_extern; i++){
    double fit = 0.0;
    for(j = 0; j < num_obs_train_extern; j++)
      fit += weights[(size_t)i*(size_t)num_obs_train_extern + (size_t)j] *
        yrow[(size_t)i*(size_t)num_obs_train_extern + (size_t)j];
    *cv -= log((fit > DBL_MIN) ? fit : DBL_MIN);
  }

  status = 0;

cleanup_cvml_shadow:
  if(weights != NULL) free(weights);
  if(yrow != NULL) free(yrow);
  if(operator_y != NULL) free(operator_y);
  return status;
}

int np_shadow_cv_con_density_ls(double *vector_scale_factor, double *cv){
  double *weights = NULL, *weights_full = NULL, *yrow = NULL, *yconv = NULL;
  int *operator_y = NULL;
  int i, j, k, status = 1;

  if((cv == NULL) || (vector_scale_factor == NULL))
    return 1;

  weights = (double *)malloc((size_t)num_obs_train_extern*(size_t)num_obs_train_extern*sizeof(double));
  weights_full = (double *)malloc((size_t)num_obs_train_extern*(size_t)num_obs_train_extern*sizeof(double));
  yrow = (double *)malloc((size_t)num_obs_train_extern*(size_t)num_obs_train_extern*sizeof(double));
  yconv = (double *)malloc((size_t)num_obs_train_extern*(size_t)num_obs_train_extern*sizeof(double));
  operator_y = (int *)calloc((size_t)MAX(1, num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern), sizeof(int));
  if((weights == NULL) || (weights_full == NULL) || (yrow == NULL) || (yconv == NULL) || (operator_y == NULL))
    goto cleanup_cvls_shadow;

  if(np_shadow_conditional_build_x_weights(vector_scale_factor, weights) != 0)
    goto cleanup_cvls_shadow;
  if(np_shadow_conditional_build_x_weights_full(vector_scale_factor, weights_full) != 0)
    goto cleanup_cvls_shadow;

  for(i = 0; i < (num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern); i++)
    operator_y[i] = OP_NORMAL;
  if(np_shadow_conditional_build_y_matrix(operator_y,
                                         vector_scale_factor,
                                         matrix_Y_unordered_train_extern,
                                         matrix_Y_ordered_train_extern,
                                         matrix_Y_continuous_train_extern,
                                         num_obs_train_extern,
                                         yrow) != 0)
    goto cleanup_cvls_shadow;

  for(i = 0; i < (num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern); i++)
    operator_y[i] = OP_CONVOLUTION;
  if(np_shadow_conditional_build_y_matrix(operator_y,
                                         vector_scale_factor,
                                         matrix_Y_unordered_train_extern,
                                         matrix_Y_ordered_train_extern,
                                         matrix_Y_continuous_train_extern,
                                         num_obs_train_extern,
                                         yconv) != 0)
    goto cleanup_cvls_shadow;

  *cv = 0.0;
  for(i = 0; i < num_obs_train_extern; i++){
    double quad = 0.0;
    double lin = 0.0;
    for(j = 0; j < num_obs_train_extern; j++){
      const double aij = weights[(size_t)i*(size_t)num_obs_train_extern + (size_t)j];
      const double aij_full = weights_full[(size_t)i*(size_t)num_obs_train_extern + (size_t)j];
      lin += aij*yrow[(size_t)i*(size_t)num_obs_train_extern + (size_t)j];
      if(aij_full == 0.0)
        continue;
      for(k = 0; k < num_obs_train_extern; k++)
        quad += aij_full*
          weights_full[(size_t)i*(size_t)num_obs_train_extern + (size_t)k]*
          yconv[(size_t)j*(size_t)num_obs_train_extern + (size_t)k];
    }
    *cv += quad - 2.0*lin;
  }
  *cv /= (double)num_obs_train_extern;
  status = 0;

cleanup_cvls_shadow:
  if(weights != NULL) free(weights);
  if(weights_full != NULL) free(weights_full);
  if(yrow != NULL) free(yrow);
  if(yconv != NULL) free(yconv);
  if(operator_y != NULL) free(operator_y);
  return status;
}

int np_shadow_cv_con_distribution_ls(double *vector_scale_factor, double *cv){
  const int num_eval = num_obs_eval_extern;
  double *weights = NULL, *ycdf = NULL;
  int *operator_y = NULL;
  int i, j, l, status = 1;

  if((cv == NULL) || (vector_scale_factor == NULL))
    return 1;

  weights = (double *)malloc((size_t)num_obs_train_extern*(size_t)num_obs_train_extern*sizeof(double));
  ycdf = (double *)malloc((size_t)num_eval*(size_t)num_obs_train_extern*sizeof(double));
  operator_y = (int *)calloc((size_t)MAX(1, num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern), sizeof(int));
  if((weights == NULL) || (ycdf == NULL) || (operator_y == NULL))
    goto cleanup_cdist_shadow;

  if(np_shadow_conditional_build_x_weights(vector_scale_factor, weights) != 0)
    goto cleanup_cdist_shadow;

  for(i = 0; i < (num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern); i++)
    operator_y[i] = OP_INTEGRAL;

  if(np_shadow_conditional_build_y_matrix(operator_y,
                                         vector_scale_factor,
                                         matrix_Y_unordered_eval_extern,
                                         matrix_Y_ordered_eval_extern,
                                         matrix_Y_continuous_eval_extern,
                                         num_eval,
                                         ycdf) != 0)
    goto cleanup_cdist_shadow;

  *cv = 0.0;
  for(i = 0; i < num_obs_train_extern; i++){
    for(j = 0; j < num_eval; j++){
      double fit = 0.0;
      const int indy = np_shadow_conditional_indicator_row(i,
                                                           j,
                                                           cdfontrain_extern,
                                                           matrix_Y_ordered_train_extern,
                                                           matrix_Y_continuous_train_extern,
                                                           matrix_Y_ordered_eval_extern,
                                                           matrix_Y_continuous_eval_extern,
                                                           num_var_ordered_extern,
                                                           num_var_continuous_extern);
      if(cdfontrain_extern && (i == j))
        continue;
      for(l = 0; l < num_obs_train_extern; l++)
        fit += weights[(size_t)i*(size_t)num_obs_train_extern + (size_t)l]*
          ycdf[(size_t)j*(size_t)num_obs_train_extern + (size_t)l];
      {
        const double tvd = ((double)indy) - fit;
        *cv += tvd*tvd;
      }
    }
  }
  *cv /= ((double)num_obs_train_extern*(double)MAX(1, num_eval));
  status = 0;

cleanup_cdist_shadow:
  if(weights != NULL) free(weights);
  if(ycdf != NULL) free(ycdf);
  if(operator_y != NULL) free(operator_y);
  return status;
}


int np_kernel_estimate_density_categorical_leave_one_out_cv(int KERNEL_den,
                                                            int KERNEL_unordered_den,
                                                            int KERNEL_ordered_den,
                                                            int BANDWIDTH_den,
                                                            int num_obs,
                                                            int num_reg_unordered,
                                                            int num_reg_ordered,
                                                            int num_reg_continuous,
                                                            double **matrix_X_unordered,
                                                            double **matrix_X_ordered,
double **matrix_X_continuous,
double *vector_scale_factor,
int *num_categories,
double *cv){
  NP_GateOverrideCtx gate_ctx_local;
  np_gate_ctx_clear(&gate_ctx_local);

  const int num_reg = num_reg_continuous+num_reg_unordered+num_reg_ordered;
  const int bwmdim = (BANDWIDTH_den==BW_GEN_NN)?num_obs:
    ((BANDWIDTH_den==BW_ADAP_NN)?num_obs:1);

  int i;

  int * operator = NULL;
  int gate_override_active = 0;
  int all_large_gate = 0;
  int *ov_cont_ok = NULL, *ov_disc_uno_ok = NULL, *ov_disc_ord_ok = NULL;
  double *ov_cont_hmin = NULL, *ov_cont_k0 = NULL;
  double *ov_disc_uno_const = NULL, *ov_disc_ord_const = NULL;
  int ov_cont_from_cache = 0;
  double **matrix_bandwidth = NULL;
  double *lambda = NULL;

  int num_obs_alloc;

#ifdef MPI2
  int stride_t = MAX((int)ceil((double) num_obs / (double) iNum_Processors),1);
  
  num_obs_alloc = stride_t*iNum_Processors;
#else
  num_obs_alloc = num_obs;
#endif

	const double log_DBL_MIN = log(DBL_MIN);
  double * rho = (double *)malloc(num_obs_alloc*sizeof(double));
  
  if(rho == NULL)
    error("failed to allocate rho");

  operator = (int *)malloc(sizeof(int)*num_reg);

  for(i = 0; i < num_reg; i++)
    operator[i] = OP_NORMAL;

  int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;

  kernel_c = (int *)malloc(sizeof(int)*num_reg_continuous);

  for(i = 0; i < num_reg_continuous; i++)
    kernel_c[i] = KERNEL_den;

  kernel_u = (int *)malloc(sizeof(int)*num_reg_unordered);

  for(i = 0; i < num_reg_unordered; i++)
    kernel_u[i] = KERNEL_unordered_den;

  kernel_o = (int *)malloc(sizeof(int)*num_reg_ordered);

  for(i = 0; i < num_reg_ordered; i++)
    kernel_o[i] = KERNEL_ordered_den;

  matrix_bandwidth = alloc_matd(bwmdim,num_reg_continuous);
  lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

  if(kernel_bandwidth_mean(KERNEL_den,
                           BANDWIDTH_den,
                           num_obs,
                           num_obs,
                           0,0,0,
                           num_reg_continuous,
                           num_reg_unordered,
                           num_reg_ordered,
                           0,
                           vector_scale_factor,
                           NULL,NULL,
                           matrix_X_continuous,
                           matrix_X_continuous,
                           NULL,
                           matrix_bandwidth,
                           lambda)==1){
    error("\n** Error: invalid bandwidth.");
  }

  if((num_reg_continuous + num_reg_unordered + num_reg_ordered) > 0){
    int ok_all = 1;

    if(num_reg_continuous > 0){
      const double rel_tol = np_cont_largeh_rel_tol();
      if(np_cont_largeh_cache_get_or_build(num_obs,
                                           num_obs,
                                           num_reg_continuous,
                                           kernel_c,
                                           matrix_X_continuous,
                                           matrix_X_continuous,
                                           rel_tol,
                                           &ov_cont_ok,
                                           &ov_cont_hmin,
                                           &ov_cont_k0)){
        ov_cont_from_cache = 1;
      } else {
        ov_cont_ok = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
        ov_cont_hmin = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
        ov_cont_k0 = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
        ok_all = (ov_cont_ok != NULL) && (ov_cont_hmin != NULL) && (ov_cont_k0 != NULL);
        if(ok_all){
          for(i = 0; i < num_reg_continuous; i++){
            const int kern = kernel_c[i];
            double xmin = DBL_MAX, xmax = -DBL_MAX;
            ov_cont_ok[i] = 0; ov_cont_hmin[i] = DBL_MAX; ov_cont_k0[i] = 0.0;
            if(!np_cont_largeh_kernel_supported(kern)) continue;
            for(int j = 0; j < num_obs; j++){
              const double v = matrix_X_continuous[i][j];
              if(!isfinite(v)) continue;
              xmin = MIN(xmin, v); xmax = MAX(xmax, v);
            }
            if(xmax >= xmin){
              const double utol = np_cont_largeh_utol(kern, rel_tol);
              if(utol > 0.0 && isfinite(utol)){
                ov_cont_ok[i] = 1;
                ov_cont_hmin[i] = (xmax - xmin)/utol;
                ov_cont_k0[i] = np_cont_largeh_k0(kern);
              }
            }
          }
        }
      }
    }

    if(ok_all && num_reg_unordered > 0){
      ov_disc_uno_ok = (int *)calloc((size_t)num_reg_unordered, sizeof(int));
      ov_disc_uno_const = (double *)malloc((size_t)num_reg_unordered*sizeof(double));
      ok_all = (ov_disc_uno_ok != NULL) && (ov_disc_uno_const != NULL);
      if(ok_all){
        double (* const ukf[])(int, double, int) = {
          np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
          np_score_uaa, np_score_unli_racine
        };
        const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
        for(i = 0; i < num_reg_unordered; i++){
          const int ku = kernel_u[i];
          const int ncat = (num_categories != NULL) ? num_categories[i] : 0;
          const double lam = lambda[i];
          ov_disc_uno_ok[i] = 0; ov_disc_uno_const[i] = 0.0;
          if(ku < 0 || ku >= nuk) continue;
          if(!np_disc_near_upper(ku, lam, ncat)) continue;
          {
            const double ks = ukf[ku](1, lam, ncat);
            const double kd = ukf[ku](0, lam, ncat);
            if(np_disc_near_const_kernel(ks, kd)){
              ov_disc_uno_ok[i] = 1;
              ov_disc_uno_const[i] = 0.5*(ks + kd);
            }
          }
        }
      }
    }

    if(ok_all && num_reg_ordered > 0){
      ov_disc_ord_ok = (int *)calloc((size_t)num_reg_ordered, sizeof(int));
      ov_disc_ord_const = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
      ok_all = (ov_disc_ord_ok != NULL) && (ov_disc_ord_const != NULL);
      if(ok_all){
        double (* const okf[])(double, double, double, double, double) = {
          np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
        };
        const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
        for(i = 0; i < num_reg_ordered; i++){
          const int oi = i + num_reg_unordered;
          const int ko = kernel_o[i];
          const int ncat = (num_categories != NULL) ? num_categories[oi] : 0;
          const double lam = lambda[oi];
          ov_disc_ord_ok[i] = 0; ov_disc_ord_const[i] = 0.0;
          if(ko < 0 || ko >= nok) continue;
          if(ncat <= 0 || matrix_categorical_vals_extern == NULL) continue;
          if(!np_disc_ordered_near_upper(ko, lam)) continue;
          {
            const double cl = matrix_categorical_vals_extern[oi][0];
            const double ch = matrix_categorical_vals_extern[oi][ncat - 1];
            const double k0 = okf[ko](cl, cl, lam, cl, ch);
            const double k1 = okf[ko](cl, ch, lam, cl, ch);
            if(np_disc_near_const_kernel(k0, k1)){
              ov_disc_ord_ok[i] = 1;
              ov_disc_ord_const[i] = 0.5*(k0 + k1);
            }
          }
        }
      }
    }

    if(ok_all){
      np_gate_ctx_set(&gate_ctx_local,
                      num_reg_continuous,
                      num_reg_unordered,
                      num_reg_ordered,
                      kernel_c,
                      kernel_u,
                      kernel_o,
                      operator,
                      ov_cont_ok,
                      ov_cont_hmin,
                      ov_cont_k0,
                      ov_disc_uno_ok,
                      ov_disc_uno_const,
                      ov_disc_ord_ok,
                      ov_disc_ord_const);
      gate_override_active = 1;
    }
  }

  all_large_gate = (BANDWIDTH_den == BW_FIXED) && gate_override_active;
  if(all_large_gate){
    for(i = 0; i < num_reg_continuous; i++){
      const double bw = matrix_bandwidth[i][0];
      if((ov_cont_ok == NULL) || (!ov_cont_ok[i]) || (!isfinite(bw)) ||
         (bw <= 0.0) || (bw < ov_cont_hmin[i])){
        all_large_gate = 0;
        break;
      }
    }
  }
  if(all_large_gate){
    for(i = 0; i < num_reg_unordered; i++){
      if((ov_disc_uno_ok == NULL) || (!ov_disc_uno_ok[i])){
        all_large_gate = 0;
        break;
      }
    }
  }
  if(all_large_gate){
    for(i = 0; i < num_reg_ordered; i++){
      if((ov_disc_ord_ok == NULL) || (!ov_disc_ord_ok[i])){
        all_large_gate = 0;
        break;
      }
    }
  }
  if(all_large_gate)
    np_fastcv_alllarge_hits++;

  kernel_weighted_sum_np_ctx(kernel_c,
                         kernel_u,
                         kernel_o,
                         BANDWIDTH_den,
                         num_obs,
                         num_obs,
                         num_reg_unordered,
                         num_reg_ordered,
                         num_reg_continuous,
                         1, // leave one out 
                         0,
                         1, // kernel_pow = 1
                         1, // bandwidth_divide = FALSE when not adaptive
                         0, 
                         0, // symmetric
                         0, // gather-scatter sum
                         0, // do not drop train
                         0, // do not drop train
                         operator, // no convolution
                         OP_NOOP, // no permutations
                         0, // no score
                         0, // no ocg
                         NULL,
                         0, //  explicity suppress parallel
                         0,
                         0,
                         int_TREE_X,
                         0,
                         kdt_extern_X,
                         NULL, NULL, NULL,
                         matrix_X_unordered,
                         matrix_X_ordered,
                         matrix_X_continuous,
                         matrix_X_unordered,
                         matrix_X_ordered,
                         matrix_X_continuous,
                         NULL,
                         NULL,
                         NULL,
                         vector_scale_factor,
                         1,matrix_bandwidth,matrix_bandwidth,lambda,
                         num_categories,
                         NULL,
                         NULL,
                         rho,  // weighted sum
                         NULL, // no permutations
                         NULL, // do not return kernel weights
                         &gate_ctx_local);
  
  

  for(i = 0, *cv = 0.0; i < num_obs; i++)
    *cv -= (rho[i] < DBL_MIN) ? log_DBL_MIN : log(rho[i]/(num_obs-1.0));
    

  free(operator);
  free(kernel_c);
  free(kernel_u);
  free(kernel_o);
  free(rho);
  free(lambda);
  free_mat(matrix_bandwidth, num_reg_continuous);
  if(gate_override_active) np_gate_ctx_clear(&gate_ctx_local);
  if((ov_cont_ok != NULL) && (!ov_cont_from_cache)) free(ov_cont_ok);
  if((ov_cont_hmin != NULL) && (!ov_cont_from_cache)) free(ov_cont_hmin);
  if((ov_cont_k0 != NULL) && (!ov_cont_from_cache)) free(ov_cont_k0);
  if(ov_disc_uno_ok != NULL) free(ov_disc_uno_ok);
  if(ov_disc_uno_const != NULL) free(ov_disc_uno_const);
  if(ov_disc_ord_ok != NULL) free(ov_disc_ord_ok);
  if(ov_disc_ord_const != NULL) free(ov_disc_ord_const);
  return(0);
}

int np_kernel_estimate_density_categorical_convolution_cv(int KERNEL_den,
                                                          int KERNEL_unordered_den,
                                                          int KERNEL_ordered_den,
                                                          int BANDWIDTH_den,
                                                          int num_obs,
                                                          int num_reg_unordered,
                                                          int num_reg_ordered,
                                                          int num_reg_continuous,
                                                          double **matrix_X_unordered,
                                                          double **matrix_X_ordered,
                                                          double **matrix_X_continuous,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *cv){
  NP_GateOverrideCtx gate_x_ctx, gate_y_ctx, gate_xy_ctx;
  np_gate_ctx_clear(&gate_x_ctx);
  np_gate_ctx_clear(&gate_y_ctx);
  np_gate_ctx_clear(&gate_xy_ctx);

  const int num_reg = num_reg_continuous+num_reg_unordered+num_reg_ordered;
  const int bwmdim = (BANDWIDTH_den==BW_GEN_NN)?num_obs:
    ((BANDWIDTH_den==BW_ADAP_NN)?num_obs:1);
  const int bounded_scalar_quadrature =
    np_density_cvls_bounded_scalar_route_ok(num_reg_unordered,
                                            num_reg_ordered,
                                            num_reg_continuous);
  const int bounded_general_quadrature =
    np_density_cvls_bounded_general_route_ok(num_reg_continuous);

  int i;
  int status = 1;

  int * operator = NULL;
  int gate_override_active = 0;
  int *ov_cont_ok = NULL, *ov_disc_uno_ok = NULL, *ov_disc_ord_ok = NULL;
  double *ov_cont_hmin = NULL, *ov_cont_k0 = NULL;
  double *ov_disc_uno_const = NULL, *ov_disc_ord_const = NULL;
  double **matrix_bandwidth = NULL;
  double *lambda = NULL;

  int num_obs_alloc;

#ifdef MPI2
  int stride_t = MAX((int)ceil((double) num_obs / (double) iNum_Processors),1);
  
  num_obs_alloc = stride_t*iNum_Processors;
#else
  num_obs_alloc = num_obs;
#endif

  double * res = (double *)np_jksum_malloc_array_or_die((size_t)num_obs_alloc, sizeof(double), "np_kernel_estimate_density_categorical_convolution_cv rho");
  double cv1, cv2;

  operator = (int *)np_jksum_malloc_array_or_die((size_t)num_reg, sizeof(int), "np_kernel_estimate_density_categorical_convolution_cv operator");

  // first the convolution portion
  for(i = 0; i < num_reg; i++)
    operator[i] = OP_CONVOLUTION;

  int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;

  kernel_c = (int *)np_jksum_malloc_array_or_die((size_t)num_reg_continuous, sizeof(int), "np_kernel_estimate_density_categorical_convolution_cv kernel_c");

  for(i = 0; i < num_reg_continuous; i++)
    kernel_c[i] = KERNEL_den;

  kernel_u = (int *)np_jksum_malloc_array_or_die((size_t)num_reg_unordered, sizeof(int), "np_kernel_estimate_density_categorical_convolution_cv kernel_u");

  for(i = 0; i < num_reg_unordered; i++)
    kernel_u[i] = KERNEL_unordered_den;

  kernel_o = (int *)np_jksum_malloc_array_or_die((size_t)num_reg_ordered, sizeof(int), "np_kernel_estimate_density_categorical_convolution_cv kernel_o");

  for(i = 0; i < num_reg_ordered; i++)
    kernel_o[i] = KERNEL_ordered_den;

  matrix_bandwidth = alloc_matd(bwmdim,num_reg_continuous);
  lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

  if(kernel_bandwidth_mean(KERNEL_den,
                           BANDWIDTH_den,
                           num_obs,
                           num_obs,
                           0,0,0,
                           num_reg_continuous,
                           num_reg_unordered,
                           num_reg_ordered,
                           0,
                           vector_scale_factor,
                           NULL,NULL,
                           matrix_X_continuous,
                           matrix_X_continuous,
                           NULL,
                           matrix_bandwidth,
                           lambda)==1){
    error("\n** Error: invalid bandwidth.");
  }

  if(bounded_scalar_quadrature){
    if(np_density_cvls_bounded_i1_quadrature(KERNEL_den,
                                             BANDWIDTH_den,
                                             num_obs,
                                             num_reg_unordered,
                                             num_reg_ordered,
                                             num_reg_continuous,
                                             matrix_X_continuous,
                                             vector_scale_factor,
                                             &cv1) != 0)
      goto cleanup_density_convolution_cv;
  } else if(bounded_general_quadrature) {
    if(np_density_cvls_bounded_i1_quadrature_general(KERNEL_den,
                                                     KERNEL_unordered_den,
                                                     KERNEL_ordered_den,
                                                     BANDWIDTH_den,
                                                     num_obs,
                                                     num_reg_unordered,
                                                     num_reg_ordered,
                                                     num_reg_continuous,
                                                     matrix_X_unordered,
                                                     matrix_X_ordered,
                                                     matrix_X_continuous,
                                                     vector_scale_factor,
                                                     num_categories,
                                                     matrix_categorical_vals,
                                                     &cv1) != 0)
      goto cleanup_density_convolution_cv;
  } else {
    if(int_cker_bound_extern){
      error("bounded npudens cv.ls currently supports up to two continuous variables");
    }

    kernel_weighted_sum_np_ctx(kernel_c,
                           kernel_u,
                           kernel_o,
                           BANDWIDTH_den,
                           num_obs,
                           num_obs,
                           num_reg_unordered,
                           num_reg_ordered,
                           num_reg_continuous,
                           0, // don't leave one out
                           0,
                           1, // kernel_pow = 1
                           1, // bandwidth_divide = FALSE when not adaptive
                           0,
                           0, // symmetric
                           0, // gather-scatter sum
                           0, // do not drop train
                           0, // do not drop train
                           operator, // convolution
                           OP_NOOP, // no permutations
                           0, // no score
                           0, // no ocg
                           NULL,
                           0, //  explicity suppress parallel
                           0,
                           0,
                           int_TREE_X,
                           0,
                           kdt_extern_X,
                           NULL, NULL, NULL,
                           matrix_X_unordered,
                           matrix_X_ordered,
                           matrix_X_continuous,
                           matrix_X_unordered,
                           matrix_X_ordered,
                           matrix_X_continuous,
                           NULL, // no ys
                           NULL, // no weights
                           NULL, // no sgn
                           vector_scale_factor,
                           1,matrix_bandwidth,matrix_bandwidth,lambda,
                           num_categories,
                           matrix_categorical_vals,
                           NULL, // no ocg
                           res,  // weighted sum
                           NULL, // no permutations
                           NULL, // do not return kernel weights
                           &gate_x_ctx);

    for(i = 0, cv1 = 0.0; i < num_obs; i++) cv1 += res[i];

    cv1 /= num_obs*num_obs;
  }

  // then the cross term
  for(i = 0; i < num_reg; i++)
    operator[i] = OP_NORMAL;

  if((num_reg_continuous + num_reg_unordered + num_reg_ordered) > 0){
    int ok_all = 1;

    if(num_reg_continuous > 0){
      ov_cont_ok = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
      ov_cont_hmin = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
      ov_cont_k0 = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
      ok_all = (ov_cont_ok != NULL) && (ov_cont_hmin != NULL) && (ov_cont_k0 != NULL);
      if(ok_all){
        const double rel_tol = np_cont_largeh_rel_tol();
        for(i = 0; i < num_reg_continuous; i++){
          const int kern = kernel_c[i];
          double xmin = DBL_MAX, xmax = -DBL_MAX;
          ov_cont_ok[i] = 0; ov_cont_hmin[i] = DBL_MAX; ov_cont_k0[i] = 0.0;
          if(!np_cont_largeh_kernel_supported(kern)) continue;
          for(int j = 0; j < num_obs; j++){
            const double v = matrix_X_continuous[i][j];
            if(!isfinite(v)) continue;
            xmin = MIN(xmin, v); xmax = MAX(xmax, v);
          }
          if(xmax >= xmin){
            const double utol = np_cont_largeh_utol(kern, rel_tol);
            if(utol > 0.0 && isfinite(utol)){
              ov_cont_ok[i] = 1;
              ov_cont_hmin[i] = (xmax - xmin)/utol;
              ov_cont_k0[i] = np_cont_largeh_k0(kern);
            }
          }
        }
      }
    }

    if(ok_all && num_reg_unordered > 0){
      ov_disc_uno_ok = (int *)calloc((size_t)num_reg_unordered, sizeof(int));
      ov_disc_uno_const = (double *)malloc((size_t)num_reg_unordered*sizeof(double));
      ok_all = (ov_disc_uno_ok != NULL) && (ov_disc_uno_const != NULL);
      if(ok_all){
        double (* const ukf[])(int, double, int) = {
          np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
          np_score_uaa, np_score_unli_racine
        };
        const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
        for(i = 0; i < num_reg_unordered; i++){
          const int ku = kernel_u[i];
          const int ncat = (num_categories != NULL) ? num_categories[i] : 0;
          const double lam = lambda[i];
          ov_disc_uno_ok[i] = 0; ov_disc_uno_const[i] = 0.0;
          if(ku < 0 || ku >= nuk) continue;
          if(!np_disc_near_upper(ku, lam, ncat)) continue;
          {
            const double ks = ukf[ku](1, lam, ncat);
            const double kd = ukf[ku](0, lam, ncat);
            if(np_disc_near_const_kernel(ks, kd)){
              ov_disc_uno_ok[i] = 1;
              ov_disc_uno_const[i] = 0.5*(ks + kd);
            }
          }
        }
      }
    }

    if(ok_all && num_reg_ordered > 0){
      ov_disc_ord_ok = (int *)calloc((size_t)num_reg_ordered, sizeof(int));
      ov_disc_ord_const = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
      ok_all = (ov_disc_ord_ok != NULL) && (ov_disc_ord_const != NULL);
      if(ok_all){
        double (* const okf[])(double, double, double, double, double) = {
          np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
        np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
        np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
        np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
        };
        const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
        for(i = 0; i < num_reg_ordered; i++){
          const int oi = i + num_reg_unordered;
          const int ko = kernel_o[i];
          const int ncat = (num_categories != NULL) ? num_categories[oi] : 0;
          const double lam = lambda[oi];
          ov_disc_ord_ok[i] = 0; ov_disc_ord_const[i] = 0.0;
          if(ko < 0 || ko >= nok) continue;
          if(ncat <= 0 || matrix_categorical_vals_extern == NULL) continue;
          if(!np_disc_ordered_near_upper(ko, lam)) continue;
          {
            const double cl = matrix_categorical_vals_extern[oi][0];
            const double ch = matrix_categorical_vals_extern[oi][ncat - 1];
            const double k0 = okf[ko](cl, cl, lam, cl, ch);
            const double k1 = okf[ko](cl, ch, lam, cl, ch);
            if(np_disc_near_const_kernel(k0, k1)){
              ov_disc_ord_ok[i] = 1;
              ov_disc_ord_const[i] = 0.5*(k0 + k1);
            }
          }
        }
      }
    }

    if(ok_all){
      np_gate_ctx_set(&gate_x_ctx,
                      num_reg_continuous,
                      num_reg_unordered,
                      num_reg_ordered,
                      kernel_c,
                      kernel_u,
                      kernel_o,
                      operator,
                      ov_cont_ok,
                      ov_cont_hmin,
                      ov_cont_k0,
                      ov_disc_uno_ok,
                      ov_disc_uno_const,
                      ov_disc_ord_ok,
                      ov_disc_ord_const);
      gate_override_active = 1;
    } else {
      np_gate_ctx_clear(&gate_x_ctx);
    }
  }

  kernel_weighted_sum_np_ctx(kernel_c,
                         kernel_u,
                         kernel_o,
                         BANDWIDTH_den,
                         num_obs,
                         num_obs,
                         num_reg_unordered,
                         num_reg_ordered,
                         num_reg_continuous,
                         1, // leave one out 
                         0,
                         1, // kernel_pow = 1
                         1, // bandwidth_divide = FALSE when not adaptive
                         0, 
                         0, // symmetric
                         0, // gather-scatter sum
                         0, // do not drop train
                         0, // do not drop train
                         operator, // convolution
                         OP_NOOP, // no permutations
                         0, // no score
                         0, // no ocg
                         NULL,
                         0, //  explicity suppress parallel
                         0,
                         0,
                         int_TREE_X,
                         0,
                         kdt_extern_X,
                         NULL, NULL, NULL,
                         matrix_X_unordered,
                         matrix_X_ordered,
                         matrix_X_continuous,
                         matrix_X_unordered,
                         matrix_X_ordered,
                         matrix_X_continuous,
                         NULL, // no ys
                         NULL, // no weights
                         NULL, // no sgn
                         vector_scale_factor,
                         1,matrix_bandwidth,matrix_bandwidth,lambda,
                         num_categories,
                         NULL,
                         NULL, // no ocg
                         res,  // weighted sum
                         NULL, // no permutations
                         NULL, // do not return kernel weights
                         &gate_x_ctx);


  for(i = 0, cv2 = 0.0; i < num_obs; i++) cv2 += res[i];

  cv2 /= num_obs*(num_obs-1.0);

  *cv = cv1 - 2.0*cv2;
  status = 0;

cleanup_density_convolution_cv:
  free(operator);
  free(kernel_c);
  free(kernel_u);
  free(kernel_o);
  free(res);
  free(lambda);
  free_mat(matrix_bandwidth, num_reg_continuous);
  if(gate_override_active) np_gate_ctx_clear(&gate_x_ctx);
  if(ov_cont_ok != NULL) free(ov_cont_ok);
  if(ov_cont_hmin != NULL) free(ov_cont_hmin);
  if(ov_cont_k0 != NULL) free(ov_cont_k0);
  if(ov_disc_uno_ok != NULL) free(ov_disc_uno_ok);
  if(ov_disc_uno_const != NULL) free(ov_disc_uno_const);
  if(ov_disc_ord_ok != NULL) free(ov_disc_ord_ok);
  if(ov_disc_ord_const != NULL) free(ov_disc_ord_const);
  return(status);

}

void kernel_estimate_dens_dist_categorical_np(int KERNEL_den,
                                              int KERNEL_unordered_den,
                                              int KERNEL_ordered_den,
                                              int BANDWIDTH_den,
                                              int num_obs_train,
                                              int num_obs_eval,
                                              int num_reg_unordered,
                                              int num_reg_ordered,
                                              int num_reg_continuous,
                                              int dop,
                                              double **matrix_X_unordered_train,
                                              double **matrix_X_ordered_train,
                                              double **matrix_X_continuous_train,
                                              double **matrix_X_unordered_eval,
                                              double **matrix_X_ordered_eval,
                                              double **matrix_X_continuous_eval,
                                              double *vector_scale_factor,
                                              int *num_categories,
                                              double ** matrix_categorical_vals,
                                              double *pdf,
                                              double *pdf_stderr,
                                              double *log_likelihood){
  NP_GateOverrideCtx gate_ctx_local;
  np_gate_ctx_clear(&gate_ctx_local);

  const int num_reg = num_reg_continuous+num_reg_unordered+num_reg_ordered;

  int i, l;

  const int bwmdim = (BANDWIDTH_den==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_den==BW_ADAP_NN)?num_obs_train:1);

  int * operator = NULL;
  int gate_override_active = 0;
  int *ov_cont_ok = NULL, *ov_disc_uno_ok = NULL, *ov_disc_ord_ok = NULL;
  double *ov_cont_hmin = NULL, *ov_cont_k0 = NULL;
  double *ov_disc_uno_const = NULL, *ov_disc_ord_const = NULL;

	double INT_KERNEL_P;					 /* Integral of K(z)^2 */
	double K_INT_KERNEL_P;				 /*  K^p */
	/* Integral of K(z-0.5)*K(z+0.5) */
	double INT_KERNEL_PM_HALF = 0.0;
	double DIFF_KER_PPM = 0.0;		 /* Difference between int K(z)^p and int K(z-.5)K(z+.5) */

  double ** matrix_bandwidth = NULL, * lambda = NULL;

  double pnh = (double)num_obs_train;

	const double log_DBL_MIN = log(DBL_MIN);

  operator = (int *)malloc(sizeof(int)*num_reg);

  for(i = 0; i < num_reg; i++)
    operator[i] = dop;

  int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;

  kernel_c = (int *)malloc(sizeof(int)*num_reg_continuous);

  for(i = 0; i < num_reg_continuous; i++)
    kernel_c[i] = KERNEL_den;

  kernel_u = (int *)malloc(sizeof(int)*num_reg_unordered);

  for(i = 0; i < num_reg_unordered; i++)
    kernel_u[i] = KERNEL_unordered_den;

  kernel_o = (int *)malloc(sizeof(int)*num_reg_ordered);

  for(i = 0; i < num_reg_ordered; i++)
    kernel_o[i] = KERNEL_ordered_den;


  if(num_reg_continuous != 0) {
    initialize_kernel_regression_asymptotic_constants(KERNEL_den,
                                                      num_reg_continuous,
                                                      &INT_KERNEL_P,
                                                      &K_INT_KERNEL_P,
                                                      &INT_KERNEL_PM_HALF,
                                                      &DIFF_KER_PPM);
  } else {
    INT_KERNEL_P = 1.0;
    K_INT_KERNEL_P = 1.0;
  }

  matrix_bandwidth = alloc_matd(bwmdim,num_reg_continuous);
  lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

  if(kernel_bandwidth_mean(KERNEL_den,
                           BANDWIDTH_den,
                           num_obs_train,
                           num_obs_eval,
                           0,0,0,
                           num_reg_continuous,
                           num_reg_unordered,
                           num_reg_ordered,
                           0,
                           vector_scale_factor,
                           NULL, NULL,
                           matrix_X_continuous_train,
                           matrix_X_continuous_eval,
	                           NULL,
	                           matrix_bandwidth,
	                           lambda)==1){
	    error("\n** Error: invalid bandwidth.");
	  }

  if((num_reg_continuous + num_reg_unordered + num_reg_ordered) > 0){
    int ok_all = 1;

    if(num_reg_continuous > 0){
      const double rel_tol = np_cont_largeh_rel_tol();
      ov_cont_ok = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
      ov_cont_hmin = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
      ov_cont_k0 = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
      ok_all = (ov_cont_ok != NULL) && (ov_cont_hmin != NULL) && (ov_cont_k0 != NULL);

      if(ok_all){
        for(i = 0; i < num_reg_continuous; i++){
          const int kern = kernel_c[i];
          double xmin = DBL_MAX, xmax = -DBL_MAX;

          ov_cont_ok[i] = 0;
          ov_cont_hmin[i] = DBL_MAX;
          ov_cont_k0[i] = 0.0;
          if(!np_cont_largeh_kernel_supported(kern)) continue;

          for(int j = 0; j < num_obs_train; j++){
            const double v = matrix_X_continuous_train[i][j];
            if(!isfinite(v)) continue;
            xmin = MIN(xmin, v);
            xmax = MAX(xmax, v);
          }

          for(int j = 0; j < num_obs_eval; j++){
            const double v = matrix_X_continuous_eval[i][j];
            if(!isfinite(v)) continue;
            xmin = MIN(xmin, v);
            xmax = MAX(xmax, v);
          }

          if(xmax >= xmin){
            const double utol = np_cont_largeh_utol(kern, rel_tol);
            if(utol > 0.0 && isfinite(utol)){
              ov_cont_ok[i] = 1;
              ov_cont_hmin[i] = (xmax - xmin)/utol;
              ov_cont_k0[i] = np_cont_largeh_k0(kern);
            }
          }
        }
      }
    }

    if(ok_all && num_reg_unordered > 0){
      ov_disc_uno_ok = (int *)calloc((size_t)num_reg_unordered, sizeof(int));
      ov_disc_uno_const = (double *)malloc((size_t)num_reg_unordered*sizeof(double));
      ok_all = (ov_disc_uno_ok != NULL) && (ov_disc_uno_const != NULL);
      if(ok_all){
        double (* const ukf[])(int, double, int) = {
          np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
          np_score_uaa, np_score_unli_racine
        };
        const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
        for(i = 0; i < num_reg_unordered; i++){
          const int ku = kernel_u[i];
          const int ncat = (num_categories != NULL) ? num_categories[i] : 0;
          const double lam = lambda[i];
          ov_disc_uno_ok[i] = 0;
          ov_disc_uno_const[i] = 0.0;
          if(ku < 0 || ku >= nuk) continue;
          if(!np_disc_near_upper(ku, lam, ncat)) continue;
          {
            const double ks = ukf[ku](1, lam, ncat);
            const double kd = ukf[ku](0, lam, ncat);
            if(np_disc_near_const_kernel(ks, kd)){
              ov_disc_uno_ok[i] = 1;
              ov_disc_uno_const[i] = 0.5*(ks + kd);
            }
          }
        }
      }
    }

    if(ok_all && num_reg_ordered > 0){
      ov_disc_ord_ok = (int *)calloc((size_t)num_reg_ordered, sizeof(int));
      ov_disc_ord_const = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
      ok_all = (ov_disc_ord_ok != NULL) && (ov_disc_ord_const != NULL);
      if(ok_all){
        double (* const okf[])(double, double, double, double, double) = {
          np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
          np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
          np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
          np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
        };
        const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
        for(i = 0; i < num_reg_ordered; i++){
          const int oi = i + num_reg_unordered;
          const int ko = kernel_o[i];
          const int ncat = (num_categories != NULL) ? num_categories[oi] : 0;
          const double lam = lambda[oi];
          ov_disc_ord_ok[i] = 0;
          ov_disc_ord_const[i] = 0.0;
          if(ko < 0 || ko >= nok) continue;
          if(ncat <= 0 || matrix_categorical_vals == NULL) continue;
          if(!np_disc_ordered_near_upper(ko, lam)) continue;
          {
            const double cl = matrix_categorical_vals[oi][0];
            const double ch = matrix_categorical_vals[oi][ncat - 1];
            const double k0 = okf[ko](cl, cl, lam, cl, ch);
            const double k1 = okf[ko](cl, ch, lam, cl, ch);
            if(np_disc_near_const_kernel(k0, k1)){
              ov_disc_ord_ok[i] = 1;
              ov_disc_ord_const[i] = 0.5*(k0 + k1);
            }
          }
        }
      }
    }

    if(ok_all){
      np_gate_ctx_set(&gate_ctx_local,
                      num_reg_continuous,
                      num_reg_unordered,
                      num_reg_ordered,
                      kernel_c,
                      kernel_u,
                      kernel_o,
                      operator,
                      ov_cont_ok,
                      ov_cont_hmin,
                      ov_cont_k0,
                      ov_disc_uno_ok,
                      ov_disc_uno_const,
                      ov_disc_ord_ok,
                      ov_disc_ord_const);
      gate_override_active = 1;
    }
  }

  kernel_weighted_sum_np_ctx(kernel_c,
                             kernel_u,
                             kernel_o,
                             BANDWIDTH_den,
                             num_obs_train,
                             num_obs_eval,
                             num_reg_unordered,
                             num_reg_ordered,
                             num_reg_continuous,
                             0, // don't leave one out
                             0,
                             1, // kernel_pow = 1
                             1, // bandwidth_divide = FALSE when not adaptive
                             0,
                             0, // symmetric
                             0, // gather-scatter sum
                             0, // do not drop train
                             0, // do not drop train
                             operator, // dens or dist
                             OP_NOOP, // no permutations
                             0, // no score
                             0, // no ocg
                             NULL,
                             0, // explicity suppress parallel
                             0,
                             0,
                             int_TREE_X,
                             0,
                             kdt_extern_X,
                             NULL, NULL, NULL,
                             matrix_X_unordered_train,
                             matrix_X_ordered_train,
                             matrix_X_continuous_train,
                             matrix_X_unordered_eval,
                             matrix_X_ordered_eval,
                             matrix_X_continuous_eval,
                             NULL, // no ys
                             NULL, // no weights
                             NULL, // no sgn
                             vector_scale_factor,
                             1,matrix_bandwidth,matrix_bandwidth,lambda,
                             num_categories,
                             matrix_categorical_vals, // if dist mcv (possibly) necessary
                             NULL, // no ocg
                             pdf,  // weighted sum
                             NULL, // no permutations
                             NULL, // do not return kernel weights
                             &gate_ctx_local); // no permutation kernel weights

  if((BANDWIDTH_den == BW_FIXED) && (dop == OP_NORMAL)){
    for(l = 0, pnh = num_obs_train; l < num_reg_continuous; l++){      
      pnh *= matrix_bandwidth[l][0];
    }
  }


  if (dop == OP_NORMAL) {
    for(i = 0, *log_likelihood = 0.0; i < num_obs_eval; i++){
      pdf[i] /= (double)num_obs_train;
      *log_likelihood += (pdf[i] < DBL_MIN) ? log_DBL_MIN : log(pdf[i]);

      if((BANDWIDTH_den == BW_GEN_NN) && (dop == OP_NORMAL)){
        for(l = 0, pnh = num_obs_train; l < num_reg_continuous; l++){
          pnh *= matrix_bandwidth[l][i];
        }
      }

      pdf_stderr[i] = sqrt(pdf[i]*K_INT_KERNEL_P/pnh);
    }
  } else {
    for(i = 0, *log_likelihood = 0.0; i < num_obs_eval; i++){
      pdf[i] /= (double)num_obs_train;
      pdf_stderr[i] = sqrt(pdf[i]*(1.0-pdf[i])/(double)num_obs_train);
    }
  }

  free(operator);
  free(kernel_c);
  free(kernel_u);
  free(kernel_o);
  free(lambda);
  free_mat(matrix_bandwidth, num_reg_continuous);
  if(gate_override_active) np_gate_ctx_clear(&gate_ctx_local);
  if(ov_cont_ok != NULL) free(ov_cont_ok);
  if(ov_cont_hmin != NULL) free(ov_cont_hmin);
  if(ov_cont_k0 != NULL) free(ov_cont_k0);
  if(ov_disc_uno_ok != NULL) free(ov_disc_uno_ok);
  if(ov_disc_uno_const != NULL) free(ov_disc_uno_const);
  if(ov_disc_ord_ok != NULL) free(ov_disc_ord_ok);
  if(ov_disc_ord_const != NULL) free(ov_disc_ord_const);

}

int np_kernel_estimate_con_density_categorical_leave_one_out_cv(int KERNEL_den,
                                                                int KERNEL_unordered_den,
                                                                int KERNEL_ordered_den,
                                                                int KERNEL_reg,
                                                                int KERNEL_unordered_reg,
                                                                int KERNEL_ordered_reg,
                                                                int BANDWIDTH_den,
                                                                int num_obs,
                                                                int num_var_unordered,
                                                                int num_var_ordered,
                                                                int num_var_continuous,
                                                                int num_reg_unordered,
                                                                int num_reg_ordered,
                                                                int num_reg_continuous,
                                                                double **matrix_Y_unordered,
                                                                double **matrix_Y_ordered,
                                                                double **matrix_Y_continuous,
                                                                double **matrix_X_unordered,
                                                                double **matrix_X_ordered,
                                                                double **matrix_X_continuous,
                                                                double **matrix_XY_unordered, 
                                                                double **matrix_XY_ordered, 
                                                                double **matrix_XY_continuous, 
                                                                double *vector_scale_factor,
                                                                int *num_categories,
                                                                double *cv){
  np_gate_override_clear();

  if(((BANDWIDTH_den == BW_FIXED) || (BANDWIDTH_den == BW_GEN_NN) || (BANDWIDTH_den == BW_ADAP_NN)) &&
     (int_ll_extern == LL_LP))
    return np_conditional_density_cvml_lp_stream(vector_scale_factor, cv);

  const int num_reg = num_reg_continuous+num_reg_unordered+num_reg_ordered;
  const int num_cvar = num_reg_continuous + num_var_continuous;
  const int num_uvar = num_reg_unordered + num_var_unordered;
  const int num_ovar = num_reg_ordered + num_var_ordered;
  const int num_all_var = num_reg + num_var_continuous + num_var_unordered + num_var_ordered;
  const int bwmdim = (BANDWIDTH_den==BW_GEN_NN)?num_obs:
    ((BANDWIDTH_den==BW_ADAP_NN)?num_obs:1);

	//const double log_DBL_MIN = log(DBL_MIN);

  int i; 
  int ret = 0;

  int * operator = NULL;

  int num_obs_alloc;

#ifdef MPI2
  int stride_t = MAX((int)ceil((double) num_obs / (double) iNum_Processors),1);
  
  num_obs_alloc = stride_t*iNum_Processors;
#else
  num_obs_alloc = num_obs;
#endif

  double * vsf_xy = NULL, * vsf_x = NULL;
  double ** matrix_bandwidth_x = NULL, ** matrix_bandwidth_xy = NULL;
  double * lambda_x = NULL, * lambda_xy = NULL;

  vsf_xy = (double *)malloc(num_all_var*sizeof(double));

  if(vsf_xy == NULL)
    error("failed to allocate vsf_xy");

  vsf_x = (double *)malloc(num_reg*sizeof(double));

  if(vsf_x == NULL)
    error("failed to allocate vsf_x");

  double * rhod = (double *)malloc(num_obs_alloc*sizeof(double));
  
  if(rhod == NULL)
    error("failed to allocate rhod");

  double * rhon = (double *)malloc(num_obs_alloc*sizeof(double));

  if(rhon == NULL)
    error("failed to allocate rhon");

  operator = (int *)malloc(sizeof(int)*num_all_var);

  for(i = 0; i < num_all_var; i++)
    operator[i] = OP_NORMAL;


  int * kernel_cx = NULL, * kernel_ux = NULL, * kernel_ox = NULL;
  int * kernel_cxy = NULL, * kernel_uxy = NULL, * kernel_oxy = NULL;

  // x data
  kernel_cx = (int *)malloc(sizeof(int)*num_reg_continuous);

  for(i = 0; i < num_reg_continuous; i++)
    kernel_cx[i] = KERNEL_reg;

  kernel_ux = (int *)malloc(sizeof(int)*num_reg_unordered);

  for(i = 0; i < num_reg_unordered; i++)
    kernel_ux[i] = KERNEL_unordered_reg;

  kernel_ox = (int *)malloc(sizeof(int)*num_reg_ordered);

  for(i = 0; i < num_reg_ordered; i++)
    kernel_ox[i] = KERNEL_ordered_reg;

  
  // xy data
  kernel_cxy = (int *)malloc(sizeof(int)*num_cvar);

  for(i = 0; i < num_reg_continuous; i++)
    kernel_cxy[i] = KERNEL_reg;

  for(i = num_reg_continuous; i < num_cvar; i++)
    kernel_cxy[i] = KERNEL_den;

  kernel_uxy = (int *)malloc(sizeof(int)*num_uvar);

  for(i = 0; i < num_reg_unordered; i++)
    kernel_uxy[i] = KERNEL_unordered_reg;

  for(i = num_reg_unordered; i < num_uvar; i++)
    kernel_uxy[i] = KERNEL_unordered_den;

  kernel_oxy = (int *)malloc(sizeof(int)*num_ovar);

  for(i = 0; i < num_reg_ordered; i++)
    kernel_oxy[i] = KERNEL_ordered_reg;

  for(i = num_reg_ordered; i < num_ovar; i++)
    kernel_oxy[i] = KERNEL_ordered_den;

  // put the correct bws in vsf_x, and vsf_xy

  np_splitxy_vsf_mcv_nc(num_var_unordered, num_var_ordered, num_var_continuous,
                        num_reg_unordered, num_reg_ordered, num_reg_continuous,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        vsf_x,
                        NULL,
                        vsf_xy,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  matrix_bandwidth_x = alloc_matd(bwmdim,num_reg_continuous);
  matrix_bandwidth_xy = alloc_matd(bwmdim,num_cvar);
  lambda_x = alloc_vecd(num_reg_unordered+num_reg_ordered);
  lambda_xy = alloc_vecd(num_uvar+num_ovar);

  if(kernel_bandwidth_mean(KERNEL_reg,
                           BANDWIDTH_den,
                           num_obs,
                           num_obs,
                           0,0,0,
                           num_reg_continuous,
                           num_reg_unordered,
                           num_reg_ordered,
                           0,
                           vsf_x,
                           NULL,NULL,
                           matrix_X_continuous,
                           matrix_X_continuous,
                           NULL,
                           matrix_bandwidth_x,
                           lambda_x)==1){
    error("\n** Error: invalid bandwidth.");
  }

  if(kernel_bandwidth_mean(KERNEL_reg,
                           BANDWIDTH_den,
                           num_obs,
                           num_obs,
                           0,0,0,
                           num_cvar,
                           num_uvar,
                           num_ovar,
                           0,
                           vsf_xy,
                           NULL,NULL,
                           matrix_XY_continuous,
                           matrix_XY_continuous,
                           NULL,
                           matrix_bandwidth_xy,
                           lambda_xy)==1){
    error("\n** Error: invalid bandwidth.");
  }

  /* All-large-X shortcut for conditional-density ML CV:
     if X kernels are effectively constant, CV reduces to unconditional
     leave-one-out density likelihood on Y only. */
  if(BANDWIDTH_den == BW_FIXED){
    int all_large_gate = 1;
    int ok_all = 1;
    int *x_cont_ok = NULL, *x_disc_uno_ok = NULL, *x_disc_ord_ok = NULL;
    double *x_cont_hmin = NULL, *x_cont_k0 = NULL;
    double *x_disc_uno_const = NULL, *x_disc_ord_const = NULL;

    if(num_reg_continuous > 0){
      const double rel_tol = np_cont_largeh_rel_tol();
      x_cont_ok = (int *)calloc((size_t)num_reg_continuous, sizeof(int));
      x_cont_hmin = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
      x_cont_k0 = (double *)malloc((size_t)num_reg_continuous*sizeof(double));
      ok_all = (x_cont_ok != NULL) && (x_cont_hmin != NULL) && (x_cont_k0 != NULL);
      if(ok_all){
        for(i = 0; i < num_reg_continuous; i++){
          const int kern = kernel_cx[i];
          double xmin = DBL_MAX, xmax = -DBL_MAX;
          x_cont_ok[i] = 0; x_cont_hmin[i] = DBL_MAX; x_cont_k0[i] = 0.0;
          if(!np_cont_largeh_kernel_supported(kern)) continue;
          for(int j = 0; j < num_obs; j++){
            const double v = matrix_X_continuous[i][j];
            if(!isfinite(v)) continue;
            xmin = MIN(xmin, v); xmax = MAX(xmax, v);
          }
          if(xmax >= xmin){
            const double utol = np_cont_largeh_utol(kern, rel_tol);
            if(utol > 0.0 && isfinite(utol)){
              x_cont_ok[i] = 1;
              x_cont_hmin[i] = (xmax - xmin)/utol;
              x_cont_k0[i] = np_cont_largeh_k0(kern);
            }
          }
        }
      }
    }

    if(ok_all && num_reg_unordered > 0){
      x_disc_uno_ok = (int *)calloc((size_t)num_reg_unordered, sizeof(int));
      x_disc_uno_const = (double *)malloc((size_t)num_reg_unordered*sizeof(double));
      ok_all = (x_disc_uno_ok != NULL) && (x_disc_uno_const != NULL);
      if(ok_all){
        double (* const ukf[])(int, double, int) = {
          np_uaa, np_unli_racine, np_econvol_uaa, np_econvol_unli_racine,
          np_score_uaa, np_score_unli_racine
        };
        const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
        for(i = 0; i < num_reg_unordered; i++){
          const int ku = kernel_ux[i];
          const int ncat = (num_categories_extern_X != NULL) ? num_categories_extern_X[i] : 0;
          const double lam = lambda_x[i];
          x_disc_uno_ok[i] = 0; x_disc_uno_const[i] = 0.0;
          if(ku < 0 || ku >= nuk) continue;
          if(!np_disc_near_upper(ku, lam, ncat)) continue;
          {
            const double ks = ukf[ku](1, lam, ncat);
            const double kd = ukf[ku](0, lam, ncat);
            if(np_disc_near_const_kernel(ks, kd)){
              x_disc_uno_ok[i] = 1;
              x_disc_uno_const[i] = 0.5*(ks + kd);
            }
          }
        }
      }
    }

    if(ok_all && num_reg_ordered > 0){
      x_disc_ord_ok = (int *)calloc((size_t)num_reg_ordered, sizeof(int));
      x_disc_ord_const = (double *)malloc((size_t)num_reg_ordered*sizeof(double));
      ok_all = (x_disc_ord_ok != NULL) && (x_disc_ord_const != NULL);
      if(ok_all){
        double (* const okf[])(double, double, double, double, double) = {
          np_owang_van_ryzin, np_oli_racine, np_onli_racine, np_oracine_li_yan,
          np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine, np_econvol_oracine_li_yan,
          np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine, np_score_oracine_li_yan,
          np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine, np_cdf_oracine_li_yan
        };
        const int nok = (int)(sizeof(okf)/sizeof(okf[0]));
        for(i = 0; i < num_reg_ordered; i++){
          const int oi = num_reg_unordered + i;
          const int ko = kernel_ox[i];
          const int ncat = (num_categories_extern_X != NULL) ? num_categories_extern_X[oi] : 0;
          const double lam = lambda_x[oi];
          x_disc_ord_ok[i] = 0; x_disc_ord_const[i] = 0.0;
          if((ko < 0) || (ko >= nok) || (!isfinite(lam)) || (!np_disc_ordered_near_upper(ko, lam)) ||
             (ncat <= 0) || (matrix_categorical_vals_extern_X == NULL))
            continue;
          {
            const double cl = matrix_categorical_vals_extern_X[oi][0];
            const double ch = matrix_categorical_vals_extern_X[oi][ncat - 1];
            const double k0 = okf[ko](cl, cl, lam, cl, ch);
            const double k1 = okf[ko](cl, ch, lam, cl, ch);
            if(np_disc_near_const_kernel(k0, k1)){
              x_disc_ord_ok[i] = 1;
              x_disc_ord_const[i] = 0.5*(k0 + k1);
            }
          }
        }
      }
    }

    if(!ok_all){
      all_large_gate = 0;
    } else {
      for(i = 0; i < num_reg_continuous; i++){
        const double bw = matrix_bandwidth_x[i][0];
        if((x_cont_ok == NULL) || (!x_cont_ok[i]) || (!isfinite(bw)) ||
           (bw <= 0.0) || (bw < x_cont_hmin[i])){
          all_large_gate = 0;
          break;
        }
      }

      if(all_large_gate){
        for(i = 0; i < num_reg_unordered; i++){
          if((x_disc_uno_ok == NULL) || (!x_disc_uno_ok[i])){
            all_large_gate = 0;
            break;
          }
        }
      }

      if(all_large_gate){
        for(i = 0; i < num_reg_ordered; i++){
          if((x_disc_ord_ok == NULL) || (!x_disc_ord_ok[i])){
            all_large_gate = 0;
            break;
          }
        }
      }
    }

    if(x_cont_ok != NULL) free(x_cont_ok);
    if(x_cont_hmin != NULL) free(x_cont_hmin);
    if(x_cont_k0 != NULL) free(x_cont_k0);
    if(x_disc_uno_ok != NULL) free(x_disc_uno_ok);
    if(x_disc_uno_const != NULL) free(x_disc_uno_const);
    if(x_disc_ord_ok != NULL) free(x_disc_ord_ok);
    if(x_disc_ord_const != NULL) free(x_disc_ord_const);

    if(all_large_gate){
      int fast_failed = 0;
      int *operator_y = NULL, *kernel_cy = NULL, *kernel_uy = NULL, *kernel_oy = NULL;
      double *vsf_y = NULL, *rho_y = NULL;
      double **matrix_bandwidth_y = NULL;
      double *lambda_y = NULL;
      const int num_var_tot = num_var_continuous + num_var_unordered + num_var_ordered;
      const int y_bwmdim = 1;
      const double log_DBL_MIN = log(DBL_MIN);

      operator_y = (int *)malloc(sizeof(int)*MAX(1,num_var_tot));
      kernel_cy = (int *)malloc(sizeof(int)*MAX(1,num_var_continuous));
      kernel_uy = (int *)malloc(sizeof(int)*MAX(1,num_var_unordered));
      kernel_oy = (int *)malloc(sizeof(int)*MAX(1,num_var_ordered));
      vsf_y = (double *)malloc(sizeof(double)*MAX(1,num_var_tot));
      rho_y = (double *)malloc(sizeof(double)*MAX(1,num_obs_alloc));
      matrix_bandwidth_y = alloc_matd(y_bwmdim, num_var_continuous);
      lambda_y = alloc_vecd(num_var_unordered + num_var_ordered);

      if(operator_y == NULL || kernel_cy == NULL || kernel_uy == NULL || kernel_oy == NULL ||
         vsf_y == NULL || rho_y == NULL ||
         ((num_var_continuous > 0) && (matrix_bandwidth_y == NULL)) ||
         (((num_var_unordered + num_var_ordered) > 0) && (lambda_y == NULL))){
        fast_failed = 1;
      }

      if(!fast_failed){
        for(i = 0; i < num_var_tot; i++) operator_y[i] = OP_NORMAL;
        for(i = 0; i < num_var_continuous; i++) kernel_cy[i] = KERNEL_den;
        for(i = 0; i < num_var_unordered; i++) kernel_uy[i] = KERNEL_unordered_den;
        for(i = 0; i < num_var_ordered; i++) kernel_oy[i] = KERNEL_ordered_den;

        np_splitxy_vsf_mcv_nc(num_var_unordered, num_var_ordered, num_var_continuous,
                              num_reg_unordered, num_reg_ordered, num_reg_continuous,
                              vector_scale_factor, NULL, NULL,
                              NULL, vsf_y, NULL,
                              NULL, NULL, NULL,
                              NULL, NULL, NULL);

        if(kernel_bandwidth_mean(KERNEL_den,
                                 BANDWIDTH_den,
                                 num_obs,
                                 num_obs,
                                 0, 0, 0,
                                 num_var_continuous,
                                 num_var_unordered,
                                 num_var_ordered,
                                 0,
                                 vsf_y,
                                 NULL, NULL,
                                 matrix_Y_continuous,
                                 matrix_Y_continuous,
                                 NULL,
                                 matrix_bandwidth_y,
                                 lambda_y) == 1){
          fast_failed = 1;
        }
      }

      if(!fast_failed){
        np_activate_bounds_y();
        kernel_weighted_sum_np(kernel_cy,
                               kernel_uy,
                               kernel_oy,
                               BANDWIDTH_den,
                               num_obs,
                               num_obs,
                               num_var_unordered,
                               num_var_ordered,
                               num_var_continuous,
                               1, // leave one out
                               0,
                               1,
                               1,
                               0,
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
                               int_TREE_Y,
                               0,
                               kdt_extern_Y,
                               NULL, NULL, NULL,
                               matrix_Y_unordered,
                               matrix_Y_ordered,
                               matrix_Y_continuous,
                               matrix_Y_unordered,
                               matrix_Y_ordered,
                               matrix_Y_continuous,
                               NULL,
                               NULL,
                               NULL,
                               vsf_y,
                               1,
                               matrix_bandwidth_y,
                               matrix_bandwidth_y,
                               lambda_y,
                               num_categories,
                               matrix_categorical_vals_extern,
                               NULL,
                               rho_y,
                               NULL,
                               NULL,
                               NULL);

        *cv = 0.0;
        ret = 0;
        for(i = 0; i < num_obs; i++){
          if(rho_y[i] < DBL_MIN){
            *cv -= log_DBL_MIN;
            ret = 1;
          } else {
            *cv -= log(rho_y[i]/(num_obs - 1.0));
          }
        }
      }

      if(operator_y != NULL) free(operator_y);
      if(kernel_cy != NULL) free(kernel_cy);
      if(kernel_uy != NULL) free(kernel_uy);
      if(kernel_oy != NULL) free(kernel_oy);
      if(vsf_y != NULL) free(vsf_y);
      if(rho_y != NULL) free(rho_y);
      if(lambda_y != NULL) free(lambda_y);
      if(matrix_bandwidth_y != NULL) free_mat(matrix_bandwidth_y, num_var_continuous);

      if(!fast_failed){
        np_fastcv_alllarge_hits++;
        goto cleanup_cvml_return;
      }
    }
  }

  // xy
  np_activate_bounds_xy();
  kernel_weighted_sum_np(kernel_cxy,
                         kernel_uxy,
                         kernel_oxy,
                         BANDWIDTH_den,
                         num_obs,
                         num_obs,
                         num_uvar,
                         num_ovar,
                         num_cvar,
                         1, // leave one out 
                         0,
                         1, // kernel_pow = 1
                         1, // bandwidth_divide = FALSE when not adaptive
                         0, 
                         0, // symmetric
                         0, // gather-scatter sum
                         0, // do not drop train
                         0, // do not drop train
                         operator, // no convolution
                         OP_NOOP, // no permutations
                         0, // no score
                         0, // no ocg
                         NULL,
                         0, //  do not explicity suppress parallel
                         0,
                         0,
                         int_TREE_XY,
                         0,
                         kdt_extern_XY,
                         NULL, NULL, NULL,
                         matrix_XY_unordered,
                         matrix_XY_ordered,
                         matrix_XY_continuous,
                         matrix_XY_unordered,
                         matrix_XY_ordered,
                         matrix_XY_continuous,
                         NULL,
                         NULL,
                         NULL,
                         vsf_xy,
                         1,matrix_bandwidth_xy,matrix_bandwidth_xy,lambda_xy,
                         num_categories_extern_XY,
                         matrix_categorical_vals_extern_XY,
                         NULL,
                         rhon,  // weighted sum
                         NULL, // no permutations
                         NULL, // do not return kernel weights
                         NULL);

  //x
  np_activate_bounds_x();
  kernel_weighted_sum_np(kernel_cx,
                         kernel_ux,
                         kernel_ox,
                         BANDWIDTH_den,
                         num_obs,
                         num_obs,
                         num_reg_unordered,
                         num_reg_ordered,
                         num_reg_continuous,
                         1, // leave one out 
                         0,
                         1, // kernel_pow = 1
                         1, // bandwidth_divide = FALSE when not adaptive
                         0, 
                         0, // symmetric
                         0, // gather-scatter sum
                         0, // do not drop train
                         0, // do not drop train
                         operator, // no convolution
                         OP_NOOP, // no permutations
                         0, // no score
                         0, // no ocg
                         NULL,
                         0, //  do not explicity suppress parallel
                         0,
                         0,
                         int_TREE_X,
                         0,
                         kdt_extern_X,
                         NULL, NULL, NULL,
                         matrix_X_unordered,
                         matrix_X_ordered,
                         matrix_X_continuous,
                         matrix_X_unordered,
                         matrix_X_ordered,
                         matrix_X_continuous,
                         NULL,
                         NULL,
                         NULL,
                         vsf_x,
                         1,matrix_bandwidth_x,matrix_bandwidth_x,lambda_x,
                         num_categories_extern_X,
                         matrix_categorical_vals_extern_X,
                         NULL,
                         rhod,  // weighted sum
                         NULL, // no permutations
                         NULL, // do not return kernel weights
                         NULL);
  
  for(i = 0, *cv = 0.0; i < num_obs; i++){
    if((rhon[i] == 0.0) || (rhod[i] == 0.0)){
      ret = 1;
      break;
      //*cv -= log_DBL_MIN;
    } else {
    *cv -= log(rhon[i]) - log(rhod[i]);
    }
  }

cleanup_cvml_return:
  free(operator);
  free(kernel_cx);
  free(kernel_ux);
  free(kernel_ox);

  free(kernel_cxy);
  free(kernel_uxy);
  free(kernel_oxy);

  free(rhon);
  free(rhod);

  free(lambda_x);
  free(lambda_xy);
  free_mat(matrix_bandwidth_x, num_reg_continuous);
  free_mat(matrix_bandwidth_xy, num_cvar);

  free(vsf_xy);
  free(vsf_x);
  np_gate_override_clear();
  return(ret);

}

void np_kernel_estimate_con_dens_dist_categorical(
int KERNEL_Y,
int KERNEL_unordered_Y,
int KERNEL_ordered_Y,
int KERNEL_X,
int KERNEL_unordered_X,
int KERNEL_ordered_X,
int BANDWIDTH_den,
int yop,
int num_obs_train,
int num_obs_eval,
int num_Y_unordered,
int num_Y_ordered,
int num_Y_continuous,
int num_X_unordered,
int num_X_ordered,
int num_X_continuous,
double **matrix_XY_unordered_train, 
double **matrix_XY_ordered_train, 
double **matrix_XY_continuous_train, 
double **matrix_XY_unordered_eval, 
double **matrix_XY_ordered_eval, 
double **matrix_XY_continuous_eval, 
double *vector_scale_factor,
int *num_categories,
int *num_categories_XY,
double ** matrix_categorical_vals,
double ** matrix_categorical_vals_XY,
double * kdf,
double * kdf_stderr,
double ** kdf_deriv,
double ** kdf_deriv_stderr,
double * log_likelihood
){


  const int num_X = num_X_continuous+num_X_unordered+num_X_ordered;
  const int num_Y = num_Y_continuous+num_Y_unordered+num_Y_ordered;
  const int num_cXY = num_X_continuous + num_Y_continuous;
  const int num_uXY = num_X_unordered + num_Y_unordered;
  const int num_oXY = num_X_ordered + num_Y_ordered;
  const int num_XY = num_X + num_Y;

  const int is_cpdf = (yop == OP_NORMAL);

  int i,l,k;
  
  double INT_KERNEL_P;					 /* Integral of K(z)^2 */
  double K_INT_KERNEL_P;				 /*  K^p */
  /* Integral of K(z-0.5)*K(z+0.5) */
	double INT_KERNEL_PM_HALF = 0.0;
	double DIFF_KER_PPM = 0.0;		 /* Difference between int K(z)^p and int K(z-.5)K(z+.5) */


  double ** matrix_bandwidth_Y = NULL, * lambda = NULL;
  double ** matrix_bandwidth_X = NULL;

  const int do_grad = (kdf_deriv != NULL); 
  const int do_gerr = (kdf_deriv_stderr != NULL);

  struct th_table * otabs = NULL;
  struct th_entry * ret = NULL;
  int ** matrix_ordered_indices = NULL;

  const int bwmdim = (BANDWIDTH_den==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_den==BW_ADAP_NN)?num_obs_train:1);

  double pnh = 1.0;

  const double log_DBL_MIN = log(DBL_MIN);

  int num_obs_eval_alloc;

#ifdef MPI2
  int stride_t = MAX((int)ceil((double) num_obs_eval / (double) iNum_Processors),1);
  
  num_obs_eval_alloc = stride_t*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval;
#endif

  int * kernel_cXY = NULL, * kernel_uXY = NULL, * kernel_oXY = NULL;
  int * operator_XY = NULL, * operator_X = NULL;

  double * vsf_XY = NULL, * vsf_X = NULL;

  // used for partial node search
  int * icX = NULL;
  NL nls = {.node = NULL, .n = 0, .nalloc = 0};

  double * ksd = NULL, * ksn = NULL;

  double * permn = NULL, * permd = NULL;

  int * bpso = NULL;

  const int p = is_cpdf ? num_cXY : num_X_continuous;

  if(p != 0) {
    initialize_kernel_regression_asymptotic_constants(KERNEL_Y,
                                                      p,
                                                      &INT_KERNEL_P,
                                                      &K_INT_KERNEL_P,
                                                      &INT_KERNEL_PM_HALF,
                                                      &DIFF_KER_PPM);
  } else {
    INT_KERNEL_P = 1.0;
    K_INT_KERNEL_P = 1.0;
  }


  const double gfac = sqrt(DIFF_KER_PPM/INT_KERNEL_P);

  if(do_grad && (num_X_ordered > 0)){
    otabs = (struct th_table *)np_jksum_malloc_array_or_die((size_t)num_X_ordered, sizeof(struct th_table), "np_kernel_estimate_con_dens_dist_categorical otabs");
    matrix_ordered_indices = (int **)np_jksum_malloc_array_or_die((size_t)num_oXY, sizeof(int *), "np_kernel_estimate_con_dens_dist_categorical matrix_ordered_indices");
    int * tc = (int *)np_jksum_malloc_array3_or_die((size_t)num_X_ordered, (size_t)num_obs_eval, sizeof(int), "np_kernel_estimate_con_dens_dist_categorical ordered index buffer");
    for(l = 0; l < num_X_ordered; l++)
      matrix_ordered_indices[l] = tc + l*num_obs_eval;

    //  pad out the moo with some fake entries at the end to stop kernel_weighted_sum_np from
    //  doing anything undefined
    
    for(; l < num_oXY; l++)
      matrix_ordered_indices[l] = tc;

    for(l = 0; l < num_X_ordered; l++){
      if(thcreate_r((size_t)ceil(1.6*num_categories[l+num_X_unordered+num_Y_unordered+num_Y_ordered]), otabs + l) == TH_ERROR)
        error("hash table creation failed");

      for(i = 0; i < num_categories[l+num_X_unordered+num_Y_unordered+num_Y_ordered]; i++){
        struct th_entry centry;
        centry.key.dkey = matrix_categorical_vals[l+num_X_unordered+num_Y_unordered+num_Y_ordered][i];
        centry.data = i;

        if(thsearch_r(&centry, TH_ENTER, &ret, otabs+l) == TH_FAILURE)
          error("insertion into hash table failed");
      }

      // now do lookups
      struct th_entry te;
      te.key.dkey = 0.0;
      te.data = -1;
      ret = &te;

      for(i = 0; i < num_obs_eval; i++){
        if(ret->key.dkey != matrix_XY_ordered_eval[l][i]){
          te.key.dkey = matrix_XY_ordered_eval[l][i];
          if(thsearch_r(&te, TH_SEARCH, &ret, otabs+l) == TH_FAILURE)
            error("hash table lookup failed (which should be impossible)");
        } 

        matrix_ordered_indices[l][i] = ret->data;
      }
    }
  }


  operator_XY = (int *)np_jksum_malloc_array_or_die((size_t)num_XY, sizeof(int), "np_kernel_estimate_con_dens_dist_categorical operator_XY");
  operator_X = (int *)np_jksum_malloc_array_or_die((size_t)num_X, sizeof(int), "np_kernel_estimate_con_dens_dist_categorical operator_X");

  kernel_cXY = (int *)np_jksum_malloc_array_or_die((size_t)num_cXY, sizeof(int), "np_kernel_estimate_con_dens_dist_categorical kernel_cXY");
  kernel_uXY = (int *)np_jksum_malloc_array_or_die((size_t)num_uXY, sizeof(int), "np_kernel_estimate_con_dens_dist_categorical kernel_uXY");
  kernel_oXY = (int *)np_jksum_malloc_array_or_die((size_t)num_oXY, sizeof(int), "np_kernel_estimate_con_dens_dist_categorical kernel_oXY");

  vsf_XY = (double *)np_jksum_malloc_array_or_die((size_t)num_XY, sizeof(double), "np_kernel_estimate_con_dens_dist_categorical vsf_XY");
  vsf_X = (double *)np_jksum_malloc_array_or_die((size_t)num_X, sizeof(double), "np_kernel_estimate_con_dens_dist_categorical vsf_X");

  ksd = (double *)np_jksum_malloc_array_or_die((size_t)num_obs_eval_alloc, sizeof(double), "np_kernel_estimate_con_dens_dist_categorical ksd");
  ksn = (double *)np_jksum_malloc_array_or_die((size_t)num_obs_eval_alloc, sizeof(double), "np_kernel_estimate_con_dens_dist_categorical ksn");

  icX = (int *)np_jksum_malloc_array_or_die((size_t)num_X_continuous, sizeof(int), "np_kernel_estimate_con_dens_dist_categorical icX");

  if(do_grad){
    permn = (double *)np_jksum_malloc_array3_or_die((size_t)num_X, (size_t)num_obs_eval_alloc, sizeof(double), "np_kernel_estimate_con_dens_dist_categorical permn");
    permd = (double *)np_jksum_malloc_array3_or_die((size_t)num_X, (size_t)num_obs_eval_alloc, sizeof(double), "np_kernel_estimate_con_dens_dist_categorical permd");
    bpso = (int *)np_jksum_malloc_array_or_die((size_t)num_XY, sizeof(int), "np_kernel_estimate_con_dens_dist_categorical bpso");

    // only enable gradients for 'X' variables
    for(i = 0; i < num_X_continuous; i++)
      bpso[i] = 1;

    for(; i < num_cXY; i++)
      bpso[i] = 0;

    for(; i < num_cXY + num_X_unordered; i++)
      bpso[i] = 1;

    for(; i < num_cXY + num_uXY; i++)
      bpso[i] = 0;

    for(; i < num_cXY + num_uXY + num_X_ordered; i++)
      bpso[i] = 1;

    for(; i < num_XY; i++)
      bpso[i] = 0;

  }

  matrix_bandwidth_Y = alloc_matd(bwmdim,num_Y_continuous);
  matrix_bandwidth_X = alloc_matd(bwmdim,num_X_continuous);
  lambda = alloc_vecd(num_uXY + num_oXY);

  if(kernel_bandwidth_mean(KERNEL_Y,
                           BANDWIDTH_den,
                           num_obs_train,
                           num_obs_eval,
                           num_Y_continuous,
                           num_Y_unordered,
                           num_Y_ordered,
                           num_X_continuous,
                           num_X_unordered,
                           num_X_ordered,
                           0,
                           vector_scale_factor,
                           matrix_XY_continuous_train + num_X_continuous,
                           matrix_XY_continuous_eval + num_X_continuous,
                           matrix_XY_continuous_train,
	                           matrix_XY_continuous_eval,
	                           matrix_bandwidth_Y,
	                           matrix_bandwidth_X,
	                           lambda)==1){
	    error("\n** Error: invalid bandwidth.");
	  }


  // relevant dimensions for partial tree search
  for(i = 0; i < num_X_continuous; i++)
    icX[i] = i;

  // nodes for pts
  nls.node = (int *)malloc(sizeof(int));
  nls.nalloc = 1;

  nls.node[0] = 0;
  nls.n = 1;

  // merge the kernels and operators

  np_kernelop_xy(KERNEL_Y, KERNEL_unordered_Y, KERNEL_ordered_Y,
		 KERNEL_X, KERNEL_unordered_X, KERNEL_ordered_X,
		 yop, // y operator
		 OP_NORMAL, // x operator
		 num_Y_unordered, num_Y_ordered, num_Y_continuous,
		 num_X_unordered, num_X_ordered, num_X_continuous,
		 NULL, NULL, NULL, NULL, NULL, NULL,
		 kernel_cXY, kernel_uXY, kernel_oXY,
		 operator_X, NULL, operator_XY);


  // put the correct bws in vsf_x, and vsf_xy

  np_splitxy_vsf_mcv_nc(num_Y_unordered, num_Y_ordered, num_Y_continuous,
                        num_X_unordered, num_X_ordered, num_X_continuous,
                        vector_scale_factor,
                        NULL,
                        NULL,
                        vsf_X,
                        NULL,
                        vsf_XY,
                        NULL, NULL, NULL,
                        NULL, NULL, NULL);

  // xy
  np_progress_fit_set_offset(0);
  np_activate_bounds_xy();
  kernel_weighted_sum_np(kernel_cXY,
                         kernel_uXY,
                         kernel_oXY,
                         BANDWIDTH_den,
                         num_obs_train,
                         num_obs_eval,
                         num_uXY,
                         num_oXY,
                         num_cXY,
                         0, // leave one out 
                         0, // '' offset
                         1, // kernel_pow = 1
                         1, // bandwidth_divide 
                         0, // '' weights
                         0, // symmetric
                         0, // gather-scatter sum
                         0, // drop train
                         0, // drop which train
                         operator_XY, 
                         do_grad ? OP_DERIVATIVE : OP_NOOP, // no permutations
                         0, // no score
                         do_grad, // no ocg
                         bpso,
                         1, // do not explicity suppress parallel
                         0, // ncol y
                         0, // ncol w
                         int_TREE_XY, // do tree
                         0, // do partial tree
                         kdt_extern_XY, // which tree
                         NULL, NULL, NULL, // partial tree data
                         matrix_XY_unordered_train,
                         matrix_XY_ordered_train,
                         matrix_XY_continuous_train,
                         matrix_XY_unordered_eval,
                         matrix_XY_ordered_eval,
                         matrix_XY_continuous_eval,
                         NULL, // matrix y
                         NULL, // matrix w
                         NULL, // sgn
                         vsf_XY,
                         0,NULL,NULL,NULL, 
                         num_categories_XY,
                         matrix_categorical_vals_XY,
                         matrix_ordered_indices, 
                         ksn,  // weighted sum
                         permn, // permutations
                         NULL, // do not return kernel weights
                         NULL);

  //x - we assume x is in xy tree order
  //  - we also reuse kernels and operators because we can

  np_progress_fit_set_offset((BANDWIDTH_den == BW_ADAP_NN) ? num_obs_train : num_obs_eval);
  np_activate_bounds_xy();
  kernel_weighted_sum_np(kernel_cXY,
                         kernel_uXY,
                         kernel_oXY,
                         BANDWIDTH_den,
                         num_obs_train,
                         num_obs_eval,
                         num_X_unordered,
                         num_X_ordered,
                         num_X_continuous,
                         0, // leave one out 
                         0,
                         1, // kernel_pow = 1
                         1, // bandwidth_divide 
                         0, 
                         0, // symmetric
                         0, // gather-scatter sum
                         0, // do not drop train
                         0, // do not drop train
                         operator_X, // no convolution
                         do_grad ? OP_DERIVATIVE : OP_NOOP, // no permutations
                         0, // no score
                         do_grad, // no ocg
                         NULL, // bpso = TRUE for all xi
                         1, //  do not explicity suppress parallel
                         0,
                         0,
                         int_TREE_XY,
                         1, // in this case we do want a partial tree
                         kdt_extern_XY,
                         &nls, icX, NULL,
                         matrix_XY_unordered_train,
                         matrix_XY_ordered_train,
                         matrix_XY_continuous_train,
                         matrix_XY_unordered_eval,
                         matrix_XY_ordered_eval,
                         matrix_XY_continuous_eval,
                         NULL,
                         NULL,
                         NULL,
                         vsf_X,
                         0,NULL,NULL,NULL,
                         num_categories + num_Y_unordered + num_Y_ordered,
                         matrix_categorical_vals + num_Y_unordered + num_Y_ordered,
                         matrix_ordered_indices, // moo
                         ksd,  // weighted sum
                         permd, // no permutations
                         NULL, // do not return kernel weights
                         NULL);

  
  if (is_cpdf) {
    if(BANDWIDTH_den == BW_FIXED){
      for(l = 0, pnh = 1.0; l < num_X_continuous; l++){      
        pnh *= matrix_bandwidth_X[l][0];
      }

      for(l = 0; l < num_Y_continuous; l++){      
        pnh *= matrix_bandwidth_Y[l][0];
      }
    }

    for(i = 0; i < num_obs_eval; i++){
      const double sk = copysign(DBL_MIN, ksd[i]) + ksd[i];
      kdf[i] = ksn[i]/sk;
      *log_likelihood += (kdf[i] < DBL_MIN) ? log_DBL_MIN : log(kdf[i]);

      if(BANDWIDTH_den == BW_GEN_NN){
        for(l = 0, pnh = 1.0; l < num_X_continuous; l++){      
          pnh *= matrix_bandwidth_X[l][i];
        }

        for(l = 0; l < num_Y_continuous; l++){
          pnh *= matrix_bandwidth_Y[l][i];
        }
      }

      kdf_stderr[i] = sqrt(kdf[i]*K_INT_KERNEL_P/(pnh*sk));
   
    }
  } else {

    if(BANDWIDTH_den == BW_FIXED){
      for(l = 0, pnh = 1.0; l < num_X_continuous; l++){      
        pnh *= matrix_bandwidth_X[l][0];
      }
    }

    for(i = 0, *log_likelihood = 0.0; i < num_obs_eval; i++){
      const double sk = copysign(DBL_MIN, ksd[i]) + ksd[i];
      kdf[i] = ksn[i]/sk;

      if(BANDWIDTH_den == BW_GEN_NN){
        for(l = 0, pnh = 1.0; l < num_X_continuous; l++){
          pnh *= matrix_bandwidth_X[l][i];
        }
      }

      kdf_stderr[i] = sqrt(kdf[i]*(1.0-kdf[i])*K_INT_KERNEL_P/(pnh*sk));
    }

  }

  if(do_grad) {
    for(l = 0; l < num_X_continuous; l++){
      for(i = 0; i < num_obs_eval; i++){
        const double sk = copysign(DBL_MIN, ksd[i]) + ksd[i];

        kdf_deriv[l][i] = (permn[l*num_obs_eval + i]-kdf[i]*permd[l*num_obs_eval + i])/sk;

        if(do_gerr){
          const double hfac = ((BANDWIDTH_den == BW_ADAP_NN) ? 1.0 : ((BANDWIDTH_den == BW_GEN_NN) ? matrix_bandwidth_X[l][i]:matrix_bandwidth_X[l][0]));
          kdf_deriv_stderr[l][i] = gfac*kdf_stderr[i]/hfac;
        }
      }
    }

    for(; l < num_X_continuous + num_X_unordered; l++){
      for(i = 0; i < num_obs_eval; i++){
        const int li = l*num_obs_eval + i;
        const double sk = copysign(DBL_MIN, permd[li]) + permd[li];
        const double s1 = permn[li]/sk;

        kdf_deriv[l][i] = kdf[i] - s1;
        
        // covariance is missing
        if(do_gerr){
          if(is_cpdf){
            if(BANDWIDTH_den == BW_GEN_NN){
             for(k = 0, pnh = 1.0; k < num_X_continuous; k++){
                pnh *= matrix_bandwidth_X[k][i];
              }

              for(k = 0; k < num_Y_continuous; k++){
                pnh *= matrix_bandwidth_Y[k][i];
              }
            }

            kdf_deriv_stderr[l][i] = sqrt(kdf_stderr[i]*kdf_stderr[i] + s1*K_INT_KERNEL_P/(pnh*sk));
          }
          else {
            if(BANDWIDTH_den == BW_GEN_NN){
             for(k = 0, pnh = 1.0; k < num_X_continuous; k++){
                pnh *= matrix_bandwidth_X[k][i];
              }
            }

            kdf_deriv_stderr[l][i] = sqrt(kdf_stderr[i]*kdf_stderr[i] + s1*(1.0-s1)*K_INT_KERNEL_P/(pnh*sk));
          }
        }
      }
    }

    for(; l < num_X; l++){
      for(i = 0; i < num_obs_eval; i++){
        const int li = l*num_obs_eval + i;
        const double sk = copysign(DBL_MIN, permd[li]) + permd[li];
        const double s1 = permn[li]/sk;

        kdf_deriv[l][i] = (kdf[i] - s1)*((matrix_ordered_indices[l - num_X_continuous - num_X_unordered][i] != 0) ? 1.0 : -1.0);

        // covariance is missing
        if(do_gerr){
          if(is_cpdf){
            if(BANDWIDTH_den == BW_GEN_NN){
             for(k = 0, pnh = 1.0; k < num_X_continuous; k++){
                pnh *= matrix_bandwidth_X[k][i];
              }

              for(k = 0; k < num_Y_continuous; k++){
                pnh *= matrix_bandwidth_Y[k][i];
              }
            }

            kdf_deriv_stderr[l][i] = sqrt(kdf_stderr[i]*kdf_stderr[i] + s1*K_INT_KERNEL_P/(pnh*sk));
          } else {
            if(BANDWIDTH_den == BW_GEN_NN){
             for(k = 0, pnh = 1.0; k < num_X_continuous; k++){
                pnh *= matrix_bandwidth_X[k][i];
              }
            }

            kdf_deriv_stderr[l][i] = sqrt(kdf_stderr[i]*kdf_stderr[i] + s1*(1.0-s1)*K_INT_KERNEL_P/(pnh*sk));
          }
        }
      }
    }


  }

  free(operator_XY);
  free(operator_X);

  free(kernel_cXY);
  free(kernel_uXY);
  free(kernel_oXY);

  free(ksn);
  free(ksd);

  free(vsf_XY);
  free(vsf_X);

  free(icX);
  clean_nl(&nls);

  if(do_grad){
    free(permn);
    free(permd);
    free(bpso);
  }


  // clean up hash stuff
  if(do_grad && (num_X_ordered > 0)){
    for(l = 0; l < num_X_ordered; l++)
      thdestroy_r(otabs+l);
    free(otabs);
    free(matrix_ordered_indices[0]);
    free(matrix_ordered_indices);
  }

  free(lambda);
  free_mat(matrix_bandwidth_Y, num_Y_continuous);
  free_mat(matrix_bandwidth_X, num_X_continuous);
}

void np_splitxy_vsf_mcv_nc(const int num_var_unordered,
                           const int num_var_ordered,
                           const int num_var_continuous,
                           const int num_reg_unordered,
                           const int num_reg_ordered,
                           const int num_reg_continuous,
                           const double * const vector_scale_factor,
                           const int * const num_categories,
                           double ** matrix_categorical_vals,
                           double * vsf_x,
                           double * vsf_y,
                           double * vsf_xy,
                           int * nc_x,
                           int * nc_y,
                           int * nc_xy,
                           double ** mcv_x,
                           double ** mcv_y,
                           double ** mcv_xy){

  int i, j, l;

  const int num_cvar = num_var_continuous + num_reg_continuous;
  const int num_uvar = num_var_unordered + num_reg_unordered;
  const int num_ovar = num_var_ordered + num_reg_ordered;

  const int num_catvar = num_uvar + num_ovar;

  const int num_all_var = num_var_continuous + num_reg_continuous + num_var_unordered + num_reg_unordered + num_var_ordered + num_reg_ordered;

  // set up xy bws
  
  if(vsf_xy != NULL){
    for(i = 0, l = 0; i < num_cvar; i++, l++){
      vsf_xy[l] = vector_scale_factor[i];
    }

    // copy vsf runo -> vsf_xy runo

    for(i = num_cvar + num_var_unordered + num_var_ordered; i < (num_cvar + num_uvar + num_var_ordered); i++, l++){
      vsf_xy[l] = vector_scale_factor[i];
    }

    // copy vsf vuno -> vsf_xy vuno

    for(i = num_cvar; i < (num_cvar + num_var_unordered); i++, l++){
      vsf_xy[l] = vector_scale_factor[i];
    }

    // copy vsf rord -> vsf_xy rord

    for(i = num_cvar + num_uvar + num_var_ordered; i < num_all_var; i++, l++){
      vsf_xy[l] = vector_scale_factor[i];
    }

    // copy vsf vord -> vsf_xy vord

    for(i = num_cvar + num_var_unordered; i < (num_cvar + num_var_unordered + num_var_ordered); i++, l++){
      vsf_xy[l] = vector_scale_factor[i];
    }

    // vsf xy DONE
  }

  if(vsf_x != NULL){
    // set up x bws

    // vsf rcon -> vsf_x rcon
    for(i = 0, l = 0; i < num_reg_continuous; i++, l++){
      vsf_x[l] = vector_scale_factor[i];
    }

    // copy vsf runo -> vsf_x runo

    for(i = num_cvar + num_var_unordered + num_var_ordered; i < num_all_var; i++, l++){
      vsf_x[l] = vector_scale_factor[i];
    }

  }

  if(vsf_y != NULL){
    // set up y bws
    for(i = num_reg_continuous, l = 0; i < (num_cvar + num_var_unordered + num_var_ordered); i++, l++){
      vsf_y[l] = vector_scale_factor[i];
    }
  }

  // copy num_categories arrays
  if(nc_xy != NULL){
    for(i = num_var_unordered + num_var_ordered, l = 0; i < (num_var_unordered + num_var_ordered + num_reg_unordered); i++, l++){
      nc_xy[l] = num_categories[i];
    }

    for(i = 0; i < num_var_unordered; i++, l++){
      nc_xy[l] = num_categories[i];
    }

    for(i = num_var_unordered + num_var_ordered + num_reg_unordered; i < num_catvar; i++, l++){
      nc_xy[l] = num_categories[i];
    }

    for(i = num_var_unordered; i < (num_var_unordered + num_var_ordered); i++, l++){
      nc_xy[l] = num_categories[i];
    }
  }

  if(nc_x != NULL){
    for(i = num_var_unordered + num_var_ordered, l = 0; i < num_catvar; i++, l++){
      nc_x[l] = num_categories[i];
    }    
  }

  if(nc_y != NULL){
    for(i = 0, l = 0; i < (num_var_unordered + num_var_ordered); i++, l++){
      nc_y[l] = num_categories[i];
    }    
  }

  // copy/fix mcv arrays

  if(mcv_xy != NULL){
    for(i = num_var_unordered + num_var_ordered, l = 0; i < (num_var_unordered + num_var_ordered + num_reg_unordered); i++, l++){
      for(j = 0; j < num_categories[i]; j++)
        mcv_xy[l][j] = matrix_categorical_vals[i][j];
    }

    for(i = 0; i < num_var_unordered; i++, l++){
      for(j = 0; j < num_categories[i]; j++)
        mcv_xy[l][j] = matrix_categorical_vals[i][j];
    }

    for(i = num_var_unordered + num_var_ordered + num_reg_unordered; i < num_catvar; i++, l++){
      for(j = 0; j < num_categories[i]; j++)
        mcv_xy[l][j] = matrix_categorical_vals[i][j];
    }

    for(i = num_var_unordered; i < (num_var_unordered + num_var_ordered); i++, l++){
      for(j = 0; j < num_categories[i]; j++)
        mcv_xy[l][j] = matrix_categorical_vals[i][j];
    }
  }

  if(mcv_x != NULL){
    for(i = num_var_unordered + num_var_ordered, l = 0; i < num_catvar; i++, l++){
      for(j = 0; j < num_categories[i]; j++)
        mcv_x[l][j] = matrix_categorical_vals[i][j];
    }    
  }

  if(mcv_y != NULL){
    for(i = 0, l = 0; i < (num_var_unordered + num_var_ordered); i++, l++){
      for(j = 0; j < num_categories[i]; j++)
        mcv_y[l][j] = matrix_categorical_vals[i][j];

    }    
  }
}

void np_kernelop_xy(const int kernel_var_continuous,
		    const int kernel_var_unordered,
		    const int kernel_var_ordered,
		    const int kernel_reg_continuous,
		    const int kernel_reg_unordered,
		    const int kernel_reg_ordered,
		    const int operator_var,
		    const int operator_reg,
		    const int num_var_unordered,
		    const int num_var_ordered,
		    const int num_var_continuous,
		    const int num_reg_unordered,
		    const int num_reg_ordered,
		    const int num_reg_continuous,
		    int * kernel_cx,
		    int * kernel_ux,
		    int * kernel_ox,
		    int * kernel_cy,
		    int * kernel_uy,
		    int * kernel_oy,
		    int * kernel_cxy,
		    int * kernel_uxy,
		    int * kernel_oxy,
		    int * operator_x,
		    int * operator_y,
		    int * operator_xy){

  const int num_reg = num_reg_continuous + num_reg_unordered + num_reg_ordered;
  const int num_var = num_var_continuous + num_var_unordered + num_var_ordered;
  const int num_all = num_reg+num_var;

  const int num_cvar = num_var_continuous + num_reg_continuous;
  const int num_uvar = num_var_unordered + num_reg_unordered;
  const int num_ovar = num_var_ordered + num_reg_ordered;

  int i;
  // x data
  
  if(kernel_cx != NULL)
    for(i = 0; i < num_reg_continuous; i++)
      kernel_cx[i] = kernel_reg_continuous;

  if(kernel_ux != NULL)
    for(i = 0; i < num_reg_unordered; i++)
      kernel_ux[i] = kernel_reg_unordered;

  if(kernel_ox != NULL)
    for(i = 0; i < num_reg_ordered; i++)
      kernel_ox[i] = kernel_reg_ordered;

  // y data
  if(kernel_cy != NULL)
    for(i = 0; i < num_var_continuous; i++)
      kernel_cy[i] = kernel_var_continuous;

  if(kernel_uy != NULL)
    for(i = 0; i < num_var_unordered; i++)
      kernel_uy[i] = kernel_var_unordered;

  if(kernel_oy != NULL)
    for(i = 0; i < num_var_ordered; i++)
      kernel_oy[i] = kernel_var_ordered;
  
  // xy data
  if(kernel_cxy != NULL){
    for(i = 0; i < num_reg_continuous; i++)
      kernel_cxy[i] = kernel_reg_continuous;

    for(i = num_reg_continuous; i < num_cvar; i++)
      kernel_cxy[i] = kernel_var_continuous;
  }

  if(kernel_uxy != NULL){
    for(i = 0; i < num_reg_unordered; i++)
      kernel_uxy[i] = kernel_reg_unordered;

    for(i = num_reg_unordered; i < num_uvar; i++)
      kernel_uxy[i] = kernel_var_unordered;
  }

  if(kernel_oxy != NULL){
    for(i = 0; i < num_reg_ordered; i++)
      kernel_oxy[i] = kernel_reg_ordered;

    for(i = num_reg_ordered; i < num_ovar; i++)
      kernel_oxy[i] = kernel_var_ordered;
  }

  // operators
  if(operator_x != NULL){
    for(i = 0; i < num_reg; i++)
      operator_x[i] = operator_reg;
  }

  if(operator_y != NULL){
    for(i = 0; i < num_var; i++)
      operator_y[i] = operator_var;
  }

  if(operator_xy != NULL){
    for(i = 0; i < num_reg_continuous; i++)
      operator_xy[i] = operator_reg;

    for(; i < num_cvar; i++)
      operator_xy[i] = operator_var;

    for(; i < num_cvar+num_reg_unordered; i++)
      operator_xy[i] = operator_reg;

    for(; i < num_cvar+num_uvar; i++)
      operator_xy[i] = operator_var;

    for(; i < num_cvar+num_uvar+num_reg_ordered; i++)
      operator_xy[i] = operator_reg;

    for(; i < num_all; i++)
      operator_xy[i] = operator_var;
  }

}
