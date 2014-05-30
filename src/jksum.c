/* Copyright (C) J. Racine, 1995-2001 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <errno.h>

#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

#include "headers.h"
#include "matrix.h"

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

#ifdef RCSID
static char rcsid[] = "$Id: jksum.c,v 1.16 2006/11/02 16:56:49 tristen Exp $";
#endif

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

extern int * ipt_lookup_extern_X;
extern int * ipt_lookup_extern_Y;
extern int * ipt_lookup_extern_XY;

extern int *num_categories_extern_XY;
extern int *num_categories_extern_X;
extern int *num_categories_extern_Y;

extern double ** matrix_categorical_vals_extern_X;
extern double ** matrix_categorical_vals_extern_Y;
extern double ** matrix_categorical_vals_extern_XY;

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
#ifdef MPI2
		MPI_Barrier(comm[1]);
		MPI_Finalize();
#endif
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

	MPI_Gather(kernel_sum, stride, MPI_DOUBLE, kernel_sum, stride, MPI_DOUBLE, 0, comm[1]);
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
  return (same_cat)?(1.0-lambda):lambda/((double)c-1.0);
}

double np_score_uaa(const int same_cat,const double lambda, const int c){
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
  return (x == y) ? -1.0 : (0.5*ipow(lambda, (int)fabs(x-y))*(fabs(x-y)/lambda - 2.0));
}

double np_oli_racine(const double x, const double y, const double lambda, const double cl, const double ch){
  return R_pow_di(lambda, (int)fabs(x-y));
}

double np_score_oli_racine(const double x, const double y, const double lambda, const double cl, const double ch){
  return (fabs(x-y)*R_pow_di(lambda, (int)fabs(x-y)-1));
}

double np_onli_racine(const double x, const double y, const double lambda, const double cl, const double ch){
  return R_pow_di(lambda, (int)fabs(x-y))*(1.0 - lambda)/(1.0 + lambda);
}

double np_score_onli_racine(const double x, const double y, const double lambda, const double cl, const double ch){
  const int cxy = (int)fabs(x-y);
  return ((cxy != 0) || (lambda != 0.0)) ? R_pow_di(lambda, cxy - 1)*(cxy*(1.0 - lambda*lambda) - 2.0 *lambda) : -2.0;
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
  const int cxy = (int)fabs(x-y);
  const double lnorm = (1.0 - lambda)/(1.0 + lambda);
  const double l2 = lambda*lambda;
  return lnorm*lnorm*R_pow_di(lambda, cxy)*((1.0 + l2)/(1.0 - l2) + cxy);

}

double np_econvol_owang_van_ryzin(const double x, const double y, const double lambda, const double cl, const double ch){
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

double np_cdf_owang_van_ryzin(const double y, const double x, const double lambda, const double cl, const double ch){
  if(x == y) return 1.0 - 0.5*lambda;
  const int cxy = (int)fabs(x-y);
  const double gee = R_pow_di(lambda, cxy);
  return (x < y) ? 0.5*gee : (1.0 - gee);
}

double np_cdf_oli_racine(const double y, const double x, const double lambda, const double cl, const double ch){
  const int xh = (x > ch) ? ch : x;
  const int cxy = (int)fabs(xh-y);
  const double gee = R_pow_di(lambda, cxy)/(1.0-lambda);
  if(x < y){
    return (x < cl) ? 0.0 : gee*(1.0-R_pow_di(lambda,(int)(x-cl+1)));
  } else {
    return (1.0 + lambda - R_pow_di(lambda,(int)(y-cl+1)))/(1.0 - lambda) - lambda*gee;
  }
}

double np_cdf_onli_racine(const double y, const double x, const double lambda, const double cl, const double ch){
  const int cxy = (int)fabs(x-y);
  const double gee = R_pow_di(lambda, cxy)/(1.0+lambda);
  return (x < y) ? gee : 1.0 - lambda*gee;
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
  const double a = sqrt(2);
  const double hx2 = hx*hx;
  const double hy2 = hy*hy;
  const double hxy2 = hx2+hy2;
  const double x2 = x*x;
  const double y2 = y*x;
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
  const double hs6 = hs4*hs2;
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
  const double hx7 = hx4*hx3;
  const double hx8 = hx4*hx4;
  const double hx9 = hx5*hx4;
  const double hx10 = hx6*hx4;
  const double hx11 = hx7*hx4;
  const double hx12 = hx8*hx4;
  const double hx13 = hx9*hx4;
  const double hx14 = hx10*hx4;
  const double hx15 = hx11*hx4;
  const double hx16 = hx12*hx4;
  const double hx17 = hx13*hx4;

  const double hy2 = hy*hy;
  const double hy3 = hy2*hy;
  const double hy4 = hy2*hy2;
  const double hy5 = hy3*hy2;
  const double hy6 = hy3*hy3;
  const double hy7 = hy4*hy3;
  const double hy8 = hy4*hy4;
  const double hy9 = hy5*hy4;
  const double hy10 = hy6*hy4;
  const double hy11 = hy7*hy4;
  const double hy12 = hy8*hy4;
  const double hy13 = hy9*hy4;
  const double hy14 = hy10*hy4;
  const double hy15 = hy11*hy4;
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
  const double hx7 = hx4*hx3;
  const double hx8 = hx4*hx4;
  const double hx9 = hx5*hx4;
  const double hx10 = hx6*hx4;
  const double hx11 = hx7*hx4;
  const double hx12 = hx8*hx4;
  const double hx13 = hx9*hx4;
  const double hx14 = hx10*hx4;
  const double hx15 = hx11*hx4;
  const double hx16 = hx12*hx4;
  const double hx17 = hx13*hx4;

  const double hy2 = hy*hy;
  const double hy3 = hy2*hy;
  const double hy4 = hy2*hy2;
  const double hy5 = hy3*hy2;
  const double hy6 = hy3*hy3;
  const double hy7 = hy4*hy3;
  const double hy8 = hy4*hy4;
  const double hy9 = hy5*hy4;
  const double hy10 = hy6*hy4;
  const double hy11 = hy7*hy4;
  const double hy12 = hy8*hy4;
  const double hy13 = hy9*hy4;
  const double hy14 = hy10*hy4;
  const double hy15 = hy11*hy4;
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
                   const int do_score){

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

  double * kbuf = NULL;

  kbuf = (double *)malloc(num_xt*sizeof(double));

  assert(kbuf != NULL);

  if(xl == NULL){
    for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
      const double kn = k[KERNEL]((x-xt[i])*sgn/h);

      result[i] = xw[j]*kn;
      kbuf[i] = kn;
      
      if(do_perm)
        p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*k[P_KERNEL]((x-xt[i])*sgn/h)*(do_score ? ((xt[i]-x)*sgn/h) : 1.0);

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
        const double kn = k[KERNEL]((x-xt[i])*sgn/h);

        result[i] = xw[j]*kn;
        kbuf[i] = kn;
      }
    }

    if(do_perm){
      for (int m = 0; m < p_xl->n; m++){
        const int istart = p_xl->istart[m];
        const int nlev = p_xl->nlev[m];
        for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
          p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*k[P_KERNEL]((x-xt[i])*sgn/h)*(do_score ? ((xt[i]-x)*sgn/h) : 1.0);
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

  free(kbuf);
}

void np_ckernelv(const int KERNEL, 
                 const double * const xt, const int num_xt, 
                 const int do_xw,
                 const double x, const double h, 
                 double * const result,
                 const XL * const xl,
                 const int swap_xxt){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i,j; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  const double sgn = swap_xxt ? -1.0 : 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

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

  if(xl == NULL)
    for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
      if(xw[j] == 0.0) continue;
      result[i] = xw[j]*k[KERNEL]((x-xt[i])*sgn/h);
    }
  else{
    for (int m = 0; m < xl->n; m++){
      const int istart = xl->istart[m];
      const int nlev = xl->nlev[m];
      for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
        if(xw[j] == 0.0) continue;
        result[i] = xw[j]*k[KERNEL]((x-xt[i])*sgn/h);
      }
    }
  }

}

void np_convol_ckernelv(const int KERNEL, 
                        const double * const xt, const int num_xt, 
                        const int do_xw,
                        const double x, 
                        double * xt_h, 
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
    if(xw[j] == 0.0) continue;
    result[i] = xw[j]*k[KERNEL](x,xt[i],h,xt_h[i])/ipow(xt_h[i], power);
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
                   const int do_ocg){

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

  double * kbuf = NULL;

  kbuf = (double *)malloc(num_xt*sizeof(double));

  assert(kbuf != NULL);

  if(xl == NULL){
    for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
      const double kn = k[KERNEL]((xt[i]==x), lambda, ncat);

      result[i] = xw[j]*kn;
      kbuf[i] = kn;

      const int iscat = (swap_xxt && do_ocg) ? (cat == x) : (xt[i] == ex);

      p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*k[P_KERNEL](iscat, lambda, ncat);
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
        const double kn = k[KERNEL]((xt[i]==x), lambda, ncat);

        result[i] = xw[j]*kn;
        kbuf[i] = kn;
      }
    }

    for (int m = 0; m < p_xl->n; m++){
      const int istart = p_xl->istart[m];
      const int nlev = p_xl->nlev[m];
      for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
        const int iscat = (swap_xxt && do_ocg) ? (cat == x) : (xt[i] == ex);
        p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*k[P_KERNEL](iscat, lambda, ncat);
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

  free(kbuf);
}

void np_ukernelv(const int KERNEL, 
                 const double * const xt, const int num_xt, 
                 const int do_xw,
                 const double x, const double lambda, const int ncat,
                 double * const result,
                 const XL * const xl){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i; 
  int j, bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double (* const k[])(int, double, int) = { np_uaa, np_unli_racine,
                                             np_econvol_uaa, np_econvol_unli_racine };

  if(xl == NULL){
    for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
      if(xw[j] == 0.0) continue;
      result[i] = xw[j]*k[KERNEL]((xt[i]==x), lambda, ncat);
    }
  } else {
    for (int m = 0; m < xl->n; m++){
      const int istart = xl->istart[m];
      const int nlev = xl->nlev[m];
      for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
        if(xw[j] == 0.0) continue;
        result[i] = xw[j]*k[KERNEL]((xt[i]==x), lambda, ncat);
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
  int j, bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  if(!swap_xxt){
    for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
      if(xw[j] == 0.0) continue;
      result[i] = xw[j]*kernel_ordered_convolution(KERNEL, xt[i], x, lambda, ncat, cat);
    }
  } else {
    for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
      if(xw[j] == 0.0) continue;
      result[i] = xw[j]*kernel_ordered_convolution(KERNEL, x, xt[i], lambda, ncat, cat);
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
                   const int swapped_index){

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
    np_owang_van_ryzin, np_oli_racine, np_onli_racine, 
    np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine,
    np_score_owang_van_ryzin, np_score_oli_racine, np_score_onli_racine,
    np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine
  };

  double * kbuf = NULL;

  kbuf = (double *)malloc(num_xt*sizeof(double));

  assert(kbuf != NULL);

  double s_cat = 0.0;

  if((!swap_xxt) && do_ocg){
    s_cat = cats[abs(swapped_index - 1)];
  }
  
  const double cl = (cats != NULL)? cats[0] : 0.0;
  const double ch = (cats != NULL)? cats[ncat - 1] : 0.0;

    if(xl == NULL){
      for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
        const double cat = do_ocg ? (swap_xxt ? cats[abs(ordered_indices[i] - 1)] : s_cat) : 0.0;
        const double c1 = swap_xxt ? x : xt[i];
        const double c2 = swap_xxt ? xt[i] : x;
        const double c3 = do_ocg ? cat : (swap_xxt ? xt[i] : x);

        const double kn = k[KERNEL](c1, c2, lambda, cl, ch);

        result[i] = xw[j]*kn;
        kbuf[i] = kn;

        p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*k[P_KERNEL](c1, c3, lambda, cl, ch);
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

          const double kn = k[KERNEL](c1, c2, lambda, cl, ch);

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

          p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*k[P_KERNEL](c1, c3, lambda, cl, ch);
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

  free(kbuf);
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
  int j, bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double (* const k[])(double, double, double, double, double) = { 
    np_owang_van_ryzin, np_oli_racine, np_onli_racine, 
    np_econvol_owang_van_ryzin, np_onull, np_econvol_onli_racine,
    np_onull, np_onull, np_onull,
    np_cdf_owang_van_ryzin, np_cdf_oli_racine, np_cdf_onli_racine
  };

  const double cl = (cats != NULL)? cats[0] : 0.0;
  const double ch = (cats != NULL)? cats[ncat - 1] : 0.0;

  if(!swap_xxt){
    if(xl == NULL){
      for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
        if(xw[j] == 0.0) continue;
        result[i] = xw[j]*k[KERNEL](xt[i], x, lambda, cl, ch);
      }
    } else {
      for (int m = 0; m < xl->n; m++){
        const int istart = xl->istart[m];
        const int nlev = xl->nlev[m];
        for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
          if(xw[j] == 0.0) continue;
          result[i] = xw[j]*k[KERNEL](xt[i], x, lambda, cl, ch);
        }
      }
    }
  } else {
    if(xl == NULL){
      for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
        if(xw[j] == 0.0) continue;
        result[i] = xw[j]*k[KERNEL](x, xt[i], lambda, cl, ch);
      }
    } else {
      for (int m = 0; m < xl->n; m++){
        const int istart = xl->istart[m];
        const int nlev = xl->nlev[m];
        for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
          if(xw[j] == 0.0) continue;
          result[i] = xw[j]*k[KERNEL](x, xt[i], lambda, cl, ch);
        }
      }
    }
  }
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

  if (do_leave_one_out) {
    temp = weights[which_k];
    weights[which_k] = 0.0;
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
                result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*ipow(weights[k]/db, kpow);
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
                result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*ipow(weights[k]/db, kpow);
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
                result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*ipow(weights[k]/db, kpow);
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
                result[k*kstride+j*max_B+i] += pmat_A[j][k*have_A]*pmat_B[i][k*have_B]*ipow(weights[k]/db, kpow);
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
                  result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*ipow(weights[k]/db, kpow);
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
                  result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*ipow(weights[k]/db, kpow);
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
                  result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*ipow(weights[k]/db, kpow);
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
                  result[k*kstride+j*max_B+i] += pmat_A[j][k*have_A]*pmat_B[i][k*have_B]*ipow(weights[k]/db, kpow);
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
}



//Warning: the MPI operations used by this function assume that the
//receiving buffers, i.e. weighted_sum have enough space allocated such that a write
//consisting of nproc*stride will not segfault. It is up to the caller
//to ensure that this is the case.

// this will be fixed by using the mpi*v functions

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
double * const kw){
  
  /* This function takes a vector Y and returns a kernel weighted
     leave-one-out sum. By default Y should be a vector of ones
     (simply compute the kernel sum). This function will allow users
     to `roll their own' with mixed data leave-one out kernel sums. */

  /* Declarations */

  int i, ii, j, kk, k, l, mstep, js, je, num_obs_eval_alloc, sum_element_length, ip;
  int do_psum, swap_xxt;

  // USED TO default to -1
  int * permutation_kernel = NULL;
  int doscoreocg = do_score || do_ocg;
  int do_perm = permutation_operator != OP_NOOP; 

  
  const int no_bpso = (NULL == bpso);

  int p_nvar;

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
  assert(igatherv != NULL);
  idisplsv = (int *)malloc(iNum_Processors*sizeof(int));
  assert(idisplsv != NULL);
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
  const int ws_step = is_adaptive? 0 :
                                 (MAX(ncol_Y, 1) * MAX(ncol_W, 1));

  double *lambda, **matrix_bandwidth, **matrix_alt_bandwidth = NULL, *m = NULL;
  double *tprod, dband, *ws, * p_ws, * tprod_mp = NULL, * p_dband = NULL;

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

  if (num_obs_eval == 0) {
    return(KWSNP_ERR_NOEVAL);
  }

  do_psum = BANDWIDTH_reg == BW_ADAP_NN;
  swap_xxt = BANDWIDTH_reg == BW_ADAP_NN;
  /* Allocate memory for objects */

  mstep = (BANDWIDTH_reg==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_train:1);


  if(bandwidth_provided){
    if(BANDWIDTH_reg == BW_GEN_NN)
      matrix_bandwidth = matrix_bw_eval;
    else if (is_adaptive){
      if (any_convolution){
        matrix_alt_bandwidth = matrix_bw_eval;        
      } 
      matrix_bandwidth = matrix_bw_train;
    } else {
      assert(0);
    }
    lambda = lambda_pre;
  } else {
    matrix_bandwidth = alloc_tmatd(mstep, num_reg_continuous);  
    lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);
  } 


  tprod = alloc_vecd((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_eval:num_obs_train);

  sum_element_length = MAX(ncol_Y, 1) * 
    MAX(ncol_W, 1);

  /* assert(!(BANDWIDTH_reg == BW_ADAP_NN)); */
  /* Conduct the estimation */

  /* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

  bpow = (int *) malloc(num_reg_continuous*sizeof(int));
  assert(bpow != NULL);

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

    assert(tprod_mp != NULL);

    p_dband = (double *)malloc(p_nvar*sizeof(double));

    assert(p_dband != NULL);

    p_ipow = (bandwidth_divide ? 1 : 0) + ((permutation_operator == OP_DERIVATIVE) ? 1 : ((permutation_operator == OP_INTEGRAL) ? -1 : 0));

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

      free(lambda);
      free_tmat(matrix_bandwidth);
      free(tprod);

      return(KWSNP_ERR_BADBW);
    }
  }

  if(!bandwidth_provided && (is_adaptive && any_convolution)){ // need additional bandwidths 
    matrix_alt_bandwidth = alloc_tmatd(num_obs_eval, num_reg_continuous);  

    // this is a bug
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

      free(lambda);
      free_tmat(matrix_bandwidth);
      free_tmat(matrix_alt_bandwidth);
      free(tprod);

      return(KWSNP_ERR_BADBW);
    }

  }
  
  if (leave_one_out && drop_one_train) {
    REprintf("\n error, leave one out estimator and drop-one estimator can't be enabled simultaneously");
    free(lambda);
    free_tmat(matrix_bandwidth);
    free(tprod);
    return(KWSNP_ERR_BADINVOC);

  }

  if ((num_obs_train < (num_obs_eval + leave_one_out_offset)) && leave_one_out){
    
    REprintf("\nnumber of training points must be >= number of evaluation points to use the 'leave one out' estimator");
    free(lambda);
    free_tmat(matrix_bandwidth);
    free(tprod);
    return(KWSNP_ERR_BADINVOC);
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

  const int leave_or_drop = leave_one_out || (drop_one_train && (BANDWIDTH_reg != BW_ADAP_NN));
  if(drop_one_train && (BANDWIDTH_reg != BW_ADAP_NN)) lod = drop_which_train;

    /* do sums */
  for(j=js; j <= je; j++, ws += ws_step, p_ws += ws_step){
    R_CheckUserInterrupt();

    dband = 1.0;

    for (ii = 0; ii < p_nvar; ii++)
      p_dband[ii] = 1.0;

    // if we are consistently dropping one obs from training data, and we are doing a parallel sum, then that means
    // we need to skip here

    if(leave_one_out) lod = j + leave_one_out_offset;

    if (num_reg_continuous > 0){
      m = matrix_bandwidth[0];
      if (BANDWIDTH_reg != BW_FIXED)
        m += j;
    }

    // do a hail mary, then generate the interaction list
    // anything but a fixed bandwidth is not yet supported
    // that includes convolutions 

    if(np_ks_tree_use){

      if(kw != NULL)
        for(i = 0; i < num_xt; i++)
          tprod[i] = 0.0;

      double bb[kdt->ndim*2];

      // reset the interaction node list
      xl.n = 0;
      if(!do_partial_tree){
        for(i = 0; i < num_reg_continuous; i++){
          const double sf = (BANDWIDTH_reg != BW_FIXED) ? (*(m+i*mstep)):m[i];
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
      } else {
        for(i = 0; i < num_reg_continuous; i++){
          const double sf = (BANDWIDTH_reg != BW_FIXED) ? (*(m+i*mstep)):m[i];
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

      if(do_perm && (permutation_operator == OP_INTEGRAL)){
        for(ii = 0, k = 0; ii < num_reg_continuous; ii++){
          if(bpso[ii]){
            // reset the interaction node list
            p_pxl[k].n = 0;

            if(!do_partial_tree){
              for(i = 0; i < num_reg_continuous; i++){
                const int knp = (i == ii) ? permutation_kernel[i] : KERNEL_reg_np[i];
                if(!is_adaptive){
                  bb[2*i] = -cksup[knp][1];
                  bb[2*i+1] = -cksup[knp][0];
                }else{
                  bb[2*i] = cksup[knp][0];
                  bb[2*i+1] = cksup[knp][1];
                }
                const double sf = (BANDWIDTH_reg != BW_FIXED) ? (*(m+i*mstep)):m[i];
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
                const double sf = (BANDWIDTH_reg != BW_FIXED) ? (*(m+i*mstep)):m[i];
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

    }
    /* continuous first */

    /* for the first iteration, no weights */
    /* for the rest, the accumulated products are the weights */
    for(i = 0, l = 0, ip = 0, k = 0; i < num_reg_continuous; i++, l++, m += mstep, ip += do_perm){
      if((BANDWIDTH_reg != BW_ADAP_NN) || (operator[l] != OP_CONVOLUTION)){
        if(p_nvar == 0){
          np_ckernelv(KERNEL_reg_np[i], xtc[i], num_xt, l, xc[i][j], *m, tprod, pxl, swap_xxt);
        } else {
          np_p_ckernelv(KERNEL_reg_np[i], (do_perm ? permutation_kernel[i] : KERNEL_reg_np[i]), k, p_nvar, xtc[i], num_xt, l, xc[i][j], *m, tprod, tprod_mp, pxl, p_pxl+k, swap_xxt, bpso[l], do_score);
        }
      }
      else
        np_convol_ckernelv(KERNEL_reg[i], xtc[i], num_xt, l, xc[i][j], 
                           matrix_alt_bandwidth[i], *m, tprod, bpow[i]);
      dband *= ipow(*m, bpow[i]);

      if(do_perm){
        for(ii = 0, kk = 0; ii < num_reg_continuous; ii++){
          if(bpso[ii]){
            if (i != ii){
              p_dband[kk] *= ipow(*m, bpow[i]);              
            } else {
              p_dband[kk] *= ipow(*m, p_ipow);
              if(((BANDWIDTH_reg == BW_FIXED) && (int_LARGE_SF == 0) && do_score)){
                p_dband[kk] *= vector_scale_factor[ii]/(*m);
              }
            }
            kk++;
          }
        }
      }
      k += bpso[l];
    }


    for(ii = k; ii < p_nvar; ii++){
      p_dband[ii] = dband;
    }

    /* unordered second */

    for(i=0; i < num_reg_unordered; i++, l++, ip += doscoreocg){
      if(doscoreocg){
        np_p_ukernelv(KERNEL_unordered_reg_np[i], ps_ukernel[i], k, p_nvar, xtu[i], num_xt, l, xu[i][j], 
                      lambda[i], num_categories[i], matrix_categorical_vals[i][0], tprod, tprod_mp, pxl, p_pxl + k, swap_xxt, (bpso[l] ? do_ocg : 0));
      } else {
        np_ukernelv(KERNEL_unordered_reg_np[i], xtu[i], num_xt, l, xu[i][j], 
                    lambda[i], num_categories[i], tprod, pxl);
      }
      k += bpso[l];
    }

    /* ordered third */
    for(i=0; i < num_reg_ordered; i++, l++, ip += doscoreocg){
      if(!doscoreocg){
        if(ps_ok_nli || (operator[l] != OP_CONVOLUTION)){
          np_okernelv(KERNEL_ordered_reg_np[i], xto[i], num_xt, l,
                      xo[i][j], lambda[num_reg_unordered+i], 
                      (matrix_categorical_vals != NULL) ? matrix_categorical_vals[i+num_reg_unordered] : NULL, 
                      (num_categories != NULL) ? num_categories[i+num_reg_unordered] : 0,
                      tprod, pxl, swap_xxt);      
        } else {
          np_convol_okernelv(KERNEL_ordered_reg[i], xto[i], num_xt, l,
                             xo[i][j], lambda[num_reg_unordered+i], 
                             num_categories[i+num_reg_unordered],
                             matrix_categorical_vals[i+num_reg_unordered],
                             tprod, swap_xxt);
        }
      } else {
        np_p_okernelv(KERNEL_ordered_reg_np[i], ps_okernel[i], k, p_nvar, xto[i], num_xt, l,
                      xo[i][j], lambda[num_reg_unordered+i], 
                      (matrix_categorical_vals != NULL) ? matrix_categorical_vals[i+num_reg_unordered] : NULL, 
                      (num_categories != NULL) ? num_categories[i+num_reg_unordered] : 0,
                      tprod, tprod_mp, pxl, p_pxl + k, swap_xxt, (bpso[l] ? do_ocg : 0), matrix_ordered_indices[i], (swap_xxt ? 0 : matrix_ordered_indices[i][j]));
      }
      k += bpso[l];
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

    if(kw != NULL){ 
      // if using adaptive bandwidths, kw is returned transposed
      if(bandwidth_divide_weights)
        for(i = 0; i < num_xt; i++)
          kw[j*num_xt + i] = tprod[i]/dband;
      else
        for(i = 0; i < num_xt; i++)
          kw[j*num_xt + i] = tprod[i];
    }
    
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

    if(kw != NULL){
      MPI_Allgather(MPI_IN_PLACE, stride * num_xt, MPI_DOUBLE, kw, stride * num_xt, MPI_DOUBLE, comm[1]);          
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

  if(no_bpso)
    free(bpso);
  
  return(0);
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

  k = (int)uki;
  j = (int)uji;
  i = (int)uii;

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
  if(!(blk_xj != NULL)) 
    error("!(blk_xj != NULL)");

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
    sum_ker_marginalf[j] =  NZD(sum_ker_marginalf[j]);
    *cv += (sum_ker_convolf[j]/sum_ker_marginalf[j]-2.0*sum_kerf[j])/sum_ker_marginalf[j];
  }
#else
  for(j = 0; j < num_obs; j++){
    /*    if(sum_ker_marginal[j] <= 0.0){
      *cv = DBL_MAX;
      break;
      } jracine 16/05/10 */
    sum_ker_marginal[j] =  NZD(sum_ker_marginal[j]);
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

// consider the leave one out local linear estimator used in cross-validation:
// g_{-i} = (sum_{j!=i}(k_{ji}*q_{ji}*t(q_{ji})))^(-1)*(sum_{j!=i}(k_{ji}*q_{ji}*y_{j})
// the expression k_{ji}*q_{ji}*t(q_{ji}) is nearly symmetric with respect to interchange of i and j:
// k_{ji}*q_{ji}*t(q_{ji}) = k_{ij}*(q_{ij}*t(q_{ij}) - 2*(q_{ij}*t(l) + l*t(q_{ij})) + 4*l*t(l))
// where t(q_{ij}) = (1 (x_{i}-x_{j})) and t(l) = (1 0) 
// we take advantage of both the quasi-parity, and 
// of the algebraic transpose symmetry of the aforementioned expression to gain a factor of ~ 4 speed-up


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

  // note that mean has 2*num_obs allocated for npksum
  int i, j, l, sf_flag = 0, num_obs_eval_alloc, tsf;

  double cv = 0.0;
  double * lambda = NULL, * vsf = NULL;
  double ** matrix_bandwidth = NULL;

  double aicc = 0.0;
  double traceH = 0.0;

  int * operator = NULL;
  int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;

  const int leave_one_out = (bwm == RBWM_CVLS)?1:0;

  operator = (int *)malloc(sizeof(int)*(num_reg_continuous+num_reg_unordered+num_reg_ordered));

  for(i = 0; i < (num_reg_continuous+num_reg_unordered+num_reg_ordered); i++)
    operator[i] = OP_NORMAL;

  kernel_c = (int *)malloc(sizeof(int)*num_reg_continuous);

  for(i = 0; i < num_reg_continuous; i++)
    kernel_c[i] = KERNEL_reg;

  kernel_u = (int *)malloc(sizeof(int)*num_reg_unordered);

  for(i = 0; i < num_reg_unordered; i++)
    kernel_u[i] = KERNEL_unordered_reg;

  kernel_o = (int *)malloc(sizeof(int)*num_reg_ordered);

  for(i = 0; i < num_reg_ordered; i++)
    kernel_o[i] = KERNEL_ordered_reg;

#ifdef MPI2
    int stride = MAX((int)ceil((double) num_obs / (double) iNum_Processors),1);
    num_obs_eval_alloc = stride*iNum_Processors;
#else
    num_obs_eval_alloc = num_obs;
#endif

    int ks_tree_use = (int_TREE_X == NP_TREE_TRUE) && (!((BANDWIDTH_reg == BW_ADAP_NN) && (int_ll == LL_LL)));

  // Allocate memory for objects 

  lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);
  matrix_bandwidth = alloc_matd(num_obs,num_reg_continuous);

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
    
    free(lambda);
    free_mat(matrix_bandwidth,num_reg_continuous);
    
    return(DBL_MAX);
  }


  if(bwm == RBWM_CVAIC){
    // compute normalisation constant

    // workaround for bwscaling = TRUE
    // really just want to get the full product kernel evaluated at zero
    tsf = int_LARGE_SF;

    int_LARGE_SF = 1;

    kernel_weighted_sum_np(kernel_c,
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
                           (BANDWIDTH_reg == BW_ADAP_NN)?1:0, // bandwidth_divide = FALSE when not adaptive
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
                           0,NULL,NULL,NULL,
                           num_categories,
                           NULL,
                           NULL,
                           &aicc,
                           NULL, // no permutations
                           NULL); // do not return kernel weights
    int_LARGE_SF = tsf;

    //fprintf(stderr,"\n%e\n",aicc);
  }

  // Conduct the estimation 

  if(int_ll == LL_LC) { // local constant
    // Nadaraya-Watson
    // Generate bandwidth vector given scale factors, nearest neighbors, or lambda 

    double * lc_Y[2];
    double * mean = (double *)malloc(2*num_obs_eval_alloc*sizeof(double));

    lc_Y[0] = vector_Y;
      
    lc_Y[1] = (double *)malloc(num_obs*sizeof(double));
    for(int ii = 0; ii < num_obs; ii++)
      lc_Y[1][ii] = 1.0;

    kernel_weighted_sum_np(kernel_c,
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
                           0,NULL,NULL,NULL,
                           num_categories,
                           NULL,
                           NULL,
                           mean,
                           NULL, // no permutations
                           NULL); // do not return kernel weights

    // every even entry in mean is sum(y*kij)
    // every odd is sum(kij)

    for(int ii = 0; ii < num_obs; ii++){
      const int ii2 = 2*ii;
      const double sk = copysign(DBL_MIN, mean[ii2+1]) + mean[ii2+1];
      const double dy = vector_Y[ii]-mean[ii2]/sk;
      cv += dy*dy;
      if(bwm == RBWM_CVAIC)
        traceH += aicc/sk;

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

    double * kwm = (double *)malloc(nrcc22*num_obs_eval_alloc*sizeof(double));

    for(int ii = 0; ii < nrcc22*num_obs; ii++)
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
      if(ks_tree_use){
        if((j % iNum_Processors) == 0){
          if((j+my_rank) < (num_obs)){
            for(l = 0; l < num_reg_continuous; l++){
          
              for(i = 0; i < num_obs; i++){
                XTKX[l+2][i] = matrix_X_continuous[l][i]-matrix_X_continuous[l][j+my_rank];
              }
              TCON[l][0] = matrix_X_continuous[l][j+my_rank]; // temporary storage
            }


            for(l = 0; l < num_reg_unordered; l++)
              TUNO[l][0] = matrix_X_unordered[l][j+my_rank];

            for(l = 0; l < num_reg_ordered; l++)
              TORD[l][0] = matrix_X_ordered[l][j+my_rank];

            kernel_weighted_sum_np(kernel_c,
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
                                   0,NULL,NULL,NULL,
                                   num_categories,
                                   NULL,
                                   NULL,
                                   kwm+(j+my_rank)*nrcc22,  // weighted sum
                                   NULL, // no permutations
                                   NULL); // do not return kernel weights

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
            }

            for(l = 0; l < num_reg_unordered; l++)
              TUNO[l][0] = matrix_X_unordered[l][j+my_rank];

            for(l = 0; l < num_reg_ordered; l++)
              TORD[l][0] = matrix_X_ordered[l][j+my_rank];

            kernel_weighted_sum_np(kernel_c,
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
                                   0,NULL,NULL,NULL,
                                   num_categories,
                                   NULL,
                                   NULL,
                                   kwm+(j+my_rank)*nrcc22,  // weighted sum
                                   NULL, // no permutations
                                   NULL); // do not return kernel weights

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
      if(ks_tree_use || (BANDWIDTH_reg == BW_ADAP_NN)){

        for(l = 0; l < num_reg_continuous; l++){
          
          for(i = 0; i < num_obs; i++){
            XTKX[l+2][i] = matrix_X_continuous[l][i]-matrix_X_continuous[l][j];
          }
          TCON[l][0] = matrix_X_continuous[l][j]; // temporary storage
        }


        for(l = 0; l < num_reg_unordered; l++)
          TUNO[l][0] = matrix_X_unordered[l][j];

        for(l = 0; l < num_reg_ordered; l++)
          TORD[l][0] = matrix_X_ordered[l][j];

        kernel_weighted_sum_np(kernel_c,
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
                               0,NULL,NULL,NULL,
                               num_categories,
                               NULL,
                               NULL,
                               kwm+j*nrcc22,  // weighted sum
                               NULL, // no permutations
                               NULL); // do not return kernel weights

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
          }

      
          for(l = 0; l < num_reg_unordered; l++)
            TUNO[l][0] = matrix_X_unordered[l][j];

          for(l = 0; l < num_reg_ordered; l++)
            TORD[l][0] = matrix_X_ordered[l][j];
      
          kernel_weighted_sum_np(kernel_c,
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
                                 0,NULL,NULL,NULL,
                                 num_categories,
                                 NULL,
                                 NULL,
                                 kwm+j*nrcc22, // weighted sum
                                 NULL, // no permutations
                                 NULL);  // no kernel weights

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
      
      // need to manipulate KWM pointers and XTKY - done

      if(bwm == RBWM_CVAIC){
        KWM[0][0] += aicc;
        XTKY[0][0] += aicc*vector_Y[j];
      }

      while(mat_inv(KWM, XTKXINV) == NULL){ // singular = ridge about
        for(int ii = 0; ii < (nrc1); ii++)
          KWM[ii][ii] += epsilon;
        nepsilon += epsilon;
      }
      
      if(bwm == RBWM_CVAIC)
        traceH += XTKXINV[0][0]*aicc;
   
      XTKY[0][0] += nepsilon*XTKY[0][0]/NZD(KWM[0][0]);

      DELTA = mat_mul(XTKXINV, XTKY, DELTA);
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
  }

  free(operator);
  free(kernel_c);
  free(kernel_u);
  free(kernel_o);
  free(lambda);
  free_mat(matrix_bandwidth,num_reg_continuous);

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
  int indy;

  int64_t i,j,l,iwx;

  int * operator = NULL;

  double **matrix_wX_unordered_eval=NULL;
  double **matrix_wX_ordered_eval=NULL;
  double **matrix_wX_continuous_eval=NULL;

  int64_t is,ie;

  int64_t N, num_obs_eval_alloc, num_obs_train_alloc, num_obs_wx_alloc;
  int64_t wx, nwx;

  size_t Nm = MIN((size_t)ceil(memfac*300000.0), (size_t)SIZE_MAX/10);

#ifdef MPI2
  int64_t stride_t = MAX((int64_t)ceil((double) num_obs_train / (double) iNum_Processors),1);

  is = stride_t * my_rank;
  ie = MIN(num_obs_train - 1, is + stride_t - 1);
  int64_t stride_e = MAX((int64_t)ceil((double) num_obs_eval / (double) iNum_Processors),1);
  
  num_obs_train_alloc = stride_t*iNum_Processors;
  num_obs_eval_alloc = stride_e*iNum_Processors;
#else
  is = 0;
  ie = num_obs_train - 1;

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
  matrix_wX_continuous_eval = (double **)malloc(num_reg_continuous*sizeof(double *));
  matrix_wX_unordered_eval = (double **)malloc(num_reg_unordered*sizeof(double *));
  matrix_wX_ordered_eval = (double **)malloc(num_reg_ordered*sizeof(double *));
 
  double * mean = (double *)malloc(num_obs_wx_alloc*sizeof(double));
  
  if(mean == NULL)
    error("failed to allocate mean");

  double ofac = num_obs_train - 1.0;

  operator = (int *)malloc(sizeof(int)*(num_reg_continuous+num_reg_unordered+num_reg_ordered));

  if(operator == NULL)
    error("failed to allocate operator");

  for(i = 0; i < (num_reg_continuous+num_reg_unordered+num_reg_ordered); i++)
    operator[i] = OP_INTEGRAL;

  int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;

  kernel_c = (int *)malloc(sizeof(int)*num_reg_continuous);

  for(i = 0; i < num_reg_continuous; i++)
    kernel_c[i] = KERNEL_den;

  kernel_u = (int *)malloc(sizeof(int)*num_reg_unordered);

  for(i = 0; i < num_reg_unordered; i++)
    kernel_u[i] = KERNEL_den_unordered;

  kernel_o = (int *)malloc(sizeof(int)*num_reg_ordered);

  for(i = 0; i < num_reg_ordered; i++)
    kernel_o[i] = KERNEL_den_ordered;

  
  *cv = 0;

  double * kwx = (double *)malloc(num_obs_train_alloc*num_obs_wx_alloc*sizeof(double));

  if(kwx == NULL)
    error("failed to allocate kwx, try reducing num_obs_eval");

  for(iwx = 0; iwx < nwx; iwx++){
    const int64_t wxo = iwx*wx;
    const int64_t dwx = (iwx != (nwx - 1)) ? wx : num_obs_eval - (nwx - 1)*wx;

    for(l = 0; l < num_reg_continuous; l++)
      matrix_wX_continuous_eval[l] = matrix_X_continuous_eval[l] + wxo;

    for(l = 0; l < num_reg_unordered; l++)
      matrix_wX_unordered_eval[l] = matrix_X_unordered_eval[l] + wxo;

    for(l = 0; l < num_reg_ordered; l++)
      matrix_wX_ordered_eval[l] = matrix_X_ordered_eval[l] + wxo;


    kernel_weighted_sum_np(kernel_c,
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
                           0,NULL,NULL,NULL,
                           num_categories,
                           matrix_categorical_vals,
                           NULL,
                           mean,
                           NULL, // no permutations
                           kwx);
    
    for(i = is; i <= ie; i++){
      for(j = wxo; j < (wxo + dwx); j++){             
        const int64_t jo = j - wxo;
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

  free(matrix_wX_continuous_eval);
  free(matrix_wX_unordered_eval);
  free(matrix_wX_ordered_eval);

  return(0);
}

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

  double vsfx[num_reg_tot];
  double vsfy[num_var_tot];
  double vsfxy[num_var_tot+num_reg_tot];
  double lambdax[num_reg_unordered+num_reg_ordered];
  double lambday[num_var_unordered+num_var_ordered];
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

  
  *cv = 0;

  if(!int_TREE_XY){
    double * kwx = (double *)malloc(num_obs_wx_alloc*num_obs_train_alloc*sizeof(double));

    if(kwx == NULL)
      error("failed to allocate kwx, tried to allocate: %" PRIi64 "bytes\n", num_obs_wx_alloc*num_obs_train_alloc*sizeof(double));

    double * kwy = (double *)malloc(num_obs_train_alloc*num_obs_wy_alloc*sizeof(double));

    if(kwy == NULL)
      error("failed to allocate kwy, try reducing num_obs_eval, tried to allocate: %" PRIi64 "bytes\n", num_obs_train_alloc*num_obs_wy_alloc*sizeof(double));
    
    for(iwx = 0; iwx < nwx; iwx++){
      const int64_t wxo = iwx*wx;
      const int64_t dwx = (iwx != (nwx - 1)) ? wx : num_obs_train - (nwx - 1)*wx;

      for(l = 0; l < num_reg_continuous; l++)
        matrix_wX_continuous_train[l] = matrix_X_continuous_train[l] + wxo;

      for(l = 0; l < num_reg_unordered; l++)
        matrix_wX_unordered_train[l] = matrix_X_unordered_train[l] + wxo;

      for(l = 0; l < num_reg_ordered; l++)
        matrix_wX_ordered_train[l] = matrix_X_ordered_train[l] + wxo;


      kernel_weighted_sum_np(kernel_cx,
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
                             0,NULL,NULL,NULL,
                             num_categories_extern_X,
                             matrix_categorical_vals_extern_X,
                             NULL, // moo
                             mean,
                             NULL, // no permutations
                             kwx);

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
        kernel_weighted_sum_np(kernel_cy,
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
                               0,NULL,NULL,NULL,
                               num_categories_extern_Y,
                               matrix_categorical_vals_extern_Y,
                               NULL,
                               NULL,
                               NULL, // no permutations
                               kwy);

        const int64_t je_dwy = MIN(je,dwy);

        for(i = wxo; i < (wxo + dwx); i++){     
          const int64_t io = i - wxo;

          for(j = (wyo + js); j < (wyo + je_dwy); j++){
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

#ifdef MPI2
    MPI_Allreduce(MPI_IN_PLACE, cv, 1, MPI_DOUBLE, MPI_SUM, comm[1]);
#endif
    *cv /= (double) num_obs_train*num_obs_eval;

    free(kwx);
    free(kwy);
  } else {
    NL nls = {.node = NULL, .n = 0, .nalloc = 0};
    NL nlps = {.node = NULL, .n = 0, .nalloc = 0};

    XL xl = {.istart = NULL, .nlev = NULL, .n = 0, .nalloc = 0};

    int icx[num_reg_continuous], icy[num_var_continuous];

    double bb[2*num_all_cvar];

    int KERNEL_XY[num_all_cvar], m;

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


      kernel_weighted_sum_np(kernel_cx,
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
                             0,NULL,NULL,NULL,
                             num_categories_extern_X,
                             matrix_categorical_vals_extern_X,
                             NULL,
                             mean,
                             NULL, // no permutations
                             kwx);

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
        kernel_weighted_sum_np(kernel_cy,
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
                               0,NULL,NULL,NULL,
                               num_categories_extern_Y,
                               matrix_categorical_vals_extern_Y,
                               NULL,
                               NULL,
                               NULL, // no permutations
                               kwy);

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

  return(0);

}

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
  
  int xyd[num_all_cvar];
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
    matrix_bandwidth_x = alloc_matd(num_obs_train, num_reg_continuous);
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



    matrix_bandwidth_y = alloc_matd(num_obs_train, num_var_continuous);
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

    matrix_bandwidth_xy = alloc_matd(num_obs_train, num_all_cvar);
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

  // extra kernel bookkeeping for trees
  int KERNEL_XY[num_all_cvar];
  double bb[2*num_all_cvar];

  for(l = 0; l < num_reg_continuous; l++)
    KERNEL_XY[l] = KERNEL_reg + OP_CFUN_OFFSETS[x_operator[l]];

  for(l = num_reg_continuous; l < num_all_cvar; l++)
    KERNEL_XY[l] = KERNEL_var + OP_CFUN_OFFSETS[y_operator[l-num_reg_continuous]];
  
  *cv = 0;

  // joint density
  kernel_weighted_sum_np(kernel_cxy,
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
                         NULL);

  // X density
  kernel_weighted_sum_np(kernel_cx,
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
                         NULL);

  if((!int_TREE_XY) && (!int_TREE_X)){
    for(i = is_i2n; i < ie_i2n; i++)
      *cv -= 2.0*jmean[i]/(mean[i] + DBL_MIN);
  } else {
    for(i = is_i2n; i < ie_i2n; i++)
      *cv -= 2.0*jmean[ipt_lookup_extern_XY[ipt_extern_X[i]]]/(mean[i] + DBL_MIN);
  }

  free(jmean);

  double * kx_ij = (double *)malloc(num_obs_wi_alloc*num_obs_wj_alloc*sizeof(double));
  pkx_ij = kx_ij;

  if(kx_ij == NULL)
    error("failed to allocate kx_ij, tried to allocate: %" PRIi64 "bytes\n", num_obs_wi_alloc*num_obs_wj_alloc*sizeof(double));

  double * kx_ik = (double *)malloc(num_obs_wi_alloc*num_obs_wk_alloc*sizeof(double));
  pkx_ik = kx_ik;

  if(kx_ik == NULL)
    error("failed to allocate kx_ik, tried to allocate: %" PRIi64 "bytes\n", num_obs_wi_alloc*num_obs_wk_alloc*sizeof(double));

  double * ky_jk = (double *)malloc(num_obs_wj_alloc*num_obs_wk_alloc*sizeof(double));
  pky_jk = ky_jk;

  if(ky_jk == NULL)
    error("failed to allocate ky_jk, tried to allocate: %" PRIi64 "bytes\n", num_obs_wj_alloc*num_obs_wk_alloc*sizeof(double));
    
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

      kernel_weighted_sum_np(kernel_cx,
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
                             kx_ij);

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
          kernel_weighted_sum_np(kernel_cx,
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
                                 kx_ik);

        } else {
          kx_ik = kx_ij;
        }
        // compute block ky_jk

        kernel_weighted_sum_np(kernel_cy,
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
                               ky_jk);

        if(!int_TREE_XY || (BANDWIDTH_den == BW_ADAP_NN)){
          const int64_t ie_dwi = MIN(ie,dwi);
          if(BANDWIDTH_den != BW_ADAP_NN){
            for(i = is+wio; i < (wio+ie_dwi); i++){
              for(j = wjo, tcvj = 0.0; j < (wjo + dwj); j++){              
                tcvk = 0.0;
                if(j == i){
                  continue;
                }

                const double tkxij = kx_ij[(i-wio)*dwj + j-wjo];

                if (tkxij != 0.0) {
                  for(k = wko; k < (wko + dwk); k++){
                    if(k == i) continue;
                    tcvk += kx_ik[(i-wio)*dwk + k-wko]*ky_jk[(j-wjo)*dwk + k-wko];
                  }               
                  tcvj += tkxij*tcvk;
                }
              }
              *cv += tcvj/(mean[i]*mean[i] + DBL_MIN);
            }
          } else {
            for(i = is+wio; i < (wio+ie_dwi); i++){
              for(j = wjo, tcvj = 0.0; j < (wjo + dwj); j++){              
                tcvk = 0.0;
                if(j == i){
                  continue;
                }

                const double tkxij = kx_ij[(j-wjo)*dwi + i-wio];

                if (tkxij != 0.0) {
                  for(k = wko; k < (wko + dwk); k++){
                    if(k == i) continue;
                    tcvk += kx_ik[(k-wko)*dwi + i-wio]*ky_jk[(j-wjo)*dwk + k-wko];
                  }               
                  tcvj += tkxij*tcvk;
                }
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

      free_mat(matrix_bandwidth_x, num_reg_continuous);
      free_mat(matrix_bandwidth_xy, num_all_cvar);
      free_mat(matrix_bandwidth_y, num_var_continuous);
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

  operator = (int *)malloc(sizeof(int)*(num_reg_continuous+num_reg_unordered+num_reg_ordered));

  for(i = 0; i < (num_reg_continuous+num_reg_unordered+num_reg_ordered); i++)
    operator[i] = OP_NORMAL;

#ifdef MPI2
  int stride_t = MAX((int)ceil((double) num_obs_train / (double) iNum_Processors),1);

  int stride_e = MAX((int)ceil((double) num_obs_eval / (double) iNum_Processors),1);
  int num_obs_eval_alloc = stride_e*iNum_Processors;

#else
  int num_obs_train_alloc = num_obs_train;
  int num_obs_eval_alloc = num_obs_eval;
#endif

  const int do_grad = (gradient != NULL); 
  const int do_gerr = (gradient_stderr != NULL);

  struct th_table * otabs = NULL;
  struct th_entry * ret = NULL;
  int ** matrix_ordered_indices = NULL;

  const int bwmdim = (BANDWIDTH_reg==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_train:1);

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
  matrix_bandwidth = alloc_matd(bwmdim,num_reg_continuous);

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
#ifdef MPI2
		MPI_Barrier(comm[1]);
		MPI_Finalize();
#endif
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
      te.key.dkey = ret->key.dkey;
      te.data = ret->data;

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

  // Conduct the estimation 

  if(int_ll == LL_LC) { // local constant
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
      assert(permy != NULL);
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
    
    kernel_weighted_sum_np(kernel_c,
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
                           0,NULL,NULL,NULL,
                           num_categories,
                           matrix_categorical_vals,
                           matrix_ordered_indices, 
                           meany,
                           permy, // permutations used for gradients
                           NULL); // do not return kernel weights

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
    MATRIX XTKXINV = mat_creat( num_reg_continuous + 1, num_reg_continuous + 1, UNDEFINED );
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

    double * kwm = (double *)malloc(nrcc33*num_obs_eval_alloc*sizeof(double));

    for(int ii = 0; ii < nrcc33*num_obs_eval_alloc; ii++)
      kwm[ii] = 0.0;

    // with local linear, we already have the gradients of the continuous components
    // so we only need to worry about unordered + ordered comps

    double * permy = NULL;
    int ** moo = NULL;

    int p_nvar = do_grad ? (num_reg_unordered + num_reg_ordered) : 0;
    int do_ocg = do_grad && (p_nvar > 0);

    if(do_ocg){
      permy = (double *)malloc(nrcc33*num_obs_eval_alloc*p_nvar*sizeof(double));
      assert(permy != NULL);
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
          }


          for(l = 0; l < num_reg_unordered; l++)
            TUNO[l][0] = matrix_X_unordered_eval[l][j+my_rank];

          for(l = 0; l < num_reg_ordered; l++)
            TORD[l][0] = matrix_X_ordered_eval[l][j+my_rank];

          kernel_weighted_sum_np(kernel_c,
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
                                 0,NULL,NULL,NULL,
                                 num_categories,
                                 matrix_categorical_vals,
                                 moo,
                                 kwm+(j+my_rank)*nrcc33,  // weighted sum
                                 do_ocg ? (permy+(j+my_rank)*nrcc33*p_nvar) : NULL, // ocg
                                 NULL); // do not return kernel weights

        }
        // synchro step
        MPI_Allgather(MPI_IN_PLACE, nrcc33, MPI_DOUBLE, kwm+j*nrcc33, nrcc33, MPI_DOUBLE, comm[1]);          
        MPI_Allgather(MPI_IN_PLACE, nrcc33*p_nvar, MPI_DOUBLE, permy+j*nrcc33*p_nvar, nrcc33*p_nvar, MPI_DOUBLE, comm[1]);    
      }
      
#else


      for(l = 0; l < num_reg_continuous; l++){
          
        for(i = 0; i < num_obs_train; i++){
          XTKX[l+3][i] = matrix_X_continuous_train[l][i]-matrix_X_continuous_eval[l][j];
        }
        TCON[l][0] = matrix_X_continuous_eval[l][j]; // temporary storage
      }


      for(l = 0; l < num_reg_unordered; l++)
        TUNO[l][0] = matrix_X_unordered_eval[l][j];

      for(l = 0; l < num_reg_ordered; l++)
        TORD[l][0] = matrix_X_ordered_eval[l][j];

      kernel_weighted_sum_np(kernel_c,
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
                             0,NULL,NULL,NULL,
                             num_categories,
                             matrix_categorical_vals,
                             moo,
                             kwm+j*nrcc33,  // weighted sum
                             do_ocg ? (permy+j*nrcc33*p_nvar) : NULL, // no permutations
                             NULL); // do not return kernel weights

#endif
      
      if(do_ocg){
        for(l = 0; l < num_reg_ordered; l++){
          moo[l]++;
        }
      }

      while(mat_inv(KWM, XTKXINV) == NULL){ // singular = ridge about
        for(int ii = 0; ii < (nrc1); ii++)
          KWM[ii][ii] += epsilon;
        nepsilon += epsilon;
      }
      
      XTKY[0][0] += nepsilon*XTKY[0][0]/NZD(KWM[0][0]);

      DELTA = mat_mul(XTKXINV, XTKY, DELTA);
      mean[j] = DELTA[0][0];

      const double sk = copysign(DBL_MIN, (kwm+j*nrcc33)[2*nrc3+2]) + (kwm+j*nrcc33)[2*nrc3+2];
      const double ey = (kwm+j*nrcc33)[nrc3+2]/sk;
      const double ey2 = (kwm+j*nrcc33)[nrc3+1]/sk;

      mean_stderr[j] = sqrt((ey2 - ey*ey)*K_INT_KERNEL_P / (sk*hprod));

      if(do_grad){
        for(int ii = 0; ii < num_reg_continuous; ii++){
          gradient[ii][j] = DELTA[ii+1][0];
          if(do_gerr)
            gradient_stderr[ii][j] = gfac*mean_stderr[j]/((BANDWIDTH_reg == BW_ADAP_NN) ? 1.0 : ((BANDWIDTH_reg == BW_GEN_NN) ? matrix_bandwidth[ii][j]:matrix_bandwidth[ii][0]));
        }
        
        // we need to do new matrix inversions here for the unordered + ordered data
        // we can safely taint KWM , XTKXINV, DELTA, and XTKY here
        for(l = num_reg_continuous; l < (num_reg_continuous + num_reg_unordered); l++){
          const int dl = l - num_reg_continuous;
          const int ojp = j*nrcc33*p_nvar + dl*nrcc33;

          for(int ii = 0; ii < nrc1; ii++){
            KWM[ii] = &permy[ojp +(ii+2)*(nrc3)+2];
            XTKY[ii] = &permy[ojp + ii + nrc3 + 2];
          }

          nepsilon = 0.0;
          while(mat_inv(KWM, XTKXINV) == NULL){ // singular = ridge about
            for(int ii = 0; ii < (nrc1); ii++)
              KWM[ii][ii] += epsilon;
            nepsilon += epsilon;
          }

          XTKY[0][0] += nepsilon*XTKY[0][0]/NZD(KWM[0][0]);

          DELTA = mat_mul(XTKXINV, XTKY, DELTA);

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
        // we can safely taint KWM , XTKXINV, DELTA, and XTKY here
        for(l = num_reg_continuous + num_reg_unordered; l < (num_reg_continuous + num_reg_unordered + num_reg_ordered); l++){

          const int dl = l - num_reg_continuous;
          const int ojp = j*nrcc33*p_nvar + dl*nrcc33;

          for(int ii = 0; ii < nrc1; ii++){
            KWM[ii] = &permy[ojp +(ii+2)*(nrc3)+2];
            XTKY[ii] = &permy[ojp + ii + nrc3 + 2];
          }

          nepsilon = 0.0;
          while(mat_inv(KWM, XTKXINV) == NULL){ // singular = ridge about
            for(int ii = 0; ii < (nrc1); ii++)
              KWM[ii][ii] += epsilon;
            nepsilon += epsilon;
          }

          XTKY[0][0] += nepsilon*XTKY[0][0]/NZD(KWM[0][0]);

          DELTA = mat_mul(XTKXINV, XTKY, DELTA);

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
    mat_free(XTKXINV);
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

  }

  // clean up hash stuff
  if(do_grad && (num_reg_ordered > 0)){
    for(l = 0; l < num_reg_ordered; l++)
      thdestroy_r(otabs+l);
    free(otabs);
    free(matrix_ordered_indices[0]);
    free(matrix_ordered_indices);
  }

  free(operator);
  free(kernel_c);
  free(kernel_u);
  free(kernel_o);
  free(lambda);
  free_mat(matrix_bandwidth,num_reg_continuous);
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

  const int num_reg = num_reg_continuous+num_reg_unordered+num_reg_ordered;

  int i;

  int * operator = NULL;

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

  kernel_weighted_sum_np(kernel_c,
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
                         0,NULL,NULL,NULL,
                         num_categories,
                         NULL,
                         NULL,
                         rho,  // weighted sum
                         NULL, // no permutations
                         NULL); // do not return kernel weights
  
  

  for(i = 0, *cv = 0.0; i < num_obs; i++)
    *cv -= (rho[i] < DBL_MIN) ? log_DBL_MIN : log(rho[i]/(num_obs-1.0));
    

  free(operator);
  free(kernel_c);
  free(kernel_u);
  free(kernel_o);
  free(rho);
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

  const int num_reg = num_reg_continuous+num_reg_unordered+num_reg_ordered;

  int i;

  int * operator = NULL;

  int num_obs_alloc;

#ifdef MPI2
  int stride_t = MAX((int)ceil((double) num_obs / (double) iNum_Processors),1);
  
  num_obs_alloc = stride_t*iNum_Processors;
#else
  num_obs_alloc = num_obs;
#endif

  double * res = (double *)malloc(num_obs_alloc*sizeof(double));
  double cv1, cv2;
  
  if(res == NULL)
    error("failed to allocate rho");

  operator = (int *)malloc(sizeof(int)*num_reg);

  // first the convolution portion
  for(i = 0; i < num_reg; i++)
    operator[i] = OP_CONVOLUTION;

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


  kernel_weighted_sum_np(kernel_c,
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
                         0,NULL,NULL,NULL,
                         num_categories,
                         matrix_categorical_vals,
                         NULL, // no ocg
                         res,  // weighted sum
                         NULL, // no permutations
                         NULL); // do not return kernel weights
  
  

  for(i = 0, cv1 = 0.0; i < num_obs; i++) cv1 += res[i];

  cv1 /= num_obs*num_obs;

  // then the cross term
  for(i = 0; i < num_reg; i++)
    operator[i] = OP_NORMAL;

  kernel_weighted_sum_np(kernel_c,
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
                         0,NULL,NULL,NULL,
                         num_categories,
                         NULL,
                         NULL, // no ocg
                         res,  // weighted sum
                         NULL, // no permutations
                         NULL); // do not return kernel weights


  for(i = 0, cv2 = 0.0; i < num_obs; i++) cv2 += res[i];

  cv2 /= num_obs*(num_obs-1.0);

  *cv = cv1 - 2.0*cv2;

  free(operator);
  free(kernel_c);
  free(kernel_u);
  free(kernel_o);
  free(res);
  return(0);

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

  const int num_reg = num_reg_continuous+num_reg_unordered+num_reg_ordered;

  int i, l;

  const int bwmdim = (BANDWIDTH_den==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_den==BW_ADAP_NN)?num_obs_train:1);

  int * operator = NULL;

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
#ifdef MPI2
		MPI_Barrier(comm[1]);
		MPI_Finalize();
#endif
    error("\n** Error: invalid bandwidth.");
  }

  kernel_weighted_sum_np(kernel_c,
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
                         0, //  explicity suppress parallel
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
                         0,NULL,NULL,NULL,
                         num_categories,
                         matrix_categorical_vals, // if dist mcv (possibly) necessary
                         NULL, // no ocg
                         pdf,  // weighted sum
                         NULL, // no permutations
                         NULL); // do not return kernel weights

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


  const int num_reg = num_reg_continuous+num_reg_unordered+num_reg_ordered;
  const int num_cvar = num_reg_continuous + num_var_continuous;
  const int num_uvar = num_reg_unordered + num_var_unordered;
  const int num_ovar = num_reg_ordered + num_var_ordered;
  const int num_all_var = num_reg + num_var_continuous + num_var_unordered + num_var_ordered;

	const double log_DBL_MIN = log(DBL_MIN);

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

  // xy
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
                         0,NULL,NULL,NULL,
                         num_categories_extern_XY,
                         matrix_categorical_vals_extern_XY,
                         NULL,
                         rhon,  // weighted sum
                         NULL, // no permutations
                         NULL); // do not return kernel weights

  //x
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
                         0,NULL,NULL,NULL,
                         num_categories_extern_X,
                         matrix_categorical_vals_extern_X,
                         NULL,
                         rhod,  // weighted sum
                         NULL, // no permutations
                         NULL); // do not return kernel weights
  
  for(i = 0, *cv = 0.0; i < num_obs; i++){
    if((rhon[i] == 0.0) || (rhod[i] == 0.0)){
      ret = 1;
      break;
      //*cv -= log_DBL_MIN;
    } else {
    *cv -= log(rhon[i]) - log(rhod[i]);
    }
  }

  free(operator);
  free(kernel_cx);
  free(kernel_ux);
  free(kernel_ox);

  free(kernel_cxy);
  free(kernel_uxy);
  free(kernel_oxy);

  free(rhon);
  free(rhod);

  free(vsf_xy);
  free(vsf_x);
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
    otabs = (struct th_table *)malloc(num_X_ordered*sizeof(struct th_table));
    matrix_ordered_indices = (int **)malloc(num_oXY*sizeof(int *));
    int * tc = (int *)malloc(num_X_ordered*num_obs_eval*sizeof(int));
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
      te.key.dkey = ret->key.dkey;
      te.data = ret->data;

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


  operator_XY = (int *)malloc(sizeof(int)*num_XY);
  operator_X = (int *)malloc(sizeof(int)*num_X);

  kernel_cXY = (int *)malloc(sizeof(int)*num_cXY);
  kernel_uXY = (int *)malloc(sizeof(int)*num_uXY);
  kernel_oXY = (int *)malloc(sizeof(int)*num_oXY);

  vsf_XY = (double *)malloc(num_XY*sizeof(double));
  vsf_X = (double *)malloc(num_X*sizeof(double));

  ksd = (double *)malloc(num_obs_eval_alloc*sizeof(double));
  ksn = (double *)malloc(num_obs_eval_alloc*sizeof(double));

  icX = (int *)malloc(sizeof(int)*num_X_continuous);

  if(do_grad){
    permn = (double *)malloc(num_X*num_obs_eval_alloc*sizeof(double));
    permd = (double *)malloc(num_X*num_obs_eval_alloc*sizeof(double));
    bpso = (int *)malloc(num_XY*sizeof(int));

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
#ifdef MPI2
    MPI_Barrier(comm[1]);
    MPI_Finalize();
#endif
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
                         0, // do not explicity suppress parallel
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
                         NULL); // do not return kernel weights

  //x - we assume x is in xy tree order
  //  - we also reuse kernels and operators because we can

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
                         0, //  do not explicity suppress parallel
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
                         NULL); // do not return kernel weights

  
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
