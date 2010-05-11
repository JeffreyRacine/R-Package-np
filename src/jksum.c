/* Copyright (C) J. Racine, 1995-2001 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <assert.h>

#include <R_ext/Utils.h>

#include "headers.h"
#include "matrix.h"



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

extern int int_DEBUG;
extern int int_VERBOSE;
extern int int_TAYLOR;
extern int int_WEIGHTS;
extern int int_LARGE_SF;

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

	double *psum;
	double *py;

#ifdef MPI2
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
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
		vector_scale_factor,
		matrix_X_continuous_train,	 /* Not used */
		matrix_X_continuous_eval,		 /* Not used */
		matrix_X_continuous_train,
		matrix_X_continuous_eval,
		matrix_bandwidth,						 /* Not used */
		matrix_bandwidth,
		lambda) == 1)
	{
#ifndef MPI2
		printf("\n** Error: invalid bandwidth.");
		printf("\nProgram Terminated.\n");
		exit(EXIT_FAILURE);
#endif
#ifdef MPI2
		if(my_rank == 0)
		{
			printf("\n** Error: invalid bandwidth.");
			printf("\nProgram Terminated.\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		exit(EXIT_FAILURE);
#endif
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

	MPI_Gather(kernel_sum, stride, MPI_DOUBLE, kernel_sum, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(kernel_sum, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

	free(lambda);

	free_mat(matrix_bandwidth,num_reg_continuous);

	return(0);

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

double np_uli_racine(const int same_cat, const double lambda, const int c){
  return (same_cat)?1.0:lambda;
}

double np_owang_van_ryzin(const double x, const double y, const double lambda){
  return (x == y)?(1.0-lambda):ipow(lambda, (int)fabs(x-y))*(1.0-lambda)*0.5;
}

double np_oli_racine(const double x, const double y, const double lambda){
  return (x == y)?1.0:ipow(lambda, (int)fabs(x-y));

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

// derivative kernels

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

// end kernels

double (* const allck[])(double) = { np_gauss2, np_gauss4, np_gauss6, np_gauss8, 
                                  np_epan2, np_epan4, np_epan6, np_epan8, 
                                  np_rect };
double (* const allok[])(double, double, double) = { np_owang_van_ryzin, np_oli_racine };
double (* const alluk[])(int, double, int) = { np_uaa, np_uli_racine };


/* 
   np_kernelv does weighted products of vectors - this is useful for 
   product kernels, where each kernel in each dimension acts as a weight.
*/

/* xt = training data */
/* xw = x weights */

void np_ckernelv(const int KERNEL, 
                 const double * const xt, const int num_xt, 
                 const int do_xw,
                 const double x, const double h, 
                 double * const result){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i,j; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double (* const k[])(double) = { np_gauss2, np_gauss4, np_gauss6, np_gauss8, //ordinary kernels
                                   np_epan2, np_epan4, np_epan6, np_epan8, 
                                   np_rect, 
                                   np_econvol_gauss2, np_econvol_gauss4, np_econvol_gauss6, np_econvol_gauss8, // convolution kernels
                                   np_econvol_epan2, np_econvol_epan4, np_econvol_epan6, np_econvol_epan8,
                                   np_deriv_gauss2, np_deriv_gauss4, np_deriv_gauss6, np_deriv_gauss8, // derivative kernels
                                   np_deriv_epan2, np_deriv_epan4, np_deriv_epan6, np_deriv_epan8, 
                                   np_deriv_rect,
                                   np_cdf_gauss2, np_cdf_gauss4, np_cdf_gauss6, np_cdf_gauss8, // cdfative kernels
                                   np_cdf_epan2, np_cdf_epan4, np_cdf_epan6, np_cdf_epan8, 
                                   np_cdf_rect };

  for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw)
    result[i] = xw[j]*k[KERNEL]((x-xt[i])/h);

}

void np_convol_ckernelv(const int KERNEL, 
                        const double * const xt, const int num_xt, 
                        const int do_xw,
                        const double x, 
                        double * xt_h, 
                        const double h, 
                        double * const result){

  int i,j; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw)
    result[i] = xw[j]*kernel_convol(KERNEL, BW_ADAP_NN, 
                                    (x-xt[i])/xt_h[i], xt_h[i], h);

}

void np_ukernelv(const int KERNEL, 
                 const double * const xt, const int num_xt, 
                 const int do_xw,
                 const double x, const double lambda, const int ncat,
                 double * const result){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i; 
  int j, bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double (* const k[])(int, double, int) = { np_uaa, np_uli_racine,
                                             np_econvol_uaa, np_econvol_uli_racine };

  for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw)
    result[i] = xw[j]*k[KERNEL]((xt[i]==x), lambda, ncat);

}

void np_convol_okernelv(const int KERNEL, 
                        const double * const xt, const int num_xt, 
                        const int do_xw,
                        const double x, const double lambda,
                        int ncat, double * cat,
                        double * const result){

  int i; 
  int j, bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw)
    result[i] = xw[j]*kernel_ordered_convolution(KERNEL, xt[i], x, lambda, ncat, cat);
}

void np_okernelv(const int KERNEL, 
                 const double * const xt, const int num_xt, 
                 const int do_xw,
                 const double x, const double lambda,
                 double * const result){
  
  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i; 
  int j, bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double (* const k[])(double, double, double) = { np_owang_van_ryzin, np_oli_racine };

  for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw)
    result[i] = xw[j]*k[KERNEL](xt[i], x, lambda);
}

// W = A
// Y = B
// outer product = AB' (assuming column vector convention)
void np_outer_weighted_sum(double * const * const mat_A, double * const sgn_A, const int ncol_A, 
                           double * const * const mat_B, const int ncol_B,
                           double * const weights, const int num_weights,
                           const int do_leave_one_out, const int which_k,
                           const int kpow,
                           const int parallel_sum,
                           const int symmetric,
                           const int gather_scatter,
                           const int bandwidth_divide, const double dband,
                           double * const result){

  int i,j,k, l = parallel_sum?which_k:0;
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

  if(!gather_scatter){
    if(!symmetric){
      if(kpow == 1){
        for (k = 0; k < num_weights; k++, l+=linc)
          for (j = 0; j < max_A; j++)
            for (i = 0; i < max_B; i++)
              result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*weights[k]/db;
      } else { // kpow != 1
        for (k = 0; k < num_weights; k++, l+=linc)
          for (j = 0; j < max_A; j++)
            for (i = 0; i < max_B; i++)
              result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*ipow(weights[k]/db, kpow);
      }
    } else { // symmetric
      if(kpow == 1){
        for (k = 0; k < num_weights; k++, l+=linc){
          for (j = 0; j < max_A; j++){
            for (i = 0; i <= j; i++){
              result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*weights[k]/db;
            }
          }
        }
      } else { // kpow != 1
        for (k = 0; k < num_weights; k++, l+=linc){
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
        for (k = 0; k < num_weights; k++, l+=linc)
          for (j = 0; j < max_A; j++)
            for (i = 0; i < max_B; i++)
              result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*weights[k]/db;
      } else { // kpow != 1
        for (k = 0; k < num_weights; k++, l+=linc)
          for (j = 0; j < max_A; j++)
            for (i = 0; i < max_B; i++)
              result[k*kstride+j*max_B+i] += pmat_A[j][l*have_A]*pmat_B[i][l*have_B]*ipow(weights[k]/db, kpow);
      }
    } else { // symmetric
      if(kpow == 1){
        for (k = 0; k < num_weights; k++){
          for (j = 0; j < max_A; j++){
            for (i = 0; i <= j; i++){
              const double tp = pmat_A[j][k*have_A]*pmat_B[i][k*have_B]*weights[k]/db;
              const double sgnp = tp*p_sgn[j*have_sgn]*p_sgn[i*have_sgn];
              result[k*kstride+j*max_B+i] += tp;
              result[(k+1)*max_AB+j*max_B+i] += sgnp;
              //fprintf(stderr,"\nmax ab %d\tsgnp %e\tres %e",max_AB,sgnp,result[(k+1)*max_AB+j*max_B+i]);
            }
          }
          // insert reference weight in top right corner
          result[(k+1)*max_AB+max_B-1] = weights[k]/db;
        }
      } else { // kpow != 1
        for (k = 0; k < num_weights; k++){
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

  if (do_leave_one_out)
    weights[which_k] = temp;
}



int kernel_weighted_sum_np(
const int KERNEL_reg,
const int KERNEL_unordered_reg,
const int KERNEL_ordered_reg,
const int BANDWIDTH_reg,
const int num_obs_train,
const int num_obs_eval,
const int num_reg_unordered,
const int num_reg_ordered,
const int num_reg_continuous,
const int leave_one_out,
const int kernel_pow,
const int bandwidth_divide,
const int do_smooth_coef_weights,
const int symmetric,
const int gather_scatter,
const int operator,
double * const * const matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double **matrix_Y,
double **matrix_W,
double * sgn,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *weighted_sum){
  
  /* This function takes a vector Y and returns a kernel weighted
     leave-one-out sum. By default Y should be a vector of ones
     (simply compute the kernel sum). This function will allow users
     to `roll their own' with mixed data leave-one out kernel sums. */

  /* num_var_ordered_extern contains number of columns of the weight matrix */


  /* Declarations */

  int i,j,l, mstep, js, je, num_obs_eval_alloc, sum_element_length;
  int do_psum; 

  const int KERNEL_reg_np = KERNEL_reg + OP_CFUN_OFFSETS[operator];
  const int KERNEL_unordered_reg_np = KERNEL_unordered_reg + OP_UFUN_OFFSETS[operator];
  const int num_xt = (BANDWIDTH_reg == BW_ADAP_NN)?num_obs_eval:num_obs_train;
  const int ws_step = (BANDWIDTH_reg == BW_ADAP_NN)? 0 :
    (MAX(num_var_continuous_extern, 1) * MAX(num_var_ordered_extern, 1));

  double *lambda, **matrix_bandwidth, **matrix_eval_bandwidth = NULL, *m = NULL;
  double *tprod, dband, *buf, *ws;

  double * const * const xtc = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_continuous_eval:matrix_X_continuous_train;
  double * const * const xtu = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_unordered_eval:matrix_X_unordered_train;
  double * const * const xto = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_ordered_eval:matrix_X_ordered_train;

  double * const * const xc = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_continuous_train:matrix_X_continuous_eval;
  double * const * const xu = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_unordered_train:matrix_X_unordered_eval;
  double * const * const xo = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_ordered_train:matrix_X_ordered_eval;

#ifdef MPI2
  int stride = MAX(ceil((double) num_obs_eval / (double) iNum_Processors),1);
  num_obs_eval_alloc = stride*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval;
#endif

  if (num_obs_eval == 0) {
    return(1);
  }

  do_psum = BANDWIDTH_reg == BW_ADAP_NN;
  /* Allocate memory for objects */

  mstep = (BANDWIDTH_reg==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_train:1);

  lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);
  matrix_bandwidth = alloc_tmatd(mstep, num_reg_continuous);  

  tprod = alloc_vecd((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_eval:num_obs_train);

  sum_element_length = MAX(num_var_continuous_extern, 1) * 
    MAX(num_var_ordered_extern, 1);

  buf = alloc_vecd(num_obs_eval_alloc*sum_element_length);

  /* clear memory in temporary buffer and in result */
  if(gather_scatter){
    for(i = 0; i < num_obs_eval_alloc*sum_element_length; i++){
      buf[i] = weighted_sum[i];
    }
  }else{
    for(i = 0; i < num_obs_eval_alloc*sum_element_length; i++){
      buf[i] = 0.0;
      weighted_sum[i] = 0.0;
    }    
  }

  /* assert(!(BANDWIDTH_reg == BW_ADAP_NN)); */
  /* Conduct the estimation */

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

    return(1);
  }

  if((BANDWIDTH_reg == BW_ADAP_NN) && (operator == OP_CONVOLUTION)){ // need additional bandwidths 
    matrix_eval_bandwidth = alloc_tmatd(num_obs_eval, num_reg_continuous);  

    if(kernel_bandwidth_mean(
                             KERNEL_reg,
                             BW_GEN_NN, // this is not an error!
                             num_obs_train,
                             num_obs_eval,
                             0,
                             0,
                             0,
                             num_reg_continuous,
                             num_reg_unordered,
                             num_reg_ordered,
                             vector_scale_factor,
                             NULL,				 /* Not used */
                             NULL,				 /* Not used */
                             matrix_X_continuous_train,
                             matrix_X_continuous_eval,
                             NULL,						 /* Not used */
                             matrix_eval_bandwidth,
                             lambda)==1){

      free(lambda);
      free_tmat(matrix_bandwidth);
      free_tmat(matrix_eval_bandwidth);
      free(tprod);

      return(1);
    }

  }
  
  if ((num_obs_train != num_obs_eval) && leave_one_out){
    
    printf("\ntraining and evaluation data must be the same to use leave one out estimator");
    free(lambda);
    free_tmat(matrix_bandwidth);
    free(tprod);
    return(1);
  }

  if (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN){
#ifdef MPI2
    js = stride * my_rank;
    je = MIN(num_obs_eval - 1, js + stride - 1);
#else
    js = 0;
    je = num_obs_eval - 1;
#endif
    
    ws = weighted_sum + js * sum_element_length;
  } else {
#ifdef MPI2
    js = stride * my_rank;
    je = MIN(num_obs_train - 1, js + stride - 1);
    ws = buf;
#else
    js = 0;
    je = num_obs_train - 1;
    ws = weighted_sum;
#endif
  }
  
    /* do sums */
  for(j=js; j <= je; j++, ws += ws_step){
    R_CheckUserInterrupt();

    dband = 1.0;

    if (num_reg_continuous > 0){
      m = matrix_bandwidth[0];
      if (BANDWIDTH_reg != BW_FIXED)
        m += j;
    }

    l = 0;

    /* continuous first */

    /* for the first iteration, no weights */
    /* for the rest, the accumulated products are the weights */
    if((BANDWIDTH_reg != BW_ADAP_NN) || (operator != OP_CONVOLUTION)){
      for(i=0; i < num_reg_continuous; i++, l++, m += mstep){
        np_ckernelv(KERNEL_reg_np, xtc[i], num_xt, l, xc[i][j], *m, tprod);
        dband *= *m;
      }
    } else {
      for(i=0; i < num_reg_continuous; i++, l++, m += mstep){
        np_convol_ckernelv(KERNEL_reg, xtc[i], num_xt, l, xc[i][j], 
                           matrix_eval_bandwidth[i], *m, tprod);
        dband *= *m;
      }      
    }

    /* unordered second */

    for(i=0; i < num_reg_unordered; i++, l++){
      np_ukernelv(KERNEL_unordered_reg_np, xtu[i], num_xt, l, xu[i][j], 
                  lambda[i], num_categories[i], tprod);
    }

    /* ordered third */
    if(operator != OP_CONVOLUTION){
      for(i=0; i < num_reg_ordered; i++, l++){
        np_okernelv(KERNEL_ordered_reg, xto[i], num_xt, l,
                    xo[i][j], lambda[num_reg_unordered+i], 
                    tprod);
      }
    } else {
      for(i=0; i < num_reg_ordered; i++, l++){
        np_convol_okernelv(KERNEL_ordered_reg, xto[i], num_xt, l,
                           xo[i][j], lambda[num_reg_unordered+i], 
                           num_categories[i+num_reg_unordered],
                           matrix_categorical_vals[i+num_reg_unordered],
                           tprod);
      }
    }

    /* expand matrix outer product, multiply by kernel weights, etc, do sum */

    np_outer_weighted_sum(matrix_W, sgn, num_var_ordered_extern, 
                          matrix_Y, num_var_continuous_extern,
                          tprod, num_obs_train,
                          leave_one_out, j,
                          kernel_pow,
                          do_psum,
                          symmetric,
                          gather_scatter,
                          bandwidth_divide, dband,
                          ws);

  }

#ifdef MPI2
  if (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN){
    MPI_Allgather(weighted_sum + js * sum_element_length, stride * sum_element_length, MPI_DOUBLE, weighted_sum, stride * sum_element_length, MPI_DOUBLE, MPI_COMM_WORLD);
  } else if(BANDWIDTH_reg == BW_ADAP_NN){
    MPI_Allreduce(ws, weighted_sum, num_obs_eval*sum_element_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
#endif

  free(lambda);
  free_tmat(matrix_bandwidth);

  if((BANDWIDTH_reg == BW_ADAP_NN) && (operator == OP_CONVOLUTION))
    free_tmat(matrix_eval_bandwidth);

  free(tprod);
  free(buf);
  
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
  double (* const yok)(double, double, double) = allok[KERNEL_ordered_den];
  double (* const xok)(double, double, double) = allok[KERNEL_ordered_reg];
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

  const uint64_t Wi = P*Q + (P > Nr)?Nr:P;
  const uint64_t Wf = Wi + Q + (P < Nr);
  const uint64_t uki = Wi/(Nb*Nb);
  const uint64_t uji = (Wi%(Nb*Nb))/Nb;
  const uint64_t uii = Wi%Nb;

  uint64_t W = Wi;

  k = (int)uki;
  j = (int)uji;
  i = (int)uii;

  sum_kerf = (double *)malloc(num_obs*sizeof(double));
  assert(sum_kerf != NULL);

  sum_ker_convolf = (double *)malloc(num_obs*sizeof(double));
  assert(sum_ker_convolf != NULL);

  sum_ker_marginalf = (double *)malloc(num_obs*sizeof(double));
  assert(sum_ker_marginalf != NULL);
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
  assert(blk_xi != NULL);

  blk_xj = (double *)malloc(sizeof(double)*blklen*blklen);
  assert(blk_xj != NULL);

  blk_yij = (double *)malloc(sizeof(double)*blklen*blklen);
  assert(blk_xj != NULL);

  sum_ker = (double *)malloc(num_obs*sizeof(double));
  assert(sum_ker != NULL);

  sum_ker_convol = (double *)malloc(num_obs*sizeof(double));
  assert(sum_ker_convol != NULL);

  sum_ker_marginal = (double *)malloc(num_obs*sizeof(double));
  assert(sum_ker_marginal != NULL);

  assert(BANDWIDTH_den == BW_FIXED);

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
            blk_xj[ib] *= xok(matrix_X_ordered[l][blk+n],matrix_X_ordered[l][blj+m],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
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
            ts *= yok(matrix_Y_ordered[l][blk+n],matrix_Y_ordered[l][blj+m],lambda[l+num_var_unordered]);
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
              blk_xi[ib] *= xok(matrix_X_ordered[l][blk+n],matrix_X_ordered[l][bli+m],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
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
  MPI_Allreduce(sum_ker, sum_kerf, num_obs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(sum_ker_convol, sum_ker_convolf, num_obs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(sum_ker_marginal, sum_ker_marginalf, num_obs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(j = 0; j < num_obs; j++){
    if(sum_ker_marginalf[j] <= 0.0){
      *cv = DBL_MAX;
      break;
    }
    *cv += (sum_ker_convolf[j]/sum_ker_marginalf[j]-2.0*sum_kerf[j])/sum_ker_marginalf[j];
  }
#else
  for(j = 0; j < num_obs; j++){
    if(sum_ker_marginal[j] <= 0.0){
      *cv = DBL_MAX;
      break;
    }
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
  int i, j, l, sf_flag = 0;

  double cv = 0.0;
  double * lambda = NULL, * vsf = NULL;
  double ** matrix_bandwidth = NULL;

  double aicc = 0.0;
  double traceH = 0.0;

  const int leave_one_out = (bwm == RBWM_CVLS)?1:0;

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

    num_var_continuous_extern = 0; // rows in the y_mat
    num_var_ordered_extern = 0; // rows in weights

    kernel_weighted_sum_np(KERNEL_reg,
                           KERNEL_unordered_reg,
                           KERNEL_ordered_reg,
                           BANDWIDTH_reg,
                           1,
                           1,
                           num_reg_unordered,
                           num_reg_ordered,
                           num_reg_continuous,
                           0, // do not leave out 
                           1, // kernel_pow = 1
                           (BANDWIDTH_reg == BW_ADAP_NN)?1:0, // bandwidth_divide = FALSE when not adaptive
                           0, // do_smooth_coef_weights = FALSE (not implemented)
                           0, // not symmetric
                           0, // do not gather-scatter
                           0, // no convolution
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
                           num_categories,
                           NULL,
                           &aicc);
    //fprintf(stderr,"\n%e\n",aicc);
  }

  // Conduct the estimation 

  if(int_ll == LL_LC) { // local constant
    // Nadaraya-Watson
    // Generate bandwidth vector given scale factors, nearest neighbors, or lambda 

    double * lc_Y[2];
    double * mean = (double *)malloc(2*num_obs*sizeof(double));

    num_var_continuous_extern = 2; // columns in the y_mat
    num_var_ordered_extern = 0; // columns in weights
    lc_Y[0] = vector_Y;
      
    lc_Y[1] = (double *)malloc(num_obs*sizeof(double));
    for(int ii = 0; ii < num_obs; ii++)
      lc_Y[1][ii] = 1.0;

    kernel_weighted_sum_np(KERNEL_reg,
                           KERNEL_unordered_reg,
                           KERNEL_ordered_reg,
                           BANDWIDTH_reg,
                           num_obs,
                           num_obs,
                           num_reg_unordered,
                           num_reg_ordered,
                           num_reg_continuous,
                           leave_one_out, 
                           1, // kernel_pow = 1
                           (BANDWIDTH_reg == BW_ADAP_NN)?1:0, // bandwidth_divide = FALSE when not adaptive
                           0, // do_smooth_coef_weights = FALSE (not implemented)
                           0, // not symmetric
                           0, // do not gather-scatter
                           0, // no convolution
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
                           num_categories,
                           NULL,
                           mean);
    
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


    

    //exit(1);

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

    double * PXC[num_reg_continuous]; 
    double * PXU[num_reg_unordered];
    double * PXO[num_reg_ordered];

    for(l = 0; l < num_reg_continuous; l++)
      PXC[l] = matrix_X_continuous[l];

    for(l = 0; l < num_reg_unordered; l++)
      PXU[l] = matrix_X_unordered[l];

    for(l = 0; l < num_reg_ordered; l++)
      PXO[l] = matrix_X_ordered[l];

    double * kwm = (double *)malloc(nrcc22*num_obs*sizeof(double));

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


    for(j = 0; j < num_obs; j++){
      nepsilon = 0.0;

      for(l = 0; l < (nrc1); l++){
        KWM[l] = &kwm[j*nrcc22+(l+1)*(nrc2)+1];
        XTKY[l] = &kwm[j*nrcc22+l+1];
      }

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
      

        num_var_continuous_extern = nrc2; // rows in the y_mat
        num_var_ordered_extern = nrc2; // rows in weights

        kernel_weighted_sum_np(KERNEL_reg,
                               KERNEL_unordered_reg,
                               KERNEL_ordered_reg,
                               BANDWIDTH_reg,
                               num_obs-j-1,
                               1,
                               num_reg_unordered,
                               num_reg_ordered,
                               num_reg_continuous,
                               0, // we leave one out via the weight matrix
                               1, // kernel_pow = 1
                               (BANDWIDTH_reg == BW_ADAP_NN)?1:0, // bandwidth_divide = FALSE when not adaptive
                               0, // do_smooth_coef_weights = FALSE (not implemented)
                               1, // symmetric
                               1, // gather-scatter sum
                               0, // no convolution
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
                               num_categories,
                               NULL,
                               kwm+j*nrcc22); // weighted sum
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
   
      XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,KWM[0][0]));

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

  free(lambda);
  free_mat(matrix_bandwidth,num_reg_continuous);

  cv /= (double)num_obs;
  if(bwm == RBWM_CVAIC){
    cv = log(cv) + (1.0+traceH/((double)num_obs))/(1.0-(traceH+2.0)/((double)num_obs));
  }
  return(cv);
}
