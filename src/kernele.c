/* Copyright (C) J. Racine, 1995-2001 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <errno.h>

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

#define IO_MIN_TRUE  1
#define IO_MIN_FALSE 0

extern int int_DEBUG;
extern int int_VERBOSE;
extern int int_MINIMIZE_IO;
extern int int_TAYLOR;
extern int int_WEIGHTS;

#ifdef RCSID
static char rcsid[] = "$Id: kernele.c,v 1.12 2006/11/02 19:50:13 tristen Exp $";
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

int kernel_estimate_density_categorical(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int BANDWIDTH_den,
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
double *vector_scale_factor,
int *num_categories,
double *pdf,
double *pdf_stderr,
double *log_likelihood)
{

	/* This function estimates a density function using both continuous */
	/* and categorical covariates with three estimation techniques and an */
	/* assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel;

	double sum_ker;
	double sum_ker_temp;

	double prod_nh;

	double *lambda;
	double **matrix_bandwidth = NULL;
	double **matrix_bandwidth_deriv = NULL;

	double INT_KERNEL_P;					 /* Integral of K(z)^p */
	double K_INT_KERNEL_P;				 /* Number of regressors times integral of K(z)^p */

	double log_DBL_MIN = log(DBL_MIN);

	#ifdef MPI2
	double log_likelihood_MPI;
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

	if((BANDWIDTH_den == 0)||(BANDWIDTH_den == 1))
	{
		matrix_bandwidth = alloc_matd(num_obs_eval,num_reg_continuous);
		matrix_bandwidth_deriv = alloc_matd(num_obs_eval,num_reg_continuous);
	}
	else if(BANDWIDTH_den == 2)
	{
		matrix_bandwidth = alloc_matd(num_obs_train,num_reg_continuous);
		matrix_bandwidth_deriv = alloc_matd(num_obs_train,num_reg_continuous);
	}

	/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

	if(kernel_bandwidth(
		KERNEL_den,
		BANDWIDTH_den,
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
		lambda,
		matrix_bandwidth_deriv) == 1)
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

	/* Initialize constants for various kernels required for asymptotic standard errors */

	initialize_kernel_density_asymptotic_constants(
		KERNEL_den,
		num_reg_continuous,
		&INT_KERNEL_P,
		&K_INT_KERNEL_P);

	#ifndef MPI2

	/* Initialize log likelihood */

	*log_likelihood = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			prod_nh = (double) num_obs_train;

			for(l = 0; l < num_reg_continuous; l++)
			{
				prod_nh *= matrix_bandwidth[l][0];
			}

			sum_ker = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel;

			}

			pdf[j] = sum_ker/prod_nh;

			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j] = ((num_reg_continuous != 0) ? sqrt(pdf[j]*K_INT_KERNEL_P/prod_nh) : sqrt(pdf[j]/prod_nh));

			if(pdf[j] > DBL_MIN)
			{
				*log_likelihood += log(pdf[j]);
			}
			else
			{
				*log_likelihood += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_density_categorical()");
				}
			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			prod_nh = (double) num_obs_train;

			for(l = 0; l < num_reg_continuous; l++)
			{
				prod_nh *= matrix_bandwidth[l][j];
			}

			sum_ker = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel;

			}

			pdf[j] = sum_ker/prod_nh;

			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j] = ((num_reg_continuous != 0) ? sqrt(pdf[j]*K_INT_KERNEL_P/prod_nh) : sqrt(pdf[j]/prod_nh));

			if(pdf[j] > DBL_MIN)
			{
				*log_likelihood += log(pdf[j]);
			}
			else
			{
				*log_likelihood += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_density_categorical()");
				}
			}

		}

	}
	else
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = sum_ker_temp = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_nh = (double) num_obs_train;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_nh *= matrix_bandwidth[l][i];
				}

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel/prod_nh;
				sum_ker_temp += prod_kernel/ipow(prod_nh, 2);

			}

			pdf[j] = sum_ker;

			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j] = ((num_reg_continuous != 0) ? sqrt(sum_ker_temp*K_INT_KERNEL_P) : sqrt(sum_ker_temp));

			if(pdf[j] > DBL_MIN)
			{
				*log_likelihood += log(pdf[j]);
			}
			else
			{
				*log_likelihood += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_density_categorical()");
				}
			}

		}

	}
	#endif

	#ifdef MPI2

	/* Initialize log likelihood */

	log_likelihood_MPI = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			prod_nh = (double) num_obs_train;

			for(l = 0; l < num_reg_continuous; l++)
			{
				prod_nh *= matrix_bandwidth[l][0];
			}

			sum_ker = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel;

			}

			pdf[j-my_rank*stride] = sum_ker/prod_nh;

			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(pdf[j-my_rank*stride]*K_INT_KERNEL_P/prod_nh) : sqrt(pdf[j-my_rank*stride]/prod_nh));

			if(pdf[j-my_rank*stride] > DBL_MIN)
			{
				log_likelihood_MPI += log(pdf[j-my_rank*stride]);
			}
			else
			{
				log_likelihood_MPI += log_DBL_MIN;
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_density_categorical()");
				}
			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			prod_nh = (double) num_obs_train;

			for(l = 0; l < num_reg_continuous; l++)
			{
				prod_nh *= matrix_bandwidth[l][j];
			}

			sum_ker = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel;

			}

			pdf[j-my_rank*stride] = sum_ker/prod_nh;

			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(pdf[j-my_rank*stride]*K_INT_KERNEL_P/prod_nh) : sqrt(pdf[j-my_rank*stride]/prod_nh));

			if(pdf[j-my_rank*stride] > DBL_MIN)
			{
				log_likelihood_MPI += log(pdf[j-my_rank*stride]);
			}
			else
			{
				log_likelihood_MPI += log_DBL_MIN;
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_density_categorical()");
				}
			}

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = sum_ker_temp = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_nh = (double) num_obs_train;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_nh *= matrix_bandwidth[l][i];
				}

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel/prod_nh;
				sum_ker_temp += prod_kernel/ipow(prod_nh,2);

			}

			pdf[j-my_rank*stride] = sum_ker;

			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(sum_ker_temp*K_INT_KERNEL_P) : sqrt(sum_ker_temp));

			if(pdf[j-my_rank*stride] > DBL_MIN)
			{
				log_likelihood_MPI += log(pdf[j-my_rank*stride]);
			}
			else
			{
				log_likelihood_MPI += log_DBL_MIN;
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_density_categorical()");
				}
			}

		}

	}

	MPI_Gather(pdf, stride, MPI_DOUBLE, pdf, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(pdf, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Gather(pdf_stderr, stride, MPI_DOUBLE, pdf_stderr, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(pdf_stderr, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Reduce(&log_likelihood_MPI, log_likelihood, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(log_likelihood, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	free(lambda);
	free_mat(matrix_bandwidth,num_reg_continuous);
	free_mat(matrix_bandwidth_deriv,num_reg_continuous);

	return(0);

}


int kernel_estimate_density_categorical_leave_one_out_cv(
int KERNEL_den,
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
double *cv)
{

	/* This function estimates a leave one out density function using both */
	/* continuous and categorical covariates with three estimation techniques */
	/* and an assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel;

	double sum_ker;

	double prod_nh;
	double temp_bw1;
	double temp_bw2;

	double pdf;

	double *lambda;
	double **matrix_bandwidth;

	double log_DBL_MIN = log(DBL_MIN);

	double *p_xj1;
	double *p_xi1;
	double *p_xj2;
	double *p_xi2;

	#ifdef MPI2
	double cv_MPI;
	int stride = ceil((double) num_obs / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

	#ifndef MPI2
	matrix_bandwidth = alloc_matd(num_obs,num_reg_continuous);
	#endif

	#ifdef MPI2
	matrix_bandwidth = alloc_matd(stride*iNum_Processors,num_reg_continuous);
	#endif

	/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

	if(kernel_bandwidth_mean(
		KERNEL_den,
		BANDWIDTH_den,
		num_obs,
		num_obs,
		0,
		0,
		0,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
		vector_scale_factor,
		matrix_X_continuous,				 /* Not used */
		matrix_X_continuous,				 /* Not used */
		matrix_X_continuous,
		matrix_X_continuous,
		matrix_bandwidth,						 /* Not used */
		matrix_bandwidth,
		lambda)==1)
	{
		free(lambda);
		free_mat(matrix_bandwidth,num_reg_continuous);
		return(1);
	}

	#ifndef MPI2

	*cv = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		prod_nh = (double) num_obs;

		for(l = 0; l < num_reg_continuous; l++)
		{
			/* Fixed bandwidth */
			prod_nh *= matrix_bandwidth[l][0];
		}

		if((num_reg_continuous == 1)&&((num_reg_unordered+num_reg_ordered) == 0))
		{

			/* For special case of fixed bandwidth for PDF/CDF for test stat */

			temp_bw1 = matrix_bandwidth[0][0];
			p_xj1 = &matrix_X_continuous[0][0];

			for(j=0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();

				sum_ker = 0.0;

				p_xi1 = &matrix_X_continuous[0][0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{
						sum_ker += kernel(KERNEL_den, (*p_xj1 - *p_xi1)/temp_bw1);
					}

					p_xi1++;

				}

				p_xj1++;

				pdf = sum_ker/prod_nh;

				if(pdf > DBL_MIN)
				{
					*cv -= log(pdf);
				}
				else
				{
					*cv -= log_DBL_MIN;
					if(int_VERBOSE == 1)
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_density_categorical_leave_one_out_cv()");
					}
				}

			}

		}
		else if((num_reg_continuous == 2)&&((num_reg_unordered+num_reg_ordered) == 0))
		{

			/* For special case of fixed bandwidth for PDF/CDF for test stat */

			temp_bw1 = matrix_bandwidth[0][0];
			temp_bw2 = matrix_bandwidth[1][0];

			p_xj1 = &matrix_X_continuous[0][0];
			p_xj2 = &matrix_X_continuous[1][0];

			for(j=0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();

				sum_ker = 0.0;

				p_xi1 = &matrix_X_continuous[0][0];
				p_xi2 = &matrix_X_continuous[1][0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{
						sum_ker += kernel(KERNEL_den, (*p_xj1 - *p_xi1)/temp_bw1)*kernel(KERNEL_den, (*p_xj2 - *p_xi2)/temp_bw2);
					}

					p_xi1++;
					p_xi2++;

				}

				p_xj1++;
				p_xj2++;

				pdf = sum_ker/prod_nh;

				if(pdf > DBL_MIN)
				{
					*cv -= log(pdf);
				}
				else
				{
					*cv -= log_DBL_MIN;
					if(int_VERBOSE == 1)
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_density_categorical_leave_one_out_cv()");
					}
				}

			}

		}
		else
		{

			for(j=0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();

				sum_ker = 0.0;

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						prod_kernel = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][0]);
						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
						}

						sum_ker += prod_kernel;

					}

				}

				pdf = sum_ker/prod_nh;

				if(pdf > DBL_MIN)
				{
					*cv -= log(pdf);
				}
				else
				{
					*cv -= log_DBL_MIN;
					if(int_VERBOSE == 1)
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_density_categorical_leave_one_out_cv()");
					}
				}

			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs; j++)
		{
		  R_CheckUserInterrupt();

			prod_nh = (double) num_obs;

			for(l = 0; l < num_reg_continuous; l++)
			{
				prod_nh *= matrix_bandwidth[l][j];
			}

			sum_ker = 0.0;

			for(i=0; i < num_obs; i++)
			{

				if(i != j)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{

						prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][j]);

					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;

				}

			}

			pdf = sum_ker/prod_nh;

			if(pdf > DBL_MIN)
			{
				*cv -= log(pdf);
			}
			else
			{
				*cv -= log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_density_categorical_leave_one_out_cv()");
				}
			}

		}

	}
	else
	{

		for(j=0; j < num_obs; j++)
		{
		  R_CheckUserInterrupt();

			sum_ker = 0.0;

			for(i=0; i < num_obs; i++)
			{

				prod_nh = (double) num_obs;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_nh *= matrix_bandwidth[l][i];
				}

				if(i != j)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][i]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;

				}

			}

			pdf = sum_ker;

			if(pdf > DBL_MIN)
			{
				*cv -= log(pdf);
			}
			else
			{
				*cv -= log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_density_categorical_leave_one_out_cv()");
				}
			}

		}
	}

	*cv /= (double) num_obs;
	#endif

	#ifdef MPI2

	cv_MPI = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		prod_nh = (double) num_obs;

		for(l = 0; l < num_reg_continuous; l++)
		{
			/* Fixed bandwidth */
			prod_nh *= matrix_bandwidth[l][0];
		}

		if((num_reg_continuous == 1)&&((num_reg_unordered+num_reg_ordered) == 0))
		{

			/* For special case of fixed bandwidth for PDF/CDF for test stat */

			temp_bw1 = matrix_bandwidth[0][0];

			p_xj1 = &matrix_X_continuous[0][my_rank*stride];

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				sum_ker = 0.0;

				p_xi1 = &matrix_X_continuous[0][0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{
						sum_ker += kernel(KERNEL_den, (*p_xj1 - *p_xi1)/temp_bw1);
					}

					p_xi1++;

				}

				p_xj1++;

				pdf = sum_ker/prod_nh;

				if(pdf > DBL_MIN)
				{
					cv_MPI -= log(pdf);
				}
				else
				{
					cv_MPI -= log_DBL_MIN;
					if((int_VERBOSE == 1)&&(my_rank == 0))
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_density_categorical_leave_one_out_cv()");
					}
				}

			}

		}
		else if((num_reg_continuous == 2)&&((num_reg_unordered+num_reg_ordered) == 0))
		{

			/* For special case of fixed bandwidth for PDF/CDF for test stat */

			temp_bw1 = matrix_bandwidth[0][0];
			temp_bw2 = matrix_bandwidth[1][0];

			p_xj1 = &matrix_X_continuous[0][my_rank*stride];
			p_xj2 = &matrix_X_continuous[1][my_rank*stride];

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				sum_ker = 0.0;

				p_xi1 = &matrix_X_continuous[0][0];
				p_xi2 = &matrix_X_continuous[1][0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{
						sum_ker += kernel(KERNEL_den, (*p_xj1 - *p_xi1)/temp_bw1)*kernel(KERNEL_den, (*p_xj2 - *p_xi2)/temp_bw2);
					}

					p_xi1++;
					p_xi2++;

				}

				p_xj1++;
				p_xj2++;

				pdf = sum_ker/prod_nh;

				if(pdf > DBL_MIN)
				{
					cv_MPI -= log(pdf);
				}
				else
				{
					cv_MPI -= log_DBL_MIN;
					if((int_VERBOSE == 1)&&(my_rank == 0))
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_density_categorical_leave_one_out_cv()");
					}
				}

			}

		}
		else
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				sum_ker = 0.0;

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						prod_kernel = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][0]);
						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
						}

						sum_ker += prod_kernel;

					}

				}

				pdf = sum_ker/prod_nh;

				if(pdf > DBL_MIN)
				{
					cv_MPI -= log(pdf);
				}
				else
				{
					cv_MPI -= log_DBL_MIN;
					if((int_VERBOSE == 1)&&(my_rank == 0))
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_density_categorical_leave_one_out_cv()");
					}
				}

			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
		{

			prod_nh = (double) num_obs;

			for(l = 0; l < num_reg_continuous; l++)
			{
				prod_nh *= matrix_bandwidth[l][j];
			}

			sum_ker = 0.0;

			for(i=0; i < num_obs; i++)
			{

				if(i != j)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{

						prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][j]);

					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;

				}

			}

			pdf = sum_ker/prod_nh;

			if(pdf > DBL_MIN)
			{
				cv_MPI -= log(pdf);
			}
			else
			{
				cv_MPI -= log_DBL_MIN;
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_density_categorical_leave_one_out_cv()");
				}
			}

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;

			for(i=0; i < num_obs; i++)
			{

				prod_nh = (double) num_obs;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_nh *= matrix_bandwidth[l][i];
				}

				if(i != j)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][i]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;

				}

			}

			pdf = sum_ker;

			if(pdf > DBL_MIN)
			{
				cv_MPI -= log(pdf);
			}
			else
			{
				cv_MPI -= log_DBL_MIN;
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_density_categorical_leave_one_out_cv()");
				}
			}

		}
	}

	/* Now reduce */

	cv_MPI /= (double) num_obs;
	MPI_Reduce(&cv_MPI, cv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(cv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	free(lambda);
	free_mat(matrix_bandwidth,num_reg_continuous);

	return(0);

}


int kernel_estimate_con_density_categorical_leave_one_out_cv(
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
double *cv)
{

	/* This function estimates a leave one out density function using both*/
	/* continuous and categorical covariates with three estimation techniques*/
	/* and an assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel_cat;
	double prod_kernel_cont;

	double prod_kernel_marginal_cat;
	double prod_kernel_marginal_cont;

	double sum_ker;
	double sum_ker_marginal;

	double prod_h;
	double prod_h_marginal;

	double pdf;

	double *lambda;
	double **matrix_bandwidth_var;
	double **matrix_bandwidth_reg;

	double log_DBL_MIN = log(DBL_MIN);

	double temp;

	#ifdef MPI2
	double cv_MPI;
	int stride = ceil((double) num_obs / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered);
	matrix_bandwidth_var = alloc_matd(num_obs,num_var_continuous);
	matrix_bandwidth_reg = alloc_matd(num_obs,num_reg_continuous);

	/* Generate bandwidth vector for continuous dep vars */

	if(kernel_bandwidth_mean(
		KERNEL_den,
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
		matrix_bandwidth_reg,				 /* Not used */
		lambda)==1)
	{
		free(lambda);
		free_mat(matrix_bandwidth_var,num_var_continuous);
		free_mat(matrix_bandwidth_reg,num_reg_continuous);
		return(1);
	}

	#ifndef MPI2

	*cv = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		prod_h = 1.0;

		for(l = 0; l < num_var_continuous; l++)
		{
			/* Fixed bandwidth */
			prod_h *= matrix_bandwidth_var[l][0];
		}

		for(j=0; j < num_obs; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs; i++)
			{

				if(i != j)
				{

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][0]);
					}

					prod_kernel_marginal_cont = prod_kernel_cont;

					for(l = 0; l < num_var_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][0]);
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
					}

					prod_kernel_marginal_cat = prod_kernel_cat;

					for(l = 0; l < num_var_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_var_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
					}

					sum_ker += prod_kernel_cont*prod_kernel_cat;
					sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				}

			}

			temp = prod_h*sum_ker_marginal;

			pdf = sum_ker/temp;

			if(pdf > DBL_MIN)
			{
				*cv -= log(pdf);
			}
			else
			{
				*cv -= log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_leave_one_out_cv()");
				}
			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs; j++)
		{
		  R_CheckUserInterrupt();
			prod_h = 1.0;

			for(l = 0; l < num_var_continuous; l++)
			{
				prod_h *= matrix_bandwidth_var[l][j];
			}

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs; i++)
			{

				if(i != j)
				{

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{

						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][j]);

					}

					prod_kernel_marginal_cont = prod_kernel_cont;

					for(l = 0; l < num_var_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][j]);
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
					}

					prod_kernel_marginal_cat = prod_kernel_cat;

					for(l = 0; l < num_var_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_var_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
					}

					sum_ker += prod_kernel_cont*prod_kernel_cat;
					sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				}

			}

			pdf = sum_ker/(prod_h*sum_ker_marginal);

			if(pdf > DBL_MIN)
			{
				*cv -= log(pdf);
			}
			else
			{
				*cv -= log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_leave_one_out_cv()");
				}
			}

		}

	}
	else
	{

		for(j=0; j < num_obs; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs; i++)
			{

				prod_h = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					/* Fixed bandwidth */
					prod_h *= matrix_bandwidth_reg[l][i];
				}

				prod_h_marginal = prod_h;

				for(l = 0; l < num_var_continuous; l++)
				{
					/* Fixed bandwidth */
					prod_h *= matrix_bandwidth_var[l][i];
				}

				if(i != j)
				{

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][i]);
					}

					prod_kernel_marginal_cont = prod_kernel_cont;

					for(l = 0; l < num_var_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][i]);
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
					}

					prod_kernel_marginal_cat = prod_kernel_cat;

					for(l = 0; l < num_var_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_var_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
					}

					sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
					sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h_marginal;

				}

			}

			pdf = sum_ker/sum_ker_marginal;

			if(pdf > DBL_MIN)
			{
				*cv -= log(pdf);
			}
			else
			{
				*cv -= log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_leave_one_out_cv()");
				}
			}

		}
	}

	*cv /= (double) num_obs;
	#endif

	#ifdef MPI2

	cv_MPI = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		prod_h = 1.0;

		for(l = 0; l < num_var_continuous; l++)
		{
			/* Fixed bandwidth */
			prod_h *= matrix_bandwidth_var[l][0];
		}

		for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs; i++)
			{

				if(i != j)
				{

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][0]);
					}

					prod_kernel_marginal_cont = prod_kernel_cont;

					for(l = 0; l < num_var_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][0]);
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
					}

					prod_kernel_marginal_cat = prod_kernel_cat;

					for(l = 0; l < num_var_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_var_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
					}

					sum_ker += prod_kernel_cont*prod_kernel_cat;
					sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				}

			}

			temp = prod_h*sum_ker_marginal;

			pdf = sum_ker/temp;

			if(pdf > DBL_MIN)
			{
				cv_MPI -= log(pdf);
			}
			else
			{
				cv_MPI -= log_DBL_MIN;
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_leave_one_out_cv()");
				}
			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
		{

			prod_h = 1.0;

			for(l = 0; l < num_var_continuous; l++)
			{
				prod_h *= matrix_bandwidth_var[l][j];
			}

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs; i++)
			{

				if(i != j)
				{

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{

						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][j]);

					}

					prod_kernel_marginal_cont = prod_kernel_cont;

					for(l = 0; l < num_var_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][j]);
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
					}

					prod_kernel_marginal_cat = prod_kernel_cat;

					for(l = 0; l < num_var_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_var_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
					}

					sum_ker += prod_kernel_cont*prod_kernel_cat;
					sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				}

			}

			pdf = sum_ker/(prod_h*sum_ker_marginal);

			if(pdf > DBL_MIN)
			{
				cv_MPI -= log(pdf);
			}
			else
			{
				cv_MPI -= log_DBL_MIN;
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_leave_one_out_cv()");
				}
			}

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs; i++)
			{

				prod_h = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					/* Fixed bandwidth */
					prod_h *= matrix_bandwidth_reg[l][i];
				}

				prod_h_marginal = prod_h;

				for(l = 0; l < num_var_continuous; l++)
				{
					/* Fixed bandwidth */
					prod_h *= matrix_bandwidth_var[l][i];
				}

				if(i != j)
				{

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][i]);
					}

					prod_kernel_marginal_cont = prod_kernel_cont;

					for(l = 0; l < num_var_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][i]);
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
					}

					prod_kernel_marginal_cat = prod_kernel_cat;

					for(l = 0; l < num_var_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_var_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
					}

					sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
					sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h_marginal;

				}

			}

			pdf = sum_ker/sum_ker_marginal;

			if(pdf > DBL_MIN)
			{
				cv_MPI -= log(pdf);
			}
			else
			{
				cv_MPI -= log_DBL_MIN;
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_leave_one_out_cv()");
				}
			}

		}
	}

	/* Now reduce */

	cv_MPI /= (double) num_obs;
	MPI_Reduce(&cv_MPI, cv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(cv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	free(lambda);
	free_mat(matrix_bandwidth_var,num_var_continuous);
	free_mat(matrix_bandwidth_reg,num_reg_continuous);

	return(0);

}


int kernel_estimate_distribution_categorical(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int BANDWIDTH_den,
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
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *cdf,
double *cdf_stderr,
double small,
int itmax)
{

	/* This function estimates a distribution function using both continuous */
	/* and categorical covariates with three estimation techniques and an */
	/* assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel;

	double sum_ker;

	double *lambda;
	double **matrix_bandwidth = NULL;
	double **matrix_bandwidth_deriv = NULL;

	#ifdef MPI2
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

	if((BANDWIDTH_den == 0)||(BANDWIDTH_den == 1))
	{
		matrix_bandwidth = alloc_matd(num_obs_eval,num_reg_continuous);
		matrix_bandwidth_deriv = alloc_matd(num_obs_eval,num_reg_continuous);
	}
	else if(BANDWIDTH_den == 2)
	{
		matrix_bandwidth = alloc_matd(num_obs_train,num_reg_continuous);
		matrix_bandwidth_deriv = alloc_matd(num_obs_train,num_reg_continuous);
	}

	/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

	if(kernel_bandwidth(
		KERNEL_den,
		BANDWIDTH_den,
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
		lambda,
		matrix_bandwidth_deriv) == 1)
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

	/* Conduct the estimation */

	#ifndef MPI2

	if(BANDWIDTH_den == 0)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= cdf_kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered],matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel;

			}

			cdf[j] = sum_ker/(double)num_obs_train;
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= cdf_kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered],matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel;

			}

			cdf[j] = sum_ker/(double)num_obs_train;
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}
	else
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= cdf_kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered],matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel;

			}

			cdf[j] = sum_ker/(double)num_obs_train;
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}
	#endif

	#ifdef MPI2

	if(BANDWIDTH_den == 0)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= cdf_kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered],matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel;

			}

			cdf[j-my_rank*stride] = sum_ker/(double)num_obs_train;
			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= cdf_kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered],matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel;

			}

			cdf[j-my_rank*stride] = sum_ker/(double)num_obs_train;
			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= cdf_kernel(KERNEL_den, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered],matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_ker += prod_kernel;

			}

			cdf[j-my_rank*stride] = sum_ker/(double)num_obs_train;
			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

		}

	}

	MPI_Gather(cdf, stride, MPI_DOUBLE, cdf, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(cdf, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(cdf_stderr, stride, MPI_DOUBLE, cdf_stderr, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(cdf_stderr, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	free(lambda);

	free_mat(matrix_bandwidth,num_reg_continuous);
	free_mat(matrix_bandwidth_deriv,num_reg_continuous);

	return(0);

}


int kernel_estimate_regression_categorical(
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
double **matrix_X_continuous_bandwidth,
double *vector_Y,
double *vector_Y_eval,
double *vector_scale_factor,
int *num_categories,
double *mean,
double **gradient,
double *mean_stderr,
double **gradient_stderr,
double *R_squared,
double *MSE,
double *MAE,
double *MAPE,
double *CORR,
double *SIGN)
{

	/* This function estimates a Nadaraya-Watson or Local Linear regression */
	/* function using both continuous and categorical covariates with three */
	/* estimation techniques and an assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int k;
	int l = INT_MAX;

	const double epsilon = 1.0/num_obs_train;
  double nepsilon;

	double prod_kernel;

	double prod_kernel_cat;
	double prod_kernel_cont;

	double sum_ker;
	double sum_y_ker;
	double sum_y_sq_ker;

	double *prod_kernel_deriv;
	double *sum_ker_deriv;
	double *sum_y_ker_deriv;

	double *lambda;
	double **matrix_bandwidth = NULL;
	double **matrix_bandwidth_deriv = NULL;

	double temp_var;
	double temp_mean_y;
	double *pointer_yi;
	double *pointer_m;
	double temp;
	double temp1 = DBL_MAX;

	MATRIX  XTKX;
	MATRIX  XTKXINV;
	MATRIX  XTKY;
	MATRIX  XTKYSQ;
	MATRIX  DELTA;

	double INT_KERNEL_P;					 /* Integral of K(z)^p */
	double K_INT_KERNEL_P;				 /* Number of regressors times integral of K(z)^p */
	/* Integral of K(z-0.5)*K(z+0.5) */
	double INT_KERNEL_PM_HALF = 0.0;
	double DIFF_KER_PPM = 0.0;		 /* Difference between int K(z)^p and int K(z-.5)K(z+.5) */

	int num_reg_cat_cont;

	#ifdef MPI2
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	if(int_TAYLOR == 1)
	{
		num_reg_cat_cont = num_reg_unordered + num_reg_ordered + num_reg_continuous;
	}
	else
	{
		num_reg_cat_cont = num_reg_continuous;
	}

	/* Initialize constants for various kernels required for asymptotic standard errors */

	initialize_kernel_regression_asymptotic_constants(
		KERNEL_reg,
		num_reg_continuous,
		&INT_KERNEL_P,
		&K_INT_KERNEL_P,
		&INT_KERNEL_PM_HALF,
		&DIFF_KER_PPM);

	/* Allocate memory for objects */

	prod_kernel_deriv = alloc_vecd(num_reg_continuous);
	sum_ker_deriv = alloc_vecd(num_reg_continuous);
	sum_y_ker_deriv = alloc_vecd(num_reg_continuous);

	lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

	if((BANDWIDTH_reg == 0)||(BANDWIDTH_reg == 1))
	{
		matrix_bandwidth = alloc_matd(num_obs_eval,num_reg_continuous);
		matrix_bandwidth_deriv = alloc_matd(num_obs_eval,num_reg_continuous);
	}
	else if(BANDWIDTH_reg == 2)
	{
		matrix_bandwidth = alloc_matd(num_obs_train,num_reg_continuous);
		matrix_bandwidth_deriv = alloc_matd(num_obs_train,num_reg_continuous);
	}

	/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

	if(kernel_bandwidth(
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
																 /* Not used */
		matrix_X_continuous_bandwidth,
		matrix_X_continuous_eval,		 /* Not used */
		matrix_X_continuous_bandwidth,
		matrix_X_continuous_eval,
		matrix_bandwidth,						 /* Not used */
		matrix_bandwidth,
		lambda,
		matrix_bandwidth_deriv) == 1)
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

	if(int_ll == 0)
	{

		/* Conduct Nadaraya-Watson estimation */

		if(BANDWIDTH_reg == 0)
		{

			pointer_m = &mean[0];

			for(j=0; j < num_obs_eval; j++)
			{
			  R_CheckUserInterrupt();
				sum_y_ker = sum_y_sq_ker = 0.0;
				sum_ker = DBL_MIN;

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_y_ker_deriv[l] = sum_ker_deriv[l] = 0.0;
				}

				for(i=0; i < num_obs_train; i++)
				{

					/* Kernel for mean */

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg,matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Kernels for derivatives */

					for (k = 0; k < num_reg_continuous; k++)
					{

						prod_kernel_deriv[k] = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][0]);

						for (l = 0; l < num_reg_continuous; l++)
						{

							if(l != k)
							{
								prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
							}

						}

					}

					/* Mean */

					sum_ker += prod_kernel;
					sum_y_ker += vector_Y[i]*prod_kernel;
					sum_y_sq_ker += ipow(vector_Y[i],2)*prod_kernel;

					/* Gradient */

					for(l = 0; l < num_reg_continuous; l++)
					{
						sum_ker_deriv[l] += prod_kernel_deriv[l];
						sum_y_ker_deriv[l] += vector_Y[i]*prod_kernel_deriv[l];
					}

				}

				*pointer_m = sum_y_ker/sum_ker;

				temp_var = (sum_y_sq_ker/sum_ker) - ipow(*pointer_m++, 2);

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / sum_ker) : sqrt(temp_var / sum_ker));

				/* gradient[0][] is that for _first_ continuous variable */

				for(l = 0; l < num_reg_continuous; l++)
				{
					gradient[l][j] = (sum_y_ker_deriv[l]-(sum_y_ker/sum_ker)*sum_ker_deriv[l])/(matrix_bandwidth_deriv[l][0]*sum_ker);

					gradient_stderr[l][j] = sqrt(temp_var * DIFF_KER_PPM /
						(sum_ker * ipow(matrix_bandwidth_deriv[l][0],2)));
				}

			}
		}
		else if(BANDWIDTH_reg == 1)
		{

			pointer_m = &mean[0];

			for(j=0; j < num_obs_eval; j++)
			{
			  R_CheckUserInterrupt();
				sum_y_ker = sum_y_sq_ker = 0.0;
				sum_ker = DBL_MIN;

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_y_ker_deriv[l] = sum_ker_deriv[l] = 0.0;
				}

				for(i=0; i < num_obs_train; i++)
				{

					/* Kernel for mean */

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg,matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Kernels for derivatives */

					for (k = 0; k < num_reg_continuous; k++)
					{

						prod_kernel_deriv[k] = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][j]);

						for (l = 0; l < num_reg_continuous; l++)
						{

							if(l != k)
							{
								prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
							}

						}

					}

					/* Mean */

					sum_ker += prod_kernel;
					sum_y_ker += vector_Y[i]*prod_kernel;
					sum_y_sq_ker += ipow(vector_Y[i],2)*prod_kernel;

					/* Gradient */

					for(l = 0; l < num_reg_continuous; l++)
					{
						sum_ker_deriv[l] += prod_kernel_deriv[l];
						sum_y_ker_deriv[l] += vector_Y[i]*prod_kernel_deriv[l];
					}

				}

				*pointer_m = sum_y_ker/sum_ker;

				temp_var = (sum_y_sq_ker/sum_ker) - ipow(*pointer_m++, 2);

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / sum_ker) : sqrt(temp_var / sum_ker));

				/* gradient[0][] is that for _first_ continuous variable */

				for(l = 0; l < num_reg_continuous; l++)
				{
					gradient[l][j] = (sum_y_ker_deriv[l]-(sum_y_ker/sum_ker)*sum_ker_deriv[l])/(matrix_bandwidth_deriv[l][j]*sum_ker);

					gradient_stderr[l][j] = sqrt(temp_var * DIFF_KER_PPM /
						(sum_ker * ipow(matrix_bandwidth_deriv[l][j],2)));
				}

			}

		}
		else
		{

			pointer_m = &mean[0];

			for(j=0; j < num_obs_eval; j++)
			{
			  R_CheckUserInterrupt();
				sum_y_ker = sum_y_sq_ker = 0.0;
				sum_ker = DBL_MIN;

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_y_ker_deriv[l] = sum_ker_deriv[l] = 0.0;
				}

				for(i=0; i < num_obs_train; i++)
				{

					/* Kernel for mean */

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg,matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg,matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Kernels for derivatives */

					for (k = 0; k < num_reg_continuous; k++)
					{

						prod_kernel_deriv[k] = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][i])/ipow(matrix_bandwidth_deriv[k][i],2);

						for (l = 0; l < num_reg_continuous; l++)
						{

							if(l != k)
							{
								prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
							}

						}

					}

					/* Mean */

					sum_ker += prod_kernel;
					sum_y_ker += vector_Y[i]*prod_kernel;
					sum_y_sq_ker += ipow(vector_Y[i],2)*prod_kernel;

					/* Gradient */

					for(l = 0; l < num_reg_continuous; l++)
					{
						sum_ker_deriv[l] += prod_kernel_deriv[l];
						sum_y_ker_deriv[l] += vector_Y[i]*prod_kernel_deriv[l];
					}

				}

				*pointer_m = sum_y_ker/sum_ker;

				temp_var = (sum_y_sq_ker/sum_ker) - ipow(*pointer_m++, 2);

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / sum_ker) : sqrt(temp_var / sum_ker));

				/* gradient[0][] is that for _first_ continuous variable */
				/* Using bandwidth defined for j not i */

				for(l = 0; l < num_reg_continuous; l++)
				{
					gradient[l][j] = (sum_y_ker_deriv[l]-(sum_y_ker/sum_ker)*sum_ker_deriv[l])/sum_ker;

					/* 11/28/01 - removed bw from denom -alloc probs */

					gradient_stderr[l][j] = sqrt(temp_var * DIFF_KER_PPM / sum_ker);

				}

			}

		}

	}
	else
	{

		/* Local linear */

		temp_mean_y = meand(num_obs_train, vector_Y);

		XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
		XTKYSQ = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
		DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

		/* Conduct the estimation */

		if(BANDWIDTH_reg == 0)
		{

			for (j = 0; j < num_obs_eval; j++)
			{
			  R_CheckUserInterrupt();
				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					XTKYSQ[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs_train; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel_cont = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][0]);
					}

					prod_kernel_cat = 1.0;

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));
					XTKYSQ[0][0] += ipow(*pointer_yi, 2) * prod_kernel;

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;
						XTKYSQ[k+1][0] += *pointer_yi * temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical()", j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

				mean[j] = DELTA[0][0];

				for(k = 0; k < num_reg_cat_cont; k++)
				{
					gradient[k][j] = - DELTA[k+1][0];
				}

				/*				DELTA =  mat_mul( XTKXINV, XTKYSQ, DELTA);

				temp_var = DELTA[0][0] - mean[j]*mean[j]; 12/9/03 - local
				linear E[Y^2|]-(E[Y|X])^2 for conditional variance is
				flawed. Now using local constant.*/

				temp_var = XTKYSQ[0][0]/XTKX[0][0] - ipow(XTKY[0][0]/XTKX[0][0],2);

				if(temp_var <= 0.0)
				{
					temp_var = DBL_MIN;
				}

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / XTKX[0][0]) : sqrt(temp_var / XTKX[0][0]));

				for(k=0; k < num_reg_cat_cont; k++)
				{
					if(k < num_reg_continuous)
					{
						gradient_stderr[k][j] = sqrt(temp_var * DIFF_KER_PPM /
							(XTKX[0][0] * ipow(matrix_bandwidth[k][0],2)));
					}
					else
					{
						gradient_stderr[k][j] = sqrt(temp_var * DIFF_KER_PPM /
							XTKX[0][0]);
					}
				}

			}

		}
		else if(BANDWIDTH_reg == 1)
		{

			for (j = 0; j < num_obs_eval; j++)
			{
			  R_CheckUserInterrupt();
				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					XTKYSQ[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs_train; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel_cont = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][j]);
					}

					prod_kernel_cat = 1.0;

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));
					XTKYSQ[0][0] += *pointer_yi * temp;

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;
						XTKYSQ[k+1][0] += *pointer_yi * temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical()", j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

				mean[j] = DELTA[0][0];

				for(k = 0; k < num_reg_cat_cont; k++)
				{
					gradient[k][j] = - DELTA[k+1][0];
				}

				/*				DELTA =  mat_mul( XTKXINV, XTKYSQ, DELTA);

				temp_var = DELTA[0][0] - mean[j]*mean[j]; 12/9/03 - local
				linear E[Y^2|]-(E[Y|X])^2 for conditional variance is
				flawed. Now using local constant.*/

				temp_var = XTKYSQ[0][0]/XTKX[0][0] - ipow(XTKY[0][0]/XTKX[0][0],2);

				if(temp_var <= 0.0)
				{
					temp_var = DBL_MIN;
				}

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / XTKX[0][0]) : sqrt(temp_var / XTKX[0][0]));

				for(k=0; k < num_reg_cat_cont; k++)
				{
					if(k < num_reg_continuous)
					{
						gradient_stderr[k][j] = sqrt(temp_var * DIFF_KER_PPM /
							(XTKX[0][0] * ipow(matrix_bandwidth[k][j],2)));
					}
					else
					{
						gradient_stderr[k][j] = sqrt(temp_var * DIFF_KER_PPM /
							XTKX[0][0]);
					}
				}

			}

		}
		else
		{

			for (j = 0; j < num_obs_eval; j++)
			{
			  R_CheckUserInterrupt();
				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					XTKYSQ[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs_train; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel_cont = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{

						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][i]);

					}

					prod_kernel_cat = 1.0;

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));
					XTKYSQ[0][0] += *pointer_yi * temp;

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
								* prod_kernel;
						}
						else if(l < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;
						XTKYSQ[k+1][0] += *pointer_yi * temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical()", j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

				mean[j] = DELTA[0][0];

				for(k = 0; k < num_reg_cat_cont; k++)
				{
					gradient[k][j] = - DELTA[k+1][0];
				}

				/*				DELTA =  mat_mul( XTKXINV, XTKYSQ, DELTA);

				temp_var = DELTA[0][0] - mean[j]*mean[j]; 12/9/03 - local
				linear E[Y^2|]-(E[Y|X])^2 for conditional variance is
				flawed. Now using local constant.*/

				temp_var = XTKYSQ[0][0]/XTKX[0][0] - ipow(XTKY[0][0]/XTKX[0][0],2);

				if(temp_var <= 0.0)
				{
					temp_var = DBL_MIN;
				}

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / XTKX[0][0]) : sqrt(temp_var / XTKX[0][0]));

				for(k=0; k < num_reg_cat_cont; k++)
				{
					/* Divide by k[j] since k[i] not avail... check N for correct way */
					if(k < num_reg_continuous)
					{
						gradient_stderr[k][j] = sqrt(temp_var * DIFF_KER_PPM /
							(XTKX[0][0] * ipow(matrix_bandwidth[k][j],2)));
					}
					else
					{
						gradient_stderr[k][j] = sqrt(temp_var * DIFF_KER_PPM /
							XTKX[0][0]);
					}
				}
			}

		}

		mat_free( XTKX );
		mat_free( XTKXINV );
		mat_free( XTKY );
		mat_free( XTKYSQ );
		mat_free( DELTA );

	}
	#endif

	#ifdef MPI2

	if(int_ll == 0)
	{

		/* Conduct Nadaraya-Watson estimation */

		if(BANDWIDTH_reg == 0)
		{

			for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
			{

				sum_y_ker = sum_y_sq_ker = 0.0;
				sum_ker = DBL_MIN;

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_y_ker_deriv[l] = sum_ker_deriv[l] = 0.0;
				}

				for(i=0; i < num_obs_train; i++)
				{

					/* Kernel for mean */

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg,matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Kernels for derivatives */

					for (k = 0; k < num_reg_continuous; k++)
					{

						prod_kernel_deriv[k] = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][0]);

						for (l = 0; l < num_reg_continuous; l++)
						{

							if(l != k)
							{
								prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
							}

						}

					}

					/* Mean */

					sum_ker += prod_kernel;
					sum_y_ker += vector_Y[i]*prod_kernel;
					sum_y_sq_ker += ipow(vector_Y[i],2)*prod_kernel;

					/* Gradient */

					for(l = 0; l < num_reg_continuous; l++)
					{
						sum_ker_deriv[l] += prod_kernel_deriv[l];
						sum_y_ker_deriv[l] += vector_Y[i]*prod_kernel_deriv[l];
					}

				}

				mean[j-my_rank*stride] = sum_y_ker/sum_ker;
				temp_var = (sum_y_sq_ker/sum_ker) - ipow(mean[j-my_rank*stride],2);

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / sum_ker) : sqrt(temp_var / sum_ker));

				/* gradient[0][] is that for _first_ continuous variable */

				for(l = 0; l < num_reg_continuous; l++)
				{
					gradient[l][j-my_rank*stride] = (sum_y_ker_deriv[l]-(sum_y_ker/sum_ker)*sum_ker_deriv[l])/(matrix_bandwidth_deriv[l][0]*sum_ker);

					gradient_stderr[l][j-my_rank*stride] = sqrt(temp_var * DIFF_KER_PPM /
						(sum_ker * ipow(matrix_bandwidth_deriv[l][0],2)));
				}

			}

		}
		else if(BANDWIDTH_reg == 1)
		{

			for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
			{

				sum_y_ker = sum_y_sq_ker = 0.0;
				sum_ker = DBL_MIN;

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_y_ker_deriv[l] = sum_ker_deriv[l] = 0.0;
				}

				for(i=0; i < num_obs_train; i++)
				{

					/* Kernel for mean */

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg,matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Kernels for derivatives */

					for (k = 0; k < num_reg_continuous; k++)
					{

						prod_kernel_deriv[k] = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][j]);

						for (l = 0; l < num_reg_continuous; l++)
						{

							if(l != k)
							{
								prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
							}

						}

					}

					/* Mean */

					sum_ker += prod_kernel;
					sum_y_ker += vector_Y[i]*prod_kernel;
					sum_y_sq_ker += ipow(vector_Y[i],2)*prod_kernel;

					/* Gradient */

					for(l = 0; l < num_reg_continuous; l++)
					{
						sum_ker_deriv[l] += prod_kernel_deriv[l];
						sum_y_ker_deriv[l] += vector_Y[i]*prod_kernel_deriv[l];
					}

				}

				mean[j-my_rank*stride] = sum_y_ker/sum_ker;

				temp_var = (sum_y_sq_ker/sum_ker) - ipow(*pointer_m++, 2);

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / sum_ker) : sqrt(temp_var / sum_ker));

				/* gradient[0][] is that for _first_ continuous variable */

				for(l = 0; l < num_reg_continuous; l++)
				{
					gradient[l][j-my_rank*stride] = (sum_y_ker_deriv[l]-(sum_y_ker/sum_ker)*sum_ker_deriv[l])/(matrix_bandwidth_deriv[l][j]*sum_ker);

					gradient_stderr[l][j-my_rank*stride] = sqrt(temp_var * DIFF_KER_PPM /
						(sum_ker * ipow(matrix_bandwidth_deriv[l][j],2)));
				}

			}

		}
		else
		{

			for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
			{

				sum_y_ker = sum_y_sq_ker = 0.0;
				sum_ker = DBL_MIN;

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_y_ker_deriv[l] = sum_ker_deriv[l] = 0.0;
				}

				for(i=0; i < num_obs_train; i++)
				{

					/* Kernel for mean */

					prod_kernel_cont = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
					}

					prod_kernel_cat = 1.0;

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg,matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg,matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Kernels for derivatives */

					for (k = 0; k < num_reg_continuous; k++)
					{

						prod_kernel_deriv[k] = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][i])/ipow(matrix_bandwidth_deriv[k][i],2);

						for (l = 0; l < num_reg_continuous; l++)
						{

							if(l != k)
							{
								prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
							}

						}

					}

					/* Mean */

					sum_ker += prod_kernel;
					sum_y_ker += vector_Y[i]*prod_kernel;
					sum_y_sq_ker += ipow(vector_Y[i],2)*prod_kernel;

					/* Gradient */

					for(l = 0; l < num_reg_continuous; l++)
					{
						sum_ker_deriv[l] += prod_kernel_deriv[l];
						sum_y_ker_deriv[l] += vector_Y[i]*prod_kernel_deriv[l];
					}

				}

				mean[j-my_rank*stride] = sum_y_ker/sum_ker;

				temp_var = (sum_y_sq_ker/sum_ker) - ipow(*pointer_m++, 2);

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / sum_ker) : sqrt(temp_var / sum_ker));

				/* gradient[0][] is that for _first_ continuous variable */
				/* Using bandwidth defined for j not i */

				for(l = 0; l < num_reg_continuous; l++)
				{
					gradient[l][j-my_rank*stride] = (sum_y_ker_deriv[l]-(sum_y_ker/sum_ker)*sum_ker_deriv[l])/sum_ker;

					/* 11/28/01 - removed bw from denom -alloc probs */

					gradient_stderr[l][j-my_rank*stride] = sqrt(temp_var * DIFF_KER_PPM / sum_ker);

				}

			}

		}

	}
	else
	{

		/* Local linear */

		temp_mean_y = meand(num_obs_train, vector_Y);

		XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
		XTKYSQ = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
		DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

		/* Conduct the estimation */

		if(BANDWIDTH_reg == 0)
		{

			for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
			{
				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					XTKYSQ[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs_train; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel_cont = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][0]);
					}

					prod_kernel_cat = 1.0;

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));
					XTKYSQ[0][0] += *pointer_yi * temp;

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;
						XTKYSQ[k+1][0] += *pointer_yi * temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical()", j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

				mean[j-my_rank*stride] = DELTA[0][0];

				for(k = 0; k < num_reg_cat_cont; k++)
				{
					gradient[k][j-my_rank*stride] = - DELTA[k+1][0];
				}

				/*				DELTA =  mat_mul( XTKXINV, XTKYSQ, DELTA);

				temp_var = DELTA[0][0] - mean[j]*mean[j]; 12/9/03 - local
				linear E[Y^2|]-(E[Y|X])^2 for conditional variance is
				flawed. Now using local constant.*/

				temp_var = XTKYSQ[0][0]/XTKX[0][0] - ipow(XTKY[0][0]/XTKX[0][0],2);

				if(temp_var <= 0.0)
				{
					temp_var = DBL_MIN;
				}

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / XTKX[0][0]) : sqrt(temp_var / XTKX[0][0]));

				for(k=0; k < num_reg_cat_cont; k++)
				{
					if(k < num_reg_continuous)
					{
						gradient_stderr[k][j-my_rank*stride] = sqrt(temp_var * DIFF_KER_PPM /
							(XTKX[0][0] * ipow(matrix_bandwidth[k][0],2)));
					}
					else
					{
						gradient_stderr[k][j-my_rank*stride] = sqrt(temp_var * DIFF_KER_PPM /
							XTKX[0][0]);
					}
				}

			}

		}
		else if(BANDWIDTH_reg == 1)
		{

			for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
			{

				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					XTKYSQ[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs_train; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel_cont = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{
						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][j]);
					}

					prod_kernel_cat = 1.0;

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));
					XTKYSQ[0][0] += *pointer_yi * temp;

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;
						XTKYSQ[k+1][0] += *pointer_yi * temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical()", j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

				mean[j-my_rank*stride] = DELTA[0][0];

				for(k = 0; k < num_reg_cat_cont; k++)
				{
					gradient[k][j-my_rank*stride] = - DELTA[k+1][0];
				}

				/*				DELTA =  mat_mul( XTKXINV, XTKYSQ, DELTA);

				temp_var = DELTA[0][0] - mean[j]*mean[j]; 12/9/03 - local
				linear E[Y^2|]-(E[Y|X])^2 for conditional variance is
				flawed. Now using local constant.*/

				temp_var = XTKYSQ[0][0]/XTKX[0][0] - ipow(XTKY[0][0]/XTKX[0][0],2);

				if(temp_var <= 0.0)
				{
					temp_var = DBL_MIN;
				}

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / XTKX[0][0]) : sqrt(temp_var / XTKX[0][0]));

				for(k=0; k < num_reg_cat_cont; k++)
				{
					if(k < num_reg_continuous)
					{
						gradient_stderr[k][j-my_rank*stride] = sqrt(temp_var * DIFF_KER_PPM /
							(XTKX[0][0] * ipow(matrix_bandwidth[k][j],2)));
					}
					else
					{
						gradient_stderr[k][j-my_rank*stride] = sqrt(temp_var * DIFF_KER_PPM /
							XTKX[0][0]);
					}
				}

			}

		}
		else
		{

			for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
			{

				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					XTKYSQ[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs_train; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel_cont = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{

						prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][i]);

					}

					prod_kernel_cat = 1.0;

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel_cat *= kernel_unordered_ratio(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
					}

					prod_kernel = prod_kernel_cont*prod_kernel_cat;

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));
					XTKYSQ[0][0] += *pointer_yi * temp;

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
								* prod_kernel;
						}
						else if(l < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;
						XTKYSQ[k+1][0] += *pointer_yi * temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical()", j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

				mean[j-my_rank*stride] = DELTA[0][0];

				for(k = 0; k < num_reg_cat_cont; k++)
				{
					gradient[k][j-my_rank*stride] = - DELTA[k+1][0];
				}

				/*				DELTA =  mat_mul( XTKXINV, XTKYSQ, DELTA);

				temp_var = DELTA[0][0] - mean[j]*mean[j]; 12/9/03 - local
				linear E[Y^2|]-(E[Y|X])^2 for conditional variance is
				flawed. Now using local constant.*/

				temp_var = XTKYSQ[0][0]/XTKX[0][0] - ipow(XTKY[0][0]/XTKX[0][0],2);

				if(temp_var <= 0.0)
				{
					temp_var = DBL_MIN;
				}

				/* With no continuous variables, need to drop K_INT_KERNEL_P */

				mean_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(temp_var * K_INT_KERNEL_P / XTKX[0][0]) : sqrt(temp_var / XTKX[0][0]));

				for(k=0; k < num_reg_cat_cont; k++)
				{
					/* Divide by k[j] since k[i] not avail... check N for correct way */
					if(k < num_reg_continuous)
					{
						gradient_stderr[k][j-my_rank*stride] = sqrt(temp_var * DIFF_KER_PPM /
							(XTKX[0][0] * ipow(matrix_bandwidth[k][j],2)));
					}
					else
					{
						gradient_stderr[k][j-my_rank*stride] = sqrt(temp_var * DIFF_KER_PPM /
							XTKX[0][0]);
					}
				}
			}

		}

		mat_free( XTKX );
		mat_free( XTKXINV );
		mat_free( XTKY );
		mat_free( XTKYSQ );
		mat_free( DELTA );

	}

	/* Important - only one gather permitted */

	/* Gather */

	MPI_Gather(mean, stride, MPI_DOUBLE, mean, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(mean, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(mean_stderr, stride, MPI_DOUBLE, mean_stderr, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(mean_stderr, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for(l = 0; l < num_reg_continuous; l++)
	{

		MPI_Gather(&gradient[l][0], stride, MPI_DOUBLE, &gradient[l][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&gradient[l][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&gradient_stderr[l][0], stride, MPI_DOUBLE, &gradient_stderr[l][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&gradient_stderr[l][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	}
	#endif

	if(vector_Y_eval)
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

	free(prod_kernel_deriv);
	free(sum_ker_deriv);
	free(sum_y_ker_deriv);

	free(lambda);

	free_mat(matrix_bandwidth,num_reg_continuous);
	free_mat(matrix_bandwidth_deriv,num_reg_continuous);

	return(0);

}


int kernel_estimate_regression_categorical_leave_one_out(
int int_ll,
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
int *num_categories,
double *mean)
{

	/* This function estimates a leave-one-out Nadaraya-Watson regression */
	/* function using both continuous and categorical covariates with three */
	/* estimation techniques and an assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int k;
	int l;

	const double epsilon = 1.0/num_obs;
  double nepsilon;

	double prod_kernel;

	double sum_ker;
	double sum_y_ker;

	double *lambda;
	double **matrix_bandwidth;

	double temp;
	double temp1 = DBL_MAX;

	MATRIX  XTKX;
	MATRIX  XTKXINV;
	MATRIX  XTKY;
	MATRIX  DELTA;

	double *pointer_yi;
	double *pointer_m;

	int num_reg_cat_cont;

	#ifdef MPI2
	int stride = ceil((double) num_obs / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	if(int_TAYLOR == 1)
	{
		num_reg_cat_cont = num_reg_unordered + num_reg_ordered + num_reg_continuous;
	}
	else
	{
		num_reg_cat_cont = num_reg_continuous;
	}

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);
	matrix_bandwidth = alloc_matd(num_obs,num_reg_continuous);

	#ifndef MPI2

	/* Conduct the estimation */

	if(int_ll == 0)
	{

		/* Nadaraya-Watson */

		/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

		if(kernel_bandwidth_mean(
			KERNEL_reg,
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
			matrix_X_continuous,			 /* Not used */
			matrix_X_continuous,			 /* Not used */
			matrix_X_continuous,
			matrix_X_continuous,
			matrix_bandwidth,					 /* Not used */
			matrix_bandwidth,
			lambda)==1)
		{

			free(lambda);
			free_mat(matrix_bandwidth,num_reg_continuous);

			return(1);
		}

		if(BANDWIDTH_reg == 0)
		{

			pointer_m = &mean[0];

			for(j=0; j < num_obs; j++)
			{

				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						prod_kernel = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][0]);
						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
						}

						sum_ker += prod_kernel;
						sum_y_ker += *pointer_yi*prod_kernel;

					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					*pointer_m++ = sum_y_ker/sum_ker;
				}
				else
				{
					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_leave_one_out()",j);
					}
					return(1);
				}

			}

		}
		else if(BANDWIDTH_reg == 1)
		{

			pointer_m = &mean[0];

			for(j=0; j < num_obs; j++)
			{

				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						prod_kernel = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][j]);
						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
						}

						sum_ker += prod_kernel;
						sum_y_ker += *pointer_yi*prod_kernel;

					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					*pointer_m++ = sum_y_ker/sum_ker;
				}
				else
				{
					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_leave_one_out()",j);
					}
					return(1);
				}

			}

		}
		else
		{

			pointer_m = &mean[0];

			for(j=0; j < num_obs; j++)
			{

				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						prod_kernel = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
						}

						sum_ker += prod_kernel;
						sum_y_ker += *pointer_yi*prod_kernel;

					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					*pointer_m++ = sum_y_ker/sum_ker;
				}
				else
				{
					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_leave_one_out()",j);
					}
					return(1);
				}

			}

		}

	}
	else
	{

		/* Local Linear */

		XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
		DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

		/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

		if(kernel_bandwidth_mean(
			KERNEL_reg,
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
			matrix_X_continuous,			 /* Not used */
			matrix_X_continuous,			 /* Not used */
			matrix_X_continuous,
			matrix_X_continuous,
			matrix_bandwidth,					 /* Not used */
			matrix_bandwidth,
			lambda)==1)
		{

			mat_free( XTKX );
			mat_free( XTKXINV );
			mat_free( XTKY );
			mat_free( DELTA );

			free(lambda);
			free_mat(matrix_bandwidth,num_reg_continuous);

			return(1);

		}

		/* Conduct the estimation */

		if(BANDWIDTH_reg == 0)
		{

			pointer_m = &mean[0];

			for (j = 0; j < num_obs; j++)
			{
				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						/* Matrix is rows/cols... not C convention... */

						prod_kernel = 1.0;

						for(k = 0; k < num_reg_continuous; k++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][0]);
						}

						for(k = 0; k < num_reg_unordered; k++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
						}

						for(k = 0; k < num_reg_ordered; k++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
						}

						/* Upper left block */

						XTKX[0][0] += prod_kernel;

						/* First element of XTKY */

						XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* First lower column of XTKX */

							if(k < num_reg_continuous)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
									* prod_kernel;
							}

							/* Diagonal of lower block of XTKX */

							XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

							/* Remaining elements of XTKY */

							XTKY[k+1][0] += temp1 * temp;

							/* Take advantage of symmetric nature of XTKX */

							for(l=0; l < k; l++)
							{
								if(l < num_reg_continuous)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
										* prod_kernel;
								}
							}

						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_leave_one_out()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
				*pointer_m++ = DELTA[0][0];

			}

		}
		else if(BANDWIDTH_reg == 1)
		{

			pointer_m = &mean[0];

			for (j = 0; j < num_obs; j++)
			{
				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						/* Matrix is rows/cols... not C convention... */

						prod_kernel = 1.0;

						for(k = 0; k < num_reg_continuous; k++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][j]);
						}

						for(k = 0; k < num_reg_unordered; k++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
						}

						for(k = 0; k < num_reg_ordered; k++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
						}

						/* Upper left block */

						XTKX[0][0] += prod_kernel;

						/* First element of XTKY */

						XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* First lower column of XTKX */

							if(k < num_reg_continuous)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
									* prod_kernel;
							}

							/* Diagonal of lower block of XTKX */

							XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

							/* Remaining elements of XTKY */

							XTKY[k+1][0] += temp1 * temp;

							/* Take advantage of symmetric nature of XTKX */

							for(l=0; l < k; l++)
							{
								if(l < num_reg_continuous)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
										* prod_kernel;
								}
							}

						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_leave_one_out()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
				*pointer_m++ = DELTA[0][0];

			}

		}
		else
		{

			pointer_m = &mean[0];

			for (j = 0; j < num_obs; j++)
			{
				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						/* Matrix is rows/cols... not C convention... */

						prod_kernel = 1.0;

						for(k = 0; k < num_reg_continuous; k++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][i]);
						}

						for(k = 0; k < num_reg_unordered; k++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
						}

						for(k = 0; k < num_reg_ordered; k++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
						}

						/* Upper left block */

						XTKX[0][0] += prod_kernel;

						/* First element of XTKY */

						XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* First lower column of XTKX */

							if(k < num_reg_continuous)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
									* prod_kernel;
							}

							/* Diagonal of lower block of XTKX */

							XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

							/* Remaining elements of XTKY */

							XTKY[k+1][0] += temp1 * temp;

							/* Take advantage of symmetric nature of XTKX */

							for(l=0; l < k; l++)
							{
								if(l < num_reg_continuous)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
										* prod_kernel;
								}
							}

						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_leave_one_out()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
				*pointer_m++ = DELTA[0][0];

			}

		}

		mat_free( XTKX );
		mat_free( XTKXINV );
		mat_free( XTKY );
		mat_free( DELTA );

	}
	#endif

	#ifdef MPI2

	/* Conduct the estimation - MPI-enables */

	if(int_ll == 0)
	{

		/* Nadaraya-Watson */

		/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

		/* May be redundant */
		if(kernel_bandwidth_mean(
			KERNEL_reg,
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
			/* Not used */
			matrix_X_continuous,
			/* Not used */
			matrix_X_continuous,
			matrix_X_continuous,
			matrix_X_continuous,
			/* Not used */
			matrix_bandwidth,
			matrix_bandwidth,
			lambda)==1)
		{

			free(lambda);
			free_mat(matrix_bandwidth,num_reg_continuous);

			return(1);
		}

		if(BANDWIDTH_reg == 0)
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						prod_kernel = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][0]);
						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
						}

						sum_ker += prod_kernel;
						sum_y_ker += *pointer_yi*prod_kernel;

					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					mean[j-my_rank*stride] = sum_y_ker/sum_ker;
				}
				else
				{
					/* Don't print if you are using MPI */
					if((int_DEBUG == 1)&&(my_rank == 0))
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_leave_one_out()",j);
					}
					return(1);
				}

			}

		}
		else if(BANDWIDTH_reg == 1)
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						prod_kernel = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][j]);
						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
						}

						sum_ker += prod_kernel;
						sum_y_ker += *pointer_yi*prod_kernel;

					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					mean[j-my_rank*stride] = sum_y_ker/sum_ker;
				}
				else
				{
					/* Don't print if you are using MPI */
					if((my_rank == 0)&&(int_DEBUG == 1))
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_leave_one_out()",j);
					}
					return(1);
				}

			}

		}
		else
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{
				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						prod_kernel = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
						}

						sum_ker += prod_kernel;
						sum_y_ker += *pointer_yi*prod_kernel;

					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					mean[j-my_rank*stride] = sum_y_ker/sum_ker;
				}
				else
				{
					/* Don't print if you are using MPI */
					if((my_rank == 0)&&(int_DEBUG == 1))
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_leave_one_out()",j);
					}
					return(1);
				}

			}

		}

	}
	else
	{

		/* Local Linear */

		XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
		DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

		/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

		if(kernel_bandwidth_mean(
			KERNEL_reg,
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
			/* Not used */
			matrix_X_continuous,
			/* Not used */
			matrix_X_continuous,
			matrix_X_continuous,
			matrix_X_continuous,
			/* Not used */
			matrix_bandwidth,
			matrix_bandwidth,
			lambda)==1)
		{

			mat_free( XTKX );
			mat_free( XTKXINV );
			mat_free( XTKY );
			mat_free( DELTA );

			free(lambda);
			free_mat(matrix_bandwidth,num_reg_continuous);

			return(1);

		}

		/* Conduct the estimation */

		if(BANDWIDTH_reg == 0)
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						/* Matrix is rows/cols... not C convention... */

						prod_kernel = 1.0;

						for(k = 0; k < num_reg_continuous; k++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][0]);
						}

						for(k = 0; k < num_reg_unordered; k++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
						}

						for(k = 0; k < num_reg_ordered; k++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
						}

						/* Upper left block */

						XTKX[0][0] += prod_kernel;

						/* First element of XTKY */

						XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* First lower column of XTKX */

							if(k < num_reg_continuous)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
									* prod_kernel;
							}

							/* Diagonal of lower block of XTKX */

							XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

							/* Remaining elements of XTKY */

							XTKY[k+1][0] += temp1 * temp;

							/* Take advantage of symmetric nature of XTKX */

							for(l=0; l < k; l++)
							{
								if(l < num_reg_continuous)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
										* prod_kernel;
								}
							}

						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if((my_rank == 0)&&(int_DEBUG == 1))
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_leave_one_out()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/


					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

				mean[j-my_rank*stride] =  DELTA[0][0];

			}

		}
		else if(BANDWIDTH_reg == 1)
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						/* Matrix is rows/cols... not C convention... */

						prod_kernel = 1.0;

						for(k = 0; k < num_reg_continuous; k++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][j]);
						}

						for(k = 0; k < num_reg_unordered; k++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
						}

						for(k = 0; k < num_reg_ordered; k++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
						}

						/* Upper left block */

						XTKX[0][0] += prod_kernel;

						/* First element of XTKY */

						XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* First lower column of XTKX */

							if(k < num_reg_continuous)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
									* prod_kernel;
							}

							/* Diagonal of lower block of XTKX */

							XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

							/* Remaining elements of XTKY */

							XTKY[k+1][0] += temp1 * temp;

							/* Take advantage of symmetric nature of XTKX */

							for(l=0; l < k; l++)
							{
								if(l < num_reg_continuous)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
										* prod_kernel;
								}
							}

						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if((my_rank == 0)&&(int_DEBUG == 1))
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_leave_one_out()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

				mean[j-my_rank*stride] =  DELTA[0][0];

			}

		}
		else
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					if(i != j)
					{

						/* Matrix is rows/cols... not C convention... */

						prod_kernel = 1.0;

						for(k = 0; k < num_reg_continuous; k++)
						{
							prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][i]);
						}

						for(k = 0; k < num_reg_unordered; k++)
						{
							prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
						}

						for(k = 0; k < num_reg_ordered; k++)
						{
							prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
						}

						/* Upper left block */

						XTKX[0][0] += prod_kernel;

						/* First element of XTKY */

						XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* First lower column of XTKX */

							if(k < num_reg_continuous)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
									* prod_kernel;
							}
							else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
									* prod_kernel;
							}

							/* Diagonal of lower block of XTKX */

							XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

							/* Remaining elements of XTKY */

							XTKY[k+1][0] += temp1 * temp;

							/* Take advantage of symmetric nature of XTKX */

							for(l=0; l < k; l++)
							{
								if(l < num_reg_continuous)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
										* prod_kernel;
								}
								else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
										* prod_kernel;
								}
							}

						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if((my_rank == 0)&&(int_DEBUG == 1))
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_leave_one_out()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

				mean[j-my_rank*stride] =  DELTA[0][0];

			}

		}

		mat_free( XTKX );
		mat_free( XTKXINV );
		mat_free( XTKY );
		mat_free( DELTA );

	}

	/* Gather */

	MPI_Gather(mean, stride, MPI_DOUBLE, mean, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(mean, num_obs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	free(lambda);
	free_mat(matrix_bandwidth,num_reg_continuous);

	return(0);

}


int kernel_estimate_regression_categorical_no_stderr(
int int_compute_gradient,
int int_ll,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_reg,
int int_WEIGHTS,
int *var_index_int,
int num_var_test_int,
double **matrix_weights_K,
double ***matrix_weights_K_deriv,
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
double **matrix_bandwidth,
double **matrix_bandwidth_deriv,
double *vector_Y,
double *lambda,
int *num_categories,
double *mean,
double **gradient)
{

	/* This function estimates a Nadaraya-Watson or Local Linear regression */
	/* function using both continuous and categorical covariates with three */
	/* estimation techniques and an assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int k;
	int tmp_k;
	int l;

	const double epsilon = 1.0/num_obs_train;
  double nepsilon;

	double prod_kernel;

	double prod_kernel_cat;
	double prod_kernel_cont;

	double sum_ker;
	double sum_y_ker;

	double *prod_kernel_deriv;

	double *p_prod_kernel_deriv;

	double *sum_ker_deriv;
	double *sum_y_ker_deriv;

	double *p_sum_ker_deriv;
	double *p_sum_y_ker_deriv;

	double sum_ker_deriv_scalar;
	double sum_y_ker_deriv_scalar;

	double temp_mean_y;

	double temp;
	double temp1 = DBL_MAX;

	double *pointer_yi;
	double *pointer_m;

	double *pointer_matrix_weights_K;
	double *pointer_matrix_weights_K_deriv;

	double *pointer_gradient;
	double *pointer_matrix_bandwidth;

	MATRIX  XTKX;
	MATRIX  XTKXINV;
	MATRIX  XTKY;
	MATRIX  DELTA;

	int num_reg_cat_cont;

	#ifdef MPI2
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	if(int_TAYLOR == 1)
	{
		num_reg_cat_cont = num_reg_unordered + num_reg_continuous;
	}
	else
	{
		num_reg_cat_cont = num_reg_continuous;
	}

	#ifndef MPI2

	if(int_compute_gradient == 1)
	{

		/* Gradient computed */

		if(int_WEIGHTS == 0)
		{

			/* Do not use weights */

			if(int_ll == 0)
			{

				prod_kernel_deriv = alloc_vecd(num_reg_continuous);

				sum_ker_deriv = alloc_vecd(num_reg_continuous);
				sum_y_ker_deriv = alloc_vecd(num_reg_continuous);

				/* Conduct Nadaraya-Watson estimation */

				if(BANDWIDTH_reg == 0)
				{

					pointer_m = &mean[0];

					for(j=0; j < num_obs_eval; j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							*p_sum_y_ker_deriv++ = *p_sum_ker_deriv++ = 0.0;
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi * prod_kernel;

							/* Kernels for derivatives */

							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for (k = 0; k < num_reg_continuous; k++)
							{

								*p_prod_kernel_deriv = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][0]);

								for (l = 0; l < num_reg_continuous; l++)
								{

									if(l != k)
									{
										*p_prod_kernel_deriv *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
									}

								}

								p_prod_kernel_deriv++;

							}

							/* Gradient */

							p_sum_ker_deriv = &sum_ker_deriv[0];
							p_sum_y_ker_deriv = &sum_y_ker_deriv[0];
							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for(l = 0; l < num_reg_continuous; l++)
							{
								*p_sum_ker_deriv++ += *p_prod_kernel_deriv;
								*p_sum_y_ker_deriv++ += *pointer_yi * *p_prod_kernel_deriv++;
							}

							pointer_yi++;

						}

						*pointer_m = sum_y_ker/sum_ker;

						/* gradient[0][] is that for _first_ continuous variable */

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							gradient[l][j] = (*p_sum_y_ker_deriv++ - *pointer_m * *p_sum_ker_deriv++)/(matrix_bandwidth_deriv[l][0]*sum_ker);
						}

						pointer_m++;

					}

				}
				else if(BANDWIDTH_reg == 1)
				{

					pointer_m = &mean[0];

					for(j=0; j < num_obs_eval; j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							*p_sum_y_ker_deriv++ = *p_sum_ker_deriv++ = 0.0;
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi * prod_kernel;

							/* Kernels for derivatives */

							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for (k = 0; k < num_reg_continuous; k++)
							{

								*p_prod_kernel_deriv = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][j]);

								for (l = 0; l < num_reg_continuous; l++)
								{

									if(l != k)
									{
										*p_prod_kernel_deriv *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
									}

								}

								p_prod_kernel_deriv++;

							}

							/* Gradient */

							p_sum_ker_deriv = &sum_ker_deriv[0];
							p_sum_y_ker_deriv = &sum_y_ker_deriv[0];
							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for(l = 0; l < num_reg_continuous; l++)
							{
								*p_sum_ker_deriv++ += *p_prod_kernel_deriv;
								*p_sum_y_ker_deriv++ += *pointer_yi * *p_prod_kernel_deriv++;
							}

							pointer_yi++;

						}

						*pointer_m = sum_y_ker/sum_ker;

						/* gradient[0][] is that for _first_ continuous variable */

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							gradient[l][j] = (*p_sum_y_ker_deriv++ - *pointer_m * *p_sum_ker_deriv++)/(matrix_bandwidth_deriv[l][j]*sum_ker);
						}

						pointer_m++;

					}

				}
				else
				{

					pointer_m = &mean[0];

					for(j=0; j < num_obs_eval; j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							*p_sum_y_ker_deriv++ = *p_sum_ker_deriv++ = 0.0;
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Kernels for derivatives */

							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for (k = 0; k < num_reg_continuous; k++)
							{

								*p_prod_kernel_deriv = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][i])/ipow(matrix_bandwidth_deriv[k][i],2);

								for (l = 0; l < num_reg_continuous; l++)
								{

									if(l != k)
									{
										*p_prod_kernel_deriv *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
									}

									p_prod_kernel_deriv++;

								}

							}

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi++*prod_kernel;

							/* Gradient */

							p_sum_ker_deriv = &sum_ker_deriv[0];
							p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for(l = 0; l < num_reg_continuous; l++)
							{
								*p_sum_ker_deriv++ += *p_prod_kernel_deriv;
								*p_sum_y_ker_deriv++ += *pointer_yi * *p_prod_kernel_deriv++;
							}

							pointer_yi++;

						}

						*pointer_m = sum_y_ker/sum_ker;

						/* gradient[0][] is that for _first_ continuous variable */

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							gradient[l][j] = (*p_sum_y_ker_deriv++ - *pointer_m * *p_sum_ker_deriv++)/sum_ker;
						}

						pointer_m++;
					}

				}

				free(prod_kernel_deriv);
				free(sum_ker_deriv);
				free(sum_y_ker_deriv);

			}
			else
			{

				/* Local linear */

				temp_mean_y = meand(num_obs_train, vector_Y);

				XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
				DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

				/* Conduct the estimation */

				if(BANDWIDTH_reg == 0)
				{

					pointer_m = &mean[0];

					for (j = 0; j < num_obs_eval; j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;
							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][0]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */

								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if(int_DEBUG == 1)
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}

							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
              XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
						*pointer_m++ = DELTA[0][0];

						for(k = 0; k < num_reg_cat_cont; k++)
						{
							gradient[k][j] = - DELTA[k+1][0];
						}

					}

				}
				else if(BANDWIDTH_reg == 1)
				{

					pointer_m = &mean[0];

					for (j = 0; j < num_obs_eval; j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;
							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][j]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */

								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if(int_DEBUG == 1)
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}
							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
              XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
						*pointer_m++ = DELTA[0][0];

						for(k = 0; k < num_reg_cat_cont; k++)
						{
							gradient[k][j] = - DELTA[k+1][0];
						}

					}

				}
				else
				{

					pointer_m = &mean[0];

					for (j = 0; j < num_obs_eval; j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;

							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][i]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(k = 0; k < num_reg_ordered; k++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */

								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if(int_DEBUG == 1)
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}

							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
              XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

						*pointer_m++ = DELTA[0][0];

						for(k = 0; k < num_reg_cat_cont; k++)
						{
							gradient[k][j] = - DELTA[k+1][0];
						}

					}

				}

				mat_free( XTKX );
				mat_free( XTKXINV );
				mat_free( XTKY );
				mat_free( DELTA );

			}

		}
		else
		{

			/* Use weights, compute gradients */

			if(int_ll == 0)
			{

				/* Nadaraya Watson */

				if(BANDWIDTH_reg == 0)
				{

					/* Mean using weight pointers */

					pointer_m = &mean[0];

					for(j=0; j < num_obs_eval; j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						pointer_yi = &vector_Y[0];
						pointer_matrix_weights_K = &matrix_weights_K[j][0];

						for (i =  0; i <  num_obs_train; i++)
						{

							sum_ker += *pointer_matrix_weights_K;
							sum_y_ker += *pointer_yi++ * *pointer_matrix_weights_K++;

						}

						*pointer_m++ = sum_y_ker/sum_ker;

					}

					/* Gradient using weight pointers - only compute necessary */

					for(k=0; k < num_var_test_int; k++)
					{
						/* Continuous variable */
						tmp_k = var_index_int[k] - num_reg_unordered;
						if(var_index_int[k] >= num_reg_unordered)
						{
							pointer_m = &mean[0];
							pointer_gradient = &gradient[tmp_k][0];
							pointer_matrix_bandwidth = &matrix_bandwidth_deriv[tmp_k][0];

							for(j=0; j < num_obs_eval; j++)
							{

								sum_ker_deriv_scalar = sum_y_ker_deriv_scalar = 0.0;
								sum_ker = DBL_MIN;

								pointer_yi = &vector_Y[0];
								pointer_matrix_weights_K = &matrix_weights_K[j][0];
								pointer_matrix_weights_K_deriv = &matrix_weights_K_deriv[tmp_k][j][0];

								for (i =  0; i <  num_obs_train; i++)
								{
									sum_ker += *pointer_matrix_weights_K++;
									sum_ker_deriv_scalar += *pointer_matrix_weights_K_deriv;
									sum_y_ker_deriv_scalar += *pointer_yi++ * *pointer_matrix_weights_K_deriv++;
								}

								*pointer_gradient++ = (sum_y_ker_deriv_scalar - *pointer_m++ * sum_ker_deriv_scalar)/(*pointer_matrix_bandwidth * sum_ker);

							}
						}
					}

				}
				else if(BANDWIDTH_reg == 1)
				{

					/* Mean using weight pointers */

					pointer_m = &mean[0];

					for(j=0; j < num_obs_eval; j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						pointer_yi = &vector_Y[0];
						pointer_matrix_weights_K = &matrix_weights_K[j][0];

						for (i =  0; i <  num_obs_train; i++)
						{

							sum_ker += *pointer_matrix_weights_K;
							sum_y_ker += *pointer_yi++ * *pointer_matrix_weights_K++;

						}

						*pointer_m++ = sum_y_ker/sum_ker;

					}

					/* Gradient using weight pointers - only compute necessary */

					for(k=0; k < num_var_test_int; k++)
					{
						/* Continuous variable */
						tmp_k = var_index_int[k] - num_reg_unordered;
						if(var_index_int[k] >= num_reg_unordered)
						{
							pointer_m = &mean[0];
							pointer_gradient = &gradient[tmp_k][0];
							pointer_matrix_bandwidth = &matrix_bandwidth_deriv[tmp_k][0];

							for(j=0; j < num_obs_eval; j++)
							{

								sum_ker_deriv_scalar = sum_y_ker_deriv_scalar = 0.0;
								sum_ker = DBL_MIN;

								pointer_yi = &vector_Y[0];
								pointer_matrix_weights_K = &matrix_weights_K[j][0];
								pointer_matrix_weights_K_deriv = &matrix_weights_K_deriv[tmp_k][j][0];

								for (i =  0; i <  num_obs_train; i++)
								{
									sum_ker += *pointer_matrix_weights_K++;
									sum_ker_deriv_scalar += *pointer_matrix_weights_K_deriv;
									sum_y_ker_deriv_scalar += *pointer_yi++ * *pointer_matrix_weights_K_deriv++;
								}

								*pointer_gradient++ = (sum_y_ker_deriv_scalar - *pointer_m++ * sum_ker_deriv_scalar)/(*pointer_matrix_bandwidth++ * sum_ker);

							}
						}
					}

				}
				else
				{

					/* Adaptive */

					pointer_m = &mean[0];

					for(j=0; j < num_obs_eval; j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						pointer_yi = &vector_Y[0];
						pointer_matrix_weights_K = &matrix_weights_K[j][0];

						for (i =  0; i <  num_obs_train; i++)
						{
							sum_ker += *pointer_matrix_weights_K;
							sum_y_ker += *pointer_yi++ * *pointer_matrix_weights_K++;
						}

						*pointer_m++ = sum_y_ker/sum_ker;

					}

					/* Gradient using weight pointers - only compute necessary */

					for(k=0; k < num_var_test_int; k++)
					{
						tmp_k = var_index_int[k] - num_reg_unordered;
						if(var_index_int[k] >= num_reg_unordered)
						{
							pointer_m = &mean[0];
							pointer_gradient = &gradient[tmp_k][0];
							for(j=0; j < num_obs_eval; j++)
							{

								sum_ker_deriv_scalar = sum_y_ker_deriv_scalar = 0.0;
								sum_ker = DBL_MIN;

								pointer_yi = &vector_Y[0];
								pointer_matrix_weights_K = &matrix_weights_K[j][0];
								pointer_matrix_weights_K_deriv = &matrix_weights_K_deriv[tmp_k][j][0];
								pointer_matrix_bandwidth = &matrix_bandwidth_deriv[tmp_k][0];

								for (i =  0; i <  num_obs_train; i++)
								{
									sum_ker += *pointer_matrix_weights_K++ * *pointer_matrix_bandwidth++;
									sum_ker_deriv_scalar += *pointer_matrix_weights_K_deriv;
									sum_y_ker_deriv_scalar += *pointer_yi++ * *pointer_matrix_weights_K_deriv++;
								}

								/* Note - sum_ker is already divided by bw at i... */

								*pointer_gradient++ =(sum_y_ker_deriv_scalar - *pointer_m++ * sum_ker_deriv_scalar)/sum_ker;

							}
						}
					}

				}

			}
			else
			{

				/* Local linear using weights (also need data) */

				temp_mean_y = meand(num_obs_train, vector_Y);

				XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
				DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

				/* Conduct the estimation */

				for (j = 0; j < num_obs_eval; j++)
				{
					/* Initialize values to zero for a given evaluation point */
					for(k=0; k <= num_reg_cat_cont; k++)
					{
						XTKY[k][0] = 0.0;
						for(l=0; l <= num_reg_cat_cont; l++)
						{
							XTKX[k][l] = 0.0;
						}
					}

					pointer_yi = &vector_Y[0];
					pointer_matrix_weights_K = &matrix_weights_K[j][0];

					for(i=0; i < num_obs_train; i++)
					{

						/* Matrix is rows/cols... not C convention... */

						/* Upper left block */

						XTKX[0][0] += *pointer_matrix_weights_K;

						/* First element of XTKY */

						XTKY[0][0] += (temp = (*pointer_yi * *pointer_matrix_weights_K));

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* First lower column of XTKX */

							if(k < num_reg_continuous)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
									* *pointer_matrix_weights_K;
							}
							else if(k < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
									* *pointer_matrix_weights_K;
							}
							else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
									* *pointer_matrix_weights_K;
							}

							/* Diagonal of lower block of XTKX */

							XTKX[k+1][k+1] += ipow(temp1, 2) * *pointer_matrix_weights_K;

							/* Remaining elements of XTKY */

							XTKY[k+1][0] += temp1 * temp;

							/* Take advantage of symmetric nature of XTKX */

							for(l=0; l < k; l++)
							{
								if(l < num_reg_continuous)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
										* *pointer_matrix_weights_K;
								}
								else if(l < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
										* *pointer_matrix_weights_K;
								}
								else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
										* *pointer_matrix_weights_K;
								}
							}

						}

						pointer_yi++;
						pointer_matrix_weights_K++;

					}

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* Take advantage of symmetric nature */

						XTKX[0][k+1] = XTKX[k+1][0];

						for(l=0; l < k; l++)
						{
							XTKX[l+1][k+1] = XTKX[k+1][l+1];
						}

					}

					/* Now compute the beast... */

					if(fabs(mat_det(XTKX)) > 0.0 )
					{

						XTKXINV = mat_inv( XTKX, XTKXINV );

					}
					else
					{

						if(int_DEBUG == 1)
						{
							printf("\r                                                                                        ");
							printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
							printf("\n");
							mat_dumpf( XTKX, "%g ");
						}

						/* Add ridge factor - epsilon goes from zero to one/n*/

						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
						}

						/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

						do
						{
							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
								nepsilon += epsilon;
							}
						} while (fabs(mat_det(XTKX)) == 0.0);

						XTKXINV = mat_inv( XTKX, XTKXINV );
						/* Add epsilon times local constant estimator to first element of XTKY */
            XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

					}

					DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

					mean[j] = DELTA[0][0];

					for(k = 0; k < num_reg_cat_cont; k++)
					{
						gradient[k][j] = - DELTA[k+1][0];
					}

				}

				mat_free( XTKX );
				mat_free( XTKXINV );
				mat_free( XTKY );
				mat_free( DELTA );

			}

		}

	}
	else
	{

		/* No gradient computed */

		if(int_WEIGHTS == 0)
		{

			/* Do not use weights */

			if(int_ll == 0)
			{

				/* Conduct Nadaraya-Watson estimation */

				if(BANDWIDTH_reg == 0)
				{

					pointer_m = &mean[0];

					for(j=0; j < num_obs_eval; j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;
						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi++ * prod_kernel;

						}

						*pointer_m++ = sum_y_ker/sum_ker;

					}

				}
				else if(BANDWIDTH_reg == 1)
				{

					pointer_m = &mean[0];

					for(j=0; j < num_obs_eval; j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;
						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi++ * prod_kernel;

						}

						*pointer_m++ = sum_y_ker/sum_ker;

					}

				}
				else
				{

					pointer_m = &mean[0];

					for(j=0; j < num_obs_eval; j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;
						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi++ * prod_kernel;

						}

						*pointer_m++ = sum_y_ker/sum_ker;

					}

				}

			}
			else
			{

				/* Local linear */

				temp_mean_y = meand(num_obs_train, vector_Y);

				XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
				DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

				/* Conduct the estimation */

				if(BANDWIDTH_reg == 0)
				{

					for (j = 0; j < num_obs_eval; j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;
							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][0]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(k = 0; k < num_reg_ordered; k++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */

								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if(int_DEBUG == 1)
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}

							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
							XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

						mean[j] = DELTA[0][0];

					}

				}
				else if(BANDWIDTH_reg == 1)
				{

					for (j = 0; j < num_obs_eval; j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;
							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][j]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(k = 0; k < num_reg_ordered; k++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */

								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if(int_DEBUG == 1)
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}

							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
							XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

						mean[j] = DELTA[0][0];

					}

				}
				else
				{

					for (j = 0; j < num_obs_eval; j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;

							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][i]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(k = 0; k < num_reg_ordered; k++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */
								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if(int_DEBUG == 1)
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}

							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
							XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

						mean[j] = DELTA[0][0];

					}

				}

				mat_free( XTKX );
				mat_free( XTKXINV );
				mat_free( XTKY );
				mat_free( DELTA );

			}

		}
		else
		{

			/* Use weights, mean only */

			if(int_ll == 0)
			{

				/* Same for NW, ADA, GEN */

				pointer_m = &mean[0];

				for(j=0; j < num_obs_eval; j++)
				{

					sum_y_ker = 0.0;
					sum_ker = DBL_MIN;

					pointer_yi = &vector_Y[0];
					pointer_matrix_weights_K = &matrix_weights_K[j][0];

					for (i =  0; i <  num_obs_train; i++)
					{
						sum_ker += *pointer_matrix_weights_K;
						sum_y_ker += *pointer_yi++ * *pointer_matrix_weights_K++;
					}

					*pointer_m++ = sum_y_ker/sum_ker;

				}

			}
			else
			{

				/* Local linear using weights (also need data) */

				temp_mean_y = meand(num_obs_train, vector_Y);

				XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
				DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

				/* Conduct the estimation */

				for (j = 0; j < num_obs_eval; j++)
				{
					/* Initialize values to zero for a given evaluation point */
					for(k=0; k <= num_reg_cat_cont; k++)
					{
						XTKY[k][0] = 0.0;
						for(l=0; l <= num_reg_cat_cont; l++)
						{
							XTKX[k][l] = 0.0;
						}
					}

					pointer_yi = &vector_Y[0];
					pointer_matrix_weights_K = &matrix_weights_K[j][0];

					for(i=0; i < num_obs_train; i++)
					{

						/* Matrix is rows/cols... not C convention... */

						/* Upper left block */

						XTKX[0][0] += *pointer_matrix_weights_K;

						/* First element of XTKY */

						XTKY[0][0] += (temp = (*pointer_yi * *pointer_matrix_weights_K));

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* First lower column of XTKX */
							if(k < num_reg_continuous)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
									* *pointer_matrix_weights_K;
							}
							else if(k < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
									* *pointer_matrix_weights_K;
							}
							else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
									* *pointer_matrix_weights_K;
							}

							/* Diagonal of lower block of XTKX */

							XTKX[k+1][k+1] += ipow(temp1, 2) * *pointer_matrix_weights_K;

							/* Remaining elements of XTKY */

							XTKY[k+1][0] += temp1 * temp;

							/* Take advantage of symmetric nature of XTKX */

							for(l=0; l < k; l++)
							{
								if(l < num_reg_continuous)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
										* *pointer_matrix_weights_K;
								}
								else if(l < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
										* *pointer_matrix_weights_K;
								}
								else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
										* *pointer_matrix_weights_K;
								}

							}

						}

						pointer_yi++;
						pointer_matrix_weights_K++;

					}

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* Take advantage of symmetric nature */

						XTKX[0][k+1] = XTKX[k+1][0];

						for(l=0; l < k; l++)
						{

							XTKX[l+1][k+1] = XTKX[k+1][l+1];

						}

					}

					/* Now compute the beast... */

					if(fabs(mat_det(XTKX)) > 0.0 )
					{

						XTKXINV = mat_inv( XTKX, XTKXINV );

					}
					else
					{

						if(int_DEBUG == 1)
						{
							printf("\r                                                                                        ");
							printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
							printf("\n");
							mat_dumpf( XTKX, "%g ");
						}

						/* Add ridge factor - epsilon goes from zero to one/n*/

						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
						}

						/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

						do
						{
							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
								nepsilon += epsilon;
							}
						} while (fabs(mat_det(XTKX)) == 0.0);

						XTKXINV = mat_inv( XTKX, XTKXINV );
						/* Add epsilon times local constant estimator to first element of XTKY */
						XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

					}

					DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

					mean[j] = DELTA[0][0];

				}

				mat_free( XTKX );
				mat_free( XTKXINV );
				mat_free( XTKY );
				mat_free( DELTA );

			}

		}

	}
	#endif

	#ifdef MPI2

	if(int_compute_gradient == 1)
	{

		/* Gradient computed */

		if(int_WEIGHTS == 0)
		{

			/* Do not use weights */

			if(int_ll == 0)
			{

				prod_kernel_deriv = alloc_vecd(num_reg_continuous);

				sum_ker_deriv = alloc_vecd(num_reg_continuous);
				sum_y_ker_deriv = alloc_vecd(num_reg_continuous);

				/* Conduct Nadaraya-Watson estimation */

				if(BANDWIDTH_reg == 0)
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							*p_sum_y_ker_deriv++ = *p_sum_ker_deriv++ = 0.0;
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi * prod_kernel;

							/* Kernels for derivatives */

							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for (k = 0; k < num_reg_continuous; k++)
							{

								*p_prod_kernel_deriv = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][0]);

								for (l = 0; l < num_reg_continuous; l++)
								{

									if(l != k)
									{
										*p_prod_kernel_deriv *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
									}

								}

								p_prod_kernel_deriv++;

							}

							/* Gradient */

							p_sum_ker_deriv = &sum_ker_deriv[0];
							p_sum_y_ker_deriv = &sum_y_ker_deriv[0];
							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for(l = 0; l < num_reg_continuous; l++)
							{
								*p_sum_ker_deriv++ += *p_prod_kernel_deriv;
								*p_sum_y_ker_deriv++ += *pointer_yi * *p_prod_kernel_deriv++;
							}

							pointer_yi++;

						}

						mean[j-my_rank*stride] = sum_y_ker/sum_ker;

						/* gradient[0][] is that for _first_ continuous variable */

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							gradient[l][j-my_rank*stride] = (*p_sum_y_ker_deriv++ - mean[j-my_rank*stride] * *p_sum_ker_deriv++)/(matrix_bandwidth_deriv[l][0]*sum_ker);
						}

					}

				}
				else if(BANDWIDTH_reg == 1)
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							*p_sum_y_ker_deriv++ = *p_sum_ker_deriv++ = 0.0;
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi * prod_kernel;

							/* Kernels for derivatives */

							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for (k = 0; k < num_reg_continuous; k++)
							{

								*p_prod_kernel_deriv = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][j]);

								for (l = 0; l < num_reg_continuous; l++)
								{

									if(l != k)
									{
										*p_prod_kernel_deriv *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
									}

								}

								p_prod_kernel_deriv++;

							}

							/* Gradient */

							p_sum_ker_deriv = &sum_ker_deriv[0];
							p_sum_y_ker_deriv = &sum_y_ker_deriv[0];
							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for(l = 0; l < num_reg_continuous; l++)
							{
								*p_sum_ker_deriv++ += *p_prod_kernel_deriv;
								*p_sum_y_ker_deriv++ += *pointer_yi * *p_prod_kernel_deriv++;
							}

							pointer_yi++;

						}

						mean[j-my_rank*stride] = sum_y_ker/sum_ker;

						/* gradient[0][] is that for _first_ continuous variable */

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							gradient[l][j-my_rank*stride] = (*p_sum_y_ker_deriv++ - mean[j-my_rank*stride] * *p_sum_ker_deriv++)/(matrix_bandwidth_deriv[l][j]*sum_ker);
						}

					}

				}
				else
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							*p_sum_y_ker_deriv++ = *p_sum_ker_deriv++ = 0.0;
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Kernels for derivatives */

							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for (k = 0; k < num_reg_continuous; k++)
							{

								*p_prod_kernel_deriv = prod_kernel_cat * kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_deriv[k][i])/ipow(matrix_bandwidth_deriv[k][i],2);

								for (l = 0; l < num_reg_continuous; l++)
								{

									if(l != k)
									{
										*p_prod_kernel_deriv *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
									}

									p_prod_kernel_deriv++;

								}

							}

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi++*prod_kernel;

							/* Gradient */

							p_sum_ker_deriv = &sum_ker_deriv[0];
							p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

							p_prod_kernel_deriv = &prod_kernel_deriv[0];

							for(l = 0; l < num_reg_continuous; l++)
							{
								*p_sum_ker_deriv++ += *p_prod_kernel_deriv;
								*p_sum_y_ker_deriv++ += *pointer_yi * *p_prod_kernel_deriv++;
							}

							pointer_yi++;

						}

						mean[j-my_rank*stride] = sum_y_ker/sum_ker;

						/* gradient[0][] is that for _first_ continuous variable */

						p_sum_ker_deriv = &sum_ker_deriv[0];
						p_sum_y_ker_deriv = &sum_y_ker_deriv[0];

						for(l = 0; l < num_reg_continuous; l++)
						{
							gradient[l][j-my_rank*stride] = (*p_sum_y_ker_deriv++ - mean[j-my_rank*stride] * *p_sum_ker_deriv++)/sum_ker;
						}

					}

				}

				free(prod_kernel_deriv);
				free(sum_ker_deriv);
				free(sum_y_ker_deriv);

			}
			else
			{

				/* Local linear */

				temp_mean_y = meand(num_obs_train, vector_Y);

				XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
				DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

				/* Conduct the estimation */

				if(BANDWIDTH_reg == 0)
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;
							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][0]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */

								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if((int_DEBUG == 1)&&(my_rank == 0))
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}

							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
							XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
						mean[j-my_rank*stride] = DELTA[0][0];

						for(k = 0; k < num_reg_cat_cont; k++)
						{
							gradient[k][j] = - DELTA[k+1][0];
						}

					}

				}
				else if(BANDWIDTH_reg == 1)
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;
							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][j]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */

								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if((int_DEBUG == 1)&&(my_rank == 0))
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}

							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
							XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
						mean[j-my_rank*stride] = DELTA[0][0];

						for(k = 0; k < num_reg_cat_cont; k++)
						{
							gradient[k][j] = - DELTA[k+1][0];
						}

					}

				}
				else
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;

							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][i]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(k = 0; k < num_reg_ordered; k++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */

								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if((int_DEBUG == 1)&&(my_rank == 0))
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}

							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
							XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

						mean[j-my_rank*stride] = DELTA[0][0];

						for(k = 0; k < num_reg_cat_cont; k++)
						{
							gradient[k][j] = - DELTA[k+1][0];
						}

					}

				}

				mat_free( XTKX );
				mat_free( XTKXINV );
				mat_free( XTKY );
				mat_free( DELTA );

			}

		}
		else
		{

			/* Use weights, compute gradients */

			if(int_ll == 0)
			{

				/* Nadaraya Watson */

				if(BANDWIDTH_reg == 0)
				{

					/* Mean using weight pointers */

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						pointer_yi = &vector_Y[0];
						pointer_matrix_weights_K = &matrix_weights_K[j][0];

						for (i =  0; i <  num_obs_train; i++)
						{

							sum_ker += *pointer_matrix_weights_K;
							sum_y_ker += *pointer_yi++ * *pointer_matrix_weights_K++;

						}

						mean[j-my_rank*stride] = sum_y_ker/sum_ker;

					}

					/* Gradient using weight pointers - only compute necessary */

					for(k=0; k < num_var_test_int; k++)
					{
						/* Continuous variable */
						tmp_k = var_index_int[k] - num_reg_unordered;
						if(var_index_int[k] >= num_reg_unordered)
						{
							pointer_gradient = &gradient[tmp_k][0];
							pointer_matrix_bandwidth = &matrix_bandwidth_deriv[tmp_k][0];

							for(j=0; j < num_obs_eval; j++)
							{

								sum_ker_deriv_scalar = sum_y_ker_deriv_scalar = 0.0;
								sum_ker = DBL_MIN;

								pointer_yi = &vector_Y[0];
								pointer_matrix_weights_K = &matrix_weights_K[j][0];
								pointer_matrix_weights_K_deriv = &matrix_weights_K_deriv[tmp_k][j][0];

								for (i =  0; i <  num_obs_train; i++)
								{
									sum_ker += *pointer_matrix_weights_K++;
									sum_ker_deriv_scalar += *pointer_matrix_weights_K_deriv;
									sum_y_ker_deriv_scalar += *pointer_yi++ * *pointer_matrix_weights_K_deriv++;
								}

								*pointer_gradient++ = (sum_y_ker_deriv_scalar - mean[j-my_rank*stride] * sum_ker_deriv_scalar)/(*pointer_matrix_bandwidth * sum_ker);

							}
						}
					}

				}
				else if(BANDWIDTH_reg == 1)
				{

					/* Mean using weight pointers */

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						pointer_yi = &vector_Y[0];
						pointer_matrix_weights_K = &matrix_weights_K[j][0];

						for (i =  0; i <  num_obs_train; i++)
						{

							sum_ker += *pointer_matrix_weights_K;
							sum_y_ker += *pointer_yi++ * *pointer_matrix_weights_K++;

						}

						mean[j-my_rank*stride] = sum_y_ker/sum_ker;

					}

					/* Gradient using weight pointers - only compute necessary */

					for(k=0; k < num_var_test_int; k++)
					{
						/* Continuous variable */
						tmp_k = var_index_int[k] - num_reg_unordered;
						if(var_index_int[k] >= num_reg_unordered)
						{
							pointer_gradient = &gradient[tmp_k][0];
							pointer_matrix_bandwidth = &matrix_bandwidth_deriv[tmp_k][0];

							for(j=0; j < num_obs_eval; j++)
							{

								sum_ker_deriv_scalar = sum_y_ker_deriv_scalar = 0.0;
								sum_ker = DBL_MIN;

								pointer_yi = &vector_Y[0];
								pointer_matrix_weights_K = &matrix_weights_K[j][0];
								pointer_matrix_weights_K_deriv = &matrix_weights_K_deriv[tmp_k][j][0];

								for (i =  0; i <  num_obs_train; i++)
								{
									sum_ker += *pointer_matrix_weights_K++;
									sum_ker_deriv_scalar += *pointer_matrix_weights_K_deriv;
									sum_y_ker_deriv_scalar += *pointer_yi++ * *pointer_matrix_weights_K_deriv++;
								}

								*pointer_gradient++ = (sum_y_ker_deriv_scalar - mean[j-my_rank*stride] * sum_ker_deriv_scalar)/(*pointer_matrix_bandwidth++ * sum_ker);

							}
						}
					}

				}
				else
				{

					/* Adaptive */

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;

						pointer_yi = &vector_Y[0];
						pointer_matrix_weights_K = &matrix_weights_K[j][0];

						for (i =  0; i <  num_obs_train; i++)
						{
							sum_ker += *pointer_matrix_weights_K;
							sum_y_ker += *pointer_yi++ * *pointer_matrix_weights_K++;
						}

						mean[j-my_rank*stride] = sum_y_ker/sum_ker;

					}

					/* Gradient using weight pointers - only compute necessary */

					for(k=0; k < num_var_test_int; k++)
					{
						tmp_k = var_index_int[k] - num_reg_unordered;
						if(var_index_int[k] >= num_reg_unordered)
						{
							pointer_gradient = &gradient[tmp_k][0];
							for(j=0; j < num_obs_eval; j++)
							{

								sum_ker_deriv_scalar = sum_y_ker_deriv_scalar = 0.0;
								sum_ker = DBL_MIN;

								pointer_yi = &vector_Y[0];
								pointer_matrix_weights_K = &matrix_weights_K[j][0];
								pointer_matrix_weights_K_deriv = &matrix_weights_K_deriv[tmp_k][j][0];
								pointer_matrix_bandwidth = &matrix_bandwidth_deriv[tmp_k][0];

								for (i =  0; i <  num_obs_train; i++)
								{
									sum_ker += *pointer_matrix_weights_K++ * *pointer_matrix_bandwidth++;
									sum_ker_deriv_scalar += *pointer_matrix_weights_K_deriv;
									sum_y_ker_deriv_scalar += *pointer_yi++ * *pointer_matrix_weights_K_deriv++;
								}

								/* Note - sum_ker is already divided by bw at i... */

								*pointer_gradient++ =(sum_y_ker_deriv_scalar - mean[j-my_rank*stride] * sum_ker_deriv_scalar)/sum_ker;

							}
						}
					}

				}

			}
			else
			{

				/* Local linear using weights (also need data) */

				temp_mean_y = meand(num_obs_train, vector_Y);

				XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
				DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

				/* Conduct the estimation */

				for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
				{
					/* Initialize values to zero for a given evaluation point */
					for(k=0; k <= num_reg_cat_cont; k++)
					{
						XTKY[k][0] = 0.0;
						for(l=0; l <= num_reg_cat_cont; l++)
						{
							XTKX[k][l] = 0.0;
						}
					}

					pointer_yi = &vector_Y[0];
					pointer_matrix_weights_K = &matrix_weights_K[j][0];

					for(i=0; i < num_obs_train; i++)
					{

						/* Matrix is rows/cols... not C convention... */

						/* Upper left block */

						XTKX[0][0] += *pointer_matrix_weights_K;

						/* First element of XTKY */

						XTKY[0][0] += (temp = (*pointer_yi * *pointer_matrix_weights_K));

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* First lower column of XTKX */

							if(k < num_reg_continuous)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
									* *pointer_matrix_weights_K;
							}
							else if(k < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
									* *pointer_matrix_weights_K;
							}
							else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
									* *pointer_matrix_weights_K;
							}

							/* Diagonal of lower block of XTKX */

							XTKX[k+1][k+1] += ipow(temp1, 2) * *pointer_matrix_weights_K;

							/* Remaining elements of XTKY */

							XTKY[k+1][0] += temp1 * temp;

							/* Take advantage of symmetric nature of XTKX */

							for(l=0; l < k; l++)
							{
								if(l < num_reg_continuous)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
										* *pointer_matrix_weights_K;
								}
								else if(l < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
										* *pointer_matrix_weights_K;
								}
								else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
										* *pointer_matrix_weights_K;
								}
							}

						}

						pointer_yi++;
						pointer_matrix_weights_K++;

					}

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* Take advantage of symmetric nature */

						XTKX[0][k+1] = XTKX[k+1][0];

						for(l=0; l < k; l++)
						{
							XTKX[l+1][k+1] = XTKX[k+1][l+1];
						}

					}

					/* Now compute the beast... */

					if(fabs(mat_det(XTKX)) > 0.0 )
					{

						XTKXINV = mat_inv( XTKX, XTKXINV );

					}
					else
					{

						if((int_DEBUG == 1)&&(my_rank == 0))
						{
							printf("\r                                                                                        ");
							printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
							printf("\n");
							mat_dumpf( XTKX, "%g ");
						}

						/* Add ridge factor - epsilon goes from zero to one/n*/

						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
						}

						/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

						do
						{
							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
								nepsilon += epsilon;
							}
						} while (fabs(mat_det(XTKX)) == 0.0);

						XTKXINV = mat_inv( XTKX, XTKXINV );
						/* Add epsilon times local constant estimator to first element of XTKY */
						XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

					}

					DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

					mean[j-my_rank*stride] = DELTA[0][0];

					for(k = 0; k < num_reg_cat_cont; k++)
					{
						gradient[k][j] = - DELTA[k+1][0];
					}

				}

				mat_free( XTKX );
				mat_free( XTKXINV );
				mat_free( XTKY );
				mat_free( DELTA );

			}

		}

	}
	else
	{

		/* No gradient computed */

		if(int_WEIGHTS == 0)
		{

			/* Do not use weights */

			if(int_ll == 0)
			{

				/* Conduct Nadaraya-Watson estimation */

				if(BANDWIDTH_reg == 0)
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;
						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi++ * prod_kernel;

						}

						mean[j-my_rank*stride] = sum_y_ker/sum_ker;

					}

				}
				else if(BANDWIDTH_reg == 1)
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;
						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi++ * prod_kernel;

						}

						mean[j-my_rank*stride] = sum_y_ker/sum_ker;

					}

				}
				else
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{

						sum_y_ker = 0.0;
						sum_ker = DBL_MIN;
						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Kernel for mean */

							prod_kernel_cont = 1.0;

							for(l = 0; l < num_reg_continuous; l++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
							}

							prod_kernel_cat = 1.0;

							for(l = 0; l < num_reg_unordered; l++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
							}

							for(l = 0; l < num_reg_ordered; l++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Mean */

							sum_ker += prod_kernel;
							sum_y_ker += *pointer_yi++ * prod_kernel;

						}

						mean[j-my_rank*stride] = sum_y_ker/sum_ker;

					}

				}

			}
			else
			{

				/* Local linear */

				temp_mean_y = meand(num_obs_train, vector_Y);

				XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
				DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

				/* Conduct the estimation */

				if(BANDWIDTH_reg == 0)
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;
							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][0]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(k = 0; k < num_reg_ordered; k++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */

								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if((int_DEBUG == 1)&&(my_rank == 0))
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}

							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
							XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

						mean[j-my_rank*stride] = DELTA[0][0];

					}

				}
				else if(BANDWIDTH_reg == 1)
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;
							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][j]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(k = 0; k < num_reg_ordered; k++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */

								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if((int_DEBUG == 1)&&(my_rank == 0))
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}

							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
							XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

						mean[j-my_rank*stride] = DELTA[0][0];

					}

				}
				else
				{

					for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
					{
						/* Initialize values to zero for a given evaluation point */
						for(k=0; k <= num_reg_cat_cont; k++)
						{
							XTKY[k][0] = 0.0;

							for(l=0; l <= num_reg_cat_cont; l++)
							{
								XTKX[k][l] = 0.0;
							}
						}

						pointer_yi = &vector_Y[0];

						for(i=0; i < num_obs_train; i++)
						{

							/* Matrix is rows/cols... not C convention... */

							prod_kernel_cont = 1.0;

							for(k = 0; k < num_reg_continuous; k++)
							{
								prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth[k][i]);
							}

							prod_kernel_cat = 1.0;

							for(k = 0; k < num_reg_unordered; k++)
							{
								prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[k][j],matrix_X_unordered_train[k][i],lambda[k],num_categories[k]);
							}

							for(k = 0; k < num_reg_ordered; k++)
							{
								prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[k][j],matrix_X_ordered_train[k][i],lambda[k+num_reg_unordered]);
							}

							prod_kernel = prod_kernel_cont*prod_kernel_cat;

							/* Upper left block */

							XTKX[0][0] += prod_kernel;

							/* First element of XTKY */

							XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

							for(k=0; k < num_reg_cat_cont; k++)
							{

								/* First lower column of XTKX */
								if(k < num_reg_continuous)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
										* prod_kernel;
								}
								else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
										* prod_kernel;
								}

								/* Diagonal of lower block of XTKX */

								XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

								/* Remaining elements of XTKY */

								XTKY[k+1][0] += temp1 * temp;

								/* Take advantage of symmetric nature of XTKX */

								for(l=0; l < k; l++)
								{
									if(l < num_reg_continuous)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
											* prod_kernel;
									}
									else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
									{
										XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
											* prod_kernel;
									}
								}

							}

							pointer_yi++;

						}

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* Take advantage of symmetric nature */

							XTKX[0][k+1] = XTKX[k+1][0];

							for(l=0; l < k; l++)
							{
								XTKX[l+1][k+1] = XTKX[k+1][l+1];
							}

						}

						/* Now compute the beast... */

						if(fabs(mat_det(XTKX)) > 0.0 )
						{

							XTKXINV = mat_inv( XTKX, XTKXINV );

						}
						else
						{

							if((int_DEBUG == 1)&&(my_rank == 0))
							{
								printf("\r                                                                                        ");
								printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
								printf("\n");
								mat_dumpf( XTKX, "%g ");
							}

							/* Add ridge factor - epsilon goes from zero to one/n*/

							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
							}

							/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

							do
							{
								for(k=0; k < num_reg_cat_cont + 1; k++)
								{
									XTKX[k][k] += epsilon;
									nepsilon += epsilon;
								}
							} while (fabs(mat_det(XTKX)) == 0.0);

							XTKXINV = mat_inv( XTKX, XTKXINV );
							/* Add epsilon times local constant estimator to first element of XTKY */
							XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

						}

						DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

						mean[j-my_rank*stride] = DELTA[0][0];

					}

				}

				mat_free( XTKX );
				mat_free( XTKXINV );
				mat_free( XTKY );
				mat_free( DELTA );

			}

		}
		else
		{

			/* Use weights, mean only */

			if(int_ll == 0)
			{

				/* Same for NW, ADA, GEN */

				for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
				{

					sum_y_ker = 0.0;
					sum_ker = DBL_MIN;

					pointer_yi = &vector_Y[0];
					pointer_matrix_weights_K = &matrix_weights_K[j][0];

					for (i =  0; i <  num_obs_train; i++)
					{
						sum_ker += *pointer_matrix_weights_K;
						sum_y_ker += *pointer_yi++ * *pointer_matrix_weights_K++;
					}

					mean[j-my_rank*stride] = sum_y_ker/sum_ker;

				}

			}
			else
			{

				/* Local linear using weights (also need data) */

				temp_mean_y = meand(num_obs_train, vector_Y);

				XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
				XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
				DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

				/* Conduct the estimation */

				for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
				{
					/* Initialize values to zero for a given evaluation point */
					for(k=0; k <= num_reg_cat_cont; k++)
					{
						XTKY[k][0] = 0.0;
						for(l=0; l <= num_reg_cat_cont; l++)
						{
							XTKX[k][l] = 0.0;
						}
					}

					pointer_yi = &vector_Y[0];
					pointer_matrix_weights_K = &matrix_weights_K[j][0];

					for(i=0; i < num_obs_train; i++)
					{

						/* Matrix is rows/cols... not C convention... */

						/* Upper left block */

						XTKX[0][0] += *pointer_matrix_weights_K;

						/* First element of XTKY */

						XTKY[0][0] += (temp = (*pointer_yi * *pointer_matrix_weights_K));

						for(k=0; k < num_reg_cat_cont; k++)
						{

							/* First lower column of XTKX */
							if(k < num_reg_continuous)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_continuous_eval[k][j] - matrix_X_continuous_train[k][i]))
									* *pointer_matrix_weights_K;
							}
							else if(k < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_unordered_eval[k-num_reg_continuous][j] - matrix_X_unordered_train[k-num_reg_continuous][i]))
									* *pointer_matrix_weights_K;
							}
							else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][0] += (temp1 = (matrix_X_ordered_eval[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[k-num_reg_continuous-num_reg_unordered][i]))
									* *pointer_matrix_weights_K;
							}

							/* Diagonal of lower block of XTKX */

							XTKX[k+1][k+1] += ipow(temp1, 2) * *pointer_matrix_weights_K;

							/* Remaining elements of XTKY */

							XTKY[k+1][0] += temp1 * temp;

							/* Take advantage of symmetric nature of XTKX */

							for(l=0; l < k; l++)
							{
								if(l < num_reg_continuous)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_continuous_eval[l][j] - matrix_X_continuous_train[l][i])
										* *pointer_matrix_weights_K;
								}
								else if(l < num_reg_continuous+num_reg_unordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_unordered_eval[l-num_reg_continuous][j] - matrix_X_unordered_train[l-num_reg_continuous][i])
										* *pointer_matrix_weights_K;
								}
								else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
								{
									XTKX[k+1][l+1] += temp1 * (matrix_X_ordered_eval[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered_train[l-num_reg_continuous-num_reg_unordered][i])
										* *pointer_matrix_weights_K;
								}

							}

						}

						pointer_yi++;
						pointer_matrix_weights_K++;

					}

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* Take advantage of symmetric nature */

						XTKX[0][k+1] = XTKX[k+1][0];

						for(l=0; l < k; l++)
						{

							XTKX[l+1][k+1] = XTKX[k+1][l+1];

						}

					}

					/* Now compute the beast... */

					if(fabs(mat_det(XTKX)) > 0.0 )
					{

						XTKXINV = mat_inv( XTKX, XTKXINV );

					}
					else
					{

						if((int_DEBUG == 1)&&(my_rank == 0))
						{
							printf("\r                                                                                        ");
							printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_no_stderr()",j);
							printf("\n");
							mat_dumpf( XTKX, "%g ");
						}

						/* Add ridge factor - epsilon goes from zero to one/n*/

						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
						}

						/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

						do
						{
							for(k=0; k < num_reg_cat_cont + 1; k++)
							{
								XTKX[k][k] += epsilon;
								nepsilon += epsilon;
							}
						} while (fabs(mat_det(XTKX)) == 0.0);

						XTKXINV = mat_inv( XTKX, XTKXINV );
						/* Add epsilon times local constant estimator to first element of XTKY */
						XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

					}

					DELTA =  mat_mul( XTKXINV, XTKY, DELTA);

					mean[j-my_rank*stride] = DELTA[0][0];

				}

				mat_free( XTKX );
				mat_free( XTKXINV );
				mat_free( XTKY );
				mat_free( DELTA );

			}

		}

	}

	/* Important - only one gather per module */

	MPI_Gather(mean, stride, MPI_DOUBLE, mean, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(mean, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(int_compute_gradient == 1)
	{

		for(l = 0; l < num_reg_continuous; l++)
		{

			MPI_Gather(&gradient[l][0], stride, MPI_DOUBLE, &gradient[l][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&gradient[l][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		}

	}
	#endif

	return(0);

}


int kernel_estimate_con_density_categorical(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int num_obs_train,
int num_obs_eval,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_Y_unordered_eval,
double **matrix_Y_ordered_eval,
double **matrix_Y_continuous_eval,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_scale_factor,
int *num_categories,
double *pdf,
double *pdf_stderr,
double *log_likelihood)
{

	/* This function estimates a density function using both continuous */
	/* and categorical covariates with three estimation techniques and an */
	/* assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel_cat;
	double prod_kernel_cont;

	double prod_kernel_marginal_cat;
	double prod_kernel_marginal_cont;

	double sum_ker;
	double sum_ker_marginal;
	double sum_ker_temp;

	double prod_h;
	double prod_h_marginal;

	double *lambda;
	double **matrix_bandwidth_var = NULL;
	double **matrix_bandwidth_reg = NULL;

	/* Integral of K(z)^p */
	double INT_KERNEL_P;
	/* Number of regressors times integral of K(z)^p */
	double K_INT_KERNEL_P;

	double log_DBL_MIN = log(DBL_MIN);

	/* Allocate memory for objects */

	#ifdef MPI2
	double log_likelihood_MPI;
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	lambda = alloc_vecd(num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered);

	if((BANDWIDTH_den == 0)||(BANDWIDTH_den == 1))
	{
		matrix_bandwidth_var = alloc_matd(num_obs_eval,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_eval,num_reg_continuous);
	}
	else if(BANDWIDTH_den == 2)
	{
		matrix_bandwidth_var = alloc_matd(num_obs_train,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_train,num_reg_continuous);
	}

	/* Bandwidths for `dependent' variables */

	if(kernel_bandwidth_mean(
		KERNEL_den,
		BANDWIDTH_den,
		num_obs_train,
		num_obs_eval,
		num_var_continuous,
		num_var_unordered,
		num_var_ordered,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
		vector_scale_factor,
		matrix_Y_continuous_train,
		matrix_Y_continuous_eval,
		matrix_X_continuous_train,
		matrix_X_continuous_eval,
		matrix_bandwidth_var,
		matrix_bandwidth_reg,
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

	/* Initialize constants for various kernels required for asymptotic standard errors */

	initialize_kernel_density_asymptotic_constants(
		KERNEL_den,
		num_reg_continuous,
		&INT_KERNEL_P,
		&K_INT_KERNEL_P);

	#ifndef MPI2

	/* Initialize log likelihood */

	*log_likelihood = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			prod_h = 1.0;

			for(l = 0; l < num_var_continuous; l++)
			{
				prod_h *= matrix_bandwidth_var[l][0];
			}

			sum_ker = 0.0;
			sum_ker_marginal = 0.0;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][0]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;

				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

			}

			pdf[j] = sum_ker/(prod_h*sum_ker_marginal);

			/* With no continuous variables, need to drop K_INT_KERNEL_P, prod_h ==1 */

			pdf_stderr[j] = ((num_reg_continuous != 0) ? sqrt(pdf[j]*K_INT_KERNEL_P/prod_h) : sqrt(pdf[j]/prod_h));

			if(pdf[j] > DBL_MIN)
			{
				*log_likelihood += log(pdf[j]);
			}
			else
			{
				*log_likelihood += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			prod_h = 1.0;

			for(l = 0; l < num_var_continuous; l++)
			{
				prod_h *= matrix_bandwidth_var[l][j];
			}

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][j]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;

				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

			}

			pdf[j] = sum_ker/(prod_h*sum_ker_marginal);

			/* With no continuous variables, need to drop K_INT_KERNEL_P, prod_h ==1 */

			pdf_stderr[j] = ((num_reg_continuous != 0) ? sqrt(pdf[j]*K_INT_KERNEL_P/prod_h) : sqrt(pdf[j]/prod_h));

			if(pdf[j] > DBL_MIN)
			{
				*log_likelihood += log(pdf[j]);
			}
			else
			{
				*log_likelihood += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

		}

	}
	else
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = sum_ker_temp = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_h = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_h *= matrix_bandwidth_var[l][i];
				}

				prod_h_marginal = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_h_marginal *= matrix_bandwidth_reg[l][i];
				}

				prod_h *= prod_h_marginal;

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i]);
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][i]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h_marginal;
				sum_ker_temp += prod_kernel_cont*prod_kernel_cat/ipow(prod_h, 2);

			}

			pdf[j] = sum_ker/sum_ker_marginal;

			/* stderr needs extra tmp I think */
			/* With no continuous variables, need to drop K_INT_KERNEL_P, prod_h ==1 */

			pdf_stderr[j] = ((num_reg_continuous != 0) ? sqrt(sum_ker_temp*K_INT_KERNEL_P/sum_ker_marginal) : sqrt(sum_ker_temp/sum_ker_marginal));

			if(pdf[j] > DBL_MIN)
			{
				*log_likelihood += log(pdf[j]);
			}
			else
			{
				*log_likelihood += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

		}

	}
	#endif

	#ifdef MPI2

	/* Initialize log likelihood */

	log_likelihood_MPI = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			prod_h = 1.0;

			for(l = 0; l < num_var_continuous; l++)
			{
				prod_h *= matrix_bandwidth_var[l][0];
			}

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][0]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;

				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

			}

			pdf[j-my_rank*stride] = sum_ker/(prod_h*sum_ker_marginal);

			/* With no continuous variables, need to drop K_INT_KERNEL_P, prod_h ==1 */

			pdf_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(pdf[j-my_rank*stride]*K_INT_KERNEL_P/prod_h) : sqrt(pdf[j-my_rank*stride]/prod_h));

			if(pdf[j-my_rank*stride] > DBL_MIN)
			{
				log_likelihood_MPI += log(pdf[j-my_rank*stride]);
			}
			else
			{
				log_likelihood_MPI += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			prod_h = 1.0;

			for(l = 0; l < num_var_continuous; l++)
			{
				prod_h *= matrix_bandwidth_var[l][j];
			}

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][j]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;

				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

			}

			pdf[j-my_rank*stride] = sum_ker/(prod_h*sum_ker_marginal);

			/* With no continuous variables, need to drop K_INT_KERNEL_P, prod_h ==1 */

			pdf_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(pdf[j-my_rank*stride]*K_INT_KERNEL_P/prod_h) : sqrt(pdf[j-my_rank*stride]/prod_h));

			if(pdf[j-my_rank*stride] > DBL_MIN)
			{
				log_likelihood_MPI += log(pdf[j-my_rank*stride]);
			}
			else
			{
				log_likelihood_MPI += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = sum_ker_temp = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_h = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_h *= matrix_bandwidth_var[l][i];
				}

				prod_h_marginal = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_h_marginal *= matrix_bandwidth_reg[l][i];
				}

				prod_h *= prod_h_marginal;

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i]);
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][i]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h_marginal;
				sum_ker_temp += prod_kernel_cont*prod_kernel_cat/ipow(prod_h, 2);

			}

			pdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;
			/* stderr needs extra tmp I think */
			/* With no continuous variables, need to drop K_INT_KERNEL_P, prod_h ==1 */

			pdf_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(sum_ker_temp*K_INT_KERNEL_P/sum_ker_marginal) : sqrt(sum_ker_temp/sum_ker_marginal));

			if(pdf[j-my_rank*stride] > DBL_MIN)
			{
				log_likelihood_MPI += log(pdf[j-my_rank*stride]);
			}
			else
			{
				log_likelihood_MPI += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

		}

	}

	MPI_Gather(pdf, stride, MPI_DOUBLE, pdf, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(pdf, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Gather(pdf_stderr, stride, MPI_DOUBLE, pdf_stderr, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(pdf_stderr, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Reduce(&log_likelihood_MPI, log_likelihood, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(log_likelihood, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	free(lambda);

	free_mat(matrix_bandwidth_var,num_var_continuous);
	free_mat(matrix_bandwidth_reg,num_reg_continuous);

	return(0);

}


int kernel_estimate_con_distribution_categorical(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int num_obs_train,
int num_obs_eval,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_Y_unordered_eval,
double **matrix_Y_ordered_eval,
double **matrix_Y_continuous_eval,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *cdf,
double *cdf_stderr,
double small,
int itmax)
{

	/* This function estimates a density function using both continuous */
	/* and categorical covariates with three estimation techniques and an */
	/* assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel_cat;
	double prod_kernel_cont;

	double prod_kernel_marginal_cat;
	double prod_kernel_marginal_cont;

	double sum_ker;
	double sum_ker_marginal;

	double prod_h;

	double *lambda;
	double **matrix_bandwidth_var = NULL;
	double **matrix_bandwidth_reg = NULL;

	#ifdef MPI2
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered);

	if((BANDWIDTH_den == 0)||(BANDWIDTH_den == 1))
	{
		matrix_bandwidth_var = alloc_matd(num_obs_eval,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_eval,num_reg_continuous);
	}
	else if(BANDWIDTH_den == 2)
	{
		matrix_bandwidth_var = alloc_matd(num_obs_train,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_train,num_reg_continuous);
	}

	/* Bandwidths for `dependent' variables */

	if(kernel_bandwidth_mean(
		KERNEL_den,
		BANDWIDTH_den,
		num_obs_train,
		num_obs_eval,
		num_var_continuous,
		num_var_unordered,
		num_var_ordered,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
		vector_scale_factor,
		matrix_Y_continuous_train,
		matrix_Y_continuous_eval,
		matrix_X_continuous_train,
		matrix_X_continuous_eval,
		matrix_bandwidth_var,
		matrix_bandwidth_reg,
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

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{

					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);

				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][0]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

			}

			cdf[j] = sum_ker/sum_ker_marginal;
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{

					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);

				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][j]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

			}

			cdf[j] = sum_ker/sum_ker_marginal;
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}
	else
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_h = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_h *= matrix_bandwidth_reg[l][i];
				}

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i]);
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][i]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h;

			}

			cdf[j] = sum_ker/sum_ker_marginal;
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}
	#endif

	#ifdef MPI2

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{

					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);

				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][0]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

			}

			cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;
			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{

					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);

				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][j]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

			}

			cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;
			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_h = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_h *= matrix_bandwidth_reg[l][i];
				}

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i]);
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][i]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h;

			}

			cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;
			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

		}

	}

	MPI_Gather(cdf, stride, MPI_DOUBLE, cdf, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(cdf, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(cdf_stderr, stride, MPI_DOUBLE, cdf_stderr, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(cdf_stderr, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	free(lambda);

	free_mat(matrix_bandwidth_var,num_var_continuous);
	free_mat(matrix_bandwidth_reg,num_reg_continuous);

	return(0);

}

/* 2/5/2010 - happy birthday! */

int kernel_estimate_con_distribution_categorical_leave_one_out(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int num_obs_train,
int num_obs_eval,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_Y_unordered_eval,
double **matrix_Y_ordered_eval,
double **matrix_Y_continuous_eval,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *cdf,
double small,
int itmax)
{

	/* This function estimates a leave-one-out distribution function
     using both continuous and categorical covariates with three
     estimation techniques and an assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel_cat;
	double prod_kernel_cont;

	double prod_kernel_marginal_cat;
	double prod_kernel_marginal_cont;

	double sum_ker;
	double sum_ker_marginal;

	double prod_h;

	double *lambda;
	double **matrix_bandwidth_var = NULL;
	double **matrix_bandwidth_reg = NULL;

	#ifdef MPI2
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered);

	if((BANDWIDTH_den == 0)||(BANDWIDTH_den == 1))
	{
		matrix_bandwidth_var = alloc_matd(num_obs_eval,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_eval,num_reg_continuous);
	}
	else if(BANDWIDTH_den == 2)
	{
		matrix_bandwidth_var = alloc_matd(num_obs_train,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_train,num_reg_continuous);
	}

	/* Bandwidths for `dependent' variables */

	if(kernel_bandwidth_mean(
    KERNEL_den,
		BANDWIDTH_den,
		num_obs_train,
		num_obs_eval,
		num_var_continuous,
		num_var_unordered,
		num_var_ordered,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
		vector_scale_factor,
		matrix_Y_continuous_train,
		matrix_Y_continuous_eval,
		matrix_X_continuous_train,
		matrix_X_continuous_eval,
		matrix_bandwidth_var,
		matrix_bandwidth_reg,
		lambda) == 1)
	{
    free(lambda);
    free_mat(matrix_bandwidth_var,num_var_continuous);
    free_mat(matrix_bandwidth_reg,num_reg_continuous);
		return(1);
	}

	#ifndef MPI2

  /* First stab could be brute force copy no saving */

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{

					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);

				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][0]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				if(i != j) {
          sum_ker += prod_kernel_cont*prod_kernel_cat;
          sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;
        }

			}

			cdf[j] = sum_ker/sum_ker_marginal;

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{

					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);

				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][j]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

			}

      if(i != j) {
        cdf[j] = sum_ker/sum_ker_marginal;
      }

		}

	}
	else
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_h = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_h *= matrix_bandwidth_reg[l][i];
				}

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i]);
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][i]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				if(i != j) {
          sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
          sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h;
        }

			}

			cdf[j] = sum_ker/sum_ker_marginal;

		}

	}
	#endif

	#ifdef MPI2

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{

					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);

				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][0]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				if(i != j) {
          sum_ker += prod_kernel_cont*prod_kernel_cat;
          sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;
        }

			}

			cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{

					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);

				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][j]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				if(i != j) {
          sum_ker += prod_kernel_cont*prod_kernel_cat;
          sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;
        }

			}

			cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_h = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_h *= matrix_bandwidth_reg[l][i];
				}

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i]);
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][i]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				if(i != j) {
          sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
          sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h;
        }

			}

			cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;

		}

	}

	MPI_Gather(cdf, stride, MPI_DOUBLE, cdf, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(cdf, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	#endif

	free(lambda);

	free_mat(matrix_bandwidth_var,num_var_continuous);
	free_mat(matrix_bandwidth_reg,num_reg_continuous);

	return(0);

}

/* Feb 6 2010 - this will call
   kernel_estimate_con_distribution_categorical_leave_one_out and
   return the cv function */

int indfunc(double a) { return(a <= 0.0 ? 1 : 0); }

int kernel_estimate_con_distribution_categorical_leave_one_out_ccdf(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int num_obs_train,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *cv, /* Returned by function */
double small,
int itmax)
{

	/* This function constructs the objective function for conditional
     distribution bandwidth selection */

  /* NOTE - currently supports multivariate continuous variables only */
  /* NOTE - sloowwww */

  /* For fully general implementation, need to feed in unique values
     for all ordered and unordered variables */

	/* Declarations and initializations  */

	int i;
  int l;
  int j;
  double indicator;
  double *cdf_loo;
	*cv = 0.0;

  /* Must declare and free evaluation */

	double **matrix_Y_unordered_eval;
	double **matrix_Y_ordered_eval;
	double **matrix_Y_continuous_eval;

  /* Must declare and free cdf_eval */

  cdf_loo = alloc_vecd(num_obs_train);

	matrix_Y_unordered_eval = alloc_matd(num_obs_train, num_reg_unordered);
	matrix_Y_ordered_eval = alloc_matd(num_obs_train, num_reg_ordered);
	matrix_Y_continuous_eval = alloc_matd(num_obs_train, num_reg_continuous);

  /* XXX */

    /*    for(l=0; l < num_reg_unordered; l++)
      for(j=0; j < num_obs_train; j++)
        matrix_Y_unordered_eval[l][j] = matrix_Y_unordered_train[l][i];
    
    for(l=0; l < num_reg_ordered; l++)
      for(j=0; j < num_obs_train; j++)
      matrix_Y_ordered_eval[l][j] = matrix_Y_ordered_train[l][i];*/
    
  for(i=0; i < num_obs_train; i++) {

    /* Brute force `copy' Y eval could of course be improved via
       pointer */

    /* Require nested loops for ordered, unordered, and continuous -
       not implemented yet (Feb 7 2010) */

    for(l=0; l < num_reg_continuous; l++)
      for(j=0; j < num_obs_train; j++)
        matrix_Y_continuous_eval[l][j] = matrix_Y_continuous_train[l][i];

    if(kernel_estimate_con_distribution_categorical_leave_one_out(
 			KERNEL_den,
			KERNEL_unordered_den,
			KERNEL_ordered_den,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_den,
			num_obs_train,
			num_obs_train,
			num_var_unordered,
			num_var_ordered,
			num_var_continuous,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_Y_unordered_train,
			matrix_Y_ordered_train,
			matrix_Y_continuous_train,
			matrix_Y_unordered_eval,
			matrix_Y_ordered_eval,
			matrix_Y_continuous_eval,
			matrix_X_unordered_train,
			matrix_X_ordered_train,
			matrix_X_continuous_train,
			matrix_X_unordered_train,
			matrix_X_ordered_train,
			matrix_X_continuous_train,
			vector_scale_factor,
			num_categories,
			matrix_categorical_vals,
			cdf_loo,
			small,
			itmax) == 1) 
    {
      return(1);
    }

    /*  If y is a discrete scalar we would replace the integral with
       the sum over all unique realizations... in a multivariate
       real-valued setting a multivariate integral becomes a
       multivariate sum so we would need to loop over all
       dimensions... for mixed data ditto - sum is certainly easier */

    /* Must use multivariate product indicator function */

    for(j=0; j < num_obs_train; j++) {
      indicator = 1.0;
      for(l=0; l < num_reg_continuous; l++) {
        indicator *= indfunc(matrix_Y_continuous_train[l][j]-matrix_Y_continuous_eval[l][j]);
      }
      *cv += ipow(indicator - cdf_loo[j],2);
    }

  }

  /* Sum over all variables - need to be careful about denominator */

  *cv /= (double) ipow(num_obs_train,1+num_reg_continuous);

  free(cdf_loo);
	free_mat(matrix_Y_unordered_eval, num_reg_unordered);
	free_mat(matrix_Y_ordered_eval, num_reg_ordered);
	free_mat(matrix_Y_continuous_eval, num_reg_continuous);

	return(0);

}


/* 11/16/04 */

int kernel_estimate_con_distribution_categorical_no_mpi(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int num_obs_train,
int num_obs_eval,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_Y_unordered_eval,
double **matrix_Y_ordered_eval,
double **matrix_Y_continuous_eval,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *cdf,
double *cdf_stderr,
double small,
int itmax)
{

	/* This function estimates a density function using both continuous */
	/* and categorical covariates with three estimation techniques and an */
	/* assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel_cat;
	double prod_kernel_cont;

	double prod_kernel_marginal_cat;
	double prod_kernel_marginal_cont;

	double sum_ker;
	double sum_ker_marginal;

	double prod_h;

	double *lambda;
	double **matrix_bandwidth_var = NULL;
	double **matrix_bandwidth_reg = NULL;

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered);

	if((BANDWIDTH_den == 0)||(BANDWIDTH_den == 1))
	{
		matrix_bandwidth_var = alloc_matd(num_obs_eval,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_eval,num_reg_continuous);
	}
	else if(BANDWIDTH_den == 2)
	{
		matrix_bandwidth_var = alloc_matd(num_obs_train,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_train,num_reg_continuous);
	}

	/* Bandwidths for `dependent' variables */

	if(kernel_bandwidth_mean(
		KERNEL_den,
		BANDWIDTH_den,
		num_obs_train,
		num_obs_eval,
		num_var_continuous,
		num_var_unordered,
		num_var_ordered,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
		vector_scale_factor,
		matrix_Y_continuous_train,
		matrix_Y_continuous_eval,
		matrix_X_continuous_train,
		matrix_X_continuous_eval,
		matrix_bandwidth_var,
		matrix_bandwidth_reg,
		lambda) == 1)
	{

		printf("\n** Error: invalid bandwidth.");
		printf("\nProgram Terminated.\n");
		exit(EXIT_FAILURE);
	}

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=0; j < num_obs_eval; j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{

					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);

				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][0]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

			}

			cdf[j] = sum_ker/sum_ker_marginal;
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs_eval; j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{

					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);

				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][j]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

			}

			cdf[j] = sum_ker/sum_ker_marginal;
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}
	else
	{

		for(j=0; j < num_obs_eval; j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(i=0; i < num_obs_train; i++)
			{

				prod_h = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_h *= matrix_bandwidth_reg[l][i];
				}

				prod_kernel_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i]);
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][i]);
				}

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h;

			}

			cdf[j] = sum_ker/sum_ker_marginal;
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}

	free(lambda);

	free_mat(matrix_bandwidth_var,num_var_continuous);
	free_mat(matrix_bandwidth_reg,num_reg_continuous);

	return(0);

}


int kernel_estimate_con_density_categorical_gradient(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int num_obs_train,
int num_obs_eval,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_Y_unordered_eval,
double **matrix_Y_ordered_eval,
double **matrix_Y_continuous_eval,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_scale_factor,
int *num_categories,
double *pdf,
double *pdf_stderr,
double **pdf_deriv,
double **pdf_deriv_stderr,
double *log_likelihood)
{

	/* This function estimates a density function using both continuous */
	/* and categorical covariates with three estimation techniques and an */
	/* assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int k;
	int l;

	double prod_kernel_cat;
	double prod_kernel_cont;

	double prod_kernel_marginal_cat;
	double prod_kernel_marginal_cont;
	double prod_kernel_var_cont;

	double *prod_kernel_deriv;
	double *sum_ker_marginal_deriv;
	double *sum_ker_deriv;

	double sum_ker;
	double sum_ker_marginal;
	double sum_ker_temp;

	double prod_h;
	double prod_h_marginal;

	double *lambda;
	double **matrix_bandwidth_var = NULL;
	double **matrix_bandwidth_reg = NULL;

	/* Initialize constants for various kernels required for asymptotic standard errors */

	/* Integral of K(z)^p */
	double INT_KERNEL_P;
	/* Number of regressors times integral of K(z)^p */
	double K_INT_KERNEL_P;
	/* Integral of K(z-0.5)*K(z+0.5) */
	double INT_KERNEL_PM_HALF;
	/* Difference between int K(z)^p and int K(z-.5)K(z+.5) */
	double DIFF_KER_PPM;

	double log_DBL_MIN = log(DBL_MIN);

	#ifdef MPI2
	double log_likelihood_MPI;
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	/* Allocate memory for objects */

	prod_kernel_deriv = alloc_vecd(num_reg_continuous);
	sum_ker_marginal_deriv = alloc_vecd(num_reg_continuous);
	sum_ker_deriv = alloc_vecd(num_reg_continuous);

	lambda = alloc_vecd(num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered);

	if((BANDWIDTH_den == 0)||(BANDWIDTH_den == 1))
	{
		matrix_bandwidth_var = alloc_matd(num_obs_eval,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_eval,num_reg_continuous);
	}
	else if(BANDWIDTH_den == 2)
	{
		matrix_bandwidth_var = alloc_matd(num_obs_train,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_train,num_reg_continuous);
	}

	/* Bandwidths for `dependent' variables */

	if(kernel_bandwidth_mean(
		KERNEL_den,
		BANDWIDTH_den,
		num_obs_train,
		num_obs_eval,
		num_var_continuous,
		num_var_unordered,
		num_var_ordered,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
		vector_scale_factor,
		matrix_Y_continuous_train,
		matrix_Y_continuous_eval,
		matrix_X_continuous_train,
		matrix_X_continuous_eval,
		matrix_bandwidth_var,
		matrix_bandwidth_reg,
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

	/* Initialize constants for various kernels required for asymptotic standard errors */

	initialize_kernel_regression_asymptotic_constants(
		KERNEL_den,
		num_reg_continuous,
		&INT_KERNEL_P,
		&K_INT_KERNEL_P,
		&INT_KERNEL_PM_HALF,
		&DIFF_KER_PPM);

	#ifndef MPI2

	/* Initialize log likelihood */

	*log_likelihood = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			prod_h = 1.0;

			for(l = 0; l < num_var_continuous; l++)
			{
				prod_h *= matrix_bandwidth_var[l][0];
			}

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_marginal_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][0]);
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][0]);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			/* If sum_ker_marginal = 0.0 then so does sum_ker since it is the elements of sum_ker_marginal
																													 multiplied by other kernels */

			if(prod_h*sum_ker_marginal > 0.0)
			{
				pdf[j] = sum_ker/(prod_h*sum_ker_marginal);
			}
			else
			{
				pdf[j] = 0.0;
			}

			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j] = ((num_reg_continuous != 0) ? sqrt(pdf[j]*K_INT_KERNEL_P/prod_h) : sqrt(pdf[j]/prod_h));

			if(pdf[j] > DBL_MIN)
			{
				*log_likelihood += log(pdf[j]);
			}
			else
			{
				*log_likelihood += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_gradient()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(prod_h*sum_ker_marginal > 0.0)
				{
					pdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
						/(prod_h*sum_ker_marginal*matrix_bandwidth_reg[l][0]);
				}
				else
				{
					pdf_deriv[l][j] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */

				pdf_deriv_stderr[l][j] = sqrt(DIFF_KER_PPM /
					(sum_ker_marginal * ipow(matrix_bandwidth_reg[l][0],2)));
			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			prod_h = 1.0;

			for(l = 0; l < num_var_continuous; l++)
			{
				prod_h *= matrix_bandwidth_var[l][j];
			}

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_marginal_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][j]);
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][j]);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			/* If sum_ker_marginal = 0.0 then so does sum_ker since it is the elements of sum_ker_marginal
																													 multiplied by other kernels */

			if(prod_h*sum_ker_marginal > 0.0)
			{
				pdf[j] = sum_ker/(prod_h*sum_ker_marginal);
			}
			else
			{
				pdf[j] = 0.0;
			}

			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j] = ((num_reg_continuous != 0) ? sqrt(pdf[j]*K_INT_KERNEL_P/prod_h) : sqrt(pdf[j]/prod_h));

			if(pdf[j] > DBL_MIN)
			{
				*log_likelihood += log(pdf[j]);
			}
			else
			{
				*log_likelihood += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_gradient()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(prod_h*sum_ker_marginal > 0.0)
				{
					pdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
						/(prod_h*sum_ker_marginal*matrix_bandwidth_reg[l][j]);
				}
				else
				{
					pdf_deriv[l][j] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */

				pdf_deriv_stderr[l][j] = sqrt(DIFF_KER_PPM /
					(sum_ker_marginal * ipow(matrix_bandwidth_reg[l][j],2)));
			}

		}

	}
	else
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = sum_ker_temp = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_h = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_h *= matrix_bandwidth_var[l][i];
				}

				prod_h_marginal = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_h_marginal *= matrix_bandwidth_reg[l][i];
				}

				prod_h *= prod_h_marginal;

				prod_kernel_marginal_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i])/matrix_bandwidth_reg[l][i];
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][i])/matrix_bandwidth_var[l][i];
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h_marginal;
				sum_ker_temp += prod_kernel_cont*prod_kernel_cat/ipow(prod_h,2);

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][i])/ipow(matrix_bandwidth_reg[k][i],2);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i])/matrix_bandwidth_reg[l][i];
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			pdf[j] = sum_ker/sum_ker_marginal;
			/* stderr needs extra tmp I think */
			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j] = ((num_reg_continuous != 0) ? sqrt(sum_ker_temp*K_INT_KERNEL_P/sum_ker_marginal) : sqrt(sum_ker_temp/sum_ker_marginal));

			if(pdf[j] > DBL_MIN)
			{
				*log_likelihood += log(pdf[j]);
			}
			else
			{
				*log_likelihood += log_DBL_MIN;
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */
			/* 11/8/01 - removed prod_h from denom */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(sum_ker_marginal > 0.0)
				{
					pdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])/sum_ker_marginal;
				}
				else
				{
					pdf_deriv[l][j] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */
				/* 11/26/01 - dropped division by prod of bws, but need to revisit */

				pdf_deriv_stderr[l][j] = sqrt(DIFF_KER_PPM / sum_ker_marginal);
			}

		}

	}
	#endif

	#ifdef MPI2

	/* Initialize log likelihood */

	log_likelihood_MPI = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			prod_h = 1.0;

			for(l = 0; l < num_var_continuous; l++)
			{
				prod_h *= matrix_bandwidth_var[l][0];
			}

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_marginal_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][0]);
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][0]);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			/* If sum_ker_marginal = 0.0 then so does sum_ker since it is the elements of sum_ker_marginal
																													 multiplied by other kernels */

			if(prod_h*sum_ker_marginal > 0.0)
			{
				pdf[j-my_rank*stride] = sum_ker/(prod_h*sum_ker_marginal);
			}
			else
			{
				pdf[j-my_rank*stride] = 0.0;
			}

			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(pdf[j-my_rank*stride]*K_INT_KERNEL_P/prod_h) : sqrt(pdf[j-my_rank*stride]/prod_h));

			if(pdf[j-my_rank*stride] > DBL_MIN)
			{
				log_likelihood_MPI += log(pdf[j-my_rank*stride]);
			}
			else
			{
				log_likelihood_MPI += log_DBL_MIN;
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_gradient()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(prod_h*sum_ker_marginal > 0.0)
				{
					pdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
						/(prod_h*sum_ker_marginal*matrix_bandwidth_reg[l][0]);
				}
				else
				{
					pdf_deriv[l][j-my_rank*stride] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */

				pdf_deriv_stderr[l][j-my_rank*stride] = sqrt(DIFF_KER_PPM /
					(sum_ker_marginal * ipow(matrix_bandwidth_reg[l][0],2)));
			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			prod_h = 1.0;

			for(l = 0; l < num_var_continuous; l++)
			{
				prod_h *= matrix_bandwidth_var[l][j];
			}

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_marginal_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][j]);
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][j]);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			/* If sum_ker_marginal = 0.0 then so does sum_ker since it is the elements of sum_ker_marginal
																													 multiplied by other kernels */

			if(prod_h*sum_ker_marginal > 0.0)
			{
				pdf[j-my_rank*stride] = sum_ker/(prod_h*sum_ker_marginal);
			}
			else
			{
				pdf[j-my_rank*stride] = 0.0;
			}

			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(pdf[j-my_rank*stride]*K_INT_KERNEL_P/prod_h) : sqrt(pdf[j-my_rank*stride]/prod_h));

			if(pdf[j-my_rank*stride] > DBL_MIN)
			{
				log_likelihood_MPI += log(pdf[j-my_rank*stride]);
			}
			else
			{
				log_likelihood_MPI += log_DBL_MIN;
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_gradient()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(prod_h*sum_ker_marginal > 0.0)
				{
					pdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
						/(prod_h*sum_ker_marginal*matrix_bandwidth_reg[l][j]);
				}
				else
				{
					pdf_deriv[l][j-my_rank*stride] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */

				pdf_deriv_stderr[l][j-my_rank*stride] = sqrt(DIFF_KER_PPM /
					(sum_ker_marginal * ipow(matrix_bandwidth_reg[l][j],2)));
			}

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = sum_ker_temp = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_h = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_h *= matrix_bandwidth_var[l][i];
				}

				prod_h_marginal = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_h_marginal *= matrix_bandwidth_reg[l][i];
				}

				prod_h *= prod_h_marginal;

				prod_kernel_marginal_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i])/matrix_bandwidth_reg[l][i];
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][i])/matrix_bandwidth_var[l][i];
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h_marginal;
				sum_ker_temp += prod_kernel_cont*prod_kernel_cat/ipow(prod_h,2);

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][i])/ipow(matrix_bandwidth_reg[k][i],2);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i])/matrix_bandwidth_reg[l][i];
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			pdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;
			/* stderr needs extra tmp I think */
			/* With no continuous variables, need to drop K_INT_KERNEL_P */

			pdf_stderr[j-my_rank*stride] = ((num_reg_continuous != 0) ? sqrt(sum_ker_temp*K_INT_KERNEL_P/sum_ker_marginal) : sqrt(sum_ker_temp/sum_ker_marginal));

			if(pdf[j-my_rank*stride] > DBL_MIN)
			{
				log_likelihood_MPI += log(pdf[j-my_rank*stride]);
			}
			else
			{
				log_likelihood_MPI += log_DBL_MIN;
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(prod_h*sum_ker_marginal > 0.0)
				{
					pdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
						/sum_ker_marginal;
				}
				else
				{
					pdf_deriv[l][j-my_rank*stride] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */
				/* 11/26/01 - dropped prod of bws in denom, need to revisit */

				pdf_deriv_stderr[l][j-my_rank*stride] = sqrt(DIFF_KER_PPM / sum_ker_marginal);
			}

		}

	}

	MPI_Gather(pdf, stride, MPI_DOUBLE, pdf, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(pdf, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(pdf_stderr, stride, MPI_DOUBLE, pdf_stderr, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(pdf_stderr, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for(l = 0; l < num_reg_continuous; l++)
	{

		MPI_Gather(&pdf_deriv[l][0], stride, MPI_DOUBLE, &pdf_deriv[l][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&pdf_deriv[l][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&pdf_deriv_stderr[l][0], stride, MPI_DOUBLE, &pdf_deriv_stderr[l][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&pdf_deriv_stderr[l][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	}

	MPI_Reduce(&log_likelihood_MPI, log_likelihood, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(log_likelihood, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	free(prod_kernel_deriv);
	free(sum_ker_marginal_deriv);
	free(sum_ker_deriv);

	free(lambda);

	free_mat(matrix_bandwidth_var,num_var_continuous);
	free_mat(matrix_bandwidth_reg,num_reg_continuous);

	return(0);

}


int kernel_estimate_con_density_categorical_gradient_categorical(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int num_obs_train,
int num_obs_eval,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
int int_ordered_categorical_gradient,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_Y_unordered_eval,
double **matrix_Y_ordered_eval,
double **matrix_Y_continuous_eval,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_scale_factor,
double **matrix_categorical_vals,
int *num_categories,
double *pdf,
double **pdf_deriv,
double **pdf_deriv_stderr)
{

	/* This function computes the gradient for the unordered variables */

	int i;
	int j;
	int l;

	double *pdf_eval;
	double *pdf_stderr;

	double log_likelihood;

	double **matrix_X_unordered_temp;
	double **matrix_X_ordered_temp;

	double *pointer_m;
	double *pointer_me;
	double *pointer_g;

	#ifdef MPI2
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	#ifndef MPI2

	pdf_eval = alloc_vecd(num_obs_eval);
	pdf_stderr = alloc_vecd(num_obs_eval);
	matrix_X_unordered_temp = alloc_matd(num_obs_eval, num_reg_unordered);
	matrix_X_ordered_temp = alloc_matd(num_obs_eval, num_reg_ordered);

	for(i=0; i < num_reg_unordered; i++)
	{

		/* For each categorical variable */

		for(j=0; j < num_obs_eval; j++)
		{

			for(l=0; l < num_reg_unordered; l++)
			{
				matrix_X_unordered_temp[l][j] = matrix_X_unordered_eval[l][j];
			}

			/* We give the user a choice here... is gradient with respect to lower
																													 adjacent class or with respect to minimum value for categorical
																													 variable? One is appropriate for ordered categorical data, the other for
																													 unordered. Of course, the two can be inferred from the other. */

			if(int_ordered_categorical_gradient == 0)
			{

				/* Compute gradient with respect to minimum class */

				if(matrix_X_unordered_eval[i][j] != matrix_categorical_vals[i][0])
				{
					matrix_X_unordered_temp[i][j] = matrix_categorical_vals[i][0];
				}

			}
			else
			{

				/* Compute gradient as ordered categorical difference between adjacent class */

				if(matrix_X_unordered_eval[i][j] > matrix_categorical_vals[i][0])
				{
					for(l = 0; l < num_categories[i]; l++)
					{
						if(matrix_X_unordered_eval[i][j] == matrix_categorical_vals[i][l])
						{
							matrix_X_unordered_temp[i][j] = matrix_categorical_vals[i][l-1];
							l += num_categories[i];
						}
					}
				}

			}

		}

		kernel_estimate_con_density_categorical(
			KERNEL_den,
			KERNEL_unordered_den,
			KERNEL_ordered_den,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_den,
			num_obs_train,
			num_obs_eval,
			num_var_unordered,
			num_var_ordered,
			num_var_continuous,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_Y_unordered_train,
			matrix_Y_ordered_train,
			matrix_Y_continuous_train,
			matrix_Y_unordered_eval,
			matrix_Y_ordered_eval,
			matrix_Y_continuous_eval,
			matrix_X_unordered_train,
			matrix_X_ordered_train,
			matrix_X_continuous_train,
		/* Only difference... evaluation data is for categorical_eval */
			matrix_X_unordered_temp,
			matrix_X_ordered_eval,
			matrix_X_continuous_eval,
			vector_scale_factor,
			num_categories,
			pdf_eval,
			pdf_stderr,
			&log_likelihood);

		/* For ith categorical variable, gradient is discrete difference */
		/* between mean for sample observation and mean for minimum category */

		pointer_m = &pdf[0];
		pointer_me = &pdf_eval[0];
		pointer_g = &pdf_deriv[i][0];

		for(j=0; j < num_obs_eval; j++)
		{
			*pointer_g++ = *pointer_m++ - *pointer_me++;
			pdf_deriv_stderr[i][j]=0.0;
		}

	}

	for(i=0; i < num_reg_ordered; i++)
	{

		/* For each categorical variable */

		for(j=0; j < num_obs_eval; j++)
		{

			for(l=0; l < num_reg_ordered; l++)
			{
				matrix_X_ordered_temp[l][j] = matrix_X_ordered_eval[l][j];
			}

			/* We give the user a choice here... is gradient with respect to lower */
			/* adjacent class or with respect to minimum value for categorical */
			/* variable? One is appropriate for ordered categorical data, the other for */
			/* unordered. Of course, the two can be inferred from the other. */

			if(int_ordered_categorical_gradient == 0)
			{

				/* Compute gradient with respect to minimum class */

				if(matrix_X_ordered_eval[i][j] != matrix_categorical_vals[i+num_reg_unordered][0])
				{
					matrix_X_ordered_temp[i][j] = matrix_categorical_vals[i+num_reg_unordered][0];
				}

			}
			else
			{

				/* Compute gradient as ordered categorical difference between adjacent class */

				if(matrix_X_ordered_eval[i][j] > matrix_categorical_vals[i+num_reg_unordered][0])
				{
					for(l = 0; l < num_categories[i+num_reg_unordered]; l++)
					{
						if(matrix_X_ordered_eval[i][j] == matrix_categorical_vals[i+num_reg_unordered][l])
						{
							matrix_X_ordered_temp[i][j] = matrix_categorical_vals[i+num_reg_unordered][l-1];
							l += num_categories[i+num_reg_unordered];
						}
					}
				}

			}

		}

		kernel_estimate_con_density_categorical(
			KERNEL_den,
			KERNEL_unordered_den,
			KERNEL_ordered_den,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_den,
			num_obs_train,
			num_obs_eval,
			num_var_unordered,
			num_var_ordered,
			num_var_continuous,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_Y_unordered_train,
			matrix_Y_ordered_train,
			matrix_Y_continuous_train,
			matrix_Y_unordered_eval,
			matrix_Y_ordered_eval,
			matrix_Y_continuous_eval,
			matrix_X_unordered_train,
			matrix_X_ordered_train,
			matrix_X_continuous_train,
			matrix_X_unordered_eval,
		/* Only difference... evaluation data is for categorical_eval */
			matrix_X_ordered_temp,
			matrix_X_continuous_eval,
			vector_scale_factor,
			num_categories,
			pdf_eval,
			pdf_stderr,
			&log_likelihood);

		/* For ith categorical variable, gradient is discrete difference
																				 between mean for sample observation and mean for minimum category */

		pointer_m = &pdf[0];
		pointer_me = &pdf_eval[0];
		pointer_g = &pdf_deriv[i+num_reg_unordered][0];

		for(j=0; j < num_obs_eval; j++)
		{
			*pointer_g++ = *pointer_m++ - *pointer_me++;
			pdf_deriv_stderr[i+num_reg_unordered][j]=0.0;
		}

	}
	#endif

	#ifdef MPI2

	pdf_eval = alloc_vecd(stride*iNum_Processors);
	pdf_stderr = alloc_vecd(stride*iNum_Processors);
	matrix_X_unordered_temp = alloc_matd(stride*iNum_Processors, num_reg_unordered);
	matrix_X_ordered_temp = alloc_matd(stride*iNum_Processors, num_reg_ordered);

	for(i=0; i < num_reg_unordered; i++)
	{

		/* For each categorical variable */

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			for(l=0; l < num_reg_unordered; l++)
			{
				matrix_X_unordered_temp[l][j] = matrix_X_unordered_eval[l][j];
			}

			/* We give the user a choice here... is gradient with respect to lower
																													 adjacent class or with respect to minimum value for categorical
																													 variable? One is appropriate for ordered categorical data, the other for
																													 unordered. Of course, the two can be inferred from the other. */

			if(int_ordered_categorical_gradient == 0)
			{

				/* Compute gradient with respect to minimum class */

				if(matrix_X_unordered_eval[i][j] != matrix_categorical_vals[i][0])
				{
					matrix_X_unordered_temp[i][j] = matrix_categorical_vals[i][0];
				}

			}
			else
			{

				/* Compute gradient as ordered categorical difference between adjacent class */

				if(matrix_X_unordered_eval[i][j] > matrix_categorical_vals[i][0])
				{
					for(l = 0; l < num_categories[i]; l++)
					{
						if(matrix_X_unordered_eval[i][j] == matrix_categorical_vals[i][l])
						{
							matrix_X_unordered_temp[i][j] = matrix_categorical_vals[i][l-1];
							l += num_categories[i];
						}
					}
				}

			}

		}

		kernel_estimate_con_density_categorical(
			KERNEL_den,
			KERNEL_unordered_den,
			KERNEL_ordered_den,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_den,
			num_obs_train,
			num_obs_eval,
			num_var_unordered,
			num_var_ordered,
			num_var_continuous,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_Y_unordered_train,
			matrix_Y_ordered_train,
			matrix_Y_continuous_train,
			matrix_Y_unordered_eval,
			matrix_Y_ordered_eval,
			matrix_Y_continuous_eval,
			matrix_X_unordered_train,
			matrix_X_ordered_train,
			matrix_X_continuous_train,
		/* Only difference... evaluation data is for categorical_eval */
			matrix_X_unordered_temp,
			matrix_X_ordered_eval,
			matrix_X_continuous_eval,
			vector_scale_factor,
			num_categories,
			pdf_eval,
			pdf_stderr,
			&log_likelihood);

		/* For ith categorical variable, gradient is discrete difference */
		/* between mean for sample observation and mean for minimum category */

		/* Fed location beginning with num_cont */

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{
			pdf_deriv[i][j-my_rank*stride] = pdf[j] - pdf_eval[j];;
			pdf_deriv_stderr[i][j-my_rank*stride]=0.0;
		}

	}

	for(i=0; i < num_reg_ordered; i++)
	{

		/* For each categorical variable */

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			for(l=0; l < num_reg_ordered; l++)
			{
				matrix_X_ordered_temp[l][j] = matrix_X_ordered_eval[l][j];
			}

			/* We give the user a choice here... is gradient with respect to lower */
			/* adjacent class or with respect to minimum value for categorical */
			/* variable? One is appropriate for ordered categorical data, the other for */
			/* unordered. Of course, the two can be inferred from the other. */

			if(int_ordered_categorical_gradient == 0)
			{

				/* Compute gradient with respect to minimum class */

				if(matrix_X_ordered_eval[i][j] != matrix_categorical_vals[i+num_reg_unordered][0])
				{
					matrix_X_ordered_temp[i][j] = matrix_categorical_vals[i+num_reg_unordered][0];
				}

			}
			else
			{

				/* Compute gradient as ordered categorical difference between adjacent class */

				if(matrix_X_ordered_eval[i][j] > matrix_categorical_vals[i+num_reg_unordered][0])
				{
					for(l = 0; l < num_categories[i+num_reg_unordered]; l++)
					{
						if(matrix_X_ordered_eval[i][j] == matrix_categorical_vals[i+num_reg_unordered][l])
						{
							matrix_X_ordered_temp[i][j] = matrix_categorical_vals[i+num_reg_unordered][l-1];
							l += num_categories[i+num_reg_unordered];
						}
					}
				}

			}

		}

		kernel_estimate_con_density_categorical(
			KERNEL_den,
			KERNEL_unordered_den,
			KERNEL_ordered_den,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_den,
			num_obs_train,
			num_obs_eval,
			num_var_unordered,
			num_var_ordered,
			num_var_continuous,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_Y_unordered_train,
			matrix_Y_ordered_train,
			matrix_Y_continuous_train,
			matrix_Y_unordered_eval,
			matrix_Y_ordered_eval,
			matrix_Y_continuous_eval,
			matrix_X_unordered_train,
			matrix_X_ordered_train,
			matrix_X_continuous_train,
			matrix_X_unordered_eval,
		/* Only difference... evaluation data is for categorical_eval */
			matrix_X_ordered_temp,
			matrix_X_continuous_eval,
			vector_scale_factor,
			num_categories,
			pdf_eval,
			pdf_stderr,
			&log_likelihood);

		/* For ith categorical variable, gradient is discrete difference */
		/* between mean for sample observation and mean for minimum category */

		/* Fed location beginning with num_cont */

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{
			pdf_deriv[i+num_reg_unordered][j-my_rank*stride] = pdf[j] - pdf_eval[j];
			pdf_deriv_stderr[i+num_reg_unordered][j-my_rank*stride]=0.0;
		}

	}

	for(l = 0; l < num_reg_unordered+num_reg_ordered; l++)
	{

		MPI_Gather(&pdf_deriv[l][0], stride, MPI_DOUBLE, &pdf_deriv[l][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&pdf_deriv[l][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&pdf_deriv_stderr[l][0], stride, MPI_DOUBLE, &pdf_deriv_stderr[l][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&pdf_deriv_stderr[l][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	}
	#endif

	free(pdf_eval);
	free(pdf_stderr);
	free_mat(matrix_X_unordered_temp, num_reg_unordered);
	free_mat(matrix_X_ordered_temp, num_reg_ordered);

	return(0);

}


/* 11/24/99 - new, definitely not debugged - need to fix con_density_cat_grad() */
/* 11/8/01 - grads finally debugged */

int kernel_estimate_con_distribution_categorical_gradient(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int num_obs_train,
int num_obs_eval,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_Y_unordered_eval,
double **matrix_Y_ordered_eval,
double **matrix_Y_continuous_eval,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *cdf,
double *cdf_stderr,
double **cdf_deriv,
double **cdf_deriv_stderr,
double small,
int itmax)
{

	/* This function estimates a density function using both continuous */
	/* and categorical covariates with three estimation techniques and an */
	/* assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int k;
	int l;

	double prod_kernel_cat;
	double prod_kernel_cont;

	double prod_kernel_marginal_cat;
	double prod_kernel_marginal_cont;
	double prod_kernel_var_cont;

	double *prod_kernel_deriv;
	double *sum_ker_marginal_deriv;
	double *sum_ker_deriv;

	double sum_ker;
	double sum_ker_marginal;

	double *lambda;
	double **matrix_bandwidth_var = NULL;
	double **matrix_bandwidth_reg = NULL;

	/* Initialize constants for various kernels required for asymptotic standard errors */

	/* Integral of K(z)^p */
	double INT_KERNEL_P;
	/* Number of regressors times integral of K(z)^p */
	double K_INT_KERNEL_P;
	/* Integral of K(z-0.5)*K(z+0.5) */
	double INT_KERNEL_PM_HALF;
	/* Difference between int K(z)^p and int K(z-.5)K(z+.5) */
	double DIFF_KER_PPM;

	#ifdef MPI2
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	/* Allocate memory for objects */

	prod_kernel_deriv = alloc_vecd(num_reg_continuous);
	sum_ker_marginal_deriv = alloc_vecd(num_reg_continuous);
	sum_ker_deriv = alloc_vecd(num_reg_continuous);

	lambda = alloc_vecd(num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered);

	if((BANDWIDTH_den == 0)||(BANDWIDTH_den == 1))
	{
		matrix_bandwidth_var = alloc_matd(num_obs_eval,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_eval,num_reg_continuous);
	}
	else if(BANDWIDTH_den == 2)
	{
		matrix_bandwidth_var = alloc_matd(num_obs_train,num_var_continuous);
		matrix_bandwidth_reg = alloc_matd(num_obs_train,num_reg_continuous);
	}

	/* Bandwidths for `dependent' variables */

	if(kernel_bandwidth_mean(
		KERNEL_den,
		BANDWIDTH_den,
		num_obs_train,
		num_obs_eval,
		num_var_continuous,
		num_var_unordered,
		num_var_ordered,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
		vector_scale_factor,
		matrix_Y_continuous_train,
		matrix_Y_continuous_eval,
		matrix_X_continuous_train,
		matrix_X_continuous_eval,
		matrix_bandwidth_var,
		matrix_bandwidth_reg,
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

	/* Initialize constants for various kernels required for asymptotic standard errors */

	initialize_kernel_regression_asymptotic_constants(
		KERNEL_den,
		num_reg_continuous,
		&INT_KERNEL_P,
		&K_INT_KERNEL_P,
		&INT_KERNEL_PM_HALF,
		&DIFF_KER_PPM);

	#ifndef MPI2

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_marginal_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][0]);
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][0]);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			/* If sum_ker_marginal = 0.0 then so does sum_ker since it is the elements of */
			/* sum_ker_marginal multiplied by other kernels */

			if(sum_ker_marginal > 0.0)
			{
				cdf[j] = sum_ker/sum_ker_marginal;
			}
			else
			{
				cdf[j] = 0.0;
			}
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(sum_ker_marginal > 0.0)
				{
					cdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
						/(sum_ker_marginal*matrix_bandwidth_reg[l][0]);
				}
				else
				{
					cdf_deriv[l][j] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */

				cdf_deriv_stderr[l][j] = sqrt(DIFF_KER_PPM /
					(sum_ker_marginal * ipow(matrix_bandwidth_reg[l][0],2)));
			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_marginal_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][j]);
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][j]);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			/* If sum_ker_marginal = 0.0 then so does sum_ker since it is the elements of */
			/* sum_ker_marginal multiplied by other kernels */

			if(sum_ker_marginal > 0.0)
			{
				cdf[j] = sum_ker/sum_ker_marginal;
			}
			else
			{
				cdf[j] = 0.0;
			}
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(sum_ker_marginal > 0.0)
				{
					cdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
						/(sum_ker_marginal*matrix_bandwidth_reg[l][j]);
				}
				else
				{
					cdf_deriv[l][j] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */

				cdf_deriv_stderr[l][j] = sqrt(DIFF_KER_PPM /
					(sum_ker_marginal * ipow(matrix_bandwidth_reg[l][j],2)));
			}

		}

	}
	else
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_marginal_cont = 1.0;

				/* Only reg gets did by h... for var integral kernel cancels */

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i])/matrix_bandwidth_reg[l][i];
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][i]);
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][i])/ipow(matrix_bandwidth_reg[k][i],2);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i])/matrix_bandwidth_reg[l][i];
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			/* If sum_ker_marginal = 0.0 then so does sum_ker since it is the elements of */
			/* sum_ker_marginal multiplied by other kernels */

			if(sum_ker_marginal > 0.0)
			{
				cdf[j] = sum_ker/sum_ker_marginal;
			}
			else
			{
				cdf[j] = 0.0;
			}

			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(sum_ker_marginal > 0.0)
				{
					cdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])/sum_ker_marginal;
				}
				else
				{
					cdf_deriv[l][j] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */
				/* 11/28/01 - removed h^2 for adaptive due to alloc issues */

				cdf_deriv_stderr[l][j] = sqrt(DIFF_KER_PPM / sum_ker_marginal);

			}

		}

	}
	#endif

	#ifdef MPI2

	/* Conduct the estimation */

	/* 4/22/2002 - code was not working for np > 1. Cut code from above, re-implemented */

	if(BANDWIDTH_den == 0)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_marginal_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][0]);
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][0]);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][0]);
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			/* If sum_ker_marginal = 0.0 then so does sum_ker since it is the elements of */
			/* sum_ker_marginal multiplied by other kernels */

			if(sum_ker_marginal > 0.0)
			{
				cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;
			}
			else
			{
				cdf[j-my_rank*stride] = 0.0;
			}
			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(sum_ker_marginal > 0.0)
				{
					cdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
						/(sum_ker_marginal*matrix_bandwidth_reg[l][0]);
				}
				else
				{
					cdf_deriv[l][j-my_rank*stride] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */

				cdf_deriv_stderr[l][j-my_rank*stride] = sqrt(DIFF_KER_PPM /
					(sum_ker_marginal * ipow(matrix_bandwidth_reg[l][0],2)));
			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_marginal_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][j]);
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][j]);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][j]);
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			/* If sum_ker_marginal = 0.0 then so does sum_ker since it is the elements of */
			/* sum_ker_marginal multiplied by other kernels */

			if(sum_ker_marginal > 0.0)
			{
				cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;
			}
			else
			{
				cdf[j-my_rank*stride] = 0.0;
			}
			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(sum_ker_marginal > 0.0)
				{
					cdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
						/(sum_ker_marginal*matrix_bandwidth_reg[l][j]);
				}
				else
				{
					cdf_deriv[l][j-my_rank*stride] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */

				cdf_deriv_stderr[l][j-my_rank*stride] = sqrt(DIFF_KER_PPM /
					(sum_ker_marginal * ipow(matrix_bandwidth_reg[l][j],2)));
			}

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = DBL_MIN;

			for(l = 0; l < num_reg_continuous; l++)
			{
				sum_ker_deriv[l] = sum_ker_marginal_deriv[l] = 0.0;
			}

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel_marginal_cont = 1.0;

				/* Only reg gets did by h... for var integral kernel cancels */

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_marginal_cont *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i])/matrix_bandwidth_reg[l][i];
				}

				prod_kernel_var_cont = 1.0;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_var_cont *= cdf_kernel(KERNEL_den, (matrix_Y_continuous_eval[l][j]-matrix_Y_continuous_train[l][i])/matrix_bandwidth_var[l][i]);
				}

				prod_kernel_cont = prod_kernel_marginal_cont*prod_kernel_var_cont;

				prod_kernel_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_eval[l][j],matrix_X_unordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_eval[l][j],matrix_X_ordered_train[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered_eval[l][j],matrix_Y_unordered_train[l][i],lambda[l],num_categories[l],matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= cdf_kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered_eval[l][j],matrix_Y_ordered_train[l][i],lambda[l+num_var_unordered],num_categories[l+num_var_unordered],matrix_categorical_vals[l+num_var_unordered]);
				}

				sum_ker += prod_kernel_cont*prod_kernel_cat;
				sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

				/* Kernels for derivatives */

				for (k = 0; k < num_reg_continuous; k++)
				{
					prod_kernel_deriv[k] = kernel_deriv(KERNEL_reg, (matrix_X_continuous_eval[k][j]-matrix_X_continuous_train[k][i])/matrix_bandwidth_reg[k][i])/ipow(matrix_bandwidth_reg[k][i],2);
					for (l = 0; l < num_reg_continuous; l++)
					{
						if(l != k)
						{
							prod_kernel_deriv[k] *= kernel(KERNEL_reg, (matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth_reg[l][i])/matrix_bandwidth_reg[l][i];
						}
					}
				}

				for(l = 0; l < num_reg_continuous; l++)
				{
					sum_ker_deriv[l] += prod_kernel_cat * prod_kernel_var_cont * prod_kernel_deriv[l];
					sum_ker_marginal_deriv[l] += prod_kernel_marginal_cat * prod_kernel_deriv[l];
				}

			}

			/* If sum_ker_marginal = 0.0 then so does sum_ker since it is the elements of */
			/* sum_ker_marginal multiplied by other kernels */

			if(sum_ker_marginal > 0.0)
			{
				cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;
			}
			else
			{
				cdf[j-my_rank*stride] = 0.0;
			}

			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

				if(sum_ker_marginal > 0.0)
				{
					cdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])/sum_ker_marginal;
				}
				else
				{
					cdf_deriv[l][j-my_rank*stride] = 0.0;
				}

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */
				/* 11/28/01 - removed h^2 for adaptive due to alloc issues */

				cdf_deriv_stderr[l][j-my_rank*stride] = sqrt(DIFF_KER_PPM / sum_ker_marginal);

			}

		}

	}

	MPI_Gather(cdf, stride, MPI_DOUBLE, cdf, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(cdf, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(cdf_stderr, stride, MPI_DOUBLE, cdf_stderr, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(cdf_stderr, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for(l = 0; l < num_reg_continuous; l++)
	{

		MPI_Gather(&cdf_deriv[l][0], stride, MPI_DOUBLE, &cdf_deriv[l][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&cdf_deriv[l][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&cdf_deriv_stderr[l][0], stride, MPI_DOUBLE, &cdf_deriv_stderr[l][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&cdf_deriv_stderr[l][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	}
	#endif

	free(prod_kernel_deriv);
	free(sum_ker_marginal_deriv);
	free(sum_ker_deriv);

	free(lambda);

	free_mat(matrix_bandwidth_var,num_var_continuous);
	free_mat(matrix_bandwidth_reg,num_reg_continuous);

	return(0);

}


int kernel_estimate_con_distribution_categorical_gradient_categorical(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int num_obs_train,
int num_obs_eval,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
int int_ordered_categorical_gradient,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_Y_unordered_eval,
double **matrix_Y_ordered_eval,
double **matrix_Y_continuous_eval,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_scale_factor,
double **matrix_categorical_vals,
int *num_categories,
double *cdf,
double **cdf_deriv,
double **cdf_deriv_stderr,
double small,
int itmax)
{

	/* This function computes the gradient for the unordered variables */

	int i;
	int j;
	int l;

	double *cdf_eval;
	double *cdf_stderr;

	double **matrix_X_unordered_temp;
	double **matrix_X_ordered_temp;

	double *pointer_m;
	double *pointer_me;
	double *pointer_g;

	#ifdef MPI2
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	#ifndef MPI2

	cdf_eval = alloc_vecd(num_obs_eval);
	cdf_stderr = alloc_vecd(num_obs_eval);
	matrix_X_unordered_temp = alloc_matd(num_obs_eval, num_reg_unordered);
	matrix_X_ordered_temp = alloc_matd(num_obs_eval, num_reg_ordered);

	for(i=0; i < num_reg_unordered; i++)
	{

		/* For each categorical variable */

		for(j=0; j < num_obs_eval; j++)
		{

			for(l=0; l < num_reg_unordered; l++)
			{
				matrix_X_unordered_temp[l][j] = matrix_X_unordered_eval[l][j];
			}

			/* We give the user a choice here... is gradient with respect to lower */
			/* adjacent class or with respect to minimum value for categorical */
			/* variable? One is appropriate for ordered categorical data, the other for */
			/* unordered. Of course, the two can be inferred from the other. */

			if(int_ordered_categorical_gradient == 0)
			{

				/* Compute gradient with respect to minimum class */

				if(matrix_X_unordered_eval[i][j] != matrix_categorical_vals[i+num_var_unordered][0])
				{
					matrix_X_unordered_temp[i][j] = matrix_categorical_vals[i+num_var_unordered][0];
				}

			}
			else
			{

				/* Compute gradient as ordered categorical difference between adjacent class */

				if(matrix_X_unordered_eval[i][j] > matrix_categorical_vals[i+num_var_unordered][0])
				{
					for(l = 0; l < num_categories[i+num_var_unordered]; l++)
					{
						if(matrix_X_unordered_eval[i][j] == matrix_categorical_vals[i+num_var_unordered][l])
						{
							matrix_X_unordered_temp[i][j] = matrix_categorical_vals[i+num_var_unordered][l-1];
							l += num_categories[i+num_var_unordered];
						}
					}
				}

			}

		}

		kernel_estimate_con_distribution_categorical(
			KERNEL_den,
			KERNEL_unordered_den,
			KERNEL_ordered_den,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_den,
			num_obs_train,
			num_obs_eval,
			num_var_unordered,
			num_var_ordered,
			num_var_continuous,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_Y_unordered_train,
			matrix_Y_ordered_train,
			matrix_Y_continuous_train,
			matrix_Y_unordered_eval,
			matrix_Y_ordered_eval,
			matrix_Y_continuous_eval,
			matrix_X_unordered_train,
			matrix_X_ordered_train,
			matrix_X_continuous_train,
		/* Only difference... evaluation data is for categorical_eval */
			matrix_X_unordered_temp,
			matrix_X_ordered_eval,
			matrix_X_continuous_eval,
			vector_scale_factor,
			num_categories,
			matrix_categorical_vals,
			cdf_eval,
			cdf_stderr,
			small,
			itmax);

		/* For ith categorical variable, gradient is discrete difference */
		/* between mean for sample observation and mean for minimum category */

		pointer_m = &cdf[0];
		pointer_me = &cdf_eval[0];
		pointer_g = &cdf_deriv[i][0];

		for(j=0; j < num_obs_eval; j++)
		{
			*pointer_g++ = *pointer_m++ - *pointer_me++;
			cdf_deriv_stderr[i][j]=0.0;
		}

	}

	for(i=0; i < num_reg_ordered; i++)
	{

		/* For each categorical variable */

		for(j=0; j < num_obs_eval; j++)
		{

			for(l=0; l < num_reg_ordered; l++)
			{
				matrix_X_ordered_temp[l][j] = matrix_X_ordered_eval[l][j];
			}

			/* We give the user a choice here... is gradient with respect to lower */
			/* adjacent class or with respect to minimum value for categorical */
			/* variable? One is appropriate for ordered categorical data, the other for */
			/* unordered. Of course, the two can be inferred from the other. */

			if(int_ordered_categorical_gradient == 0)
			{

				/* Compute gradient with respect to minimum class */

				if(matrix_X_ordered_eval[i][j] != matrix_categorical_vals[i+num_var_unordered+num_var_ordered][0])
				{
					matrix_X_ordered_temp[i][j] = matrix_categorical_vals[i+num_var_unordered+num_var_ordered][0];
				}

			}
			else
			{

				/* Compute gradient as ordered categorical difference between adjacent class */

				if(matrix_X_ordered_eval[i][j] > matrix_categorical_vals[i+num_var_unordered+num_var_ordered][0])
				{
					for(l = 0; l < num_categories[i+num_var_unordered+num_var_ordered]; l++)
					{
						if(matrix_X_ordered_eval[i][j] == matrix_categorical_vals[i+num_var_unordered+num_var_ordered][l])
						{
							matrix_X_ordered_temp[i][j] = matrix_categorical_vals[i+num_var_unordered+num_var_ordered][l-1];
							l += num_categories[i+num_var_unordered+num_var_ordered];
						}
					}
				}

			}

		}

		kernel_estimate_con_distribution_categorical(
			KERNEL_den,
			KERNEL_unordered_den,
			KERNEL_ordered_den,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_den,
			num_obs_train,
			num_obs_eval,
			num_var_unordered,
			num_var_ordered,
			num_var_continuous,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_Y_unordered_train,
			matrix_Y_ordered_train,
			matrix_Y_continuous_train,
			matrix_Y_unordered_eval,
			matrix_Y_ordered_eval,
			matrix_Y_continuous_eval,
			matrix_X_unordered_train,
			matrix_X_ordered_train,
			matrix_X_continuous_train,
			matrix_X_unordered_eval,
		/* Only difference... evaluation data is for categorical_eval */
			matrix_X_ordered_temp,
			matrix_X_continuous_eval,
			vector_scale_factor,
			num_categories,
			matrix_categorical_vals,
			cdf_eval,
			cdf_stderr,
			small,
			itmax);

		/* For ith categorical variable, gradient is discrete difference */
		/* between mean for sample observation and mean for minimum category */

		pointer_m = &cdf[0];
		pointer_me = &cdf_eval[0];
		pointer_g = &cdf_deriv[i+num_reg_unordered][0];

		for(j=0; j < num_obs_eval; j++)
		{
			*pointer_g++ = *pointer_m++ - *pointer_me++;
			cdf_deriv_stderr[i+num_reg_unordered][j]=0.0;
		}

	}
	#endif

	#ifdef MPI2

	cdf_eval = alloc_vecd(stride*iNum_Processors);
	cdf_stderr = alloc_vecd(stride*iNum_Processors);
	matrix_X_unordered_temp = alloc_matd(stride*iNum_Processors, num_reg_unordered);
	matrix_X_ordered_temp = alloc_matd(stride*iNum_Processors, num_reg_ordered);

	for(i=0; i < num_reg_unordered; i++)
	{

		/* For each categorical variable */

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			for(l=0; l < num_reg_unordered; l++)
			{
				matrix_X_unordered_temp[l][j] = matrix_X_unordered_eval[l][j];
			}

			/* We give the user a choice here... is gradient with respect to lower */
			/* adjacent class or with respect to minimum value for categorical */
			/* variable? One is appropriate for ordered categorical data, the other for */
			/* unordered. Of course, the two can be inferred from the other. */

			if(int_ordered_categorical_gradient == 0)
			{

				/* Compute gradient with respect to minimum class */

				if(matrix_X_unordered_eval[i][j] != matrix_categorical_vals[i+num_var_unordered][0])
				{
					matrix_X_unordered_temp[i][j] = matrix_categorical_vals[i+num_var_unordered][0];
				}

			}
			else
			{

				/* Compute gradient as ordered categorical difference between adjacent class */

				if(matrix_X_unordered_eval[i][j] > matrix_categorical_vals[i+num_var_unordered][0])
				{
					for(l = 0; l < num_categories[i+num_var_unordered]; l++)
					{
						if(matrix_X_unordered_eval[i][j] == matrix_categorical_vals[i+num_var_unordered][l])
						{
							matrix_X_unordered_temp[i][j] = matrix_categorical_vals[i+num_var_unordered][l-1];
							l += num_categories[i+num_var_unordered];
						}
					}
				}

			}

		}

		kernel_estimate_con_distribution_categorical(
			KERNEL_den,
			KERNEL_unordered_den,
			KERNEL_ordered_den,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_den,
			num_obs_train,
			num_obs_eval,
			num_var_unordered,
			num_var_ordered,
			num_var_continuous,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_Y_unordered_train,
			matrix_Y_ordered_train,
			matrix_Y_continuous_train,
			matrix_Y_unordered_eval,
			matrix_Y_ordered_eval,
			matrix_Y_continuous_eval,
			matrix_X_unordered_train,
			matrix_X_ordered_train,
			matrix_X_continuous_train,
		/* Only difference... evaluation data is for categorical_eval */
			matrix_X_unordered_temp,
			matrix_X_ordered_eval,
			matrix_X_continuous_eval,
			vector_scale_factor,
			num_categories,
			matrix_categorical_vals,
			cdf_eval,
			cdf_stderr,
			small,
			itmax);

		/* For ith categorical variable, gradient is discrete difference */
		/* between mean for sample observation and mean for minimum category */

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{
			cdf_deriv[i][j-my_rank*stride] = cdf[j] - cdf_eval[j];;
			cdf_deriv_stderr[i][j-my_rank*stride]=0.0;
		}

	}

	for(i=0; i < num_reg_ordered; i++)
	{

		/* For each categorical variable */

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			for(l=0; l < num_reg_ordered; l++)
			{
				matrix_X_ordered_temp[l][j] = matrix_X_ordered_eval[l][j];
			}

			/* We give the user a choice here... is gradient with respect to lower */
			/* adjacent class or with respect to minimum value for categorical */
			/* variable? One is appropriate for ordered categorical data, the other for */
			/* unordered. Of course, the two can be inferred from the other. */

			if(int_ordered_categorical_gradient == 0)
			{

				/* Compute gradient with respect to minimum class */

				if(matrix_X_ordered_eval[i][j] != matrix_categorical_vals[i+num_var_unordered+num_var_ordered][0])
				{
					matrix_X_ordered_temp[i][j] = matrix_categorical_vals[i+num_var_unordered+num_var_ordered][0];
				}

			}
			else
			{

				/* Compute gradient as ordered categorical difference between adjacent class */

				if(matrix_X_ordered_eval[i][j] > matrix_categorical_vals[i+num_var_unordered+num_var_ordered][0])
				{
					for(l = 0; l < num_categories[i+num_var_unordered+num_var_ordered]; l++)
					{
						if(matrix_X_ordered_eval[i][j] == matrix_categorical_vals[i+num_var_unordered+num_var_ordered][l])
						{
							matrix_X_ordered_temp[i][j] = matrix_categorical_vals[i+num_var_unordered+num_var_ordered][l-1];
							l += num_categories[i+num_var_unordered+num_var_ordered];
						}
					}
				}

			}

		}

		kernel_estimate_con_distribution_categorical(
			KERNEL_den,
			KERNEL_unordered_den,
			KERNEL_ordered_den,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_den,
			num_obs_train,
			num_obs_eval,
			num_var_unordered,
			num_var_ordered,
			num_var_continuous,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_Y_unordered_train,
			matrix_Y_ordered_train,
			matrix_Y_continuous_train,
			matrix_Y_unordered_eval,
			matrix_Y_ordered_eval,
			matrix_Y_continuous_eval,
			matrix_X_unordered_train,
			matrix_X_ordered_train,
			matrix_X_continuous_train,
			matrix_X_unordered_eval,
		/* Only difference... evaluation data is for categorical_eval */
			matrix_X_ordered_temp,
			matrix_X_continuous_eval,
			vector_scale_factor,
			num_categories,
			matrix_categorical_vals,
			cdf_eval,
			cdf_stderr,
			small,
			itmax);

		/* For ith categorical variable, gradient is discrete difference */
		/* between mean for sample observation and mean for minimum category */

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{
			cdf_deriv[i+num_reg_unordered][j-my_rank*stride] = cdf[j] - cdf_eval[j];;
			cdf_deriv_stderr[i+num_reg_unordered][j-my_rank*stride]=0.0;
		}

	}

	for(l = 0; l < num_reg_unordered+num_reg_ordered; l++)
	{

		MPI_Gather(&cdf_deriv[l][0], stride, MPI_DOUBLE, &cdf_deriv[l][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&cdf_deriv[l][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&cdf_deriv_stderr[l][0], stride, MPI_DOUBLE, &cdf_deriv_stderr[l][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&cdf_deriv_stderr[l][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	}
	#endif

	free(cdf_eval);
	free(cdf_stderr);
	free_mat(matrix_X_unordered_temp, num_reg_unordered);
	free_mat(matrix_X_ordered_temp, num_reg_ordered);

	return(0);

}


/* Convolution kernel implemented...*/

int kernel_estimate_density_categorical_convolution_cv(
int KERNEL_den,
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
double *cv)
{

	/* This function estimates a leave one out density function using both */
	/* continuous and categorical covariates with three estimation techniques */
	/* and an assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel;
	double prod_kernel_convol;

	double sum_ker;
	double sum_ker_convol;

	double temp_bw1;
	double temp;

	double *lambda;
	double **matrix_bandwidth;

	double *p_xj1;
	double *p_xi1;
	double temp_xj;

	double n_sq_inv = 1.0/ipow(num_obs,2);
	double n_times_n_minus_1_inv = 1.0/((double)num_obs*(num_obs-1.0));

	double nh_sq_inv;
	double nh_times_n_minus_1_inv;

	/* Allocate memory for objects */

	#ifdef MPI2
	double cv_MPI;
	int stride = ceil((double) num_obs / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);
	matrix_bandwidth = alloc_matd(num_obs,num_reg_continuous);

	/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

	if(kernel_bandwidth_mean(
		KERNEL_den,
		BANDWIDTH_den,
		num_obs,
		num_obs,
		0,
		0,
		0,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
		vector_scale_factor,
		/* Not used */
		matrix_X_continuous,
		/* Not used */
		matrix_X_continuous,
		matrix_X_continuous,
		matrix_X_continuous,
		/* Not used */
		matrix_bandwidth,
		matrix_bandwidth,
		lambda)==1)
	{
		free(lambda);
		free_mat(matrix_bandwidth,num_reg_continuous);
		return(1);
	}

	#ifndef MPI2

	/* Initialize cv function */

	*cv = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		if((num_reg_continuous == 1)&&(num_reg_unordered + num_reg_ordered == 0))
		{

			/* For special case of fixed bandwidth for PDF/CDF for test stat */

			temp_bw1 = matrix_bandwidth[0][0];

			nh_sq_inv = n_sq_inv/temp_bw1;
			nh_times_n_minus_1_inv = n_times_n_minus_1_inv/temp_bw1;

			/* For all i and j */

			p_xj1 = &matrix_X_continuous[0][0];

			for (j = 0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();

				sum_ker = 0.0;
				sum_ker_convol = 0.0;

				/* Move unnecessary pointer out of loop (VTune) */

				temp_xj = *p_xj1++;
				p_xi1 = &matrix_X_continuous[0][0];

				for (i =  0; i <  num_obs; i++)
				{

					temp = (temp_xj - *p_xi1++)/temp_bw1;

					sum_ker_convol += kernel_convol(KERNEL_den,BANDWIDTH_den,temp,temp_bw1,temp_bw1);

					if(i != j)
					{
						sum_ker += kernel(KERNEL_den,temp);
					}

				}

				*cv += sum_ker_convol*nh_sq_inv-2.0*sum_ker*nh_times_n_minus_1_inv;

			}

		}
		else
		{

			for(j=0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();

				sum_ker = 0.0;
				sum_ker_convol = 0.0;

				for(i=0; i < num_obs; i++)
				{

					prod_kernel = 1.0;
					prod_kernel_convol = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][0])/matrix_bandwidth[l][0];
						prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
							(matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][0],matrix_bandwidth[l][0],matrix_bandwidth[l][0])
							/matrix_bandwidth[l][0];
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
						prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l], matrix_categorical_vals[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
						prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered], matrix_categorical_vals[l+num_reg_unordered]);
					}

					if(i != j)
					{
						sum_ker += prod_kernel;
					}

					sum_ker_convol += prod_kernel_convol;

				}

				*cv += sum_ker_convol*n_sq_inv-2.0*sum_ker*n_times_n_minus_1_inv;

			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs; j++)
		{
		  R_CheckUserInterrupt();

			sum_ker = 0.0;
			sum_ker_convol = 0.0;

			for(i=0; i < num_obs; i++)
			{

				prod_kernel = 1.0;
				prod_kernel_convol = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][j])/matrix_bandwidth[l][j];
					prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
						(matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][j],matrix_bandwidth[l][i],matrix_bandwidth[l][j])
						/matrix_bandwidth[l][j];
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l], matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered], matrix_categorical_vals[l+num_reg_unordered]);
				}

				if(i != j)
				{
					sum_ker += prod_kernel;
				}

				sum_ker_convol += prod_kernel_convol;

			}

			*cv += sum_ker_convol*n_sq_inv-2.0*sum_ker*n_times_n_minus_1_inv;

		}

	}
	else
	{

		for(j=0; j < num_obs; j++)
		{
		  R_CheckUserInterrupt();

			sum_ker = 0.0;
			sum_ker_convol = 0.0;

			for(i=0; i < num_obs; i++)
			{

				prod_kernel = 1.0;
				prod_kernel_convol = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
					prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
						(matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][i],matrix_bandwidth[l][j],matrix_bandwidth[l][i])
						/matrix_bandwidth[l][i];
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l], matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered], matrix_categorical_vals[l+num_reg_unordered]);
				}

				if(i != j)
				{
					sum_ker += prod_kernel;
				}

				sum_ker_convol += prod_kernel_convol;

			}

			*cv += sum_ker_convol*n_sq_inv-2.0*sum_ker*n_times_n_minus_1_inv;

		}

	}
	#endif

	#ifdef MPI2

	/* Initialize cv function */

	cv_MPI = 0.0;

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		if((num_reg_continuous == 1)&&(num_reg_unordered + num_reg_ordered == 0))
		{

			/* For special case of fixed bandwidth for PDF/CDF for test stat */

			temp_bw1 = matrix_bandwidth[0][0];

			nh_sq_inv = n_sq_inv/temp_bw1;
			nh_times_n_minus_1_inv = n_times_n_minus_1_inv/temp_bw1;

			/* For all i and j */

			p_xj1 = &matrix_X_continuous[0][my_rank*stride];

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				sum_ker = 0.0;
				sum_ker_convol = 0.0;

				/* Move unnecessary pointer out of loop (VTune) */

				temp_xj = *p_xj1++;
				p_xi1 = &matrix_X_continuous[0][0];

				for (i =  0; i <  num_obs; i++)
				{

					temp = (temp_xj - *p_xi1++)/temp_bw1;

					sum_ker_convol += kernel_convol(KERNEL_den,BANDWIDTH_den,temp,temp_bw1,temp_bw1);

					if(i != j)
					{
						sum_ker += kernel(KERNEL_den,temp);
					}

				}

				cv_MPI += sum_ker_convol*nh_sq_inv-2.0*sum_ker*nh_times_n_minus_1_inv;

			}

		}
		else
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				sum_ker = 0.0;
				sum_ker_convol = 0.0;

				for(i=0; i < num_obs; i++)
				{

					prod_kernel = 1.0;
					prod_kernel_convol = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][0])/matrix_bandwidth[l][0];
						prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
							(matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][0],matrix_bandwidth[l][0],matrix_bandwidth[l][0])
							/matrix_bandwidth[l][0];
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
						prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l], matrix_categorical_vals[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
						prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered], matrix_categorical_vals[l+num_reg_unordered]);
					}

					if(i != j)
					{
						sum_ker += prod_kernel;
					}

					sum_ker_convol += prod_kernel_convol;

				}

				cv_MPI += sum_ker_convol*n_sq_inv-2.0*sum_ker*n_times_n_minus_1_inv;

			}

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_convol = 0.0;

			for(i=0; i < num_obs; i++)
			{

				prod_kernel = 1.0;
				prod_kernel_convol = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][j])/matrix_bandwidth[l][j];
					prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
						(matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][j],matrix_bandwidth[l][i],matrix_bandwidth[l][j])
						/matrix_bandwidth[l][j];
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l], matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered], matrix_categorical_vals[l+num_reg_unordered]);
				}

				if(i != j)
				{
					sum_ker += prod_kernel;
				}

				sum_ker_convol += prod_kernel_convol;

			}

			cv_MPI += sum_ker_convol*n_sq_inv-2.0*sum_ker*n_times_n_minus_1_inv;

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_convol = 0.0;

			for(i=0; i < num_obs; i++)
			{

				prod_kernel = 1.0;
				prod_kernel_convol = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel(KERNEL_den, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
					prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
						(matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][i],matrix_bandwidth[l][j],matrix_bandwidth[l][i])
						/matrix_bandwidth[l][i];
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l], matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered],num_categories[l+num_reg_unordered], matrix_categorical_vals[l+num_reg_unordered]);
				}

				if(i != j)
				{
					sum_ker += prod_kernel;
				}

				sum_ker_convol += prod_kernel_convol;

			}

			cv_MPI += sum_ker_convol*n_sq_inv-2.0*sum_ker*n_times_n_minus_1_inv;

		}

	}

	MPI_Reduce(&cv_MPI, cv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(cv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	free(lambda);
	free_mat(matrix_bandwidth,num_reg_continuous);

	return(0);

}


int kernel_estimate_con_density_categorical_convolution_cv(
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
double *cv)
{

	/* This function estimates a leave one out density function using both */
	/* continuous and categorical covariates with three estimation techniques */
	/* and an assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int k;
	int l;

	double prod_kernel_convol;

	double prod_kernel_cat;
	double prod_kernel_cont;

	double prod_kernel_marginal_cat;
	double prod_kernel_marginal_cont;

	double sum_ker_convol;
	double sum_ker_marginal;

	double sum_ker;

	double *lambda;
	double **matrix_bandwidth_var;
	double **matrix_bandwidth_reg;

	double **matrix_weights_K_x;
	double **matrix_weights_K_xy;
	double **matrix_weights_K_convol_y;

	double *pointer_k_xy;
	double *pointer_k_x;
	double *pointer_k_x_kj;
	double *pointer_k_convol_y;

	#ifdef MPI2
	double cv_MPI;
	int stride = ceil((double) num_obs / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered);
	matrix_bandwidth_var = alloc_matd(num_obs,num_var_continuous);
	matrix_bandwidth_reg = alloc_matd(num_obs,num_reg_continuous);

	/* Generate bandwidth vector for continuous dep vars */

	if(kernel_bandwidth_mean(
		KERNEL_den,
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
		lambda)==1)
	{
		free(lambda);
		free_mat(matrix_bandwidth_var,num_var_continuous);
		free_mat(matrix_bandwidth_reg,num_reg_continuous);
		return(1);
	}

	#ifndef MPI2

	/* Initialize cv function */

	*cv = 0.0;

	/* Conduct the estimation */

	if(int_WEIGHTS == 0)
	{

		if(BANDWIDTH_den == 0)
		{

			for(j=0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();
				/* Evaluate at point j */

				sum_ker = 0.0;
				sum_ker_convol = 0.0;
				sum_ker_marginal = DBL_MIN;

				for(i=0; i < num_obs; i++)
				{

					/* Estimate using points i */

					if(i != j)
					{

						/* K(i,j) i != j */

						prod_kernel_cont = 1.0;
						prod_kernel_cat = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][0])/matrix_bandwidth_reg[l][0];
						}

						prod_kernel_marginal_cont = prod_kernel_cont;

						for(l = 0; l < num_var_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];

						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
						}

						prod_kernel_marginal_cat = prod_kernel_cat;

						for(l = 0; l < num_var_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_var_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
						}

						sum_ker += prod_kernel_cont*prod_kernel_cat;
						sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

						/* The dreaded third sum  - fixed j */

						for(k=0; k < num_obs; k++)
						{

							if(k != j)
							{
								/* Initialize to K(i,j), i!=j for X's */

								prod_kernel_convol = prod_kernel_marginal_cont*prod_kernel_marginal_cat;

								/* Multiply by K(k,j), k!=j */

								for(l = 0; l < num_reg_continuous; l++)
								{
									prod_kernel_convol *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][k])/matrix_bandwidth_reg[l][0])/matrix_bandwidth_reg[l][0];
								}

								for(l = 0; l < num_reg_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][k],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
								}

								for(l = 0; l < num_reg_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][k],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
								}

								/* Multiply by K^(2)(i,k) */

								for(l = 0; l < num_var_continuous; l++)
								{
									prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
										(matrix_Y_continuous[l][k]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][0],matrix_bandwidth_var[l][0],matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];
								}

								for(l = 0; l < num_var_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_Y_unordered[l][k],matrix_Y_unordered[l][i],lambda[l], num_categories[l], matrix_categorical_vals[l]);
								}

								for(l = 0; l < num_var_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_Y_ordered[l][k],matrix_Y_ordered[l][i],lambda[l+num_var_unordered], num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
								}

								/* Sum over k,i,j, i != j i != k */

								sum_ker_convol += prod_kernel_convol;

							}

						}

					}

				}

				/* First term has marginal density squared, second does not */

				if(sum_ker_marginal > 0.0)
				{
					*cv += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
				}
				else
				{
					if(int_VERBOSE == 1)
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_con_density_categorical_convolution_cv()");
					}
					*cv += DBL_MAX;
				}

			}

			/* Don't forget!!! */

			*cv /= (double) num_obs;

		}
		else if(BANDWIDTH_den == 1)
		{

			for(j=0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();
				/* Evaluate at point j */

				sum_ker = 0.0;
				sum_ker_convol = 0.0;
				sum_ker_marginal = DBL_MIN;

				for(i=0; i < num_obs; i++)
				{

					/* Estimate using points i */

					if(i != j)
					{

						/* K(i,j) i != j */

						prod_kernel_cont = 1.0;
						prod_kernel_cat = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][j])/matrix_bandwidth_reg[l][j];
						}

						prod_kernel_marginal_cont = prod_kernel_cont;

						for(l = 0; l < num_var_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][j])/matrix_bandwidth_var[l][j];

						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
						}

						prod_kernel_marginal_cat = prod_kernel_cat;

						for(l = 0; l < num_var_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_var_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
						}

						sum_ker += prod_kernel_cont*prod_kernel_cat;
						sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

						/* The dreaded third sum  - fixed j */

						for(k=0; k < num_obs; k++)
						{

							if(k != j)
							{
								/* Initialize to K(i,j), i!=j for X's */

								prod_kernel_convol = prod_kernel_marginal_cont*prod_kernel_marginal_cat;

								/* Multiply by K(k,j), k!=j */

								for(l = 0; l < num_reg_continuous; l++)
								{
									prod_kernel_convol *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][k])/matrix_bandwidth_reg[l][j])/matrix_bandwidth_reg[l][j];
								}

								for(l = 0; l < num_reg_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][k],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
								}

								for(l = 0; l < num_reg_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][k],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
								}

								/* Multiply by K^(2)(i,k) */

								for(l = 0; l < num_var_continuous; l++)
								{
									prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
										(matrix_Y_continuous[l][k]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][k],matrix_bandwidth_var[l][i],matrix_bandwidth_var[l][k])/matrix_bandwidth_var[l][k];
								}

								for(l = 0; l < num_var_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_Y_unordered[l][k],matrix_Y_unordered[l][i],lambda[l], num_categories[l], matrix_categorical_vals[l]);
								}

								for(l = 0; l < num_var_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_Y_ordered[l][k],matrix_Y_ordered[l][i],lambda[l+num_var_unordered], num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
								}

								/* Sum over k,i,j, i != j i != k */

								sum_ker_convol += prod_kernel_convol;

							}

						}

					}

				}

				/* First term has marginal density squared, second does not */

				if(sum_ker_marginal > 0.0)
				{
					*cv += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
				}
				else
				{
					if(int_VERBOSE == 1)
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_con_density_categorical_convolution_cv()");
					}
					*cv += DBL_MAX;
				}

			}

			/* Don't forget!!! */

			*cv /= (double) num_obs;

		}
		else
		{

			for(j=0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();
				/* Evaluate at point j */

				sum_ker = 0.0;
				sum_ker_convol = 0.0;
				sum_ker_marginal = DBL_MIN;

				for(i=0; i < num_obs; i++)
				{

					/* Estimate using points i */

					if(i != j)
					{

						/* K(i,j) i != j */

						prod_kernel_cont = 1.0;
						prod_kernel_cat = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][j])/matrix_bandwidth_reg[l][i];
						}

						prod_kernel_marginal_cont = prod_kernel_cont;

						for(l = 0; l < num_var_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][j])/matrix_bandwidth_var[l][i];

						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
						}

						prod_kernel_marginal_cat = prod_kernel_cat;

						for(l = 0; l < num_var_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_var_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
						}

						sum_ker += prod_kernel_cont*prod_kernel_cat;
						sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

						/* The dreaded third sum  - fixed j */

						for(k=0; k < num_obs; k++)
						{

							if(k != j)
							{
								/* Initialize to K(i,j), i!=j for X's */

								prod_kernel_convol = prod_kernel_marginal_cont*prod_kernel_marginal_cat;

								/* Multiply by K(k,j), k!=j */

								for(l = 0; l < num_reg_continuous; l++)
								{
									prod_kernel_convol *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][k])/matrix_bandwidth_reg[l][j])/matrix_bandwidth_reg[l][k];
								}

								for(l = 0; l < num_reg_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][k],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
								}

								for(l = 0; l < num_reg_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][k],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
								}

								/* Multiply by K^(2)(i,k) */

								for(l = 0; l < num_var_continuous; l++)
								{
									prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
										(matrix_Y_continuous[l][k]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][i],matrix_bandwidth_var[l][k],matrix_bandwidth_var[l][i])/matrix_bandwidth_var[l][i];
								}

								for(l = 0; l < num_var_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_Y_unordered[l][k],matrix_Y_unordered[l][i],lambda[l], num_categories[l], matrix_categorical_vals[l]);
								}

								for(l = 0; l < num_var_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_Y_ordered[l][k],matrix_Y_ordered[l][i],lambda[l+num_var_unordered], num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
								}

								/* Sum over k,i,j, i != j i != k */

								sum_ker_convol += prod_kernel_convol;

							}

						}

					}

				}

				/* First term has marginal density squared, second does not */

				if(sum_ker_marginal > 0.0)
				{
					*cv += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
				}
				else
				{
					if(int_VERBOSE == 1)
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_con_density_categorical_convolution_cv()");
					}
					*cv += DBL_MAX;
				}

			}

			*cv /= (double) num_obs;

		}

	}
	else
	{

		/* Allocate memory */

		matrix_weights_K_x = alloc_matd(num_obs,num_obs);
		matrix_weights_K_xy = alloc_matd(num_obs,num_obs);
		matrix_weights_K_convol_y = alloc_matd(num_obs,num_obs);

		/* Compute kernel weights */

		kernel_weights_conditional_convolution_cv(
			int_WEIGHTS,
			KERNEL_den,
			KERNEL_unordered_den,
			KERNEL_ordered_den,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_den,
			num_obs,
			num_var_unordered,
			num_var_ordered,
			num_var_continuous,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_Y_unordered,
			matrix_Y_ordered,
			matrix_Y_continuous,
			matrix_X_unordered,
			matrix_X_ordered,
			matrix_X_continuous,
			lambda,
			matrix_bandwidth_var,
			matrix_bandwidth_reg,
			num_categories,
			matrix_categorical_vals,
			matrix_weights_K_x,
			matrix_weights_K_xy,
			matrix_weights_K_convol_y);

		for(j=0; j < num_obs; j++)
		{
		  R_CheckUserInterrupt();
			/* Evaluate at point j */

			sum_ker = 0.0;
			sum_ker_convol = 0.0;
			sum_ker_marginal = DBL_MIN;

			pointer_k_xy = &matrix_weights_K_xy[j][0];
			pointer_k_x = &matrix_weights_K_x[j][0];

			for(i=0; i < num_obs; i++)
			{

				/* Estimate using points i */

				if(i != j)
				{

					/* K(i,j) i != j */

					sum_ker += *pointer_k_xy;
					sum_ker_marginal += *pointer_k_x;

					/* The dreaded third sum  - fixed j */

					pointer_k_x_kj = &matrix_weights_K_x[j][0];
					pointer_k_convol_y = &matrix_weights_K_convol_y[i][0];

					for(k=0; k < num_obs; k++)
					{

						if(k != j)
						{
							/* K(i,j)*K(k,j)*K^(2)(i,k) */
							sum_ker_convol += *pointer_k_x_kj * *pointer_k_x * *pointer_k_convol_y;
						}

						pointer_k_x_kj++;
						pointer_k_convol_y++;

					}

				}

				pointer_k_xy++;
				pointer_k_x++;

			}

			/* First term has marginal density squared, second does not */

			if(sum_ker_marginal > 0.0)
			{
				*cv += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
			}
			else
			{
				if(int_VERBOSE == 1)
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_convolution_cv()");
				}
				*cv += DBL_MAX;
			}

		}

		/* Don't forget!!! */

		*cv /= (double) num_obs;

		free_mat(matrix_weights_K_x,num_obs);
		free_mat(matrix_weights_K_xy,num_obs);
		free_mat(matrix_weights_K_convol_y,num_obs);

	}
	#endif

	#ifdef MPI2

	/* Initialize cv function */

	cv_MPI = 0.0;

	/* Conduct the estimation */

	if(int_WEIGHTS == 0)
	{

		if(BANDWIDTH_den == 0)
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				/* Evaluate at point j */

				sum_ker = 0.0;
				sum_ker_convol = 0.0;
				sum_ker_marginal = DBL_MIN;

				for(i=0; i < num_obs; i++)
				{

					/* Estimate using points i */

					if(i != j)
					{

						/* K(i,j) i != j */

						prod_kernel_cont = 1.0;
						prod_kernel_cat = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][0])/matrix_bandwidth_reg[l][0];
						}

						prod_kernel_marginal_cont = prod_kernel_cont;

						for(l = 0; l < num_var_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];

						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
						}

						prod_kernel_marginal_cat = prod_kernel_cat;

						for(l = 0; l < num_var_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_var_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
						}

						sum_ker += prod_kernel_cont*prod_kernel_cat;
						sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

						/* The dreaded third sum  - fixed j */

						for(k=0; k < num_obs; k++)
						{

							if(k != j)
							{
								/* Initialize to K(i,j), i!=j for X's */

								prod_kernel_convol = prod_kernel_marginal_cont*prod_kernel_marginal_cat;

								/* Multiply by K(k,j), k!=j */

								for(l = 0; l < num_reg_continuous; l++)
								{
									prod_kernel_convol *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][k])/matrix_bandwidth_reg[l][0])/matrix_bandwidth_reg[l][0];
								}

								for(l = 0; l < num_reg_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][k],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
								}

								for(l = 0; l < num_reg_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][k],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
								}

								/* Multiply by K^(2)(i,k) */

								for(l = 0; l < num_var_continuous; l++)
								{
									prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
										(matrix_Y_continuous[l][k]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][0],matrix_bandwidth_var[l][0],matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];
								}

								for(l = 0; l < num_var_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_Y_unordered[l][k],matrix_Y_unordered[l][i],lambda[l], num_categories[l], matrix_categorical_vals[l]);
								}

								for(l = 0; l < num_var_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_Y_ordered[l][k],matrix_Y_ordered[l][i],lambda[l+num_var_unordered], num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
								}

								/* Sum over k,i,j, i != j i != k */

								sum_ker_convol += prod_kernel_convol;

							}

						}

					}

				}

				/* First term has marginal density squared, second does not */

				if(sum_ker_marginal > 0.0)
				{
					cv_MPI += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
				}
				else
				{
					if((int_VERBOSE == 1)&&(my_rank == 0))
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_con_density_categorical_convolution_cv()");
					}
					cv_MPI += DBL_MAX;
				}

			}

			/* Don't forget!!! */

			cv_MPI /= (double) num_obs;

		}
		else if(BANDWIDTH_den == 1)
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				/* Evaluate at point j */

				sum_ker = 0.0;
				sum_ker_convol = 0.0;
				sum_ker_marginal = DBL_MIN;

				for(i=0; i < num_obs; i++)
				{

					/* Estimate using points i */

					if(i != j)
					{

						/* K(i,j) i != j */

						prod_kernel_cont = 1.0;
						prod_kernel_cat = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][j])/matrix_bandwidth_reg[l][j];
						}

						prod_kernel_marginal_cont = prod_kernel_cont;

						for(l = 0; l < num_var_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][j])/matrix_bandwidth_var[l][j];

						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
						}

						prod_kernel_marginal_cat = prod_kernel_cat;

						for(l = 0; l < num_var_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_var_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
						}

						sum_ker += prod_kernel_cont*prod_kernel_cat;
						sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

						/* The dreaded third sum  - fixed j */

						for(k=0; k < num_obs; k++)
						{

							if(k != j)
							{
								/* Initialize to K(i,j), i!=j for X's */

								prod_kernel_convol = prod_kernel_marginal_cont*prod_kernel_marginal_cat;

								/* Multiply by K(k,j), k!=j */

								for(l = 0; l < num_reg_continuous; l++)
								{
									prod_kernel_convol *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][k])/matrix_bandwidth_reg[l][j])/matrix_bandwidth_reg[l][j];
								}

								for(l = 0; l < num_reg_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][k],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
								}

								for(l = 0; l < num_reg_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][k],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
								}

								/* Multiply by K^(2)(i,k) */

								for(l = 0; l < num_var_continuous; l++)
								{
									prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
										(matrix_Y_continuous[l][k]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][k],matrix_bandwidth_var[l][i],matrix_bandwidth_var[l][k])/matrix_bandwidth_var[l][k];
								}

								for(l = 0; l < num_var_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_Y_unordered[l][k],matrix_Y_unordered[l][i],lambda[l], num_categories[l], matrix_categorical_vals[l]);
								}

								for(l = 0; l < num_var_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_Y_ordered[l][k],matrix_Y_ordered[l][i],lambda[l+num_var_unordered], num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
								}

								/* Sum over k,i,j, i != j i != k */

								sum_ker_convol += prod_kernel_convol;

							}

						}

					}

				}

				/* First term has marginal density squared, second does not */

				if(sum_ker_marginal > 0.0)
				{
					cv_MPI += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
				}
				else
				{
					if((int_VERBOSE == 1)&&(my_rank == 0))
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_con_density_categorical_convolution_cv()");
					}
					cv_MPI += DBL_MAX;
				}

			}

			/* Don't forget!!! */

			cv_MPI /= (double) num_obs;

		}
		else
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				/* Evaluate at point j */

				sum_ker = 0.0;
				sum_ker_convol = 0.0;
				sum_ker_marginal = DBL_MIN;

				for(i=0; i < num_obs; i++)
				{

					/* Estimate using points i */

					if(i != j)
					{

						/* K(i,j) i != j */

						prod_kernel_cont = 1.0;
						prod_kernel_cat = 1.0;

						for(l = 0; l < num_reg_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][j])/matrix_bandwidth_reg[l][i];
						}

						prod_kernel_marginal_cont = prod_kernel_cont;

						for(l = 0; l < num_var_continuous; l++)
						{
							prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][j])/matrix_bandwidth_var[l][i];

						}

						for(l = 0; l < num_reg_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
						}

						for(l = 0; l < num_reg_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
						}

						prod_kernel_marginal_cat = prod_kernel_cat;

						for(l = 0; l < num_var_unordered; l++)
						{
							prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
						}

						for(l = 0; l < num_var_ordered; l++)
						{
							prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
						}

						sum_ker += prod_kernel_cont*prod_kernel_cat;
						sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;

						/* The dreaded third sum  - fixed j */

						for(k=0; k < num_obs; k++)
						{

							if(k != j)
							{
								/* Initialize to K(i,j), i!=j for X's */

								prod_kernel_convol = prod_kernel_marginal_cont*prod_kernel_marginal_cat;

								/* Multiply by K(k,j), k!=j */

								for(l = 0; l < num_reg_continuous; l++)
								{
									prod_kernel_convol *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][k])/matrix_bandwidth_reg[l][j])/matrix_bandwidth_reg[l][k];
								}

								for(l = 0; l < num_reg_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][k],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
								}

								for(l = 0; l < num_reg_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][k],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
								}

								/* Multiply by K^(2)(i,k) */

								for(l = 0; l < num_var_continuous; l++)
								{
									prod_kernel_convol *= kernel_convol(KERNEL_den,BANDWIDTH_den,
										(matrix_Y_continuous[l][k]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][i],matrix_bandwidth_var[l][k],matrix_bandwidth_var[l][i])/matrix_bandwidth_var[l][i];
								}

								for(l = 0; l < num_var_unordered; l++)
								{
									prod_kernel_convol *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_Y_unordered[l][k],matrix_Y_unordered[l][i],lambda[l], num_categories[l], matrix_categorical_vals[l]);
								}

								for(l = 0; l < num_var_ordered; l++)
								{
									prod_kernel_convol *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_Y_ordered[l][k],matrix_Y_ordered[l][i],lambda[l+num_var_unordered], num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
								}

								/* Sum over k,i,j, i != j i != k */

								sum_ker_convol += prod_kernel_convol;

							}

						}

					}

				}

				/* First term has marginal density squared, second does not */

				if(sum_ker_marginal > 0.0)
				{
					cv_MPI += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
				}
				else
				{
					if((int_VERBOSE == 1)&&(my_rank == 0))
					{
						printf("\r                                                                           ");
						printf("\r** Trimming binding in kernel_estimate_con_density_categorical_convolution_cv()");
					}
					cv_MPI += DBL_MAX;
				}

			}

			cv_MPI /= (double) num_obs;

		}

	}
	else
	{

		/* Allocate memory */

		matrix_weights_K_x = alloc_matd(stride*iNum_Processors,stride*iNum_Processors);
		matrix_weights_K_xy = alloc_matd(stride*iNum_Processors,stride*iNum_Processors);
		matrix_weights_K_convol_y = alloc_matd(stride*iNum_Processors,stride*iNum_Processors);

		/* Compute kernel weights */

		kernel_weights_conditional_convolution_cv(
			int_WEIGHTS,
			KERNEL_den,
			KERNEL_unordered_den,
			KERNEL_ordered_den,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_den,
			num_obs,
			num_var_unordered,
			num_var_ordered,
			num_var_continuous,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_Y_unordered,
			matrix_Y_ordered,
			matrix_Y_continuous,
			matrix_X_unordered,
			matrix_X_ordered,
			matrix_X_continuous,
			lambda,
			matrix_bandwidth_var,
			matrix_bandwidth_reg,
			num_categories,
			matrix_categorical_vals,
			matrix_weights_K_x,
			matrix_weights_K_xy,
			matrix_weights_K_convol_y);

		for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
		{

			/* Evaluate at point j */

			sum_ker = 0.0;
			sum_ker_convol = 0.0;
			sum_ker_marginal = DBL_MIN;

			pointer_k_xy = &matrix_weights_K_xy[j][0];
			pointer_k_x = &matrix_weights_K_x[j][0];

			for(i=0; i < num_obs; i++)
			{

				/* Estimate using points i */

				if(i != j)
				{

					/* K(i,j) i != j */

					sum_ker += *pointer_k_xy;
					sum_ker_marginal += *pointer_k_x;

					/* The dreaded third sum  - fixed j */

					pointer_k_x_kj = &matrix_weights_K_x[j][0];
					pointer_k_convol_y = &matrix_weights_K_convol_y[i][0];

					for(k=0; k < num_obs; k++)
					{

						if(k != j)
						{
							/* K(i,j)*K(k,j)*K^(2)(i,k) */
							sum_ker_convol += *pointer_k_x_kj * *pointer_k_x * *pointer_k_convol_y;
						}

						pointer_k_x_kj++;
						pointer_k_convol_y++;

					}

				}

				pointer_k_xy++;
				pointer_k_x++;

			}

			/* First term has marginal density squared, second does not */

			if(sum_ker_marginal > 0.0)
			{
				cv_MPI += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
			}
			else
			{
				if((int_VERBOSE == 1)&&(my_rank == 0))
				{
					printf("\r                                                                           ");
					printf("\r** Trimming binding in kernel_estimate_con_density_categorical_convolution_cv()");
				}
				cv_MPI += DBL_MAX;
			}

		}

		/* Don't forget!!! */

		cv_MPI /= (double) num_obs;

		free_mat(matrix_weights_K_x,stride*iNum_Processors);
		free_mat(matrix_weights_K_xy,stride*iNum_Processors);
		free_mat(matrix_weights_K_convol_y,stride*iNum_Processors);

	}
	MPI_Reduce(&cv_MPI, cv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(cv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	free(lambda);
	free_mat(matrix_bandwidth_var,num_var_continuous);
	free_mat(matrix_bandwidth_reg,num_reg_continuous);

	return(0);

}


int kernel_estimate_quantile(
int gradient_compute,
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int BANDWIDTH_den,
int num_obs_train,
int num_obs_eval,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_Y_unordered_train,
double **matrix_Y_ordered_train,
double **matrix_Y_continuous_train,
double **matrix_Y_unordered_eval,
double **matrix_Y_ordered_eval,
double **matrix_Y_continuous_eval,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_scale_factor,
double *quan,
double *quan_stderr,
double **quan_gradient,
int seed,
double ftol,
double tol,
double small,
int itmax,
int iMax_Num_Multistart,
double zero)
{

	int i;
	int j;
	int k;
	double fret;
	double fret_best;
	double quantile[2];
	double quantile_multistart[2];
	int iImproved;
	int iMs_counter;
	int iter;
	int iNum_Ms;
	double **matrix_y;

	double quantile_l;
	double quantile_u;

	double *lambda = NULL;
	double **matrix_bandwidth_var = NULL;
	double **matrix_bandwidth_reg = NULL;

	#ifdef MPI2
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	if(gradient_compute == 1)
	{

		lambda = alloc_vecd(num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered);

		if((BANDWIDTH_den == 0)||(BANDWIDTH_den == 1))
		{
			matrix_bandwidth_var = alloc_matd(num_obs_eval,num_var_continuous);
			matrix_bandwidth_reg = alloc_matd(num_obs_eval,num_reg_continuous);
		}
		else if(BANDWIDTH_den == 2)
		{
			matrix_bandwidth_var = alloc_matd(num_obs_train,num_var_continuous);
			matrix_bandwidth_reg = alloc_matd(num_obs_train,num_reg_continuous);
		}

		/* Bandwidths for `dependent' variables */

		if(kernel_bandwidth_mean(
			KERNEL_den,
			BANDWIDTH_den,
			num_obs_train,
			num_obs_eval,
			num_var_continuous,
			num_var_unordered,
			num_var_ordered,
			num_reg_continuous,
			num_reg_unordered,
			num_reg_ordered,
			vector_scale_factor,
			matrix_Y_continuous_train,
			matrix_Y_continuous_train, /* Same Y for training and evaluation */
			matrix_X_continuous_train,
			matrix_X_continuous_eval,
			matrix_bandwidth_var,
			matrix_bandwidth_reg,
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

	}

	#ifndef MPI2

	matrix_y = alloc_matd(2,2);

	y_min_extern = y_max_extern = matrix_Y_continuous_train[0][0];

	/* Was zero for erfun() */

	itmax_extern = itmax;
	small_extern = small;

	for(i=0; i < num_obs_train; i++)
	{
		if(matrix_Y_continuous_train[0][i] < y_min_extern)
		{
			y_min_extern = matrix_Y_continuous_train[0][i];
		}
		if(matrix_Y_continuous_train[0][i] > y_max_extern)
		{
			y_max_extern = matrix_Y_continuous_train[0][i];
		}
	}

	for(i=0; i < num_obs_eval; i++)
	{
	  R_CheckUserInterrupt();
		for(j = 0; j < num_reg_unordered; j++)
		{
			matrix_X_unordered_quantile_extern[j][0] = matrix_X_unordered_eval[j][i];
		}

		for(j = 0; j < num_reg_ordered; j++)
		{
			matrix_X_ordered_quantile_extern[j][0] = matrix_X_ordered_eval[j][i];
		}

		for(j = 0; j < num_reg_continuous; j++)
		{
			matrix_X_continuous_quantile_extern[j][0] = matrix_X_continuous_eval[j][i];
		}

		quantile[1] = (y_max_extern-y_min_extern)/2.0;

		initialize_nr_hessian(1, matrix_y);

		powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

		if(fret > zero)
		{
			fret_best = fret;
			quantile_multistart[1] = quantile[1];

			for(iMs_counter = 1, iNum_Ms = 1; iMs_counter < iMax_Num_Multistart; iMs_counter++, iNum_Ms++)
			{

				quantile[1] = y_min_extern + ran3(&seed)*(y_max_extern-y_min_extern);

				initialize_nr_hessian(1, matrix_y);

				powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

				if(fret < fret_best)
				{
					fret_best = fret;
					iImproved = 1;
					quantile_multistart[1] = quantile[1];
					if(fret <= zero)
					{
						iMs_counter = iMax_Num_Multistart;
					}
				}
				else
				{
					iImproved = 0;
				}

			}

			fret = fret_best;
			quantile[1] = quantile_multistart[1];

      if(int_MINIMIZE_IO != IO_MIN_TRUE){
        printf("\r                                                                             ");
        printf("\rWorking... (observation %d/%d required %d restarts to attain %g)", i+1, num_obs_eval, iNum_Ms, zero);
        fflush(stdout);
      }

		}

		quan[i] = quantile[1];
		quan_stderr[i] = 0.0;

		/* Need to correct for bw0 */

		if(gradient_compute == 1)
		{

			for(k = 0; k < num_reg_continuous; k++)
			{

				/* Gradient for continuous regressors - quantile evaluated at x-h */

				for(j = 0; j < num_reg_continuous; j++)
				{
					matrix_X_continuous_quantile_extern[j][0] = matrix_X_continuous_eval[j][i];
				}

				matrix_X_continuous_quantile_extern[k][0] = matrix_X_continuous_eval[k][i]  - matrix_bandwidth_reg[k][i];

				quantile[1] = (y_max_extern-y_min_extern)/2.0;

				initialize_nr_hessian(1, matrix_y);

				powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

				if(fret > zero)
				{
					fret_best = fret;
					quantile_multistart[1] = quantile[1];

					for(iMs_counter = 1, iNum_Ms = 1; iMs_counter < iMax_Num_Multistart; iMs_counter++, iNum_Ms++)
					{

						quantile[1] = y_min_extern + ran3(&seed)*(y_max_extern-y_min_extern);

						initialize_nr_hessian(1, matrix_y);

						powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

						if(fret < fret_best)
						{
							fret_best = fret;
							iImproved = 1;
							quantile_multistart[1] = quantile[1];
							if(fret <= zero)
							{
								iMs_counter = iMax_Num_Multistart;
							}
						}
						else
						{
							iImproved = 0;
						}

					}

					fret = fret_best;
					quantile[1] = quantile_multistart[1];
          if(int_MINIMIZE_IO != IO_MIN_TRUE){
            printf("\r                                                                             ");
            printf("\rWorking... (observation %d/%d required %d restarts to attain %g)", i+1, num_obs_eval, iNum_Ms, zero);
            fflush(stdout);
          }
				}

				quantile_l = quantile[1];

				/* Gradient for continuous regressors - quantile evaluated at x+h */

				for(j = 0; j < num_reg_continuous; j++)
				{
					matrix_X_continuous_quantile_extern[j][0] = matrix_X_continuous_eval[j][i];
				}

				matrix_X_continuous_quantile_extern[k][0] = matrix_X_continuous_eval[k][i] + matrix_bandwidth_reg[k][i]/2.0;

				quantile[1] = (y_max_extern-y_min_extern)/2.0;

				initialize_nr_hessian(1, matrix_y);

				powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

				if(fret > zero)
				{
					fret_best = fret;
					quantile_multistart[1] = quantile[1];

					for(iMs_counter = 1, iNum_Ms = 1; iMs_counter < iMax_Num_Multistart; iMs_counter++, iNum_Ms++)
					{

						quantile[1] = y_min_extern + ran3(&seed)*(y_max_extern-y_min_extern);

						initialize_nr_hessian(1, matrix_y);

						powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

						if(fret < fret_best)
						{
							fret_best = fret;
							iImproved = 1;
							quantile_multistart[1] = quantile[1];
							if(fret <= zero)
							{
								iMs_counter = iMax_Num_Multistart;
							}
						}
						else
						{
							iImproved = 0;
						}

					}

					fret = fret_best;
					quantile[1] = quantile_multistart[1];
          if(int_MINIMIZE_IO != IO_MIN_TRUE){
            printf("\r                                                                             ");
            printf("\rWorking... (observation %d/%d required %d restarts to attain %g)", i+1, num_obs_eval, iNum_Ms, zero);
            fflush(stdout);
          }

				}

				quantile_u = quantile[1];

				quan_gradient[k][i] = (quantile_u-quantile_l)/(2.0*matrix_bandwidth_reg[k][i]);

			}

		}

		/* End gradient */

	}

	free_mat(matrix_y, 2);
	#endif

	#ifdef MPI2

	matrix_y = alloc_matd(2,2);

	y_min_extern = y_max_extern = matrix_Y_continuous_train[0][0];

	/* Was zero for erfun() */

	itmax_extern = itmax;
	small_extern = small;

	for(i=0; i < num_obs_train; i++)
	{
		if(matrix_Y_continuous_train[0][i] < y_min_extern)
		{
			y_min_extern = matrix_Y_continuous_train[0][i];
		}
		if(matrix_Y_continuous_train[0][i] > y_max_extern)
		{
			y_max_extern = matrix_Y_continuous_train[0][i];
		}
	}

	/* Converting to  MPI aware 11/16/04 */

	for(i=my_rank*stride; (i < num_obs_eval) && (i < (my_rank+1)*stride); i++)
	{

		for(j = 0; j < num_reg_unordered; j++)
		{
			matrix_X_unordered_quantile_extern[j][0] = matrix_X_unordered_eval[j][i];
		}

		for(j = 0; j < num_reg_ordered; j++)
		{
			matrix_X_ordered_quantile_extern[j][0] = matrix_X_ordered_eval[j][i];
		}

		for(j = 0; j < num_reg_continuous; j++)
		{
			matrix_X_continuous_quantile_extern[j][0] = matrix_X_continuous_eval[j][i];
		}

		quantile[1] = (y_max_extern-y_min_extern)/2.0;

		initialize_nr_hessian(1, matrix_y);

		powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

		if(fret > zero)
		{
			fret_best = fret;
			quantile_multistart[1] = quantile[1];

			for(iMs_counter = 1, iNum_Ms = 1; iMs_counter < iMax_Num_Multistart; iMs_counter++, iNum_Ms++)
			{

				quantile[1] = y_min_extern + ran3(&seed)*(y_max_extern-y_min_extern);

				initialize_nr_hessian(1, matrix_y);

				powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

				if(fret < fret_best)
				{
					fret_best = fret;
					iImproved = 1;
					quantile_multistart[1] = quantile[1];
					if(fret <= zero)
					{
						iMs_counter = iMax_Num_Multistart;
					}
				}
				else
				{
					iImproved = 0;
				}

			}

			fret = fret_best;
			quantile[1] = quantile_multistart[1];

			if((my_rank == 0) && (int_MINIMIZE_IO != IO_MIN_TRUE))
			{
				printf("\r                                                                             ");
				printf("\rWorking... (observation %d/%d required %d restarts to attain %g)", i+1, stride, iNum_Ms, zero);
				fflush(stdout);

			}

		}

		quan[i-my_rank*stride] = quantile[1];
		quan_stderr[i-my_rank*stride] = 0.0;

		if(gradient_compute == 1)
		{

			for(k = 0; k < num_reg_continuous; k++)
			{

				/* Gradient for continuous regressors - quantile evaluated at x-h */

				for(j = 0; j < num_reg_continuous; j++)
				{
					matrix_X_continuous_quantile_extern[j][0] = matrix_X_continuous_eval[j][i];
				}

				matrix_X_continuous_quantile_extern[k][0] = matrix_X_continuous_eval[k][i]  - matrix_bandwidth_reg[k][i];

				quantile[1] = (y_max_extern-y_min_extern)/2.0;

				initialize_nr_hessian(1, matrix_y);

				powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

				if(fret > zero)
				{
					fret_best = fret;
					quantile_multistart[1] = quantile[1];

					for(iMs_counter = 1, iNum_Ms = 1; iMs_counter < iMax_Num_Multistart; iMs_counter++, iNum_Ms++)
					{

						quantile[1] = y_min_extern + ran3(&seed)*(y_max_extern-y_min_extern);

						initialize_nr_hessian(1, matrix_y);

						powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

						if(fret < fret_best)
						{
							fret_best = fret;
							iImproved = 1;
							quantile_multistart[1] = quantile[1];
							if(fret <= zero)
							{
								iMs_counter = iMax_Num_Multistart;
							}
						}
						else
						{
							iImproved = 0;
						}

					}

					fret = fret_best;
					quantile[1] = quantile_multistart[1];

          if((my_rank == 0) && (int_MINIMIZE_IO != IO_MIN_TRUE))
					{

						printf("\r                                                                             ");
						printf("\rWorking... (observation %d/%d required %d restarts to attain %g)", i+1, num_obs_eval, iNum_Ms, zero);
						fflush(stdout);

					}

				}

				quantile_l = quantile[1];

				/* Gradient for continuous regressors - quantile evaluated at x+h */

				for(j = 0; j < num_reg_continuous; j++)
				{
					matrix_X_continuous_quantile_extern[j][0] = matrix_X_continuous_eval[j][i];
				}

				matrix_X_continuous_quantile_extern[k][0] = matrix_X_continuous_eval[k][i] + matrix_bandwidth_reg[k][i]/2.0;

				quantile[1] = (y_max_extern-y_min_extern)/2.0;

				initialize_nr_hessian(1, matrix_y);

				powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

				if(fret > zero)
				{
					fret_best = fret;
					quantile_multistart[1] = quantile[1];

					for(iMs_counter = 1, iNum_Ms = 1; iMs_counter < iMax_Num_Multistart; iMs_counter++, iNum_Ms++)
					{

						quantile[1] = y_min_extern + ran3(&seed)*(y_max_extern-y_min_extern);

						initialize_nr_hessian(1, matrix_y);

						powell(0, 0, quantile, quantile, matrix_y, 1, ftol, tol, small, itmax, &iter, &fret, func_con_density_quantile);

						if(fret < fret_best)
						{
							fret_best = fret;
							iImproved = 1;
							quantile_multistart[1] = quantile[1];
							if(fret <= zero)
							{
								iMs_counter = iMax_Num_Multistart;
							}
						}
						else
						{
							iImproved = 0;
						}

					}

					fret = fret_best;
					quantile[1] = quantile_multistart[1];

          if((my_rank == 0) && (int_MINIMIZE_IO != IO_MIN_TRUE))
					{
						printf("\r                                                                             ");
						printf("\rWorking... (observation %d/%d required %d restarts to attain %g)", i+1, num_obs_eval, iNum_Ms, zero);
						fflush(stdout);

					}

				}

				quantile_u = quantile[1];

				quan_gradient[k][i-my_rank*stride] = (quantile_u-quantile_l)/(2.0*matrix_bandwidth_reg[k][i]);

			}

		}

		/* End gradient */

	}

	/* Collect */

	MPI_Gather(quan, stride, MPI_DOUBLE, quan, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(quan, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Gather(quan_stderr, stride, MPI_DOUBLE, quan_stderr, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(quan_stderr, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(gradient_compute == 1)
	{

		for(k = 0; k < num_reg_continuous; k++)
		{

			MPI_Gather(&quan_gradient[k][0], stride, MPI_DOUBLE, &quan_gradient[k][0], stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&quan_gradient[k][0], num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		}

	}

	free_mat(matrix_y, 2);
	#endif

	if(gradient_compute == 1)
	{
		free(lambda);
		free_mat(matrix_bandwidth_var,num_var_continuous);
		free_mat(matrix_bandwidth_reg,num_reg_continuous);
	}

	return(0);

}


double kernel_estimate_regression_categorical_aic_c(
int int_ll,
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
int *num_categories)
{

	/* This function estimates a leave-one-out Nadaraya-Watson regression */
	/* function using both continuous and categorical covariates with three */
	/* estimation techniques and an assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int k;
	int l;

	const double epsilon = 1.0/num_obs;
  double nepsilon;

	double prod_kernel;

	double sum_ker;
	double sum_y_ker;

	double *lambda;
	double **matrix_bandwidth;

	double temp;
	double temp1 = DBL_MAX;

	MATRIX  XTKX;
	MATRIX  XTKXINV;
	MATRIX  XTKY;
	MATRIX  DELTA;

	double *pointer_yi;
	double *pointer_m;

	int num_reg_cat_cont;

	double *mean;
	double prod_kernel_i_eq_j = DBL_MAX;
	double trace_H = 0.0;
	double aic_c = 0.0;
	double sigmasq = 0.0;

	#ifdef MPI2
	double trace_H_MPI = 0.0;
	int stride = ceil((double) num_obs / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	mean = alloc_vecd(stride*iNum_Processors);
	#endif
	#ifndef MPI2
	mean = alloc_vecd(num_obs);
	#endif

	if(int_TAYLOR == 1)
	{
		num_reg_cat_cont = num_reg_unordered + num_reg_ordered + num_reg_continuous;
	}
	else
	{
		num_reg_cat_cont = num_reg_continuous;
	}

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);
	matrix_bandwidth = alloc_matd(num_obs,num_reg_continuous);

	#ifndef MPI2

	/* Conduct the estimation */

	if(int_ll == 0)
	{

		/* Nadaraya-Watson */

		/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

		if(kernel_bandwidth_mean(
			KERNEL_reg,
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
			matrix_X_continuous,			 /* Not used */
			matrix_X_continuous,			 /* Not used */
			matrix_X_continuous,
			matrix_X_continuous,
			matrix_bandwidth,					 /* Not used */
			matrix_bandwidth,
			lambda)==1)
		{

			free(lambda);
			free_mat(matrix_bandwidth,num_reg_continuous);

			return(DBL_MAX);
		}

		if(BANDWIDTH_reg == 0)
		{

			pointer_m = &mean[0];

			for(j=0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();
				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][0]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_y_ker += *pointer_yi*prod_kernel;

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					*pointer_m++ = sum_y_ker/sum_ker;
					trace_H += prod_kernel_i_eq_j/sum_ker;
				}
				else
				{
					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_aic()",j);
					}
					return(DBL_MAX);
				}

			}

		}
		else if(BANDWIDTH_reg == 1)
		{

			pointer_m = &mean[0];

			for(j=0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();
				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][j]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_y_ker += *pointer_yi*prod_kernel;

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					*pointer_m++ = sum_y_ker/sum_ker;
					trace_H += prod_kernel_i_eq_j/sum_ker;
				}
				else
				{
					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_aic()",j);
					}
					return(DBL_MAX);
				}

			}

		}
		else
		{

			pointer_m = &mean[0];

			for(j=0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();
				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_y_ker += *pointer_yi*prod_kernel;

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					*pointer_m++ = sum_y_ker/sum_ker;
					trace_H += prod_kernel_i_eq_j/sum_ker;
				}
				else
				{
					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_aic()",j);
					}
					return(DBL_MAX);
				}

			}

		}

	}
	else
	{

		/* Local Linear */

		XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
		DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

		/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

		if(kernel_bandwidth_mean(
			KERNEL_reg,
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
			matrix_X_continuous,			 /* Not used */
			matrix_X_continuous,			 /* Not used */
			matrix_X_continuous,
			matrix_X_continuous,
			matrix_bandwidth,					 /* Not used */
			matrix_bandwidth,
			lambda)==1)
		{

			mat_free( XTKX );
			mat_free( XTKXINV );
			mat_free( XTKY );
			mat_free( DELTA );

			free(lambda);
			free_mat(matrix_bandwidth,num_reg_continuous);

			return(DBL_MAX);

		}

		/* Conduct the estimation */

		if(BANDWIDTH_reg == 0)
		{

			pointer_m = &mean[0];

			for (j = 0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();
				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][0]);
					}

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
					}

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_aic()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
				trace_H += XTKXINV[0][0]*prod_kernel_i_eq_j;
				*pointer_m++ = DELTA[0][0];

			}

		}
		else if(BANDWIDTH_reg == 1)
		{

			pointer_m = &mean[0];

			for (j = 0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();
				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][j]);
					}

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
					}

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_aic()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
				trace_H += XTKXINV[0][0]*prod_kernel_i_eq_j;
				*pointer_m++ = DELTA[0][0];

			}

		}
		else
		{

			pointer_m = &mean[0];

			for (j = 0; j < num_obs; j++)
			{
			  R_CheckUserInterrupt();
				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][i]);
					}

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
					}

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if(int_DEBUG == 1)
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_aic()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
				trace_H += XTKXINV[0][0]*prod_kernel_i_eq_j;
				*pointer_m++ = DELTA[0][0];

			}

		}

		mat_free( XTKX );
		mat_free( XTKXINV );
		mat_free( XTKY );
		mat_free( DELTA );

	}
	#endif

	#ifdef MPI2

	/* Conduct the estimation - MPI-enabled */

	if(int_ll == 0)
	{

		/* Nadaraya-Watson */

		/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

		/* May be redundant */

		if(kernel_bandwidth_mean(
			KERNEL_reg,
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
			/* Not used */
			matrix_X_continuous,
			/* Not used */
			matrix_X_continuous,
			matrix_X_continuous,
			matrix_X_continuous,
			/* Not used */
			matrix_bandwidth,
			matrix_bandwidth,
			lambda)==1)
		{

			free(lambda);
			free_mat(matrix_bandwidth,num_reg_continuous);

			return(DBL_MAX);
		}

		if(BANDWIDTH_reg == 0)
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][0]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_y_ker += *pointer_yi*prod_kernel;

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					mean[j-my_rank*stride] = sum_y_ker/sum_ker;
					trace_H_MPI += prod_kernel_i_eq_j/sum_ker;
				}
				else
				{
					/* Don't print if you are using MPI */
					if((int_DEBUG == 1)&&(my_rank == 0))
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_aic()",j);
					}
					return(DBL_MAX);
				}

			}

		}
		else if(BANDWIDTH_reg == 1)
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][j]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_y_ker += *pointer_yi*prod_kernel;

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					mean[j-my_rank*stride] = sum_y_ker/sum_ker;
					trace_H_MPI += prod_kernel_i_eq_j/sum_ker;
				}
				else
				{
					/* Don't print if you are using MPI */
					if((my_rank == 0)&&(int_DEBUG == 1))
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_aic()",j);
					}
					return(DBL_MAX);
				}

			}

		}
		else
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{
				sum_y_ker = 0.0;
				sum_ker = DBL_MIN;

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth[l][i])/matrix_bandwidth[l][i];
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_y_ker += *pointer_yi*prod_kernel;

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					pointer_yi++;

				}

				if(sum_ker > 0.0)
				{
					mean[j-my_rank*stride] = sum_y_ker/sum_ker;
					trace_H_MPI += prod_kernel_i_eq_j/sum_ker;
				}
				else
				{
					/* Don't print if you are using MPI */
					if((my_rank == 0)&&(int_DEBUG == 1))
					{
						printf("\r                                                                                        ");
						printf("\r** sum_ker[%d]==0.0 in kernel_regression_categorical_aic()",j);
					}
					return(DBL_MAX);
				}

			}

		}

	}
	else
	{

		/* Local Linear */

		XTKX = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKXINV = mat_creat( num_reg_cat_cont + 1, num_reg_cat_cont + 1, UNDEFINED );
		XTKY = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );
		DELTA = mat_creat( num_reg_cat_cont + 1, 1, UNDEFINED );

		/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

		if(kernel_bandwidth_mean(
			KERNEL_reg,
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
			/* Not used */
			matrix_X_continuous,
			/* Not used */
			matrix_X_continuous,
			matrix_X_continuous,
			matrix_X_continuous,
			/* Not used */
			matrix_bandwidth,
			matrix_bandwidth,
			lambda)==1)
		{

			mat_free( XTKX );
			mat_free( XTKXINV );
			mat_free( XTKY );
			mat_free( DELTA );

			free(lambda);
			free_mat(matrix_bandwidth,num_reg_continuous);

			return(DBL_MAX);

		}

		/* Conduct the estimation */

		if(BANDWIDTH_reg == 0)
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][0]);
					}

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
					}

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if((my_rank == 0)&&(int_DEBUG == 1))
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_aic()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
				trace_H_MPI += XTKXINV[0][0]*prod_kernel_i_eq_j;
				mean[j-my_rank*stride] =  DELTA[0][0];

			}

		}
		else if(BANDWIDTH_reg == 1)
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][j]);
					}

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
					}

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if((my_rank == 0)&&(int_DEBUG == 1))
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_aic()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
				trace_H_MPI += XTKXINV[0][0]*prod_kernel_i_eq_j;
				mean[j-my_rank*stride] =  DELTA[0][0];

			}

		}
		else
		{

			for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
			{

				/* Initialize values to zero for a given evaluation point */
				for(k=0; k <= num_reg_cat_cont; k++)
				{
					XTKY[k][0] = 0.0;
					for(l=0; l <= num_reg_cat_cont; l++)
					{
						XTKX[k][l] = 0.0;
					}
				}

				pointer_yi = &vector_Y[0];

				for(i=0; i < num_obs; i++)
				{

					/* Matrix is rows/cols... not C convention... */

					prod_kernel = 1.0;

					for(k = 0; k < num_reg_continuous; k++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous[k][j]-matrix_X_continuous[k][i])/matrix_bandwidth[k][i]);
					}

					for(k = 0; k < num_reg_unordered; k++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[k][j],matrix_X_unordered[k][i],lambda[k],num_categories[k]);
					}

					for(k = 0; k < num_reg_ordered; k++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered[k][j],matrix_X_ordered[k][i],lambda[k+num_reg_unordered]);
					}

					if(i == j)
					{
						prod_kernel_i_eq_j = prod_kernel;
					}

					/* Upper left block */

					XTKX[0][0] += prod_kernel;

					/* First element of XTKY */

					XTKY[0][0] += (temp = (*pointer_yi * prod_kernel));

					for(k=0; k < num_reg_cat_cont; k++)
					{

						/* First lower column of XTKX */

						if(k < num_reg_continuous)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_continuous[k][j] - matrix_X_continuous[k][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_unordered[k-num_reg_continuous][j] - matrix_X_unordered[k-num_reg_continuous][i]))
								* prod_kernel;
						}
						else if(k < num_reg_continuous+num_reg_unordered+num_reg_ordered)
						{
							XTKX[k+1][0] += (temp1 = (matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[k-num_reg_continuous-num_reg_unordered][i]))
								* prod_kernel;
						}

						/* Diagonal of lower block of XTKX */

						XTKX[k+1][k+1] += ipow(temp1, 2) * prod_kernel;

						/* Remaining elements of XTKY */

						XTKY[k+1][0] += temp1 * temp;

						/* Take advantage of symmetric nature of XTKX */

						for(l=0; l < k; l++)
						{
							if(l < num_reg_continuous)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_continuous[l][j] - matrix_X_continuous[l][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_unordered[l-num_reg_continuous][j] - matrix_X_unordered[l-num_reg_continuous][i])
									* prod_kernel;
							}
							else if(l < num_reg_continuous+num_reg_unordered+num_reg_ordered)
							{
								XTKX[k+1][l+1] += temp1 * (matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][j] - matrix_X_ordered[l-num_reg_continuous-num_reg_unordered][i])
									* prod_kernel;
							}
						}

					}

					pointer_yi++;

				}

				for(k=0; k < num_reg_cat_cont; k++)
				{

					/* Take advantage of symmetric nature */

					XTKX[0][k+1] = XTKX[k+1][0];

					for(l=0; l < k; l++)
					{
						XTKX[l+1][k+1] = XTKX[k+1][l+1];
					}

				}

				/* Now compute the beast... */

				if(fabs(mat_det(XTKX)) > 0.0 )
				{

					XTKXINV = mat_inv( XTKX, XTKXINV );

				}
				else
				{

					if((my_rank == 0)&&(int_DEBUG == 1))
					{
						printf("\r                                                                                        ");
						printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical_aic()",j);
						printf("\n");
						mat_dumpf( XTKX, "%g ");
					}

					/* Add ridge factor - epsilon goes from zero to one/n*/

					for(k=0; k < num_reg_cat_cont + 1; k++)
					{
						XTKX[k][k] += epsilon;
					}

					/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

					do
					{
						for(k=0; k < num_reg_cat_cont + 1; k++)
						{
							XTKX[k][k] += epsilon;
							nepsilon += epsilon;
						}
					} while (fabs(mat_det(XTKX)) == 0.0);

					XTKXINV = mat_inv( XTKX, XTKXINV );
					/* Add epsilon times local constant estimator to first element of XTKY */
					XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

				}

				DELTA =  mat_mul( XTKXINV, XTKY, DELTA);
				trace_H_MPI += XTKXINV[0][0]*prod_kernel_i_eq_j;
				mean[j-my_rank*stride] =  DELTA[0][0];

			}

		}

		mat_free( XTKX );
		mat_free( XTKXINV );
		mat_free( XTKY );
		mat_free( DELTA );

	}

	/* Gather */

	MPI_Gather(mean, stride, MPI_DOUBLE, mean, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(mean, num_obs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Reduce(&trace_H_MPI, &trace_H, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&trace_H, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	/* Now compute Hurvich, Simonoff, and Tsai's AICc */

	for(i=0; i < num_obs;i++)
	{
		sigmasq += ipow((vector_Y[i] - mean[i]),2);
	}

	sigmasq /= (double) num_obs;

	aic_c = log(sigmasq) + (1.0+trace_H/((double) num_obs))/(1.0-(trace_H+2.0)/((double) num_obs));

	free(mean);
	free(lambda);
	free_mat(matrix_bandwidth,num_reg_continuous);

	/* Negative penalties are treated as infinite: Hurvich et al pg 277 */

	if((1.0+trace_H/((double) num_obs))/(1.0-(trace_H+2.0)/((double) num_obs))>0.0)
	{
		return(aic_c);
	}
	else
	{
		return(DBL_MAX);
	}

}

int kernel_estimate_ate_categorical_leave_one_out(
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_reg,
int num_obs_train,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double *vector_Y,
double *vector_T,
double *vector_scale_factor,
int *num_categories,
double *mean,
double *tau)
{

	/* This function estimates a density function using both continuous */
	/* and categorical covariates with three estimation techniques and an */
	/* assortment of kernels. */

	/* Declarations */

	int i;
	int j;
	int k;
	int l;

	const double epsilon = 1.0/num_obs_train;
  double nepsilon;

	double prod_kernel;

	double sum_ker;
	double sum_ker_t;
	double sum_ker_t_sq;
	double sum_ker_y;
	double sum_ker_ty;

	double *lambda;
	double **matrix_bandwidth = NULL;

	MATRIX  XTKX;
	MATRIX  XTKXINV;
	MATRIX  XTKY;
	MATRIX  DELTA;

	#ifdef MPI2
	int stride = ceil((double) num_obs_train / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

	/* Allocate memory for objects */

	XTKX = mat_creat( 2, 2, UNDEFINED );
	XTKXINV = mat_creat( 2, 2, UNDEFINED );
	XTKY = mat_creat( 2, 1, UNDEFINED );
	DELTA = mat_creat( 2, 1, UNDEFINED );

	#ifndef MPI2

	lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

	if((BANDWIDTH_reg == 0)||(BANDWIDTH_reg == 1))
	{
		matrix_bandwidth = alloc_matd(num_obs_train,num_reg_continuous);
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
		num_obs_train,
		0,
		0,
		0,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
		vector_scale_factor,
		matrix_X_continuous_train,	 /* Not used */
		matrix_X_continuous_train,	 /* Not used */
		matrix_X_continuous_train,
		matrix_X_continuous_train,
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

	/* Conduct the estimation */

	if(BANDWIDTH_reg == 0)
	{

		for(j=0; j < num_obs_train; j++)
		{

			sum_ker = sum_ker_t = sum_ker_t_sq = sum_ker_y = sum_ker_ty = 0.0;

			for(i=0; i < num_obs_train; i++)
			{
				if(i != j)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous_train[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_train[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_train[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_ker_t += vector_T[i]*prod_kernel;
					sum_ker_t_sq += ipow(vector_T[i],2)*prod_kernel;
					sum_ker_y += vector_Y[i]*prod_kernel;
					sum_ker_ty += vector_T[i]*vector_Y[i]*prod_kernel;

				}

			}

			XTKX[0][0] = sum_ker;
			XTKX[0][1] = sum_ker_t;
			XTKX[1][0] = sum_ker_t;
			XTKX[1][1] = sum_ker_t_sq;
			XTKY[0][0] = sum_ker_y;
			XTKY[1][0]= sum_ker_ty;

			if(fabs(mat_det(XTKX)) > 0.0 )
			{

				XTKXINV = mat_inv( XTKX, XTKXINV );

			}
			else
			{

				if(int_DEBUG == 1)
				{
					printf("\r                                                                                        ");
					printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical()", j);
					printf("\n");
					mat_dumpf( XTKX, "%g ");
				}

				/* Add ridge factor - epsilon goes from zero to one/n*/

				for(k=0; k < 2; k++)
				{
					XTKX[k][k] += epsilon;
				}

				/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

				do
				{
					for(k=0; k <2; k++)
					{
						XTKX[k][k] += epsilon;
						nepsilon += epsilon;
					}
				} while (fabs(mat_det(XTKX)) == 0.0);

				XTKXINV = mat_inv( XTKX, XTKXINV );
				/* Add epsilon times local constant estimator to first element of XTKY */
				XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));

			}

			/*			XTKXINV = mat_inv( XTKX, XTKXINV );*/
			DELTA =  mat_mul( XTKXINV, XTKY, DELTA );
			mean[j] = DELTA[0][0];
			tau[j] = DELTA[1][0];

		}

	}
	else if(BANDWIDTH_reg == 1)
	{

		for(j=0; j < num_obs_train; j++)
		{

			sum_ker = sum_ker_t = sum_ker_t_sq = sum_ker_y = sum_ker_ty = 0.0;

			for(i=0; i < num_obs_train; i++)
			{
				if(i != j)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous_train[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_train[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_train[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_ker_t += vector_T[i]*prod_kernel;
					sum_ker_t_sq += ipow(vector_T[i],2)*prod_kernel;
					sum_ker_y += vector_Y[i]*prod_kernel;
					sum_ker_ty += vector_T[i]*vector_Y[i]*prod_kernel;

				}

			}

			XTKX[0][0] = sum_ker;
			XTKX[0][1] = sum_ker_t;
			XTKX[1][0] = sum_ker_t;
			XTKX[1][1] = sum_ker_t_sq;
			XTKY[0][0] = sum_ker_y;
			XTKY[1][0]= sum_ker_ty;

			XTKXINV = mat_inv( XTKX, XTKXINV );
			DELTA =  mat_mul( XTKXINV, XTKY, DELTA );

			mean[j] = DELTA[0][0];
			tau[j] = DELTA[1][0];

		}

	}
	else
	{

		for(j=0; j < num_obs_train; j++)
		{

			sum_ker = sum_ker_t = sum_ker_t_sq = sum_ker_y = sum_ker_ty = 0.0;

			for(i=0; i < num_obs_train; i++)
			{
				if(i != j)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous_train[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_train[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_train[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_ker_t += vector_T[i]*prod_kernel;
					sum_ker_t_sq += ipow(vector_T[i],2)*prod_kernel;
					sum_ker_y += vector_Y[i]*prod_kernel;
					sum_ker_ty += vector_T[i]*vector_Y[i]*prod_kernel;

				}

			}

			XTKX[0][0] = sum_ker;
			XTKX[0][1] = sum_ker_t;
			XTKX[1][0] = sum_ker_t;
			XTKX[1][1] = sum_ker_t_sq;
			XTKY[0][0] = sum_ker_y;
			XTKY[1][0]= sum_ker_ty;

			XTKXINV = mat_inv( XTKX, XTKXINV );
			DELTA =  mat_mul( XTKXINV, XTKY, DELTA );

			mean[j] = DELTA[0][0];
			tau[j] = DELTA[1][0];

		}

	}
	#endif

	#ifdef MPI2

	lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

	if((BANDWIDTH_reg == 0)||(BANDWIDTH_reg == 1))
	{
		matrix_bandwidth = alloc_matd(stride*iNum_Processors,num_reg_continuous);
	}
	else if(BANDWIDTH_reg == 2)
	{
		matrix_bandwidth = alloc_matd(stride*iNum_Processors,num_reg_continuous);
	}

	/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

	kernel_bandwidth_mean(
		KERNEL_reg,
		BANDWIDTH_reg,
		num_obs_train,
		num_obs_train,
		0,
		0,
		0,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
		vector_scale_factor,
		matrix_X_continuous_train,	 /* Not used */
		matrix_X_continuous_train,	 /* Not used */
		matrix_X_continuous_train,
		matrix_X_continuous_train,
		matrix_bandwidth,						 /* Not used */
		matrix_bandwidth,
		lambda);

	/* Conduct the estimation */

	if(BANDWIDTH_reg == 0)
	{

		for(j=my_rank*stride; (j < num_obs_train) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = sum_ker_t = sum_ker_t_sq = sum_ker_y = sum_ker_ty = 0.0;

			for(i=0; i < num_obs_train; i++)
			{
				if(i != j)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous_train[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_train[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_train[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_ker_t += vector_T[i]*prod_kernel;
					sum_ker_t_sq += ipow(vector_T[i],2)*prod_kernel;
					sum_ker_y += vector_Y[i]*prod_kernel;
					sum_ker_ty += vector_T[i]*vector_Y[i]*prod_kernel;

				}

			}

			XTKX[0][0] = sum_ker;
			XTKX[0][1] = sum_ker_t;
			XTKX[1][0] = sum_ker_t;
			XTKX[1][1] = sum_ker_t_sq;
			XTKY[0][0] = sum_ker_y;
			XTKY[1][0]= sum_ker_ty;

			if(fabs(mat_det(XTKX)) > 0.0 )
			{

				XTKXINV = mat_inv( XTKX, XTKXINV );

			}
			else
			{

				if(int_DEBUG == 1)
				{
					printf("\r                                                                                        ");
					printf("\r** XTKX[%d] is singular: adding ridge factor in kernel_estimate_regression_categorical()", j);
					printf("\n");
					mat_dumpf( XTKX, "%g ");
				}

				/* Add ridge factor - epsilon goes from zero to one/n*/

				for(k=0; k < 2; k++)
				{
					XTKX[k][k] += epsilon;
				}

				/* Keep increasing ridge factor until non-singular - machine and architecture dependent */

				do
				{
					for(k=0; k < 2; k++)
					{
						XTKX[k][k] += epsilon;
						nepsilon += epsilon;
					}
				} while (fabs(mat_det(XTKX)) == 0.0);

				XTKXINV = mat_inv( XTKX, XTKXINV );
				/* Add epsilon times local constant estimator to first element of XTKY */
				XTKY[0][0] += nepsilon*XTKY[0][0]/(MAX(DBL_MIN,XTKX[0][0]));
			}

			/*			XTKXINV = mat_inv( XTKX, XTKXINV );*/
			DELTA =  mat_mul( XTKXINV, XTKY, DELTA );
			mean[j-my_rank*stride] = DELTA[0][0];
			tau[j-my_rank*stride] = DELTA[1][0];

		}

	}
	else if(BANDWIDTH_reg == 1)
	{

		for(j=my_rank*stride; (j < num_obs_train) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = sum_ker_t = sum_ker_t_sq = sum_ker_y = sum_ker_ty = 0.0;

			for(i=0; i < num_obs_train; i++)
			{
				if(i != j)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous_train[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_train[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_train[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_ker_t += vector_T[i]*prod_kernel;
					sum_ker_t_sq += ipow(vector_T[i],2)*prod_kernel;
					sum_ker_y += vector_Y[i]*prod_kernel;
					sum_ker_ty += vector_T[i]*vector_Y[i]*prod_kernel;

				}

			}

			XTKX[0][0] = sum_ker;
			XTKX[0][1] = sum_ker_t;
			XTKX[1][0] = sum_ker_t;
			XTKX[1][1] = sum_ker_t_sq;
			XTKY[0][0] = sum_ker_y;
			XTKY[1][0]= sum_ker_ty;

			XTKXINV = mat_inv( XTKX, XTKXINV );
			DELTA =  mat_mul( XTKXINV, XTKY, DELTA );

			mean[j-my_rank*stride] = DELTA[0][0];
			tau[j-my_rank*stride] = DELTA[1][0];

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_train) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = sum_ker_t = sum_ker_t_sq = sum_ker_y = sum_ker_ty = 0.0;

			for(i=0; i < num_obs_train; i++)
			{
				if(i != j)
				{

					prod_kernel = 1.0;

					for(l = 0; l < num_reg_continuous; l++)
					{
						prod_kernel *= kernel(KERNEL_reg, (matrix_X_continuous_train[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i]);
					}

					for(l = 0; l < num_reg_unordered; l++)
					{
						prod_kernel *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered_train[l][j],matrix_X_unordered_train[l][i],lambda[l],num_categories[l]);
					}

					for(l = 0; l < num_reg_ordered; l++)
					{
						prod_kernel *= kernel_ordered(KERNEL_ordered_reg, matrix_X_ordered_train[l][j],matrix_X_ordered_train[l][i],lambda[l+num_reg_unordered]);
					}

					sum_ker += prod_kernel;
					sum_ker_t += vector_T[i]*prod_kernel;
					sum_ker_t_sq += ipow(vector_T[i],2)*prod_kernel;
					sum_ker_y += vector_Y[i]*prod_kernel;
					sum_ker_ty += vector_T[i]*vector_Y[i]*prod_kernel;

				}

			}

			XTKX[0][0] = sum_ker;
			XTKX[0][1] = sum_ker_t;
			XTKX[1][0] = sum_ker_t;
			XTKX[1][1] = sum_ker_t_sq;
			XTKY[0][0] = sum_ker_y;
			XTKY[1][0]= sum_ker_ty;

			XTKXINV = mat_inv( XTKX, XTKXINV );
			DELTA =  mat_mul( XTKXINV, XTKY, DELTA );

			mean[j-my_rank*stride] = DELTA[0][0];
			tau[j-my_rank*stride] = DELTA[1][0];

		}

	}

	MPI_Gather(mean, stride, MPI_DOUBLE, mean, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(mean, num_obs_train, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Gather(tau, stride, MPI_DOUBLE, tau, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(tau, num_obs_train, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	mat_free( XTKX );
	mat_free( XTKXINV );
	mat_free( XTKY );
	mat_free( DELTA );

	free(lambda);
	free_mat(matrix_bandwidth,num_reg_continuous);

	return(0);

}


/* tristen's modifications to certain functions */

int kecg_compare(const void * a, const void * b)
{
	return ((*(double *)a) == (*(double *)b))?0:(((*(double *)a) > (*(double *)b))*2-1);
}


/* not parallel friendly */

int kernel_estimate_categorical_gradient_ocg_fast(
int int_COMPUTE,
int *var_index_int,
int num_var_test_int,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_reg,
int int_ll,
int int_ordered_categorical_gradient,
int num_obs_train,
int num_obs_eval,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double *vector_Y,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double **matrix_bandwidth,
double **matrix_bandwidth_deriv,
double *lambda,
int *num_categories,
double **matrix_categorical_vals,
double *mean,
double **gradient_categorical)
{

	/* Having a separate module for gradient wrt unordered variables */
	/* is most natural - in particular, local linear cannot be handled with */
	/* a derivative approach, though the Nadaraya-Watson can. Since the local linear */
	/* is so important, we choose to have a separate module.  */

	/*   1) For each categorical variable, create evaluation data set having */
	/*   minimum value for categorical variable. */

	/*   2) Take as argument mean for sample observations, then compute that for*/
	/*   minimum val.*/

	/*   3) Take difference between the two. If var is at min, diff is zero.*/
	/*   Otherwise, get discrete diff between min category and whatever the*/
	/*   variable is for the sample.*/

	int i, j, l, num_obs_eval_alloc;

	double * mean_eval;
	double ** gradient = NULL;

	double **matrix_weights_K = NULL;
	double ***matrix_weights_K_deriv = NULL;

	double **matrix_X_unordered_eval_temp;
	double **matrix_X_ordered_eval_temp;

	double *pointer_m;
	double *pointer_me;
	double *pointer_g;

	double *iord;

	#ifdef MPI2
	num_obs_eval_alloc = MAX(ceil((double) num_obs_eval / (double) iNum_Processors),1)*iNum_Processors;
	#else
	num_obs_eval_alloc = num_obs_eval;
	#endif

	if(!int_COMPUTE)
		return(0);

	mean_eval = alloc_vecd(num_obs_eval_alloc);

	matrix_X_unordered_eval_temp = alloc_matd(num_obs_eval, num_reg_unordered);
	matrix_X_ordered_eval_temp = alloc_matd(num_obs_eval, num_reg_ordered);

	/* unordered variables */
	/* first, copy the evaluation matrix */

	for(l=0; l < num_reg_unordered; l++)
		for(j=0; j < num_obs_eval; j++)
			matrix_X_unordered_eval_temp[l][j] = matrix_X_unordered_eval[l][j];

	/* loop over the matrix, changing one column at a time */
	/* then evaluate, and change the column back */

	for(i=0; i < num_reg_unordered; i++)
	{

		/* change column*/

		for(j=0; j < num_obs_eval; j++)
			matrix_X_unordered_eval_temp[i][j] = matrix_categorical_vals[i][0];

		/* Evaluate on the evaluation data */

		kernel_estimate_regression_categorical_no_stderr(
			0,												 /* Don't compute gradient */
			int_ll,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_reg,
			0,												 /* Don't use weights */
			var_index_int,
			num_var_test_int,
			matrix_weights_K,
			matrix_weights_K_deriv,
			num_obs_train,
			num_obs_eval,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_X_unordered_train,	 /* Train */
			matrix_X_ordered_train,		 /* Train */
			matrix_X_continuous_train,
			matrix_X_unordered_eval_temp,/* Eval */
			matrix_X_ordered_eval,		 /* Eval */
			matrix_X_continuous_eval,
			matrix_bandwidth,
			matrix_bandwidth_deriv,
			vector_Y,
			lambda,
			num_categories,
			mean_eval,
			gradient);

		/* For ith categorical variable, gradient is discrete difference */
		/* between mean for sample observation and mean for minimum category */

		pointer_m = &mean[0];
		pointer_me = &mean_eval[0];
		pointer_g = &gradient_categorical[i][0];

		for(j=0; j < num_obs_eval; j++)
			*pointer_g++ = *pointer_m++ - *pointer_me++;

		/* change back */

		for(j=0; j < num_obs_eval; j++)
			matrix_X_unordered_eval_temp[i][j] = matrix_X_unordered_eval[i][j];

	}

	/* Next for the ordered variables */

	/* use the same procedure for the ordered categories, except we do */
	/* gradients wrt adjacent categories */

	/* copy first */

	for(l=0; l < num_reg_ordered; l++)
		for(j=0; j < num_obs_eval; j++)
			matrix_X_ordered_eval_temp[l][j] = matrix_X_ordered_eval[l][j];

	for(i=0; i < num_reg_ordered; i++)
	{

		/* ideally, matrix_X_ordered_eval would be a matrix of indices of */
		/* matrix_categorical_vals -- then you could make a truly blazing algo */
		for(j=0; j < num_obs_eval; j++)
		{
			iord = (double*) bsearch((const void *)(&matrix_X_ordered_eval[i][j]),
				(const void *)(matrix_categorical_vals[i+num_reg_unordered]),
				(size_t) num_categories[i+num_reg_unordered], sizeof(double),
				kecg_compare);
			if(iord != NULL)
				matrix_X_ordered_eval_temp[i][j] =
					(matrix_X_ordered_eval_temp[i][j] == matrix_categorical_vals[i+num_reg_unordered][0])?
					matrix_categorical_vals[i+num_reg_unordered][1]:iord[-1];
		}

		/* Evaluate on the evaluation data */

		kernel_estimate_regression_categorical_no_stderr(
			0,												 /* Don't compute gradient */
			int_ll,
			KERNEL_reg,
			KERNEL_unordered_reg,
			KERNEL_ordered_reg,
			BANDWIDTH_reg,
			0,												 /* Don't use weights */
			var_index_int,
			num_var_test_int,
			matrix_weights_K,
			matrix_weights_K_deriv,
			num_obs_train,
			num_obs_eval,
			num_reg_unordered,
			num_reg_ordered,
			num_reg_continuous,
			matrix_X_unordered_train,	 /* Train */
			matrix_X_ordered_train,		 /* Train */
			matrix_X_continuous_train,
			matrix_X_unordered_eval,	 /* Eval */
			matrix_X_ordered_eval_temp,/* Eval */
			matrix_X_continuous_eval,
			matrix_bandwidth,
			matrix_bandwidth_deriv,
			vector_Y,
			lambda,
			num_categories,
			mean_eval,
			gradient);

		/* For ith categorical variable, gradient is discrete difference */
		/* between mean for sample observation and mean for minimum category */

		pointer_m = &mean[0];
		pointer_me = &mean_eval[0];
		pointer_g = &gradient_categorical[i+num_reg_unordered][0];

		/* takes the difference and negates it for observations of the first category */

		for(j=0; j < num_obs_eval; j++)
			*pointer_g++ = (*pointer_m++ - *pointer_me++)*
				((matrix_X_ordered_eval[i][j] != matrix_categorical_vals[i+num_reg_unordered][0])?1:-1);

		/* restore the ith column */

		for(j=0; j < num_obs_eval; j++)
			matrix_X_ordered_eval_temp[i][j] = matrix_X_ordered_eval[i][j];

	}

	free(mean_eval);
	free_mat(matrix_X_unordered_eval_temp, num_reg_unordered);
	free_mat(matrix_X_ordered_eval_temp, num_reg_ordered);

	return(0);

}
