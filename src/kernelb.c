/* This module contains the functions for the kernel bandwidth function. */

/* Copyright (C) J. Racine, 1995-2000 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <errno.h>

#include "headers.h"

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

extern int int_LARGE_SF;
/*
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
int int_WEIGHTS;
*/

#ifdef RCSID
static char rcsid[] = "$Id: kernelb.c,v 1.4 2006/11/02 16:56:49 tristen Exp $";
#endif

/* Overloaded these modules to handle conditional distributions */

int kernel_bandwidth(int KERNEL,
int BANDWIDTH,
int num_obs_train,
int num_obs_eval,
int num_var_cont,
int num_var_un,
int num_var_or,
int num_reg_cont,
int num_reg_un,
int num_reg_or,
double *vector_scale_factor,
double **matrix_Y_train,
double **matrix_Y_eval,
double **matrix_X_train,
double **matrix_X_eval,
double **matrix_bandwidth_Y,
double **matrix_bandwidth_X,
double *vector_lambda,
double **matrix_bandwidth_deriv)
{

/* This computes a matrix of bandwidths for fixed, generalized nearest */
/* neighbor, or adaptive nearest neighbor estimation for a density or */
/* regression function as well as for derivatives of the beasts. We permit */
/* two sets of continuous variables, X, and Y, and permit categorical */
/* variables for X and Y as well. Finally, the user has the option of */
/* using/writing the `raw' bandwidth (which approaches zero as n increases) */
/* for all variables, or to use a `normalized' scaling factor which is */
/* constant for all sample sizes. */

	int i;
	int j;

	double temp_inv;

	double temp_pow1 = DBL_MAX;
	double temp_pow2 = DBL_MAX;

	double *vec_sdev_x = NULL;
	double *vec_sdev_y = NULL;

	double *nn_distance = NULL;

	double *pointer_bw;
	double *pointer_bwd;
	double *pointer_nn;

#ifdef MPI2 
	int stride;
#endif

	if(num_obs_train == 0) return(1);

/* Don't compute unnecessary standard deviations  */

	if(int_LARGE_SF == 0)
	{

/* Continuous variables */

		vec_sdev_x = alloc_vecd(num_reg_cont);
		vec_sdev_y = alloc_vecd(num_var_cont);

/* Compute the standard deviation to test for variables which are in */
/* fact constant. */

		for(i=0; i < num_reg_cont; i++)
		{

			vec_sdev_x[i] = standerrd(num_obs_train, matrix_X_train[i]);

			if(vec_sdev_x[i] <= DBL_MIN)
			{
#ifdef MPI2
				if(my_rank == 0) {
#endif
				printf("\r ** Fatal Error in routine kernel_bandwidth() ** The variable appears to be constant!");
				printf("\n ** Program terminated abnormally!\n");
#ifdef MPI2
				}
				/* Since we are exiting and this will terminate the program, clean up */
				MPI_Finalize();
#endif
				exit(EXIT_SUCCESS);
			}

		}

		for(i=0; i < num_var_cont; i++)
		{

			vec_sdev_y[i] = standerrd(num_obs_train, matrix_Y_train[i]);

			if(vec_sdev_y[i] <= DBL_MIN)
			{
#ifdef MPI2
				if(my_rank == 0) {
#endif
				printf("\r ** Fatal Error in routine kernel_bandwidth() ** The variable appears to be constant!");
				printf("\n ** Program terminated abnormally!\n");
#ifdef MPI2
				}
				/* Since we are exiting and this will terminate the program, clean up */
				MPI_Finalize();
#endif
				exit(EXIT_SUCCESS);
			}

		}

	}

#ifndef MPI2

	if(BANDWIDTH == 1)
	{
		nn_distance = alloc_vecd(num_obs_eval);
	}

	if(BANDWIDTH == 2)
	{
		nn_distance = alloc_vecd(num_obs_train);
	}

#endif

#ifdef MPI2

	if(BANDWIDTH == 1)
	{
		stride = ceil((double) num_obs_eval / (double) iNum_Processors);
		if(stride < 1) stride = 1;
		nn_distance = alloc_vecd(stride*iNum_Processors);
	}

	if(BANDWIDTH == 2)
	{
		stride = ceil((double) num_obs_train / (double) iNum_Processors);
		if(stride < 1) stride = 1;
		nn_distance = alloc_vecd(stride*iNum_Processors);
	}

#endif

/* Set appropriate constants for scaling factor */

	switch(KERNEL)
	{

		case 0:

/* Gaussian Kernel */

			temp_pow1 = 1.0/pow((double)num_obs_train, (1.0/(4.0 + (double) num_reg_cont + num_var_cont)));
			temp_pow2 = 1.0/pow((double)num_obs_train, (1.0/(6.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 1:

/* Fourth Order Gaussian Kernel */

			temp_pow1 = 1.0/pow((double)num_obs_train, (1.0/(8.0 + (double) num_reg_cont + num_var_cont)));
			temp_pow2 = 1.0/pow((double)num_obs_train, (1.0/(10.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 2:

/* Sixth Order Gaussian Kernel */

			temp_pow1 = 1.0/pow((double)num_obs_train, (1.0/(12.0 + (double) num_reg_cont + num_var_cont)));
			temp_pow2 = 1.0/pow((double)num_obs_train, (1.0/(14.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 3:

/* Eighth Order Gaussian Kernel */

			temp_pow1 = 1.0/pow((double)num_obs_train, (1.0/(16.0 + (double) num_reg_cont + num_var_cont)));
			temp_pow2 = 1.0/pow((double)num_obs_train, (1.0/(18.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 4:

/* Second Order Epanechnikov Kernel */

			temp_pow1 = 1.0/pow((double)num_obs_train, (1.0/(4.0 + (double) num_reg_cont + num_var_cont)));
			temp_pow2 = 1.0/pow((double)num_obs_train, (1.0/(6.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 5:

/* Fourth Order Epanechnikov Kernel */

			temp_pow1 = 1.0/pow((double)num_obs_train, (1.0/(8.0 + (double) num_reg_cont + num_var_cont)));
			temp_pow2 = 1.0/pow((double)num_obs_train, (1.0/(10.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 6:

/* Sixth Order Epanechnikov Kernel */

			temp_pow1 = 1.0/pow((double)num_obs_train, (1.0/(12.0 + (double) num_reg_cont + num_var_cont)));
			temp_pow2 = 1.0/pow((double)num_obs_train, (1.0/(14.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 7:

/* Eighth Order Epanechnikov Kernel */

			temp_pow1 = 1.0/pow((double)num_obs_train, (1.0/(16.0 + (double) num_reg_cont + num_var_cont)));
			temp_pow2 = 1.0/pow((double)num_obs_train, (1.0/(18.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 8:

/* Rectangular kernel - using second order Epanechnikov for now */

/* Second Order Epanechnikov Kernel */

			temp_pow1 = 1.0/pow((double)num_obs_train, (1.0/(4.0 + (double) num_reg_cont + num_var_cont)));
			temp_pow2 = 1.0/pow((double)num_obs_train, (1.0/(6.0 + (double) num_reg_cont + num_var_cont)));

			break;

	}

	if(BANDWIDTH == 0)
	{

/* fixed */

		for(i=0; i < num_reg_cont; i++)
		{

/* Save on computation since bandwidth is fixed */

			if(int_LARGE_SF == 0)
			{
				matrix_bandwidth_X[i][0] = vector_scale_factor[i] * vec_sdev_x[i] * temp_pow1;
				matrix_bandwidth_deriv[i][0] = vector_scale_factor[i] * vec_sdev_x[i] * temp_pow2;
			}
			else
			{
				matrix_bandwidth_X[i][0] = vector_scale_factor[i];
				matrix_bandwidth_deriv[i][0] = vector_scale_factor[i];
			}

		}

		for(i=0; i < num_var_cont; i++)
		{

/* Save on computation since bandwidth is fixed */

			if(int_LARGE_SF == 0)
			{
				matrix_bandwidth_Y[i][0] = vector_scale_factor[i+num_reg_cont] * vec_sdev_y[i] * temp_pow1;
			}
			else
			{
				matrix_bandwidth_Y[i][0] = vector_scale_factor[i+num_reg_cont];
			}

		}

	}
	else if(BANDWIDTH == 1)
	{

/* Generalized NN */

		for(i=0; i < num_reg_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */

			if(compute_nn_distance_train_eval(num_obs_train,num_obs_eval, matrix_X_train[i], matrix_X_eval[i], fround(vector_scale_factor[i]), nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_X[i][0];
			pointer_bwd = &matrix_bandwidth_deriv[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_eval; j++)
			{

				*pointer_bw++ = *pointer_nn;
				*pointer_bwd++ = *pointer_nn++;

			}

		}

		for(i=0; i < num_var_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */

			if(compute_nn_distance_train_eval(num_obs_train,num_obs_eval, matrix_Y_train[i], matrix_Y_eval[i], fround(vector_scale_factor[i+num_reg_cont]), nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_Y[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_eval; j++)
			{

				*pointer_bw++ = *pointer_nn++;

			}

		}

	}
	else if(BANDWIDTH == 2)
	{

/* Adaptive */

		for(i=0; i < num_reg_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */
			if(compute_nn_distance(num_obs_train, matrix_X_train[i], fround(vector_scale_factor[i]), nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_X[i][0];
			pointer_bwd = &matrix_bandwidth_deriv[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_train; j++)
			{

				*pointer_bw++ = *pointer_nn;
				*pointer_bwd++ = *pointer_nn++;

			}

		}

		for(i=0; i < num_var_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */
			if(compute_nn_distance(num_obs_train, matrix_Y_train[i], fround(vector_scale_factor[i+num_reg_cont]), nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_Y[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_train; j++)
			{

				*pointer_bw++ = *pointer_nn++;

			}

		}

	}                                               /* End generalized NN or adaptive */

/* In vector_scale_factor, order is continuous reg, continuous var, */
/* unordered variables, ordered variables, unordered regressors, ordered regressors */

/* Unordered categorical variables */

	temp_inv = ipow(temp_pow1, 2);

/* Unordered categorical variables */

	for(i=0; i < num_var_un; i++)
	{
		if(int_LARGE_SF == 0)
		{
			vector_lambda[i] = vector_scale_factor[i+num_reg_cont+num_var_cont]*temp_inv;
		}
		else
		{
			vector_lambda[i] = vector_scale_factor[i+num_reg_cont+num_var_cont];
		}
	}

/* Ordered categorical variables */

	for(i=0; i < num_var_or; i++)
	{
		if(int_LARGE_SF == 0)
		{
			vector_lambda[i+num_var_un] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un]*temp_inv;
		}
		else
		{
			vector_lambda[i+num_var_un] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un];
		}
	}

/* Unordered categorical regressors */

	for(i=0; i < num_reg_un; i++)
	{
		if(int_LARGE_SF == 0)
		{
			vector_lambda[i+num_var_un+num_var_or] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un+num_var_or]*temp_inv;
		}
		else
		{
			vector_lambda[i+num_var_un+num_var_or] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un+num_var_or];
		}
	}

/* Ordered categorical regressors */

	for(i=0; i < num_reg_or; i++)
	{
		if(int_LARGE_SF == 0)
		{
			vector_lambda[i+num_var_un+num_var_or+num_reg_un] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un+num_var_or+num_reg_un]*temp_inv;
		}
		else
		{
			vector_lambda[i+num_var_un+num_var_or+num_reg_un] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un+num_var_or+num_reg_un];
		}
	}

	if((BANDWIDTH == 1) || (BANDWIDTH == 2))
	{
		free(nn_distance);
	}

/* Don't compute unnecessary standard deviations  */

	if(int_LARGE_SF == 0)
	{

		free(vec_sdev_x);
		free(vec_sdev_y);

	}

	return(0);

}

int kernel_bandwidth_mean(int KERNEL,
int BANDWIDTH,
int num_obs_train,
int num_obs_eval,
int num_var_cont,
int num_var_un,
int num_var_or,
int num_reg_cont,
int num_reg_un,
int num_reg_or,
double *vector_scale_factor,
double **matrix_Y_train,
double **matrix_Y_eval,
double **matrix_X_train,
double **matrix_X_eval,
double **matrix_bandwidth_Y,
double **matrix_bandwidth_X,
double *vector_lambda)
{

/* This computes a matrix of bandwidths for fixed, generalized nearest */
/* neighbor, or adaptive nearest neighbor estimation for a density or */
/* regression function. We permit two sets of continuous variables, X, and */
/* Y, and permit unordered variables for X and Y as well. Finally, the */
/* user has the option of using/writing the `raw' bandwidth (which */
/* approaches zero as n increases) for all variables, or to use a */
/* `normalized' scaling factor which is constant for all sample sizes. */

	int i;
	int j;

	double temp_pow = DBL_MAX;

	double temp_inv;

	double *vec_sdev_x = NULL;
	double *vec_sdev_y = NULL;

	double *nn_distance = NULL;

	double *pointer_bw;
	double *pointer_nn;

#ifdef MPI2
	int stride;
#endif

	if(num_obs_train == 0) return(1);

/* Don't compute unnecessary standard deviations  */

	if(int_LARGE_SF == 0)
	{

/* Continuous variables */

		vec_sdev_x = alloc_vecd(num_reg_cont);
		vec_sdev_y = alloc_vecd(num_var_cont);


/* Compute the standard deviation to test for variables which are in
fact constant. */

		for(i=0; i < num_reg_cont; i++)
		{

			vec_sdev_x[i] = standerrd(num_obs_train, matrix_X_train[i]);

			if(vec_sdev_x[i] <= DBL_MIN)
			{
#ifdef MPI2
				if(my_rank == 0) {
#endif
				printf("\r ** Fatal Error in routine kernel_bandwidth() ** The variable appears to be constant!");
				printf("\n ** Program terminated abnormally!\n");
#ifdef MPI2
				}
				/* Since we are exiting and this will terminate the program, clean up */
				MPI_Finalize();
#endif
				exit(EXIT_SUCCESS);
			}

		}

		for(i=0; i < num_var_cont; i++)
		{

			vec_sdev_y[i] = standerrd(num_obs_train, matrix_Y_train[i]);

			if(vec_sdev_y[i] <= DBL_MIN)
			{
#ifdef MPI2
				if(my_rank == 0) {
#endif
				printf("\r ** Fatal Error in routine kernel_bandwidth() ** The variable appears to be constant!");
				printf("\n ** Program terminated abnormally!\n");
#ifdef MPI2
				}
				/* Since we are exiting and this will terminate the program, clean up */
				MPI_Finalize();
#endif
				exit(EXIT_SUCCESS);
			}

		}

	}

#ifndef MPI2

	if(BANDWIDTH == 1)
	{
		nn_distance = alloc_vecd(num_obs_eval);
	}

	if(BANDWIDTH == 2)
	{
		nn_distance = alloc_vecd(num_obs_train);
	}

#endif

#ifdef MPI2

	if(BANDWIDTH == 1)
	{
		stride = ceil((double) num_obs_eval / (double) iNum_Processors);
		if(stride < 1) stride = 1;
		nn_distance = alloc_vecd(stride*iNum_Processors);
	}

	if(BANDWIDTH == 2)
	{
		stride = ceil((double) num_obs_train / (double) iNum_Processors);
		if(stride < 1) stride = 1;
		nn_distance = alloc_vecd(stride*iNum_Processors);
	}

#endif

/* Set appropriate constants for scaling factor */

	switch(KERNEL)
	{

		case 0:

/* Gaussian Kernel */

			temp_pow = 1.0/pow((double)num_obs_train, (1.0/(4.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 1:

/* Fourth Order Gaussian Kernel */

			temp_pow = 1.0/pow((double)num_obs_train, (1.0/(8.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 2:

/* Sixth Order Gaussian Kernel */

			temp_pow = 1.0/pow((double)num_obs_train, (1.0/(12.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 3:

/* Eighth Order Gaussian Kernel */

			temp_pow = 1.0/pow((double)num_obs_train, (1.0/(16.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 4:

/* Second Order Epanechnikov Kernel */

			temp_pow = 1.0/pow((double)num_obs_train, (1.0/(4.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 5:

/* Fourth Order Epanechnikov Kernel */

			temp_pow = 1.0/pow((double)num_obs_train, (1.0/(8.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 6:

/* Sixth Order Epanechnikov Kernel */

			temp_pow = 1.0/pow((double)num_obs_train, (1.0/(12.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 7:

/* Eighth Order Epanechnikov Kernel */

			temp_pow = 1.0/pow((double)num_obs_train, (1.0/(16.0 + (double) num_reg_cont + num_var_cont)));

			break;

		case 8:

/* Rectangular kernel - using second order Epanechnikov for now */

/* Second Order Epanechnikov Kernel */

			temp_pow = 1.0/pow((double)num_obs_train, (1.0/(4.0 + (double) num_reg_cont + num_var_cont)));

			break;

	}

	if(BANDWIDTH == 0)
	{

/* fixed */

		for(i=0; i < num_reg_cont; i++)
		{

/* Save on computation since bandwidth is fixed */

			if(int_LARGE_SF == 0)
			{
				matrix_bandwidth_X[i][0] = vector_scale_factor[i] * vec_sdev_x[i] * temp_pow;
			}
			else
			{
				matrix_bandwidth_X[i][0] = vector_scale_factor[i];
			}

		}

		for(i=0; i < num_var_cont; i++)
		{

/* Save on computation since bandwidth is fixed */

			if(int_LARGE_SF == 0)
			{
				matrix_bandwidth_Y[i][0] = vector_scale_factor[i+num_reg_cont] * vec_sdev_y[i] * temp_pow;
			}
			else
			{
				matrix_bandwidth_Y[i][0] = vector_scale_factor[i+num_reg_cont];
			}

		}

	}
	else if(BANDWIDTH == 1)
	{

/* Generalized NN */

		for(i=0; i < num_reg_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */

			if(compute_nn_distance_train_eval(num_obs_train,num_obs_eval, matrix_X_train[i], matrix_X_eval[i], fround(vector_scale_factor[i]), nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_X[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_eval; j++)
			{

				*pointer_bw++ = *pointer_nn++;

			}

		}

		for(i=0; i < num_var_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */

			if(compute_nn_distance_train_eval(num_obs_train,num_obs_eval, matrix_Y_train[i], matrix_Y_eval[i], fround(vector_scale_factor[i+num_reg_cont]), nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_Y[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_eval; j++)
			{

				*pointer_bw++ = *pointer_nn++;

			}

		}

	}
	else if(BANDWIDTH == 2)
	{

/* Adaptive */

		for(i=0; i < num_reg_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */
			if(compute_nn_distance(num_obs_train, matrix_X_train[i], fround(vector_scale_factor[i]), nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_X[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_train; j++)
			{

				*pointer_bw++ = *pointer_nn++;

			}

		}

		for(i=0; i < num_var_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */
			if(compute_nn_distance(num_obs_train, matrix_Y_train[i], fround(vector_scale_factor[i+num_reg_cont]), nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_Y[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_train; j++)
			{

				*pointer_bw++ = *pointer_nn++;

			}

		}

	}                                               /* End generalized NN or adaptive */

/* In vector_scale_factor, order is continuous reg, continuous var, */
/* unordered variables, ordered variables, unordered regressors, ordered regressors */

	temp_inv = ipow(temp_pow, 2);

/* Unordered categorical variables */

	for(i=0; i < num_var_un; i++)
	{
		if(int_LARGE_SF == 0)
		{
			vector_lambda[i] = vector_scale_factor[i+num_reg_cont+num_var_cont]*temp_inv;
		}
		else
		{
			vector_lambda[i] = vector_scale_factor[i+num_reg_cont+num_var_cont];
		}
	}

/* Ordered categorical variables */

	for(i=0; i < num_var_or; i++)
	{
		if(int_LARGE_SF == 0)
		{
			vector_lambda[i+num_var_un] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un]*temp_inv;
		}
		else
		{
			vector_lambda[i+num_var_un] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un];
		}
	}

/* Unordered categorical regressors */

	for(i=0; i < num_reg_un; i++)
	{
		if(int_LARGE_SF == 0)
		{
			vector_lambda[i+num_var_un+num_var_or] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un+num_var_or]*temp_inv;
		}
		else
		{
			vector_lambda[i+num_var_un+num_var_or] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un+num_var_or];
		}
	}

/* Ordered categorical regressors */

	for(i=0; i < num_reg_or; i++)
	{
		if(int_LARGE_SF == 0)
		{
			vector_lambda[i+num_var_un+num_var_or+num_reg_un] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un+num_var_or+num_reg_un]*temp_inv;
		}
		else
		{
			vector_lambda[i+num_var_un+num_var_or+num_reg_un] = vector_scale_factor[i+num_reg_cont+num_var_cont+num_var_un+num_var_or+num_reg_un];
		}
	}

	if((BANDWIDTH == 1) || (BANDWIDTH == 2))
	{
		free(nn_distance);
	}

	if(int_LARGE_SF == 0)
	{

		free(vec_sdev_x);
		free(vec_sdev_y);

	}

	return(0);

}


