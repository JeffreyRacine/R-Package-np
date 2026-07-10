/* This module contains the functions for the kernel bandwidth function. */

/* Copyright (C) J. Racine, 1995-2000 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>
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
extern double nconfac_extern;
extern double ncatfac_extern;
extern double * vector_continuous_stddev_extern;

typedef struct {
	int valid;
	int num_obs_train;
	int num_obs_eval;
	int suppress_parallel;
	int lookup_k;
	const double *train;
	const double *eval;
	uint64_t train_hash;
	uint64_t eval_hash;
	double *distance;
} np_nn_distance_cache_entry;

static np_nn_distance_cache_entry *np_nn_distance_cache = NULL;
static int np_nn_distance_cache_size = 0;
static int np_nn_distance_cache_capacity = 0;

static void np_nn_distance_cache_clear(void)
{
	int i;
	for(i=0; i < np_nn_distance_cache_size; i++)
	{
		safe_free(np_nn_distance_cache[i].distance);
	}
	safe_free(np_nn_distance_cache);
	np_nn_distance_cache = NULL;
	np_nn_distance_cache_size = 0;
	np_nn_distance_cache_capacity = 0;
}

static int np_nn_distance_cache_entry_matches(const np_nn_distance_cache_entry *entry,
const int num_obs_train,
const int num_obs_eval,
const int suppress_parallel,
const double *train,
const double *eval,
const int lookup_k,
const uint64_t train_hash,
const uint64_t eval_hash)
{
	return (entry != NULL) &&
		entry->valid &&
		(entry->num_obs_train == num_obs_train) &&
		(entry->num_obs_eval == num_obs_eval) &&
		(entry->suppress_parallel == suppress_parallel) &&
		(entry->lookup_k == lookup_k) &&
		(entry->train == train) &&
		(entry->eval == eval) &&
		(entry->train_hash == train_hash) &&
		(entry->eval_hash == eval_hash) &&
		(entry->distance != NULL);
}

static uint64_t np_nn_distance_hash_vector(const double *x, const int n)
{
	int i;
	uint64_t h = UINT64_C(1469598103934665603);
	if((x == NULL) || (n <= 0))
	{
		return(h);
	}

	for(i=0; i < n; i++)
	{
		const unsigned char *p = (const unsigned char *)(const void *)&x[i];
		size_t b;
		for(b=0; b < sizeof(double); b++)
		{
			h ^= (uint64_t)p[b];
			h *= UINT64_C(1099511628211);
		}
	}
	return(h);
}

static int np_nn_distance_cache_find(const int num_obs_train,
const int num_obs_eval,
const int suppress_parallel,
const double *train,
const double *eval,
const int lookup_k,
const uint64_t train_hash,
const uint64_t eval_hash)
{
	int i;
	for(i=0; i < np_nn_distance_cache_size; i++)
	{
		if(np_nn_distance_cache_entry_matches(&np_nn_distance_cache[i],
		                                       num_obs_train,
		                                       num_obs_eval,
		                                       suppress_parallel,
		                                       train,
		                                       eval,
		                                       lookup_k,
		                                       train_hash,
		                                       eval_hash))
		{
			return(i);
		}
	}
	return(-1);
}

static np_nn_distance_cache_entry *np_nn_distance_cache_add(const int num_obs_train,
const int num_obs_eval,
const int suppress_parallel,
const double *train,
const double *eval,
const int lookup_k,
const uint64_t train_hash,
const uint64_t eval_hash,
const double *distance)
{
	np_nn_distance_cache_entry *tmp;
	np_nn_distance_cache_entry *entry;
	double *copy;

	if((num_obs_eval < 1) || (distance == NULL))
	{
		return(NULL);
	}

	if(np_nn_distance_cache_size >= 64)
	{
		np_nn_distance_cache_clear();
	}

	if(np_nn_distance_cache_size >= np_nn_distance_cache_capacity)
	{
		const int new_capacity = (np_nn_distance_cache_capacity == 0) ? 8 : 2*np_nn_distance_cache_capacity;
		tmp = (np_nn_distance_cache_entry *)realloc(np_nn_distance_cache,
		                                             (size_t)new_capacity * sizeof(np_nn_distance_cache_entry));
		if(tmp == NULL)
		{
			np_nn_distance_cache_clear();
			return(NULL);
		}
		np_nn_distance_cache = tmp;
		memset(np_nn_distance_cache + np_nn_distance_cache_capacity,
		       0,
		       (size_t)(new_capacity - np_nn_distance_cache_capacity) * sizeof(np_nn_distance_cache_entry));
		np_nn_distance_cache_capacity = new_capacity;
	}

	copy = (double *)malloc((size_t)num_obs_eval * sizeof(double));
	if(copy == NULL)
	{
		return(NULL);
	}
	memcpy(copy, distance, (size_t)num_obs_eval * sizeof(double));

	entry = &np_nn_distance_cache[np_nn_distance_cache_size++];
	entry->valid = 1;
	entry->num_obs_train = num_obs_train;
	entry->num_obs_eval = num_obs_eval;
	entry->suppress_parallel = suppress_parallel;
	entry->lookup_k = lookup_k;
	entry->train = train;
	entry->eval = eval;
	entry->train_hash = train_hash;
	entry->eval_hash = eval_hash;
	entry->distance = copy;
	return(entry);
}

static int np_compute_nn_distance_train_eval_cached(const int num_obs_train,
const int num_obs_eval,
const int suppress_parallel,
double *vector_data_train,
double *vector_data_eval,
const int lookup_k,
const int use_cache,
double *nn_distance)
{
	int idx;
	uint64_t train_hash = 0;
	uint64_t eval_hash = 0;
	if(!use_cache)
	{
		return(compute_nn_distance_train_eval(num_obs_train,
		                                      num_obs_eval,
		                                      suppress_parallel,
		                                      vector_data_train,
		                                      vector_data_eval,
		                                      lookup_k,
		                                      nn_distance));
	}

	train_hash = np_nn_distance_hash_vector(vector_data_train, num_obs_train);
	eval_hash = np_nn_distance_hash_vector(vector_data_eval, num_obs_eval);
	idx = np_nn_distance_cache_find(num_obs_train,
	                                 num_obs_eval,
	                                 suppress_parallel,
	                                 vector_data_train,
	                                 vector_data_eval,
	                                 lookup_k,
	                                 train_hash,
	                                 eval_hash);
	if(idx >= 0)
	{
		memcpy(nn_distance,
		       np_nn_distance_cache[idx].distance,
		       (size_t)num_obs_eval * sizeof(double));
		return(0);
	}

	if(compute_nn_distance_train_eval(num_obs_train,
	                                  num_obs_eval,
	                                  suppress_parallel,
	                                  vector_data_train,
	                                  vector_data_eval,
	                                  lookup_k,
	                                  nn_distance)==1)
	{
		return(1);
	}

	np_nn_distance_cache_add(num_obs_train,
	                          num_obs_eval,
	                          suppress_parallel,
	                          vector_data_train,
	                          vector_data_eval,
	                          lookup_k,
	                          train_hash,
	                          eval_hash,
	                          nn_distance);
	return(0);
}

static int np_compute_nn_distance_cached(const int num_obs,
const int suppress_parallel,
double *vector_data,
const int lookup_k,
const int use_cache,
double *nn_distance)
{
	int idx;
	uint64_t train_hash = 0;
	if(!use_cache)
	{
		return(compute_nn_distance(num_obs,
		                           suppress_parallel,
		                           vector_data,
		                           lookup_k,
		                           nn_distance));
	}

	train_hash = np_nn_distance_hash_vector(vector_data, num_obs);
	idx = np_nn_distance_cache_find(num_obs,
	                                num_obs,
	                                suppress_parallel,
	                                vector_data,
	                                NULL,
	                                lookup_k,
	                                train_hash,
	                                UINT64_C(0));
	if(idx >= 0)
	{
		memcpy(nn_distance,
		       np_nn_distance_cache[idx].distance,
		       (size_t)num_obs * sizeof(double));
		return(0);
	}

	if(compute_nn_distance(num_obs,
	                       suppress_parallel,
	                       vector_data,
	                       lookup_k,
	                       nn_distance)==1)
	{
		return(1);
	}

	np_nn_distance_cache_add(num_obs,
	                         num_obs,
	                         suppress_parallel,
	                         vector_data,
	                         NULL,
	                         lookup_k,
	                         train_hash,
	                         UINT64_C(0),
	                         nn_distance);
	return(0);
}

static int np_extendednn_enabled(void)
{
	const SEXP val = Rf_GetOption1(Rf_install("np.extendednn"));
	const int flag = Rf_asLogical(val);
	return flag == TRUE;
}

static int np_nn_lookup_from_scale(const int num_obs_train,
const int allow_extended,
const double scale_factor,
int *lookup_k,
double *distance_scale,
int *is_extended)
{
	const int max_k = num_obs_train - 1;
	int rounded_k;

	if((lookup_k == NULL) || (distance_scale == NULL) || (max_k < 1))
	{
		return(1);
	}

	if(is_extended != NULL)
	{
		*is_extended = 0;
	}

	if(!isfinite(scale_factor) || (scale_factor < 1.0) || (scale_factor > ((double)INT_MAX / 2.0)))
	{
		return(1);
	}

	rounded_k = np_fround(scale_factor);

	if(rounded_k < 1)
	{
		return(1);
	}

	if(rounded_k > max_k)
	{
		if(!allow_extended || !np_extendednn_enabled())
		{
			return(1);
		}

		*lookup_k = max_k;
		*distance_scale = ((double)rounded_k)/((double)max_k);
		if(is_extended != NULL)
		{
			*is_extended = 1;
		}
		return(0);
	}

	*lookup_k = rounded_k;
	*distance_scale = 1.0;
	return(0);
}

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
	double nn_scale;
	int nn_extended;
	int int_nn_k;

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
				/* Since we are exiting and this will terminate the program, clean up */
				MPI_Finalize();
#endif
				error("\r ** Fatal Error in routine kernel_bandwidth() ** The variable appears to be constant!");
			}

		}

		for(i=0; i < num_var_cont; i++)
		{

			vec_sdev_y[i] = standerrd(num_obs_train, matrix_Y_train[i]);

			if(vec_sdev_y[i] <= DBL_MIN)
			{
#ifdef MPI2
				/* Since we are exiting and this will terminate the program, clean up */
				MPI_Finalize();
#endif
				error("\r ** Fatal Error in routine kernel_bandwidth() ** The variable appears to be constant!");
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
		stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
		if(stride < 1) stride = 1;
		nn_distance = alloc_vecd(stride*iNum_Processors);
	}

	if(BANDWIDTH == 2)
	{
		stride = (int)ceil((double) num_obs_train / (double) iNum_Processors);
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

		case 9:

/* Gaussian Kernel */

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

			if(np_nn_lookup_from_scale(num_obs_train, 1, vector_scale_factor[i], &int_nn_k, &nn_scale, &nn_extended)==1)
			{
				return(1);
			}

			if(np_compute_nn_distance_train_eval_cached(num_obs_train,num_obs_eval, 0,matrix_X_train[i], matrix_X_eval[i], int_nn_k, 1, nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_X[i][0];
			pointer_bwd = &matrix_bandwidth_deriv[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_eval; j++)
			{

				*pointer_bw++ = nn_scale * *pointer_nn;
				*pointer_bwd++ = nn_scale * *pointer_nn++;

			}

		}

		for(i=0; i < num_var_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */

			if(np_nn_lookup_from_scale(num_obs_train, 1, vector_scale_factor[i+num_reg_cont], &int_nn_k, &nn_scale, &nn_extended)==1)
			{
				return(1);
			}

			if(np_compute_nn_distance_train_eval_cached(num_obs_train,num_obs_eval, 0, matrix_Y_train[i], matrix_Y_eval[i], int_nn_k, 1, nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_Y[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_eval; j++)
			{

				*pointer_bw++ = nn_scale * *pointer_nn++;

			}

		}

	}
	else if(BANDWIDTH == 2)
	{

/* Adaptive */

		for(i=0; i < num_reg_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */
			if(np_nn_lookup_from_scale(num_obs_train, 1, vector_scale_factor[i], &int_nn_k, &nn_scale, &nn_extended)==1)
			{
				return(1);
			}

			if(np_compute_nn_distance_cached(num_obs_train, 0, matrix_X_train[i], int_nn_k, 1, nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_X[i][0];
			pointer_bwd = &matrix_bandwidth_deriv[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_train; j++)
			{

				*pointer_bw++ = nn_scale * *pointer_nn;
				*pointer_bwd++ = nn_scale * *pointer_nn++;

			}

		}

		for(i=0; i < num_var_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */
			if(np_nn_lookup_from_scale(num_obs_train, 1, vector_scale_factor[i+num_reg_cont], &int_nn_k, &nn_scale, &nn_extended)==1)
			{
				return(1);
			}

			if(np_compute_nn_distance_cached(num_obs_train, 0, matrix_Y_train[i], int_nn_k, 1, nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_Y[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_train; j++)
			{

				*pointer_bw++ = nn_scale * *pointer_nn++;

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
                          int suppress_parallel,
                          double *vector_scale_factor,
                          double **matrix_Y_train,
                          double **matrix_Y_eval,
                          double **matrix_X_train,
                          double **matrix_X_eval,
                          double **matrix_bandwidth_Y,
                          double **matrix_bandwidth_X,
                          double *vector_lambda){

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
	double nn_scale;
	int nn_extended;
	int int_nn_k;

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

			//vec_sdev_x[i] = standerrd(num_obs_train, matrix_X_train[i]);
      vec_sdev_x[i] = vector_continuous_stddev_extern[i];

			if(vec_sdev_x[i] <= DBL_MIN)
			{
#ifdef MPI2
				/* Since we are exiting and this will terminate the program, clean up */
				MPI_Finalize();
#endif
				error("\r ** Fatal Error in routine kernel_bandwidth() ** The variable appears to be constant!");
			}

		}

		for(i=0; i < num_var_cont; i++)
		{

			//vec_sdev_y[i] = standerrd(num_obs_train, matrix_Y_train[i]);
      vec_sdev_y[i] = vector_continuous_stddev_extern[i+num_reg_cont];

			if(vec_sdev_y[i] <= DBL_MIN)
			{
#ifdef MPI2
				/* Since we are exiting and this will terminate the program, clean up */
				MPI_Finalize();
#endif
				error("\r ** Fatal Error in routine kernel_bandwidth() ** The variable appears to be constant!");
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
		stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
		if(stride < 1) stride = 1;
		nn_distance = alloc_vecd(stride*iNum_Processors);
	}

	if(BANDWIDTH == 2)
	{
		stride = (int)ceil((double) num_obs_train / (double) iNum_Processors);
		if(stride < 1) stride = 1;
		nn_distance = alloc_vecd(stride*iNum_Processors);
	}

#endif

/* Set appropriate constants for scaling factor */
  temp_pow = nconfac_extern;
	

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

			if(np_nn_lookup_from_scale(num_obs_train, 1, vector_scale_factor[i], &int_nn_k, &nn_scale, &nn_extended)==1)
			{
				return(1);
			}

			if(np_compute_nn_distance_train_eval_cached(num_obs_train,num_obs_eval, suppress_parallel, matrix_X_train[i], matrix_X_eval[i], int_nn_k, 1, nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_X[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_eval; j++)
			{

				*pointer_bw++ = nn_scale * *pointer_nn++;

			}

		}

		for(i=0; i < num_var_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */

			if(np_nn_lookup_from_scale(num_obs_train, 1, vector_scale_factor[i+num_reg_cont], &int_nn_k, &nn_scale, &nn_extended)==1)
			{
				return(1);
			}

			if(np_compute_nn_distance_train_eval_cached(num_obs_train,num_obs_eval, suppress_parallel, matrix_Y_train[i], matrix_Y_eval[i], int_nn_k, 1, nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_Y[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_eval; j++)
			{

				*pointer_bw++ = nn_scale * *pointer_nn++;

			}

		}

	}
	else if(BANDWIDTH == 2)
	{

/* Adaptive */

		for(i=0; i < num_reg_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */
			if(np_nn_lookup_from_scale(num_obs_train, 1, vector_scale_factor[i], &int_nn_k, &nn_scale, &nn_extended)==1)
			{
				return(1);
			}

			if(np_compute_nn_distance_cached(num_obs_train, suppress_parallel, matrix_X_train[i], int_nn_k, 1, nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_X[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_train; j++)
			{

				*pointer_bw++ = nn_scale * *pointer_nn++;

			}

		}

		for(i=0; i < num_var_cont; i++)
		{

/* Return 1 for nearest-neighbor which is zero */
			if(np_nn_lookup_from_scale(num_obs_train, 1, vector_scale_factor[i+num_reg_cont], &int_nn_k, &nn_scale, &nn_extended)==1)
			{
				return(1);
			}

			if(np_compute_nn_distance_cached(num_obs_train, suppress_parallel, matrix_Y_train[i], int_nn_k, 1, nn_distance)==1)
			{
				return(1);
			}

/* Compute the nearest neighbor distances */

			pointer_bw = &matrix_bandwidth_Y[i][0];
			pointer_nn = &nn_distance[0];

			for(j=0; j < num_obs_train; j++)
			{

				*pointer_bw++ = nn_scale * *pointer_nn++;

			}

		}

	}                                               /* End generalized NN or adaptive */

/* In vector_scale_factor, order is continuous reg, continuous var, */
/* unordered variables, ordered variables, unordered regressors, ordered regressors */

	temp_inv = ncatfac_extern;

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

void np_nn_distance_cache_clear_extern(void)
{
	np_nn_distance_cache_clear();
}
