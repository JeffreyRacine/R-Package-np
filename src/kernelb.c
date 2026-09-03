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
#include "np_native_safety.h"

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

#ifndef MPI2
typedef struct {
	int valid;
	int geometry_contract;
	NPNNQueryMode query_mode;
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

enum {
	NP_NN_CACHE_INCUMBENT = 0,
	NP_NN_CACHE_CANONICAL = 1
};

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
const int geometry_contract,
const NPNNQueryMode query_mode,
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
		(entry->geometry_contract == geometry_contract) &&
		(entry->query_mode == query_mode) &&
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
const int geometry_contract,
const NPNNQueryMode query_mode,
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
		                                       geometry_contract,
		                                       query_mode,
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
const int geometry_contract,
const NPNNQueryMode query_mode,
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
	entry->geometry_contract = geometry_contract;
	entry->query_mode = query_mode;
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

static NPNNGeometryStatus np_compute_nn_distance_train_eval_cached(const int num_obs_train,
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
		return compute_nn_distance_train_eval_status(
		  num_obs_train, num_obs_eval, suppress_parallel,
		  vector_data_train, vector_data_eval, lookup_k, nn_distance);
	}

	train_hash = np_nn_distance_hash_vector(vector_data_train, num_obs_train);
	eval_hash = np_nn_distance_hash_vector(vector_data_eval, num_obs_eval);
	idx = np_nn_distance_cache_find(num_obs_train,
	                                 num_obs_eval,
	                                 suppress_parallel,
	                                 NP_NN_CACHE_INCUMBENT,
	                                 NP_NN_QUERY_EXTERNAL,
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
		return NP_NN_GEOMETRY_OK;
	}

	{
		const NPNNGeometryStatus status = compute_nn_distance_train_eval_status(
		  num_obs_train, num_obs_eval, suppress_parallel,
		  vector_data_train, vector_data_eval, lookup_k, nn_distance);
		if(status != NP_NN_GEOMETRY_OK)
			return status;
	}

	np_nn_distance_cache_add(num_obs_train,
	                          num_obs_eval,
	                          suppress_parallel,
	                          NP_NN_CACHE_INCUMBENT,
	                          NP_NN_QUERY_EXTERNAL,
	                          vector_data_train,
	                          vector_data_eval,
	                          lookup_k,
	                          train_hash,
	                          eval_hash,
	                          nn_distance);
	return NP_NN_GEOMETRY_OK;
}

static NPNNGeometryStatus np_compute_nn_distance_train_eval_context_cached(
const int num_obs_train,
const int num_obs_eval,
const int suppress_parallel,
const double *vector_data_train,
const double *vector_data_eval,
const int lookup_k,
const NPNNGeometryContext *geometry_context,
const int use_cache,
double *nn_distance)
{
	int idx;
	uint64_t train_hash;
	uint64_t eval_hash;
	NPNNGeometryStatus status;

	if(geometry_context == NULL)
		return NP_NN_GEOMETRY_INVALID_ARGUMENT;
	/* A map is caller-owned mutable identity state; bypass until a bounded,
	 * content-complete map key is introduced with a measured need. */
	if(geometry_context->mode == NP_NN_QUERY_TRAINING_MAP)
		return compute_nn_distance_train_eval_ctx(
			num_obs_train, num_obs_eval, suppress_parallel,
			vector_data_train, vector_data_eval, lookup_k,
			geometry_context, nn_distance);
	if(!use_cache)
		return compute_nn_distance_train_eval_ctx(
			num_obs_train, num_obs_eval, suppress_parallel,
			vector_data_train, vector_data_eval, lookup_k,
			geometry_context, nn_distance);

	train_hash = np_nn_distance_hash_vector(vector_data_train, num_obs_train);
	eval_hash = np_nn_distance_hash_vector(vector_data_eval, num_obs_eval);

	idx = np_nn_distance_cache_find(
		num_obs_train, num_obs_eval, suppress_parallel,
		NP_NN_CACHE_CANONICAL, geometry_context->mode,
		vector_data_train, vector_data_eval,
		lookup_k, train_hash, eval_hash);
	if(idx >= 0)
	{
		memcpy(nn_distance, np_nn_distance_cache[idx].distance,
		       (size_t)num_obs_eval*sizeof(double));
		return NP_NN_GEOMETRY_OK;
	}

	status = compute_nn_distance_train_eval_ctx(
		num_obs_train, num_obs_eval, suppress_parallel,
		vector_data_train, vector_data_eval, lookup_k,
		geometry_context, nn_distance);
	if(status != NP_NN_GEOMETRY_OK)
		return status;

	np_nn_distance_cache_add(
		num_obs_train, num_obs_eval, suppress_parallel,
		NP_NN_CACHE_CANONICAL, geometry_context->mode,
		vector_data_train, vector_data_eval,
		lookup_k, train_hash, eval_hash, nn_distance);
	return NP_NN_GEOMETRY_OK;
}

static NPNNGeometryStatus np_compute_nn_distance_cached(const int num_obs,
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
		return compute_nn_distance_status(
		  num_obs, suppress_parallel, vector_data, lookup_k, nn_distance);
	}

	train_hash = np_nn_distance_hash_vector(vector_data, num_obs);
	idx = np_nn_distance_cache_find(num_obs,
	                                num_obs,
	                                suppress_parallel,
	                                NP_NN_CACHE_INCUMBENT,
	                                NP_NN_QUERY_EXTERNAL,
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
		return NP_NN_GEOMETRY_OK;
	}

	{
		const NPNNGeometryStatus status = compute_nn_distance_status(
		  num_obs, suppress_parallel, vector_data, lookup_k, nn_distance);
		if(status != NP_NN_GEOMETRY_OK)
			return status;
	}

	np_nn_distance_cache_add(num_obs,
	                         num_obs,
	                         suppress_parallel,
	                         NP_NN_CACHE_INCUMBENT,
	                         NP_NN_QUERY_EXTERNAL,
	                         vector_data,
	                         NULL,
	                         lookup_k,
	                         train_hash,
	                         UINT64_C(0),
	                         nn_distance);
	return NP_NN_GEOMETRY_OK;
}
#endif

static int np_extendednn_enabled(void)
{
	const SEXP val = Rf_GetOption1(Rf_install("np.extendednn"));
	const int flag = Rf_asLogical(val);
	return flag == TRUE;
}

int np_nn_lookup_from_scale(const int num_obs_train,
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

/* Exact adaptive delete-one rows have one fewer donor-neighbour than the
 * full sample. Resolve the requested observation count against that fold
 * cardinality once, using the same extended-NN linear radius convention as
 * the ordinary bandwidth owner. */
static int np_adaptive_fold_lookup_from_scale(const int num_obs_train,
                                              const double scale_factor,
                                              int *lookup_k,
                                              double *distance_scale)
{
	const int max_k = num_obs_train - 2;
	int requested_k;

	if(lookup_k == NULL || distance_scale == NULL || max_k < 1 ||
	   !isfinite(scale_factor) || scale_factor < 1.0 ||
	   scale_factor > ((double)INT_MAX / 2.0))
		return(1);

	requested_k = np_fround(scale_factor);
	if(requested_k < 1)
		return(1);
	if(requested_k > max_k && !np_extendednn_enabled())
		return(1);

	*lookup_k = MIN(requested_k, max_k);
	*distance_scale = ((double)requested_k)/((double)*lookup_k);
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

int kernel_bandwidth_ctx(int KERNEL,
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
double **matrix_bandwidth_deriv,
const NPNNGeometryContext *x_geometry_context,
const NPNNGeometryContext *y_geometry_context,
NPNNGeometryStatus *geometry_status)
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
	int status = 0;

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
	int int_nn_k;
	int nn_extended;
	NPNNGeometryStatus nn_geometry_status = NP_NN_GEOMETRY_OK;

#ifdef MPI2 
	int nn_distance_alloc = 0;
#endif

	if(geometry_status != NULL)
		*geometry_status = NP_NN_GEOMETRY_OK;
	if(num_obs_train == 0) return(1);

#ifdef MPI2
	if((BANDWIDTH == 1) || (BANDWIDTH == 2))
	{
		int stride_ignored = 0;
		const int nn_count = (BANDWIDTH == 1) ? num_obs_eval : num_obs_train;

		if(!np_int_padded_count_nonnegative(nn_count,
		                                    iNum_Processors,
		                                    1,
		                                    &stride_ignored,
		                                    &nn_distance_alloc))
			return(1);
	}
#endif

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
				error("\r ** Fatal Error in routine kernel_bandwidth() ** The variable appears to be constant!");
			}

		}

		for(i=0; i < num_var_cont; i++)
		{

			vec_sdev_y[i] = standerrd(num_obs_train, matrix_Y_train[i]);

			if(vec_sdev_y[i] <= DBL_MIN)
			{
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
		nn_distance = alloc_vecd(nn_distance_alloc);
	}

	if(BANDWIDTH == 2)
	{
		nn_distance = alloc_vecd(nn_distance_alloc);
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
				status = 1;
				goto cleanup;
			}


#ifndef MPI2
			if(x_geometry_context != NULL)
				nn_geometry_status = np_compute_nn_distance_train_eval_context_cached(
					num_obs_train, num_obs_eval, 0,
					matrix_X_train[i], matrix_X_eval[i], int_nn_k,
					x_geometry_context, 1, nn_distance);
			else
				nn_geometry_status = np_compute_nn_distance_train_eval_cached(
					num_obs_train, num_obs_eval, 0,
					matrix_X_train[i], matrix_X_eval[i], int_nn_k,
					1, nn_distance);
#else
			if(x_geometry_context != NULL)
				nn_geometry_status = compute_nn_distance_train_eval_ctx(
					num_obs_train, num_obs_eval, 0,
					matrix_X_train[i], matrix_X_eval[i], int_nn_k,
					x_geometry_context, nn_distance);
			else
				nn_geometry_status = compute_nn_distance_train_eval_status(
					num_obs_train, num_obs_eval, 0,
					matrix_X_train[i], matrix_X_eval[i], int_nn_k,
					nn_distance);
#endif
			if(nn_geometry_status != NP_NN_GEOMETRY_OK)
			{
				if(geometry_status != NULL)
					*geometry_status = nn_geometry_status;
				status = 1;
				goto cleanup;
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
				status = 1;
				goto cleanup;
			}

			nn_geometry_status = NP_NN_GEOMETRY_OK;
#ifndef MPI2
			if(y_geometry_context != NULL)
				nn_geometry_status = np_compute_nn_distance_train_eval_context_cached(
					num_obs_train, num_obs_eval, 0,
					matrix_Y_train[i], matrix_Y_eval[i], int_nn_k,
					y_geometry_context, 1, nn_distance);
			else
				nn_geometry_status = np_compute_nn_distance_train_eval_cached(
					num_obs_train, num_obs_eval, 0,
					matrix_Y_train[i], matrix_Y_eval[i], int_nn_k,
					1, nn_distance);
#else
			if(y_geometry_context != NULL)
				nn_geometry_status = compute_nn_distance_train_eval_ctx(
					num_obs_train, num_obs_eval, 0,
					matrix_Y_train[i], matrix_Y_eval[i], int_nn_k,
					y_geometry_context, nn_distance);
			else
				nn_geometry_status = compute_nn_distance_train_eval_status(
					num_obs_train, num_obs_eval, 0,
					matrix_Y_train[i], matrix_Y_eval[i], int_nn_k,
					nn_distance);
#endif
			if(nn_geometry_status != NP_NN_GEOMETRY_OK)
			{
				if(geometry_status != NULL)
					*geometry_status = nn_geometry_status;
				status = 1;
				goto cleanup;
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
				status = 1;
				goto cleanup;
			}

			nn_geometry_status = compute_nn_distance_status(
				num_obs_train, 0, matrix_X_train[i], int_nn_k, nn_distance);
			if(nn_geometry_status != NP_NN_GEOMETRY_OK)
			{
				if(geometry_status != NULL)
					*geometry_status = nn_geometry_status;
				status = 1;
				goto cleanup;
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
				status = 1;
				goto cleanup;
			}

			nn_geometry_status = compute_nn_distance_status(
				num_obs_train, 0, matrix_Y_train[i], int_nn_k, nn_distance);
			if(nn_geometry_status != NP_NN_GEOMETRY_OK)
			{
				if(geometry_status != NULL)
					*geometry_status = nn_geometry_status;
				status = 1;
				goto cleanup;
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

cleanup:
	free(nn_distance);
	free(vec_sdev_x);
	free(vec_sdev_y);

	return(status);

	}

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
	return kernel_bandwidth_ctx(
		KERNEL, BANDWIDTH, num_obs_train, num_obs_eval,
		num_var_cont, num_var_un, num_var_or,
		num_reg_cont, num_reg_un, num_reg_or,
		vector_scale_factor,
		matrix_Y_train, matrix_Y_eval,
		matrix_X_train, matrix_X_eval,
		matrix_bandwidth_Y, matrix_bandwidth_X,
		vector_lambda, matrix_bandwidth_deriv,
		NULL, NULL, NULL);
}

static int np_kernel_bandwidth_continuous_nn_into_ctx(
  int BANDWIDTH,
  int num_obs_train,
  int num_obs_eval,
  int num_cont,
  int suppress_parallel,
  double *vector_scale_factor,
  double **matrix_train,
  double **matrix_eval,
  double **matrix_bandwidth,
  double *nn_distance,
  const NPNNGeometryContext *geometry_context,
  NPNNGeometryStatus *geometry_status)
{
  int dimension;
  int observation;
  int int_nn_k;
  double nn_scale;
  int nn_extended;

  if(geometry_status != NULL)
    *geometry_status = NP_NN_GEOMETRY_OK;

  for(dimension = 0; dimension < num_cont; ++dimension) {
    if(np_nn_lookup_from_scale(num_obs_train, 1,
                               vector_scale_factor[dimension],
                               &int_nn_k, &nn_scale, &nn_extended) == 1)
      return 1;

    if(BANDWIDTH == BW_GEN_NN) {
      NPNNGeometryStatus query_status = NP_NN_GEOMETRY_OK;
      if(geometry_context != NULL) {
#ifndef MPI2
        query_status = np_compute_nn_distance_train_eval_context_cached(
          num_obs_train, num_obs_eval, suppress_parallel,
          matrix_train[dimension], matrix_eval[dimension],
          int_nn_k, geometry_context, 1, nn_distance);
#else
        query_status = compute_nn_distance_train_eval_ctx(
          num_obs_train, num_obs_eval, suppress_parallel,
          matrix_train[dimension], matrix_eval[dimension],
          int_nn_k, geometry_context, nn_distance);
#endif
      } else {
#ifndef MPI2
        query_status = np_compute_nn_distance_train_eval_cached(
          num_obs_train, num_obs_eval, suppress_parallel,
          matrix_train[dimension], matrix_eval[dimension],
          int_nn_k, 1, nn_distance);
#else
        query_status = compute_nn_distance_train_eval_status(
          num_obs_train, num_obs_eval, suppress_parallel,
          matrix_train[dimension], matrix_eval[dimension],
          int_nn_k, nn_distance);
#endif
      }
      if(query_status != NP_NN_GEOMETRY_OK) {
        if(geometry_status != NULL)
          *geometry_status = query_status;
        return 1;
      }
      for(observation = 0; observation < num_obs_eval; ++observation) {
        const double bandwidth = nn_scale*nn_distance[observation];
        if(!isfinite(bandwidth) || bandwidth <= 0.0) {
          if(geometry_status != NULL)
            *geometry_status = !isfinite(bandwidth) ?
              NP_NN_GEOMETRY_NONFINITE_RADIUS : NP_NN_GEOMETRY_ZERO_RADIUS;
          return 1;
        }
        matrix_bandwidth[dimension][observation] = bandwidth;
      }
    } else {
      if(geometry_context != NULL &&
         geometry_context->mode == NP_NN_QUERY_ADAPTIVE_FOLD_PREPARE){
        int fold_lookup_k;
        double fold_scale;
        NPNNGeometryStatus pair_status;

        if(num_obs_train < 3 || num_obs_eval != num_obs_train ||
           geometry_context->eval_to_train != NULL ||
           geometry_context->adaptive_successor == NULL ||
           geometry_context->adaptive_fold_scale == NULL ||
           geometry_context->adaptive_successor[dimension] == NULL ||
           (geometry_context->adaptive_full != NULL &&
            geometry_context->adaptive_full[dimension] == NULL)){
          if(geometry_status != NULL)
            *geometry_status = NP_NN_GEOMETRY_INVALID_ARGUMENT;
          return 1;
        }
        /* A full-sample NN count can be outside the delete-one domain.
           Keep candidate rejection distinct from a malformed fold owner. */
        if(np_adaptive_fold_lookup_from_scale(
             num_obs_train, vector_scale_factor[dimension],
             &fold_lookup_k, &fold_scale) != 0){
          if(geometry_status != NULL)
            *geometry_status = NP_NN_GEOMETRY_INVALID_SCALE;
          return 1;
        }
        pair_status = compute_nn_adaptive_distance_pair(
          num_obs_train, matrix_train[dimension], fold_lookup_k,
          matrix_bandwidth[dimension],
          geometry_context->adaptive_successor[dimension]);
        if(pair_status != NP_NN_GEOMETRY_OK){
          if(geometry_status != NULL)
            *geometry_status = pair_status;
          return 1;
        }
        geometry_context->adaptive_fold_scale[dimension] = fold_scale;
        if(geometry_context->adaptive_full != NULL){
          const int full_uses_successor = int_nn_k != fold_lookup_k;

          if((!full_uses_successor && nn_scale != 1.0) ||
             (full_uses_successor && int_nn_k != fold_lookup_k + 1)){
            if(geometry_status != NULL)
              *geometry_status = NP_NN_GEOMETRY_INVALID_ARGUMENT;
            return 1;
          }
          for(observation = 0; observation < num_obs_train; ++observation)
            geometry_context->adaptive_full[dimension][observation] =
              nn_scale*(full_uses_successor ?
                geometry_context->adaptive_successor[dimension][observation] :
                matrix_bandwidth[dimension][observation]);
        }
        if(fold_scale != 1.0)
          for(observation = 0; observation < num_obs_train; ++observation){
            matrix_bandwidth[dimension][observation] *= fold_scale;
            geometry_context->adaptive_successor[dimension][observation] *=
              fold_scale;
          }
        continue;
      }
#ifndef MPI2
      {
        const NPNNGeometryStatus query_status =
          np_compute_nn_distance_cached(
            num_obs_train, suppress_parallel, matrix_train[dimension],
            int_nn_k, 1, nn_distance);
        if(query_status != NP_NN_GEOMETRY_OK) {
          if(geometry_status != NULL)
            *geometry_status = query_status;
          return 1;
        }
      }
#else
      {
        const NPNNGeometryStatus query_status =
          compute_nn_distance_status(
            num_obs_train, suppress_parallel, matrix_train[dimension],
            int_nn_k, nn_distance);
        if(query_status != NP_NN_GEOMETRY_OK) {
          if(geometry_status != NULL)
            *geometry_status = query_status;
          return 1;
        }
      }
#endif
      for(observation = 0; observation < num_obs_train; ++observation) {
        const double bandwidth = nn_scale * nn_distance[observation];
        if(!isfinite(bandwidth) || bandwidth <= 0.0) {
          if(geometry_status != NULL)
            *geometry_status = !isfinite(bandwidth) ?
              NP_NN_GEOMETRY_NONFINITE_RADIUS : NP_NN_GEOMETRY_ZERO_RADIUS;
          return 1;
        }
        matrix_bandwidth[dimension][observation] = bandwidth;
      }
    }
  }

  return 0;
}

int np_kernel_bandwidth_continuous_nn(
  int BANDWIDTH,
  int num_obs_train,
  int num_obs_eval,
  int num_cont,
  int suppress_parallel,
  double *vector_scale_factor,
  double **matrix_train,
  double **matrix_eval,
  double **matrix_bandwidth,
  NPNNGeometryStatus *geometry_status)
{
  int dimension;
  int allocation_count;
  int status;
  double *nn_distance;
#ifdef MPI2
  int stride_ignored;
#endif

  if((BANDWIDTH != BW_GEN_NN && BANDWIDTH != BW_ADAP_NN) ||
     num_obs_train <= 0 || num_obs_eval <= 0 || num_cont <= 0 ||
     vector_scale_factor == NULL || matrix_train == NULL ||
     matrix_eval == NULL || matrix_bandwidth == NULL)
    return 1;

  for(dimension = 0; dimension < num_cont; ++dimension)
    if(matrix_train[dimension] == NULL || matrix_eval[dimension] == NULL ||
       matrix_bandwidth[dimension] == NULL)
      return 1;

  allocation_count =
    (BANDWIDTH == BW_GEN_NN) ? num_obs_eval : num_obs_train;
#ifdef MPI2
  if(!np_int_padded_count_nonnegative(allocation_count,
                                      iNum_Processors,
                                      1,
                                      &stride_ignored,
                                      &allocation_count))
    return 1;
#endif

  nn_distance = alloc_vecd(allocation_count);
  status = np_kernel_bandwidth_continuous_nn_into_ctx(
    BANDWIDTH, num_obs_train, num_obs_eval, num_cont, suppress_parallel,
    vector_scale_factor, matrix_train, matrix_eval, matrix_bandwidth,
    nn_distance, NULL, geometry_status);
  free(nn_distance);
  return status;
}


int kernel_bandwidth_mean_ctx(int KERNEL,
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
                          double *vector_lambda,
                          const NPNNGeometryContext *x_geometry_context,
                          const NPNNGeometryContext *y_geometry_context,
                          NPNNGeometryStatus *geometry_status){

/* This computes a matrix of bandwidths for fixed, generalized nearest */
/* neighbor, or adaptive nearest neighbor estimation for a density or */
/* regression function. We permit two sets of continuous variables, X, and */
/* Y, and permit unordered variables for X and Y as well. Finally, the */
/* user has the option of using/writing the `raw' bandwidth (which */
/* approaches zero as n increases) for all variables, or to use a */
/* `normalized' scaling factor which is constant for all sample sizes. */

	int i;
	int status = 0;

	double temp_pow = DBL_MAX;

	double temp_inv;

	double *vec_sdev_x = NULL;
	double *vec_sdev_y = NULL;

	double *nn_distance = NULL;

#ifdef MPI2
	int nn_distance_alloc = 0;
#endif

		if(geometry_status != NULL)
			*geometry_status = NP_NN_GEOMETRY_OK;
		if(num_obs_train == 0) return(1);

#ifdef MPI2
	if((BANDWIDTH == 1) || (BANDWIDTH == 2))
	{
		int stride_ignored = 0;
		const int nn_count = (BANDWIDTH == 1) ? num_obs_eval : num_obs_train;

		if(!np_int_padded_count_nonnegative(nn_count,
		                                    iNum_Processors,
		                                    1,
		                                    &stride_ignored,
		                                    &nn_distance_alloc))
			return(1);
	}
#endif

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
				error("\r ** Fatal Error in routine kernel_bandwidth() ** The variable appears to be constant!");
			}

		}

		for(i=0; i < num_var_cont; i++)
		{

			//vec_sdev_y[i] = standerrd(num_obs_train, matrix_Y_train[i]);
      vec_sdev_y[i] = vector_continuous_stddev_extern[i+num_reg_cont];

			if(vec_sdev_y[i] <= DBL_MIN)
			{
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
		nn_distance = alloc_vecd(nn_distance_alloc);
	}

	if(BANDWIDTH == 2)
	{
		nn_distance = alloc_vecd(nn_distance_alloc);
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
	else if(BANDWIDTH == BW_GEN_NN || BANDWIDTH == BW_ADAP_NN)
	{
			status = np_kernel_bandwidth_continuous_nn_into_ctx(
				BANDWIDTH, num_obs_train, num_obs_eval, num_reg_cont,
				suppress_parallel, vector_scale_factor,
				matrix_X_train, matrix_X_eval, matrix_bandwidth_X,
				nn_distance, x_geometry_context, geometry_status);
		if(status != 0) {
			goto cleanup;
		}

			status = np_kernel_bandwidth_continuous_nn_into_ctx(
				BANDWIDTH, num_obs_train, num_obs_eval, num_var_cont,
				suppress_parallel, vector_scale_factor + num_reg_cont,
				matrix_Y_train, matrix_Y_eval, matrix_bandwidth_Y,
				nn_distance, y_geometry_context, geometry_status);
		if(status != 0) {
			goto cleanup;
		}
	}

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

cleanup:
	free(nn_distance);
	free(vec_sdev_x);
	free(vec_sdev_y);

	return(status);

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
                          double *vector_lambda)
{
	return kernel_bandwidth_mean_ctx(
		KERNEL, BANDWIDTH, num_obs_train, num_obs_eval,
		num_var_cont, num_var_un, num_var_or,
		num_reg_cont, num_reg_un, num_reg_or,
		suppress_parallel, vector_scale_factor,
		matrix_Y_train, matrix_Y_eval,
		matrix_X_train, matrix_X_eval,
		matrix_bandwidth_Y, matrix_bandwidth_X, vector_lambda,
		NULL, NULL, NULL);
}
