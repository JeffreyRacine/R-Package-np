/* Copyright (C) J. Racine, 1995-2001 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <string.h>

#include <R.h>
#include <R_ext/Utils.h>

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
extern MPI_Comm	*comm;
#endif

#define IO_MIN_TRUE  1
#define IO_MIN_FALSE 0

extern int int_DEBUG;
extern int int_VERBOSE;
extern int int_MINIMIZE_IO;
extern int int_TAYLOR;
extern int int_WEIGHTS;

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

extern int np_lp_engine_extern;

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

// so that quantile stuff gives a more sensible multistart message
extern int imstot;

#define NP_QREG_GRID_POINTS 33

static double np_qreg_quantile_objective_scalar(double y)
{
  double quantile[2];

  quantile[0] = 0.0;
  quantile[1] = y;

  return func_con_density_quantile(quantile);
}

static int np_qreg_build_grid_window(
  double *left,
  double *right,
  double *best_x,
  double *best_f)
{
  double grid_x[NP_QREG_GRID_POINTS];
  double grid_f[NP_QREG_GRID_POINTS];
  double span;
  int best_idx;
  int i;

  span = y_max_extern - y_min_extern;
  if(span <= 0.0)
  {
    *left = y_min_extern;
    *right = y_min_extern;
    *best_x = y_min_extern;
    *best_f = np_qreg_quantile_objective_scalar(y_min_extern);
    return 1;
  }

  best_idx = 0;
  for(i = 0; i < NP_QREG_GRID_POINTS; i++)
  {
    grid_x[i] = y_min_extern + span * ((double)i / (double)(NP_QREG_GRID_POINTS - 1));
    grid_f[i] = np_qreg_quantile_objective_scalar(grid_x[i]);
    if(grid_f[i] < grid_f[best_idx])
    {
      best_idx = i;
    }
  }

  *best_x = grid_x[best_idx];
  *best_f = grid_f[best_idx];

  if((best_idx == 0) || (best_idx == NP_QREG_GRID_POINTS - 1))
  {
    if(best_idx == 0)
    {
      *left = grid_x[0];
      *right = grid_x[1];
    }
    else
    {
      *left = grid_x[NP_QREG_GRID_POINTS-2];
      *right = grid_x[NP_QREG_GRID_POINTS-1];
    }
    return 1;
  }

  *left = grid_x[best_idx-1];
  *right = grid_x[best_idx+1];
  return 1;
}

static double np_qreg_refine_golden_1d(
  double left,
  double right,
  double tol,
  double small,
  int itmax,
  double *quantile_best)
{
  const double phi = 0.6180339887498948482;
  double a;
  double b;
  double c;
  double d;
  double fa;
  double fb;
  double fc;
  double fd;
  int iter;
  int maxiter;

  a = left;
  b = right;
  if(b < a)
  {
    double tmp = a;
    a = b;
    b = tmp;
  }

  if((b - a) <= small)
  {
    fa = np_qreg_quantile_objective_scalar(a);
    fb = np_qreg_quantile_objective_scalar(b);
    if(fa <= fb)
    {
      *quantile_best = a;
      return fa;
    }
    *quantile_best = b;
    return fb;
  }

  c = b - phi * (b - a);
  d = a + phi * (b - a);
  fc = np_qreg_quantile_objective_scalar(c);
  fd = np_qreg_quantile_objective_scalar(d);

  maxiter = (itmax < 64) ? itmax : 64;
  for(iter = 0; iter < maxiter; iter++)
  {
    if(fabs(b - a) <= tol * (fabs(c) + fabs(d)) + small)
    {
      break;
    }

    if(fc <= fd)
    {
      b = d;
      d = c;
      fd = fc;
      c = b - phi * (b - a);
      fc = np_qreg_quantile_objective_scalar(c);
    }
    else
    {
      a = c;
      c = d;
      fc = fd;
      d = a + phi * (b - a);
      fd = np_qreg_quantile_objective_scalar(d);
    }
  }

  fa = np_qreg_quantile_objective_scalar(a);
  fb = np_qreg_quantile_objective_scalar(b);

  *quantile_best = a;
  if(fc < fa)
  {
    *quantile_best = c;
    fa = fc;
  }
  if(fd < fa)
  {
    *quantile_best = d;
    fa = fd;
  }
  if(fb < fa)
  {
    *quantile_best = b;
    fa = fb;
  }

  return fa;
}

static double np_qreg_extract_quantile_1d(
  double tol,
  double small,
  int itmax,
  double *quantile_best)
{
  double left;
  double right;
  double xmin;
  double fret;
  double grid_best_x;
  double grid_best_f;
  int window_status;

  window_status = np_qreg_build_grid_window(&left, &right, &grid_best_x, &grid_best_f);

  if(window_status == 0)
  {
    error("C_np_quantile_conditional: canonical one-dimensional quantile extraction failed to bracket a finite support window");
  }

  fret = np_qreg_refine_golden_1d(left, right, tol, small, itmax, &xmin);

  if((!R_FINITE(fret)) ||
     (!R_FINITE(xmin)) ||
     (xmin < y_min_extern) ||
     (xmin > y_max_extern))
  {
    error("C_np_quantile_conditional: canonical one-dimensional quantile extraction failed to produce a finite in-support candidate");
  }

  if(grid_best_f < fret)
  {
    *quantile_best = grid_best_x;
    return grid_best_f;
  }

  *quantile_best = xmin;

  return fret;
}

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
    0, // do not suppress_parallel
		vector_scale_factor,
		matrix_Y_continuous_train,
		matrix_Y_continuous_eval,
		matrix_X_continuous_train,
		matrix_X_continuous_eval,
		matrix_bandwidth_var,
		matrix_bandwidth_reg,
		lambda) == 1)
	{
#ifdef MPI2
		MPI_Barrier(comm[1]);
		MPI_Finalize();
#endif
		error("\n** Error: invalid bandwidth.");
	}

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=0; j < num_obs_eval; j++)
		{

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

			cdf[j] = sum_ker/NZD(sum_ker_marginal);
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs_eval; j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = 0.0;

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

			cdf[j] = sum_ker/NZD(sum_ker_marginal);
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}
	else
	{

		for(j=0; j < num_obs_eval; j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = 0.0;

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

			cdf[j] = sum_ker/NZD(sum_ker_marginal);
			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

		}

	}

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
double zero,
double lbc_dir,
int dfc_dir,
double c_dir,
double initc_dir,
double lbd_dir,
double  hbd_dir,
double  d_dir,
double  initd_dir)
{

	int i;
	int j;
	int k;
	double quantile[2];
	double **matrix_y;

	double quantile_l;
	double quantile_u;

	double *lambda = NULL;
	double **matrix_bandwidth_var = NULL;
	double **matrix_bandwidth_reg = NULL;

	#ifdef MPI2
	int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
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
      0, // do not suppress_parallel
			vector_scale_factor,
			matrix_Y_continuous_train,
			matrix_Y_continuous_train, /* Same Y for training and evaluation */
			matrix_X_continuous_train,
			matrix_X_continuous_eval,
			matrix_bandwidth_var,
			matrix_bandwidth_reg,
			lambda) == 1)
		{
#ifdef MPI2
			MPI_Barrier(comm[1]);
			MPI_Finalize();
#endif
      error("\n** Error: invalid bandwidth.");
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

			(void) np_qreg_extract_quantile_1d(tol,
	                                     small,
	                                     itmax,
	                                     &quantile[1]);

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

					(void) np_qreg_extract_quantile_1d(tol,
	                                         small,
	                                         itmax,
	                                         &quantile[1]);

				quantile_l = quantile[1];

				/* Gradient for continuous regressors - quantile evaluated at x+h */

				for(j = 0; j < num_reg_continuous; j++)
				{
					matrix_X_continuous_quantile_extern[j][0] = matrix_X_continuous_eval[j][i];
				}

				matrix_X_continuous_quantile_extern[k][0] = matrix_X_continuous_eval[k][i] + matrix_bandwidth_reg[k][i]/2.0;

					(void) np_qreg_extract_quantile_1d(tol,
	                                         small,
	                                         itmax,
	                                         &quantile[1]);

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

			(void) np_qreg_extract_quantile_1d(tol,
	                                     small,
	                                     itmax,
	                                     &quantile[1]);

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

					(void) np_qreg_extract_quantile_1d(tol,
	                                         small,
	                                         itmax,
	                                         &quantile[1]);

				quantile_l = quantile[1];

				/* Gradient for continuous regressors - quantile evaluated at x+h */

				for(j = 0; j < num_reg_continuous; j++)
				{
					matrix_X_continuous_quantile_extern[j][0] = matrix_X_continuous_eval[j][i];
				}

				matrix_X_continuous_quantile_extern[k][0] = matrix_X_continuous_eval[k][i] + matrix_bandwidth_reg[k][i]/2.0;

					(void) np_qreg_extract_quantile_1d(tol,
	                                         small,
	                                         itmax,
	                                         &quantile[1]);

				quantile_u = quantile[1];

				quan_gradient[k][i-my_rank*stride] = (quantile_u-quantile_l)/(2.0*matrix_bandwidth_reg[k][i]);

			}

		}

		/* End gradient */

	}

	/* Collect */

	MPI_Gather(quan, stride, MPI_DOUBLE, quan, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(quan, num_obs_eval, MPI_DOUBLE, 0, comm[1]);

	MPI_Gather(quan_stderr, stride, MPI_DOUBLE, quan_stderr, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(quan_stderr, num_obs_eval, MPI_DOUBLE, 0, comm[1]);

	if(gradient_compute == 1)
	{

		for(k = 0; k < num_reg_continuous; k++)
		{

			MPI_Gather(&quan_gradient[k][0], stride, MPI_DOUBLE, &quan_gradient[k][0], stride, MPI_DOUBLE, 0, comm[1]);
			MPI_Bcast(&quan_gradient[k][0], num_obs_eval, MPI_DOUBLE, 0, comm[1]);

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
