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
	int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
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
#ifdef MPI2
		MPI_Barrier(comm[1]);
		MPI_Finalize();
#endif
    error("\n** Error: invalid bandwidth.");
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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_density_categorical()");
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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_density_categorical()");
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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_density_categorical()");
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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_density_categorical()");
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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_density_categorical()");
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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_density_categorical()");
				}
			}

		}

	}

	MPI_Gather(pdf, stride, MPI_DOUBLE, pdf, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(pdf, num_obs_eval, MPI_DOUBLE, 0, comm[1]);

	MPI_Gather(pdf_stderr, stride, MPI_DOUBLE, pdf_stderr, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(pdf_stderr, num_obs_eval, MPI_DOUBLE, 0, comm[1]);

	MPI_Reduce(&log_likelihood_MPI, log_likelihood, 1, MPI_DOUBLE, MPI_SUM, 0, comm[1]);
	MPI_Bcast(log_likelihood, 1, MPI_DOUBLE, 0, comm[1]);
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

	double *p_xj1;
	double *p_xi1;
	double *p_xj2;
	double *p_xi2;

	#ifdef MPI2
	double cv_MPI;
	int stride = (int)ceil((double) num_obs / (double) iNum_Processors);
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
    0, // do not suppress_parallel
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

				if(!(pdf > DBL_MIN) && (int_VERBOSE == 1))
				{
					REprintf("\r                                                                           ");
					REprintf("\r** Guarded CVML contribution in kernel_estimate_density_categorical_leave_one_out_cv()");
				}
				*cv += np_guarded_cvml_contribution(pdf);

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

				if(!(pdf > DBL_MIN) && (int_VERBOSE == 1))
				{
					REprintf("\r                                                                           ");
					REprintf("\r** Guarded CVML contribution in kernel_estimate_density_categorical_leave_one_out_cv()");
				}
				*cv += np_guarded_cvml_contribution(pdf);

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

				if(!(pdf > DBL_MIN) && (int_VERBOSE == 1))
				{
					REprintf("\r                                                                           ");
					REprintf("\r** Guarded CVML contribution in kernel_estimate_density_categorical_leave_one_out_cv()");
				}
				*cv += np_guarded_cvml_contribution(pdf);

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

			if(!(pdf > DBL_MIN) && (int_VERBOSE == 1))
			{
				REprintf("\r                                                                           ");
				REprintf("\r** Guarded CVML contribution in kernel_estimate_density_categorical_leave_one_out_cv()");
			}
			*cv += np_guarded_cvml_contribution(pdf);

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

			if(!(pdf > DBL_MIN) && (int_VERBOSE == 1))
			{
				REprintf("\r                                                                           ");
				REprintf("\r** Guarded CVML contribution in kernel_estimate_density_categorical_leave_one_out_cv()");
			}
			*cv += np_guarded_cvml_contribution(pdf);

		}
	}


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

				if(!(pdf > DBL_MIN) && (int_VERBOSE == 1) && (my_rank == 0))
				{
					REprintf("\r                                                                           ");
					REprintf("\r** Guarded CVML contribution in kernel_estimate_density_categorical_leave_one_out_cv()");
				}
				cv_MPI += np_guarded_cvml_contribution(pdf);

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

				if(!(pdf > DBL_MIN) && (int_VERBOSE == 1) && (my_rank == 0))
				{
					REprintf("\r                                                                           ");
					REprintf("\r** Guarded CVML contribution in kernel_estimate_density_categorical_leave_one_out_cv()");
				}
				cv_MPI += np_guarded_cvml_contribution(pdf);

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

				if(!(pdf > DBL_MIN) && (int_VERBOSE == 1) && (my_rank == 0))
				{
					REprintf("\r                                                                           ");
					REprintf("\r** Guarded CVML contribution in kernel_estimate_density_categorical_leave_one_out_cv()");
				}
				cv_MPI += np_guarded_cvml_contribution(pdf);

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

			if(!(pdf > DBL_MIN) && (int_VERBOSE == 1) && (my_rank == 0))
			{
				REprintf("\r                                                                           ");
				REprintf("\r** Guarded CVML contribution in kernel_estimate_density_categorical_leave_one_out_cv()");
			}
			cv_MPI += np_guarded_cvml_contribution(pdf);

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

			if(!(pdf > DBL_MIN) && (int_VERBOSE == 1) && (my_rank == 0))
			{
				REprintf("\r                                                                           ");
				REprintf("\r** Guarded CVML contribution in kernel_estimate_density_categorical_leave_one_out_cv()");
			}
			cv_MPI += np_guarded_cvml_contribution(pdf);

		}
	}

	/* Now reduce */

	MPI_Reduce(&cv_MPI, cv, 1, MPI_DOUBLE, MPI_SUM, 0, comm[1]);
	MPI_Bcast(cv, 1, MPI_DOUBLE, 0, comm[1]);
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

	#ifdef MPI2
	double cv_MPI;
	int stride = (int)ceil((double) num_obs / (double) iNum_Processors);
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
    0, // do not suppress_parallel
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
			sum_ker_marginal = 0.0;

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

      pdf = sum_ker/(prod_h*NZD(sum_ker_marginal));

			if(!(pdf > DBL_MIN) && (int_VERBOSE == 1))
			{
				REprintf("\r                                                                           ");
				REprintf("\r** Guarded CVML contribution in kernel_estimate_con_density_categorical_leave_one_out_cv()");
			}
			*cv += np_guarded_cvml_contribution(pdf);

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
			sum_ker_marginal = 0.0;

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

			pdf = sum_ker/(prod_h*NZD(sum_ker_marginal));

			if(!(pdf > DBL_MIN) && (int_VERBOSE == 1))
			{
				REprintf("\r                                                                           ");
				REprintf("\r** Guarded CVML contribution in kernel_estimate_con_density_categorical_leave_one_out_cv()");
			}
			*cv += np_guarded_cvml_contribution(pdf);

		}

	}
	else
	{

		for(j=0; j < num_obs; j++)
		{
		  R_CheckUserInterrupt();
			sum_ker = 0.0;
			sum_ker_marginal = 0.0;

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

			pdf = sum_ker/NZD(sum_ker_marginal);

			if(!(pdf > DBL_MIN) && (int_VERBOSE == 1))
			{
				REprintf("\r                                                                           ");
				REprintf("\r** Guarded CVML contribution in kernel_estimate_con_density_categorical_leave_one_out_cv()");
			}
			*cv += np_guarded_cvml_contribution(pdf);

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
			sum_ker_marginal = 0.0;

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

      pdf = sum_ker/(prod_h*NZD(sum_ker_marginal));

			if(!(pdf > DBL_MIN) && (int_VERBOSE == 1) && (my_rank == 0))
			{
				REprintf("\r                                                                           ");
				REprintf("\r** Guarded CVML contribution in kernel_estimate_con_density_categorical_leave_one_out_cv()");
			}
			cv_MPI += np_guarded_cvml_contribution(pdf);

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
			sum_ker_marginal = 0.0;

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

			pdf = sum_ker/(prod_h*NZD(sum_ker_marginal));

			if(!(pdf > DBL_MIN) && (int_VERBOSE == 1) && (my_rank == 0))
			{
				REprintf("\r                                                                           ");
				REprintf("\r** Guarded CVML contribution in kernel_estimate_con_density_categorical_leave_one_out_cv()");
			}
			cv_MPI += np_guarded_cvml_contribution(pdf);

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = 0.0;
			sum_ker_marginal = 0.0;

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

			pdf = sum_ker/NZD(sum_ker_marginal);

			if(!(pdf > DBL_MIN) && (int_VERBOSE == 1) && (my_rank == 0))
			{
				REprintf("\r                                                                           ");
				REprintf("\r** Guarded CVML contribution in kernel_estimate_con_density_categorical_leave_one_out_cv()");
			}
			cv_MPI += np_guarded_cvml_contribution(pdf);

		}
	}

	/* Now reduce */

	cv_MPI /= (double) num_obs;
	MPI_Reduce(&cv_MPI, cv, 1, MPI_DOUBLE, MPI_SUM, 0, comm[1]);
	MPI_Bcast(cv, 1, MPI_DOUBLE, 0, comm[1]);
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
	int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
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
#ifdef MPI2
		MPI_Barrier(comm[1]);
		MPI_Finalize();
#endif
    error("\n** Error: invalid bandwidth.");
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

	MPI_Gather(cdf, stride, MPI_DOUBLE, cdf, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(cdf, num_obs_eval, MPI_DOUBLE, 0, comm[1]);
	MPI_Gather(cdf_stderr, stride, MPI_DOUBLE, cdf_stderr, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(cdf_stderr, num_obs_eval, MPI_DOUBLE, 0, comm[1]);
	#endif

	free(lambda);

	free_mat(matrix_bandwidth,num_reg_continuous);
	free_mat(matrix_bandwidth_deriv,num_reg_continuous);

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
	int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
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

			pdf[j] = sum_ker/(prod_h*NZD(sum_ker_marginal));

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
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

			pdf[j] = sum_ker/(prod_h*NZD(sum_ker_marginal));

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
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
			sum_ker_marginal = 0.0;

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

      /* Don't keep calling NZD for successive divides */

      sum_ker_marginal = NZD(sum_ker_marginal);

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
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

			pdf[j-my_rank*stride] = sum_ker/(prod_h*NZD(sum_ker_marginal));

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
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

			pdf[j-my_rank*stride] = sum_ker/(prod_h*NZD(sum_ker_marginal));

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_ker = sum_ker_temp = 0.0;
			sum_ker_marginal = 0.0;

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

      /* Don't keep calling NZD for successive divides */

      sum_ker_marginal = NZD(sum_ker_marginal);

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

		}

	}

	MPI_Gather(pdf, stride, MPI_DOUBLE, pdf, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(pdf, num_obs_eval, MPI_DOUBLE, 0, comm[1]);

	MPI_Gather(pdf_stderr, stride, MPI_DOUBLE, pdf_stderr, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(pdf_stderr, num_obs_eval, MPI_DOUBLE, 0, comm[1]);

	MPI_Reduce(&log_likelihood_MPI, log_likelihood, 1, MPI_DOUBLE, MPI_SUM, 0, comm[1]);
	MPI_Bcast(log_likelihood, 1, MPI_DOUBLE, 0, comm[1]);
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
	int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
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

	#ifndef MPI2

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
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
		  R_CheckUserInterrupt();
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
		  R_CheckUserInterrupt();
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
	#endif

	#ifdef MPI2

	/* Conduct the estimation */

	if(BANDWIDTH_den == 0)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
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

			cdf[j-my_rank*stride] = sum_ker/NZD(sum_ker_marginal);
			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
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

			cdf[j-my_rank*stride] = sum_ker/NZD(sum_ker_marginal);
			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
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

			cdf[j-my_rank*stride] = sum_ker/NZD(sum_ker_marginal);
			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

		}

	}

	MPI_Gather(cdf, stride, MPI_DOUBLE, cdf, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(cdf, num_obs_eval, MPI_DOUBLE, 0, comm[1]);
	MPI_Gather(cdf_stderr, stride, MPI_DOUBLE, cdf_stderr, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(cdf_stderr, num_obs_eval, MPI_DOUBLE, 0, comm[1]);
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
	int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
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

				if(i != j) {
          sum_ker += prod_kernel_cont*prod_kernel_cat;
          sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;
        }

			}

			cdf[j] = sum_ker/NZD(sum_ker_marginal);

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
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

      if(i != j) {
        cdf[j] = sum_ker/NZD(sum_ker_marginal);
      }

		}

	}
	else
	{

		for(j=0; j < num_obs_eval; j++)
		{
		  R_CheckUserInterrupt();
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

				if(i != j) {
          sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
          sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h;
        }

			}

			cdf[j] = sum_ker/NZD(sum_ker_marginal);

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

				if(i != j) {
          sum_ker += prod_kernel_cont*prod_kernel_cat;
          sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;
        }

			}

			cdf[j-my_rank*stride] = sum_ker/NZD(sum_ker_marginal);

		}

	}
	else if(BANDWIDTH_den == 1)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
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

				if(i != j) {
          sum_ker += prod_kernel_cont*prod_kernel_cat;
          sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat;
        }

			}

			cdf[j-my_rank*stride] = sum_ker/NZD(sum_ker_marginal);

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
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

				if(i != j) {
          sum_ker += prod_kernel_cont*prod_kernel_cat/prod_h;
          sum_ker_marginal += prod_kernel_marginal_cont*prod_kernel_marginal_cat/prod_h;
        }

			}

			cdf[j-my_rank*stride] = sum_ker/NZD(sum_ker_marginal);

		}

	}

	MPI_Gather(cdf, stride, MPI_DOUBLE, cdf, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(cdf, num_obs_eval, MPI_DOUBLE, 0, comm[1]);

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

	#ifdef MPI2
	int stride = (int)ceil((double) num_obs_train / (double) iNum_Processors);
	if(stride < 1) stride = 1;
	#endif

  /* Must declare and free cdf_eval */

	#ifndef MPI2
  cdf_loo = alloc_vecd(num_obs_train);

	matrix_Y_unordered_eval = alloc_matd(num_obs_train, num_var_unordered);
	matrix_Y_ordered_eval = alloc_matd(num_obs_train, num_var_ordered);
	matrix_Y_continuous_eval = alloc_matd(num_obs_train, num_var_continuous);
	#endif

	#ifdef MPI2
  cdf_loo = alloc_vecd(stride*iNum_Processors);

	matrix_Y_unordered_eval = alloc_matd(stride*iNum_Processors, num_var_unordered);
	matrix_Y_ordered_eval = alloc_matd(stride*iNum_Processors, num_var_ordered);
	matrix_Y_continuous_eval = alloc_matd(stride*iNum_Processors, num_var_continuous);
	#endif

  /* XXX */

    /*    for(l=0; l < num_var_unordered; l++)
      for(j=0; j < num_obs_train; j++)
        matrix_Y_unordered_eval[l][j] = matrix_Y_unordered_train[l][i];

    for(l=0; l < num_var_ordered; l++)
      for(j=0; j < num_obs_train; j++)
	      matrix_Y_ordered_eval[l][j] = matrix_Y_ordered_train[l][i];*/

  for(i=0; i < num_obs_train; i++) {

    /* Brute force `copy' Y eval could of course be improved via
       pointer */

    /* Require nested loops for ordered, unordered, and continuous -
       not implemented yet (Feb 7 2010) */

    for(l=0; l < num_var_continuous; l++)
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
      for(l=0; l < num_var_continuous; l++) {
        indicator *= indfunc(matrix_Y_continuous_train[l][j]-matrix_Y_continuous_eval[l][j]);
      }
      *cv += ipow(indicator - cdf_loo[j],2);
    }

  }

  /* Sum over all variables - need to be careful about denominator */
  /* Note - we are only set up for 1 continuous Y at the moment */

  *cv /= (double) ipow(num_obs_train,1+num_reg_continuous);

  free(cdf_loo);
	free_mat(matrix_Y_unordered_eval, num_var_unordered);
	free_mat(matrix_Y_ordered_eval, num_var_ordered);
	free_mat(matrix_Y_continuous_eval, num_var_continuous);

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
	int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
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
			sum_ker_marginal = 0.0;

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

      sum_ker_marginal = NZD(sum_ker_marginal);

      pdf[j] = sum_ker/(prod_h*sum_ker_marginal);

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical_gradient()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

        pdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
          /(prod_h*sum_ker_marginal*matrix_bandwidth_reg[l][0]);
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
			sum_ker_marginal = 0.0;

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

      sum_ker_marginal = NZD(sum_ker_marginal);

      pdf[j] = sum_ker/(prod_h*sum_ker_marginal);

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical_gradient()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

        pdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
          /(prod_h*sum_ker_marginal*matrix_bandwidth_reg[l][j]);

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
			sum_ker_marginal = 0.0;

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

      sum_ker_marginal = NZD(sum_ker_marginal);

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */
			/* 11/8/01 - removed prod_h from denom */

			for(l = 0; l < num_reg_continuous; l++)
			{

        pdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])/sum_ker_marginal;


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
			sum_ker_marginal = 0.0;

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

			sum_ker_marginal = 	NZD(sum_ker_marginal);

      pdf[j-my_rank*stride] = sum_ker/(prod_h*sum_ker_marginal);

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical_gradient()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

        pdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
        /(prod_h*sum_ker_marginal*matrix_bandwidth_reg[l][0]);

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
			sum_ker_marginal = 0.0;

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

			sum_ker_marginal = 	NZD(sum_ker_marginal);

      pdf[j-my_rank*stride] = sum_ker/(prod_h*sum_ker_marginal);

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical_gradient()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

        pdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
          /(prod_h*sum_ker_marginal*matrix_bandwidth_reg[l][j]);

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
			sum_ker_marginal = 0.0;

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

			sum_ker_marginal = 	NZD(sum_ker_marginal);

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
					REprintf("\r                                                                           ");
					REprintf("\r** Trimming binding in kernel_estimate_con_density_categorical()");
				}
			}

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

        pdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
          /sum_ker_marginal;

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */
				/* 11/26/01 - dropped prod of bws in denom, need to revisit */

				pdf_deriv_stderr[l][j-my_rank*stride] = sqrt(DIFF_KER_PPM / sum_ker_marginal);
			}

		}

	}

	MPI_Gather(pdf, stride, MPI_DOUBLE, pdf, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(pdf, num_obs_eval, MPI_DOUBLE, 0, comm[1]);
	MPI_Gather(pdf_stderr, stride, MPI_DOUBLE, pdf_stderr, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(pdf_stderr, num_obs_eval, MPI_DOUBLE, 0, comm[1]);

	for(l = 0; l < num_reg_continuous; l++)
	{

		MPI_Gather(&pdf_deriv[l][0], stride, MPI_DOUBLE, &pdf_deriv[l][0], stride, MPI_DOUBLE, 0, comm[1]);
		MPI_Bcast(&pdf_deriv[l][0], num_obs_eval, MPI_DOUBLE, 0, comm[1]);
		MPI_Gather(&pdf_deriv_stderr[l][0], stride, MPI_DOUBLE, &pdf_deriv_stderr[l][0], stride, MPI_DOUBLE, 0, comm[1]);
		MPI_Bcast(&pdf_deriv_stderr[l][0], num_obs_eval, MPI_DOUBLE, 0, comm[1]);

	}

	MPI_Reduce(&log_likelihood_MPI, log_likelihood, 1, MPI_DOUBLE, MPI_SUM, 0, comm[1]);
	MPI_Bcast(log_likelihood, 1, MPI_DOUBLE, 0, comm[1]);
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

#ifndef MPI2
	double *pointer_m;
	double *pointer_me;
	double *pointer_g;
#endif

	#ifdef MPI2
	int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
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

		MPI_Gather(&pdf_deriv[l][0], stride, MPI_DOUBLE, &pdf_deriv[l][0], stride, MPI_DOUBLE, 0, comm[1]);
		MPI_Bcast(&pdf_deriv[l][0], num_obs_eval, MPI_DOUBLE, 0, comm[1]);
		MPI_Gather(&pdf_deriv_stderr[l][0], stride, MPI_DOUBLE, &pdf_deriv_stderr[l][0], stride, MPI_DOUBLE, 0, comm[1]);
		MPI_Bcast(&pdf_deriv_stderr[l][0], num_obs_eval, MPI_DOUBLE, 0, comm[1]);

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
	int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
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
			sum_ker_marginal = 0.0;

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

			sum_ker_marginal = 	NZD(sum_ker_marginal);

      cdf[j] = sum_ker/sum_ker_marginal;

			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

        cdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
          /(sum_ker_marginal*matrix_bandwidth_reg[l][0]);

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
			sum_ker_marginal = 0.0;

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

			sum_ker_marginal = 	NZD(sum_ker_marginal);

      cdf[j] = sum_ker/sum_ker_marginal;

			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

        cdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
          /(sum_ker_marginal*matrix_bandwidth_reg[l][j]);

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
			sum_ker_marginal = 0.0;

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

			sum_ker_marginal = 	NZD(sum_ker_marginal);

      cdf[j] = sum_ker/sum_ker_marginal;

			cdf_stderr[j] = sqrt(cdf[j]*(1.0-cdf[j])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

        cdf_deriv[l][j] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])/sum_ker_marginal;

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
			sum_ker_marginal = 0.0;

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

			sum_ker_marginal = 	NZD(sum_ker_marginal);

      cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;

			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

        cdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
          /(sum_ker_marginal*matrix_bandwidth_reg[l][0]);

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
			sum_ker_marginal = 0.0;

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

			sum_ker_marginal = 	NZD(sum_ker_marginal);

      cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;

			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			sum_ker_marginal = 	NZD(sum_ker_marginal);

			for(l = 0; l < num_reg_continuous; l++)
			{

        cdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])
          /(sum_ker_marginal*matrix_bandwidth_reg[l][j]);

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
			sum_ker_marginal = 0.0;

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

			sum_ker_marginal = 	NZD(sum_ker_marginal);

      cdf[j-my_rank*stride] = sum_ker/sum_ker_marginal;


			cdf_stderr[j-my_rank*stride] = sqrt(cdf[j-my_rank*stride]*(1.0-cdf[j-my_rank*stride])/(double)num_obs_train);

			/* gradient[0][] is that for _first_ continuous variable */

			for(l = 0; l < num_reg_continuous; l++)
			{

        cdf_deriv[l][j-my_rank*stride] = (sum_ker_deriv[l]-(sum_ker/sum_ker_marginal)*sum_ker_marginal_deriv[l])/sum_ker_marginal;

				/* grads definitely incorrect... dropped tmp_var after sqrt(, and using formula for regression */
				/* 11/28/01 - removed h^2 for adaptive due to alloc issues */

				cdf_deriv_stderr[l][j-my_rank*stride] = sqrt(DIFF_KER_PPM / sum_ker_marginal);

			}

		}

	}

	MPI_Gather(cdf, stride, MPI_DOUBLE, cdf, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(cdf, num_obs_eval, MPI_DOUBLE, 0, comm[1]);
	MPI_Gather(cdf_stderr, stride, MPI_DOUBLE, cdf_stderr, stride, MPI_DOUBLE, 0, comm[1]);
	MPI_Bcast(cdf_stderr, num_obs_eval, MPI_DOUBLE, 0, comm[1]);

	for(l = 0; l < num_reg_continuous; l++)
	{

		MPI_Gather(&cdf_deriv[l][0], stride, MPI_DOUBLE, &cdf_deriv[l][0], stride, MPI_DOUBLE, 0, comm[1]);
		MPI_Bcast(&cdf_deriv[l][0], num_obs_eval, MPI_DOUBLE, 0, comm[1]);
		MPI_Gather(&cdf_deriv_stderr[l][0], stride, MPI_DOUBLE, &cdf_deriv_stderr[l][0], stride, MPI_DOUBLE, 0, comm[1]);
		MPI_Bcast(&cdf_deriv_stderr[l][0], num_obs_eval, MPI_DOUBLE, 0, comm[1]);

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

#ifndef MPI2
	double *pointer_m;
	double *pointer_me;
	double *pointer_g;
#endif

	#ifdef MPI2
	int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
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

		MPI_Gather(&cdf_deriv[l][0], stride, MPI_DOUBLE, &cdf_deriv[l][0], stride, MPI_DOUBLE, 0, comm[1]);
		MPI_Bcast(&cdf_deriv[l][0], num_obs_eval, MPI_DOUBLE, 0, comm[1]);
		MPI_Gather(&cdf_deriv_stderr[l][0], stride, MPI_DOUBLE, &cdf_deriv_stderr[l][0], stride, MPI_DOUBLE, 0, comm[1]);
		MPI_Bcast(&cdf_deriv_stderr[l][0], num_obs_eval, MPI_DOUBLE, 0, comm[1]);

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
	int stride = (int)ceil((double) num_obs / (double) iNum_Processors);
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
    0, // do not suppress_parallel
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

	MPI_Reduce(&cv_MPI, cv, 1, MPI_DOUBLE, MPI_SUM, 0, comm[1]);
	MPI_Bcast(cv, 1, MPI_DOUBLE, 0, comm[1]);
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
	int stride = (int)ceil((double) num_obs / (double) iNum_Processors);
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
    0, // do not suppress_parallel
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
				sum_ker_marginal = 0.0;

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

        sum_ker_marginal = 	NZD(sum_ker_marginal);

        *cv += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
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
				sum_ker_marginal = 0.0;

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

        sum_ker_marginal = 	NZD(sum_ker_marginal);

        *cv += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;

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
				sum_ker_marginal = 0.0;

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

        sum_ker_marginal = 	NZD(sum_ker_marginal);

        *cv += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;

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
			sum_ker_marginal = 0.0;

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

      sum_ker_marginal = 	NZD(sum_ker_marginal);

      *cv += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;

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
				sum_ker_marginal = 0.0;

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

        sum_ker_marginal = 	NZD(sum_ker_marginal);

        cv_MPI += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
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
				sum_ker_marginal = 0.0;

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

        sum_ker_marginal = 	NZD(sum_ker_marginal);

        cv_MPI += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
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
				sum_ker_marginal = 0.0;

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

        sum_ker_marginal = 	NZD(sum_ker_marginal);

        cv_MPI += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
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
			sum_ker_marginal = 0.0;

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

      sum_ker_marginal = 	NZD(sum_ker_marginal);

      cv_MPI += (sum_ker_convol/sum_ker_marginal-2.0*sum_ker)/sum_ker_marginal;
		}

		/* Don't forget!!! */

		cv_MPI /= (double) num_obs;

		free_mat(matrix_weights_K_x,stride*iNum_Processors);
		free_mat(matrix_weights_K_xy,stride*iNum_Processors);
		free_mat(matrix_weights_K_convol_y,stride*iNum_Processors);

	}
	MPI_Reduce(&cv_MPI, cv, 1, MPI_DOUBLE, MPI_SUM, 0, comm[1]);
	MPI_Bcast(cv, 1, MPI_DOUBLE, 0, comm[1]);
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
