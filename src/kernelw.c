/* This module contains the functions for computing the kernel weights. */

/* Copyright (C) J. Racine, 1995-2000 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "headers.h"

#ifdef RCSID
static char rcsid[] = "$Id: kernelw.c,v 1.4 2006/11/02 16:56:49 tristen Exp $";
#endif

int kernel_weights_conditional_convolution_cv(
int int_WEIGHTS,
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
double *lambda,
double **matrix_bandwidth_var,
double **matrix_bandwidth_reg,
int *num_categories,
double **matrix_categorical_vals,
double **matrix_weights_K_x,
double **matrix_weights_K_xy,
double **matrix_weights_K_convol_y)
{

/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel_convol_cat;
	double prod_kernel_convol_cont;

	double prod_kernel_cat;
	double prod_kernel_cont;

	double prod_kernel_marginal_cat;
	double prod_kernel_marginal_cont;

	double *pointer_k_x;
	double *pointer_k_xy;
	double *pointer_k_convol_y;

	double K_xy_jj;
	double K_x_jj;
	double K_convol_y_jj;

/* Only compute if weights requested */

	if(int_WEIGHTS == 0)
	{
		return(0);
	}

	if(BANDWIDTH_den == 0)
	{

/* Symmetry */

/* i == j - don't compute inside of loop */

		prod_kernel_cont = 1.0;
		prod_kernel_convol_cont = 1.0;

		for(l = 0; l < num_reg_continuous; l++)
		{
			prod_kernel_cont *= kernel(KERNEL_reg, 0.0/matrix_bandwidth_reg[l][0])/matrix_bandwidth_reg[l][0];
		}

		prod_kernel_marginal_cont = prod_kernel_cont;

		for(l = 0; l < num_var_continuous; l++)
		{
			prod_kernel_cont *= kernel(KERNEL_den, 0.0/matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];
			prod_kernel_convol_cont *= kernel_convol(KERNEL_den,BANDWIDTH_den,
				0.0/matrix_bandwidth_var[l][0],matrix_bandwidth_var[l][0],matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];
		}

		prod_kernel_cat = 1.0;
		prod_kernel_convol_cat = 1.0;

		for(l = 0; l < num_reg_unordered; l++)
		{
			prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, 0.0,0.0,lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
		}

		for(l = 0; l < num_reg_ordered; l++)
		{
			prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg,0.0,0.0,lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
		}

		prod_kernel_marginal_cat = prod_kernel_cat;

		for(l = 0; l < num_var_unordered; l++)
		{
			prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, 0.0,0.0,lambda[l],num_categories[l]);
			prod_kernel_convol_cat *= kernel_unordered_convolution(KERNEL_unordered_den, 0.0,0.0,lambda[l], num_categories[l], matrix_categorical_vals[l]);
		}

		for(l = 0; l < num_var_ordered; l++)
		{
			prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, 0.0,0.0,lambda[l]);
			prod_kernel_convol_cat *= kernel_ordered_convolution(KERNEL_ordered_den, 0.0,0.0,lambda[l+num_var_unordered], num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
		}

		K_xy_jj = prod_kernel_cont*prod_kernel_cat;
		K_x_jj = prod_kernel_marginal_cont*prod_kernel_marginal_cat;
		K_convol_y_jj = prod_kernel_convol_cont*prod_kernel_convol_cat;

		for(j=0; j < num_obs; j++)
		{

/* Diagonal */

			matrix_weights_K_xy[j][j] = K_xy_jj;
			matrix_weights_K_x[j][j] = K_x_jj;
			matrix_weights_K_convol_y[j][j] = K_convol_y_jj;

			pointer_k_xy = &matrix_weights_K_xy[j][0];
			pointer_k_x = &matrix_weights_K_x[j][0];
			pointer_k_convol_y = &matrix_weights_K_convol_y[j][0];

			for (i =  0; i <  j; i++)
			{

				prod_kernel_cont = 1.0;
				prod_kernel_convol_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][0])/matrix_bandwidth_reg[l][0];
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];
					prod_kernel_convol_cont *= kernel_convol(KERNEL_den,BANDWIDTH_den,
						(matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][0],matrix_bandwidth_var[l][0],matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];
				}

				prod_kernel_cat = 1.0;
				prod_kernel_convol_cat = 1.0;

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
					prod_kernel_convol_cat *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l], num_categories[l], matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
					prod_kernel_convol_cat *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered], num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
				}

				matrix_weights_K_xy[j][i] = matrix_weights_K_xy[i][j] = prod_kernel_cont*prod_kernel_cat;
				matrix_weights_K_x[j][i] = matrix_weights_K_x[i][j] = prod_kernel_marginal_cont*prod_kernel_marginal_cat;
				matrix_weights_K_convol_y[j][i] = matrix_weights_K_convol_y[i][j] = prod_kernel_convol_cont*prod_kernel_convol_cat;

			}

		}

	} else if(BANDWIDTH_den == 1)
	{

		for(j=0; j < num_obs; j++)
		{

			pointer_k_xy = &matrix_weights_K_xy[j][0];
			pointer_k_x = &matrix_weights_K_x[j][0];
			pointer_k_convol_y = &matrix_weights_K_convol_y[j][0];

			for(i=0; i < num_obs; i++)
			{

				prod_kernel_cont = 1.0;
				prod_kernel_convol_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][j])/matrix_bandwidth_reg[l][j];
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][j])/matrix_bandwidth_var[l][j];
					prod_kernel_convol_cont *= kernel_convol(KERNEL_den,BANDWIDTH_den,
						(matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][j],matrix_bandwidth_var[l][i],matrix_bandwidth_var[l][j])/matrix_bandwidth_var[l][j];
				}

				prod_kernel_cat = 1.0;
				prod_kernel_convol_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg,matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
					prod_kernel_convol_cat *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l], num_categories[l], matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
					prod_kernel_convol_cat *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered], num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
				}

				*pointer_k_xy++ = prod_kernel_cont*prod_kernel_cat;
				*pointer_k_x++ = prod_kernel_marginal_cont*prod_kernel_marginal_cat;
				*pointer_k_convol_y++ = prod_kernel_convol_cont*prod_kernel_convol_cat;

			}

		}

	}
	else
	{

		for(j=0; j < num_obs; j++)
		{

			pointer_k_xy = &matrix_weights_K_xy[j][0];
			pointer_k_x = &matrix_weights_K_x[j][0];
			pointer_k_convol_y = &matrix_weights_K_convol_y[j][0];

			for(i=0; i < num_obs; i++)
			{

				prod_kernel_cont = 1.0;
				prod_kernel_convol_cont = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_reg, (matrix_X_continuous[l][j]-matrix_X_continuous[l][i])/matrix_bandwidth_reg[l][i])/matrix_bandwidth_reg[l][i];
				}

				prod_kernel_marginal_cont = prod_kernel_cont;

				for(l = 0; l < num_var_continuous; l++)
				{
					prod_kernel_cont *= kernel(KERNEL_den, (matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][i])/matrix_bandwidth_var[l][i];
					prod_kernel_convol_cont *= kernel_convol(KERNEL_den,BANDWIDTH_den,
						(matrix_Y_continuous[l][j]-matrix_Y_continuous[l][i])/matrix_bandwidth_var[l][i],matrix_bandwidth_var[l][j],matrix_bandwidth_var[l][i])/matrix_bandwidth_var[l][i];
				}

				prod_kernel_cat = 1.0;
				prod_kernel_convol_cat = 1.0;

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_reg, matrix_X_unordered[l][j],matrix_X_unordered[l][i],lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_reg,matrix_X_ordered[l][j],matrix_X_ordered[l][i],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
				}

				prod_kernel_marginal_cat = prod_kernel_cat;

				for(l = 0; l < num_var_unordered; l++)
				{
					prod_kernel_cat *= kernel_unordered(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l],num_categories[l]);
					prod_kernel_convol_cat *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_Y_unordered[l][j],matrix_Y_unordered[l][i],lambda[l], num_categories[l], matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_var_ordered; l++)
				{
					prod_kernel_cat *= kernel_ordered(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered]);
					prod_kernel_convol_cat *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_Y_ordered[l][j],matrix_Y_ordered[l][i],lambda[l+num_var_unordered], num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
				}

				*pointer_k_xy++ = prod_kernel_cont*prod_kernel_cat;
				*pointer_k_x++ = prod_kernel_marginal_cont*prod_kernel_marginal_cat;
				*pointer_k_convol_y++ = prod_kernel_convol_cont*prod_kernel_convol_cat;

			}

		}

	}

	return(0);

}


