/* Copyright (C) Jeff Racine, 1995-2004 */

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include <R.h>

#ifdef MPI2
#include "mpi.h"
int my_rank;
int iNum_Processors;
int source;
int dest;
int tag;
int iSeed_my_rank;
MPI_Status status;
extern MPI_Comm	*comm;
#endif

/*
-f to enable Fast resampling (memory intensive)
-n to use Nonparametric measures of central tendency and dispersion
-o to use Ordered categorical gradients
-p to compute univariate categorical conditional Predictions
-r to enable a different Random seed with each program invocation

*/

/* headers.h has all definitions of routines used by main() and related modules */

#include "headers.h"
#include "matrix.h"

#ifdef RCSID
static char rcsid[] = "$Id: np.c,v 1.35 2006/11/02 16:56:49 tristen Exp $";
#endif

int int_DEBUG;
int int_VERBOSE;
int int_TAYLOR;
int int_WEIGHTS = 0;
int int_LARGE_SF;
int int_NOKEYPRESS;
int int_DISPLAY_CV;
int int_RANDOM_SEED = 42;
int int_MINIMIZE_IO;
int int_ORDERED_CATEGORICAL_GRADIENT;
int int_PREDICT;
int int_ROBUST;
int int_SIMULATION;

int int_RESTART_FROM_MIN;

/* Some externals for numerical routines */
/* Some externals for numerical routines */

int num_obs_train_extern=0;
int num_obs_eval_extern=0;
int num_var_continuous_extern=0;
int num_var_unordered_extern=0;
int num_var_ordered_extern=0;
int num_reg_continuous_extern=0;
int num_reg_unordered_extern=0;
int num_reg_ordered_extern=0;


int *num_categories_extern;
double **matrix_categorical_vals_extern;

double **matrix_X_continuous_train_extern;
double **matrix_X_unordered_train_extern;
double **matrix_X_ordered_train_extern;
double **matrix_X_continuous_eval_extern;
double **matrix_X_unordered_eval_extern;
double **matrix_X_ordered_eval_extern;

double **matrix_Y_continuous_train_extern;
double **matrix_Y_unordered_train_extern;
double **matrix_Y_ordered_train_extern;

double **matrix_Y_continuous_eval_extern;
double **matrix_Y_unordered_eval_extern;
double **matrix_Y_ordered_eval_extern;

/* Quantile - no Y ordered or unordered used, but defined anyways */

double **matrix_Y_continuous_quantile_extern;
double **matrix_Y_unordered_quantile_extern;
double **matrix_Y_ordered_quantile_extern;
double **matrix_X_continuous_quantile_extern;
double **matrix_X_unordered_quantile_extern;
double **matrix_X_ordered_quantile_extern;

double *vector_Y_extern;
double *vector_T_extern;
double *vector_T_resample;
double *vector_Y_eval_extern;
double *vector_Y_null;


int int_ll_extern=0;

int KERNEL_reg_extern=0;
int KERNEL_reg_unordered_extern=0;
int KERNEL_reg_ordered_extern=0;
int KERNEL_den_extern=0;
int KERNEL_den_unordered_extern=0;
int KERNEL_den_ordered_extern=0;
int BANDWIDTH_reg_extern;
int BANDWIDTH_den_extern;

/* Statics for dependence metric */

int num_lag_extern;
int int_lag_extern;
int int_iter_extern;

int itmax_extern;
double small_extern;

double *vector_scale_factor_dep_met_bivar_extern;
double *vector_scale_factor_dep_met_univar_extern;
double *vector_scale_factor_dep_met_univar_lag_extern;

double *vector_scale_factor_extern;
double gamma_extern = 0.5;
double y_min_extern;
double y_max_extern;

int imsnum = 0;
int imstot = 0;

extern int iff;

void spinner(int num) {
  if(int_MINIMIZE_IO == IO_MIN_FALSE){
    const char spinney[] = { '|', '/', '-', '\\' };
    Rprintf("\rMultistart %d of %d %c", imsnum+1, imstot, spinney[num%4]);
    R_FlushConsole();
  }
}

void np_set_seed(int * num){
  int_RANDOM_SEED = *num;
  iff = 0;
}

void np_mpi_init(int * mpi_status){
#ifdef MPI2 
  MPI_Comm_rank(comm[1], &my_rank);
  MPI_Comm_size(comm[1], &iNum_Processors);
  mpi_status[MPI_RANKI] = my_rank;
  mpi_status[MPI_NUMPI] = iNum_Processors;
#else
  mpi_status[MPI_RANKI] = -1;
  mpi_status[MPI_NUMPI] = -1;
#endif
}


void np_density_bw(double * myuno, double * myord, double * mycon, 
                   int * myopti, double * myoptd, double * myans, double * fval){
  /* Likelihood bandwidth selection for density estimation */

  double **matrix_y;

  double *vector_continuous_stddev;
  double *vector_scale_factor, *vector_scale_factor_multistart;

  double fret, fret_best;
  double ftol, tol;
  double (* bwmfunc)(double *);

  double small;
  
  int i,j;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, iImproved;
  int itmax, iter;
  int int_use_starting_values;

  num_reg_unordered_extern = myopti[BW_NUNOI];
  num_reg_ordered_extern = myopti[BW_NORDI];
  num_reg_continuous_extern = myopti[BW_NCONI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[BW_NOBSI];
  iMultistart = myopti[BW_IMULTII];
  iNum_Multistart = myopti[BW_NMULTII];

  KERNEL_den_extern = myopti[BW_CKRNEVI];
  KERNEL_den_unordered_extern = 0;
  KERNEL_den_ordered_extern = 0;

  int_use_starting_values= myopti[BW_USTARTI];
  int_LARGE_SF=myopti[BW_LSFI];
  BANDWIDTH_den_extern=myopti[BW_DENI];
  int_RESTART_FROM_MIN = myopti[BW_REMINI];
  int_MINIMIZE_IO = myopti[BW_MINIOI];

  itmax=myopti[BW_ITMAXI];

  ftol=myoptd[BW_FTOLD];
  tol=myoptd[BW_TOLD];
  small=myoptd[BW_SMALLD];

/* Allocate memory for objects */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);


  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_y = alloc_matd(num_var + 1, num_var +1);
  vector_scale_factor = alloc_vecd(num_var + 1);
  matrix_categorical_vals_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  if (int_use_starting_values)
    for( i=0;i<num_var; i++ )
      vector_scale_factor[i+1] = myans[i];

/* Parse data */

  for( j=0;j<num_reg_unordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_unordered_train_extern[j][i]=myuno[j*num_obs_train_extern+i];
    

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=myord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=mycon[j*num_obs_train_extern+i];


  determine_categorical_vals(
                             num_obs_train_extern,
                             0,
                             0,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             matrix_Y_unordered_train_extern,
                             matrix_Y_ordered_train_extern,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern,
                             matrix_categorical_vals_extern);

  vector_continuous_stddev = alloc_vecd(num_reg_continuous_extern);

  compute_continuous_stddev(
                            int_LARGE_SF,
                            num_obs_train_extern,
                            0,
                            num_reg_continuous_extern,
                            matrix_Y_continuous_train_extern,
                            matrix_X_continuous_train_extern,
                            vector_continuous_stddev);

  /* Initialize scale factors and Hessian for NR modules */

  initialize_nr_vector_scale_factor(
                                    BANDWIDTH_den_extern,
                                    BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    0,                /* regression (0) regression ml (1) */
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    matrix_Y_continuous_train_extern,
                                    matrix_X_continuous_train_extern,
                                    int_use_starting_values,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor);

  initialize_nr_hessian(num_var, matrix_y);

  /* When multistarting, set counter */

  imsnum = iMs_counter = 0;
  imstot = iNum_Multistart;
  
  /* Conduct direction set search */

  /* assign the function to be optimized */
  switch(myopti[BW_MI]){
  case BWM_CVML : bwmfunc = cv_func_density_categorical_ml; break;
  case BWM_CVLS : bwmfunc = cv_func_density_categorical_ls; break;
    //case BWM_CVML_NP : bwmfunc = cv_func_np_density_categorical_ml; break;
  default : REprintf("np.c: invalid bandwidth selection method.");
    exit(0); break;
  }

  spinner(0);

  fret_best = bwmfunc(vector_scale_factor);
  iImproved = 0;

  powell(0,
         0,
         vector_scale_factor,
         vector_scale_factor,
         matrix_y,
         num_var,
         ftol,
         tol,
         small,
         itmax,
         &iter,
         &fret,
         bwmfunc);

  /* int_RESTART_FROM_MIN needs to be set */

  if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

    initialize_nr_hessian(num_var, matrix_y);

    powell(0,
           0,
           vector_scale_factor,
           vector_scale_factor,
           matrix_y,
           num_var,
           ftol,
           tol,
           small,
           itmax,
           &iter,
           &fret,
           bwmfunc);
  }

  iImproved = (fret < fret_best);

  /* When multistarting save initial minimum of objective function and scale factors */

  if(iMultistart == IMULTI_TRUE){
    fret_best = fret;
    vector_scale_factor_multistart = alloc_vecd(num_var + 1);
    for(i = 1; i <= num_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
    		
    /* Conduct search from new random values of the search parameters */
       	
    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and hessian for NR modules */

      initialize_nr_vector_scale_factor(
                                        BANDWIDTH_den_extern,
                                        BANDWIDTH_den_extern,
                                        1,        /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        0,        /* regression (0) regression ml (1) */
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        0,
                                        0,
                                        0,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        matrix_Y_continuous_train_extern,
                                        matrix_X_continuous_train_extern,
                                        int_use_starting_values,
                                        pow((double)4.0/(double)3.0,0.2),     /* Init for continuous vars */
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor);

      initialize_nr_hessian(num_var, matrix_y);

      /* Conduct direction set search */

      powell(0,
             0,
             vector_scale_factor,
             vector_scale_factor,
             matrix_y,
             num_var,
             ftol,
             tol,
             small,
             itmax,
             &iter,
             &fret,
             bwmfunc);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

        initialize_nr_hessian(num_var, matrix_y);

        powell(0,
               0,
               vector_scale_factor,
               vector_scale_factor,
               matrix_y,
               num_var,
               ftol,
               tol,
               small,
               itmax,
               &iter,
               &fret,
               bwmfunc);
      }
       			
      /* If this run resulted in an improved minimum save information */
       		
      if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter;
       
        for(i = 1; i <= num_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }

    }

    /* Save best for estimation */

    fret = fret_best;

    for(i = 1; i <= num_var; i++)
      vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];

    free(vector_scale_factor_multistart);

  }

  /* return data to R */
  if (BANDWIDTH_den_extern == BW_GEN_NN || 
      BANDWIDTH_den_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern; i++ )
      vector_scale_factor[i+1]=fround(vector_scale_factor[i+1]);
  }

  for( i=0; i<num_var; i++ )
    myans[i]=vector_scale_factor[i+1];

  fval[0] = fret;
  fval[1] = iImproved;
  /* end return data */

  /* Free data objects */

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  free_mat(matrix_y, num_var + 1);
  free(vector_scale_factor);
  free(num_categories_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  free(vector_continuous_stddev);

  if(int_MINIMIZE_IO != IO_MIN_TRUE)
    Rprintf("\r                   \r");

  return ;
  
}


void np_density_conditional_bw(double * c_uno, double * c_ord, double * c_con, 
                               double * u_uno, double * u_ord, double * u_con,
                               int * myopti, double * myoptd, double * myans, double * fval){
/* Likelihood bandwidth selection for density estimation */

  double **matrix_y;

  double *vector_continuous_stddev;
  double *vector_scale_factor, *vector_scale_factor_multistart;

  double fret, fret_best;
  double ftol, tol;
  double (* bwmfunc)(double *);

  double small;
  
  int i,j;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, num_all_var, num_var_var, iImproved;
  int itmax, iter;
  int int_use_starting_values, autoSelectCVLS, ibwmfunc;

#ifdef MPI2
  int work_np;
#endif // MPI

  num_var_unordered_extern = myopti[CBW_CNUNOI];
  num_var_ordered_extern = myopti[CBW_CNORDI];
  num_var_continuous_extern = myopti[CBW_CNCONI];

  num_reg_unordered_extern = myopti[CBW_UNUNOI];
  num_reg_ordered_extern = myopti[CBW_UNORDI];
  num_reg_continuous_extern = myopti[CBW_UNCONI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;
  num_var_var = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  num_all_var = num_var+num_var_var;

  num_obs_train_extern = myopti[CBW_NOBSI];
  iMultistart = myopti[CBW_IMULTII];
  iNum_Multistart = myopti[CBW_NMULTII];

  KERNEL_reg_extern = myopti[CBW_CXKRNEVI];
  KERNEL_den_extern = myopti[CBW_CYKRNEVI];

  KERNEL_reg_unordered_extern = myopti[CBW_UXKRNEVI];
  KERNEL_den_unordered_extern = myopti[CBW_UYKRNEVI];

  KERNEL_reg_ordered_extern = myopti[CBW_OXKRNEVI];
  KERNEL_den_ordered_extern = myopti[CBW_OYKRNEVI];

  int_use_starting_values= myopti[CBW_USTARTI];
  int_LARGE_SF=myopti[CBW_LSFI];
  BANDWIDTH_den_extern=myopti[CBW_DENI];
  int_RESTART_FROM_MIN = myopti[CBW_REMINI];
  int_MINIMIZE_IO = myopti[CBW_MINIOI];

  itmax=myopti[CBW_ITMAXI];
  int_WEIGHTS = myopti[CBW_FASTI];
  autoSelectCVLS = myopti[CBW_AUTOI];

  ftol=myoptd[CBW_FTOLD];
  tol=myoptd[CBW_TOLD];
  small=myoptd[CBW_SMALLD];

/* Allocate memory for objects */

  matrix_Y_unordered_train_extern = alloc_matd(num_obs_train_extern, num_var_unordered_extern);
  matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);
  matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);
	
  num_categories_extern = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                     num_reg_unordered_extern + num_reg_ordered_extern);
  matrix_y = alloc_matd(num_all_var + 1, num_all_var + 1);
  vector_scale_factor = alloc_vecd(num_all_var + 1);
  
  matrix_categorical_vals_extern = 
    alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern + 
               num_reg_unordered_extern + num_reg_ordered_extern);

  /* in v_s_f order is creg, cvar, uvar, ovar, ureg, oreg  */

  if (int_use_starting_values)
    for( i=0;i<num_all_var; i++ )
      vector_scale_factor[i+1] = myans[i];

/* Parse data */

  for(j=0;j<num_var_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_unordered_train_extern[j][i]=c_uno[j*num_obs_train_extern+i];

  for(j=0;j<num_var_ordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_ordered_train_extern[j][i]=c_ord[j*num_obs_train_extern+i];

  for(j=0;j<num_var_continuous_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_continuous_train_extern[j][i]=c_con[j*num_obs_train_extern+i];


  for(j=0;j<num_reg_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_X_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+i];


  determine_categorical_vals(
                             num_obs_train_extern,
                             num_var_unordered_extern,
                             num_var_ordered_extern,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             matrix_Y_unordered_train_extern,
                             matrix_Y_ordered_train_extern,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern,
                             matrix_categorical_vals_extern);


  vector_continuous_stddev = alloc_vecd(num_var_continuous_extern + num_reg_continuous_extern);

  compute_continuous_stddev(
                            int_LARGE_SF,
                            num_obs_train_extern,
                            num_var_continuous_extern,
                            num_reg_continuous_extern,
                            matrix_Y_continuous_train_extern,
                            matrix_X_continuous_train_extern,
                            vector_continuous_stddev);

  /* Initialize scale factors and Hessian for NR modules */

  initialize_nr_vector_scale_factor(
                                    BANDWIDTH_den_extern,
                                    BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    0,                /* regression (0) regression ml (1) */
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    num_var_continuous_extern,
                                    num_var_unordered_extern,
                                    num_var_ordered_extern,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    matrix_Y_continuous_train_extern,
                                    matrix_X_continuous_train_extern,
                                    int_use_starting_values,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor);

  initialize_nr_hessian(num_all_var, matrix_y);

  /* When multistarting, set counter */

  imsnum = iMs_counter = 0;
  imstot = iNum_Multistart;

  /* Conduct direction set search */

  /* assign the function to be optimized */

  ibwmfunc = myopti[CBW_MI];

  /* 7/2/2010 */
  
  /*  if((ibwmfunc != CBWM_CVML) && autoSelectCVLS){*/
  if((ibwmfunc != CBWM_CVML && ibwmfunc != CBWM_CCDF) && autoSelectCVLS){
#ifdef MPI2
    int nobs_proc = (num_obs_train_extern / iNum_Processors) + 
      ((num_obs_train_extern%iNum_Processors) != 0);
    
    ibwmfunc = (nobs_proc < CBW_MINOBS) ? CBWM_CVLS : CBWM_NPLS;
    int_WEIGHTS = (nobs_proc < CBW_MINOBS);
#else //MPI
    ibwmfunc = (num_obs_train_extern < CBW_MINOBS) ? CBWM_CVLS : CBWM_NPLS;
    int_WEIGHTS = (num_obs_train_extern < CBW_MINOBS);
#endif //MPI
  }

  switch(ibwmfunc){
  case CBWM_CVML : bwmfunc = cv_func_con_density_categorical_ml; break;
  case CBWM_CVLS : bwmfunc = cv_func_con_density_categorical_ls; break;
  case CBWM_NPLS : bwmfunc = np_cv_func_con_density_categorical_ls;break;
  case CBWM_CCDF : bwmfunc = cv_func_con_distribution_categorical_ccdf; break;
  default : REprintf("np.c: invalid bandwidth selection method.");
    exit(0); break;
  }

  spinner(0);

  fret_best = bwmfunc(vector_scale_factor);
  iImproved = 0;

  powell(0,
         0,
         vector_scale_factor,
         vector_scale_factor,
         matrix_y,
         num_all_var,
         ftol,
         tol,
         small,
         itmax,
         &iter,
         &fret,
         bwmfunc);

  if(int_RESTART_FROM_MIN == RE_MIN_TRUE){
    initialize_nr_hessian(num_all_var, matrix_y);

    powell(0,
           0,
           vector_scale_factor,
           vector_scale_factor,
           matrix_y,
           num_all_var,
           ftol,
           tol,
           small,
           itmax,
           &iter,
           &fret,
           bwmfunc);

  }

  iImproved = (fret < fret_best);

  /* When multistarting save initial minimum of objective function and scale factors */


  if(iMultistart == IMULTI_TRUE){
    fret_best = fret;
    vector_scale_factor_multistart = alloc_vecd(num_all_var + 1);
    for(i = 1; i <= num_all_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
			

    /* Conduct search from new random values of the search parameters */
		
    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and hessian for NR modules */
      initialize_nr_vector_scale_factor(
                                        BANDWIDTH_den_extern,
                                        BANDWIDTH_den_extern,
                                        1,                /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        0,                /* regression (0) regression ml (1) */
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        num_var_continuous_extern,
                                        num_var_unordered_extern,
                                        num_var_ordered_extern,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        matrix_Y_continuous_train_extern,
                                        matrix_X_continuous_train_extern,
                                        int_use_starting_values,
                                        pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor);

      initialize_nr_hessian(num_all_var, matrix_y);

      /* Conduct direction set search */
      
      powell(0,
             0,
             vector_scale_factor,
             vector_scale_factor,
             matrix_y,
             num_all_var,
             ftol,
             tol,
             small,
             itmax,
             &iter,
             &fret,
             bwmfunc);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

        initialize_nr_hessian(num_all_var, matrix_y);

        powell(0,
               0,
               vector_scale_factor,
               vector_scale_factor,
               matrix_y,
               num_all_var,
               ftol,
               tol,
               small,
               itmax,
               &iter,
               &fret,
               bwmfunc);
      }
				
      /* If this run resulted in an improved minimum save information */
      
      if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        
        for(i = 1; i <= num_all_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }
      
    }

    /* Save best for estimation */

    fret = fret_best;
    for(i = 1; i <= num_all_var; i++)
      vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];
    free(vector_scale_factor_multistart);
  }

  /* return data to R */
  if (BANDWIDTH_den_extern == BW_GEN_NN || 
      BANDWIDTH_den_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern+num_var_continuous_extern; i++ )
      vector_scale_factor[i+1]=fround(vector_scale_factor[i+1]);
  }

  for( i=0; i<num_all_var; i++ )
    myans[i]=vector_scale_factor[i+1];

  fval[0] = fret;
  fval[1] = iImproved;
  /* end return data */

  /* Free data objects */

  free_mat(matrix_Y_unordered_train_extern, num_var_unordered_extern);
  free_mat(matrix_Y_ordered_train_extern, num_var_ordered_extern);
  free_mat(matrix_Y_continuous_train_extern, num_var_continuous_extern);

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  free_mat(matrix_y, num_all_var + 1);
  safe_free(vector_scale_factor);
  safe_free(num_categories_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);

  safe_free(vector_continuous_stddev);

  int_WEIGHTS = 0;

  if(int_MINIMIZE_IO != IO_MIN_TRUE)
    Rprintf("\r                   \r");

  return ;
}


void np_density_conditional(double * tc_uno, double * tc_ord, double * tc_con, 
                            double * tu_uno, double * tu_ord, double * tu_con,
                            double * ec_uno, double * ec_ord, double * ec_con, 
                            double * eu_uno, double * eu_ord, double * eu_con,
                            double * mybw, 
                            double * ymcv, double * ypadnum,
                            double * xmcv, double * xpadnum,
                            int * myopti, 
                            double * cdens, double * cderr, 
                            double * cg, double * cgerr,
                            double * ll){
  /* Likelihood bandwidth selection for density estimation */

  double small = 1.0e-16;
  double *vector_scale_factor, *pdf, *pdf_stderr, log_likelihood = 0.0;
  double ** pdf_deriv = NULL, ** pdf_deriv_stderr = NULL;
  double xpad_num, ypad_num;

  int itmax = 10000;
  int i,j;
  int num_var;

  int num_all_var, num_var_var, train_is_eval, do_grad, num_obs_eval_alloc;
  int xmax_lev, ymax_lev, dens_or_dist, t_num;


  num_var_unordered_extern = myopti[CD_CNUNOI];
  num_var_ordered_extern = myopti[CD_CNORDI];
  num_var_continuous_extern = myopti[CD_CNCONI];

  num_reg_unordered_extern = myopti[CD_UNUNOI];
  num_reg_ordered_extern = myopti[CD_UNORDI];
  num_reg_continuous_extern = myopti[CD_UNCONI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;
  num_var_var = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  num_all_var = num_var + num_var_var;

  num_obs_train_extern = myopti[CD_TNOBSI];
  num_obs_eval_extern = myopti[CD_ENOBSI];

  if((train_is_eval = myopti[CD_TISEI]) && 
     (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_density_conditional): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    exit(0);
  }

  KERNEL_reg_extern = myopti[CD_CXKRNEVI];
  KERNEL_den_extern = myopti[CD_CYKRNEVI];

  KERNEL_reg_unordered_extern = myopti[CD_UXKRNEVI];
  KERNEL_den_unordered_extern = myopti[CD_UYKRNEVI];

  KERNEL_reg_ordered_extern = myopti[CD_OXKRNEVI];
  KERNEL_den_ordered_extern = myopti[CD_OYKRNEVI];

  int_LARGE_SF = myopti[CD_LSFI];
  BANDWIDTH_den_extern = myopti[CD_DENI];
  int_MINIMIZE_IO = myopti[CD_MINIOI];
  do_grad = myopti[CD_GRAD];

  ymax_lev = myopti[CD_YMLEVI];
  xmax_lev = myopti[CD_XMLEVI];

  ypad_num = *ypadnum;
  xpad_num = *xpadnum;

  dens_or_dist = myopti[CD_DODENI];


#ifdef MPI2
  num_obs_eval_alloc = MAX(ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif

  /* Allocate memory for objects */

  matrix_Y_unordered_train_extern = alloc_matd(num_obs_train_extern, num_var_unordered_extern);
  matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);
  matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  if(train_is_eval) {
    matrix_Y_unordered_eval_extern = matrix_Y_unordered_train_extern;
    matrix_Y_ordered_eval_extern = matrix_Y_ordered_train_extern;
    matrix_Y_continuous_eval_extern = matrix_Y_continuous_train_extern;

    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;
  } else {
    matrix_Y_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_var_unordered_extern);
    matrix_Y_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_var_ordered_extern);
    matrix_Y_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_var_continuous_extern);

    matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
    matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
    matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);
  }
	
  num_categories_extern = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                     num_reg_unordered_extern + num_reg_ordered_extern);
  vector_scale_factor = alloc_vecd(num_all_var + 1);
  
  matrix_categorical_vals_extern = 
    alloc_matd(MAX(xmax_lev,ymax_lev), num_var_unordered_extern + num_var_ordered_extern + 
               num_reg_unordered_extern + num_reg_ordered_extern);

  /* notice use of num_obs_eval_alloc for MPI compatibility */
  pdf = alloc_vecd(num_obs_eval_alloc);
  pdf_stderr = alloc_vecd(num_obs_eval_alloc);

  if (do_grad){
    pdf_deriv = alloc_matd(num_obs_eval_alloc, num_var);
    pdf_deriv_stderr = alloc_matd(num_obs_eval_alloc, num_var);
  }

  /* in v_s_f order is creg, cvar, uvar, ovar, ureg, oreg  */

  for( i=0;i<num_all_var; i++ )
    vector_scale_factor[i+1] = mybw[i];

  /* Parse data */

  /* train */

  for(j=0;j<num_var_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_unordered_train_extern[j][i]=tc_uno[j*num_obs_train_extern+i];

  for(j=0;j<num_var_ordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_ordered_train_extern[j][i]=tc_ord[j*num_obs_train_extern+i];

  for(j=0;j<num_var_continuous_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_continuous_train_extern[j][i]=tc_con[j*num_obs_train_extern+i];


  for(j=0;j<num_reg_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_X_unordered_train_extern[j][i]=tu_uno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=tu_ord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=tu_con[j*num_obs_train_extern+i];

  /* eval */
  if(!train_is_eval){
    for(j=0;j<num_var_unordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_Y_unordered_eval_extern[j][i]=ec_uno[j*num_obs_eval_extern+i];

    for(j=0;j<num_var_ordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_Y_ordered_eval_extern[j][i]=ec_ord[j*num_obs_eval_extern+i];

    for(j=0;j<num_var_continuous_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_Y_continuous_eval_extern[j][i]=ec_con[j*num_obs_eval_extern+i];


    for(j=0;j<num_reg_unordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_X_unordered_eval_extern[j][i]=eu_uno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=eu_ord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=eu_con[j*num_obs_eval_extern+i];
  }

  /* fix up categories */
  for(j=0; j < (num_var_unordered_extern + num_var_ordered_extern); j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = ymcv[j*ymax_lev+i];
    } while(++i < ymax_lev && ymcv[j*ymax_lev+i] != ypad_num);
    num_categories_extern[j] = i;
  }

  t_num = j;

  for(j=0; j < (num_reg_unordered_extern+num_reg_ordered_extern); j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j+t_num][i] = xmcv[j*xmax_lev+i];
    } while(++i < xmax_lev && xmcv[j*xmax_lev+i] != xpad_num);
    num_categories_extern[j+t_num] = i;
  }

  if (dens_or_dist == NP_DO_DENS) {
    if (!do_grad) {
      kernel_estimate_con_density_categorical(KERNEL_den_extern,
                                              KERNEL_den_unordered_extern,
                                              KERNEL_den_ordered_extern,
                                              KERNEL_reg_extern,
                                              KERNEL_reg_unordered_extern,
                                              KERNEL_reg_ordered_extern,
                                              BANDWIDTH_den_extern,
                                              num_obs_train_extern,
                                              num_obs_eval_extern,
                                              num_var_unordered_extern,
                                              num_var_ordered_extern,
                                              num_var_continuous_extern,
                                              num_reg_unordered_extern,
                                              num_reg_ordered_extern,
                                              num_reg_continuous_extern,
                                              /* Train */
                                              matrix_Y_unordered_train_extern,
                                              matrix_Y_ordered_train_extern,
                                              matrix_Y_continuous_train_extern,
                                              /* Eval */
                                              matrix_Y_unordered_eval_extern,
                                              matrix_Y_ordered_eval_extern,
                                              matrix_Y_continuous_eval_extern,
                                              /* Train */
                                              matrix_X_unordered_train_extern,
                                              matrix_X_ordered_train_extern,
                                              matrix_X_continuous_train_extern,
                                              /* Eval */
                                              matrix_X_unordered_eval_extern,
                                              matrix_X_ordered_eval_extern,
                                              matrix_X_continuous_eval_extern,
                                              &vector_scale_factor[1],
                                              num_categories_extern,
                                              pdf,
                                              pdf_stderr,
                                              &log_likelihood);
    } else {
      kernel_estimate_con_density_categorical_gradient(KERNEL_den_extern,
                                                       KERNEL_den_unordered_extern,
                                                       KERNEL_den_ordered_extern,
                                                       KERNEL_reg_extern,
                                                       KERNEL_reg_unordered_extern,
                                                       KERNEL_reg_ordered_extern,
                                                       BANDWIDTH_den_extern,
                                                       num_obs_train_extern,
                                                       num_obs_eval_extern,
                                                       num_var_unordered_extern,
                                                       num_var_ordered_extern,
                                                       num_var_continuous_extern,
                                                       num_reg_unordered_extern,
                                                       num_reg_ordered_extern,
                                                       num_reg_continuous_extern,
                                                       /* Train */
                                                       matrix_Y_unordered_train_extern,
                                                       matrix_Y_ordered_train_extern,
                                                       matrix_Y_continuous_train_extern,
                                                       /* Eval */
                                                       matrix_Y_unordered_eval_extern,
                                                       matrix_Y_ordered_eval_extern,
                                                       matrix_Y_continuous_eval_extern,
                                                       /* Train */
                                                       matrix_X_unordered_train_extern,
                                                       matrix_X_ordered_train_extern,
                                                       matrix_X_continuous_train_extern,
                                                       /* Eval */
                                                       matrix_X_unordered_eval_extern,
                                                       matrix_X_ordered_eval_extern,
                                                       matrix_X_continuous_eval_extern,
                                                       &vector_scale_factor[1],
                                                       num_categories_extern,
                                                       pdf,
                                                       pdf_stderr,
                                                       pdf_deriv,
                                                       pdf_deriv_stderr,
                                                       &log_likelihood);

      kernel_estimate_con_density_categorical_gradient_categorical(KERNEL_den_extern,
                                                                   KERNEL_den_unordered_extern,
                                                                   KERNEL_den_ordered_extern,
                                                                   KERNEL_reg_extern,
                                                                   KERNEL_reg_unordered_extern,
                                                                   KERNEL_reg_ordered_extern,
                                                                   BANDWIDTH_den_extern,
                                                                   num_obs_train_extern,
                                                                   num_obs_eval_extern,
                                                                   num_var_unordered_extern,
                                                                   num_var_ordered_extern,
                                                                   num_var_continuous_extern,
                                                                   num_reg_unordered_extern,
                                                                   num_reg_ordered_extern,
                                                                   num_reg_continuous_extern,
                                                                   0,
                                                                   /* Train */
                                                                   matrix_Y_unordered_train_extern,
                                                                   matrix_Y_ordered_train_extern,
                                                                   matrix_Y_continuous_train_extern,
                                                                   /* Eval */
                                                                   matrix_Y_unordered_eval_extern,
                                                                   matrix_Y_ordered_eval_extern,
                                                                   matrix_Y_continuous_eval_extern,
                                                                   /* Train */
                                                                   matrix_X_unordered_train_extern,
                                                                   matrix_X_ordered_train_extern,
                                                                   matrix_X_continuous_train_extern,
                                                                   /* Eval */
                                                                   matrix_X_unordered_eval_extern,
                                                                   matrix_X_ordered_eval_extern,
                                                                   matrix_X_continuous_eval_extern,
                                                                   &vector_scale_factor[1],
                                                                   matrix_categorical_vals_extern,
                                                                   num_categories_extern,
                                                                   pdf,
                                                                   &pdf_deriv[num_reg_continuous_extern],
                                                                   &pdf_deriv_stderr[num_reg_continuous_extern]);

    }
  } else if (dens_or_dist == NP_DO_DIST) {
    if(!do_grad){
      kernel_estimate_con_distribution_categorical(KERNEL_den_extern,
                                                   KERNEL_den_unordered_extern,
                                                   KERNEL_den_ordered_extern,
                                                   KERNEL_reg_extern,
                                                   KERNEL_reg_unordered_extern,
                                                   KERNEL_reg_ordered_extern,
                                                   BANDWIDTH_den_extern,
                                                   num_obs_train_extern,
                                                   num_obs_eval_extern,
                                                   num_var_unordered_extern,
                                                   num_var_ordered_extern,
                                                   num_var_continuous_extern,
                                                   num_reg_unordered_extern,
                                                   num_reg_ordered_extern,
                                                   num_reg_continuous_extern,
                                                   /* Train */
                                                   matrix_Y_unordered_train_extern,
                                                   matrix_Y_ordered_train_extern,
                                                   matrix_Y_continuous_train_extern,
                                                   /* Eval */
                                                   matrix_Y_unordered_eval_extern,
                                                   matrix_Y_ordered_eval_extern,
                                                   matrix_Y_continuous_eval_extern,
                                                   /* Train */
                                                   matrix_X_unordered_train_extern,
                                                   matrix_X_ordered_train_extern,
                                                   matrix_X_continuous_train_extern,
                                                   /* Eval */
                                                   matrix_X_unordered_eval_extern,
                                                   matrix_X_ordered_eval_extern,
                                                   matrix_X_continuous_eval_extern,
                                                   &vector_scale_factor[1],
                                                   num_categories_extern,
                                                   matrix_categorical_vals_extern,
                                                   pdf,
                                                   pdf_stderr,
                                                   small, itmax);
    } else {
      kernel_estimate_con_distribution_categorical_gradient(KERNEL_den_extern,
                                                            KERNEL_den_unordered_extern,
                                                            KERNEL_den_ordered_extern,
                                                            KERNEL_reg_extern,
                                                            KERNEL_reg_unordered_extern,
                                                            KERNEL_reg_ordered_extern,
                                                            BANDWIDTH_den_extern,
                                                            num_obs_train_extern,
                                                            num_obs_eval_extern,
                                                            num_var_unordered_extern,
                                                            num_var_ordered_extern,
                                                            num_var_continuous_extern,
                                                            num_reg_unordered_extern,
                                                            num_reg_ordered_extern,
                                                            num_reg_continuous_extern,
                                                            /* Train */
                                                            matrix_Y_unordered_train_extern,
                                                            matrix_Y_ordered_train_extern,
                                                            matrix_Y_continuous_train_extern,
                                                            /* Eval */
                                                            matrix_Y_unordered_eval_extern,
                                                            matrix_Y_ordered_eval_extern,
                                                            matrix_Y_continuous_eval_extern,
                                                            /* Train */
                                                            matrix_X_unordered_train_extern,
                                                            matrix_X_ordered_train_extern,
                                                            matrix_X_continuous_train_extern,
                                                            /* Eval */
                                                            matrix_X_unordered_eval_extern,
                                                            matrix_X_ordered_eval_extern,
                                                            matrix_X_continuous_eval_extern,
                                                            &vector_scale_factor[1],
                                                            num_categories_extern,
                                                            matrix_categorical_vals_extern,
                                                            pdf,
                                                            pdf_stderr,
                                                            pdf_deriv,
                                                            pdf_deriv_stderr,
                                                            small, itmax);

      kernel_estimate_con_distribution_categorical_gradient_categorical(KERNEL_den_extern,
                                                                        KERNEL_den_unordered_extern,
                                                                        KERNEL_den_ordered_extern,
                                                                        KERNEL_reg_extern,
                                                                        KERNEL_reg_unordered_extern,
                                                                        KERNEL_reg_ordered_extern,
                                                                        BANDWIDTH_den_extern,
                                                                        num_obs_train_extern,
                                                                        num_obs_eval_extern,
                                                                        num_var_unordered_extern,
                                                                        num_var_ordered_extern,
                                                                        num_var_continuous_extern,
                                                                        num_reg_unordered_extern,
                                                                        num_reg_ordered_extern,
                                                                        num_reg_continuous_extern,
                                                                        0,
                                                                        /* Train */
                                                                        matrix_Y_unordered_train_extern,
                                                                        matrix_Y_ordered_train_extern,
                                                                        matrix_Y_continuous_train_extern,
                                                                        /* Eval */
                                                                        matrix_Y_unordered_eval_extern,
                                                                        matrix_Y_ordered_eval_extern,
                                                                        matrix_Y_continuous_eval_extern,
                                                                        /* Train */
                                                                        matrix_X_unordered_train_extern,
                                                                        matrix_X_ordered_train_extern,
                                                                        matrix_X_continuous_train_extern,
                                                                        /* Eval */
                                                                        matrix_X_unordered_eval_extern,
                                                                        matrix_X_ordered_eval_extern,
                                                                        matrix_X_continuous_eval_extern,
                                                                        &vector_scale_factor[1],
                                                                        matrix_categorical_vals_extern,
                                                                        num_categories_extern,
                                                                        pdf,
                                                                        &pdf_deriv[num_reg_continuous_extern],
                                                                        &pdf_deriv_stderr[num_reg_continuous_extern],
                                                                        small, itmax);
    }
  }

  /* return data to R */
  for( i=0; i<num_obs_eval_extern; i++ )
    cdens[i]=pdf[i];

  for( i=0; i<num_obs_eval_extern; i++ )
    cderr[i]=pdf_stderr[i];

  if (do_grad) {
    for(j=0;j<num_var;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        cgerr[j*num_obs_eval_extern+i]=pdf_deriv_stderr[j][i];

    for(j=0;j<num_var;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        cg[j*num_obs_eval_extern+i]=pdf_deriv[j][i];
  }



  *ll = log_likelihood;
  /* end return data */

  /* Free data objects */

  free_mat(matrix_Y_unordered_train_extern, num_var_unordered_extern);
  free_mat(matrix_Y_ordered_train_extern, num_var_ordered_extern);
  free_mat(matrix_Y_continuous_train_extern, num_var_continuous_extern);

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  if(train_is_eval){
    matrix_Y_unordered_eval_extern = NULL;
    matrix_Y_ordered_eval_extern = NULL;
    matrix_Y_continuous_eval_extern = NULL;

    matrix_X_unordered_eval_extern = NULL;
    matrix_X_ordered_eval_extern = NULL;
    matrix_X_continuous_eval_extern = NULL;
  } else {
    free_mat(matrix_Y_unordered_eval_extern, num_var_unordered_extern);
    free_mat(matrix_Y_ordered_eval_extern, num_var_ordered_extern);
    free_mat(matrix_Y_continuous_eval_extern, num_var_continuous_extern);

    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  if (do_grad){
    free_mat(pdf_deriv, num_var);
    free_mat(pdf_deriv_stderr, num_var);
  }

  safe_free(vector_scale_factor);
  safe_free(num_categories_extern);
  safe_free(pdf);
  safe_free(pdf_stderr);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);
  return ;
}


void np_density(double * tuno, double * tord, double * tcon, 
                double * euno, double * eord, double * econ, 
                double * dbw, 
                double * mcv, double * padnum, 
                int * myopti, double * mydens, double * myderr, double * ll){


  double small = 1.0e-16;
  double * vector_scale_factor, * pdf, * pdf_stderr, log_likelihood = 0.0;
  double pad_num;

  int itmax = 10000;
  int i,j;
  int num_var, num_obs_eval_alloc, max_lev, train_is_eval, dens_or_dist;
  

  /* match integer options with their globals */

  num_reg_continuous_extern = myopti[DEN_NCONI];
  num_reg_unordered_extern = myopti[DEN_NUNOI];
  num_reg_ordered_extern = myopti[DEN_NORDI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[DEN_TNOBSI];
  num_obs_eval_extern = myopti[DEN_ENOBSI];

  KERNEL_den_extern = myopti[DEN_CKRNEVI];
  KERNEL_den_unordered_extern = 0;
  KERNEL_den_ordered_extern = 0;

  int_LARGE_SF = myopti[DEN_LSFI];
  int_MINIMIZE_IO = myopti[DEN_MINIOI];
  BANDWIDTH_den_extern = myopti[DEN_DENI];

  train_is_eval = myopti[DEN_TISEI];

  max_lev = myopti[DEN_MLEVI];
  pad_num = *padnum;

  dens_or_dist = myopti[DEN_DODENI];

#ifdef MPI2
  num_obs_eval_alloc = MAX(ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif

  /* Allocate memory for objects */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  if(!train_is_eval){
    matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
    matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
    matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);
  } else {
    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;
  }

  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  vector_scale_factor = alloc_vecd(num_var + 1);
  matrix_categorical_vals_extern = alloc_matd(max_lev, num_reg_unordered_extern + num_reg_ordered_extern);

  /* note use of num_obs_eval_alloc */
  pdf = alloc_vecd(num_obs_eval_alloc);
  pdf_stderr = alloc_vecd(num_obs_eval_alloc);

  
  /* Parse data */
	
  /* train */

  for( j=0;j<num_reg_unordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+i];

  /* eval */
  if (!train_is_eval) {
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+i];
  }

  /*  bandwidths/scale factors */

  for( i=0; i<num_var; i++ )
    vector_scale_factor[i+1]=dbw[i];

  /* fix up categories */
  
  for(j=0; j < (num_reg_unordered_extern + num_reg_ordered_extern); j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;
  }

  /* Conduct estimation */
  
  if (dens_or_dist == NP_DO_DENS){
  /* nb - KERNEL_(|un)ordered_den are set to zero upon declaration 
     - they have only one kernel type each at the moment */
    kernel_estimate_density_categorical(KERNEL_den_extern,
                                        KERNEL_den_unordered_extern,
                                        KERNEL_den_ordered_extern,
                                        BANDWIDTH_den_extern,
                                        num_obs_train_extern,
                                        num_obs_eval_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        num_reg_continuous_extern,
                                        /* Train */
                                        matrix_X_unordered_train_extern,
                                        matrix_X_ordered_train_extern,
                                        matrix_X_continuous_train_extern,
                                        /* Eval */
                                        matrix_X_unordered_eval_extern,
                                        matrix_X_ordered_eval_extern,
                                        matrix_X_continuous_eval_extern,
                                        &vector_scale_factor[1],
                                        num_categories_extern,
                                        pdf,
                                        pdf_stderr,
                                        &log_likelihood);
  } else if (dens_or_dist == NP_DO_DIST) {
    kernel_estimate_distribution_categorical(KERNEL_den_extern,
                                             KERNEL_den_unordered_extern,
                                             KERNEL_den_ordered_extern,
                                             BANDWIDTH_den_extern,
                                             num_obs_train_extern,
                                             num_obs_eval_extern,
                                             num_reg_unordered_extern,
                                             num_reg_ordered_extern,
                                             num_reg_continuous_extern,
                                             /* Train */
                                             matrix_X_unordered_train_extern,
                                             matrix_X_ordered_train_extern,
                                             matrix_X_continuous_train_extern,
                                             /* Eval */
                                             matrix_X_unordered_eval_extern,
                                             matrix_X_ordered_eval_extern,
                                             matrix_X_continuous_eval_extern,
                                             &vector_scale_factor[1],
                                             num_categories_extern,
                                             matrix_categorical_vals_extern,
                                             pdf,
                                             pdf_stderr,
                                             small, itmax);

  }
  
  
  /* write the return values */

  for(i=0;i<num_obs_eval_extern;i++){
    mydens[i] = pdf[i];
    myderr[i] = pdf_stderr[i];
  }
  *ll = log_likelihood;

  /* clean up and wave goodbye */

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  if (!train_is_eval){
    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  safe_free(vector_scale_factor);
  safe_free(num_categories_extern);
  safe_free(pdf_stderr);
  safe_free(pdf);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  return;
}


void np_regression_bw(double * runo, double * rord, double * rcon, double * y,
                      int * myopti, double * myoptd, double * rbw, double * fval){
  
  double **matrix_y;

  double *vector_continuous_stddev;
  double *vector_scale_factor, *vector_scale_factor_multistart;

  double fret, fret_best;
  double ftol, tol, small;
  double (* bwmfunc)(double *);

  int i,j;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, iImproved;
  int itmax, iter;
  int int_use_starting_values;

  num_reg_continuous_extern = myopti[RBW_NCONI];
  num_reg_unordered_extern = myopti[RBW_NUNOI];
  num_reg_ordered_extern = myopti[RBW_NORDI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[RBW_NOBSI];
  iMultistart = myopti[RBW_IMULTII];
  iNum_Multistart = myopti[RBW_NMULTII];

  KERNEL_reg_extern = myopti[RBW_CKRNEVI];
  KERNEL_reg_unordered_extern = myopti[RBW_UKRNEVI];
  KERNEL_reg_ordered_extern = myopti[RBW_OKRNEVI];

  int_use_starting_values= myopti[RBW_USTARTI];
  int_LARGE_SF=myopti[RBW_LSFI];

  BANDWIDTH_reg_extern=myopti[RBW_REGI];
  BANDWIDTH_den_extern=0;

  itmax=myopti[RBW_ITMAXI];
  int_RESTART_FROM_MIN = myopti[RBW_REMINI];
  int_MINIMIZE_IO = myopti[RBW_MINIOI];

  int_ll_extern = myopti[RBW_LL];

  ftol=myoptd[RBW_FTOLD];
  tol=myoptd[RBW_TOLD];
  small=myoptd[RBW_SMALLD];

  imsnum = 0;
  imstot = iNum_Multistart;

  /* Allocate memory for objects */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  vector_Y_extern = alloc_vecd(num_obs_train_extern);
	
  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_y = alloc_matd(num_var + 1, num_var +1);
  vector_scale_factor = alloc_vecd(num_var + 1);
  matrix_categorical_vals_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  vector_continuous_stddev = alloc_vecd(num_reg_continuous_extern);

  /* Request starting values for optimization if values already exist */

  /* bandwidths */

  if (int_use_starting_values)
    for( i=0;i<num_var; i++ )
      vector_scale_factor[i+1] = rbw[i];

  /* regressors */

  for( j=0;j<num_reg_unordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_unordered_train_extern[j][i]=runo[j*num_obs_train_extern+i];
    

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=rord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=rcon[j*num_obs_train_extern+i];

  /* response variable */
  for( i=0;i<num_obs_train_extern;i++ )
    vector_Y_extern[i] = y[i];

  determine_categorical_vals(
                             num_obs_train_extern,
                             0,
                             0,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             matrix_Y_unordered_train_extern,
                             matrix_Y_ordered_train_extern,
                             matrix_X_unordered_train_extern,
                             matrix_X_ordered_train_extern,
                             num_categories_extern,
                             matrix_categorical_vals_extern);


  compute_continuous_stddev(
                            int_LARGE_SF,
                            num_obs_train_extern,
                            0,
                            num_reg_continuous_extern,
                            matrix_Y_continuous_train_extern,
                            matrix_X_continuous_train_extern,
                            vector_continuous_stddev);

  /* Initialize scale factors and Hessian for NR modules */

  initialize_nr_vector_scale_factor(
                                    BANDWIDTH_reg_extern,
                                    BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    0,                /* regression (0) regression ml (1) */
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    matrix_Y_continuous_train_extern,
                                    matrix_X_continuous_train_extern,
                                    int_use_starting_values,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor);

  initialize_nr_hessian(num_var, matrix_y);

  /* When multistarting, set counter */

  iMs_counter = 0;

  /* assign the function to be optimized */
  switch(myopti[RBW_MI]){
  case RBWM_CVAIC : bwmfunc = cv_func_regression_categorical_aic_c; break;
  case RBWM_CVLS : bwmfunc = cv_func_regression_categorical_ls; break;
  default : REprintf("np.c: invalid bandwidth selection method.");
    exit(0);break;
  }

  spinner(0);

  fret_best = bwmfunc(vector_scale_factor);
  iImproved = 0;

  powell(0,
         0,
         vector_scale_factor,
         vector_scale_factor,
         matrix_y,
         num_var,
         ftol,
         tol,
         small,
         itmax,
         &iter,
         &fret,
         bwmfunc);


  if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

    initialize_nr_hessian(num_var, matrix_y);

    powell(0,
           0,
           vector_scale_factor,
           vector_scale_factor,
           matrix_y,
           num_var,
           ftol,
           tol,
           small,
           itmax,
           &iter,
           &fret,
           bwmfunc);

  }

  iImproved = (fret < fret_best);

  /* When multistarting save initial minimum of objective function and scale factors */


  if(iMultistart == IMULTI_TRUE){
    fret_best = fret;
    vector_scale_factor_multistart = alloc_vecd(num_var + 1);

    for(i = 1; i <= num_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];

    /* Conduct search from new random values of the search parameters */

    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and hessian for NR modules */
				
      initialize_nr_vector_scale_factor(BANDWIDTH_reg_extern,
                                        BANDWIDTH_den_extern,
                                        1,        /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        0,        /* regression (0) regression ml (1) */
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        0,
                                        0,
                                        0,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        matrix_Y_continuous_train_extern,
                                        matrix_X_continuous_train_extern,
                                        int_use_starting_values,
                                        pow((double)4.0/(double)3.0,0.2),     /* Init for continuous vars */
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor);
      initialize_nr_hessian(num_var, matrix_y);

      /* Conduct direction set search */

      powell(0,
             0,
             vector_scale_factor,
             vector_scale_factor,
             matrix_y,
             num_var,
             ftol,
             tol,
             small,
             itmax,
             &iter,
             &fret,
             bwmfunc);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE)	{
						
        initialize_nr_hessian(num_var, matrix_y);
						
        powell(0,
               0,
               vector_scale_factor,
               vector_scale_factor,
               matrix_y,
               num_var,
               ftol,
               tol,
               small,
               itmax,
               &iter,
               &fret,
               bwmfunc);

      }

      /* If this run resulted in an improved minimum save information */
      
      if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        
        for(i = 1; i <= num_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }


    }

    /* Save best for estimation */

    fret = fret_best;

    for(i = 1; i <= num_var; i++)
      vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];

    free(vector_scale_factor_multistart);

  }

  /* return data to R */
  if (BANDWIDTH_reg_extern == BW_GEN_NN || 
      BANDWIDTH_reg_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern; i++ )
      vector_scale_factor[i+1]=fround(vector_scale_factor[i+1]);
  }
  for( i=0; i<num_var; i++ )
    rbw[i]=vector_scale_factor[i+1];

  fval[0] = fret;
  fval[1] = iImproved;
  /* end return data */

  /* Free data objects */

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  safe_free(vector_Y_extern);

  free_mat(matrix_y, num_var + 1);
  safe_free(vector_scale_factor);
  safe_free(num_categories_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  free(vector_continuous_stddev);

  if(int_MINIMIZE_IO != IO_MIN_TRUE)
    Rprintf("\r                   \r");

  //fprintf(stderr,"\nNP TOASTY\n");
  return ;
  
}


void np_regression(double * tuno, double * tord, double * tcon, double * ty,
                   double * euno, double * eord, double * econ, double * ey,
                   double * rbw, 
                   double * mcv, double * padnum, 
                   int * myopti, 
                   double * cm, double * cmerr, double * g, double *gerr, 
                   double * xtra){

  double * vector_scale_factor, * ecm, * ecmerr, ** eg, **egerr;
  double * lambda, ** matrix_bandwidth;
  double RS, MSE, MAE, MAPE, CORR, SIGN, pad_num;

  int i,j, num_var;
  int ey_is_ty, do_grad, train_is_eval, num_obs_eval_alloc, max_lev;

  /* match integer options with their globals */

  num_reg_continuous_extern = myopti[REG_NCONI];
  num_reg_unordered_extern = myopti[REG_NUNOI];
  num_reg_ordered_extern = myopti[REG_NORDI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  train_is_eval = myopti[REG_TISEI];
  ey_is_ty = myopti[REG_EY];

  num_obs_train_extern = myopti[REG_TNOBSI];
  num_obs_eval_extern = myopti[REG_ENOBSI];

  if(train_is_eval && (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_regression): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    exit(0);
  }

  KERNEL_reg_extern = myopti[REG_CKRNEVI];
  KERNEL_reg_unordered_extern = myopti[REG_UKRNEVI];
  KERNEL_reg_ordered_extern = myopti[REG_OKRNEVI];

  int_LARGE_SF = myopti[REG_LSFI];
  int_MINIMIZE_IO = myopti[REG_MINIOI];
  BANDWIDTH_reg_extern = myopti[REG_BWI];

  do_grad = myopti[REG_GRAD];
  int_ll_extern = myopti[REG_LL];

  max_lev = myopti[REG_MLEVI];
  pad_num = *padnum;

#ifdef MPI2
  num_obs_eval_alloc = MAX(ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif


  /* Allocate memory for objects */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  vector_Y_extern = alloc_vecd(num_obs_train_extern);

  if(!train_is_eval){
    matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
    matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
    matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);

    if(!ey_is_ty)
      vector_Y_eval_extern = alloc_vecd(num_obs_eval_extern);
    else
      vector_Y_eval_extern = NULL;

  } else {
    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;

    if(!ey_is_ty)
      vector_Y_eval_extern = alloc_vecd(num_obs_eval_extern);
    else
      vector_Y_eval_extern = vector_Y_extern;

  }

  ecm = alloc_vecd(num_obs_eval_alloc);
  ecmerr = alloc_vecd(num_obs_eval_alloc);
  
  eg = alloc_matd(num_obs_eval_alloc, num_var);
  egerr = alloc_matd(num_obs_eval_alloc, num_var);

  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  vector_scale_factor = alloc_vecd(num_var + 1);
  matrix_categorical_vals_extern = alloc_matd(max_lev, num_reg_unordered_extern + num_reg_ordered_extern);

  lambda =  alloc_vecd(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_bandwidth = alloc_matd((BANDWIDTH_reg_extern==BW_GEN_NN)?num_obs_eval_extern:
                                ((BANDWIDTH_reg_extern==BW_ADAP_NN)?num_obs_train_extern:1),num_reg_continuous_extern);  

  /* train */

  for( j=0;j<num_reg_unordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+i];

  for( i=0;i<num_obs_train_extern;i++ )
    vector_Y_extern[i] = ty[i];

  /* eval */
  if(!train_is_eval){
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+i];
  }

  if(!ey_is_ty)
    for(i=0;i<num_obs_eval_extern;i++)
      vector_Y_eval_extern[i] = ey[i];

  /*  bandwidths/scale factors */

  for( i=0; i<num_var; i++ )
    vector_scale_factor[i+1] = rbw[i];

  /* fix up categories */

  for(j=0; j < (num_reg_unordered_extern + num_reg_ordered_extern); j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;
  }


  /* Conduct estimation */
	
  /* 
     nb - KERNEL_(|un)ordered_den are set to zero upon declaration 
     - they have only one kernel type each at the moment 
  */

  kernel_estimate_regression_categorical(int_ll_extern,
                                         KERNEL_reg_extern,
                                         KERNEL_reg_unordered_extern,
                                         KERNEL_reg_ordered_extern,
                                         BANDWIDTH_reg_extern,
                                         num_obs_train_extern,
                                         num_obs_eval_extern,
                                         num_reg_unordered_extern,
                                         num_reg_ordered_extern,
                                         num_reg_continuous_extern,
                                         /* Train */
                                         matrix_X_unordered_train_extern,
                                         matrix_X_ordered_train_extern,
                                         matrix_X_continuous_train_extern,
                                         /* Eval */
                                         matrix_X_unordered_eval_extern,
                                         matrix_X_ordered_eval_extern,
                                         matrix_X_continuous_eval_extern,
                                         /* Bandwidth */
                                         matrix_X_continuous_train_extern,
                                         vector_Y_extern,
                                         vector_Y_eval_extern,
                                         &vector_scale_factor[1],
                                         num_categories_extern,
                                         ecm,
                                         eg,
                                         ecmerr,
                                         egerr,
                                         &RS,
                                         &MSE,
                                         &MAE,
                                         &MAPE,
                                         &CORR,
                                         &SIGN);

  if (do_grad){
    kernel_bandwidth_mean(KERNEL_reg_extern,
                          BANDWIDTH_reg_extern,
                          num_obs_train_extern,
                          num_obs_eval_extern,
                          0,
                          0,
                          0,
                          num_reg_continuous_extern,
                          num_reg_unordered_extern,
                          num_reg_ordered_extern,
                          &vector_scale_factor[1],
                          /* Not used */
                          matrix_Y_continuous_train_extern,
                          /* Not used */
                          matrix_Y_continuous_train_extern,
                          matrix_X_continuous_train_extern,
                          matrix_X_continuous_eval_extern,
                          matrix_bandwidth,/* Not used */
                          matrix_bandwidth,
                          lambda);
    kernel_estimate_categorical_gradient_ocg_fast(1,
                                                  NULL,
                                                  0,
                                                  KERNEL_reg_extern,
                                                  KERNEL_reg_unordered_extern,
                                                  KERNEL_reg_ordered_extern,
                                                  BANDWIDTH_reg_extern,
                                                  int_ll_extern,
                                                  0,
                                                  num_obs_train_extern,
                                                  num_obs_eval_extern,
                                                  num_reg_unordered_extern,
                                                  num_reg_ordered_extern,
                                                  num_reg_continuous_extern,
                                                  vector_Y_extern,
                                                  matrix_X_unordered_train_extern,
                                                  matrix_X_ordered_train_extern,
                                                  matrix_X_continuous_train_extern,
                                                  matrix_X_unordered_eval_extern,
                                                  matrix_X_ordered_eval_extern,
                                                  matrix_X_continuous_eval_extern,
                                                  matrix_bandwidth,
                                                  NULL,
                                                  lambda,
                                                  num_categories_extern,
                                                  matrix_categorical_vals_extern,
                                                  ecm,
                                                  &eg[num_reg_continuous_extern]);

  }

  /* write the return values */

  for(i=0;i<num_obs_eval_extern;i++)
    cm[i] = ecm[i];

  for(i=0;i<num_obs_eval_extern;i++)
    cmerr[i] = ecmerr[i];

  if(do_grad){
    for(j=0;j<num_var;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        g[j*num_obs_eval_extern+i]=eg[j][i];

    for(j=0;j<num_reg_continuous_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        gerr[j*num_obs_eval_extern+i]=egerr[j][i];
  }

  xtra[0] = RS;
  xtra[1] = MSE;
  xtra[2] = MAE;
  xtra[3] = MAPE;
  xtra[4] = CORR;
  xtra[5] = SIGN;

  /* clean up and wave goodbye */

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  if(!train_is_eval){
    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  free_mat(eg, num_var);
  free_mat(egerr, num_var);

  free_mat(matrix_bandwidth, num_reg_continuous_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  safe_free(vector_Y_extern);
  if(!ey_is_ty)
    safe_free(vector_Y_eval_extern);

  safe_free(ecm);
  safe_free(ecmerr);

  safe_free(num_categories_extern);
  safe_free(vector_scale_factor);

  safe_free(lambda);

  return;
}


void np_kernelsum(double * tuno, double * tord, double * tcon, 
                  double * ty, double * weights,
                  double * euno, double * eord, double * econ, 
                  double * bw,
                  double * mcv, double * padnum, 
                  int * myopti, double * kpow, double * weighted_sum){
      
  /* the ys are the weights */

  double * vector_scale_factor, * ksum, pad_num;
  int i,j, num_var, num_obs_eval_alloc;
  int no_y, do_ipow, leave_one_out, train_is_eval, do_divide_bw, operator;
  int max_lev, do_smooth_coef_weights, no_weights, sum_element_length;


  /* match integer options with their globals */

  num_reg_continuous_extern = myopti[KWS_NCONI];
  num_reg_unordered_extern = myopti[KWS_NUNOI];
  num_reg_ordered_extern = myopti[KWS_NORDI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[KWS_TNOBSI];
  num_obs_eval_extern = myopti[KWS_ENOBSI];

  KERNEL_reg_extern = myopti[KWS_CKRNEVI];
  KERNEL_reg_unordered_extern = myopti[KWS_UKRNEVI];
  KERNEL_reg_ordered_extern = myopti[KWS_OKRNEVI];

  int_LARGE_SF = myopti[KWS_LSFI];
  int_MINIMIZE_IO = myopti[KWS_MINIOI];
  BANDWIDTH_reg_extern = myopti[KWS_BWI];

  train_is_eval = myopti[KWS_TISEI];
  // no_y = myopti[KWS_NOYI];
  leave_one_out = myopti[KWS_LOOI];
  do_ipow = myopti[KWS_IPOWI];
  do_divide_bw = myopti[KWS_BDIVI];
  operator = myopti[KWS_OPI];

  max_lev = myopti[KWS_MLEVI];
  pad_num = *padnum;

  // do_weights = myopti[KWS_DOWI];
  do_smooth_coef_weights = myopti[KWS_SCOEFI];

  /* the y and weight matrices will be contained in these variables */
  num_var_continuous_extern = myopti[KWS_YNCOLI];
  num_var_ordered_extern = myopti[KWS_WNCOLI];

  no_y = (num_var_continuous_extern == 0);
  no_weights = (num_var_ordered_extern == 0);

  sum_element_length = (no_y ? 1 : num_var_continuous_extern)*(no_weights ? 1 : num_var_ordered_extern);

#ifdef MPI2
  num_obs_eval_alloc = MAX(ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif

  if(train_is_eval && (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_kernelsum): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    exit(0);
  }

  /* allocate */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);
  
  /* for the moment we will just allocate a vector of ones */
  /* vector_Y_extern = (no_y)?NULL:alloc_vecd(num_obs_train_extern); */

  matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);
  matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);

  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_categorical_vals_extern = alloc_matd(max_lev, num_reg_unordered_extern + num_reg_ordered_extern);

  vector_scale_factor = alloc_vecd(num_var + 1);
  ksum = alloc_vecd(num_obs_eval_alloc*sum_element_length);

  if(!train_is_eval){
    matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
    matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
    matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);
  } else {
    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;
  }

  /* train */

  for( j=0;j<num_reg_unordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+i];

  for( j = 0; j < num_var_continuous_extern; j++ )
    for( i = 0; i < num_obs_train_extern; i++ )
      matrix_Y_continuous_train_extern[j][i] = ty[j*num_obs_train_extern+i];

  for( j = 0; j < num_var_ordered_extern; j++ )
    for( i = 0; i < num_obs_train_extern; i++ )
      matrix_Y_ordered_train_extern[j][i] = weights[j*num_obs_train_extern+i];

  /* eval */

  if(!train_is_eval){
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+i];
  }


  /* bandwidths */
  for( i=0; i<num_var; i++ )
    vector_scale_factor[i+1] = bw[i];
  
  /* fix up categories */

  for(j=0; j<num_reg_unordered_extern; j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;
  }

  for(j=num_reg_unordered_extern; j < (num_reg_unordered_extern+num_reg_ordered_extern); j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;
  }

  assert(!((operator == OP_CONVOLUTION) && (BANDWIDTH_reg_extern != BW_ADAP_NN) && (KERNEL_reg_extern == 8)));

  
  kernel_weighted_sum_np(KERNEL_reg_extern,
                         KERNEL_reg_unordered_extern,
                         KERNEL_reg_ordered_extern,
                         BANDWIDTH_reg_extern,
                         num_obs_train_extern,
                         num_obs_eval_extern,
                         num_reg_unordered_extern,
                         num_reg_ordered_extern,
                         num_reg_continuous_extern,
                         leave_one_out,
                         (int)(*kpow),
                         do_divide_bw,
                         do_smooth_coef_weights,
                         0, //not symmetric
                         0, //disable 'twisting'
                         operator,
                         matrix_X_unordered_train_extern,
                         matrix_X_ordered_train_extern,
                         matrix_X_continuous_train_extern,
                         matrix_X_unordered_eval_extern,
                         matrix_X_ordered_eval_extern,
                         matrix_X_continuous_eval_extern,
                         /* ys matrix */
                         matrix_Y_continuous_train_extern,
                         /* weights matrix */
                         matrix_Y_ordered_train_extern,
                         NULL,
                         &vector_scale_factor[1],
                         num_categories_extern,
                         matrix_categorical_vals_extern,
                         ksum);
  /*
    kernel_convolution_weighted_sum(KERNEL_reg_extern,
                                    KERNEL_reg_unordered_extern,
                                    KERNEL_reg_ordered_extern,
                                    BANDWIDTH_reg_extern,
                                    num_obs_train_extern,
                                    num_obs_eval_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    num_reg_continuous_extern,
                                    matrix_X_unordered_train_extern,
                                    matrix_X_ordered_train_extern,
                                    matrix_X_continuous_train_extern,
                                    matrix_X_unordered_eval_extern,
                                    matrix_X_ordered_eval_extern,
                                    matrix_X_continuous_eval_extern,
                                    vector_Y_extern,
                                    &vector_scale_factor[1],
                                    num_categories_extern,
                                    matrix_categorical_vals_extern,
                                    ksum);
  */


  for(i = 0; i < sum_element_length * num_obs_eval_extern; i++)
    weighted_sum[i] = ksum[i];

  /* clean up */

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  
  if(!train_is_eval){
    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);
  free_mat(matrix_Y_continuous_train_extern, num_var_continuous_extern);
  free_mat(matrix_Y_ordered_train_extern, num_var_ordered_extern);

  safe_free(num_categories_extern);
  safe_free(vector_scale_factor);
  safe_free(ksum);

  return;
}


void np_quantile_conditional(double * tc_con,
                             double * tu_uno, double * tu_ord, double * tu_con,
                             double * eu_uno, double * eu_ord, double * eu_con,
                             double * quantile,
                             double * mybw, 
                             double * mcv, double *padnum,
                             int * myopti, double * myoptd,
                             double * yq, double * yqerr, double *yg){
  /* Likelihood bandwidth selection for density estimation */

  double **g = NULL, * eq, * eqerr;
  double ftol, small, tol, itmax;
  double pad_num;

  int i,j, max_lev;
  int num_var, num_obs_eval_alloc;
  int num_all_var, num_var_var, train_is_eval, do_gradients;

  imsnum = 0;
  imstot = myopti[BW_NMULTII]; /* iNum_Multistart */

  num_var_unordered_extern = 0;
  num_var_ordered_extern = 0;
  num_var_continuous_extern = 1;

  num_reg_unordered_extern = myopti[CQ_UNUNOI];
  num_reg_ordered_extern = myopti[CQ_UNORDI];
  num_reg_continuous_extern = myopti[CQ_UNCONI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;
  num_var_var = 1;
  num_all_var = num_var + num_var_var;

  num_obs_train_extern = myopti[CQ_TNOBSI];
  num_obs_eval_extern = myopti[CQ_ENOBSI];

  if((train_is_eval = myopti[CQ_TISEI]) && 
     (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_quantile_conditional): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    exit(0);
  }

  KERNEL_reg_extern = myopti[CQ_CXKRNEVI];
  KERNEL_den_extern = myopti[CQ_CYKRNEVI];

  KERNEL_reg_unordered_extern = myopti[CQ_UXKRNEVI];
  KERNEL_den_unordered_extern = myopti[CQ_UYKRNEVI];

  KERNEL_reg_ordered_extern = myopti[CQ_OXKRNEVI];
  KERNEL_den_ordered_extern = myopti[CQ_OYKRNEVI];

  int_LARGE_SF = myopti[CQ_LSFI];
  BANDWIDTH_den_extern = myopti[CQ_DENI];
  int_MINIMIZE_IO = myopti[CQ_MINIOI];
  do_gradients = myopti[CQ_GRADI];
  itmax = myopti[CQ_ITMAXI];

  max_lev = myopti[CQ_MLEVI];
  pad_num = *padnum;

  ftol = myoptd[CQ_FTOLD];
  tol = myoptd[CQ_TOLD];
  small = myoptd[CQ_SMALLD];

  gamma_extern = *quantile;

#ifdef MPI2
  num_obs_eval_alloc = MAX(ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif


  /* Allocate memory for objects */
  matrix_Y_continuous_quantile_extern = alloc_matd(1, num_var_continuous_extern);
  matrix_X_unordered_quantile_extern = alloc_matd(1, num_reg_unordered_extern);
  matrix_X_ordered_quantile_extern = alloc_matd(1, num_reg_ordered_extern);
  matrix_X_continuous_quantile_extern = alloc_matd(1, num_reg_continuous_extern);
  /* */

  matrix_Y_unordered_train_extern = alloc_matd(num_obs_train_extern, num_var_unordered_extern);
  matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);
  matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  if(train_is_eval) {
    matrix_Y_unordered_eval_extern = alloc_matd(num_obs_train_extern, num_var_unordered_extern);
    matrix_Y_ordered_eval_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);
    matrix_Y_continuous_eval_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);

    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;
  } else {
    matrix_Y_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_var_unordered_extern);
    matrix_Y_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_var_ordered_extern);
    matrix_Y_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_var_continuous_extern);

    matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
    matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
    matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);
  }
	
  num_categories_extern = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                     num_reg_unordered_extern + num_reg_ordered_extern);
  vector_scale_factor_extern = alloc_vecd(num_all_var + 1);
  
  matrix_categorical_vals_extern = 
    alloc_matd(max_lev, num_var_unordered_extern + num_var_ordered_extern + 
               num_reg_unordered_extern + num_reg_ordered_extern);

  eq = alloc_vecd(num_obs_eval_alloc);
  eqerr = alloc_vecd(num_obs_eval_alloc);

  if(do_gradients)
    g = alloc_matd(num_obs_eval_alloc, num_var);

  /* in v_s_f order is creg, cvar, uvar, ovar, ureg, oreg  */

  for( i=0;i<num_all_var; i++ )
    vector_scale_factor_extern[i+1] = mybw[i];

  /* Parse data */

  /* train */

  for(j=0;j<num_var_continuous_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_Y_continuous_train_extern[j][i]=tc_con[j*num_obs_train_extern+i];

  for(j=0;j<num_reg_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_X_unordered_train_extern[j][i]=tu_uno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_ordered_train_extern[j][i]=tu_ord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_X_continuous_train_extern[j][i]=tu_con[j*num_obs_train_extern+i];

  /* eval */
  if(!train_is_eval){
    for(j=0;j<num_reg_unordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_X_unordered_eval_extern[j][i]=eu_uno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=eu_ord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=eu_con[j*num_obs_eval_extern+i];
  }

  /* fix up categories */
  
  for(j=0; j < (num_reg_unordered_extern + num_reg_ordered_extern); j++){
    i = 0;
    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;
  }

  kernel_estimate_quantile(do_gradients,
                           KERNEL_den_extern,
                           KERNEL_den_unordered_extern,
                           KERNEL_den_ordered_extern,
                           BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_obs_eval_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           num_var_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           num_reg_continuous_extern,
                           matrix_Y_unordered_train_extern,
                           matrix_Y_ordered_train_extern,
                           matrix_Y_continuous_train_extern,
                           matrix_Y_unordered_eval_extern,
                           matrix_Y_ordered_eval_extern,
                           matrix_Y_continuous_eval_extern,
                           matrix_X_unordered_train_extern,
                           matrix_X_ordered_train_extern,
                           matrix_X_continuous_train_extern,
                           matrix_X_unordered_eval_extern,
                           matrix_X_ordered_eval_extern,
                           matrix_X_continuous_eval_extern,
                           &vector_scale_factor_extern[1],
                           eq,
                           eqerr,
                           g,
                           int_RANDOM_SEED,
                           ftol,
                           tol,
                           small,
                           itmax,
                           itmax,            /* Maximum number of multistarts */
                           1.0e-10);         /* Zero for all intents and purposes */


  /* return data to R */

  for(i=0; i < num_obs_eval_extern; i++)
    yq[i] = eq[i];

  for(i=0; i < num_obs_eval_extern; i++)
    yqerr[i] = eqerr[i];

  if(do_gradients)
    for(j=0; j < num_var; j++)
      for(i=0; i < num_obs_eval_extern; i++)
        yg[j*num_obs_eval_extern+i] = g[j][i];
  
  /* end return data */

  /* Free data objects */

  free_mat(matrix_Y_unordered_train_extern, num_var_unordered_extern);
  free_mat(matrix_Y_ordered_train_extern, num_var_ordered_extern);
  free_mat(matrix_Y_continuous_train_extern, num_var_continuous_extern);

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);

  free_mat(matrix_Y_continuous_quantile_extern, num_var_continuous_extern); 
  free_mat(matrix_X_unordered_quantile_extern, num_reg_unordered_extern); 
  free_mat(matrix_X_ordered_quantile_extern, num_reg_ordered_extern); 
  free_mat(matrix_X_continuous_quantile_extern, num_reg_continuous_extern); 

  free_mat(g, num_var);

  free_mat(matrix_Y_unordered_eval_extern, num_var_unordered_extern);
  free_mat(matrix_Y_ordered_eval_extern, num_var_ordered_extern);
  free_mat(matrix_Y_continuous_eval_extern, num_var_continuous_extern);


  if(train_is_eval){
    matrix_X_unordered_eval_extern = NULL;
    matrix_X_ordered_eval_extern = NULL;
    matrix_X_continuous_eval_extern = NULL;
  } else {
    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  safe_free(vector_scale_factor_extern);
  safe_free(num_categories_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);

  safe_free(eq);
  safe_free(eqerr);

  if(int_MINIMIZE_IO != IO_MIN_TRUE)
    Rprintf("\r                   \r");
  return ;
}

