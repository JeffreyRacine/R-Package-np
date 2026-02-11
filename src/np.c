/* Copyright (C) Jeff Racine, 1995-2004 */

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <R.h>
#include <R_ext/Arith.h>
#include <stdint.h>

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
#include "tree.h"

// categorical hashing
#include "hash.h"

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

int int_TREE_X;
int int_TREE_Y;
int int_TREE_XY;

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

/* these are data which are sorted into an 'alternate' order */
/* this allows us to support 2 trees simultaneously !*/

int * num_categories_extern_XY;
int * num_categories_extern_X;
int * num_categories_extern_Y;

double **matrix_categorical_vals_extern_XY;
double **matrix_categorical_vals_extern_X;
double **matrix_categorical_vals_extern_Y;

double **matrix_XY_continuous_train_extern;
double **matrix_XY_unordered_train_extern;
double **matrix_XY_ordered_train_extern;
double **matrix_XY_continuous_eval_extern;
double **matrix_XY_unordered_eval_extern;
double **matrix_XY_ordered_eval_extern;

int * ipt_extern_X;
int * ipt_extern_Y;
int * ipt_extern_XY;

static double (*bwmfunc_raw)(double *) = NULL;
static double bwm_eval_count = 0.0;
static double bwm_invalid_count = 0.0;
static int bwm_use_transform = 0;
static int bwm_num_reg_continuous = 0;
static int bwm_num_reg_unordered = 0;
static int bwm_num_reg_ordered = 0;
static int bwm_kernel_unordered = 0;
static int *bwm_kernel_unordered_vec = NULL;
static int bwm_kernel_unordered_len = 0;
static int *bwm_num_categories = NULL;
static double *bwm_transform_buf = NULL;
static int bwm_transform_buf_len = 0;
static int bwm_penalty_mode = 0;
static double bwm_penalty_value = DBL_MAX;

static void bwm_reset_counters(void)
{
  bwm_eval_count = 0.0;
  bwm_invalid_count = 0.0;
}

static double bwm_sigmoid(double x)
{
  if (x >= 0.0) {
    double z = exp(-x);
    return 1.0/(1.0+z);
  } else {
    double z = exp(x);
    return z/(1.0+z);
  }
}

static double bwm_safe_exp(double x)
{
  if (x > 700.0) return DBL_MAX/2.0;
  if (x < -700.0) return 0.0;
  return exp(x);
}

static double bwm_logit(double p)
{
  const double eps = 1e-12;
  if (p < eps) p = eps;
  if (p > 1.0 - eps) p = 1.0 - eps;
  return log(p/(1.0-p));
}

static void bwm_apply_transform(const double *p, double *out, int n)
{
  int i;
  int idx;

  for (i = 1; i <= n; i++)
    out[i] = p[i];

  for (i = 1; i <= bwm_num_reg_continuous; i++)
    out[i] = bwm_safe_exp(out[i]);

  for (i = 0; i < bwm_num_reg_unordered; i++) {
    idx = bwm_num_reg_continuous + 1 + i;
    if (bwm_num_categories != NULL) {
      int kern = (bwm_kernel_unordered_vec != NULL && i < bwm_kernel_unordered_len) ?
        bwm_kernel_unordered_vec[i] : bwm_kernel_unordered;
      double maxbw = max_unordered_bw(bwm_num_categories[i], kern);
      out[idx] = bwm_sigmoid(out[idx]) * maxbw;
    } else {
      out[idx] = bwm_sigmoid(out[idx]);
    }
  }

  for (i = 0; i < bwm_num_reg_ordered; i++) {
    idx = bwm_num_reg_continuous + bwm_num_reg_unordered + 1 + i;
    out[idx] = bwm_sigmoid(out[idx]);
  }
}

static void bwm_to_unconstrained(double *p, int n)
{
  int i;
  int idx;
  const double eps = 1e-12;

  if (!bwm_use_transform)
    return;

  for (i = 1; i <= n; i++)
    bwm_transform_buf[i] = p[i];

  for (i = 1; i <= bwm_num_reg_continuous; i++) {
    double v = bwm_transform_buf[i];
    if (v <= 0.0) v = eps;
    p[i] = log(v);
  }

  for (i = 0; i < bwm_num_reg_unordered; i++) {
    idx = bwm_num_reg_continuous + 1 + i;
    int kern = (bwm_kernel_unordered_vec != NULL && i < bwm_kernel_unordered_len) ?
      bwm_kernel_unordered_vec[i] : bwm_kernel_unordered;
    double maxbw = (bwm_num_categories != NULL) ? max_unordered_bw(bwm_num_categories[i], kern) : 1.0;
    double v = bwm_transform_buf[idx];
    if (maxbw <= 0.0) {
      p[idx] = 0.0;
    } else {
      double frac = v / maxbw;
      p[idx] = bwm_logit(frac);
    }
  }

  for (i = 0; i < bwm_num_reg_ordered; i++) {
    idx = bwm_num_reg_continuous + bwm_num_reg_unordered + 1 + i;
    p[idx] = bwm_logit(bwm_transform_buf[idx]);
  }
}

static void bwm_to_constrained(double *p, int n)
{
  int i;
  if (!bwm_use_transform)
    return;

  bwm_apply_transform(p, bwm_transform_buf, n);
  for (i = 1; i <= n; i++)
    p[i] = bwm_transform_buf[i];
}

static double bwmfunc_wrapper(double *p)
{
  double val;
  double *use_p = p;

  bwm_eval_count += 1.0;
  if (bwm_use_transform) {
    int n = bwm_num_reg_continuous + bwm_num_reg_unordered + bwm_num_reg_ordered;
    if (bwm_transform_buf_len < n + 1) {
      bwm_transform_buf = (double *) realloc(bwm_transform_buf, (n + 1) * sizeof(double));
      bwm_transform_buf_len = n + 1;
    }
    bwm_apply_transform(p, bwm_transform_buf, n);
    use_p = bwm_transform_buf;
  }

  val = bwmfunc_raw(use_p);

  if (!R_FINITE(val) || val == DBL_MAX) {
    bwm_invalid_count += 1.0;
    if (bwm_penalty_mode == 1 && R_FINITE(bwm_penalty_value))
      return bwm_penalty_value;
  }

  return val;
}

int * ipt_lookup_extern_X;
int * ipt_lookup_extern_Y;
int * ipt_lookup_extern_XY;

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

// cdf algorithm extern
double dbl_memfac_ccdf_extern = 1.0;
double dbl_memfac_dls_extern = 1.0;
int cdfontrain_extern = 0;

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

KDT * kdt_extern_X = NULL;
KDT * kdt_extern_Y = NULL;
KDT * kdt_extern_XY = NULL;

// to facilitate bandwidth->scale-factor conversions
 double * vector_continuous_stddev_extern = NULL;
double nconfac_extern = 0.0;
double ncatfac_extern = 0.0;

extern int iff;

double np_tgauss2_b = 3.0, np_tgauss2_alpha = 1.030174731161562;
double np_tgauss2_c0 = .004565578246317041;

double np_tgauss2_a0 = 0.2993759121518507, np_tgauss2_a1 = 2.0844504723243343E-5;
double np_tgauss2_a2 = -2.0*0.002351671671248367;

double np_tgauss2_k = 2.90113075268188e-01, np_tgauss2_k2 = 9.17819591566274e-01;
double np_tgauss2_k22 = 1.40866160472795e-01, np_tgauss2_km = 2.23983611906613e-01;

extern double cksup[OP_NCFUN][2];

double timing_extern  = -1.0;

void np_set_tgauss2(double * coefficients){
  np_tgauss2_b = coefficients[TG2_B];
  np_tgauss2_alpha = coefficients[TG2_ALPHA];
  np_tgauss2_c0 = coefficients[TG2_C0];

  np_tgauss2_a0 = coefficients[TG2_A0];
  np_tgauss2_a1 = coefficients[TG2_A1];
  np_tgauss2_a2 = coefficients[TG2_A2];

  np_tgauss2_k = coefficients[TG2_K];
  np_tgauss2_k2 = coefficients[TG2_K2];
  np_tgauss2_k22 = coefficients[TG2_K22];
  np_tgauss2_km = coefficients[TG2_KM];

  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_NORMAL]][0] = -np_tgauss2_b;
  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_NORMAL]][1] = np_tgauss2_b;

  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_CONVOLUTION]][0] = -2.0*np_tgauss2_b;
  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_CONVOLUTION]][1] = 2.0*np_tgauss2_b;

  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_DERIVATIVE]][0] = -np_tgauss2_b;
  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_DERIVATIVE]][1] = np_tgauss2_b;

  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_INTEGRAL]][0] = -np_tgauss2_b;
  cksup[CK_TGAUSS2 + OP_CFUN_OFFSETS[OP_INTEGRAL]][1] = DBL_MAX;
}

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
                   double * mysd, int * myopti, double * myoptd, double * myans, double * fval,
                   double * objective_function_values, double * objective_function_evals,
                   double * objective_function_invalid, double * timing,
                   int * penalty_mode, double * penalty_mult){
  /* Likelihood bandwidth selection for density estimation */

  double **matrix_y;

  double *vector_continuous_stddev;
  double *vsfh, *vector_scale_factor, *vector_scale_factor_multistart;

  double fret, fret_best;
  double ftol, tol;
  double (* bwmfunc)(double *) = NULL;

  double small, lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init; 
  double lbd_init, hbd_init, d_init;
  int dfc_dir;
  
  int i,j;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, iImproved;
  int itmax, iter;
  int int_use_starting_values;
  int scale_cat;

  int * ipt = NULL;  // point permutation, see tree.c
  int old_bw;


  num_reg_unordered_extern = myopti[BW_NUNOI];
  num_reg_ordered_extern = myopti[BW_NORDI];
  num_reg_continuous_extern = myopti[BW_NCONI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[BW_NOBSI];
  iMultistart = myopti[BW_IMULTII];
  iNum_Multistart = myopti[BW_NMULTII];

  KERNEL_den_extern = myopti[BW_CKRNEVI];
  KERNEL_den_unordered_extern = myopti[BW_UKRNEVI];
  KERNEL_den_ordered_extern = myopti[BW_OKRNEVI];

  int_use_starting_values= myopti[BW_USTARTI];
  int_LARGE_SF=myopti[BW_LSFI];
  BANDWIDTH_den_extern=myopti[BW_DENI];
  int_RESTART_FROM_MIN = myopti[BW_REMINI];
  int_MINIMIZE_IO = myopti[BW_MINIOI];

  itmax=myopti[BW_ITMAXI];
  old_bw=myopti[BW_OLDBW];
  int_TREE_X = myopti[BW_DOTREEI];
  scale_cat = myopti[BW_SCATI];
  bwm_use_transform = myopti[BW_TBNDI];
  if (BANDWIDTH_den_extern != BW_FIXED)
    bwm_use_transform = 0;
  if (bwm_use_transform) {
    int n = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
    if (bwm_transform_buf_len < n + 1) {
      bwm_transform_buf = (double *) realloc(bwm_transform_buf, (n + 1) * sizeof(double));
      bwm_transform_buf_len = n + 1;
    }
  }

  ftol=myoptd[BW_FTOLD];
  tol=myoptd[BW_TOLD];
  small=myoptd[BW_SMALLD];

  dfc_dir = myopti[BW_DFC_DIRI];

  lbc_dir = myoptd[BW_LBC_DIRD];
  c_dir = myoptd[BW_C_DIRD];
  initc_dir = myoptd[BW_INITC_DIRD]; 
  lbd_dir = myoptd[BW_LBD_DIRD]; 
  hbd_dir = myoptd[BW_HBD_DIRD]; 
  d_dir = myoptd[BW_D_DIRD]; 
  initd_dir = myoptd[BW_INITD_DIRD]; 

  lbc_init = myoptd[BW_LBC_INITD]; 
  hbc_init = myoptd[BW_HBC_INITD]; 
  c_init = myoptd[BW_C_INITD]; 

  lbd_init = myoptd[BW_LBD_INITD]; 
  hbd_init = myoptd[BW_HBD_INITD]; 
  d_init = myoptd[BW_D_INITD]; 

  nconfac_extern = myoptd[BW_NCONFD];
  ncatfac_extern = myoptd[BW_NCATFD];

/* Allocate memory for objects */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);


  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_y = alloc_matd(num_var + 1, num_var +1);
  vector_scale_factor = alloc_vecd(num_var + 1);
  vsfh = alloc_vecd(num_var + 1);

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

  ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  // attempt tree build, if enabled 
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                 4*num_reg_continuous_extern, ipt, &kdt_extern_X);
  

    //put training data into tree-order using the index array

    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_unordered_train_extern[j][i]=myuno[j*num_obs_train_extern+ipt[i]];
    
    
    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_ordered_train_extern[j][i]=myord[j*num_obs_train_extern+ipt[i]];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_continuous_train_extern[j][i]=mycon[j*num_obs_train_extern+ipt[i]];

  }


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

  for (j = 0; j < num_reg_continuous_extern; j++)
    vector_continuous_stddev[j] = mysd[j];

  vector_continuous_stddev_extern = vector_continuous_stddev;

  /* Initialize scale factors and Hessian for NR modules */

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0, 
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0, 
                                    KERNEL_den_unordered_extern,                                    
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0, 
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0, 
                                    KERNEL_den_unordered_extern,                                    
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vsfh,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           0,
                           0,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED, 
                           lbc_dir, dfc_dir, c_dir,initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);

  /* When multistarting, set counter */

  imsnum = iMs_counter = 0;
  imstot = iNum_Multistart;
  
  /* Conduct direction set search */

  /* assign the function to be optimized */
  if(old_bw){
    switch(myopti[BW_MI]){
    case BWM_CVML : bwmfunc = cv_func_density_categorical_ml; break;
    case BWM_CVLS : bwmfunc = cv_func_density_categorical_ls; break;
      //case BWM_CVML_NP : bwmfunc = cv_func_np_density_categorical_ml; break;
    default : REprintf("np.c: invalid bandwidth selection method.");
      error("np.c: invalid bandwidth selection method."); break;
    }
  } else {
    switch(myopti[BW_MI]){
    case BWM_CVML : bwmfunc = np_cv_func_density_categorical_ml; break;
    case BWM_CVLS : bwmfunc = np_cv_func_density_categorical_ls; break;
    default : REprintf("np.c: invalid bandwidth selection method.");
      error("np.c: invalid bandwidth selection method."); break;
    }

  }

  if (bwm_use_transform)
    bwm_to_unconstrained(vector_scale_factor, num_var);

  spinner(0);

  bwmfunc_raw = bwmfunc;
  bwm_num_reg_continuous = num_reg_continuous_extern;
  bwm_num_reg_unordered = num_reg_unordered_extern;
  bwm_num_reg_ordered = num_reg_ordered_extern;
  bwm_kernel_unordered = KERNEL_den_unordered_extern;
  bwm_kernel_unordered_vec = NULL;
  bwm_kernel_unordered_len = 0;
  bwm_num_categories = num_categories_extern;
  bwm_reset_counters();
  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  if (penalty_mode[0] == 1) {
    double pmult = penalty_mult[0];
    double baseline;
    if (pmult < 1.0) pmult = 1.0;
    baseline = bwmfunc_raw(vector_scale_factor);
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      double *tmp = alloc_vecd(num_var + 1);
      for (i = 1; i <= num_var; i++)
        tmp[i] = vector_scale_factor[i];
      for (i = 1; i <= num_reg_continuous_extern; i++)
        tmp[i] *= 2.0;
      for (i = 0; i < num_reg_unordered_extern; i++) {
        int idx = num_reg_continuous_extern + 1 + i;
        double maxbw = max_unordered_bw(num_categories_extern[i], KERNEL_den_unordered_extern);
        tmp[idx] = 0.5*maxbw;
      }
      for (i = 0; i < num_reg_ordered_extern; i++) {
        int idx = num_reg_continuous_extern + num_reg_unordered_extern + 1 + i;
        tmp[idx] = 0.5;
      }
      baseline = bwmfunc_raw(tmp);
      safe_free(tmp);
    }
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      bwm_penalty_value = pmult * 1.0e6;
    } else {
      bwm_penalty_value = baseline + (fabs(baseline) + 1.0) * pmult;
    }
    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }

  fret_best = bwmfunc_wrapper(vector_scale_factor);
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
         bwmfunc_wrapper);

  /* int_RESTART_FROM_MIN needs to be set */

  if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

    initialize_nr_directions(BANDWIDTH_den_extern,
                             num_obs_train_extern,
                             num_reg_continuous_extern,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             0,
                             0,
                             0,
                             vsfh,
                             num_categories_extern,
                             matrix_y,
                             0, int_RANDOM_SEED, 
                             lbc_dir, dfc_dir, c_dir, initc_dir,
                             lbd_dir, hbd_dir, d_dir, initd_dir,
                             matrix_X_continuous_train_extern,
                             matrix_Y_continuous_train_extern);



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
           bwmfunc_wrapper);
  }

  iImproved = (fret < fret_best);
  *timing = timing_extern;

  /* When multistarting save initial minimum of objective function and scale factors */
  objective_function_values[0]=-fret;
  objective_function_evals[0]=bwm_eval_count;
  objective_function_invalid[0]=bwm_invalid_count;

  if(iMultistart == IMULTI_TRUE){
    fret_best = fret;
    vector_scale_factor_multistart = alloc_vecd(num_var + 1);
    for(i = 1; i <= num_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
    		
    /* Conduct search from new random values of the search parameters */
       	
    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and directions for NR modules */

      initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                        1,                /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        0, 
                                        0,
                                        0,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        0, 
                                        KERNEL_den_unordered_extern,                                    
                                        0,
                                        scale_cat,
                                        pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                        nconfac_extern, ncatfac_extern,
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor,
                                        lbc_init, hbc_init, c_init, 
                                        lbd_init, hbd_init, d_init,
                                        matrix_X_continuous_train_extern,
                                        matrix_Y_continuous_train_extern);


      initialize_nr_directions(BANDWIDTH_den_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               0,
                               0,
                               0,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               1, int_RANDOM_SEED, 
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);



      /* Conduct direction set search */

      if (bwm_use_transform)
        bwm_to_unconstrained(vector_scale_factor, num_var);

      bwm_reset_counters();

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
             bwmfunc_wrapper);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE){
        initialize_nr_directions(BANDWIDTH_den_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 num_reg_unordered_extern,
                                 num_reg_ordered_extern,
                                 0,
                                 0,
                                 0,
                                 vsfh,
                                 num_categories_extern,
                                 matrix_y,
                                 0, int_RANDOM_SEED, 
                                 lbc_dir, dfc_dir, c_dir, initc_dir,
                                 lbd_dir, hbd_dir, d_dir, initd_dir,
                                 matrix_X_continuous_train_extern,
                                 matrix_Y_continuous_train_extern);



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
               bwmfunc_wrapper);
      }
       			
      /* If this run resulted in an improved minimum save information */
       		
      if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        *timing = timing_extern;
       
        for(i = 1; i <= num_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }
      objective_function_values[iMs_counter]=-fret;
      objective_function_evals[iMs_counter]=bwm_eval_count;
      objective_function_invalid[iMs_counter]=bwm_invalid_count;
    }

    /* Save best for estimation */

    fret = fret_best;

    for(i = 1; i <= num_var; i++)
      vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];

    free(vector_scale_factor_multistart);

  }

  if (bwm_use_transform)
    bwm_to_constrained(vector_scale_factor, num_var);

  /* return data to R */
  if (BANDWIDTH_den_extern == BW_GEN_NN || 
      BANDWIDTH_den_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern; i++ )
      vector_scale_factor[i+1]=np_fround(vector_scale_factor[i+1]);
  }

  for( i=0; i<num_var; i++ )
    myans[i]=vector_scale_factor[i+1];

  fval[0] = -fret;
  fval[1] = iImproved;

  /* end return data */

  /* Free data objects */

  free_mat(matrix_X_unordered_train_extern, num_reg_unordered_extern);
  free_mat(matrix_X_ordered_train_extern, num_reg_ordered_extern);
  free_mat(matrix_X_continuous_train_extern, num_reg_continuous_extern);
  free_mat(matrix_y, num_var + 1);
  free(vector_scale_factor);
  free(vsfh);
  free(num_categories_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  free(vector_continuous_stddev);

  free(ipt);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  if(int_MINIMIZE_IO != IO_MIN_TRUE)
    Rprintf("\r                   \r");

  return ;
  
}

 
// For distributions the bandwidth selection involves evaluating the CDF at a number of points.
// We allow one to specify those points, passed in by mye{uno,ord,con}.

void np_distribution_bw(double * myuno, double * myord, double * mycon, 
                        double * myeuno, double * myeord, double * myecon, double * mysd,
                        int * myopti, double * myoptd, double * myans, double * fval,
                        double * objective_function_values, double * objective_function_evals,
                        double * objective_function_invalid, double * timing,
                        int * penalty_mode, double * penalty_mult){
  /* Likelihood bandwidth selection for density estimation */

  double **matrix_y;

  double *vector_continuous_stddev;
  double *vsfh, *vector_scale_factor, *vector_scale_factor_multistart;

  double fret, fret_best;
  double ftol, tol;
  double (* bwmfunc)(double *) = NULL;

  double small, lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init; 
  double lbd_init, hbd_init, d_init;
  int dfc_dir;
  
  int i,j;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, iImproved;
  int itmax, iter;
  int int_use_starting_values, cdfontrain;

  int scale_cat;

  int * ipt = NULL, * ipe = NULL;

  cdfontrain_extern = cdfontrain =  myopti[DBW_CDFONTRAIN];

  num_reg_unordered_extern = myopti[DBW_NUNOI];
  num_reg_ordered_extern = myopti[DBW_NORDI];
  num_reg_continuous_extern = myopti[DBW_NCONI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[DBW_NOBSI];
  
  num_obs_eval_extern = cdfontrain ? num_obs_train_extern : myopti[DBW_NEVALI];

  iMultistart = myopti[DBW_IMULTII];
  iNum_Multistart = myopti[DBW_NMULTII];

  KERNEL_den_extern = myopti[DBW_CKRNEVI];
  KERNEL_den_unordered_extern = myopti[DBW_UKRNEVI];
  KERNEL_den_ordered_extern = myopti[DBW_OKRNEVI];

  int_use_starting_values= myopti[DBW_USTARTI];
  int_LARGE_SF=myopti[DBW_LSFI];
  BANDWIDTH_den_extern=myopti[DBW_DENI];
  int_RESTART_FROM_MIN = myopti[DBW_REMINI];
  int_MINIMIZE_IO = myopti[DBW_MINIOI];

  itmax=myopti[DBW_ITMAXI];

  int_TREE_X = myopti[DBW_DOTREEI];
  scale_cat = myopti[DBW_SCATI];
  bwm_use_transform = myopti[DBW_TBNDI];
  if (BANDWIDTH_den_extern != BW_FIXED)
    bwm_use_transform = 0;
  if (bwm_use_transform) {
    int n = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
    if (bwm_transform_buf_len < n + 1) {
      bwm_transform_buf = (double *) realloc(bwm_transform_buf, (n + 1) * sizeof(double));
      bwm_transform_buf_len = n + 1;
    }
  }

  ftol=myoptd[DBW_FTOLD];
  tol=myoptd[DBW_TOLD];
  small=myoptd[DBW_SMALLD];

  dfc_dir = myopti[DBW_DFC_DIRI];
  lbc_dir = myoptd[DBW_LBC_DIRD];
  c_dir = myoptd[DBW_C_DIRD];
  initc_dir = myoptd[DBW_INITC_DIRD]; 

  lbd_dir = myoptd[DBW_LBD_DIRD]; 
  hbd_dir = myoptd[DBW_HBD_DIRD]; 
  d_dir = myoptd[DBW_D_DIRD]; 
  initd_dir = myoptd[DBW_INITD_DIRD]; 

  lbc_init = myoptd[DBW_LBC_INITD]; 
  hbc_init = myoptd[DBW_HBC_INITD]; 
  c_init = myoptd[DBW_C_INITD]; 

  lbd_init = myoptd[DBW_LBD_INITD]; 
  hbd_init = myoptd[DBW_HBD_INITD]; 
  d_init = myoptd[DBW_D_INITD]; 

  nconfac_extern = myoptd[DBW_NCONFD];
  ncatfac_extern = myoptd[DBW_NCATFD];

  dbl_memfac_dls_extern = myoptd[DBW_MEMORYD];

/* Allocate memory for objects */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  if(cdfontrain){
    matrix_X_unordered_eval_extern = matrix_X_unordered_train_extern;
    matrix_X_ordered_eval_extern = matrix_X_ordered_train_extern;
    matrix_X_continuous_eval_extern = matrix_X_continuous_train_extern;
  } else {
    matrix_X_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_unordered_extern);
    matrix_X_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_ordered_extern);
    matrix_X_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_reg_continuous_extern);
  }

  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_y = alloc_matd(num_var + 1, num_var +1);
  vector_scale_factor = alloc_vecd(num_var + 1);
  vsfh = alloc_vecd(num_var + 1);
  // nb check vals
  matrix_categorical_vals_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  if(num_reg_unordered_extern > 0)
    error("np.c: distribution bw selection only works on ordered and continuous data."); 

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


  // points for evaluating the cdf
  if(!cdfontrain){
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_unordered_eval_extern[j][i]=myeuno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_ordered_eval_extern[j][i]=myeord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_X_continuous_eval_extern[j][i]=myecon[j*num_obs_eval_extern+i];
  }

  ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  if(!cdfontrain) {
    ipe = (int *)malloc(num_obs_eval_extern*sizeof(int));
    if(!(ipe != NULL))
      error("!(ipe != NULL)");

    for(i = 0; i < num_obs_eval_extern; i++){
      ipe[i] = i;
    }
  } else {
    ipe = ipt;
  }

  // attempt tree build, if enabled 
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    if((BANDWIDTH_den_extern != BW_ADAP_NN) || ((BANDWIDTH_den_extern == BW_ADAP_NN) && cdfontrain)){
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipt, &kdt_extern_X);

      //put training data into tree-order using the index array

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_unordered_train_extern[j][i]=myuno[j*num_obs_train_extern+ipt[i]];  
    
      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_ordered_train_extern[j][i]=myord[j*num_obs_train_extern+ipt[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_continuous_train_extern[j][i]=mycon[j*num_obs_train_extern+ipt[i]];

    } else {
      build_kdtree(matrix_X_continuous_eval_extern, num_obs_eval_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipe, &kdt_extern_X);      

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_unordered_eval_extern[j][i]=myeuno[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_ordered_eval_extern[j][i]=myeord[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_continuous_eval_extern[j][i]=myecon[j*num_obs_eval_extern+ipe[i]];

    }

  }

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

  for (j = 0; j < num_reg_continuous_extern; j++)
    vector_continuous_stddev[j] = mysd[j];

  vector_continuous_stddev_extern = vector_continuous_stddev;


  /* Initialize scale factors and Directions for NR modules */

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0,
                                    KERNEL_den_unordered_extern,
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_eval_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0,
                                    KERNEL_den_unordered_extern,
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vsfh,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_eval_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           0,
                           0,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED, 
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);


  /* When multistarting, set counter */

  imsnum = iMs_counter = 0;
  imstot = iNum_Multistart;
  
  /* Conduct direction set search */

  /* assign the function to be optimized */
  switch(myopti[DBW_MI]){
  case DBWM_CVLS : bwmfunc = cv_func_distribution_categorical_ls; break;
  default : REprintf("np.c: invalid bandwidth selection method.");
    error("np.c: invalid bandwidth selection method."); break;
  }

  if (bwm_use_transform)
    bwm_to_unconstrained(vector_scale_factor, num_var);

  spinner(0);

  bwmfunc_raw = bwmfunc;
  bwm_num_reg_continuous = num_reg_continuous_extern;
  bwm_num_reg_unordered = num_reg_unordered_extern;
  bwm_num_reg_ordered = num_reg_ordered_extern;
  bwm_kernel_unordered = KERNEL_den_unordered_extern;
  bwm_kernel_unordered_vec = NULL;
  bwm_kernel_unordered_len = 0;
  bwm_num_categories = num_categories_extern;
  bwm_reset_counters();
  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  if (penalty_mode[0] == 1) {
    double pmult = penalty_mult[0];
    double baseline;
    if (pmult < 1.0) pmult = 1.0;
    baseline = bwmfunc_raw(vector_scale_factor);
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      double *tmp = alloc_vecd(num_var + 1);
      for (i = 1; i <= num_var; i++)
        tmp[i] = vector_scale_factor[i];
      for (i = 1; i <= num_reg_continuous_extern; i++)
        tmp[i] *= 2.0;
      for (i = 0; i < num_reg_unordered_extern; i++) {
        int idx = num_reg_continuous_extern + 1 + i;
        double maxbw = max_unordered_bw(num_categories_extern[i], KERNEL_den_unordered_extern);
        tmp[idx] = 0.5*maxbw;
      }
      for (i = 0; i < num_reg_ordered_extern; i++) {
        int idx = num_reg_continuous_extern + num_reg_unordered_extern + 1 + i;
        tmp[idx] = 0.5;
      }
      baseline = bwmfunc_raw(tmp);
      safe_free(tmp);
    }
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      bwm_penalty_value = pmult * 1.0e6;
    } else {
      bwm_penalty_value = baseline + (fabs(baseline) + 1.0) * pmult;
    }
    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }

  fret_best = bwmfunc_wrapper(vector_scale_factor);
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
         bwmfunc_wrapper);

  /* int_RESTART_FROM_MIN needs to be set */

  if(int_RESTART_FROM_MIN == RE_MIN_TRUE){
    initialize_nr_directions(BANDWIDTH_den_extern,
                             num_obs_train_extern,
                             num_reg_continuous_extern,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             0,
                             0,
                             0,
                             vsfh,
                             num_categories_extern,
                             matrix_y,
                             0, int_RANDOM_SEED,  
                             lbc_dir, dfc_dir, c_dir, initc_dir,
                             lbd_dir, hbd_dir, d_dir, initd_dir,
                             matrix_X_continuous_train_extern,
                             matrix_Y_continuous_train_extern);



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
           bwmfunc_wrapper);
  }

  iImproved = (fret < fret_best);
  *timing = timing_extern;

  objective_function_values[0]=fret;
  objective_function_evals[0]=bwm_eval_count;
  objective_function_invalid[0]=bwm_invalid_count;
  /* When multistarting save initial minimum of objective function and scale factors */

  if(iMultistart == IMULTI_TRUE){
    fret_best = fret;
    vector_scale_factor_multistart = alloc_vecd(num_var + 1);
    for(i = 1; i <= num_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
    		
    /* Conduct search from new random values of the search parameters */
       	
    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and directions for NR modules */

      initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                        1,        /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        0,
                                        0,
                                        0,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        0,
                                        KERNEL_den_unordered_extern,
                                        0,
                                        scale_cat,
                                        pow((double)4.0/(double)3.0,0.2),     /* Init for continuous vars */
                                        nconfac_extern, ncatfac_extern,
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor,
                                        lbc_init, hbc_init, c_init, 
                                        lbd_init, hbd_init, d_init,
                                        matrix_X_continuous_eval_extern,
                                        matrix_Y_continuous_train_extern);

      initialize_nr_directions(BANDWIDTH_den_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               0,
                               0,
                               0,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               1, int_RANDOM_SEED,  
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);

      /* Conduct direction set search */

      if (bwm_use_transform)
        bwm_to_unconstrained(vector_scale_factor, num_var);

      bwm_reset_counters();

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
             bwmfunc_wrapper);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

        initialize_nr_directions(BANDWIDTH_den_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 num_reg_unordered_extern,
                                 num_reg_ordered_extern,
                                 0,
                                 0,
                                 0,
                                 vsfh,
                                 num_categories_extern,
                                 matrix_y,
                                 0, int_RANDOM_SEED,  
                                 lbc_dir, dfc_dir, c_dir, initc_dir,
                                 lbd_dir, hbd_dir, d_dir, initd_dir,
                                 matrix_X_continuous_train_extern,
                                 matrix_Y_continuous_train_extern);

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
               bwmfunc_wrapper);
      }
       			
      /* If this run resulted in an improved minimum save information */
       		
      if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        *timing = timing_extern;
       
        for(i = 1; i <= num_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }
      objective_function_values[iMs_counter]=fret;
      objective_function_evals[iMs_counter]=bwm_eval_count;
      objective_function_invalid[iMs_counter]=bwm_invalid_count;
    }

    /* Save best for estimation */

    fret = fret_best;

    for(i = 1; i <= num_var; i++)
      vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];

    free(vector_scale_factor_multistart);

  }

  if (bwm_use_transform)
    bwm_to_constrained(vector_scale_factor, num_var);

  /* return data to R */
  if (BANDWIDTH_den_extern == BW_GEN_NN || 
      BANDWIDTH_den_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern; i++ )
      vector_scale_factor[i+1]=np_fround(vector_scale_factor[i+1]);
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

  if(!cdfontrain){
    free_mat(matrix_X_unordered_eval_extern, num_reg_unordered_extern);
    free_mat(matrix_X_ordered_eval_extern, num_reg_ordered_extern);
    free_mat(matrix_X_continuous_eval_extern, num_reg_continuous_extern);
  }

  free_mat(matrix_y, num_var + 1);
  free(vector_scale_factor);
  free(vsfh);
  free(num_categories_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  free(vector_continuous_stddev);

  free(ipt);
  if(!cdfontrain)
    free(ipe);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  if(int_MINIMIZE_IO != IO_MIN_TRUE)
    Rprintf("\r                   \r");

  return ;
  
}


void np_density_conditional_bw(double * c_uno, double * c_ord, double * c_con, 
                               double * u_uno, double * u_ord, double * u_con,
                               double * mysd,
                               int * myopti, double * myoptd, double * myans, double * fval,
                               double * objective_function_values, double * objective_function_evals,
                               double * objective_function_invalid, double * timing,
                               int * penalty_mode, double * penalty_mult){
/* Likelihood bandwidth selection for density estimation */

  double **matrix_y = NULL;

  double *vector_continuous_stddev;
  double *vsfh, *vector_scale_factor, *vector_scale_factor_multistart;

  double fret, fret_best;
  double ftol, tol;
  double (* bwmfunc)(double *) = NULL;

  double small, lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init; 
  double lbd_init, hbd_init, d_init;
  int dfc_dir;
  
  int i,j;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, num_all_var, num_var_var, iImproved;
  int itmax, iter;
  int int_use_starting_values, ibwmfunc, old_cdens, scale_cat;

  int num_all_cvar, num_all_uvar, num_all_ovar;

  int * ipt_X = NULL, * ipt_XY = NULL, * ipt_Y = NULL; 
  int * ipt_lookup_XY = NULL, * ipt_lookup_Y = NULL, * ipt_lookup_X = NULL;

  /* Ensure optional Y-only categorical arrays are reset each call */
  num_categories_extern_Y = NULL;
  matrix_categorical_vals_extern_Y = NULL;

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
  old_cdens = myopti[CBW_OLDI];
  int_TREE_XY = int_TREE_Y = int_TREE_X = myopti[CBW_TREEI];
  scale_cat = myopti[CBW_SCATI];
  
  ibwmfunc = myopti[CBW_MI];
  bwm_use_transform = myopti[CBW_TBNDI];
  if (BANDWIDTH_den_extern != BW_FIXED)
    bwm_use_transform = 0;
  if (bwm_use_transform) {
    int n = num_var_continuous_extern + num_reg_continuous_extern +
      num_var_unordered_extern + num_reg_unordered_extern +
      num_var_ordered_extern + num_reg_ordered_extern;
    if (bwm_transform_buf_len < n + 1) {
      bwm_transform_buf = (double *) realloc(bwm_transform_buf, (n + 1) * sizeof(double));
      bwm_transform_buf_len = n + 1;
    }
  }

  ftol=myoptd[CBW_FTOLD];
  tol=myoptd[CBW_TOLD];
  small=myoptd[CBW_SMALLD];
  dbl_memfac_ccdf_extern = myoptd[CBW_MEMFACD];

  dfc_dir = myopti[CBW_DFC_DIRI];
  lbc_dir = myoptd[CBW_LBC_DIRD];
  c_dir = myoptd[CBW_C_DIRD];
  initc_dir = myoptd[CBW_INITC_DIRD]; 

  lbd_dir = myoptd[CBW_LBD_DIRD]; 
  hbd_dir = myoptd[CBW_HBD_DIRD]; 
  d_dir = myoptd[CBW_D_DIRD]; 
  initd_dir = myoptd[CBW_INITD_DIRD]; 

  lbc_init = myoptd[CBW_LBC_INITD]; 
  hbc_init = myoptd[CBW_HBC_INITD]; 
  c_init = myoptd[CBW_C_INITD]; 

  lbd_init = myoptd[CBW_LBD_INITD]; 
  hbd_init = myoptd[CBW_HBD_INITD]; 
  d_init = myoptd[CBW_D_INITD]; 

  nconfac_extern = myoptd[CBW_NCONFD];
  ncatfac_extern = myoptd[CBW_NCATFD];

/* Allocate memory for objects */

  //if((BANDWIDTH_den_extern != BW_FIXED) && (ibwmfunc == CBWM_CVLS))
  //old_cdens = 1;

  matrix_Y_unordered_train_extern = alloc_matd(num_obs_train_extern, num_var_unordered_extern);
  matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);
  matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);
	
  num_categories_extern = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                     num_reg_unordered_extern + num_reg_ordered_extern);

  num_categories_extern_XY = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                     num_reg_unordered_extern + num_reg_ordered_extern);

  num_categories_extern_X = alloc_vecu(num_reg_unordered_extern + num_reg_ordered_extern);
  
  if(ibwmfunc == CBWM_CVLS){
    num_categories_extern_Y = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern);
  }

  matrix_y = alloc_matd(num_all_var + 1, num_all_var + 1);
  vector_scale_factor = alloc_vecd(num_all_var + 1);
  vsfh = alloc_vecd(num_all_var + 1);
  
  matrix_categorical_vals_extern = 
    alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern + 
               num_reg_unordered_extern + num_reg_ordered_extern);

  matrix_categorical_vals_extern_X = 
    alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  if(ibwmfunc == CBWM_CVLS){
    matrix_categorical_vals_extern_Y = 
      alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern);
  }

  matrix_categorical_vals_extern_XY = 
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

  /* data has been parsed, make the joint xy dataset */

  // the main tree is x, the alt tree is xy

  num_all_cvar = num_reg_continuous_extern + num_var_continuous_extern;
  num_all_uvar = num_reg_unordered_extern + num_var_unordered_extern;
  num_all_ovar = num_reg_ordered_extern + num_var_ordered_extern;

  // we use 3 trees for cpdf ls, and 2 for cpdf ml

  ipt_X = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_X != NULL))
    error("!(ipt_X != NULL)");

  ipt_lookup_X = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_lookup_X != NULL))
    error("!(ipt_lookup_X != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt_lookup_X[i] = ipt_X[i] = i;
  }

  ipt_extern_X = ipt_X;
  ipt_lookup_extern_X = ipt_lookup_X;

  if(ibwmfunc == CBWM_CVLS){
    ipt_Y = (int *)malloc(num_obs_train_extern*sizeof(int));
    if(!(ipt_Y != NULL))
      error("!(ipt_Y != NULL)");

    ipt_lookup_Y = (int *)malloc(num_obs_train_extern*sizeof(int));
    if(!(ipt_lookup_Y != NULL))
      error("!(ipt_lookup_Y != NULL)");

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_Y[i] = ipt_Y[i] = i;
    }

  } else {
    ipt_Y = ipt_X;
    ipt_lookup_Y = ipt_lookup_X;
  }

  ipt_extern_Y = ipt_Y;
  ipt_lookup_extern_Y = ipt_lookup_Y;

  ipt_XY = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_XY != NULL))
    error("!(ipt_XY != NULL)");

  ipt_lookup_XY = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_lookup_XY != NULL))
    error("!(ipt_lookup_XY != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt_lookup_XY[i] = ipt_XY[i] = i;
  }

  ipt_extern_XY = ipt_XY;
  ipt_lookup_extern_XY = ipt_lookup_XY;

  int_TREE_XY = int_TREE_XY && (((num_all_cvar) != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);

  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);

  int_TREE_Y = int_TREE_Y && (ibwmfunc == CBWM_CVLS) && ((num_var_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);


  if(int_TREE_X == NP_TREE_TRUE){
    build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                 4*num_reg_continuous_extern, ipt_X, &kdt_extern_X);
  

    //put training data into tree-order using the index array

    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+ipt_X[i]];
    
    
    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+ipt_X[i]];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+ipt_X[i]];

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_X[ipt_X[i]] = i;
    }
  }

  if(int_TREE_Y == NP_TREE_TRUE){
    build_kdtree(matrix_Y_continuous_train_extern, num_obs_train_extern, num_var_continuous_extern, 
                 4*num_var_continuous_extern, ipt_Y, &kdt_extern_Y);

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_Y[ipt_Y[i]] = i;
    }

  }


  if((int_TREE_X == NP_TREE_TRUE) || (int_TREE_Y == NP_TREE_TRUE)){
    for(j=0;j<num_var_unordered_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_unordered_train_extern[j][i]=c_uno[j*num_obs_train_extern+ipt_Y[i]];

    for(j=0;j<num_var_ordered_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_ordered_train_extern[j][i]=c_ord[j*num_obs_train_extern+ipt_Y[i]];

    for(j=0;j<num_var_continuous_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_continuous_train_extern[j][i]=c_con[j*num_obs_train_extern+ipt_Y[i]];

  }


  matrix_XY_continuous_train_extern = alloc_matd(num_obs_train_extern, num_all_cvar);
  matrix_XY_unordered_train_extern = alloc_matd(num_obs_train_extern, num_all_uvar);
  matrix_XY_ordered_train_extern = alloc_matd(num_obs_train_extern, num_all_ovar);


  for(j = 0; j < num_reg_unordered_extern; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+i];

  for(j = num_reg_unordered_extern; j < num_all_uvar; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_unordered_train_extern[j][i]=c_uno[(j-num_reg_unordered_extern)*num_obs_train_extern+i];


  for(j = 0; j < num_reg_ordered_extern; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+i];

  for(j = num_reg_ordered_extern; j < num_all_ovar; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_ordered_train_extern[j][i]=c_ord[(j-num_reg_ordered_extern)*num_obs_train_extern+i];


  for(j = 0; j < num_reg_continuous_extern; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+i];

  for(j = num_reg_continuous_extern; j < num_all_cvar; j++)
    for(i = 0; i < num_obs_train_extern; i++)
      matrix_XY_continuous_train_extern[j][i]=c_con[(j-num_reg_continuous_extern)*num_obs_train_extern+i];

  if(int_TREE_XY == NP_TREE_TRUE){

    build_kdtree(matrix_XY_continuous_train_extern, num_obs_train_extern, num_all_cvar, 
                 4*num_all_cvar, ipt_XY, &kdt_extern_XY);

    // put data into tree-order
    for(j = 0; j < num_reg_unordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_unordered_extern; j < num_all_uvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=c_uno[(j-num_reg_unordered_extern)*num_obs_train_extern+ipt_XY[i]];


    for(j = 0; j < num_reg_ordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_ordered_extern; j < num_all_ovar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=c_ord[(j-num_reg_ordered_extern)*num_obs_train_extern+ipt_XY[i]];


    for(j = 0; j < num_reg_continuous_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_continuous_extern; j < num_all_cvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=c_con[(j-num_reg_continuous_extern)*num_obs_train_extern+ipt_XY[i]];

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_XY[ipt_XY[i]] = i;
    }

  }

  determine_categorical_vals(num_obs_train_extern,
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

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern, num_var_ordered_extern, num_var_continuous_extern,
                        num_reg_unordered_extern, num_reg_ordered_extern, num_reg_continuous_extern,
                        vector_scale_factor+1,
                        num_categories_extern,
                        matrix_categorical_vals_extern,
                        NULL, NULL, NULL,
                        num_categories_extern_X, (ibwmfunc == CBWM_CVLS) ? num_categories_extern_Y : NULL, num_categories_extern_XY,
                        matrix_categorical_vals_extern_X, (ibwmfunc == CBWM_CVLS) ? matrix_categorical_vals_extern_Y : NULL, matrix_categorical_vals_extern_XY);
  

  vector_continuous_stddev = alloc_vecd(num_var_continuous_extern + num_reg_continuous_extern);

  for(j = 0; j < (num_var_continuous_extern + num_reg_continuous_extern); j++)
    vector_continuous_stddev[j] = mysd[j];

  vector_continuous_stddev_extern = vector_continuous_stddev;

  /* Initialize scale factors and Directions for NR modules */

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    num_var_continuous_extern,
                                    num_var_unordered_extern,
                                    num_var_ordered_extern,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    KERNEL_den_unordered_extern,
                                    KERNEL_reg_unordered_extern,
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    num_var_continuous_extern,
                                    num_var_unordered_extern,
                                    num_var_ordered_extern,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    KERNEL_den_unordered_extern,
                                    KERNEL_reg_unordered_extern,
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vsfh,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           num_var_continuous_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED,  
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);

  /* When multistarting, set counter */

  imsnum = iMs_counter = 0;
  imstot = iNum_Multistart;

  /* Conduct direction set search */

  /* assign the function to be optimized */

  /* 7/2/2010 */  
  if(old_cdens){
    switch(ibwmfunc){
    case CBWM_CVML : bwmfunc = cv_func_con_density_categorical_ml; break;
    case CBWM_CVLS : bwmfunc = cv_func_con_density_categorical_ls; break;
    case CBWM_NPLS : bwmfunc = np_cv_func_con_density_categorical_ls;break;
    case CBWM_CCDF : bwmfunc = cv_func_con_distribution_categorical_ccdf; break;
    default : REprintf("np.c: invalid bandwidth selection method.");
      error("np.c: invalid bandwidth selection method."); break;
    }
  } else {
    switch(ibwmfunc){
    case CBWM_CVML : bwmfunc = np_cv_func_con_density_categorical_ml; break;
    case CBWM_CVLS : bwmfunc = np_cv_func_con_density_categorical_ls_npksum; break;
    default : REprintf("np.c: invalid bandwidth selection method.");
      error("np.c: invalid bandwidth selection method."); break;
    }
  }

  if (bwm_use_transform)
    bwm_to_unconstrained(vector_scale_factor, num_all_var);

  spinner(0);

  bwmfunc_raw = bwmfunc;
  bwm_num_reg_continuous = num_var_continuous_extern + num_reg_continuous_extern;
  bwm_num_reg_unordered = num_var_unordered_extern + num_reg_unordered_extern;
  bwm_num_reg_ordered = num_var_ordered_extern + num_reg_ordered_extern;
  bwm_kernel_unordered = KERNEL_den_unordered_extern;
  bwm_kernel_unordered_len = bwm_num_reg_unordered;
  if (bwm_kernel_unordered_len > 0) {
    bwm_kernel_unordered_vec = alloc_vecu(bwm_kernel_unordered_len);
    for (i = 0; i < num_var_unordered_extern; i++)
      bwm_kernel_unordered_vec[i] = KERNEL_den_unordered_extern;
    for (i = 0; i < num_reg_unordered_extern; i++)
      bwm_kernel_unordered_vec[num_var_unordered_extern + i] = KERNEL_reg_unordered_extern;
  } else {
    bwm_kernel_unordered_vec = NULL;
  }
  bwm_num_categories = num_categories_extern;
  bwm_reset_counters();
  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  if (penalty_mode[0] == 1) {
    double pmult = penalty_mult[0];
    double baseline;
    if (pmult < 1.0) pmult = 1.0;
    baseline = bwmfunc_raw(vector_scale_factor);
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      double *tmp = alloc_vecd(num_all_var + 1);
      for (i = 1; i <= num_all_var; i++)
        tmp[i] = vector_scale_factor[i];
      for (i = 1; i <= (num_var_continuous_extern + num_reg_continuous_extern); i++)
        tmp[i] *= 2.0;
      for (i = 0; i < (num_var_unordered_extern + num_reg_unordered_extern); i++) {
        int idx = num_var_continuous_extern + num_reg_continuous_extern + 1 + i;
        double maxbw = max_unordered_bw(num_categories_extern[i], KERNEL_den_unordered_extern);
        tmp[idx] = 0.5*maxbw;
      }
      for (i = 0; i < (num_var_ordered_extern + num_reg_ordered_extern); i++) {
        int idx = num_var_continuous_extern + num_reg_continuous_extern +
          num_var_unordered_extern + num_reg_unordered_extern + 1 + i;
        tmp[idx] = 0.5;
      }
      baseline = bwmfunc_raw(tmp);
      safe_free(tmp);
    }
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      bwm_penalty_value = pmult * 1.0e6;
    } else {
      bwm_penalty_value = baseline + (fabs(baseline) + 1.0) * pmult;
    }
    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }

  fret_best = bwmfunc_wrapper(vector_scale_factor);
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
         bwmfunc_wrapper);

  if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

    initialize_nr_directions(BANDWIDTH_den_extern,
                             num_obs_train_extern,
                             num_reg_continuous_extern,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             num_var_continuous_extern,
                             num_var_unordered_extern,
                             num_var_ordered_extern,
                             vsfh,
                             num_categories_extern,
                             matrix_y,
                             0, int_RANDOM_SEED,  
                             lbc_dir, dfc_dir, c_dir, initc_dir,
                             lbd_dir, hbd_dir, d_dir, initd_dir,
                             matrix_X_continuous_train_extern,
                             matrix_Y_continuous_train_extern);


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
           bwmfunc_wrapper);

  }

  iImproved = (fret < fret_best);
  *timing = timing_extern;

  objective_function_values[0]=-fret;
  objective_function_evals[0]=bwm_eval_count;
  objective_function_invalid[0]=bwm_invalid_count;
  /* When multistarting save initial minimum of objective function and scale factors */


  if(iMultistart == IMULTI_TRUE){
    fret_best = fret;
    vector_scale_factor_multistart = alloc_vecd(num_all_var + 1);
    for(i = 1; i <= num_all_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
			

    /* Conduct search from new random values of the search parameters */
		
    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and directions for NR modules */
      initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                        1,                /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        num_var_continuous_extern,
                                        num_var_unordered_extern,
                                        num_var_ordered_extern,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        KERNEL_den_unordered_extern,
                                        KERNEL_reg_unordered_extern,
                                        0,
                                        scale_cat,
                                        pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                        nconfac_extern, ncatfac_extern,
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor,
                                        lbc_init, hbc_init, c_init, 
                                        lbd_init, hbd_init, d_init,
                                        matrix_X_continuous_train_extern,
                                        matrix_Y_continuous_train_extern);

      initialize_nr_directions(BANDWIDTH_den_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               num_var_continuous_extern,
                               num_var_unordered_extern,
                               num_var_ordered_extern,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               1, int_RANDOM_SEED,  
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);

      /* Conduct direction set search */

      if (bwm_use_transform)
        bwm_to_unconstrained(vector_scale_factor, num_all_var);

      bwm_reset_counters();
      
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
             bwmfunc_wrapper);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

        initialize_nr_directions(BANDWIDTH_den_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 num_reg_unordered_extern,
                                 num_reg_ordered_extern,
                                 num_var_continuous_extern,
                                 num_var_unordered_extern,
                                 num_var_ordered_extern,
                                 vsfh,
                                 num_categories_extern,
                                 matrix_y,
                                 0, int_RANDOM_SEED,  
                                 lbc_dir, dfc_dir, c_dir, initc_dir,
                                 lbd_dir, hbd_dir, d_dir, initd_dir,
                                 matrix_X_continuous_train_extern,
                                 matrix_Y_continuous_train_extern);

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
               bwmfunc_wrapper);
      }
				
      /* If this run resulted in an improved minimum save information */
      
      if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        *timing = timing_extern;
        
        for(i = 1; i <= num_all_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }
      objective_function_values[iMs_counter]=-fret;
      objective_function_evals[iMs_counter]=bwm_eval_count;
      objective_function_invalid[iMs_counter]=bwm_invalid_count;
    }

    /* Save best for estimation */

    fret = fret_best;
    for(i = 1; i <= num_all_var; i++)
      vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];
    free(vector_scale_factor_multistart);
  }

  if (bwm_use_transform)
    bwm_to_constrained(vector_scale_factor, num_all_var);

  /* return data to R */
  if (BANDWIDTH_den_extern == BW_GEN_NN || 
      BANDWIDTH_den_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern+num_var_continuous_extern; i++ )
      vector_scale_factor[i+1]=np_fround(vector_scale_factor[i+1]);
  }

  for( i=0; i<num_all_var; i++ )
    myans[i]=vector_scale_factor[i+1];

  fval[0] = -fret;
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
  safe_free(vsfh);
  safe_free(num_categories_extern);
  safe_free(num_categories_extern_XY);
  safe_free(num_categories_extern_X);
  safe_free(num_categories_extern_Y);
  safe_free(bwm_kernel_unordered_vec);
  bwm_kernel_unordered_vec = NULL;
  bwm_kernel_unordered_len = 0;

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);

  free_mat(matrix_categorical_vals_extern_X, num_reg_unordered_extern + num_reg_ordered_extern);

  free_mat(matrix_categorical_vals_extern_XY, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);

  if(ibwmfunc == CBWM_CVLS)
    free_mat(matrix_categorical_vals_extern_Y, num_var_unordered_extern + num_var_ordered_extern);
  matrix_categorical_vals_extern_Y = NULL;

  safe_free(vector_continuous_stddev);

  safe_free(ipt_X);
  safe_free(ipt_XY);

  safe_free(ipt_lookup_X);
  safe_free(ipt_lookup_XY);

  if(ibwmfunc == CBWM_CVLS){
    safe_free(ipt_Y);
    safe_free(ipt_lookup_Y);
  }
  num_categories_extern_Y = NULL;

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  if(int_TREE_Y == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_Y);
    int_TREE_Y = NP_TREE_FALSE;
  }

  free_mat(matrix_XY_continuous_train_extern, num_all_cvar);
  free_mat(matrix_XY_unordered_train_extern, num_all_uvar);
  free_mat(matrix_XY_ordered_train_extern, num_all_ovar);

  if(int_TREE_XY == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_XY);
    int_TREE_XY = NP_TREE_FALSE;
  }


  int_WEIGHTS = 0;

  if(int_MINIMIZE_IO != IO_MIN_TRUE)
    Rprintf("\r                   \r");

  return ;
}

void np_distribution_conditional_bw(double * c_uno, double * c_ord, double * c_con, 
                                    double * u_uno, double * u_ord, double * u_con,
                                    double * cg_uno, double * cg_ord, double * cg_con, double * mysd,
                                    int * myopti, double * myoptd, double * myans, double * fval,
                                    double * objective_function_values, double * objective_function_evals,
                                    double * objective_function_invalid, double * timing,
                                    int * penalty_mode, double * penalty_mult){
/* Likelihood bandwidth selection for density estimation */

  double **matrix_y;

  double *vector_continuous_stddev;
  double *vsfh, *vector_scale_factor, *vector_scale_factor_multistart;

  double fret, fret_best;
  double ftol, tol;
  double (* bwmfunc)(double *) = NULL;

  double small, lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init; 
  double lbd_init, hbd_init, d_init;
  int dfc_dir;
  
  int i,j;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, num_all_var, num_var_var, iImproved;
  int itmax, iter;
  int int_use_starting_values, ibwmfunc;
  int cdfontrain;
  int scale_cat;

  int num_all_cvar, num_all_uvar, num_all_ovar;

  int * ipt_X = NULL, * ipt_XY = NULL, * ipt_Y = NULL; 
  int * ipt_lookup_XY = NULL, * ipt_lookup_Y = NULL, * ipt_lookup_X = NULL;
  int num_obs_alt;

  cdfontrain_extern = cdfontrain =  myopti[CDBW_CDFONTRAIN];

  num_var_unordered_extern = myopti[CDBW_CNUNOI];
  num_var_ordered_extern = myopti[CDBW_CNORDI];
  num_var_continuous_extern = myopti[CDBW_CNCONI];

  num_reg_unordered_extern = myopti[CDBW_UNUNOI];
  num_reg_ordered_extern = myopti[CDBW_UNORDI];
  num_reg_continuous_extern = myopti[CDBW_UNCONI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;
  num_var_var = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  num_all_var = num_var+num_var_var;

  num_obs_train_extern = myopti[CDBW_NOBSI];
  num_obs_eval_extern  = cdfontrain ? num_obs_train_extern : myopti[CDBW_NEVALI];

  iMultistart = myopti[CDBW_IMULTII];
  iNum_Multistart = myopti[CDBW_NMULTII];

  KERNEL_reg_extern = myopti[CDBW_CXKRNEVI];
  KERNEL_den_extern = myopti[CDBW_CYKRNEVI];

  KERNEL_reg_unordered_extern = myopti[CDBW_UXKRNEVI];
  KERNEL_den_unordered_extern = myopti[CDBW_UYKRNEVI];

  KERNEL_reg_ordered_extern = myopti[CDBW_OXKRNEVI];
  KERNEL_den_ordered_extern = myopti[CDBW_OYKRNEVI];

  int_use_starting_values= myopti[CDBW_USTARTI];
  int_LARGE_SF=myopti[CDBW_LSFI];
  BANDWIDTH_den_extern=myopti[CDBW_DENI];
  int_RESTART_FROM_MIN = myopti[CDBW_REMINI];
  int_MINIMIZE_IO = myopti[CDBW_MINIOI];

  itmax=myopti[CDBW_ITMAXI];

  int_TREE_XY = int_TREE_Y = int_TREE_X = myopti[CDBW_TREEI];

  scale_cat = myopti[CDBW_SCATI];
  bwm_use_transform = myopti[CDBW_TBNDI];
  if (BANDWIDTH_den_extern != BW_FIXED)
    bwm_use_transform = 0;
  if (bwm_use_transform) {
    int n = num_var_continuous_extern + num_reg_continuous_extern +
      num_var_unordered_extern + num_reg_unordered_extern +
      num_var_ordered_extern + num_reg_ordered_extern;
    if (bwm_transform_buf_len < n + 1) {
      bwm_transform_buf = (double *) realloc(bwm_transform_buf, (n + 1) * sizeof(double));
      bwm_transform_buf_len = n + 1;
    }
  }

  ftol=myoptd[CDBW_FTOLD];
  tol=myoptd[CDBW_TOLD];
  small=myoptd[CDBW_SMALLD];
  dbl_memfac_ccdf_extern = myoptd[CDBW_MEMFACD];

  dfc_dir = myopti[CDBW_DFC_DIRI];
  lbc_dir = myoptd[CDBW_LBC_DIRD];
  c_dir = myoptd[CDBW_C_DIRD];
  initc_dir = myoptd[CDBW_INITC_DIRD]; 

  lbd_dir = myoptd[CDBW_LBD_DIRD]; 
  hbd_dir = myoptd[CDBW_HBD_DIRD]; 
  d_dir = myoptd[CDBW_D_DIRD]; 
  initd_dir = myoptd[CDBW_INITD_DIRD]; 

  lbc_init = myoptd[CDBW_LBC_INITD]; 
  hbc_init = myoptd[CDBW_HBC_INITD]; 
  c_init = myoptd[CDBW_C_INITD]; 

  lbd_init = myoptd[CDBW_LBD_INITD]; 
  hbd_init = myoptd[CDBW_HBD_INITD]; 
  d_init = myoptd[CDBW_D_INITD]; 

  nconfac_extern = myoptd[CDBW_NCONFD];
  ncatfac_extern = myoptd[CDBW_NCATFD];

/* Allocate memory for objects */

  matrix_Y_unordered_train_extern = alloc_matd(num_obs_train_extern, num_var_unordered_extern);
  matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, num_var_ordered_extern);
  matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, num_var_continuous_extern);

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);

  if(cdfontrain){
    matrix_Y_unordered_eval_extern = matrix_Y_unordered_train_extern;
    matrix_Y_ordered_eval_extern = matrix_Y_ordered_train_extern;
    matrix_Y_continuous_eval_extern = matrix_Y_continuous_train_extern;
  } else {
    matrix_Y_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_var_unordered_extern);
    matrix_Y_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_var_ordered_extern);
    matrix_Y_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_var_continuous_extern);
  }

  num_categories_extern = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                     num_reg_unordered_extern + num_reg_ordered_extern);

  num_categories_extern_X = alloc_vecu(num_reg_unordered_extern + num_reg_ordered_extern);
  num_categories_extern_Y = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern);

  num_categories_extern_XY = alloc_vecu(num_var_unordered_extern + num_var_ordered_extern +
                                        num_reg_unordered_extern + num_reg_ordered_extern);
  
  matrix_y = alloc_matd(num_all_var + 1, num_all_var + 1);
  vector_scale_factor = alloc_vecd(num_all_var + 1);
  vsfh = alloc_vecd(num_all_var + 1);
  
  matrix_categorical_vals_extern = 
    alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern + 
               num_reg_unordered_extern + num_reg_ordered_extern);

  matrix_categorical_vals_extern_X = 
    alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  matrix_categorical_vals_extern_Y = 
    alloc_matd(num_obs_train_extern, num_var_unordered_extern + num_var_ordered_extern);

  matrix_categorical_vals_extern_XY = 
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

  if(!cdfontrain){
    for(j=0;j<num_var_unordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_Y_unordered_eval_extern[j][i]=cg_uno[j*num_obs_eval_extern+i];

    for(j=0;j<num_var_ordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_Y_ordered_eval_extern[j][i]=cg_ord[j*num_obs_eval_extern+i];

    for(j=0;j<num_var_continuous_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_Y_continuous_eval_extern[j][i]=cg_con[j*num_obs_eval_extern+i];

  }

  num_all_cvar = num_reg_continuous_extern + num_var_continuous_extern;
  num_all_uvar = num_reg_unordered_extern + num_var_unordered_extern;
  num_all_ovar = num_reg_ordered_extern + num_var_ordered_extern;

  // we need 3 trees to accelerate cg_concdf :)

  ipt_X = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_X != NULL))
    error("!(ipt_X != NULL)");

  ipt_lookup_X = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_lookup_X != NULL))
    error("!(ipt_lookup_X != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt_lookup_X[i] = ipt_X[i] = i;
  }

  ipt_extern_X = ipt_X;
  ipt_lookup_extern_X = ipt_lookup_X;


  ipt_Y = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_Y != NULL))
    error("!(ipt_Y != NULL)");

  ipt_lookup_Y = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_lookup_Y != NULL))
    error("!(ipt_lookup_Y != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt_lookup_Y[i] = ipt_Y[i] = i;
  }

  ipt_extern_Y = ipt_Y;
  ipt_lookup_extern_Y = ipt_lookup_Y;

  num_obs_alt = (BANDWIDTH_den_extern != BW_ADAP_NN) ? num_obs_train_extern : 0;

  ipt_XY = (int *)malloc(num_obs_alt*sizeof(int));
  if(!(ipt_XY != NULL))
    error("!(ipt_XY != NULL)");

  ipt_lookup_XY = (int *)malloc(num_obs_alt*sizeof(int));
  if(!(ipt_lookup_XY != NULL))
    error("!(ipt_lookup_XY != NULL)");

  for(i = 0; i < num_obs_alt; i++){
    ipt_lookup_XY[i] = ipt_XY[i] = i;
  }

  ipt_extern_XY = ipt_XY;
  ipt_lookup_extern_XY = ipt_lookup_XY;

  int_TREE_XY = int_TREE_XY && (((num_all_cvar) != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);

  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);

  int_TREE_Y = int_TREE_Y && ((num_var_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE) && (BANDWIDTH_den_extern != BW_ADAP_NN);

  if(int_TREE_X == NP_TREE_TRUE){
    build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                 4*num_reg_continuous_extern, ipt_X, &kdt_extern_X);
  
    // put x data into x-tree order
    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+ipt_X[i]];
    
    
    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+ipt_X[i]];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+ipt_X[i]];

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_X[ipt_X[i]] = i;
    }

  }

  if(int_TREE_Y == NP_TREE_TRUE){
    build_kdtree(matrix_Y_continuous_train_extern, num_obs_train_extern, num_var_continuous_extern, 
                 4*num_var_continuous_extern, ipt_Y, &kdt_extern_Y);
  
    // put y data into y-tree order

    for(j=0;j<num_var_unordered_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_unordered_train_extern[j][i]=c_uno[j*num_obs_train_extern+ipt_Y[i]];

    for(j=0;j<num_var_ordered_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_ordered_train_extern[j][i]=c_ord[j*num_obs_train_extern+ipt_Y[i]];

    for(j=0;j<num_var_continuous_extern;j++)
      for(i=0;i<num_obs_train_extern;i++)
        matrix_Y_continuous_train_extern[j][i]=c_con[j*num_obs_train_extern+ipt_Y[i]];

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_Y[ipt_Y[i]] = i;
    }
    
  }

  matrix_XY_continuous_train_extern = alloc_matd(num_obs_alt, num_all_cvar);
  matrix_XY_unordered_train_extern = alloc_matd(num_obs_alt, num_all_uvar);
  matrix_XY_ordered_train_extern = alloc_matd(num_obs_alt, num_all_ovar);

  if(int_TREE_XY == NP_TREE_TRUE){

    for(j = 0; j < num_reg_unordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+i];

    for(j = num_reg_unordered_extern; j < num_all_uvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=c_uno[(j-num_reg_unordered_extern)*num_obs_train_extern+i];


    for(j = 0; j < num_reg_ordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+i];

    for(j = num_reg_ordered_extern; j < num_all_ovar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=c_ord[(j-num_reg_ordered_extern)*num_obs_train_extern+i];


    for(j = 0; j < num_reg_continuous_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+i];

    for(j = num_reg_continuous_extern; j < num_all_cvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=c_con[(j-num_reg_continuous_extern)*num_obs_train_extern+i];

  // XY tree!

    build_kdtree(matrix_XY_continuous_train_extern, num_obs_train_extern, num_all_cvar, 
                 4*num_all_cvar, ipt_XY, &kdt_extern_XY);

    // put data into xy-tree order
    for(j = 0; j < num_reg_unordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=u_uno[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_unordered_extern; j < num_all_uvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_unordered_train_extern[j][i]=c_uno[(j-num_reg_unordered_extern)*num_obs_train_extern+ipt_XY[i]];


    for(j = 0; j < num_reg_ordered_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=u_ord[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_ordered_extern; j < num_all_ovar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_ordered_train_extern[j][i]=c_ord[(j-num_reg_ordered_extern)*num_obs_train_extern+ipt_XY[i]];


    for(j = 0; j < num_reg_continuous_extern; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=u_con[j*num_obs_train_extern+ipt_XY[i]];

    for(j = num_reg_continuous_extern; j < num_all_cvar; j++)
      for(i = 0; i < num_obs_train_extern; i++)
        matrix_XY_continuous_train_extern[j][i]=c_con[(j-num_reg_continuous_extern)*num_obs_train_extern+ipt_XY[i]];

    for(i = 0; i < num_obs_train_extern; i++){
      ipt_lookup_XY[ipt_XY[i]] = i;
    }
    
  }

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

  np_splitxy_vsf_mcv_nc(num_var_unordered_extern, num_var_ordered_extern, num_var_continuous_extern,
                        num_reg_unordered_extern, num_reg_ordered_extern, num_reg_continuous_extern,
                        vector_scale_factor+1,
                        num_categories_extern,
                        matrix_categorical_vals_extern,
                        NULL, NULL, NULL,
                        num_categories_extern_X, num_categories_extern_Y, num_categories_extern_XY,
                        matrix_categorical_vals_extern_X, matrix_categorical_vals_extern_Y, matrix_categorical_vals_extern_XY);


  vector_continuous_stddev = vector_continuous_stddev_extern = mysd;


  /* Initialize scale factors and Directions for NR modules */

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    num_var_continuous_extern,
                                    num_var_unordered_extern,
                                    num_var_ordered_extern,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    KERNEL_den_unordered_extern,
                                    KERNEL_reg_unordered_extern,
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    num_var_continuous_extern,
                                    num_var_unordered_extern,
                                    num_var_ordered_extern,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    KERNEL_den_unordered_extern,
                                    KERNEL_reg_unordered_extern,
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vsfh,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_den_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           num_var_continuous_extern,
                           num_var_unordered_extern,
                           num_var_ordered_extern,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED,  
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);


  /* When multistarting, set counter */

  imsnum = iMs_counter = 0;
  imstot = iNum_Multistart;

  /* Conduct direction set search */

  /* assign the function to be optimized */

  ibwmfunc = myopti[CDBW_MI];

  switch(ibwmfunc){
  case CDBWM_CVLS : bwmfunc = cv_func_con_distribution_categorical_ls; break;
  default : REprintf("np.c: invalid bandwidth selection method.");
    error("np.c: invalid bandwidth selection method."); break;
  }

  if (bwm_use_transform)
    bwm_to_unconstrained(vector_scale_factor, num_all_var);

  spinner(0);

  bwmfunc_raw = bwmfunc;
  bwm_num_reg_continuous = num_var_continuous_extern + num_reg_continuous_extern;
  bwm_num_reg_unordered = num_var_unordered_extern + num_reg_unordered_extern;
  bwm_num_reg_ordered = num_var_ordered_extern + num_reg_ordered_extern;
  bwm_kernel_unordered = KERNEL_den_unordered_extern;
  bwm_kernel_unordered_len = bwm_num_reg_unordered;
  if (bwm_kernel_unordered_len > 0) {
    bwm_kernel_unordered_vec = alloc_vecu(bwm_kernel_unordered_len);
    for (i = 0; i < num_var_unordered_extern; i++)
      bwm_kernel_unordered_vec[i] = KERNEL_den_unordered_extern;
    for (i = 0; i < num_reg_unordered_extern; i++)
      bwm_kernel_unordered_vec[num_var_unordered_extern + i] = KERNEL_reg_unordered_extern;
  } else {
    bwm_kernel_unordered_vec = NULL;
  }
  bwm_num_categories = num_categories_extern;
  bwm_reset_counters();
  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  if (penalty_mode[0] == 1) {
    double pmult = penalty_mult[0];
    double baseline;
    if (pmult < 1.0) pmult = 1.0;
    baseline = bwmfunc_raw(vector_scale_factor);
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      double *tmp = alloc_vecd(num_all_var + 1);
      for (i = 1; i <= num_all_var; i++)
        tmp[i] = vector_scale_factor[i];
      for (i = 1; i <= (num_var_continuous_extern + num_reg_continuous_extern); i++)
        tmp[i] *= 2.0;
      for (i = 0; i < (num_var_unordered_extern + num_reg_unordered_extern); i++) {
        int idx = num_var_continuous_extern + num_reg_continuous_extern + 1 + i;
        double maxbw = max_unordered_bw(num_categories_extern[i], KERNEL_den_unordered_extern);
        tmp[idx] = 0.5*maxbw;
      }
      for (i = 0; i < (num_var_ordered_extern + num_reg_ordered_extern); i++) {
        int idx = num_var_continuous_extern + num_reg_continuous_extern +
          num_var_unordered_extern + num_reg_unordered_extern + 1 + i;
        tmp[idx] = 0.5;
      }
      baseline = bwmfunc_raw(tmp);
      safe_free(tmp);
    }
    if (!R_FINITE(baseline) || baseline == DBL_MAX) {
      bwm_penalty_value = pmult * 1.0e6;
    } else {
      bwm_penalty_value = baseline + (fabs(baseline) + 1.0) * pmult;
    }
    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }

  fret_best = bwmfunc_wrapper(vector_scale_factor);
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
         bwmfunc_wrapper);

  if(int_RESTART_FROM_MIN == RE_MIN_TRUE){
    initialize_nr_directions(BANDWIDTH_den_extern,
                             num_obs_train_extern,
                             num_reg_continuous_extern,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             num_var_continuous_extern,
                             num_var_unordered_extern,
                             num_var_ordered_extern,
                             vsfh,
                             num_categories_extern,
                             matrix_y,
                             0, int_RANDOM_SEED,  
                             lbc_dir, dfc_dir, c_dir, initc_dir,
                             lbd_dir, hbd_dir, d_dir, initd_dir,
                             matrix_X_continuous_train_extern,
                             matrix_Y_continuous_train_extern);

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
           bwmfunc_wrapper);

  }

  iImproved = (fret < fret_best);
  *timing = timing_extern;

  objective_function_values[0]=fret;
  objective_function_evals[0]=bwm_eval_count;
  objective_function_invalid[0]=bwm_invalid_count;
  /* When multistarting save initial minimum of objective function and scale factors */


  if(iMultistart == IMULTI_TRUE){
    fret_best = fret;
    vector_scale_factor_multistart = alloc_vecd(num_all_var + 1);
    for(i = 1; i <= num_all_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
			

    /* Conduct search from new random values of the search parameters */
		
    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and directions for NR modules */
      initialize_nr_vector_scale_factor(BANDWIDTH_den_extern,
                                        1,                /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        num_var_continuous_extern,
                                        num_var_unordered_extern,
                                        num_var_ordered_extern,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        KERNEL_den_unordered_extern,
                                        KERNEL_reg_unordered_extern,
                                        0,
                                        scale_cat,
                                        pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                        nconfac_extern, ncatfac_extern,
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor,
                                        lbc_init, hbc_init, c_init, 
                                        lbd_init, hbd_init, d_init,
                                        matrix_X_continuous_train_extern,
                                        matrix_Y_continuous_train_extern);

      initialize_nr_directions(BANDWIDTH_den_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               num_var_continuous_extern,
                               num_var_unordered_extern,
                               num_var_ordered_extern,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               1, int_RANDOM_SEED,  
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);


      /* Conduct direction set search */

      bwm_reset_counters();
      
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
             bwmfunc_wrapper);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

        initialize_nr_directions(BANDWIDTH_den_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 num_reg_unordered_extern,
                                 num_reg_ordered_extern,
                                 num_var_continuous_extern,
                                 num_var_unordered_extern,
                                 num_var_ordered_extern,
                                 vsfh,
                                 num_categories_extern,
                                 matrix_y, 
                                 0, int_RANDOM_SEED,  
                                 lbc_dir, dfc_dir, c_dir, initc_dir,
                                 lbd_dir, hbd_dir, d_dir, initd_dir,
                                 matrix_X_continuous_train_extern,
                                 matrix_Y_continuous_train_extern);


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
               bwmfunc_wrapper);
      }
				
      /* If this run resulted in an improved minimum save information */
      
      if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        *timing = timing_extern;
        
        for(i = 1; i <= num_all_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }
      objective_function_values[iMs_counter]=fret;
      objective_function_evals[iMs_counter]=bwm_eval_count;
      objective_function_invalid[iMs_counter]=bwm_invalid_count;
    }

    /* Save best for estimation */

    fret = fret_best;
    for(i = 1; i <= num_all_var; i++)
      vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];
    free(vector_scale_factor_multistart);
  }

  if (bwm_use_transform)
    bwm_to_constrained(vector_scale_factor, num_all_var);

  /* return data to R */
  if (BANDWIDTH_den_extern == BW_GEN_NN || 
      BANDWIDTH_den_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern+num_var_continuous_extern; i++ )
      vector_scale_factor[i+1]=np_fround(vector_scale_factor[i+1]);
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

  if(!cdfontrain){
    free_mat(matrix_Y_unordered_eval_extern, num_var_unordered_extern);
    free_mat(matrix_Y_ordered_eval_extern, num_var_ordered_extern);
    free_mat(matrix_Y_continuous_eval_extern, num_var_continuous_extern);
  }

  free_mat(matrix_y, num_all_var + 1);
  safe_free(vector_scale_factor);
  safe_free(vsfh);
  safe_free(num_categories_extern);
  safe_free(num_categories_extern_X);
  safe_free(num_categories_extern_Y);
  safe_free(num_categories_extern_XY);
  safe_free(bwm_kernel_unordered_vec);
  bwm_kernel_unordered_vec = NULL;
  bwm_kernel_unordered_len = 0;

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);

  free_mat(matrix_categorical_vals_extern_X, num_reg_unordered_extern + num_reg_ordered_extern);

  free_mat(matrix_categorical_vals_extern_Y, num_var_unordered_extern + num_var_ordered_extern);
  free_mat(matrix_categorical_vals_extern_XY, num_reg_unordered_extern + num_reg_ordered_extern +
           num_var_unordered_extern + num_var_ordered_extern);

  safe_free(ipt_X);
  safe_free(ipt_Y);
  safe_free(ipt_XY);

  safe_free(ipt_lookup_X);
  safe_free(ipt_lookup_Y);
  safe_free(ipt_lookup_XY);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

 if(int_TREE_Y == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_Y);
    int_TREE_Y = NP_TREE_FALSE;
  }

  free_mat(matrix_XY_continuous_train_extern, num_all_cvar);
  free_mat(matrix_XY_unordered_train_extern, num_all_uvar);
  free_mat(matrix_XY_ordered_train_extern, num_all_ovar);

  if(int_TREE_XY == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_XY);
    int_TREE_XY = NP_TREE_FALSE;
  }


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
                            double * nconfac, double * ncatfac, double * mysd,
                            int * myopti, 
                            double * cdens, double * cderr, 
                            double * cg, double * cgerr,
                            double * ll){
  /* Likelihood bandwidth selection for density estimation */

  double *vector_scale_factor, *pdf, *pdf_stderr, log_likelihood = 0.0;
  double ** pdf_deriv = NULL, ** pdf_deriv_stderr = NULL;
  double xpad_num, ypad_num;

  int i,j;
  int num_var;

  int num_all_var, num_var_var, train_is_eval, do_grad, num_obs_eval_alloc;
  int num_all_cvar, num_all_uvar, num_all_ovar, num_all_catvar;
  int xmax_lev, ymax_lev, dens_or_dist, t_num;

  int * ipt_XY = NULL, *ipe_XY = NULL;
  int operator;


  num_var_unordered_extern = myopti[CD_CNUNOI];
  num_var_ordered_extern = myopti[CD_CNORDI];
  num_var_continuous_extern = myopti[CD_CNCONI];

  num_reg_unordered_extern = myopti[CD_UNUNOI];
  num_reg_ordered_extern = myopti[CD_UNORDI];
  num_reg_continuous_extern = myopti[CD_UNCONI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;
  num_var_var = num_var_continuous_extern + num_var_unordered_extern + num_var_ordered_extern;
  num_all_var = num_var + num_var_var;

  num_all_cvar = num_reg_continuous_extern + num_var_continuous_extern;
  num_all_uvar = num_reg_unordered_extern + num_var_unordered_extern;
  num_all_ovar = num_reg_ordered_extern + num_var_ordered_extern;

  num_all_catvar = num_all_uvar + num_all_ovar;

  num_obs_train_extern = myopti[CD_TNOBSI];
  num_obs_eval_extern = myopti[CD_ENOBSI];

  if((train_is_eval = myopti[CD_TISEI]) && 
     (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_density_conditional): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    error("\n(np_density_conditional): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
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

  nconfac_extern = *nconfac;
  ncatfac_extern = *ncatfac;

  dens_or_dist = myopti[CD_DODENI];

  int_TREE_XY = myopti[CD_TREEI]; // we just build a single xy tree
  
  int_TREE_X = int_TREE_Y = NP_TREE_FALSE;

  operator = (dens_or_dist == NP_DO_DENS) ? OP_NORMAL : OP_INTEGRAL;

#ifdef MPI2
  num_obs_eval_alloc = MAX(ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif

  // our method of evaluation uses a single joint xy data matrix

  matrix_XY_continuous_train_extern = alloc_matd(num_obs_train_extern, num_all_cvar);
  matrix_XY_unordered_train_extern = alloc_matd(num_obs_train_extern, num_all_uvar);
  matrix_XY_ordered_train_extern = alloc_matd(num_obs_train_extern, num_all_ovar);

  if(train_is_eval) {
    matrix_XY_continuous_eval_extern = matrix_XY_continuous_train_extern;
    matrix_XY_unordered_eval_extern = matrix_XY_unordered_train_extern;
    matrix_XY_ordered_eval_extern = matrix_XY_ordered_train_extern;
  } else {
    matrix_XY_continuous_eval_extern = alloc_matd(num_obs_eval_extern, num_all_cvar);
    matrix_XY_unordered_eval_extern = alloc_matd(num_obs_eval_extern, num_all_uvar);
    matrix_XY_ordered_eval_extern = alloc_matd(num_obs_eval_extern, num_all_ovar);
  }

  num_categories_extern = alloc_vecu(num_all_catvar);
	
  num_categories_extern_XY = alloc_vecu(num_all_catvar);

  vector_scale_factor = alloc_vecd(num_all_var + 1);
  
  matrix_categorical_vals_extern = alloc_matd(MAX(xmax_lev,ymax_lev), num_all_catvar);
  matrix_categorical_vals_extern_XY = alloc_matd(MAX(xmax_lev,ymax_lev), num_all_catvar);


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

  vector_continuous_stddev_extern = mysd;

  /* Parse data */

  // xy arrays contain first x, then y data

  for(j=0;j<num_reg_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_XY_unordered_train_extern[j][i]=tu_uno[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_ordered_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_XY_ordered_train_extern[j][i]=tu_ord[j*num_obs_train_extern+i];

  for( j=0;j<num_reg_continuous_extern;j++)
    for( i=0;i<num_obs_train_extern;i++ )
      matrix_XY_continuous_train_extern[j][i]=tu_con[j*num_obs_train_extern+i];

  // y data
  for(j=0;j<num_var_unordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_XY_unordered_train_extern[j+num_reg_unordered_extern][i]=tc_uno[j*num_obs_train_extern+i];

  for(j=0;j<num_var_ordered_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_XY_ordered_train_extern[j+num_reg_ordered_extern][i]=tc_ord[j*num_obs_train_extern+i];

  for(j=0;j<num_var_continuous_extern;j++)
    for(i=0;i<num_obs_train_extern;i++)
      matrix_XY_continuous_train_extern[j+num_reg_continuous_extern][i]=tc_con[j*num_obs_train_extern+i];

  /* eval */
  if(!train_is_eval){
    for(j=0;j<num_reg_unordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_XY_unordered_eval_extern[j][i]=eu_uno[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_XY_ordered_eval_extern[j][i]=eu_ord[j*num_obs_eval_extern+i];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_eval_extern;i++ )
        matrix_XY_continuous_eval_extern[j][i]=eu_con[j*num_obs_eval_extern+i];

    // y data
    for(j=0;j<num_var_unordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_XY_unordered_eval_extern[j+num_reg_unordered_extern][i]=ec_uno[j*num_obs_eval_extern+i];

    for(j=0;j<num_var_ordered_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_XY_ordered_eval_extern[j+num_reg_ordered_extern][i]=ec_ord[j*num_obs_eval_extern+i];

    for(j=0;j<num_var_continuous_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        matrix_XY_continuous_eval_extern[j+num_reg_continuous_extern][i]=ec_con[j*num_obs_eval_extern+i];
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

  // properly fill-in xy categorical data
  np_splitxy_vsf_mcv_nc(num_var_unordered_extern, num_var_ordered_extern, num_var_continuous_extern,
                        num_reg_unordered_extern, num_reg_ordered_extern, num_reg_continuous_extern,
                        vector_scale_factor+1,
                        num_categories_extern,
                        matrix_categorical_vals_extern,
                        NULL, NULL, NULL,
                        NULL, NULL, num_categories_extern_XY,
                        NULL, NULL, matrix_categorical_vals_extern_XY);

  
  // set up indexing data for tree
  ipt_XY = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt_XY != NULL))
    error("!(ipt_XY != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt_XY[i] = i;
  }

  if(!train_is_eval) {
    ipe_XY = (int *)malloc(num_obs_eval_extern*sizeof(int));
    if(!(ipe_XY != NULL))
      error("!(ipe_XY != NULL)");

    for(i = 0; i < num_obs_eval_extern; i++){
      ipe_XY[i] = i;
    }
  } else {
    ipe_XY = ipt_XY;
  }

  int_TREE_XY = int_TREE_XY && ((num_all_cvar != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_XY == NP_TREE_TRUE){
    if((BANDWIDTH_den_extern != BW_ADAP_NN) || ((BANDWIDTH_den_extern == BW_ADAP_NN) && train_is_eval)){
      build_kdtree(matrix_XY_continuous_train_extern, num_obs_train_extern, num_all_cvar, 
                   4*num_all_cvar, ipt_XY, &kdt_extern_XY);

      // x
      for(j = 0; j < num_reg_unordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_unordered_train_extern[j][i]=tu_uno[j*num_obs_train_extern+ipt_XY[i]];
        
      for(j = 0; j < num_reg_ordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_ordered_train_extern[j][i]=tu_ord[j*num_obs_train_extern+ipt_XY[i]];

      for(j = 0; j < num_reg_continuous_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_continuous_train_extern[j][i]=tu_con[j*num_obs_train_extern+ipt_XY[i]];
    
      // y
      for(j = 0; j < num_var_unordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_unordered_train_extern[j+num_reg_unordered_extern][i]=tc_uno[j*num_obs_train_extern+ipt_XY[i]];
        
      for(j = 0; j < num_var_ordered_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_ordered_train_extern[j+num_reg_ordered_extern][i]=tc_ord[j*num_obs_train_extern+ipt_XY[i]];

      for(j = 0; j < num_var_continuous_extern; j++)
        for(i = 0; i < num_obs_train_extern; i++)
          matrix_XY_continuous_train_extern[j+num_reg_continuous_extern][i]=tc_con[j*num_obs_train_extern+ipt_XY[i]];
    } else {
      build_kdtree(matrix_XY_continuous_eval_extern, num_obs_eval_extern, num_all_cvar, 
                   4*num_all_cvar, ipe_XY, &kdt_extern_XY);

      // x
      for(j = 0; j < num_reg_unordered_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_unordered_eval_extern[j][i]=eu_uno[j*num_obs_eval_extern+ipe_XY[i]];
        
      for(j = 0; j < num_reg_ordered_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_ordered_eval_extern[j][i]=eu_ord[j*num_obs_eval_extern+ipe_XY[i]];

      for(j = 0; j < num_reg_continuous_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_continuous_eval_extern[j][i]=eu_con[j*num_obs_eval_extern+ipe_XY[i]];
    
      // y
      for(j = 0; j < num_var_unordered_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_unordered_eval_extern[j+num_reg_unordered_extern][i]=ec_uno[j*num_obs_eval_extern+ipe_XY[i]];
        
      for(j = 0; j < num_var_ordered_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_ordered_eval_extern[j+num_reg_ordered_extern][i]=ec_ord[j*num_obs_eval_extern+ipe_XY[i]];

      for(j = 0; j < num_var_continuous_extern; j++)
        for(i = 0; i < num_obs_eval_extern; i++)
          matrix_XY_continuous_eval_extern[j+num_reg_continuous_extern][i]=ec_con[j*num_obs_eval_extern+ipe_XY[i]];

    }
  }

  np_kernel_estimate_con_dens_dist_categorical(KERNEL_den_extern,
                                               KERNEL_den_unordered_extern,
                                               KERNEL_den_ordered_extern,
                                               KERNEL_reg_extern,
                                               KERNEL_reg_unordered_extern,
                                               KERNEL_reg_ordered_extern,
                                               BANDWIDTH_den_extern,
                                               operator,
                                               num_obs_train_extern,
                                               num_obs_eval_extern,
                                               num_var_unordered_extern,
                                               num_var_ordered_extern,
                                               num_var_continuous_extern,
                                               num_reg_unordered_extern,
                                               num_reg_ordered_extern,
                                               num_reg_continuous_extern,
                                               matrix_XY_unordered_train_extern, 
                                               matrix_XY_ordered_train_extern, 
                                               matrix_XY_continuous_train_extern, 
                                               matrix_XY_unordered_eval_extern, 
                                               matrix_XY_ordered_eval_extern, 
                                               matrix_XY_continuous_eval_extern, 
                                               &vector_scale_factor[1],
                                               num_categories_extern,
                                               num_categories_extern_XY,
                                               matrix_categorical_vals_extern,
                                               matrix_categorical_vals_extern_XY,
                                               pdf,
                                               pdf_stderr,
                                               pdf_deriv,
                                               pdf_deriv_stderr,
                                               &log_likelihood);


  /* return data to R */
  for( i=0; i<num_obs_eval_extern; i++ )
    cdens[ipe_XY[i]]=pdf[i];

  for( i=0; i<num_obs_eval_extern; i++ )
    cderr[ipe_XY[i]]=pdf_stderr[i];
  
  if (do_grad) {
    for(j=0;j<num_var;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        cgerr[j*num_obs_eval_extern+ipe_XY[i]]=pdf_deriv_stderr[j][i];

    for(j=0;j<num_var;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        cg[j*num_obs_eval_extern+ipe_XY[i]]=pdf_deriv[j][i];
  }



  *ll = log_likelihood;
  /* end return data */

  /* Free data objects */

  free_mat(matrix_XY_unordered_train_extern, num_all_uvar);
  free_mat(matrix_XY_ordered_train_extern, num_all_ovar);
  free_mat(matrix_XY_continuous_train_extern, num_all_cvar);

  if(train_is_eval){
    matrix_XY_unordered_eval_extern = NULL;
    matrix_XY_ordered_eval_extern = NULL;
    matrix_XY_continuous_eval_extern = NULL;
  } else {
    free_mat(matrix_XY_unordered_eval_extern, num_all_uvar);
    free_mat(matrix_XY_ordered_eval_extern, num_all_ovar);
    free_mat(matrix_XY_continuous_eval_extern, num_all_cvar);
  }

  if (do_grad){
    free_mat(pdf_deriv, num_var);
    free_mat(pdf_deriv_stderr, num_var);
  }

  vector_continuous_stddev_extern = NULL;

  safe_free(vector_scale_factor);
  safe_free(num_categories_extern);
  safe_free(num_categories_extern_XY);
  safe_free(pdf);
  safe_free(pdf_stderr);

  free_mat(matrix_categorical_vals_extern, num_all_catvar);
  free_mat(matrix_categorical_vals_extern_XY, num_all_catvar);

  safe_free(ipt_XY);

  if(!train_is_eval)
    safe_free(ipe_XY);

  if(int_TREE_XY == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_XY);
    int_TREE_XY = NP_TREE_FALSE;
  }


  return;
}


void np_density(double * tuno, double * tord, double * tcon, 
                double * euno, double * eord, double * econ, 
                double * dbw, 
                double * mcv, double * padnum, 
                double * nconfac, double * ncatfac, double * mysd,
                int * myopti, double * mydens, double * myderr, double * ll){


  double small = 1.0e-16;
  double * vector_scale_factor, * pdf, * pdf_stderr, log_likelihood = 0.0;
  double pad_num;

  int itmax = 10000;
  int i,j;
  int num_var, num_obs_eval_alloc, max_lev, train_is_eval, dens_or_dist, old_dens;

  int * ipt = NULL, * ipe = NULL;
  

  /* match integer options with their globals */

  num_reg_continuous_extern = myopti[DEN_NCONI];
  num_reg_unordered_extern = myopti[DEN_NUNOI];
  num_reg_ordered_extern = myopti[DEN_NORDI];

  num_var = num_reg_ordered_extern + num_reg_continuous_extern + num_reg_unordered_extern;

  num_obs_train_extern = myopti[DEN_TNOBSI];
  num_obs_eval_extern = myopti[DEN_ENOBSI];

  KERNEL_den_extern = myopti[DEN_CKRNEVI];
  KERNEL_den_unordered_extern = myopti[DEN_UKRNEVI];
  KERNEL_den_ordered_extern = myopti[DEN_OKRNEVI];

  int_LARGE_SF = myopti[DEN_LSFI];
  int_MINIMIZE_IO = myopti[DEN_MINIOI];
  BANDWIDTH_den_extern = myopti[DEN_DENI];

  train_is_eval = myopti[DEN_TISEI];

  max_lev = myopti[DEN_MLEVI];
  pad_num = *padnum;

  nconfac_extern = *nconfac;
  ncatfac_extern = *ncatfac;

  dens_or_dist = myopti[DEN_DODENI];
  old_dens = myopti[DEN_OLDI];
  int_TREE_X = myopti[DEN_TREEI];

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
  
  vector_continuous_stddev_extern = mysd;

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

  /* data has been copied, now build tree */

  ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  if(!train_is_eval) {
    ipe = (int *)malloc(num_obs_eval_extern*sizeof(int));
    if(!(ipe != NULL))
      error("!(ipe != NULL)");

    for(i = 0; i < num_obs_eval_extern; i++){
      ipe[i] = i;
    }
  } else {
    ipe = ipt;
  }

  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    if((BANDWIDTH_den_extern != BW_ADAP_NN) || ((BANDWIDTH_den_extern == BW_ADAP_NN) && train_is_eval)){
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipt, &kdt_extern_X);

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+ipt[i]];
    
    
      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+ipt[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+ipt[i]];

    } else {
      build_kdtree(matrix_X_continuous_eval_extern, num_obs_eval_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipe, &kdt_extern_X);


      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+ipe[i]];
    }

  }


  /* Conduct estimation */
  
  if(old_dens){
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
  } else {
    const int dop = (dens_or_dist == NP_DO_DENS) ? OP_NORMAL : OP_INTEGRAL;

      kernel_estimate_dens_dist_categorical_np(KERNEL_den_extern,
                                               KERNEL_den_unordered_extern,
                                               KERNEL_den_ordered_extern,
                                               BANDWIDTH_den_extern,
                                               num_obs_train_extern,
                                               num_obs_eval_extern,
                                               num_reg_unordered_extern,
                                               num_reg_ordered_extern,
                                               num_reg_continuous_extern,
                                               dop,
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
                                               &log_likelihood);
  }
  
  
  /* write the return values */

  for(i=0;i<num_obs_eval_extern;i++){
    mydens[ipe[i]] = pdf[i];
    myderr[ipe[i]] = pdf_stderr[i];
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

  vector_continuous_stddev_extern = NULL;

  safe_free(vector_scale_factor);
  safe_free(num_categories_extern);
  safe_free(pdf_stderr);
  safe_free(pdf);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  safe_free(ipt);

  if(!train_is_eval)
    safe_free(ipe);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  return;
}


void np_regression_bw(double * runo, double * rord, double * rcon, double * y,
                      double * mysd, int * myopti, double * myoptd, double * rbw, double * fval,
                      double * objective_function_values, double * objective_function_evals,
                      double * objective_function_invalid, double * timing,
                      int * penalty_mode, double * penalty_mult){
  //KDT * kdt = NULL; // tree structure
  //NL nl = { .node = NULL, .n = 0, .nalloc = 0 };// a node list structure -- used for searching - here for testing
  //double tb[4] = {0.25, 0.5, 0.3, 0.75};
  int * ipt = NULL;  // point permutation, see tree.c

  double **matrix_y;

  double *vector_continuous_stddev;
  double *vector_scale_factor, *vector_scale_factor_multistart, * vsfh;

  double fret, fret_best;
  double ftol, tol, small;
  double (* bwmfunc)(double *) = NULL;

  double lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  double lbc_init, hbc_init, c_init; 
  double lbd_init, hbd_init, d_init;
  int dfc_dir;

  int i,j;
  int num_var;
  int iMultistart, iMs_counter, iNum_Multistart, iImproved;
  int itmax, iter;
  int int_use_starting_values;

  int scale_cat;

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
#ifdef MPI2
  {
    int ll_min = 0, ll_max = 0;
    MPI_Allreduce(&int_ll_extern, &ll_min, 1, MPI_INT, MPI_MIN, comm[1]);
    MPI_Allreduce(&int_ll_extern, &ll_max, 1, MPI_INT, MPI_MAX, comm[1]);
    if(ll_min != ll_max){
      if(my_rank == 0){
        REprintf("\n[npRmpi] Warning: inconsistent regression type across ranks (min=%d max=%d). Forcing all ranks to min.\n", ll_min, ll_max);
        R_FlushConsole();
      }
      int_ll_extern = ll_min;
    }
  }
#endif

  int_TREE_X = myopti[RBW_DOTREEI];
  scale_cat = myopti[RBW_SCATI];
  bwm_use_transform = myopti[RBW_TBNDI];
  if (BANDWIDTH_reg_extern != BW_FIXED)
    bwm_use_transform = 0;

  ftol=myoptd[RBW_FTOLD];
  tol=myoptd[RBW_TOLD];
  small=myoptd[RBW_SMALLD];

  dfc_dir = myopti[RBW_DFC_DIRI];
  lbc_dir = myoptd[RBW_LBC_DIRD];
  c_dir = myoptd[RBW_C_DIRD];
  initc_dir = myoptd[RBW_INITC_DIRD]; 

  lbd_dir = myoptd[RBW_LBD_DIRD]; 
  hbd_dir = myoptd[RBW_HBD_DIRD]; 
  d_dir = myoptd[RBW_D_DIRD]; 
  initd_dir = myoptd[RBW_INITD_DIRD]; 

  lbc_init = myoptd[RBW_LBC_INITD]; 
  hbc_init = myoptd[RBW_HBC_INITD]; 
  c_init = myoptd[RBW_C_INITD]; 

  lbd_init = myoptd[RBW_LBD_INITD]; 
  hbd_init = myoptd[RBW_HBD_INITD]; 
  d_init = myoptd[RBW_D_INITD]; 

  nconfac_extern = myoptd[RBW_NCONFD];
  ncatfac_extern = myoptd[RBW_NCATFD];

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
  vsfh = alloc_vecd(num_var + 1);
  matrix_categorical_vals_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern + num_reg_ordered_extern);

  vector_continuous_stddev = alloc_vecd(num_reg_continuous_extern);

  for(j = 0; j < num_reg_continuous_extern; j++)
    vector_continuous_stddev[j] = mysd[j];

  vector_continuous_stddev_extern = vector_continuous_stddev;

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

  bwm_penalty_mode = 0;
  bwm_penalty_value = DBL_MAX;
  if (penalty_mode[0] == 1) {
    double pmult = penalty_mult[0];
    double y_mean = 0.0;
    double mse0 = 0.0;
    if (pmult < 1.0) pmult = 1.0;
    for (i = 0; i < num_obs_train_extern; i++)
      y_mean += vector_Y_extern[i];
    y_mean /= (double) num_obs_train_extern;
    for (i = 0; i < num_obs_train_extern; i++) {
      const double dy = vector_Y_extern[i] - y_mean;
      mse0 += dy*dy;
    }
    mse0 /= (double) num_obs_train_extern;
    if (mse0 <= 0.0) mse0 = DBL_MIN;
    if (myopti[RBW_MI] == RBWM_CVAIC) {
      const double denom = 1.0 - 2.0/((double) num_obs_train_extern);
      if (denom > 0.0) {
        const double base_aic = log(mse0) + (1.0/denom);
        bwm_penalty_value = base_aic + log(pmult);
      }
    } else {
      bwm_penalty_value = pmult * mse0;
    }
    if (R_FINITE(bwm_penalty_value))
      bwm_penalty_mode = 1;
  }

  // initialize permutation arrays
  ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  // attempt tree build, if enabled 
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                 4*num_reg_continuous_extern, ipt, &kdt_extern_X);

    //put training data into tree-order using the index array

    for( j=0;j<num_reg_unordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_unordered_train_extern[j][i]=runo[j*num_obs_train_extern+ipt[i]];
    
    
    for( j=0;j<num_reg_ordered_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_ordered_train_extern[j][i]=rord[j*num_obs_train_extern+ipt[i]];

    for( j=0;j<num_reg_continuous_extern;j++)
      for( i=0;i<num_obs_train_extern;i++ )
        matrix_X_continuous_train_extern[j][i]=rcon[j*num_obs_train_extern+ipt[i]];

    /* response variable */
    for( i=0;i<num_obs_train_extern;i++ )
      vector_Y_extern[i] = y[ipt[i]];
    
    //boxSearch(kdt_extern, 0, tb, &nl);
  }

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


  /* Initialize scale factors and Directions for NR modules */

  initialize_nr_vector_scale_factor(BANDWIDTH_reg_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0,
                                    KERNEL_reg_unordered_extern,
                                    int_use_starting_values,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vector_scale_factor,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_vector_scale_factor(BANDWIDTH_reg_extern,
                                    0,                /* Not Random (0) Random (1) */
                                    int_RANDOM_SEED,
                                    int_LARGE_SF,
                                    num_obs_train_extern,
                                    0,
                                    0,
                                    0,
                                    num_reg_continuous_extern,
                                    num_reg_unordered_extern,
                                    num_reg_ordered_extern,
                                    0,
                                    KERNEL_reg_unordered_extern,
                                    0,
                                    scale_cat,
                                    pow((double)4.0/(double)3.0,0.2),             /* Init for continuous vars */
                                    nconfac_extern, ncatfac_extern,
                                    num_categories_extern,
                                    vector_continuous_stddev,
                                    vsfh,
                                    lbc_init, hbc_init, c_init, 
                                    lbd_init, hbd_init, d_init,
                                    matrix_X_continuous_train_extern,
                                    matrix_Y_continuous_train_extern);

  initialize_nr_directions(BANDWIDTH_reg_extern,
                           num_obs_train_extern,
                           num_reg_continuous_extern,
                           num_reg_unordered_extern,
                           num_reg_ordered_extern,
                           0,
                           0,
                           0,
                           vsfh,
                           num_categories_extern,
                           matrix_y,
                           0, int_RANDOM_SEED, 
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir,
                           matrix_X_continuous_train_extern,
                           matrix_Y_continuous_train_extern);


  /* When multistarting, set counter */

  iMs_counter = 0;

  /* assign the function to be optimized */
  switch(myopti[RBW_MI]){
  case RBWM_CVAIC : bwmfunc = cv_func_regression_categorical_aic_c; break;
  case RBWM_CVLS : bwmfunc = cv_func_regression_categorical_ls; break;
  default : REprintf("np.c: invalid bandwidth selection method.");
    error("np.c: invalid bandwidth selection method.");break;
  }

  spinner(0);

  bwmfunc_raw = bwmfunc;
  bwm_num_reg_continuous = num_reg_continuous_extern;
  bwm_num_reg_unordered = num_reg_unordered_extern;
  bwm_num_reg_ordered = num_reg_ordered_extern;
  bwm_kernel_unordered = KERNEL_reg_unordered_extern;
  bwm_num_categories = num_categories_extern;
  if (bwm_use_transform) {
    int n = bwm_num_reg_continuous + bwm_num_reg_unordered + bwm_num_reg_ordered;
    if (bwm_transform_buf_len < n + 1) {
      bwm_transform_buf = (double *) realloc(bwm_transform_buf, (n + 1) * sizeof(double));
      bwm_transform_buf_len = n + 1;
    }
  }
  bwm_reset_counters();

  fret_best = bwmfunc_wrapper(vector_scale_factor);
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
         bwmfunc_wrapper);


  if(int_RESTART_FROM_MIN == RE_MIN_TRUE){

    initialize_nr_directions(BANDWIDTH_reg_extern,
                             num_obs_train_extern,
                             num_reg_continuous_extern,
                             num_reg_unordered_extern,
                             num_reg_ordered_extern,
                             0,
                             0,
                             0,
                             vsfh,
                             num_categories_extern,
                             matrix_y,
                             0, int_RANDOM_SEED, 
                             lbc_dir, dfc_dir, c_dir, initc_dir,
                             lbd_dir, hbd_dir, d_dir, initd_dir,
                             matrix_X_continuous_train_extern,
                             matrix_Y_continuous_train_extern);


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
           bwmfunc_wrapper);

  }

  iImproved = (fret < fret_best);
  *timing = timing_extern;

  objective_function_values[0]=fret;
  objective_function_evals[0]=bwm_eval_count;
  objective_function_invalid[0]=bwm_invalid_count;
  /* When multistarting save initial minimum of objective function and scale factors */


  if(iMultistart == IMULTI_TRUE){
    fret_best = fret;
    vector_scale_factor_multistart = alloc_vecd(num_var + 1);

    for(i = 1; i <= num_var; i++)
      vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];

    /* Conduct search from new random values of the search parameters */

    for(imsnum = iMs_counter = 1; iMs_counter < iNum_Multistart; imsnum++,iMs_counter++){

      /* Initialize scale factors and directions for NR modules */
				
      initialize_nr_vector_scale_factor(BANDWIDTH_reg_extern,
                                        1,        /* Not Random (0) Random (1) */
                                        int_RANDOM_SEED,
                                        int_LARGE_SF,
                                        num_obs_train_extern,
                                        0,
                                        0,
                                        0,
                                        num_reg_continuous_extern,
                                        num_reg_unordered_extern,
                                        num_reg_ordered_extern,
                                        0,
                                        KERNEL_reg_unordered_extern,
                                        int_use_starting_values,
                                        scale_cat,
                                        pow((double)4.0/(double)3.0,0.2),     /* Init for continuous vars */
                                        nconfac_extern, ncatfac_extern,
                                        num_categories_extern,
                                        vector_continuous_stddev,
                                        vector_scale_factor,
                                        lbc_init, hbc_init, c_init, 
                                        lbd_init, hbd_init, d_init,
                                        matrix_X_continuous_train_extern,
                                        matrix_Y_continuous_train_extern);

      initialize_nr_directions(BANDWIDTH_reg_extern,
                               num_obs_train_extern,
                               num_reg_continuous_extern,
                               num_reg_unordered_extern,
                               num_reg_ordered_extern,
                               0,
                               0,
                               0,
                               vsfh,
                               num_categories_extern,
                               matrix_y,
                               1, int_RANDOM_SEED, 
                               lbc_dir, dfc_dir, c_dir, initc_dir,
                               lbd_dir, hbd_dir, d_dir, initd_dir,
                               matrix_X_continuous_train_extern,
                               matrix_Y_continuous_train_extern);


      /* Conduct direction set search */

      bwm_reset_counters();

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
             bwmfunc_wrapper);

      if(int_RESTART_FROM_MIN == RE_MIN_TRUE)	{
						
        initialize_nr_directions(BANDWIDTH_reg_extern,
                                 num_obs_train_extern,
                                 num_reg_continuous_extern,
                                 num_reg_unordered_extern,
                                 num_reg_ordered_extern,
                                 0,
                                 0,
                                 0,
                                 vsfh,
                                 num_categories_extern,
                                 matrix_y,
                                 0, int_RANDOM_SEED, 
                                 lbc_dir, dfc_dir, c_dir, initc_dir,
                                 lbd_dir, hbd_dir, d_dir, initd_dir,
                                 matrix_X_continuous_train_extern,
                                 matrix_Y_continuous_train_extern);

						
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
               bwmfunc_wrapper);

      }

      /* If this run resulted in an improved minimum save information */

      if(fret < fret_best){
        fret_best = fret;
        iImproved = iMs_counter+1;
        *timing = timing_extern;
        
        for(i = 1; i <= num_var; i++)	
          vector_scale_factor_multistart[i] = (double) vector_scale_factor[i];
      }
      objective_function_values[iMs_counter]=fret;
      objective_function_evals[iMs_counter]=bwm_eval_count;
      objective_function_invalid[iMs_counter]=bwm_invalid_count;

    }

    /* Save best for estimation */

    fret = fret_best;

    for(i = 1; i <= num_var; i++)
      vector_scale_factor[i] = (double) vector_scale_factor_multistart[i];

    free(vector_scale_factor_multistart);

  }

  if (bwm_use_transform)
    bwm_to_constrained(vector_scale_factor, num_var);

  /* return data to R */
  if (BANDWIDTH_reg_extern == BW_GEN_NN || 
      BANDWIDTH_reg_extern == BW_ADAP_NN){
    for( i=0; i<num_reg_continuous_extern; i++ )
      vector_scale_factor[i+1]=np_fround(vector_scale_factor[i+1]);
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
  safe_free(vsfh);
  safe_free(num_categories_extern);

  free_mat(matrix_categorical_vals_extern, num_reg_unordered_extern+num_reg_ordered_extern);

  free(vector_continuous_stddev);

  safe_free(ipt);
  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  if(int_MINIMIZE_IO != IO_MIN_TRUE)
    Rprintf("\r                   \r");

  //fprintf(stderr,"\nNP TOASTY\n");
  return ;
  
}


void np_regression(double * tuno, double * tord, double * tcon, double * ty,
                   double * euno, double * eord, double * econ, double * ey,
                   double * rbw, 
                   double * mcv, double * padnum, 
                   double * nconfac, double * ncatfac, double * mysd,
                   int * myopti, 
                   double * cm, double * cmerr, double * g, double *gerr, 
                   double * xtra){

  double * vector_scale_factor, * ecm = NULL, * ecmerr = NULL, ** eg = NULL, **egerr = NULL;
  double * lambda, ** matrix_bandwidth;
  double RS, MSE, MAE, MAPE, CORR, SIGN, pad_num;

  int i,j, num_var;
  int ey_is_ty, do_grad, train_is_eval, num_obs_eval_alloc, max_lev, old_reg;

  int * ipt = NULL, * ipe = NULL;  // point permutation, see tree.c
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
    error("\n(np_regression): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
  }

  KERNEL_reg_extern = myopti[REG_CKRNEVI];
  KERNEL_reg_unordered_extern = myopti[REG_UKRNEVI];
  KERNEL_reg_ordered_extern = myopti[REG_OKRNEVI];

  int_LARGE_SF = myopti[REG_LSFI];
  int_MINIMIZE_IO = myopti[REG_MINIOI];
  BANDWIDTH_reg_extern = myopti[REG_BWI];

  do_grad = myopti[REG_GRAD];
  int_ll_extern = myopti[REG_LL];
#ifdef MPI2
  {
    int ll_min = 0, ll_max = 0;
    MPI_Allreduce(&int_ll_extern, &ll_min, 1, MPI_INT, MPI_MIN, comm[1]);
    MPI_Allreduce(&int_ll_extern, &ll_max, 1, MPI_INT, MPI_MAX, comm[1]);
    if(ll_min != ll_max){
      if(my_rank == 0){
        REprintf("\n[npRmpi] Warning: inconsistent regression type across ranks (min=%d max=%d). Forcing all ranks to min.\n", ll_min, ll_max);
        R_FlushConsole();
      }
      int_ll_extern = ll_min;
    }
  }
#endif

  max_lev = myopti[REG_MLEVI];
  pad_num = *padnum;

  nconfac_extern = *nconfac;
  ncatfac_extern = *ncatfac;

  int_TREE_X = myopti[REG_DOTREEI];
  old_reg = myopti[REG_OLDREGI];

#ifdef MPI2
  num_obs_eval_alloc = MAX((int)ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
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

  vector_continuous_stddev_extern = mysd;
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

  if (!ey_is_ty)
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

  ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  if(!train_is_eval) {
    ipe = (int *)malloc(num_obs_eval_extern*sizeof(int));
    if(!(ipe != NULL))
      error("!(ipe != NULL)");

    for(i = 0; i < num_obs_eval_extern; i++){
      ipe[i] = i;
    }
  } else {
    ipe = ipt;
  }

  // attempt tree build, if enabled 
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    if((BANDWIDTH_reg_extern != BW_ADAP_NN) || ((BANDWIDTH_reg_extern == BW_ADAP_NN) && train_is_eval)){
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipt, &kdt_extern_X);

      //put training data into tree-order using the index array

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+ipt[i]];
    
    
      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+ipt[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+ipt[i]];

      /* response variable */
      for( i=0;i<num_obs_train_extern;i++ )
        vector_Y_extern[i] = ty[ipt[i]];

    } else {
      build_kdtree(matrix_X_continuous_eval_extern, num_obs_eval_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipe, &kdt_extern_X);

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+ipe[i]];

      if(!ey_is_ty)
        for(i=0;i<num_obs_eval_extern;i++)
          vector_Y_eval_extern[i] = ey[ipe[i]];

    }
  }


  /* Conduct estimation */
	
  /* 
     nb - KERNEL_(|un)ordered_den are set to zero upon declaration 
     - they have only one kernel type each at the moment 
  */

  if(old_reg){
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
                            0, // do not suppress_parallel
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
  } else {

    kernel_estimate_regression_categorical_tree_np(int_ll_extern,
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
                                                   vector_Y_extern,
                                                   vector_Y_eval_extern,
                                                   &vector_scale_factor[1],
                                                   num_categories_extern,
                                                   matrix_categorical_vals_extern,
                                                   ecm,
                                                   do_grad ? eg : NULL,
                                                   ecmerr,
                                                   do_grad ? egerr : NULL,
                                                   &RS,
                                                   &MSE,
                                                   &MAE,
                                                   &MAPE,
                                                   &CORR,
                                                   &SIGN);


  }

  for(i=0;i<num_obs_eval_extern;i++)
    cm[ipe[i]] = ecm[i];
      
  for(i=0;i<num_obs_eval_extern;i++)
    cmerr[ipe[i]] = ecmerr[i];


  if(do_grad){
    for(j=0;j<num_var;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        g[j*num_obs_eval_extern+ipe[i]]=eg[j][i];

    for(j=0;j<num_reg_continuous_extern;j++)
      for(i=0;i<num_obs_eval_extern;i++)
        gerr[j*num_obs_eval_extern+ipe[i]]=egerr[j][i];
  }


  /* write the return values */


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

  safe_free(ipt);

  if(!train_is_eval)
    safe_free(ipe);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
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
  vector_continuous_stddev_extern = NULL;

  safe_free(lambda);

  return;
}

void np_kernelsum(double * tuno, double * tord, double * tcon, 
                  double * ty, double * weights,
                  double * euno, double * eord, double * econ, 
                  double * bw,
                  double * mcv, double * padnum, 
                  int * operator,
                  int * myopti, double * kpow, 
                  double * weighted_sum, double * weighted_p_sum,
                  double * kernel_weights){

  int * ipt = NULL, * ipe = NULL;  // point permutation, see tree.c
      
  /* the ys are the weights */

  double * vector_scale_factor, * ksum, * p_ksum = NULL, pad_num, * kw = NULL;
  int i,j,k, num_var, num_obs_eval_alloc;
  int no_y, leave_one_out, train_is_eval, do_divide_bw;
  int max_lev, no_weights, sum_element_length, return_kernel_weights;
  int p_operator, do_score, do_ocg, p_nvar = 0;
  int ksum_is_output = 0, pksum_is_output = 0, kw_is_output = 0;

  struct th_table * otabs = NULL;
  struct th_entry * ret = NULL;
  int ** matrix_ordered_indices = NULL;

  int ncol_Y, ncol_W;

  int * kernel_c = NULL, * kernel_u = NULL, * kernel_o = NULL;

  int npks_err = 0;

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
  do_divide_bw = myopti[KWS_BDIVI];
  
  max_lev = myopti[KWS_MLEVI];
  pad_num = *padnum;

  /* the y and weight matrices will be contained in these variables */
  ncol_Y = myopti[KWS_YNCOLI];
  ncol_W = myopti[KWS_WNCOLI];

  int_TREE_X = myopti[KWS_DOTREEI];
  return_kernel_weights = myopti[KWS_RKWI];
  p_operator = myopti[KWS_POPI];
  do_score = myopti[KWS_PSCOREI];
  do_ocg = myopti[KWS_POCGI];

  nconfac_extern = ncatfac_extern = 0.0;

  no_y = (ncol_Y == 0);
  no_weights = (ncol_W == 0);

  sum_element_length = (no_y ? 1 : ncol_Y)*(no_weights ? 1 : ncol_W);

#ifdef MPI2
  num_obs_eval_alloc = MAX(ceil((double) num_obs_eval_extern / (double) iNum_Processors),1)*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval_extern;
#endif

  if(train_is_eval && (num_obs_eval_extern != num_obs_train_extern)){
    REprintf("\n(np_kernelsum): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
    error("\n(np_kernelsum): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
  }

  /* allocate */

  matrix_X_unordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_unordered_extern);
  matrix_X_ordered_train_extern = alloc_matd(num_obs_train_extern, num_reg_ordered_extern);
  matrix_X_continuous_train_extern = alloc_matd(num_obs_train_extern, num_reg_continuous_extern);
  
  /* for the moment we will just allocate a vector of ones */
  /* vector_Y_extern = (no_y)?NULL:alloc_vecd(num_obs_train_extern); */

  matrix_Y_continuous_train_extern = alloc_matd(num_obs_train_extern, ncol_Y);
  matrix_Y_ordered_train_extern = alloc_matd(num_obs_train_extern, ncol_W);

  num_categories_extern = alloc_vecu(num_reg_unordered_extern+num_reg_ordered_extern);
  matrix_categorical_vals_extern = alloc_matd(max_lev, num_reg_unordered_extern + num_reg_ordered_extern);

  vector_scale_factor = alloc_vecd(num_var + 1);
  if(!int_TREE_X && (num_obs_eval_alloc == num_obs_eval_extern)){
    ksum = weighted_sum;
    ksum_is_output = 1;
  } else {
    ksum = alloc_vecd(num_obs_eval_alloc*sum_element_length);
  }

  if((p_operator != OP_NOOP) || do_ocg){
    p_nvar = ((p_operator != OP_NOOP) ? num_reg_continuous_extern : 0) + ((do_score || do_ocg) ? num_reg_unordered_extern + num_reg_ordered_extern : 0);
    if(!int_TREE_X && (num_obs_eval_alloc == num_obs_eval_extern)){
      p_ksum = weighted_p_sum;
      pksum_is_output = 1;
    } else {
      p_ksum = alloc_vecd(num_obs_eval_alloc*sum_element_length*p_nvar);
    }
  }

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

  for( j = 0; j < ncol_Y; j++ )
    for( i = 0; i < num_obs_train_extern; i++ )
      matrix_Y_continuous_train_extern[j][i] = ty[j*num_obs_train_extern+i];

  for( j = 0; j < ncol_W; j++ )
    for( i = 0; i < num_obs_train_extern; i++ )
      matrix_Y_ordered_train_extern[j][i] = weights[j*num_obs_train_extern+i];

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

  ipt = (int *)malloc(num_obs_train_extern*sizeof(int));
  if(!(ipt != NULL))
    error("!(ipt != NULL)");

  for(i = 0; i < num_obs_train_extern; i++){
    ipt[i] = i;
  }

  if(!train_is_eval) {
    ipe = (int *)malloc(num_obs_eval_extern*sizeof(int));
    if(!(ipe != NULL))
      error("!(ipe != NULL)");

    for(i = 0; i < num_obs_eval_extern; i++){
      ipe[i] = i;
    }
  } else {
    ipe = ipt;
  }

  // attempt tree build, if enabled 
  int_TREE_X = int_TREE_X && ((num_reg_continuous_extern != 0) ? NP_TREE_TRUE : NP_TREE_FALSE);

  if(int_TREE_X == NP_TREE_TRUE){
    if((BANDWIDTH_reg_extern != BW_ADAP_NN) || ((BANDWIDTH_reg_extern == BW_ADAP_NN) && train_is_eval)){
      build_kdtree(matrix_X_continuous_train_extern, num_obs_train_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipt, &kdt_extern_X);

      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_unordered_train_extern[j][i]=tuno[j*num_obs_train_extern+ipt[i]];
    
    
      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_ordered_train_extern[j][i]=tord[j*num_obs_train_extern+ipt[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_train_extern;i++ )
          matrix_X_continuous_train_extern[j][i]=tcon[j*num_obs_train_extern+ipt[i]];

      for( j = 0; j < ncol_Y; j++ )
        for( i = 0; i < num_obs_train_extern; i++ )
          matrix_Y_continuous_train_extern[j][i] = ty[j*num_obs_train_extern+ipt[i]];

      for( j = 0; j < ncol_W; j++ )
        for( i = 0; i < num_obs_train_extern; i++ )
          matrix_Y_ordered_train_extern[j][i] = weights[j*num_obs_train_extern+ipt[i]];

    } else {
      build_kdtree(matrix_X_continuous_eval_extern, num_obs_eval_extern, num_reg_continuous_extern, 
                   4*num_reg_continuous_extern, ipe, &kdt_extern_X);


      for( j=0;j<num_reg_unordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_unordered_eval_extern[j][i]=euno[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_ordered_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_ordered_eval_extern[j][i]=eord[j*num_obs_eval_extern+ipe[i]];

      for( j=0;j<num_reg_continuous_extern;j++)
        for( i=0;i<num_obs_eval_extern;i++ )
          matrix_X_continuous_eval_extern[j][i]=econ[j*num_obs_eval_extern+ipe[i]];

    }
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

  if(do_ocg && (num_reg_ordered_extern > 0)){
    otabs = (struct th_table *)malloc(num_reg_ordered_extern*sizeof(struct th_table));
    matrix_ordered_indices = (int **)malloc(num_reg_ordered_extern*sizeof(int *));
    int * tc = (int *)malloc(num_reg_ordered_extern*num_obs_eval_extern*sizeof(int));
    for(i = 0; i < num_reg_ordered_extern; i++)
      matrix_ordered_indices[i] = tc + i*num_obs_eval_extern;
  }

  for(j=num_reg_unordered_extern, k = 0; j < (num_reg_unordered_extern+num_reg_ordered_extern); j++, k++){
    i = 0;

    do { 
      matrix_categorical_vals_extern[j][i] = mcv[j*max_lev+i];
    } while(++i < max_lev && mcv[j*max_lev+i] != pad_num);
    num_categories_extern[j] = i;

    if(do_ocg){
      if(thcreate_r((size_t)ceil(1.6*num_categories_extern[j]), otabs + k) == TH_ERROR)
        error("hash table creation failed");

      for(i = 0; i < num_categories_extern[j]; i++){
        struct th_entry centry;
        centry.key.dkey = matrix_categorical_vals_extern[j][i];
        centry.data = i;

        if(thsearch_r(&centry, TH_ENTER, &ret, otabs+k) == TH_FAILURE)
          error("insertion into hash table failed");
      }
      
      // now do lookups
      struct th_entry te;
      te.key.dkey = pad_num;
      te.data = -1;

      ret = &te;
      for(i = 0; i < num_obs_eval_extern; i++){
        if(ret->key.dkey != matrix_X_ordered_eval_extern[k][i]){
          te.key.dkey = matrix_X_ordered_eval_extern[k][i];
          if(thsearch_r(&te, TH_SEARCH, &ret, otabs+k) == TH_FAILURE)
            error("hash table lookup failed (which should be impossible)");

        } 

        matrix_ordered_indices[k][i] = ret->data;
      }
    }
  }

  if(return_kernel_weights){
    if(!int_TREE_X && (BANDWIDTH_reg_extern != BW_ADAP_NN) && (num_obs_eval_alloc == num_obs_eval_extern)){
      kw = kernel_weights;
      kw_is_output = 1;
    } else {
      kw = alloc_vecd(num_obs_train_extern*num_obs_eval_extern);
    }
  }
  
  kernel_c = (int *)malloc(sizeof(int)*num_reg_continuous_extern);

  for(i = 0; i < num_reg_continuous_extern; i++)
    kernel_c[i] = KERNEL_reg_extern;

  kernel_u = (int *)malloc(sizeof(int)*num_reg_unordered_extern);

  for(i = 0; i < num_reg_unordered_extern; i++)
    kernel_u[i] = KERNEL_reg_unordered_extern;

  kernel_o = (int *)malloc(sizeof(int)*num_reg_ordered_extern);

  for(i = 0; i < num_reg_ordered_extern; i++)
    kernel_o[i] = KERNEL_reg_ordered_extern;

  int *bpso = NULL;
  {
    const int bpso_len = num_reg_continuous_extern + num_reg_unordered_extern + num_reg_ordered_extern;
    if(bpso_len > 0){
      const int do_perm = (p_operator != OP_NOOP);
      const int doscoreocg = (do_score || do_ocg);

      bpso = (int *)R_alloc((size_t)bpso_len, sizeof(int));
      for(i = 0; i < num_reg_continuous_extern; i++)
        bpso[i] = do_perm;
      for(i = num_reg_continuous_extern; i < bpso_len; i++)
        bpso[i] = doscoreocg;
    }
  }
  
  
  if((npks_err=kernel_weighted_sum_np(kernel_c,
                                      kernel_u,
                                      kernel_o,
                                      BANDWIDTH_reg_extern,
                                      num_obs_train_extern,
                                      num_obs_eval_extern,
                                      num_reg_unordered_extern,
                                      num_reg_ordered_extern,
                                      num_reg_continuous_extern,
                                      leave_one_out,
                                      0,
                                      (int)(*kpow),
                                      do_divide_bw,
                                      0, 
                                      0, //not symmetric
                                      0, //disable 'twisting'
                                      0, // do not drop train
                                      0, // do not drop train
                                      operator,
                                      p_operator,
                                      do_score,
                                      do_ocg, // no ocg (for now)
                                      bpso,
                                      0, // don't explicity suppress parallel
                                      ncol_Y,
                                      ncol_W,
                                      int_TREE_X,
                                      0,
                                      kdt_extern_X,
                                      NULL, NULL, NULL,
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
                                      0,NULL,NULL,NULL,
                                      num_categories_extern,
                                      matrix_categorical_vals_extern,
                                      matrix_ordered_indices,
                                      ksum,
                                      p_ksum,
                                      kw)) == 1){
    Rprintf("kernel_weighted_sum_np has reported an error, probably due to invalid bandwidths\n");
  }


  if(!npks_err){
    if(int_TREE_X || !ksum_is_output){
      for(j = 0; j < num_obs_eval_extern; j++)
        for(i = 0; i < sum_element_length; i++)
          weighted_sum[ipe[j]*sum_element_length + i] = ksum[j*sum_element_length+i];
    }

    if(return_kernel_weights){
      if(!kw_is_output){
        if(BANDWIDTH_reg_extern != BW_ADAP_NN){ // adaptive weights are currently returned transposed...
          for(j = 0; j < num_obs_eval_extern; j++)
            for(i = 0; i < num_obs_train_extern; i++)
              kernel_weights[ipe[j]*num_obs_train_extern + ipt[i]] = kw[j*num_obs_train_extern + i];
        } else {
          for(j = 0; j < num_obs_train_extern; j++)
            for(i = 0; i < num_obs_eval_extern; i++)
              kernel_weights[ipe[i]*num_obs_train_extern + ipt[j]] = kw[j*num_obs_eval_extern + i];      
        }
      }
    }

    if(p_nvar > 0){
      if(int_TREE_X || !pksum_is_output){
        for(k = 0; k < p_nvar; k++){
          const int kidx = k*num_obs_eval_extern*sum_element_length;
          for(j = 0; j < num_obs_eval_extern; j++)
            for(i = 0; i < sum_element_length; i++)
              weighted_p_sum[kidx + ipe[j]*sum_element_length + i] = p_ksum[kidx + j*sum_element_length + i];
        }
      }
    }

  }
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
  free_mat(matrix_Y_continuous_train_extern, ncol_Y);
  free_mat(matrix_Y_ordered_train_extern, ncol_W);

  safe_free(num_categories_extern);
  safe_free(vector_scale_factor);
  if(!ksum_is_output)
    safe_free(ksum);

  if(!kw_is_output)
    safe_free(kw);

  safe_free(ipt);

  free(kernel_c);
  free(kernel_u);
  free(kernel_o);

  if(!train_is_eval)
    safe_free(ipe);

  if(int_TREE_X == NP_TREE_TRUE){
    free_kdtree(&kdt_extern_X);
    int_TREE_X = NP_TREE_FALSE;
  }

  if(p_nvar > 0){
    if(!pksum_is_output)
      safe_free(p_ksum);
  }

  if(do_ocg && (num_reg_ordered_extern > 0)){

    for(k = 0; k < num_reg_ordered_extern; k++){
      thdestroy_r(otabs+k);
    }
    free(otabs);
    free(matrix_ordered_indices[0]);
    free(matrix_ordered_indices);
  }

  return;
}


void np_quantile_conditional(double * tc_con,
                             double * tu_uno, double * tu_ord, double * tu_con,
                             double * eu_uno, double * eu_ord, double * eu_con,
                             double * quantile,
                             double * mybw, 
                             double * mcv, double *padnum,
                             double * nconfac, double * ncatfac, double * mysd,
                             int * myopti, double * myoptd,
                             double * yq, double * yqerr, double *yg){
  /* Likelihood bandwidth selection for density estimation */

  double **g = NULL, * eq, * eqerr;
  double ftol, small, tol;
  double pad_num;

  int i,j, max_lev;
  int num_var, num_obs_eval_alloc;
  int num_all_var, num_var_var, train_is_eval, do_gradients;
  int itmax;

  int iNum_Multistart;

  double lbc_dir, c_dir;
  double initc_dir;
  double lbd_dir, hbd_dir, d_dir, initd_dir;
  int dfc_dir;

  iNum_Multistart = myopti[CQ_NMULTII];
  imsnum = 0;
  imstot = myopti[CQ_NMULTII]; /* iNum_Multistart */

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
    error("\n(np_quantile_conditional): consistency check failed, train_is_eval but num_obs_train_extern != num_obs_eval_extern. bailing\n");
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

  dfc_dir = myopti[CQ_DFC_DIRI];
  lbc_dir = myoptd[CQ_LBC_DIRD];
  c_dir = myoptd[CQ_C_DIRD];
  initc_dir = myoptd[CQ_INITC_DIRD]; 

  lbd_dir = myoptd[CQ_LBD_DIRD]; 
  hbd_dir = myoptd[CQ_HBD_DIRD]; 
  d_dir = myoptd[CQ_D_DIRD]; 
  initd_dir = myoptd[CQ_INITD_DIRD]; 

  gamma_extern = *quantile;
  nconfac_extern = *nconfac;
  ncatfac_extern = *ncatfac;
  vector_continuous_stddev_extern = mysd;

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
                           iNum_Multistart,            /* Maximum number of multistarts */
                           1.0e-10,         /* Zero for all intents and purposes */
                           lbc_dir, dfc_dir, c_dir, initc_dir,
                           lbd_dir, hbd_dir, d_dir, initd_dir);

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
