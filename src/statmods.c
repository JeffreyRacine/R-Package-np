#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <R.h>

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

/* This function generates the index for bootstrap samples of size n */

extern int int_LARGE_SF;
extern int int_DEBUG;
extern int int_VERBOSE;
extern int int_ROBUST;
extern int int_nn_k_min_extern;
extern double nconfac_extern;
extern double *vector_continuous_stddev_extern;


#include <math.h>

/*
  Preserve the existing Numerical Recipes sort behavior without forming
  one-element-before-start pointers at call sites.
*/
static void sort_safe(int n, double *values)
{
    double *work;
    int i;

    if ((values == NULL) || (n <= 1))
        return;

    work = alloc_vecd(n + 1);
    for (i = 0; i < n; i++)
        work[i + 1] = values[i];

    sort(n, work);

    for (i = 0; i < n; i++)
        values[i] = work[i + 1];

    free(work);
}

int np_fround(double x)
{
    double intpart, fracpart;
    int origint;

    fracpart = modf(x, &intpart);
    origint = (int) intpart;
    if (fracpart < 0.5)                           /* small fraction, round down */
        return origint;
    else if (fracpart > 0.5)                      /* large, round up */
        return origint + 1;
    else if (origint % 2)                         /* halfway case, odd number */
        return origint + 1;
    else                                          /* halfway case, even number */
        return origint;
}


#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>


int simple_unique(int n, double * vector){
  int i, m;

  /* gcc 11 Found the following significant warnings: statmods.c:73:3:
     warning: 'sort' accessing 8 bytes in a region of size 0
     [-Wstringop-overflow=] */
  
  /*  double * v=NULL; */
  double *v;

  /* v=(double *)malloc(sizeof(double)*n); */

  v = alloc_vecd(n);

  for(i=0; i<n; i++)
    v[i]=vector[i];

  sort_safe(n, v);

  for(i=0, m=1; i < (n - 1); i++)
    m += v[i]!=v[i+1];
  free(v);

  return(m);
}

/* 7/24/95: Added pointer arithmetic for efficiency */

/* This will generate the mean of a vector */

double meand(int n, double *vector)
{

    int i;
    double *pi;
    double *pi_temp;
    double sum = 0.0;
    double stat;

    if(int_ROBUST == 1)
    {

        double *vector_temp;
        int int_med, int_medl, int_medh;

/* Robust employs median */
/* First create temporary vector and sort... */

        vector_temp = alloc_vecd(n);
        pi = &vector[0];
        pi_temp = &vector_temp[0];

        for (i=0; i<n; i++)
        {
            *pi_temp++ = *pi++;
        }

/* Sort... */

        sort_safe(n, vector_temp);

        int_med = np_fround(((double)n-1.0)/2.0);
        int_medl = np_fround(((double)n-2.0)/2.0);
        int_medh = np_fround(((double)n)/2.0);

/*stat = (n % 2 ? odd number of obs : even number of obs);*/

        stat = (n % 2 ? vector_temp[int_med] :
        0.5*(vector_temp[int_medl]+vector_temp[int_medh]));

        free(vector_temp);

    }
    else
    {

        pi = &vector[0];

        for(i=0; i < n; i++)
        {
            sum += (double) *pi++;
        }

        stat = sum / (double) n;

    }

    return((double) stat);

}


double standerrd(int n, double *vector)
{

    int i;

    double *pi;
    double *p;
    double sum = 0.0;
    double sumsq = 0.0;
    double temp;
    double temp1;
    double std = 0.0;
    double IQR = 0.0;

/* November 18, 2008, using adaptive measure of spread */

    double *vector_temp;
    int int_25, int_25l, int_25h, int_75, int_75l, int_75h;

/* First create temporary vector and sort... */

    vector_temp = alloc_vecd(n);

    pi = &vector_temp[0];
    p = &vector[0];
    
    for (i=0; i<n; i++)
      {
        *pi++ = *p++;
      }

/* Sort... */

    sort_safe(n, vector_temp);

/* Interquartile Range */

    int_25 = np_fround(0.25*((double)n+1.0)-1);
    int_25l = np_fround(0.25*((double)n)-1);
    int_25h = np_fround(0.25*((double)n));

    int_75 = np_fround(0.75*((double)n+1.0)-1);
    int_75l = np_fround(0.75*((double)n)-1);
    int_75h = np_fround(0.75*((double)n));

/* (n % 2 ? odd number of obs : even number of obs) */

    IQR = (n % 2 ? (vector_temp[int_75]-vector_temp[int_25])
           : (0.25*vector_temp[int_75l]+0.75*vector_temp[int_75h])
           - (0.75*vector_temp[int_25l]+0.25*vector_temp[int_25h]));

    free(vector_temp);

    pi = &vector[0];
    
    for(i=0; i < n; i++)
      {
        sum += (double) (temp = (double) *pi++);
        sumsq += (double) ipow(temp, 2);
      }

    /* Jan 15 2009, no longer using ml estimate of variance */
    
    temp1 = (sumsq - ipow(sum,2)/(double)n)/((double)(n-1));

    if(temp1 > 0.0)
      {
        std = (double) sqrt(temp1);           /* Variance */
      }
    else
      {
        if(int_VERBOSE == 1)
          {
#ifdef MPI2
            if(my_rank == 0)
              {
#endif
                REprintf("\nFunction standerrd(): invalid standard error estimate (%lg)", temp1);
                REprintf("\nsum = %lg, sumsq = %lg, n = %d", (double) sum, (double) sumsq, n);
                REprintf("\nVar 1");
#ifdef MPI2
              }
#endif
          }
        return((double)0.0);
      }

    if(IQR > 0) {
      /* Return min(std,IQR/1.348980) */
      if(std < IQR/1.348980)
	{
	  return((double)std);
	} 
      else 
	{
	  return((double)IQR/1.348980);
	}

    } else {
      /* April 19 2009 - pathological case encountered when IQR is
	 zero but std > 0. Rather than bombing out of R (we return
	 min(std,IQR/1.348), in this case we return std. */
      return((double)std);
    }

}


double max_unordered_bw(int num_categories,
                        int kernel){
  return((kernel == UKERNEL_UAA) ? (num_categories - 1.0)/num_categories : 1.0);
}

int is_valid_unordered_bw(double lambda,
                          int num_categories,
                          int kernel){
  return ((lambda >= 0.0) && (lambda <= max_unordered_bw(num_categories, kernel)));
}

static int build_sorted_unique_support(int n,
                                       const double *vector_data,
                                       double **support_out,
                                       int **support_count_out,
                                       int *support_n_out)
{
  int i, m;
  double *support;
  int *support_count;

  support = alloc_vecd(n);
  support_count = NULL;

  if (support_count_out != NULL)
    support_count = alloc_vecu(n);

  for (i = 0; i < n; i++)
    support[i] = vector_data[i];

  sort_safe(n, support);

  m = 0;
  for (i = 0; i < n; i++) {
    if ((m == 0) || (support[i] != support[m - 1])) {
      support[m] = support[i];
      if (support_count != NULL)
        support_count[m] = 1;
      m++;
    } else if (support_count != NULL) {
      support_count[m - 1]++;
    }
  }

  *support_out = support;
  if (support_count_out != NULL)
    *support_count_out = support_count;
  *support_n_out = m;

  return 0;
}

static int find_support_index(int support_n,
                              const double *support,
                              double value)
{
  int lo, hi;

  lo = 0;
  hi = support_n - 1;

  while (lo <= hi) {
    const int mid = lo + (hi - lo) / 2;

    if (support[mid] == value)
      return mid;

    if (support[mid] < value) {
      lo = mid + 1;
    } else {
      hi = mid - 1;
    }
  }

  return -1;
}

static int nearest_positive_radius_from_support(int support_n,
                                                const double *support,
                                                int center_idx,
                                                double *radius_out)
{
  const double dleft =
    (center_idx > 0) ? (support[center_idx] - support[center_idx - 1]) : DBL_MAX;
  const double dright =
    (center_idx + 1 < support_n) ? (support[center_idx + 1] - support[center_idx]) : DBL_MAX;
  const double distance = (dleft < dright) ? dleft : dright;

  if (distance <= DBL_MIN)
    return 1;

  *radius_out = distance;
  return 0;
}

static int kth_observation_radius_from_support(int support_n,
                                               const double *support,
                                               const int *support_count,
                                               int center_idx,
                                               int int_k_nn,
                                               double *radius_out)
{
  int left, right, count;

  left = center_idx - 1;
  right = center_idx + 1;
  count = support_count[center_idx] - 1;

  if (count >= int_k_nn)
    return nearest_positive_radius_from_support(support_n, support, center_idx, radius_out);

  while ((left >= 0) || (right < support_n)) {
    double dleft = DBL_MAX;
    double dright = DBL_MAX;
    double distance;

    if (left >= 0)
      dleft = support[center_idx] - support[left];
    if (right < support_n)
      dright = support[right] - support[center_idx];

    distance = (dleft < dright) ? dleft : dright;

    if (distance <= DBL_MIN)
      return 1;

    if ((left >= 0) && ((support[center_idx] - support[left]) == distance)) {
      count += support_count[left];
      left--;
    }

    if ((right < support_n) && ((support[right] - support[center_idx]) == distance)) {
      count += support_count[right];
      right++;
    }

    if (count >= int_k_nn) {
      *radius_out = distance;
      return 0;
    }
  }

  return 1;
}

static int compute_nn_distance_observation_support_subset(int num_obs,
                                                          double *vector_data,
                                                          int int_k_nn,
                                                          int query_start,
                                                          int query_end,
                                                          double *nn_distance)
{
  int i, j, support_n;
  double *support;
  int *support_count;
  double *support_radius;

  support = NULL;
  support_count = NULL;
  support_radius = NULL;

  if ((query_start < 0) || (query_end >= num_obs) || (query_start > query_end))
    return 1;

  if (build_sorted_unique_support(num_obs, vector_data, &support, &support_count, &support_n) != 0)
    return 1;

  if ((int_k_nn < 1) || (int_k_nn > num_obs - 1)) {
    free(support);
    if (support_count != NULL)
      free(support_count);
    return 1;
  }

  support_radius = alloc_vecd(support_n);

  for (i = 0; i < support_n; i++) {
    if (kth_observation_radius_from_support(
          support_n,
          support,
          support_count,
          i,
          int_k_nn,
          &support_radius[i]
        ) != 0) {
      free(support);
      free(support_count);
      free(support_radius);
      return 1;
    }
  }

  for (i = query_start, j = 0; i <= query_end; i++, j++) {
    const int idx = find_support_index(support_n, support, vector_data[i]);

    if (idx < 0) {
      free(support);
      free(support_count);
      free(support_radius);
      return 1;
    }

    nn_distance[j] = support_radius[idx];
    if (nn_distance[j] <= DBL_MIN) {
      free(support);
      free(support_count);
      free(support_radius);
      return 1;
    }
  }

  free(support);
  free(support_count);
  free(support_radius);

  return 0;
}

static int lower_bound_support(int support_n,
                               const double *support,
                               double value)
{
  int lo, hi;

  lo = 0;
  hi = support_n;

  while (lo < hi) {
    const int mid = lo + (hi - lo) / 2;

    if (support[mid] < value) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }

  return lo;
}

static int kth_observation_radius_for_eval_from_support(int support_n,
                                                        const double *support,
                                                        const int *support_count,
                                                        double eval_value,
                                                        int int_k_nn,
                                                        double *radius_out)
{
  int left, right, count;
  int insert_idx;
  int exact_match;

  insert_idx = lower_bound_support(support_n, support, eval_value);
  exact_match = ((insert_idx < support_n) && (support[insert_idx] == eval_value));

  if (exact_match) {
    left = insert_idx - 1;
    right = insert_idx + 1;
    count = support_count[insert_idx];
  } else {
    left = insert_idx - 1;
    right = insert_idx;
    count = 0;
  }

  if (count >= int_k_nn) {
    if (exact_match)
      return nearest_positive_radius_from_support(support_n, support, insert_idx, radius_out);
    return 1;
  }

  while ((left >= 0) || (right < support_n)) {
    double dleft = DBL_MAX;
    double dright = DBL_MAX;
    double distance;

    if (left >= 0)
      dleft = eval_value - support[left];
    if (right < support_n)
      dright = support[right] - eval_value;

    distance = (dleft < dright) ? dleft : dright;

    if (distance <= DBL_MIN)
      return 1;

    if ((left >= 0) && ((eval_value - support[left]) == distance)) {
      count += support_count[left];
      left--;
    }
    if ((right < support_n) && ((support[right] - eval_value) == distance)) {
      count += support_count[right];
      right++;
    }

    if (count >= int_k_nn) {
      *radius_out = distance;
      return 0;
    }
  }

  return 1;
}

static int compute_nn_distance_train_eval_observation_support_subset(int num_obs_train,
                                                                     int num_obs_eval,
                                                                     double *vector_data_train,
                                                                     double *vector_data_eval,
                                                                     int int_k_nn,
                                                                     int query_start,
                                                                     int query_end,
                                                                     double *nn_distance)
{
  int i, j, support_n;
  double *support;
  int *support_count;

  support = NULL;
  support_count = NULL;

  if ((query_start < 0) || (query_end >= num_obs_eval) || (query_start > query_end))
    return 1;

  if (build_sorted_unique_support(num_obs_train, vector_data_train, &support, &support_count, &support_n) != 0)
    return 1;

  if ((int_k_nn < 1) || (int_k_nn > num_obs_train - 1)) {
    free(support);
    free(support_count);
    return 1;
  }

  for (i = query_start, j = 0; i <= query_end; i++, j++) {
    if (kth_observation_radius_for_eval_from_support(
          support_n, support, support_count, vector_data_eval[i], int_k_nn, &nn_distance[j]
        ) != 0) {
      free(support);
      free(support_count);
      return 1;
    }
  }

  free(support);
  free(support_count);
  return 0;
}

/* Population variance, double precision */
/* Returns 0 upon success, 1 upon failure (constant most likely) */


int compute_nn_distance(int num_obs, int suppress_parallel, double *vector_data,
int int_k_nn, double *nn_distance)
{

#ifdef MPI2
    int stride = (int)ceil((double) num_obs / (double) iNum_Processors);
		int return_flag = 0;
		int return_flag_MPI = 0;
    if(stride < 1) stride = 1;

    int is, ie;
#endif

#ifndef MPI2

    if((int_k_nn < 1)||(int_k_nn > num_obs-1))
    {
			if(int_VERBOSE == 1)
        {
					REprintf("\n** Error: Invalid Kth nearest neighbor (%d).", int_k_nn);
        }
			return(1);
    }

    return(compute_nn_distance_observation_support_subset(num_obs, vector_data, int_k_nn, 0, num_obs - 1, nn_distance));

#endif

#ifdef MPI2

    if((int_k_nn < 1)||(int_k_nn > num_obs-1))
    {
        if(int_VERBOSE == 1)
        {
            if(my_rank == 0)
            {
                REprintf("\n** Error: Invalid Kth nearest neighbor (%d).", int_k_nn);
            }
        }
        return(1);
    }


    if(!suppress_parallel){
      is = my_rank*stride;
      ie = MIN(num_obs,(my_rank+1)*stride) - 1;

    } else {
      is = 0;
      ie = num_obs - 1;
    }

    if (compute_nn_distance_observation_support_subset(num_obs, vector_data, int_k_nn, is, ie, nn_distance) != 0)
        return_flag_MPI = 1;

    if(!suppress_parallel){
      MPI_Reduce(&return_flag_MPI, &return_flag, 1, MPI_INT, MPI_SUM, 0, comm[1]);
      MPI_Bcast(&return_flag, 1, MPI_INT, 0, comm[1]);
    }

		if(return_flag > 0) {
			return(1);
		}

    if(!suppress_parallel){
      if(my_rank == 0){
        MPI_Gather(MPI_IN_PLACE, stride, MPI_DOUBLE, nn_distance, stride, MPI_DOUBLE, 0, comm[1]);
      } else {
        MPI_Gather(nn_distance, stride, MPI_DOUBLE, NULL, stride, MPI_DOUBLE, 0, comm[1]);
      }
      MPI_Bcast(nn_distance, num_obs, MPI_DOUBLE, 0, comm[1]);
    }
#endif

    return(0);

}


int compute_nn_distance_train_eval(int num_obs_train,
                                   int num_obs_eval,
                                   int suppress_parallel,
                                   double *vector_data_train,
                                   double *vector_data_eval,
                                   int int_k_nn,
                                   double *nn_distance){

#ifdef MPI2
    int stride = (int)ceil((double) num_obs_eval / (double) iNum_Processors);
		int return_flag = 0;
		int return_flag_MPI = 0;
    if(stride < 1) stride = 1;
    int is, ie;
#endif

    if((int_k_nn < 1)||(int_k_nn > num_obs_train-1))
    {
        if(int_VERBOSE == 1)
        {
#ifdef MPI2
            if(my_rank == 0)
            {
#endif
                REprintf("\n** Error: Invalid Kth nearest neighbor (%d).", int_k_nn);
#ifdef MPI2
            }
#endif

        }
        return(1);
    }

#ifndef MPI2
    return(compute_nn_distance_train_eval_observation_support_subset(
      num_obs_train,
      num_obs_eval,
      vector_data_train,
      vector_data_eval,
      int_k_nn,
      0,
      num_obs_eval - 1,
      nn_distance
    ));

#endif

#ifdef MPI2

    if(!suppress_parallel){
      is = my_rank*stride;
      ie = MIN(num_obs_eval,(my_rank+1)*stride) - 1;
    } else {
      is = 0;
      ie = num_obs_eval - 1;
    }

    if (compute_nn_distance_train_eval_observation_support_subset(
          num_obs_train,
          num_obs_eval,
          vector_data_train,
          vector_data_eval,
          int_k_nn,
          is,
          ie,
          nn_distance
        ) != 0)
      return_flag_MPI = 1;

    if(!suppress_parallel){
      MPI_Reduce(&return_flag_MPI, &return_flag, 1, MPI_INT, MPI_SUM, 0, comm[1]);
      MPI_Bcast(&return_flag, 1, MPI_INT, 0, comm[1]);
    }
		if(return_flag > 0) {
			return(1);
		}

    if(!suppress_parallel){
      if(my_rank == 0){
        MPI_Gather(MPI_IN_PLACE, stride, MPI_DOUBLE, nn_distance, stride, MPI_DOUBLE, 0, comm[1]);
      } else {
        MPI_Gather(nn_distance, stride, MPI_DOUBLE, NULL, stride, MPI_DOUBLE, 0, comm[1]);
      }
      MPI_Bcast(nn_distance, num_obs_eval, MPI_DOUBLE, 0, comm[1]);
    }
#endif

    return(0);

}



int initialize_nr_directions(int BANDWIDTH, 
                             int num_obs,
                             int num_reg_continuous, 
                             int num_reg_unordered, 
                             int num_reg_ordered, 
                             int num_var_continuous, 
                             int num_var_unordered, int num_var_ordered, 
                             double * vector_scale_factor, 
                             int * num_categories, 
                             double **matrix_y, 
                             int random, 
                             int seed, 
                             double lbc_dir, 
                             int dfc_dir, 
                             double c_dir,
                             double initc_dir,
                             double lbd_dir, 
                             double hbd_dir, 
                             double d_dir,
                             double initd_dir,
                             double ** matrix_x_continuous,
                             double ** matrix_y_continuous){

  int i, j, li;
  // sfac ought to be smaller than lbd_dir, 1-hbd_dir in
  // initialize_nr_vector_scale_factor() and sfac constant here in
  // order to keep from going out of bounds (so e.g. lbd_dir=.2,hbd_dir=.8,
  // sfac constant .25 would work, or .3,.7,.375, etc.). Note the
  // 3-sqrt(5) is related to golden section search, golden ratios, and
  // golden means

  //const double sfac = 0.25*(3.0-sqrt(5)); 
  //const double csfac = 2.5*(3.0-sqrt(5));

  li =  num_reg_continuous + num_reg_unordered + num_reg_ordered + 
    num_var_continuous + num_var_unordered + num_var_ordered;

  for(i = 1; i <= li; i++)
    for(j = 1; j <= li; j++)
      matrix_y[j][i] = (j == i)? 1.0 : 0.0;

  if(vector_scale_factor == NULL) return(0);

  // nvc + nrc
  // this is only to ensure that initial cv function probes don't
  // go outside of the allowed ranges for bws

  li =  num_reg_continuous + num_var_continuous;

  if(BANDWIDTH==BW_FIXED){
    for(i = 1; i <= li; i++)
      matrix_y[i][i] = vector_scale_factor[i]*(random ? chidev(&seed, dfc_dir)  + lbc_dir: initc_dir)*c_dir;
  }else{
    for(i = 1; i <= num_reg_continuous; i++){
      const double bw_max =
        ((BANDWIDTH == BW_ADAP_NN) || (BANDWIDTH == BW_GEN_NN)) ?
        (double)(num_obs - 1) :
        (double)(simple_unique(num_obs,matrix_x_continuous[i-1]) - 1);
      matrix_y[i][i] = ceil(MIN(vector_scale_factor[i], bw_max - vector_scale_factor[i])*(random ? ran3(&seed): 1.0));
    }
    for(i = num_reg_continuous+1; i <= li; i++){
      const double bw_max =
        ((BANDWIDTH == BW_ADAP_NN) || (BANDWIDTH == BW_GEN_NN)) ?
        (double)(num_obs - 1) :
        (double)(simple_unique(num_obs,matrix_y_continuous[i-num_reg_continuous-1]) - 1);
      matrix_y[i][i] = ceil(MIN(vector_scale_factor[i], bw_max - vector_scale_factor[i])*(random ? ran3(&seed): 1.0));
    }
  }
  if(num_categories == NULL) return(0);

  // nvu
  li = num_reg_continuous + num_var_continuous;
  
  for(i = li + 1, j = 0; i <= (li + num_var_unordered); i++, j++) 
    matrix_y[i][i] = MIN(vector_scale_factor[i], 1.0 - vector_scale_factor[i])*(random ? (hbd_dir-lbd_dir)*ran3(&seed) + lbd_dir: initd_dir)*d_dir;

  // nvo
  li += num_var_unordered;

  for(; i <= (li + num_var_ordered); i++) 
    matrix_y[i][i] = MIN(vector_scale_factor[i], (1.0 - vector_scale_factor[i]))*(random ? (hbd_dir-lbd_dir)*ran3(&seed) + lbd_dir: initd_dir)*d_dir;

  //nru
  j += num_var_ordered;
  li += num_var_ordered;

  for(; i <= (li + num_reg_unordered); i++, j++)
    matrix_y[i][i] = MIN(vector_scale_factor[i], 1.0 - vector_scale_factor[i])*(random ? (hbd_dir-lbd_dir)*ran3(&seed) + lbd_dir: initd_dir)*d_dir;

  // nro
  li += num_reg_unordered;

  for(; i <= (li + num_reg_ordered); i++)
    matrix_y[i][i] = MIN(vector_scale_factor[i], (1.0 - vector_scale_factor[i]))*(random ? (hbd_dir-lbd_dir)*ran3(&seed) + lbd_dir: initd_dir)*d_dir;

  return(0);

}

void initialize_nr_vector_scale_factor(int BANDWIDTH,
                                       int RANDOM,
                                       int seed,
                                       int int_large,
                                       int num_obs,
                                       int num_var_continuous,
                                       int num_var_unordered,
                                       int num_var_ordered,
                                       int num_reg_continuous,
                                       int num_reg_unordered,
                                       int num_reg_ordered,
                                       int kernel_yu,
                                       int kernel_xu,
                                       int int_use_starting_values,
                                       int scale_cat,
                                       double init_continuous,
                                       double nconfac, 
                                       double ncatfac,
                                       int *num_categories,
                                       double *vector_continuous_stddev,
                                       double *vector_scale_factor,  
                                       double lbc_init,
                                       double hbc_init,
                                       double c_init,
                                       double lbd_init,
                                       double hbd_init,
                                       double d_init,
                                       double ** matrix_x_continuous,
                                       double ** matrix_y_continuous){
  int i, l = 0;

  // lbc and hbc and init_continuous [fed in] play a similar role to
  // lbd, jbd, and initd - provide a range of sfs and point of
  // reference for first (always non-random) attempt. If calling
  // program uses init_continuous of 0.25, this starts at
  // 0.25*EssDee()*n^(-1/(2*p+q)) for e.g. regression. If you then
  // want to vary from quite small bws to larger ones, you want hbc to
  // be, say, 2 hence restarting we get both small and large
  // bws. init_continuous=0.5 and hbc=2.0 means search will start from
  // initial bws that are 50% smaller than the previous defaults of
  // 1.06 and (implicitly) 1.0. Also, by now setting a lower bound > 0
  // (e.g. 0.1) we can avoid any start with impossibly small starting
  // values.


  const int fixed_bw = (BANDWIDTH == BW_FIXED);
  const int count_bw = ((BANDWIDTH == BW_ADAP_NN) || (BANDWIDTH == BW_GEN_NN));
  double bw_nf = 0;
  const double bw_cmin = fixed_bw ? 0.0 : (double)MAX(1, int_nn_k_min_extern);
  const double bw_cmax = fixed_bw ? DBL_MAX : num_obs-1;
  const int ncon = num_reg_continuous + num_var_continuous;

  int_use_starting_values = int_use_starting_values && (!RANDOM);

  // x continuous
  for(i = 0; i < num_reg_continuous; i++,l++){
    if(!fixed_bw){
      bw_nf = MAX((double)MAX(1, int_nn_k_min_extern),
                  ceil(sqrt(count_bw ? num_obs : simple_unique(num_obs,matrix_x_continuous[i]))));
    }
    const double bwi = fixed_bw ? (int_large ? vector_continuous_stddev[l] * nconfac : 1.0) : bw_nf;

    if(!int_use_starting_values){
      if(RANDOM){
        if(fixed_bw){
          vector_scale_factor[l+1] = bwi*((hbc_init-lbc_init)*ran3(&seed)+lbc_init);
        } else {
          vector_scale_factor[l+1] = ceil(bwi*((hbc_init-lbc_init)*ran3(&seed)+lbc_init)+1.0);
        }
      } else {
        if(fixed_bw){
          vector_scale_factor[l+1] = bwi*c_init;
        } else {
          vector_scale_factor[l+1] = ceil(bwi*c_init);
        }
      }
    } else {
      if(fixed_bw) {
        if((vector_scale_factor[l+1] < bw_cmin) || (vector_scale_factor[l+1] > bw_cmax)){
          REprintf("\n** Warning: invalid sf in init_nr_sf() [%g]\n", vector_scale_factor[l+1]);
          vector_scale_factor[l+1] = bwi*c_init;
        }
      } else {
        const double bw_kmax =
          count_bw ? (double)(num_obs - 1) : (double)(simple_unique(num_obs,matrix_x_continuous[i]) - 1);
        if((vector_scale_factor[l+1] < bw_cmin) || (vector_scale_factor[l+1] > bw_kmax)){
          REprintf("\n** Warning: invalid sf in init_nr_sf() [%g]\n", vector_scale_factor[l+1]);
          vector_scale_factor[l+1] = ceil(bwi*c_init);
        }
      }
    }
  }

  // y continuous
  for(i = 0; i < num_var_continuous; i++,l++){
    if(!fixed_bw){
      bw_nf = MAX((double)MAX(1, int_nn_k_min_extern),
                  ceil(sqrt(count_bw ? num_obs : simple_unique(num_obs,matrix_y_continuous[i]))));
    }
    const double bwi = fixed_bw ? (int_large ? vector_continuous_stddev[l] * nconfac : 1.0) : bw_nf;

    if(!int_use_starting_values){
      if(RANDOM){
        if(fixed_bw){
          vector_scale_factor[l+1] = bwi*((hbc_init-lbc_init)*ran3(&seed)+lbc_init);
        } else {
          vector_scale_factor[l+1] = ceil(bwi*((hbc_init-lbc_init)*ran3(&seed)+lbc_init)+1.0);
        }
      } else {
        if(fixed_bw){
          vector_scale_factor[l+1] = bwi*c_init;
        } else {
          vector_scale_factor[l+1] = ceil(bwi*c_init);
        }
      }
    } else {
      if(fixed_bw) {
        if((vector_scale_factor[l+1] < bw_cmin) || (vector_scale_factor[l+1] > bw_cmax)){
          REprintf("\n** Warning: invalid sf in init_nr_sf() [%g]\n", vector_scale_factor[l+1]);
          vector_scale_factor[l+1] = bwi*c_init;
        }
      } else {
        const double bw_kmax =
          count_bw ? (double)(num_obs - 1) : (double)(simple_unique(num_obs,matrix_y_continuous[i]) - 1);
        if((vector_scale_factor[l+1] < bw_cmin) || (vector_scale_factor[l+1] > bw_kmax)){
          REprintf("\n** Warning: invalid sf in init_nr_sf() [%g]\n", vector_scale_factor[l+1]);
          vector_scale_factor[l+1] = ceil(bwi*c_init);
        }
      }
    }
  }

  for(i = 0; i < num_var_unordered; i++,l++){
    const double bwi = (scale_cat ? 1.0 : 1.0/ncatfac)*(int_large ? ncatfac : 1.0)*max_unordered_bw(num_categories[l-ncon], kernel_yu);

    if(!int_use_starting_values){
      vector_scale_factor[l+1] = bwi*(RANDOM ? (hbd_init-lbd_init)*ran3(&seed)+lbd_init : d_init);
    } else {
      if(!is_valid_unordered_bw(vector_scale_factor[l+1], num_categories[l-ncon], kernel_yu)){
        REprintf("\n** Warning: invalid sf in init_nr_sf() [%g]\n", vector_scale_factor[l+1]);
        vector_scale_factor[l+1] = bwi;
      }
    }
  }

  for(i = 0; i < num_var_ordered; i++,l++){
    const double bwi = (scale_cat ? 1.0 : 1.0/ncatfac) * (int_large ? ncatfac : 1.0);

    if(!int_use_starting_values){
      vector_scale_factor[l+1] = bwi*(RANDOM ? (hbd_init-lbd_init)*ran3(&seed)+lbd_init : d_init);
    } else {
      if((vector_scale_factor[l+1] < 0.0) || (vector_scale_factor[l+1] > 1.0)){
        REprintf("\n** Warning: invalid sf in init_nr_sf() [%g]\n", vector_scale_factor[l+1]);
        vector_scale_factor[l+1] = bwi;
      }
    }
  }

  for(i = 0; i < num_reg_unordered; i++,l++){
    const double bwi = (scale_cat ? 1.0 : 1.0/ncatfac) * (int_large ? ncatfac : 1.0)*max_unordered_bw(num_categories[l-ncon], kernel_xu);

    if(!int_use_starting_values){
      vector_scale_factor[l+1] = bwi*(RANDOM ? (hbd_init-lbd_init)*ran3(&seed)+lbd_init : d_init);
    } else {
      if(!is_valid_unordered_bw(vector_scale_factor[l+1], num_categories[l-ncon], kernel_xu)){
        REprintf("\n** Warning: invalid sf in init_nr_sf() [%g]\n", vector_scale_factor[l+1]);
        vector_scale_factor[l+1] = bwi;
      }
    }
  }

  for(i = 0; i < num_reg_ordered; i++,l++){
    const double bwi = (scale_cat ? 1.0 : 1.0/ncatfac) * (int_large ? ncatfac : 1.0);

    if(!int_use_starting_values){
      vector_scale_factor[l+1] = bwi*(RANDOM ? (hbd_init-lbd_init)*ran3(&seed)+lbd_init : d_init);
    } else {
      if((vector_scale_factor[l+1] < 0.0) || (vector_scale_factor[l+1] > 1.0)){
        REprintf("\n** Warning: invalid sf in init_nr_sf() [%g]\n", vector_scale_factor[l+1]);
        vector_scale_factor[l+1] = bwi;
      }
    }
  }

}

double fGoodness_of_Fit(int iNum_Obs, double *fvector_Y, double *fkernel_fit)
{
    int i;
    double dRSQ;
    double dMean = 0.0;
    double dCov_Sum = 0.0;
    double dvector_Y_Sum_SQ = 0.0;
    double dkernel_fit_Sum_SQ = 0.0;

    for(i=0;i < iNum_Obs; i++)
    {
        dMean += fvector_Y[i];
    }

    dMean /= (double) iNum_Obs;

    for(i=0;i < iNum_Obs; i++)
    {
        dCov_Sum += ((double) fvector_Y[i] - dMean)*((double) fkernel_fit[i] - dMean);
        dvector_Y_Sum_SQ += ipow(fvector_Y[i] - dMean,2);
        dkernel_fit_Sum_SQ += ipow(fkernel_fit[i] - dMean,2);
    }

    if((dvector_Y_Sum_SQ !=0.0)&&(dkernel_fit_Sum_SQ !=0.0))
    {
        dRSQ = (dCov_Sum*dCov_Sum)/(dvector_Y_Sum_SQ*dkernel_fit_Sum_SQ);
    }
    else
    {
        dRSQ = 0.0;
    }

    return((double) dRSQ);

}


double fMSE(int iNum_Obs, double *fvector_Y, double *fkernel_fit)
{

/* Mean square error */

    int i;
    double sum = 0.0;

    for(i=0;i < iNum_Obs; i++)
    {
        sum += (fvector_Y[i]-fkernel_fit[i])*(fvector_Y[i]-fkernel_fit[i]);
    }

    return((double) sum / (double) iNum_Obs);

}


double fCORR(int iNum_Obs, double *fvector_Y, double *fkernel_fit)
{

/* Pearson's correlation coefficient */

    int i;
    double sum_act_fit = 0.0;
    double sum_actsq = 0.0;
    double sum_fitsq = 0.0;

    double mean_act = meand(iNum_Obs, fvector_Y);
    double mean_fit = meand(iNum_Obs, fkernel_fit);

    for(i=0;i < iNum_Obs; i++)
    {
        sum_act_fit += (fvector_Y[i]-mean_act)*(fkernel_fit[i]-mean_fit);
        sum_actsq += ipow(fvector_Y[i]-mean_act,2);
        sum_fitsq += ipow(fkernel_fit[i]-mean_fit,2);
    }

    if((sum_actsq !=0.0)&&(sum_fitsq !=0.0))
    {
        return((double) (sum_act_fit/(sqrt(sum_actsq)*sqrt(sum_fitsq))));
    }
    else
    {
        return((double)0.0);
    }

}


double fMAE(int iNum_Obs, double *fvector_Y, double *fkernel_fit)
{

/* Mean absolute error */

    int i;
    double sum = 0.0;

    for(i=0;i < iNum_Obs; i++)
    {
        sum += fabs(fvector_Y[i]-fkernel_fit[i]);
    }

    return((double) sum / (double) iNum_Obs);

}


double fMAPE(int iNum_Obs, double *fvector_Y, double *fkernel_fit)
{

/* Mean absolute percentage error */

    int i;
    double sum = 0.0;

    for(i=0;i < iNum_Obs; i++)
    {
        if(fvector_Y[i] != 0.0)
        {
            sum += fabs((fvector_Y[i]-fkernel_fit[i])/fvector_Y[i]);
        }
        else
        {
            sum += fabs((fvector_Y[i] - fkernel_fit[i])/(fvector_Y[i] + fkernel_fit[i])/2.0);
        }
    }

    return((double) sum / (double) iNum_Obs);

}


double fSIGN(int iNum_Obs, double *fvector_Y, double *fkernel_fit)
{

/* Mean correct sign */

    int i;
    double sum = 0.0;

    for(i=0;i < iNum_Obs; i++)
    {
        if((fvector_Y[i]*fkernel_fit[i]) >= 0.0)
        {
            sum +=1;
        }
    }

    return((double) sum / (double) iNum_Obs);

}


int determine_categorical_vals(
int num_obs,
int num_var_unordered,
int num_var_ordered,
int num_reg_unordered,
int num_reg_ordered,
double **matrix_Y_unordered,
double **matrix_Y_ordered,
double **matrix_X_unordered,
double **matrix_X_ordered,
int *num_categories,
double **matrix_categorical_vals)
{

/* We treat variables in the following order - Y_un, Y_or, X_un, X_or */

    int i;
    int k;
    int l;

    double **matrix;

    FILE *File10 = NULL;

    if((num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered) == 0)
    {
        return(0);
    }

#ifdef MPI2
    if(my_rank == 0)
    {
#endif
        if(int_DEBUG == 1)
        {
            File10 = fopen("cat_dat.dbg", "w");
        }
#ifdef MPI2
    }
#endif

/* Variables (conditioned variables in density estimation) */

    matrix=alloc_matd(num_obs,num_var_unordered);

/* Copy data to matrix to be sorted */

    for(i=0; i < num_obs;i++)
    {
        for(k = 0; k < num_var_unordered; k++)
        {
            matrix[k][i]=matrix_Y_unordered[k][i];
        }
    }

    for(k = 0; k < num_var_unordered; k++)
    {
        sort_safe(num_obs, matrix[k]);
        matrix_categorical_vals[k][0]=matrix[k][0];
        for(i=1, l=1; i < num_obs;i++)
        {
            if(matrix[k][i]!=matrix[k][i-1]) matrix_categorical_vals[k][l++]=matrix[k][i];
        }

        num_categories[k] = l;

        if(int_VERBOSE == 1)
        {

            if(num_categories[k] == num_obs)
            {
#ifdef MPI2
                if(my_rank == 0)
                {
#endif
                    REprintf("\n** Note: unordered variable %d contains strictly unique values\n** [%d out of %d are unique]", (int) k+1, (int) num_categories[k], (int) num_obs);
#ifdef MPI2
                }
#endif

            }

        }

#ifdef MPI2
        if(my_rank == 0)
        {
#endif
            if((int_DEBUG == 1) && (File10 != NULL))
            {

                fprintf(File10, "\nThere are %d unique values for unordered variable %d.", l, k+1);

                for(i=0; i < l;i++)
                {
                    fprintf(File10, "\nValue %d unique for unordered variable %d is %g", i+1, k+1, matrix_categorical_vals[k][i]);
                }

            }
#ifdef MPI2
        }
#endif

    }

    free_mat(matrix,num_var_unordered);

    matrix=alloc_matd(num_obs,num_var_ordered);

/* Copy data to matrix to be sorted */

    for(i=0; i < num_obs;i++)
    {
        for(k = 0; k < num_var_ordered; k++)
        {
            matrix[k][i]=matrix_Y_ordered[k][i];
        }
    }

    for(k = 0; k < num_var_ordered; k++)
    {
        sort_safe(num_obs, matrix[k]);
        matrix_categorical_vals[k+num_var_unordered][0]=matrix[k][0];
        for(i=1, l=1; i < num_obs;i++)
        {
            if(matrix[k][i]!=matrix[k][i-1]) matrix_categorical_vals[k+num_var_unordered][l++]=matrix[k][i];
        }

        num_categories[k+num_var_unordered] = l;

        if(int_VERBOSE == 1)
        {

            if(num_categories[k+num_var_unordered] == num_obs)
            {
#ifdef MPI2
                if(my_rank == 0)
                {
#endif
                    REprintf("\n** Note: ordered variable %d contains strictly unique values\n** [%d out of %d are unique]", (int) k+num_var_ordered+1, (int) num_categories[k+num_var_unordered], (int) num_obs);
#ifdef MPI2
                }
#endif

            }

        }

#ifdef MPI2
        if(my_rank == 0)
        {
#endif
            if((int_DEBUG == 1) && (File10 != NULL))
            {

                fprintf(File10, "\nThere are %d unique values for ordered variable %d.", l, k+1);

                for(i=0; i < l;i++)
                {
                    fprintf(File10, "\nValue %d unique for ordered variable %d is %g", i+1, k+1, matrix_categorical_vals[k+num_var_unordered][i]);
                }

            }
#ifdef MPI2
        }
#endif

    }

    free_mat(matrix,num_var_ordered);

/* Regressors */

    matrix=alloc_matd(num_obs,num_reg_unordered);

/* Copy data to matrix to be sorted */

    for(i=0; i < num_obs;i++)
    {
        for(k = 0; k < num_reg_unordered; k++)
        {
            matrix[k][i]=matrix_X_unordered[k][i];
        }
    }

    for(k = 0; k < num_reg_unordered; k++)
    {
        sort_safe(num_obs, matrix[k]);
        matrix_categorical_vals[k+num_var_unordered+num_var_ordered][0]=matrix[k][0];
        for(i=1, l=1; i < num_obs;i++)
        {
            if(matrix[k][i]!=matrix[k][i-1]) matrix_categorical_vals[k+num_var_unordered+num_var_ordered][l++]=matrix[k][i];
        }

        num_categories[k+num_var_unordered+num_var_ordered] = l;

        if(int_VERBOSE == 1)
        {

            if(num_categories[k+num_var_unordered+num_var_ordered] == num_obs)
            {
#ifdef MPI2
                if(my_rank == 0)
                {
#endif
                    REprintf("\n** Note: unordered predictor %d contains strictly unique values\n** [%d out of %d are unique]", (int) k+1, (int) num_categories[k+num_var_unordered+num_var_ordered], (int) num_obs);
#ifdef MPI2
                }
#endif

            }

        }

#ifdef MPI2
        if(my_rank == 0)
        {
#endif
            if((int_DEBUG == 1) && (File10 != NULL))
            {

                fprintf(File10, "\nThere are %d unique values for unordered predictor %d.", l, k+1);

                for(i=0; i < l;i++)
                {
                    fprintf(File10, "\nValue %d for unordered predictor %d is %g", i+1, k+1, matrix_categorical_vals[k+num_var_unordered+num_var_ordered][i]);
                }

            }
#ifdef MPI2
        }
#endif

    }

    free_mat(matrix,num_reg_unordered);

    matrix=alloc_matd(num_obs,num_reg_ordered);

/* Copy data to matrix to be sorted */

    for(i=0; i < num_obs;i++)
    {
        for(k = 0; k < num_reg_ordered; k++)
        {
            matrix[k][i]=matrix_X_ordered[k][i];
        }
    }

    for(k = 0; k < num_reg_ordered; k++)
    {
        sort_safe(num_obs, matrix[k]);
        matrix_categorical_vals[k+num_var_unordered+num_var_ordered+num_reg_unordered][0]=matrix[k][0];
        for(i=1, l=1; i < num_obs;i++)
        {
            if(matrix[k][i]!=matrix[k][i-1]) matrix_categorical_vals[k+num_var_unordered+num_var_ordered+num_reg_unordered][l++]=matrix[k][i];
        }

        num_categories[k+num_var_unordered+num_var_ordered+num_reg_unordered] = l;

        if(int_VERBOSE == 1)
        {

            if(num_categories[k+num_var_unordered+num_var_ordered+num_reg_unordered] == num_obs)
            {
#ifdef MPI2
                if(my_rank == 0)
                {
#endif
                    REprintf("\n** Note: ordered predictor %d contains strictly unique values\n** [%d out of %d are unique]", (int) k+1, (int) num_categories[k+num_var_unordered+num_var_ordered+num_reg_unordered], (int) num_obs);
#ifdef MPI2
                }
#endif

            }

        }

#ifdef MPI2
        if(my_rank == 0)
        {
#endif
            if((int_DEBUG == 1) && (File10 != NULL))
            {

                fprintf(File10, "\nThere are %d unique values for ordered predictor %d.", l, k+1);

                for(i=0; i < l;i++)
                {
                    fprintf(File10, "\nValue %d for ordered predictor %d is %g", i+1, k+1, matrix_categorical_vals[k+num_var_unordered+num_var_ordered+num_reg_unordered][i]);
                }

            }
#ifdef MPI2
        }
#endif

    }

    if(int_VERBOSE == 1)
    {
#ifdef MPI2
        if(my_rank == 0)
        {
#endif
            REprintf("\n");
#ifdef MPI2
        }
#endif
    }

#ifdef MPI2
    if(my_rank == 0)
    {
#endif
        if((int_DEBUG == 1) && (File10 != NULL))
        {
            fclose(File10);
        }
#ifdef MPI2
    }
#endif

    free_mat(matrix,num_reg_ordered);

    return(0);

}



/* For resampling conditional data under the null - `scramble' variables in test */


/* For resampling conditional data under the null - `scramble'
     (permute) variables in test, and also first resample the x & y */

int check_valid_scale_factor_cv(
int KERNEL,
int KERNEL_unordered_liracine,
int BANDWIDTH,
int BANDWIDTH_den_ml,
int REGRESSION_ML,
int num_obs,
int num_var_continuous,
int num_var_unordered,
int num_var_ordered,
int num_reg_continuous,
int num_reg_unordered,
int num_reg_ordered,
int *num_categories,
double *vector_scale_factor)
{

  /* April 30, 2008 - adding test for the LI_RACINE unordered kernel */

/* Dec 19, 2000 - changed around - now counter reflects variable order rather */
/* than in the variable index itself - can only improve speed */

    int i;
    int num_obs_m_1 = num_obs - 1;

/* In vector_scale_factor, order is continuous reg, continuous var, */
/* unordered variables, ordered variables, unordered regressors, ordered regressors */

    double temp_pow = DBL_MAX;

/* Set appropriate constants for categorical variable scaling factor O(h^2) */

    switch(KERNEL)
    {

        case 0:

/* Gaussian Kernel */

            temp_pow = 1.0/pow((double)num_obs, (2.0/(4.0 + (double) num_reg_continuous + num_var_continuous)));

            break;

        case 1:

/* Fourth Order Gaussian Kernel */

            temp_pow = 1.0/pow((double)num_obs, (2.0/(8.0 + (double) num_reg_continuous + num_var_continuous)));

            break;

        case 2:

/* Sixth Order Gaussian Kernel */

            temp_pow = 1.0/pow((double)num_obs, (2.0/(12.0 + (double) num_reg_continuous + num_var_continuous)));

            break;

        case 3:

/* Eighth Order Gaussian Kernel */

            temp_pow = 1.0/pow((double)num_obs, (2.0/(16.0 + (double) num_reg_continuous + num_var_continuous)));

            break;

        case 4:

/* Second Order Epanechnikov Kernel */

            temp_pow = 1.0/pow((double)num_obs, (2.0/(4.0 + (double) num_reg_continuous + num_var_continuous)));

            break;

        case 5:

/* Fourth Order Epanechnikov Kernel */

            temp_pow = 1.0/pow((double)num_obs, (2.0/(8.0 + (double) num_reg_continuous + num_var_continuous)));

            break;

        case 6:

/* Sixth Order Epanechnikov Kernel */

            temp_pow = 1.0/pow((double)num_obs, (2.0/(12.0 + (double) num_reg_continuous + num_var_continuous)));

            break;

        case 7:

/* Eighth Order Epanechnikov Kernel */

            temp_pow = 1.0/pow((double)num_obs, (2.0/(16.0 + (double) num_reg_continuous + num_var_continuous)));

            break;

        case 8:

/* Rectangular kernel - using second order Epanechnikov for now */

/* Second Order Epanechnikov Kernel */

            temp_pow = 1.0/pow((double)num_obs, (2.0/(4.0 + (double) num_reg_continuous + num_var_continuous)));

            break;

    }

/* Continuous Regressors */

    for(i = 1; i <= num_reg_continuous; i++)
    {
/* Check for admissible values */
        if(BANDWIDTH == 0)
        {
            if( (vector_scale_factor[i] <= 0.0) || (vector_scale_factor[i] > DBL_MAX) )
            {
                return(1);
            }
        }
        else if((BANDWIDTH == 1)||(BANDWIDTH == 2))
        {
            if( (np_fround(vector_scale_factor[i]) < MAX(1, int_nn_k_min_extern)) || (np_fround(vector_scale_factor[i]) > num_obs_m_1) )
            {
                return(1);
            }
        }
    }

/* Continuous Variables */

    for(i = num_reg_continuous+1; i <= num_reg_continuous+num_var_continuous; i++)
    {
/* Check for admissible values */
        if(BANDWIDTH == 0)
        {
            if( (vector_scale_factor[i] <= 0.0) || (vector_scale_factor[i] > DBL_MAX) )
            {
                return(1);
            }
        }
        else if((BANDWIDTH == 1)||(BANDWIDTH == 2))
        {
            if( (np_fround(vector_scale_factor[i]) < MAX(1, int_nn_k_min_extern)) || (np_fround(vector_scale_factor[i]) > num_obs_m_1 ) )
            {
                return(1);
            }
        }
    }

/* Unordered categorical Variables */

    for(i = num_var_continuous+num_reg_continuous+1; i <= num_var_continuous+num_reg_continuous+num_var_unordered; i++)
    {
/* Check for admissible value of lambda - in inadmissible, set to admissible */
      if(KERNEL_unordered_liracine == 1) {
        /* If using the unordered li_racine kernel use ordered bounds */
        if(int_LARGE_SF == 0)
        {
            if( (vector_scale_factor[i]*temp_pow > 1.0) || (vector_scale_factor[i] < 0.0) )
            {
                return(1);
            }
        }
        else
        {
            if( (vector_scale_factor[i] > 1.0) || (vector_scale_factor[i] < 0.0) )
            {
                return(1);
            }
        }
      } else {
        if(int_LARGE_SF == 0)
        {
            if( (vector_scale_factor[i]*temp_pow > (1.0 - 1.0/(double) num_categories[i-num_var_continuous-num_reg_continuous-1])) ||
                (vector_scale_factor[i] < 0.0))
            {
                return(1);
            }
        }
        else
        {
            if( (vector_scale_factor[i] > (1.0 - 1.0/(double) num_categories[i-num_var_continuous-num_reg_continuous-1])) ||
                (vector_scale_factor[i] < 0.0))
            {
                return(1);
            }
        }
      }
    }

/* Ordered categorical Variables */

    for(i = num_var_continuous+num_reg_continuous+num_var_unordered+1; i <= num_var_continuous+num_reg_continuous+num_var_unordered+num_var_ordered; i++)
    {
/* Check for admissible value of lambda - in inadmissible, set to admissible */
        if(int_LARGE_SF == 0)
        {
            if( (vector_scale_factor[i]*temp_pow > 1.0) || (vector_scale_factor[i] < 0.0) )
            {
                return(1);
            }
        }
        else
        {
            if( (vector_scale_factor[i] > 1.0) || (vector_scale_factor[i] < 0.0) )
            {
                return(1);
            }
        }
    }

/* Unordered categorical regressors */

    for(i = num_var_continuous+num_reg_continuous+num_var_unordered+num_var_ordered+1; i <= num_var_continuous+num_reg_continuous+num_var_unordered+num_var_ordered+num_reg_unordered; i++)
    {
/* Check for admissible value of lambda - in inadmissible, set to admissible */
      if(KERNEL_unordered_liracine == 1) {
        /* If using the unordered li_racine kernel use ordered bounds */
        if(int_LARGE_SF == 0)
        {
            if( (vector_scale_factor[i]*temp_pow > 1.0) || (vector_scale_factor[i] < 0.0) )
            {
                return(1);
            }
        }
        else
        {
            if( (vector_scale_factor[i] > 1.0) || (vector_scale_factor[i] < 0.0) )
            {
                return(1);
            }
        }
      } else {
        if(int_LARGE_SF == 0)
        {
            if( (vector_scale_factor[i]*temp_pow > (1.0 - 1.0/(double) num_categories[i-num_var_continuous-num_reg_continuous-1])) ||
                (vector_scale_factor[i] < 0.0))
            {
                return(1);
            }
        }
        else
        {
            if( (vector_scale_factor[i]> (1.0 - 1.0/(double) num_categories[i-num_var_continuous-num_reg_continuous-1])) ||
                (vector_scale_factor[i] < 0.0))
            {
                return(1);
            }
        }
      }
    }

/* Ordered categorical regressors */

    for(i = num_var_continuous+num_reg_continuous+num_var_unordered+num_var_ordered+num_reg_unordered+1; i <= num_var_continuous+num_reg_continuous+num_var_unordered+num_var_ordered+num_reg_unordered+num_reg_ordered; i++)
    {
/* Check for admissible value of lambda - in inadmissible, set to admissible */
        if(int_LARGE_SF == 0)
        {
            if( (vector_scale_factor[i]*temp_pow > 1.0) || (vector_scale_factor[i] < 0.0) )
            {
                return(1);
            }
        }
        else
        {
            if( (vector_scale_factor[i] > 1.0) || (vector_scale_factor[i] < 0.0) )
            {
                return(1);
            }
        }
    }

/* Density SF for regression */

    if(REGRESSION_ML == 1)
    {
        if(BANDWIDTH_den_ml == 0)
        {

/* Check for admissible value of density scale factor */
            if( (vector_scale_factor[num_var_continuous+num_reg_continuous+num_var_unordered+num_var_ordered+num_reg_unordered+num_reg_ordered+1] <= 0.0) ||
                (vector_scale_factor[num_var_continuous+num_reg_continuous+num_var_unordered+num_var_ordered+num_reg_unordered+num_reg_ordered+1] > DBL_MAX))
            {
                return(1);
            }

        }
        else if((BANDWIDTH_den_ml == 1)||(BANDWIDTH_den_ml == 2))
        {
            if( (np_fround(vector_scale_factor[num_var_continuous+num_reg_continuous+num_var_unordered+num_var_ordered+num_reg_unordered+num_reg_ordered+1]) < 1) ||
                (np_fround(vector_scale_factor[num_var_continuous+num_reg_continuous+num_var_unordered+num_var_ordered+num_reg_unordered+num_reg_ordered+1]) > (int) (num_obs-1) ))
            {
                return(1);
            }
        }

    }

    return(0);

}


int compute_continuous_stddev(
int int_large,
int num_obs,
int num_var_continuous,
int num_reg_continuous,
double **matrix_Y_continuous,
double **matrix_X_continuous,
double *vector_continuous_stddev)
{

    int i;

/* In vector_scale_factor, order is continuous reg, continuous var, */

/* Only compute if not using standardized values */

    if(int_large == 1)
    {

        for(i=0; i < num_reg_continuous; i++)
        {

            vector_continuous_stddev[i] = standerrd(num_obs, matrix_X_continuous[i]);

            if(vector_continuous_stddev[i] <= DBL_MIN)
            {
              error("\r ** Fatal Error in routine kernel_bandwidth() ** variable %d appears to be constant!", i);
            }

        }

        for(i=0; i < num_var_continuous; i++)
        {

            vector_continuous_stddev[i+num_reg_continuous] = standerrd(num_obs, matrix_Y_continuous[i]);

            if(vector_continuous_stddev[i+num_reg_continuous] <= DBL_MIN)
            {
              error("\r ** Fatal Error in routine kernel_bandwidth() ** variable %d appears to be constant!", i+num_reg_continuous);
            }
            
        }

    }

    return(0);

}





/* +++Date last modified: 05-Jul-1997 */

/*
 **  IPOW.C - Raise a number to an integer power
 **
 **  public domain by Mark Stephen with suggestions by Keiichi Nakasato
 */

double ipow(double x, int n)                      /* return x^n */
{
    double t = 1.0;

    if (!n)
        return 1.0;                               /* At the top. 0**0 = 1 */
    if (n < 0)
    {
        n = -n;
        x = 1.0/x;                                /* error if x == 0. Good                        */
    }                                             /* ZTC/SC returns inf, which is even better     */
    if (x == 0.0)
        return 0.0;
    do
    {
        if (n & 1)
            t *= x;
        n /= 2;                                   /* KN prefers if (n/=2) x*=x; This avoids an    */
        x *= x;                                   /* unnecessary but benign multiplication on     */
    } while (n);                                  /* the last pass, but the comparison is always  */
    return t;                                     /* true _except_ on the last pass.              */
}



int unique(
int num_obs,
double *x)
{

    int i;
    int unique = num_obs;
    double *dist;
    double x_max;

    dist = alloc_vecd(num_obs);

/* Obtain maximum in an unsorted list */

    x_max = x[0];

    for(i=1; i < num_obs; i++)
    {
        if(x[i] > x_max) x_max = x[i];
    }

/* Compute maximum number of unique distances (could be minimum?) */

    for(i=0; i < num_obs; i++)
    {
        dist[i] = fabs(x[i]-x_max);
    }

/* Sort distances and determine the number of unique ones */

    sort_safe(num_obs, dist);

    for(i=1; i < num_obs; i++)
    {
        if(dist[i]==dist[i-1]) unique--;
    }

		free(dist);

    return(unique);

}
