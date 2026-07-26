#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "headers.h"
#include "jksum_gaussian_density.h"

#if defined(NP_USE_ACCELERATE_GAUSS) && defined(__APPLE__) && defined(__arm64__)
#define NP_GAUSSIAN_DENSITY_PAIR_COMPILED 1
typedef long np_vDSP_Stride;
typedef unsigned long np_vDSP_Length;
extern void vDSP_vsmsaD(const double *, np_vDSP_Stride, const double *,
                        const double *, double *, np_vDSP_Stride,
                        np_vDSP_Length);
extern void vDSP_vsqD(const double *, np_vDSP_Stride, double *,
                      np_vDSP_Stride, np_vDSP_Length);
extern void vDSP_vsmulD(const double *, np_vDSP_Stride, const double *,
                        double *, np_vDSP_Stride, np_vDSP_Length);
extern void vDSP_vaddD(const double *, np_vDSP_Stride, const double *,
                       np_vDSP_Stride, double *, np_vDSP_Stride,
                       np_vDSP_Length);
extern void vDSP_sveD(const double *, np_vDSP_Stride, double *,
                      np_vDSP_Length);
extern void vvexp(double *, const double *, const int *);
#else
#define NP_GAUSSIAN_DENSITY_PAIR_COMPILED 0
#endif

#ifdef MPI2
#include "mpi.h"

extern int my_rank;
extern int iNum_Processors;
extern MPI_Comm *comm;
#endif

attribute_hidden int np_fixed_gaussian_density_cvls_pair_try(
  const int enabled,
  const int num_obs,
  const int num_reg_continuous,
  double * const * const matrix_X_continuous,
  double * const * const matrix_bandwidth,
  double * const cv)
{
#if NP_GAUSSIAN_DENSITY_PAIR_COMPILED
  const double minus_one = -1.0;
  const double zero = 0.0;
  const int vector_length = num_obs;
  const size_t scratch_length = (size_t)num_obs;
  double *difference = NULL;
  double *exponent = NULL;
  double *kernel_value = NULL;
  double convolution_scale = 1.0;
  double kernel_scale = 1.0;
  double bandwidth_product = 1.0;
  double convolution_total = 0.0;
  double kernel_total = 0.0;
  int eval_start = 0;
  int eval_end = num_obs;
#ifdef MPI2
  int stride = 0;
  size_t row_capacity = 0;
  double *convolution_rows = NULL;
  double *kernel_rows = NULL;
#endif

  if(!enabled || (num_obs < 256) || (num_reg_continuous <= 0) ||
     (matrix_X_continuous == NULL) || (matrix_bandwidth == NULL) ||
     (cv == NULL) || (scratch_length > SIZE_MAX/sizeof(double)))
    return 0;

  for(int dimension = 0; dimension < num_reg_continuous; dimension++){
    double h;
    double h2;
    double sqrt_h2;
    double dimension_scale;

    if((matrix_X_continuous[dimension] == NULL) ||
       (matrix_bandwidth[dimension] == NULL))
      return 0;

    h = matrix_bandwidth[dimension][0];
    h2 = h*h + h*h;
    sqrt_h2 = sqrt(h2);
    if((!isfinite(h)) || (h <= 0.0) || (!isfinite(h2)) ||
       (h2 <= 0.0) || (!isfinite(sqrt_h2)) || (sqrt_h2 <= 0.0))
      return 0;

    dimension_scale = (0.3989422803*h*h)/(sqrt_h2*h);
    if(!isfinite(dimension_scale))
      return 0;

    convolution_scale *= dimension_scale;
    kernel_scale *= ONE_OVER_SQRT_TWO_PI;
    bandwidth_product *= h;
  }

  if((!isfinite(convolution_scale)) || (!isfinite(kernel_scale)) ||
     (!isfinite(bandwidth_product)) || (bandwidth_product <= 0.0))
    return 0;

  difference = (double *)malloc(scratch_length*sizeof(double));
  exponent = (double *)malloc(scratch_length*sizeof(double));
  kernel_value = (double *)malloc(scratch_length*sizeof(double));
#ifdef MPI2
  stride = (num_obs + iNum_Processors - 1)/iNum_Processors;
  if((stride <= 0) || (iNum_Processors <= 0) ||
     ((size_t)stride > SIZE_MAX/(size_t)iNum_Processors)){
    free(difference);
    free(exponent);
    free(kernel_value);
    return 0;
  }
  row_capacity = (size_t)stride*(size_t)iNum_Processors;
  if(row_capacity > SIZE_MAX/sizeof(double)){
    free(difference);
    free(exponent);
    free(kernel_value);
    return 0;
  }
  convolution_rows = (double *)calloc(row_capacity, sizeof(double));
  kernel_rows = (double *)calloc(row_capacity, sizeof(double));
  {
    const int local_failure =
      (difference == NULL) || (exponent == NULL) || (kernel_value == NULL) ||
      (convolution_rows == NULL) || (kernel_rows == NULL);
    int any_failure = 0;
    MPI_Allreduce(&local_failure, &any_failure, 1, MPI_INT, MPI_MAX, comm[1]);
    if(any_failure){
      free(difference);
      free(exponent);
      free(kernel_value);
      free(convolution_rows);
      free(kernel_rows);
      return 0;
    }
  }
  eval_start = stride*my_rank;
  eval_end = eval_start + stride;
  if(eval_end > num_obs)
    eval_end = num_obs;
#else
  if((difference == NULL) || (exponent == NULL) || (kernel_value == NULL)){
    free(difference);
    free(exponent);
    free(kernel_value);
    return 0;
  }
#endif

  for(int eval_index = eval_start; eval_index < eval_end; eval_index++){
    double convolution_sum = 0.0;
    double kernel_sum = 0.0;

    for(int dimension = 0; dimension < num_reg_continuous; dimension++){
      const double h = matrix_bandwidth[dimension][0];
      const double h2 = h*h + h*h;
      const double exponent_scale = -0.5/h2;
      const double x = matrix_X_continuous[dimension][eval_index];

      vDSP_vsmsaD(matrix_X_continuous[dimension], 1, &minus_one, &x,
                  difference, 1, (np_vDSP_Length)num_obs);
      vDSP_vsqD(difference, 1, difference, 1,
                (np_vDSP_Length)num_obs);
      if(dimension == 0){
        vDSP_vsmulD(difference, 1, &exponent_scale,
                    exponent, 1, (np_vDSP_Length)num_obs);
      } else {
        vDSP_vsmulD(difference, 1, &exponent_scale,
                    difference, 1, (np_vDSP_Length)num_obs);
        vDSP_vaddD(exponent, 1, difference, 1,
                   exponent, 1, (np_vDSP_Length)num_obs);
      }
    }

    vvexp(kernel_value, exponent, &vector_length);
    vDSP_vsqD(kernel_value, 1, difference, 1,
              (np_vDSP_Length)num_obs);

    vDSP_sveD(kernel_value, 1, &convolution_sum,
              (np_vDSP_Length)num_obs);
#ifdef MPI2
    convolution_rows[eval_index] =
      convolution_scale*convolution_sum/bandwidth_product;
#else
    convolution_total +=
      convolution_scale*convolution_sum/bandwidth_product;
#endif

    difference[eval_index] = zero;
    vDSP_sveD(difference, 1, &kernel_sum,
              (np_vDSP_Length)num_obs);
#ifdef MPI2
    kernel_rows[eval_index] =
      kernel_scale*kernel_sum/bandwidth_product;
#else
    kernel_total += kernel_scale*kernel_sum/bandwidth_product;
#endif
  }

#ifdef MPI2
  MPI_Allgather(MPI_IN_PLACE, stride, MPI_DOUBLE,
                convolution_rows, stride, MPI_DOUBLE, comm[1]);
  MPI_Allgather(MPI_IN_PLACE, stride, MPI_DOUBLE,
                kernel_rows, stride, MPI_DOUBLE, comm[1]);
  for(int eval_index = 0; eval_index < num_obs; eval_index++){
    convolution_total += convolution_rows[eval_index];
    kernel_total += kernel_rows[eval_index];
  }
  free(convolution_rows);
  free(kernel_rows);
#endif
  free(difference);
  free(exponent);
  free(kernel_value);

  *cv = convolution_total/((double)num_obs*(double)num_obs) -
    2.0*kernel_total/((double)num_obs*((double)num_obs - 1.0));
  return isfinite(*cv);
#else
  (void)enabled;
  (void)num_obs;
  (void)num_reg_continuous;
  (void)matrix_X_continuous;
  (void)matrix_bandwidth;
  (void)cv;
  return 0;
#endif
}
