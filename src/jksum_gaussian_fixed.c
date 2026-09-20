#include <math.h>

#include "headers.h"
#include "jksum_gaussian_fixed.h"

attribute_hidden int np_fixed_gaussian_convolution_row_try(
  const int kernel,
  const double * const xt,
  const int num_xt,
  const int do_xw,
  const double x,
  const double hy,
  const double h,
  double * const result,
  const int power)
{
  int i;
  np_gaussian_convolution_polynomial p;
  if((xt == NULL) || (result == NULL) || (num_xt < 0) ||
     ((kernel != 1) && (kernel != 2)))
    return 0;

  p = np_gaussian_convolution_prepare(kernel, h, hy);
  p.scale /= ipow(hy, power);
  for(i = 0; i < num_xt; i++){
    const double weight = do_xw > 0 ? result[i] : 1.0;
    if(weight == 0.0)
      continue;
    result[i] = weight*np_gaussian_convolution_evaluate(&p, kernel, x-xt[i]);
  }
  return 1;
}

attribute_hidden int np_fixed_gaussian_convolution_product_try(
  const int * const kernels,
  double * const * const xt,
  double * const * const xeval,
  double * const * const bandwidth,
  double * const * const alt_bandwidth,
  const int * const bandwidth_power,
  const int ndim,
  const int n,
  const int eval_index,
  const int bandwidth_index,
  double * const result)
{
  int dimension;

  if((kernels == NULL) || (xt == NULL) || (xeval == NULL) ||
     (bandwidth == NULL) || (alt_bandwidth == NULL) ||
     (bandwidth_power == NULL) || (result == NULL) ||
     (ndim <= 0) || (n < 0))
    return 0;

  for(dimension = 0; dimension < ndim; dimension++){
    if(((kernels[dimension] != 1) && (kernels[dimension] != 2)) ||
       (xt[dimension] == NULL) || (xeval[dimension] == NULL) ||
       (bandwidth[dimension] == NULL) ||
       (alt_bandwidth[dimension] == NULL))
      return 0;
  }

  for(dimension = 0; dimension < ndim; dimension++){
    if(!np_fixed_gaussian_convolution_row_try(
         kernels[dimension], xt[dimension], n, dimension != 0,
         xeval[dimension][eval_index], alt_bandwidth[dimension][0],
         bandwidth[dimension][bandwidth_index], result,
         bandwidth_power[dimension]))
      return 0;
  }

  return 1;
}
