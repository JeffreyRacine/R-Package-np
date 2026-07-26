#ifndef NP_JKSUM_GAUSSIAN_FIXED_H
#define NP_JKSUM_GAUSSIAN_FIXED_H

#include <R_ext/Visibility.h>

attribute_hidden int np_fixed_gaussian_convolution_row_try(
  int kernel,
  const double *xt,
  int num_xt,
  int do_xw,
  double x,
  double hy,
  double h,
  double *result,
  int power);

attribute_hidden int np_fixed_gaussian_convolution_product_try(
  const int *kernels,
  double * const *xt,
  double * const *xeval,
  double * const *bandwidth,
  double * const *alt_bandwidth,
  const int *bandwidth_power,
  int ndim,
  int n,
  int eval_index,
  int bandwidth_index,
  double *result);

#endif
