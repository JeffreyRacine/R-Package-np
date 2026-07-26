#ifndef NP_JKSUM_GAUSSIAN_DENSITY_H
#define NP_JKSUM_GAUSSIAN_DENSITY_H

#include <R_ext/Visibility.h>

attribute_hidden int np_fixed_gaussian_density_cvls_pair_try(
  int enabled,
  int num_obs,
  int num_reg_continuous,
  double * const *matrix_X_continuous,
  double * const *matrix_bandwidth,
  double *cv);

attribute_hidden int np_fixed_gaussian_density_cvls_pair_dispatch_try(
  int KERNEL_den,
  int BANDWIDTH_den,
  int num_obs,
  int num_reg_unordered,
  int num_reg_ordered,
  int num_reg_continuous,
  double **matrix_X_continuous,
  double *vector_scale_factor,
  double *cv);

#endif
