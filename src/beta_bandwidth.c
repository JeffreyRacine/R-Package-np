#include <stddef.h>
#include <string.h>

#include <R.h>

#include "headers.h"
#include "beta_bandwidth.h"

extern int int_LARGE_SF;
extern double nconfac_extern;
extern double ncatfac_extern;

const char *np_beta_bandwidth_prepare_status_message(
  np_beta_bandwidth_prepare_status status)
{
  switch(status) {
  case NP_BETA_BANDWIDTH_PREPARE_OK:
    return "success";
  case NP_BETA_BANDWIDTH_PREPARE_ERR_LAYOUT:
    return "invalid beta nearest-neighbor bandwidth layout";
  case NP_BETA_BANDWIDTH_PREPARE_ERR_DISTANCE:
    return "invalid beta nearest-neighbor bandwidth or distance neighborhood";
  default:
    return "unknown beta nearest-neighbor bandwidth status";
  }
}

np_beta_bandwidth_prepare_status
np_beta_bandwidth_prepare(np_beta_bandwidth_mode bandwidth_mode,
                          const double *train_continuous,
                          const double *eval_continuous,
                          const double *nearest_neighbor,
                          int num_train,
                          int num_eval,
                          int num_continuous,
                          int train_is_eval,
                          int need_eval,
                          int need_train,
                          int suppress_parallel,
                          double *bandwidth_eval,
                          double *bandwidth_train)
{
  double **matrix_train;
  double **matrix_eval;
  double **matrix_bandwidth;
  int large_sf_save;
  double nconfac_save;
  double ncatfac_save;
  int status = 0;
  int dimension;

  if((bandwidth_mode != NP_BETA_BANDWIDTH_GENERALIZED_NN &&
      bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN) ||
     train_continuous == NULL || nearest_neighbor == NULL ||
     num_train <= 0 || num_eval <= 0 || num_continuous <= 0 ||
     (!train_is_eval && eval_continuous == NULL) ||
     (need_eval && bandwidth_eval == NULL) ||
     (need_train && bandwidth_train == NULL))
    return NP_BETA_BANDWIDTH_PREPARE_ERR_LAYOUT;

  matrix_train = (double **)R_alloc((size_t)num_continuous,
                                    sizeof(double *));
  matrix_eval = (double **)R_alloc((size_t)num_continuous,
                                   sizeof(double *));
  matrix_bandwidth = (double **)R_alloc((size_t)num_continuous,
                                        sizeof(double *));
  for(dimension = 0; dimension < num_continuous; ++dimension) {
    matrix_train[dimension] =
      (double *)train_continuous + (size_t)dimension * (size_t)num_train;
    matrix_eval[dimension] = train_is_eval ? matrix_train[dimension] :
      (double *)eval_continuous + (size_t)dimension * (size_t)num_eval;
  }

  large_sf_save = int_LARGE_SF;
  nconfac_save = nconfac_extern;
  ncatfac_save = ncatfac_extern;
  int_LARGE_SF = SF_ARB;
  nconfac_extern = 0.0;
  ncatfac_extern = 0.0;

  if(need_train) {
    for(dimension = 0; dimension < num_continuous; ++dimension)
      matrix_bandwidth[dimension] = bandwidth_train +
        (size_t)dimension * (size_t)num_train;
    status = kernel_bandwidth_mean(
      0,
      (bandwidth_mode == NP_BETA_BANDWIDTH_GENERALIZED_NN) ?
        BW_GEN_NN : BW_ADAP_NN,
      num_train, num_train,
      0, 0, 0, num_continuous, 0, 0, suppress_parallel,
      (double *)nearest_neighbor,
      matrix_train, matrix_train, matrix_train, matrix_train,
      matrix_bandwidth, matrix_bandwidth, NULL);
  }

  if(status == 0 && need_eval) {
    if(bandwidth_mode == NP_BETA_BANDWIDTH_ADAPTIVE_NN && train_is_eval) {
      if(!need_train) {
        status = 1;
      } else {
        memcpy(bandwidth_eval, bandwidth_train,
               (size_t)num_continuous * (size_t)num_train * sizeof(double));
      }
    } else {
      for(dimension = 0; dimension < num_continuous; ++dimension)
        matrix_bandwidth[dimension] = bandwidth_eval +
          (size_t)dimension * (size_t)num_eval;
      status = kernel_bandwidth_mean(
        0, BW_GEN_NN, num_train, num_eval,
        0, 0, 0, num_continuous, 0, 0, suppress_parallel,
        (double *)nearest_neighbor,
        matrix_train, matrix_eval, matrix_train, matrix_eval,
        matrix_bandwidth, matrix_bandwidth, NULL);
    }
  }

  int_LARGE_SF = large_sf_save;
  nconfac_extern = nconfac_save;
  ncatfac_extern = ncatfac_save;

  return (status == 0) ? NP_BETA_BANDWIDTH_PREPARE_OK :
    NP_BETA_BANDWIDTH_PREPARE_ERR_DISTANCE;
}
