#include <stddef.h>
#include <string.h>

#include <R.h>

#include "headers.h"
#include "beta_bandwidth.h"

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
  case NP_BETA_BANDWIDTH_PREPARE_ERR_ZERO_RADIUS:
    return "beta nearest-neighbor bandwidth has a zero literal radius";
  default:
    return "unknown beta nearest-neighbor bandwidth status";
  }
}

np_beta_bandwidth_prepare_status
np_beta_bandwidth_prepare_matrix(
  np_beta_bandwidth_mode bandwidth_mode,
  double * const *train_continuous,
  double * const *eval_continuous,
  const double *nearest_neighbor,
  int num_train,
  int num_eval,
  int num_continuous,
  int train_is_eval,
  int need_eval,
  int need_train,
  int suppress_parallel,
  double **bandwidth_eval,
  double **bandwidth_train)
{
  int status = 0;
  int dimension;
  NPNNGeometryStatus geometry_status = NP_NN_GEOMETRY_OK;

  if((bandwidth_mode != NP_BETA_BANDWIDTH_GENERALIZED_NN &&
      bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN) ||
     train_continuous == NULL || nearest_neighbor == NULL ||
     num_train <= 0 || num_eval <= 0 || num_continuous <= 0 ||
     (train_is_eval != 0 && train_is_eval != 1) ||
     (train_is_eval && num_train != num_eval) ||
     (!train_is_eval && eval_continuous == NULL) ||
     (need_eval && bandwidth_eval == NULL) ||
     (need_train && bandwidth_train == NULL))
    return NP_BETA_BANDWIDTH_PREPARE_ERR_LAYOUT;
  for(dimension = 0; dimension < num_continuous; ++dimension)
    if(train_continuous[dimension] == NULL ||
       (!train_is_eval && eval_continuous[dimension] == NULL) ||
       (need_eval && bandwidth_eval[dimension] == NULL) ||
       (need_train && bandwidth_train[dimension] == NULL))
      return NP_BETA_BANDWIDTH_PREPARE_ERR_LAYOUT;

  if(need_train)
    status = np_kernel_bandwidth_continuous_nn(
      (bandwidth_mode == NP_BETA_BANDWIDTH_GENERALIZED_NN) ?
        BW_GEN_NN : BW_ADAP_NN,
      num_train, num_train, num_continuous, suppress_parallel,
      (double *)nearest_neighbor,
      (double **)train_continuous, (double **)train_continuous,
      bandwidth_train, &geometry_status);

  if(status == 0 && need_eval) {
    if(bandwidth_mode == NP_BETA_BANDWIDTH_ADAPTIVE_NN && train_is_eval) {
      if(!need_train) {
        status = 1;
      } else {
        for(dimension = 0; dimension < num_continuous; ++dimension)
          memcpy(bandwidth_eval[dimension], bandwidth_train[dimension],
                 (size_t)num_train * sizeof(double));
      }
    } else {
      status = np_kernel_bandwidth_continuous_nn(
        BW_GEN_NN, num_train, num_eval, num_continuous, suppress_parallel,
        (double *)nearest_neighbor,
        (double **)train_continuous, (double **)eval_continuous,
        bandwidth_eval, &geometry_status);
    }
  }

  if(status == 0)
    return NP_BETA_BANDWIDTH_PREPARE_OK;
  return geometry_status == NP_NN_GEOMETRY_ZERO_RADIUS ?
    NP_BETA_BANDWIDTH_PREPARE_ERR_ZERO_RADIUS :
    NP_BETA_BANDWIDTH_PREPARE_ERR_DISTANCE;
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
  int dimension;

  if((bandwidth_mode != NP_BETA_BANDWIDTH_GENERALIZED_NN &&
      bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN) ||
     train_continuous == NULL || nearest_neighbor == NULL ||
     num_train <= 0 || num_eval <= 0 || num_continuous <= 0 ||
     (train_is_eval != 0 && train_is_eval != 1) ||
     (train_is_eval && num_train != num_eval) ||
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

  if(need_train)
    for(dimension = 0; dimension < num_continuous; ++dimension)
      matrix_bandwidth[dimension] = bandwidth_train +
        (size_t)dimension * (size_t)num_train;
  if(need_eval) {
    double **matrix_bandwidth_eval = (double **)R_alloc(
      (size_t)num_continuous, sizeof(double *));

    for(dimension = 0; dimension < num_continuous; ++dimension)
      matrix_bandwidth_eval[dimension] = bandwidth_eval +
        (size_t)dimension * (size_t)num_eval;
    return np_beta_bandwidth_prepare_matrix(
      bandwidth_mode, matrix_train, matrix_eval, nearest_neighbor,
      num_train, num_eval, num_continuous, train_is_eval,
      need_eval, need_train, suppress_parallel,
      matrix_bandwidth_eval, matrix_bandwidth);
  }
  return np_beta_bandwidth_prepare_matrix(
    bandwidth_mode, matrix_train, matrix_eval, nearest_neighbor,
    num_train, num_eval, num_continuous, train_is_eval,
    need_eval, need_train, suppress_parallel, NULL, matrix_bandwidth);
}
