#ifndef NP_BETA_BANDWIDTH_H
#define NP_BETA_BANDWIDTH_H

#include "regression_contract.h"

typedef enum {
  NP_BETA_BANDWIDTH_FIXED = 0,
  NP_BETA_BANDWIDTH_GENERALIZED_NN = 1,
  NP_BETA_BANDWIDTH_ADAPTIVE_NN = 2
} np_beta_bandwidth_mode;

typedef enum {
  NP_BETA_BANDWIDTH_PREPARE_OK = 0,
  NP_BETA_BANDWIDTH_PREPARE_ERR_LAYOUT = 1,
  NP_BETA_BANDWIDTH_PREPARE_ERR_DISTANCE = 2
} np_beta_bandwidth_prepare_status;

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
  double **bandwidth_train);

/* Prepare dimension-major metric bandwidth matrices with the package's
 * existing observation-support nearest-neighbor engine.
 *
 * bandwidth_eval has num_continuous * num_eval elements and always uses the
 * generalized-NN rule. bandwidth_train has num_continuous * num_train
 * elements and uses the public mode: generalized at each training center or
 * adaptive with the center observation excluded. For adaptive train-is-eval
 * overlap, the evaluation matrix is copied from the adaptive training matrix
 * so the same center receives the same bandwidth on both sides. */
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
                          double *bandwidth_train);

const char *np_beta_bandwidth_prepare_status_message(
  np_beta_bandwidth_prepare_status status);

/* Fill a bounded dimension-major continuous-bandwidth layout and matching
 * categorical lambda vector for any canonical beta row consumer. */
int np_beta_continuous_bandwidth_prepare_canonical(
  int bandwidth_mode,
  int num_obs_train,
  int num_obs_eval,
  int num_unordered,
  int num_ordered,
  int num_continuous,
  double **matrix_continuous_train,
  double **matrix_continuous_eval,
  double *vector_scale_factor,
  double **matrix_bandwidth,
  double **matrix_bandwidth_deriv,
  double *lambda,
  const NPContinuousPreparedBandwidthView *prepared_bandwidth);

#endif
