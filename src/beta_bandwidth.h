#ifndef NP_BETA_BANDWIDTH_H
#define NP_BETA_BANDWIDTH_H

#include "beta_kernelsum.h"

typedef enum {
  NP_BETA_BANDWIDTH_PREPARE_OK = 0,
  NP_BETA_BANDWIDTH_PREPARE_ERR_LAYOUT = 1,
  NP_BETA_BANDWIDTH_PREPARE_ERR_DISTANCE = 2
} np_beta_bandwidth_prepare_status;

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

#endif
