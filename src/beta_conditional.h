#ifndef NP_BETA_CONDITIONAL_H
#define NP_BETA_CONDITIONAL_H

#include "beta_kernelsum.h"
#include "kernel_registry.h"

typedef enum {
  NP_BETA_CONDITIONAL_OK = 0,
  NP_BETA_CONDITIONAL_ERR_LAYOUT = 1,
  NP_BETA_CONDITIONAL_ERR_KERNEL = 2,
  NP_BETA_CONDITIONAL_ERR_ZERO_WEIGHT = 3,
  NP_BETA_CONDITIONAL_ERR_NUMERIC = 4,
  NP_BETA_CONDITIONAL_ERR_MEMORY = 5
} np_beta_conditional_status;

np_beta_conditional_status
np_beta_conditional_lc(const double *train_x,
                       const double *train_y,
                       const double *eval_x,
                       const double *eval_y,
                       const double *bandwidth_eval_x,
                       const double *bandwidth_train_x,
                       const double *bandwidth_eval_y,
                       const double *bandwidth_train_y,
                       const double *lower_x,
                       const double *upper_x,
                       const double *lower_y,
                       const double *upper_y,
                       np_continuous_kernel_family family_x,
                       int kernel_code_x,
                       int order_x,
                       np_continuous_kernel_family family_y,
                       int kernel_code_y,
                       int order_y,
                       np_beta_bandwidth_mode bandwidth_mode,
                       int do_distribution,
                       int num_train,
                       int num_eval,
                       int num_x,
                       int num_y,
                       int train_is_eval,
                       double *conditional,
                       double *conditional_stderr,
                       double *log_likelihood,
                       int *bad_evaluation,
                       int *bad_dimension,
                       np_beta_status *kernel_status,
                       np_beta_kernelsum_progress_callback progress_callback);

const char *np_beta_conditional_status_message(
  np_beta_conditional_status status);

#endif
