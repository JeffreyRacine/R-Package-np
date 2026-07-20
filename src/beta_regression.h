#ifndef NP_BETA_REGRESSION_H
#define NP_BETA_REGRESSION_H

#include "beta_kernelsum.h"

typedef enum {
  NP_BETA_REGRESSION_OK = 0,
  NP_BETA_REGRESSION_ERR_LAYOUT = 1,
  NP_BETA_REGRESSION_ERR_KERNEL = 2,
  NP_BETA_REGRESSION_ERR_ZERO_WEIGHT = 3,
  NP_BETA_REGRESSION_ERR_NUMERIC = 4,
  NP_BETA_REGRESSION_ERR_MEMORY = 5
} np_beta_regression_status;

typedef void (*np_beta_regression_progress_callback)(int done,
                                                     int natural_total);

np_beta_regression_status
np_beta_regression_lc(const double *train_continuous,
                      const double *eval_continuous,
                      const double *response,
                      const double *bandwidth_eval,
                      const double *bandwidth_train,
                      const double *lower,
                      const double *upper,
                      np_beta_bandwidth_mode bandwidth_mode,
                      int order,
                      int num_train,
                      int num_eval,
                      int num_continuous,
                      int train_is_eval,
                      double *mean,
                      double *mean_stderr,
                      int *bad_evaluation,
                      int *bad_dimension,
                      np_beta_status *kernel_status,
                      np_beta_regression_progress_callback progress_callback);

np_beta_regression_status
np_beta_regression_fixed_lc(const double *train_continuous,
                            const double *eval_continuous,
                            const double *response,
                            const double *bandwidth,
                            const double *lower,
                            const double *upper,
                            int num_train,
                            int num_eval,
                            int num_continuous,
                            int train_is_eval,
                            double *mean,
                            double *mean_stderr,
                            int *bad_evaluation,
                            int *bad_dimension,
                            np_beta_status *kernel_status,
                            np_beta_regression_progress_callback progress_callback);

np_beta_regression_status
np_beta_regression_lc_gradient(const double *train_continuous,
                               const double *eval_continuous,
                               const double *response,
                               const double *bandwidth_eval,
                               const double *bandwidth_train,
                               const double *lower,
                               const double *upper,
                               np_beta_bandwidth_mode bandwidth_mode,
                               int order,
                               int num_train,
                               int num_eval,
                               int num_continuous,
                               int train_is_eval,
                               double *gradient,
                               double *gradient_stderr,
                               int *infinite_count,
                               int *undefined_count,
                               int *bad_evaluation,
                               int *bad_dimension,
                               np_beta_status *kernel_status);

const char *np_beta_regression_status_message(np_beta_regression_status status);

#endif
