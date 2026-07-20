#ifndef NP_BETA_KERNELSUM_H
#define NP_BETA_KERNELSUM_H

#include "beta_kernel.h"

typedef enum {
  NP_BETA_KERNELSUM_OK = 0,
  NP_BETA_KERNELSUM_ERR_LAYOUT = 1,
  NP_BETA_KERNELSUM_ERR_KERNEL = 2,
  NP_BETA_KERNELSUM_ERR_NUMERIC = 3
} np_beta_kernelsum_status;

typedef void (*np_beta_kernelsum_progress_callback)(int done,
                                                    int natural_total);

typedef enum {
  NP_BETA_OPERATOR_PDF = 0,
  NP_BETA_OPERATOR_CDF = 1,
  NP_BETA_OPERATOR_OVERLAP = 2,
  NP_BETA_OPERATOR_DERIVATIVE = 3
} np_beta_operator;

typedef enum {
  NP_BETA_BANDWIDTH_FIXED = 0,
  NP_BETA_BANDWIDTH_GENERALIZED_NN = 1,
  NP_BETA_BANDWIDTH_ADAPTIVE_NN = 2
} np_beta_bandwidth_mode;

np_beta_kernelsum_status
np_beta_kernelsum(const double *train_continuous,
                  const double *eval_continuous,
                  const double *response,
                  const double *weights,
                  const double *bandwidth_eval,
                  const double *bandwidth_train,
                  const double *lower,
                  const double *upper,
                  const np_beta_operator *operators,
                  np_beta_bandwidth_mode bandwidth_mode,
                  int order,
                  int num_train,
                  int num_eval,
                  int num_continuous,
                  int num_response_columns,
                  int num_weight_columns,
                  int train_is_eval,
                  int leave_one_out,
                  int return_kernel_weights,
                  double *weighted_sum,
                  double *kernel_weights,
                  double *kernel_square_sum,
                  double *kernel_centered_m2,
                  int *bad_dimension,
                  np_beta_status *kernel_status,
                  np_beta_kernelsum_progress_callback progress_callback);

np_beta_kernelsum_status
np_beta_kernelsum_fixed(const double *train_continuous,
                        const double *eval_continuous,
                        const double *response,
                        const double *weights,
                        const double *bandwidth,
                        const double *lower,
                        const double *upper,
                        const np_beta_operator *operators,
                        int num_train,
                        int num_eval,
                        int num_continuous,
                        int num_response_columns,
                        int num_weight_columns,
                        int train_is_eval,
                        int leave_one_out,
                        int return_kernel_weights,
                        double *weighted_sum,
                        double *kernel_weights,
                        double *kernel_square_sum,
                        double *kernel_centered_m2,
                        int *bad_dimension,
                        np_beta_status *kernel_status,
                        np_beta_kernelsum_progress_callback progress_callback);

np_beta_kernelsum_status
np_beta_kernelsum_derivative(const double *train_continuous,
                             const double *eval_continuous,
                             const double *response,
                             const double *weights,
                             const double *bandwidth_eval,
                             const double *bandwidth_train,
                             const double *lower,
                             const double *upper,
                             const np_beta_operator *operators,
                             int derivative_dimension,
                             np_beta_bandwidth_mode bandwidth_mode,
                             int order,
                             int num_train,
                             int num_eval,
                             int num_continuous,
                             int num_response_columns,
                             int num_weight_columns,
                             int train_is_eval,
                             int leave_one_out,
                             int return_kernel_weights,
                             double *weighted_sum,
                             double *kernel_weights,
                             int *undefined_count,
                             int *bad_dimension,
                             np_beta_status *kernel_status,
                             np_beta_kernelsum_progress_callback progress_callback);

np_beta_kernelsum_status
np_beta_kernelsum_fixed_pdf(const double *train_continuous,
                            const double *eval_continuous,
                            const double *response,
                            const double *weights,
                            const double *bandwidth,
                            const double *lower,
                            const double *upper,
                            int num_train,
                            int num_eval,
                            int num_continuous,
                            int num_response_columns,
                            int num_weight_columns,
                            int train_is_eval,
                            int leave_one_out,
                            int return_kernel_weights,
                            double *weighted_sum,
                            double *kernel_weights,
                            double *kernel_square_sum,
                            int *bad_dimension,
                            np_beta_status *kernel_status,
                            np_beta_kernelsum_progress_callback progress_callback);

const char *np_beta_kernelsum_status_message(np_beta_kernelsum_status status);

#endif
