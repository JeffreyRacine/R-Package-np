#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include <R_ext/Arith.h>

#include "beta_regression.h"

const char *np_beta_regression_status_message(np_beta_regression_status status)
{
  switch(status) {
  case NP_BETA_REGRESSION_OK:
    return "success";
  case NP_BETA_REGRESSION_ERR_LAYOUT:
    return "unsupported beta regression layout";
  case NP_BETA_REGRESSION_ERR_KERNEL:
    return "beta scalar kernel evaluation failed";
  case NP_BETA_REGRESSION_ERR_ZERO_WEIGHT:
    return "all beta regression weights are zero at an evaluation point";
  case NP_BETA_REGRESSION_ERR_NUMERIC:
    return "beta regression produced a non-finite result";
  case NP_BETA_REGRESSION_ERR_MEMORY:
    return "beta regression could not allocate its log-weight workspace";
  default:
    return "unknown beta regression status";
  }
}

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
                      np_beta_regression_progress_callback progress_callback)
{
  double *log_weights;
  double *weight_signs;
  int evaluation_index;

  if(bad_evaluation != NULL)
    *bad_evaluation = -1;
  if(bad_dimension != NULL)
    *bad_dimension = -1;
  if(kernel_status != NULL)
    *kernel_status = NP_BETA_OK;

  if(train_continuous == NULL || response == NULL ||
     lower == NULL || upper == NULL || mean == NULL || mean_stderr == NULL ||
     num_train <= 0 || num_eval <= 0 || num_continuous <= 0 ||
     (!train_is_eval && eval_continuous == NULL) ||
     (train_is_eval && num_eval != num_train) ||
     (bandwidth_mode != NP_BETA_BANDWIDTH_FIXED &&
      bandwidth_mode != NP_BETA_BANDWIDTH_GENERALIZED_NN &&
      bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN) ||
     (bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN &&
      bandwidth_eval == NULL) ||
     (bandwidth_mode == NP_BETA_BANDWIDTH_ADAPTIVE_NN &&
      bandwidth_train == NULL) ||
     !np_beta_order_supported(order))
    return NP_BETA_REGRESSION_ERR_LAYOUT;

  if((size_t)num_train > SIZE_MAX / (2U * sizeof(double)))
    return NP_BETA_REGRESSION_ERR_MEMORY;
  log_weights = (double *)malloc(2U * (size_t)num_train * sizeof(double));
  if(log_weights == NULL)
    return NP_BETA_REGRESSION_ERR_MEMORY;
  weight_signs = log_weights + num_train;

  for(evaluation_index = 0; evaluation_index < num_eval; ++evaluation_index) {
    double max_log_weight = -INFINITY;
    double total_weight = 0.0;
    double weighted_mean = 0.0;
    double weighted_m2 = 0.0;
    double squared_weight_sum = 0.0;
    int observation_index;

    for(observation_index = 0; observation_index < num_train;
        ++observation_index) {
      double log_product = 0.0;
      int product_sign = 1;
      int dimension;

      if(!R_FINITE(response[observation_index])) {
        free(log_weights);
        if(bad_evaluation != NULL)
          *bad_evaluation = evaluation_index;
        return NP_BETA_REGRESSION_ERR_NUMERIC;
      }

      for(dimension = 0; dimension < num_continuous; ++dimension) {
        const double evaluation = train_is_eval ?
          train_continuous[dimension * num_train + evaluation_index] :
          eval_continuous[dimension * num_eval + evaluation_index];
        const double observation =
          train_continuous[dimension * num_train + observation_index];
        const double bandwidth =
          (bandwidth_mode == NP_BETA_BANDWIDTH_FIXED) ?
          bandwidth_eval[dimension] :
          ((bandwidth_mode == NP_BETA_BANDWIDTH_GENERALIZED_NN) ?
           bandwidth_eval[dimension * num_eval + evaluation_index] :
           bandwidth_train[dimension * num_train + observation_index]);
        np_beta_status scalar_status = NP_BETA_OK;
        int scalar_sign = 1;
        const double log_value = np_beta_log_abs_pdf_order(
          evaluation, observation, bandwidth,
          lower[dimension], upper[dimension], order,
          &scalar_sign, &scalar_status);
        if(scalar_status != NP_BETA_OK) {
          free(log_weights);
          if(bad_evaluation != NULL)
            *bad_evaluation = evaluation_index;
          if(bad_dimension != NULL)
            *bad_dimension = dimension;
          if(kernel_status != NULL)
            *kernel_status = scalar_status;
          return NP_BETA_REGRESSION_ERR_KERNEL;
        }
        if(log_value == -INFINITY) {
          log_product = -INFINITY;
          product_sign = 0;
          break;
        }
        product_sign *= scalar_sign;
        log_product += log_value;
        if(!R_FINITE(log_product)) {
          free(log_weights);
          if(bad_evaluation != NULL)
            *bad_evaluation = evaluation_index;
          if(bad_dimension != NULL)
            *bad_dimension = dimension;
          return NP_BETA_REGRESSION_ERR_NUMERIC;
        }
      }

      log_weights[observation_index] = log_product;
      weight_signs[observation_index] = (double)product_sign;
      if(log_product > max_log_weight)
        max_log_weight = log_product;
    }

    if(max_log_weight == -INFINITY) {
      free(log_weights);
      if(bad_evaluation != NULL)
        *bad_evaluation = evaluation_index;
      return NP_BETA_REGRESSION_ERR_ZERO_WEIGHT;
    }

    if(order == 2) {
      /* Weighted Welford accumulation is invariant to the common log-sum-exp
       * scaling and avoids cancellation in the positive-weight route. */
      for(observation_index = 0; observation_index < num_train;
          ++observation_index) {
        const double log_weight = log_weights[observation_index];
        const double weight = (log_weight == -INFINITY) ?
          0.0 : exp(log_weight - max_log_weight);

        if(weight > 0.0) {
          const double response_value = response[observation_index];
          const double new_total_weight = total_weight + weight;
          const double delta = response_value - weighted_mean;
          const double new_mean = weighted_mean +
            (weight / new_total_weight) * delta;

          weighted_m2 += weight * delta * (response_value - new_mean);
          squared_weight_sum += weight * weight;
          total_weight = new_total_weight;
          weighted_mean = new_mean;
        }
      }
    } else {
      double weighted_response_sum = 0.0;

      /* Signed higher-order weights are scaled by their largest log absolute
       * value.  The common factor cancels from the Nadaraya--Watson ratio. */
      for(observation_index = 0; observation_index < num_train;
          ++observation_index) {
        const double log_weight = log_weights[observation_index];
        const double weight = (log_weight == -INFINITY) ? 0.0 :
          weight_signs[observation_index] *
          exp(log_weight - max_log_weight);

        total_weight += weight;
        weighted_response_sum += weight * response[observation_index];
      }
      if(!R_FINITE(total_weight) || total_weight == 0.0 ||
         !R_FINITE(weighted_response_sum)) {
        free(log_weights);
        if(bad_evaluation != NULL)
          *bad_evaluation = evaluation_index;
        return (total_weight == 0.0) ? NP_BETA_REGRESSION_ERR_ZERO_WEIGHT :
          NP_BETA_REGRESSION_ERR_NUMERIC;
      }
      weighted_mean = weighted_response_sum / total_weight;

      for(observation_index = 0; observation_index < num_train;
          ++observation_index) {
        const double log_weight = log_weights[observation_index];
        const double weight = (log_weight == -INFINITY) ? 0.0 :
          weight_signs[observation_index] *
          exp(log_weight - max_log_weight);
        const double residual = response[observation_index] - weighted_mean;

        weighted_m2 += weight * residual * residual;
        squared_weight_sum += weight * weight;
      }
    }

    if(!R_FINITE(total_weight) ||
       (order == 2 ? total_weight <= 0.0 : total_weight == 0.0) ||
       !R_FINITE(weighted_mean) || !R_FINITE(weighted_m2) ||
       !R_FINITE(squared_weight_sum)) {
      free(log_weights);
      if(bad_evaluation != NULL)
        *bad_evaluation = evaluation_index;
      return NP_BETA_REGRESSION_ERR_NUMERIC;
    }

    if((order == 2 && weighted_m2 < 0.0) ||
       (order != 2 && weighted_m2 / total_weight < 0.0))
      weighted_m2 = 0.0;
    mean[evaluation_index] = weighted_mean;
    mean_stderr[evaluation_index] = sqrt(
      (weighted_m2 / total_weight) *
      (squared_weight_sum / (total_weight * total_weight)));
    if(!R_FINITE(mean_stderr[evaluation_index])) {
      free(log_weights);
      if(bad_evaluation != NULL)
        *bad_evaluation = evaluation_index;
      return NP_BETA_REGRESSION_ERR_NUMERIC;
    }

    if(progress_callback != NULL)
      progress_callback(evaluation_index + 1, num_eval);
  }

  free(log_weights);
  return NP_BETA_REGRESSION_OK;
}

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
                            np_beta_regression_progress_callback progress_callback)
{
  return np_beta_regression_lc(
    train_continuous, eval_continuous, response, bandwidth, bandwidth,
    lower, upper, NP_BETA_BANDWIDTH_FIXED,
    2,
    num_train, num_eval, num_continuous, train_is_eval,
    mean, mean_stderr, bad_evaluation, bad_dimension, kernel_status,
    progress_callback);
}
