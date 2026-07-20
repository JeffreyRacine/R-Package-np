#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include <R_ext/Arith.h>

#include "beta_kernelsum.h"

const char *np_beta_kernelsum_status_message(np_beta_kernelsum_status status)
{
  switch(status) {
  case NP_BETA_KERNELSUM_OK:
    return "success";
  case NP_BETA_KERNELSUM_ERR_LAYOUT:
    return "unsupported beta kernel-sum layout";
  case NP_BETA_KERNELSUM_ERR_KERNEL:
    return "beta scalar kernel evaluation failed";
  case NP_BETA_KERNELSUM_ERR_NUMERIC:
    return "beta kernel sum produced a non-finite result";
  default:
    return "unknown beta kernel-sum status";
  }
}

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
                  np_beta_kernelsum_progress_callback progress_callback)
{
  const int response_extent = (num_response_columns > 0) ?
    num_response_columns : 1;
  const int weight_extent = (num_weight_columns > 0) ?
    num_weight_columns : 1;
  const int sum_extent = response_extent * weight_extent;
  int evaluation_index;
  int observation_index;
  int dimension;
  int response_column;
  int weight_column;

  if(bad_dimension != NULL)
    *bad_dimension = -1;
  if(kernel_status != NULL)
    *kernel_status = NP_BETA_OK;

  if(train_continuous == NULL || lower == NULL ||
     upper == NULL || weighted_sum == NULL || num_train <= 0 ||
     num_eval <= 0 || num_continuous <= 0 ||
     num_response_columns < 0 || num_weight_columns < 0 ||
     (!train_is_eval && eval_continuous == NULL) ||
     (num_response_columns > 0 && response == NULL) ||
     (num_weight_columns > 0 && weights == NULL) ||
     (return_kernel_weights && kernel_weights == NULL) ||
     (bandwidth_mode != NP_BETA_BANDWIDTH_FIXED &&
      bandwidth_mode != NP_BETA_BANDWIDTH_GENERALIZED_NN &&
      bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN) ||
     (bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN &&
      bandwidth_eval == NULL) ||
     (bandwidth_mode == NP_BETA_BANDWIDTH_ADAPTIVE_NN &&
      bandwidth_train == NULL) ||
     !np_beta_order_supported(order))
    return NP_BETA_KERNELSUM_ERR_LAYOUT;

  for(evaluation_index = 0; evaluation_index < num_eval; ++evaluation_index) {
    double kernel_running_mean = 0.0;
    double kernel_running_m2 = 0.0;
    int kernel_running_count = 0;

    if(kernel_square_sum != NULL)
      kernel_square_sum[evaluation_index] = 0.0;
    if(kernel_centered_m2 != NULL)
      kernel_centered_m2[evaluation_index] = 0.0;

    for(response_column = 0; response_column < response_extent;
        ++response_column) {
      for(weight_column = 0; weight_column < weight_extent; ++weight_column) {
        const int output_offset = evaluation_index * sum_extent +
          response_column * weight_extent + weight_column;
        weighted_sum[output_offset] = 0.0;
      }
    }

    for(observation_index = 0; observation_index < num_train;
        ++observation_index) {
      double product = 1.0;

      if(leave_one_out && train_is_eval &&
         observation_index == evaluation_index) {
        if(return_kernel_weights)
          kernel_weights[evaluation_index * num_train + observation_index] = 0.0;
        continue;
      }

      for(dimension = 0; dimension < num_continuous; ++dimension) {
        const double evaluation = train_is_eval ?
          train_continuous[dimension * num_train + evaluation_index] :
          eval_continuous[dimension * num_eval + evaluation_index];
        const double observation =
          train_continuous[dimension * num_train + observation_index];
        np_beta_status scalar_status = NP_BETA_OK;
        double evaluation_bandwidth;
        double observation_bandwidth;
        double value;

        if(bandwidth_mode == NP_BETA_BANDWIDTH_FIXED) {
          evaluation_bandwidth = bandwidth_eval[dimension];
          observation_bandwidth = evaluation_bandwidth;
        } else if(bandwidth_mode == NP_BETA_BANDWIDTH_GENERALIZED_NN) {
          evaluation_bandwidth =
            bandwidth_eval[dimension * num_eval + evaluation_index];
          observation_bandwidth = (bandwidth_train == NULL) ?
            evaluation_bandwidth :
            bandwidth_train[dimension * num_train + observation_index];
        } else {
          observation_bandwidth =
            bandwidth_train[dimension * num_train + observation_index];
          evaluation_bandwidth = observation_bandwidth;
          if(operators != NULL &&
             operators[dimension] == NP_BETA_OPERATOR_OVERLAP) {
            if(bandwidth_eval == NULL) {
              if(bad_dimension != NULL)
                *bad_dimension = dimension;
              return NP_BETA_KERNELSUM_ERR_LAYOUT;
            }
            evaluation_bandwidth =
              bandwidth_eval[dimension * num_eval + evaluation_index];
          }
        }

        if(operators == NULL ||
           operators[dimension] == NP_BETA_OPERATOR_PDF) {
          const double pdf_bandwidth =
            (bandwidth_mode == NP_BETA_BANDWIDTH_ADAPTIVE_NN) ?
            observation_bandwidth : evaluation_bandwidth;
          value = np_beta_pdf_order(evaluation,
                                    observation,
                                    pdf_bandwidth,
                                    lower[dimension],
                                    upper[dimension],
                                    order,
                                    &scalar_status);
        } else if(operators[dimension] == NP_BETA_OPERATOR_CDF) {
          const double cdf_bandwidth =
            (bandwidth_mode == NP_BETA_BANDWIDTH_ADAPTIVE_NN) ?
            observation_bandwidth : evaluation_bandwidth;
          value = np_beta_cdf_order(evaluation,
                                    observation,
                                    cdf_bandwidth,
                                    lower[dimension],
                                    upper[dimension],
                                    order,
                                    &scalar_status);
        } else if(operators[dimension] == NP_BETA_OPERATOR_OVERLAP) {
          value = np_beta_overlap_order(evaluation,
                                        evaluation_bandwidth,
                                        observation,
                                        observation_bandwidth,
                                        lower[dimension],
                                        upper[dimension],
                                        order,
                                        &scalar_status);
        } else {
          if(bad_dimension != NULL)
            *bad_dimension = dimension;
          return NP_BETA_KERNELSUM_ERR_LAYOUT;
        }
        if(scalar_status != NP_BETA_OK) {
          if(bad_dimension != NULL)
            *bad_dimension = dimension;
          if(kernel_status != NULL)
            *kernel_status = scalar_status;
          return NP_BETA_KERNELSUM_ERR_KERNEL;
        }
        product *= value;
        if(!R_FINITE(product)) {
          if(bad_dimension != NULL)
            *bad_dimension = dimension;
          return NP_BETA_KERNELSUM_ERR_NUMERIC;
        }
      }

      if(return_kernel_weights)
        kernel_weights[evaluation_index * num_train + observation_index] = product;

      if(kernel_square_sum != NULL) {
        kernel_square_sum[evaluation_index] += product * product;
        if(!R_FINITE(kernel_square_sum[evaluation_index]))
          return NP_BETA_KERNELSUM_ERR_NUMERIC;
      }

      if(kernel_centered_m2 != NULL) {
        const double delta = product - kernel_running_mean;
        ++kernel_running_count;
        kernel_running_mean += delta / (double)kernel_running_count;
        kernel_running_m2 += delta * (product - kernel_running_mean);
        if(!R_FINITE(kernel_running_mean) || !R_FINITE(kernel_running_m2))
          return NP_BETA_KERNELSUM_ERR_NUMERIC;
      }

      for(response_column = 0; response_column < response_extent;
          ++response_column) {
        const double response_value = (num_response_columns > 0) ?
          response[response_column * num_train + observation_index] : 1.0;
        for(weight_column = 0; weight_column < weight_extent; ++weight_column) {
          const double weight_value = (num_weight_columns > 0) ?
            weights[weight_column * num_train + observation_index] : 1.0;
          const int output_offset = evaluation_index * sum_extent +
            response_column * weight_extent + weight_column;
          weighted_sum[output_offset] += product * response_value * weight_value;
          if(!R_FINITE(weighted_sum[output_offset]))
            return NP_BETA_KERNELSUM_ERR_NUMERIC;
        }
      }
    }

    if(kernel_centered_m2 != NULL) {
      if(kernel_running_m2 < 0.0)
        kernel_running_m2 = 0.0;
      kernel_centered_m2[evaluation_index] = kernel_running_m2;
    }

    if(progress_callback != NULL)
      progress_callback(evaluation_index + 1, num_eval);
  }

  return NP_BETA_KERNELSUM_OK;
}

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
                        np_beta_kernelsum_progress_callback progress_callback)
{
  return np_beta_kernelsum(train_continuous,
                           eval_continuous,
                           response,
                           weights,
                           bandwidth,
                           bandwidth,
                           lower,
                           upper,
                           operators,
                           NP_BETA_BANDWIDTH_FIXED,
                           2,
                           num_train,
                           num_eval,
                           num_continuous,
                           num_response_columns,
                           num_weight_columns,
                           train_is_eval,
                           leave_one_out,
                           return_kernel_weights,
                           weighted_sum,
                           kernel_weights,
                           kernel_square_sum,
                           kernel_centered_m2,
                           bad_dimension,
                           kernel_status,
                           progress_callback);
}

static double np_beta_kernelsum_log_add(double accumulator, double term)
{
  double maximum;
  double minimum;

  if(accumulator == -INFINITY)
    return term;
  if(term == -INFINITY)
    return accumulator;
  maximum = fmax(accumulator, term);
  minimum = fmin(accumulator, term);
  return maximum + log1p(exp(minimum - maximum));
}

static void np_beta_kernelsum_signed_add(double log_absolute,
                                         int sign,
                                         double *positive_log,
                                         double *negative_log)
{
  if(sign > 0)
    *positive_log = np_beta_kernelsum_log_add(*positive_log, log_absolute);
  else if(sign < 0)
    *negative_log = np_beta_kernelsum_log_add(*negative_log, log_absolute);
}

static double np_beta_kernelsum_signed_value(double positive_log,
                                             double negative_log,
                                             int *undefined)
{
  double log_absolute = -INFINITY;
  int sign = 0;
  const np_beta_status status = np_beta_signed_log_absolute(
    positive_log, negative_log, &log_absolute, &sign);

  if(status != NP_BETA_OK || ISNAN(log_absolute)) {
    if(undefined != NULL)
      *undefined = 1;
    return NA_REAL;
  }
  if(sign == 0)
    return 0.0;
  if(log_absolute > log(DBL_MAX))
    return (sign > 0) ? INFINITY : -INFINITY;
  return (sign > 0) ? exp(log_absolute) : -exp(log_absolute);
}

static np_beta_status
np_beta_kernelsum_operator_log_absolute(double evaluation,
                                        double observation,
                                        double evaluation_bandwidth,
                                        double observation_bandwidth,
                                        double lower,
                                        double upper,
                                        np_beta_operator operator,
                                        int order,
                                        double *log_absolute,
                                        int *sign)
{
  np_beta_status status = NP_BETA_OK;
  double value;

  *log_absolute = -INFINITY;
  *sign = 0;
  if(operator == NP_BETA_OPERATOR_PDF) {
    *log_absolute = np_beta_log_abs_pdf_order(
      evaluation, observation, evaluation_bandwidth,
      lower, upper, order, sign, &status);
    return status;
  }
  if(operator == NP_BETA_OPERATOR_CDF)
    value = np_beta_cdf_order(evaluation, observation,
                              evaluation_bandwidth,
                              lower, upper, order, &status);
  else if(operator == NP_BETA_OPERATOR_OVERLAP)
    value = np_beta_overlap_order(evaluation, evaluation_bandwidth,
                                  observation, observation_bandwidth,
                                  lower, upper, order, &status);
  else
    return NP_BETA_ERR_NUMERIC;
  if(status != NP_BETA_OK)
    return status;
  if(value == 0.0)
    return NP_BETA_OK;
  if(!R_FINITE(value))
    return NP_BETA_ERR_RANGE;
  *sign = (value > 0.0) ? 1 : -1;
  *log_absolute = log(fabs(value));
  return NP_BETA_OK;
}

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
                             np_beta_kernelsum_progress_callback progress_callback)
{
  const int response_extent = (num_response_columns > 0) ?
    num_response_columns : 1;
  const int weight_extent = (num_weight_columns > 0) ?
    num_weight_columns : 1;
  const int sum_extent = response_extent * weight_extent;
  double *regular_positive;
  double *regular_negative;
  double *jump_positive;
  double *jump_negative;
  int evaluation_index;
  int observation_index;
  int dimension;
  int response_column;
  int weight_column;

  if(undefined_count != NULL)
    *undefined_count = 0;
  if(bad_dimension != NULL)
    *bad_dimension = -1;
  if(kernel_status != NULL)
    *kernel_status = NP_BETA_OK;
  if(train_continuous == NULL || lower == NULL || upper == NULL ||
     weighted_sum == NULL || num_train <= 0 || num_eval <= 0 ||
     num_continuous <= 0 || derivative_dimension < 0 ||
     derivative_dimension >= num_continuous ||
     num_response_columns < 0 || num_weight_columns < 0 ||
     (!train_is_eval && eval_continuous == NULL) ||
     (num_response_columns > 0 && response == NULL) ||
     (num_weight_columns > 0 && weights == NULL) ||
     (return_kernel_weights && kernel_weights == NULL) ||
     (bandwidth_mode != NP_BETA_BANDWIDTH_FIXED &&
      bandwidth_mode != NP_BETA_BANDWIDTH_GENERALIZED_NN &&
      bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN) ||
     (bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN &&
      bandwidth_eval == NULL) ||
     (bandwidth_mode == NP_BETA_BANDWIDTH_ADAPTIVE_NN &&
      bandwidth_train == NULL) || !np_beta_order_supported(order))
    return NP_BETA_KERNELSUM_ERR_LAYOUT;

  regular_positive = (double *)calloc((size_t)sum_extent, sizeof(double));
  regular_negative = (double *)calloc((size_t)sum_extent, sizeof(double));
  jump_positive = (double *)calloc((size_t)sum_extent, sizeof(double));
  jump_negative = (double *)calloc((size_t)sum_extent, sizeof(double));
  if(regular_positive == NULL || regular_negative == NULL ||
     jump_positive == NULL || jump_negative == NULL) {
    free(regular_positive); free(regular_negative);
    free(jump_positive); free(jump_negative);
    return NP_BETA_KERNELSUM_ERR_NUMERIC;
  }

  for(evaluation_index = 0; evaluation_index < num_eval; ++evaluation_index) {
    for(dimension = 0; dimension < sum_extent; ++dimension) {
      regular_positive[dimension] = -INFINITY;
      regular_negative[dimension] = -INFINITY;
      jump_positive[dimension] = -INFINITY;
      jump_negative[dimension] = -INFINITY;
    }

    for(observation_index = 0; observation_index < num_train;
        ++observation_index) {
      double other_log_absolute = 0.0;
      int other_sign = 1;
      np_beta_derivative derivative;

      if(leave_one_out && train_is_eval &&
         observation_index == evaluation_index) {
        if(return_kernel_weights)
          kernel_weights[evaluation_index * num_train + observation_index] = 0.0;
        continue;
      }

      for(dimension = 0; dimension < num_continuous; ++dimension) {
        const double evaluation = train_is_eval ?
          train_continuous[dimension * num_train + evaluation_index] :
          eval_continuous[dimension * num_eval + evaluation_index];
        const double observation =
          train_continuous[dimension * num_train + observation_index];
        double evaluation_bandwidth;
        double observation_bandwidth;
        np_beta_status status = NP_BETA_OK;

        if(bandwidth_mode == NP_BETA_BANDWIDTH_FIXED) {
          evaluation_bandwidth = bandwidth_eval[dimension];
          observation_bandwidth = evaluation_bandwidth;
        } else if(bandwidth_mode == NP_BETA_BANDWIDTH_GENERALIZED_NN) {
          evaluation_bandwidth =
            bandwidth_eval[dimension * num_eval + evaluation_index];
          observation_bandwidth = (bandwidth_train == NULL) ?
            evaluation_bandwidth :
            bandwidth_train[dimension * num_train + observation_index];
        } else {
          observation_bandwidth =
            bandwidth_train[dimension * num_train + observation_index];
          evaluation_bandwidth = observation_bandwidth;
          if(operators != NULL &&
             operators[dimension] == NP_BETA_OPERATOR_OVERLAP) {
            if(bandwidth_eval == NULL) {
              if(bad_dimension != NULL) *bad_dimension = dimension;
              free(regular_positive); free(regular_negative);
              free(jump_positive); free(jump_negative);
              return NP_BETA_KERNELSUM_ERR_LAYOUT;
            }
            evaluation_bandwidth =
              bandwidth_eval[dimension * num_eval + evaluation_index];
          }
        }

        if(dimension == derivative_dimension) {
          const double derivative_bandwidth =
            (bandwidth_mode == NP_BETA_BANDWIDTH_ADAPTIVE_NN) ?
            observation_bandwidth : evaluation_bandwidth;
          status = np_beta_pdf_derivative_order(
            evaluation, observation, derivative_bandwidth,
            lower[dimension], upper[dimension], order, &derivative);
        } else {
          double factor_log_absolute = -INFINITY;
          int factor_sign = 0;
          const np_beta_operator operator = (operators == NULL) ?
            NP_BETA_OPERATOR_PDF : operators[dimension];
          status = np_beta_kernelsum_operator_log_absolute(
            evaluation, observation,
            (bandwidth_mode == NP_BETA_BANDWIDTH_ADAPTIVE_NN &&
             operator != NP_BETA_OPERATOR_OVERLAP) ?
            observation_bandwidth : evaluation_bandwidth,
            observation_bandwidth, lower[dimension], upper[dimension],
            operator, order, &factor_log_absolute, &factor_sign);
          if(status == NP_BETA_OK) {
            if(factor_sign == 0)
              other_sign = 0;
            else if(other_sign != 0) {
              other_sign *= factor_sign;
              other_log_absolute += factor_log_absolute;
            }
          }
        }
        if(status != NP_BETA_OK) {
          if(bad_dimension != NULL) *bad_dimension = dimension;
          if(kernel_status != NULL) *kernel_status = status;
          free(regular_positive); free(regular_negative);
          free(jump_positive); free(jump_negative);
          return NP_BETA_KERNELSUM_ERR_KERNEL;
        }
      }

      if(return_kernel_weights) {
        double value = 0.0;
        int undefined = 0;
        if(other_sign != 0 && derivative.jump_sign != 0)
          value = (other_sign * derivative.jump_sign > 0) ?
            INFINITY : -INFINITY;
        else if(other_sign != 0 && derivative.regular_sign != 0) {
          const double log_absolute = other_log_absolute +
            derivative.regular_log_absolute;
          if(log_absolute > log(DBL_MAX))
            value = (other_sign * derivative.regular_sign > 0) ?
              INFINITY : -INFINITY;
          else
            value = (other_sign * derivative.regular_sign > 0) ?
              exp(log_absolute) : -exp(log_absolute);
        }
        if(undefined) {
          value = NA_REAL;
          if(undefined_count != NULL) ++*undefined_count;
        }
        kernel_weights[evaluation_index * num_train + observation_index] = value;
      }

      if(other_sign == 0)
        continue;
      for(response_column = 0; response_column < response_extent;
          ++response_column) {
        const double response_value = (num_response_columns > 0) ?
          response[response_column * num_train + observation_index] : 1.0;
        for(weight_column = 0; weight_column < weight_extent; ++weight_column) {
          const double weight_value = (num_weight_columns > 0) ?
            weights[weight_column * num_train + observation_index] : 1.0;
          const int output_offset = response_column * weight_extent +
            weight_column;
          double multiplier_log_absolute;
          int multiplier_sign;

          if(!R_FINITE(response_value) || !R_FINITE(weight_value)) {
            free(regular_positive); free(regular_negative);
            free(jump_positive); free(jump_negative);
            return NP_BETA_KERNELSUM_ERR_NUMERIC;
          }
          if(response_value == 0.0 || weight_value == 0.0)
            continue;
          multiplier_sign = other_sign *
            ((response_value > 0.0) ? 1 : -1) *
            ((weight_value > 0.0) ? 1 : -1);
          multiplier_log_absolute = other_log_absolute +
            log(fabs(response_value)) + log(fabs(weight_value));
          if(derivative.regular_sign != 0)
            np_beta_kernelsum_signed_add(
              multiplier_log_absolute + derivative.regular_log_absolute,
              multiplier_sign * derivative.regular_sign,
              &regular_positive[output_offset],
              &regular_negative[output_offset]);
          if(derivative.jump_sign != 0)
            np_beta_kernelsum_signed_add(
              multiplier_log_absolute + derivative.jump_log_absolute,
              multiplier_sign * derivative.jump_sign,
              &jump_positive[output_offset],
              &jump_negative[output_offset]);
        }
      }
    }

    for(response_column = 0; response_column < response_extent;
        ++response_column) {
      for(weight_column = 0; weight_column < weight_extent; ++weight_column) {
        const int local_offset = response_column * weight_extent + weight_column;
        const int output_offset = evaluation_index * sum_extent + local_offset;
        int undefined = 0;
        double jump = np_beta_kernelsum_signed_value(
          jump_positive[local_offset], jump_negative[local_offset], &undefined);
        if(undefined) {
          weighted_sum[output_offset] = NA_REAL;
          if(undefined_count != NULL) ++*undefined_count;
        } else if(jump != 0.0) {
          weighted_sum[output_offset] = (jump > 0.0) ? INFINITY : -INFINITY;
        } else {
          weighted_sum[output_offset] = np_beta_kernelsum_signed_value(
            regular_positive[local_offset], regular_negative[local_offset],
            &undefined);
          if(undefined && undefined_count != NULL) ++*undefined_count;
        }
      }
    }
    if(progress_callback != NULL)
      progress_callback(evaluation_index + 1, num_eval);
  }

  free(regular_positive); free(regular_negative);
  free(jump_positive); free(jump_negative);
  return NP_BETA_KERNELSUM_OK;
}

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
                            np_beta_kernelsum_progress_callback progress_callback)
{
  return np_beta_kernelsum_fixed(train_continuous,
                                 eval_continuous,
                                 response,
                                 weights,
                                 bandwidth,
                                 lower,
                                 upper,
                                 NULL,
                                 num_train,
                                 num_eval,
                                 num_continuous,
                                 num_response_columns,
                                 num_weight_columns,
                                 train_is_eval,
                                 leave_one_out,
                                 return_kernel_weights,
                                 weighted_sum,
                                 kernel_weights,
                                 kernel_square_sum,
                                 NULL,
                                 bad_dimension,
                                 kernel_status,
                                 progress_callback);
}
