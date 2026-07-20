#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include <R_ext/Arith.h>

#include "beta_conditional.h"

double kernel(int kernel_code, double value);
double cdf_kernel(int kernel_code, double value);
double kernel_deriv(int kernel_code, double value);

static double np_beta_conditional_log_add(double accumulator, double term)
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

static int np_beta_conditional_finite_bound(double value)
{
  return R_FINITE(value) && fabs(value) < 0.5 * DBL_MAX;
}

static double np_beta_conditional_legacy_pdf(int kernel_code,
                                             double evaluation,
                                             double observation,
                                             double bandwidth,
                                             double lower,
                                             double upper)
{
  const int finite_lower = np_beta_conditional_finite_bound(lower);
  const int finite_upper = np_beta_conditional_finite_bound(upper);
  const double base = kernel(kernel_code,
                             (evaluation - observation) / bandwidth);
  double denominator = 1.0;

  if(!R_FINITE(evaluation) || !R_FINITE(observation) ||
     !R_FINITE(bandwidth) || bandwidth <= 0.0 || !R_FINITE(base))
    return NAN;
  if(finite_lower || finite_upper) {
    const double lower_mass = finite_lower ?
      cdf_kernel(kernel_code, (lower - evaluation) / bandwidth) : 0.0;
    const double upper_mass = finite_upper ?
      cdf_kernel(kernel_code, (upper - evaluation) / bandwidth) : 1.0;
    denominator = upper_mass - lower_mass;
  }
  if(!R_FINITE(denominator) || denominator <= 0.0)
    return NAN;
  return base / (bandwidth * denominator);
}

static double np_beta_conditional_legacy_cdf(int kernel_code,
                                             double evaluation,
                                             double observation,
                                             double bandwidth,
                                             double lower,
                                             double upper)
{
  const int finite_lower = np_beta_conditional_finite_bound(lower);
  const int finite_upper = np_beta_conditional_finite_bound(upper);
  double lower_mass;
  double upper_mass;
  double evaluation_mass;
  double denominator;

  if(!R_FINITE(evaluation) || !R_FINITE(observation) ||
     !R_FINITE(bandwidth) || bandwidth <= 0.0)
    return NAN;
  if(finite_lower && evaluation <= lower)
    return 0.0;
  if(finite_upper && evaluation >= upper)
    return 1.0;

  lower_mass = finite_lower ?
    cdf_kernel(kernel_code, (lower - observation) / bandwidth) : 0.0;
  upper_mass = finite_upper ?
    cdf_kernel(kernel_code, (upper - observation) / bandwidth) : 1.0;
  evaluation_mass = cdf_kernel(
    kernel_code, (evaluation - observation) / bandwidth);
  denominator = upper_mass - lower_mass;
  if(!R_FINITE(lower_mass) || !R_FINITE(upper_mass) ||
     !R_FINITE(evaluation_mass) || !R_FINITE(denominator) ||
     denominator <= 0.0)
    return NAN;
  return (evaluation_mass - lower_mass) / denominator;
}

np_beta_conditional_status
np_beta_conditional_scalar_log(np_continuous_kernel_family family,
                               int kernel_code,
                               int order,
                               int do_cdf,
                               double evaluation,
                               double observation,
                               double bandwidth,
                               double lower,
                               double upper,
                               double *log_absolute,
                               int *sign,
                               np_beta_status *kernel_status)
{
  double value;

  if(log_absolute == NULL || sign == NULL)
    return NP_BETA_CONDITIONAL_ERR_LAYOUT;
  *log_absolute = -INFINITY;
  *sign = 0;

  if(family == NP_CKERNEL_FAMILY_BETA) {
    np_beta_status scalar_status = NP_BETA_OK;

    if(do_cdf) {
      value = np_beta_cdf_order(evaluation, observation, bandwidth,
                                lower, upper, order, &scalar_status);
      if(scalar_status == NP_BETA_OK && value != 0.0) {
        *log_absolute = log(fabs(value));
        *sign = (value > 0.0) ? 1 : -1;
      }
    } else {
      *log_absolute = np_beta_log_abs_pdf_order(
        evaluation, observation, bandwidth, lower, upper,
        order, sign, &scalar_status);
    }
    if(scalar_status != NP_BETA_OK) {
      if(kernel_status != NULL)
        *kernel_status = scalar_status;
      return NP_BETA_CONDITIONAL_ERR_KERNEL;
    }
  } else if(family == NP_CKERNEL_FAMILY_LEGACY) {
    value = do_cdf ?
      np_beta_conditional_legacy_cdf(kernel_code, evaluation, observation,
                                     bandwidth, lower, upper) :
      np_beta_conditional_legacy_pdf(kernel_code, evaluation, observation,
                                     bandwidth, lower, upper);
    if(!R_FINITE(value))
      return NP_BETA_CONDITIONAL_ERR_KERNEL;
    if(value != 0.0) {
      *log_absolute = log(fabs(value));
      *sign = (value > 0.0) ? 1 : -1;
    }
  } else {
    return NP_BETA_CONDITIONAL_ERR_LAYOUT;
  }

  if(ISNAN(*log_absolute) || *log_absolute == INFINITY)
    return NP_BETA_CONDITIONAL_ERR_NUMERIC;
  return NP_BETA_CONDITIONAL_OK;
}

static double np_beta_conditional_bandwidth(
  np_beta_bandwidth_mode bandwidth_mode,
  const double *bandwidth_eval,
  const double *bandwidth_train,
  int dimension,
  int evaluation_index,
  int observation_index,
  int num_eval,
  int num_train)
{
  if(bandwidth_mode == NP_BETA_BANDWIDTH_FIXED)
    return bandwidth_eval[dimension];
  if(bandwidth_mode == NP_BETA_BANDWIDTH_GENERALIZED_NN)
    return bandwidth_eval[dimension * num_eval + evaluation_index];
  return bandwidth_train[dimension * num_train + observation_index];
}

const char *np_beta_conditional_status_message(
  np_beta_conditional_status status)
{
  switch(status) {
  case NP_BETA_CONDITIONAL_OK:
    return "success";
  case NP_BETA_CONDITIONAL_ERR_LAYOUT:
    return "unsupported beta conditional-estimator layout";
  case NP_BETA_CONDITIONAL_ERR_KERNEL:
    return "conditional scalar kernel evaluation failed";
  case NP_BETA_CONDITIONAL_ERR_ZERO_WEIGHT:
    return "conditional explanatory-kernel denominator is zero";
  case NP_BETA_CONDITIONAL_ERR_NUMERIC:
    return "beta conditional estimator produced a non-finite result";
  case NP_BETA_CONDITIONAL_ERR_MEMORY:
    return "beta conditional estimator could not allocate its workspace";
  default:
    return "unknown beta conditional-estimator status";
  }
}

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
                       np_beta_kernelsum_progress_callback progress_callback)
{
  double *workspace;
  double *x_log;
  double *y_log;
  double *x_sign;
  double *y_sign;
  int evaluation_index;

  if(bad_evaluation != NULL)
    *bad_evaluation = -1;
  if(bad_dimension != NULL)
    *bad_dimension = -1;
  if(kernel_status != NULL)
    *kernel_status = NP_BETA_OK;

  if(train_x == NULL || train_y == NULL || lower_x == NULL ||
     upper_x == NULL || lower_y == NULL || upper_y == NULL ||
     conditional == NULL || conditional_stderr == NULL ||
     log_likelihood == NULL || num_train <= 0 || num_eval <= 0 ||
     num_x <= 0 || num_y <= 0 || (!train_is_eval &&
      (eval_x == NULL || eval_y == NULL)) ||
     (train_is_eval && num_train != num_eval) ||
     !np_beta_order_supported(order_x) ||
     !np_beta_order_supported(order_y) ||
     (bandwidth_mode != NP_BETA_BANDWIDTH_FIXED &&
      bandwidth_mode != NP_BETA_BANDWIDTH_GENERALIZED_NN &&
      bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN) ||
     (bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN &&
      (bandwidth_eval_x == NULL || bandwidth_eval_y == NULL)) ||
     (bandwidth_mode == NP_BETA_BANDWIDTH_ADAPTIVE_NN &&
      (bandwidth_train_x == NULL || bandwidth_train_y == NULL)))
    return NP_BETA_CONDITIONAL_ERR_LAYOUT;

  if((size_t)num_train > SIZE_MAX / (4U * sizeof(double)))
    return NP_BETA_CONDITIONAL_ERR_MEMORY;
  workspace = (double *)malloc(4U * (size_t)num_train * sizeof(double));
  if(workspace == NULL)
    return NP_BETA_CONDITIONAL_ERR_MEMORY;
  x_log = workspace;
  y_log = workspace + num_train;
  x_sign = workspace + 2 * num_train;
  y_sign = workspace + 3 * num_train;
  *log_likelihood = 0.0;

  for(evaluation_index = 0; evaluation_index < num_eval; ++evaluation_index) {
    double denominator_positive = -INFINITY;
    double denominator_negative = -INFINITY;
    double numerator_positive = -INFINITY;
    double numerator_negative = -INFINITY;
    double denominator_log = -INFINITY;
    double numerator_log = -INFINITY;
    double conditional_log = -INFINITY;
    double squared_influence_log = -INFINITY;
    int denominator_sign = 0;
    int numerator_sign = 0;
    int conditional_sign = 0;
    int observation_index;

    for(observation_index = 0; observation_index < num_train;
        ++observation_index) {
      double log_product_x = 0.0;
      double log_product_y = 0.0;
      int product_sign_x = 1;
      int product_sign_y = 1;
      int dimension;

      for(dimension = 0; dimension < num_x; ++dimension) {
        const double evaluation = train_is_eval ?
          train_x[dimension * num_train + evaluation_index] :
          eval_x[dimension * num_eval + evaluation_index];
        const double observation =
          train_x[dimension * num_train + observation_index];
        const double bandwidth = np_beta_conditional_bandwidth(
          bandwidth_mode, bandwidth_eval_x, bandwidth_train_x,
          dimension, evaluation_index, observation_index,
          num_eval, num_train);
        double scalar_log = -INFINITY;
        int scalar_sign = 0;
        const np_beta_conditional_status scalar_status =
          np_beta_conditional_scalar_log(
            family_x, kernel_code_x, order_x, 0,
            evaluation, observation, bandwidth,
            lower_x[dimension], upper_x[dimension],
            &scalar_log, &scalar_sign, kernel_status);

        if(scalar_status != NP_BETA_CONDITIONAL_OK) {
          free(workspace);
          if(bad_evaluation != NULL)
            *bad_evaluation = evaluation_index;
          if(bad_dimension != NULL)
            *bad_dimension = dimension;
          return scalar_status;
        }
        if(scalar_sign == 0) {
          log_product_x = -INFINITY;
          product_sign_x = 0;
          break;
        }
        log_product_x += scalar_log;
        product_sign_x *= scalar_sign;
      }

      for(dimension = 0; dimension < num_y; ++dimension) {
        const double evaluation = train_is_eval ?
          train_y[dimension * num_train + evaluation_index] :
          eval_y[dimension * num_eval + evaluation_index];
        const double observation =
          train_y[dimension * num_train + observation_index];
        const double bandwidth = np_beta_conditional_bandwidth(
          bandwidth_mode, bandwidth_eval_y, bandwidth_train_y,
          dimension, evaluation_index, observation_index,
          num_eval, num_train);
        double scalar_log = -INFINITY;
        int scalar_sign = 0;
        const np_beta_conditional_status scalar_status =
          np_beta_conditional_scalar_log(
            family_y, kernel_code_y, order_y, do_distribution,
            evaluation, observation, bandwidth,
            lower_y[dimension], upper_y[dimension],
            &scalar_log, &scalar_sign, kernel_status);

        if(scalar_status != NP_BETA_CONDITIONAL_OK) {
          free(workspace);
          if(bad_evaluation != NULL)
            *bad_evaluation = evaluation_index;
          if(bad_dimension != NULL)
            *bad_dimension = num_x + dimension;
          return scalar_status;
        }
        if(scalar_sign == 0) {
          log_product_y = -INFINITY;
          product_sign_y = 0;
          break;
        }
        log_product_y += scalar_log;
        product_sign_y *= scalar_sign;
      }

      x_log[observation_index] = log_product_x;
      y_log[observation_index] = log_product_y;
      x_sign[observation_index] = (double)product_sign_x;
      y_sign[observation_index] = (double)product_sign_y;

      if(product_sign_x > 0)
        denominator_positive = np_beta_conditional_log_add(
          denominator_positive, log_product_x);
      else if(product_sign_x < 0)
        denominator_negative = np_beta_conditional_log_add(
          denominator_negative, log_product_x);

      if(product_sign_x != 0 && product_sign_y != 0) {
        const double joint_log = log_product_x + log_product_y;
        const int joint_sign = product_sign_x * product_sign_y;
        if(joint_sign > 0)
          numerator_positive = np_beta_conditional_log_add(
            numerator_positive, joint_log);
        else
          numerator_negative = np_beta_conditional_log_add(
            numerator_negative, joint_log);
      }
    }

    if(np_beta_signed_log_absolute(
         denominator_positive, denominator_negative,
         &denominator_log, &denominator_sign) != NP_BETA_OK ||
       denominator_sign == 0) {
      free(workspace);
      if(bad_evaluation != NULL)
        *bad_evaluation = evaluation_index;
      return NP_BETA_CONDITIONAL_ERR_ZERO_WEIGHT;
    }
    if(np_beta_signed_log_absolute(
         numerator_positive, numerator_negative,
         &numerator_log, &numerator_sign) != NP_BETA_OK) {
      free(workspace);
      if(bad_evaluation != NULL)
        *bad_evaluation = evaluation_index;
      return NP_BETA_CONDITIONAL_ERR_NUMERIC;
    }

    conditional_sign = numerator_sign * denominator_sign;
    if(conditional_sign == 0) {
      conditional[evaluation_index] = 0.0;
      conditional_log = -INFINITY;
    } else {
      conditional_log = numerator_log - denominator_log;
      if(!R_FINITE(conditional_log) || conditional_log > log(DBL_MAX)) {
        free(workspace);
        if(bad_evaluation != NULL)
          *bad_evaluation = evaluation_index;
        return NP_BETA_CONDITIONAL_ERR_NUMERIC;
      }
      conditional[evaluation_index] = (conditional_sign > 0 ? 1.0 : -1.0) *
        exp(conditional_log);
    }

    if(num_train <= 1) {
      conditional_stderr[evaluation_index] = 0.0;
    } else {
      for(observation_index = 0; observation_index < num_train;
          ++observation_index) {
        double difference_positive = -INFINITY;
        double difference_negative = -INFINITY;
        double difference_log = -INFINITY;
        int difference_sign = 0;
        const int observation_y_sign = (int)y_sign[observation_index];

        if(observation_y_sign > 0)
          difference_positive = y_log[observation_index];
        else if(observation_y_sign < 0)
          difference_negative = y_log[observation_index];
        if(conditional_sign < 0)
          difference_positive = np_beta_conditional_log_add(
            difference_positive, conditional_log);
        else if(conditional_sign > 0)
          difference_negative = np_beta_conditional_log_add(
            difference_negative, conditional_log);

        if(np_beta_signed_log_absolute(
             difference_positive, difference_negative,
             &difference_log, &difference_sign) != NP_BETA_OK) {
          free(workspace);
          if(bad_evaluation != NULL)
            *bad_evaluation = evaluation_index;
          return NP_BETA_CONDITIONAL_ERR_NUMERIC;
        }
        if((int)x_sign[observation_index] != 0 && difference_sign != 0) {
          const double influence_log = x_log[observation_index] +
            difference_log;
          squared_influence_log = np_beta_conditional_log_add(
            squared_influence_log, 2.0 * influence_log);
        }
      }

      if(squared_influence_log == -INFINITY) {
        conditional_stderr[evaluation_index] = 0.0;
      } else {
        const double standard_error_log =
          0.5 * squared_influence_log - denominator_log -
          0.5 * log((double)(num_train - 1));
        if(!R_FINITE(standard_error_log) ||
           standard_error_log > log(DBL_MAX)) {
          free(workspace);
          if(bad_evaluation != NULL)
            *bad_evaluation = evaluation_index;
          return NP_BETA_CONDITIONAL_ERR_NUMERIC;
        }
        conditional_stderr[evaluation_index] = exp(standard_error_log);
      }
    }

    if(!do_distribution)
      *log_likelihood += log((conditional[evaluation_index] < DBL_MIN) ?
                             DBL_MIN : conditional[evaluation_index]);
    if(progress_callback != NULL)
      progress_callback(evaluation_index + 1, num_eval);
  }

  free(workspace);
  return NP_BETA_CONDITIONAL_OK;
}

static np_beta_conditional_status
np_beta_conditional_legacy_pdf_derivative(int kernel_code,
                                          double evaluation,
                                          double observation,
                                          double bandwidth,
                                          double lower,
                                          double upper,
                                          double *log_absolute,
                                          int *sign)
{
  const int finite_lower = np_beta_conditional_finite_bound(lower);
  const int finite_upper = np_beta_conditional_finite_bound(upper);
  double denominator = 1.0;
  double value;

  *log_absolute = -INFINITY;
  *sign = 0;
  if(!R_FINITE(evaluation) || !R_FINITE(observation) ||
     !R_FINITE(bandwidth) || bandwidth <= 0.0)
    return NP_BETA_CONDITIONAL_ERR_KERNEL;
  if(finite_lower || finite_upper) {
    const double lower_mass = finite_lower ?
      cdf_kernel(kernel_code, (lower - evaluation) / bandwidth) : 0.0;
    const double upper_mass = finite_upper ?
      cdf_kernel(kernel_code, (upper - evaluation) / bandwidth) : 1.0;
    denominator = upper_mass - lower_mass;
  }
  if(!R_FINITE(denominator) || denominator <= 0.0)
    return NP_BETA_CONDITIONAL_ERR_KERNEL;

  /* The derivative of the bounded normalizer is common to every
   * observation and therefore cancels exactly in a normalized conditional
   * ratio.  Retaining only the observation-varying term is both sufficient
   * and consistent with the established legacy local-constant route. */
  value = kernel_deriv(kernel_code,
                       (evaluation - observation) / bandwidth) /
    (bandwidth * bandwidth * denominator);
  if(!R_FINITE(value))
    return NP_BETA_CONDITIONAL_ERR_KERNEL;
  if(value != 0.0) {
    *log_absolute = log(fabs(value));
    *sign = (value > 0.0) ? 1 : -1;
  }
  return NP_BETA_CONDITIONAL_OK;
}

np_beta_conditional_status
np_beta_conditional_lc_gradient(const double *train_x,
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
                                double *gradient,
                                double *gradient_stderr,
                                int *infinite_count,
                                int *undefined_count,
                                int *bad_evaluation,
                                int *bad_dimension,
                                np_beta_status *kernel_status)
{
  double *workspace;
  double *level_log;
  double *regular_log;
  double *jump_log;
  double *y_log;
  double *level_sign;
  double *regular_sign;
  double *jump_sign;
  double *y_sign;
  int evaluation_index;
  int derivative_dimension;

  if(infinite_count != NULL) *infinite_count = 0;
  if(undefined_count != NULL) *undefined_count = 0;
  if(bad_evaluation != NULL) *bad_evaluation = -1;
  if(bad_dimension != NULL) *bad_dimension = -1;
  if(kernel_status != NULL) *kernel_status = NP_BETA_OK;
  if(train_x == NULL || train_y == NULL || lower_x == NULL ||
     upper_x == NULL || lower_y == NULL || upper_y == NULL ||
     gradient == NULL || gradient_stderr == NULL || num_train <= 0 ||
     num_eval <= 0 || num_x <= 0 || num_y <= 0 ||
     (!train_is_eval && (eval_x == NULL || eval_y == NULL)) ||
     (train_is_eval && num_train != num_eval) ||
     !np_beta_order_supported(order_x) ||
     !np_beta_order_supported(order_y) ||
     (bandwidth_mode != NP_BETA_BANDWIDTH_FIXED &&
      bandwidth_mode != NP_BETA_BANDWIDTH_GENERALIZED_NN &&
      bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN) ||
     (bandwidth_mode != NP_BETA_BANDWIDTH_ADAPTIVE_NN &&
      (bandwidth_eval_x == NULL || bandwidth_eval_y == NULL)) ||
     (bandwidth_mode == NP_BETA_BANDWIDTH_ADAPTIVE_NN &&
      (bandwidth_train_x == NULL || bandwidth_train_y == NULL)))
    return NP_BETA_CONDITIONAL_ERR_LAYOUT;
  if((size_t)num_train > SIZE_MAX / (8U * sizeof(double)))
    return NP_BETA_CONDITIONAL_ERR_MEMORY;
  workspace = (double *)malloc(8U * (size_t)num_train * sizeof(double));
  if(workspace == NULL)
    return NP_BETA_CONDITIONAL_ERR_MEMORY;
  level_log = workspace;
  regular_log = level_log + num_train;
  jump_log = regular_log + num_train;
  y_log = jump_log + num_train;
  level_sign = y_log + num_train;
  regular_sign = level_sign + num_train;
  jump_sign = regular_sign + num_train;
  y_sign = jump_sign + num_train;

  for(evaluation_index = 0; evaluation_index < num_eval; ++evaluation_index) {
    for(derivative_dimension = 0; derivative_dimension < num_x;
        ++derivative_dimension) {
      double maximum_x_log = -INFINITY;
      double total_weight = 0.0;
      double weighted_y = 0.0;
      double regular_total = 0.0;
      double regular_y = 0.0;
      double jump_total = 0.0;
      double jump_y = 0.0;
      int observation_index;

      for(observation_index = 0; observation_index < num_train;
          ++observation_index) {
        double other_x_log = 0.0;
        int other_x_sign = 1;
        double derivative_level_log = -INFINITY;
        int derivative_level_sign = 0;
        np_beta_derivative derivative;
        int dimension;

        derivative.regular_log_absolute = -INFINITY;
        derivative.jump_log_absolute = -INFINITY;
        derivative.regular_sign = 0;
        derivative.jump_sign = 0;
        for(dimension = 0; dimension < num_x; ++dimension) {
          const double evaluation = train_is_eval ?
            train_x[dimension * num_train + evaluation_index] :
            eval_x[dimension * num_eval + evaluation_index];
          const double observation =
            train_x[dimension * num_train + observation_index];
          const double bandwidth = np_beta_conditional_bandwidth(
            bandwidth_mode, bandwidth_eval_x, bandwidth_train_x,
            dimension, evaluation_index, observation_index,
            num_eval, num_train);
          double scalar_log = -INFINITY;
          int scalar_sign = 0;
          np_beta_conditional_status status =
            np_beta_conditional_scalar_log(
              family_x, kernel_code_x, order_x, 0,
              evaluation, observation, bandwidth,
              lower_x[dimension], upper_x[dimension],
              &scalar_log, &scalar_sign, kernel_status);

          if(status == NP_BETA_CONDITIONAL_OK &&
             dimension == derivative_dimension) {
            derivative_level_log = scalar_log;
            derivative_level_sign = scalar_sign;
            if(family_x == NP_CKERNEL_FAMILY_BETA) {
              const np_beta_status beta_status =
                np_beta_pdf_derivative_order(
                  evaluation, observation, bandwidth,
                  lower_x[dimension], upper_x[dimension],
                  order_x, &derivative);
              if(beta_status != NP_BETA_OK) {
                if(kernel_status != NULL) *kernel_status = beta_status;
                status = NP_BETA_CONDITIONAL_ERR_KERNEL;
              }
            } else {
              status = np_beta_conditional_legacy_pdf_derivative(
                kernel_code_x, evaluation, observation, bandwidth,
                lower_x[dimension], upper_x[dimension],
                &derivative.regular_log_absolute,
                &derivative.regular_sign);
            }
          } else if(status == NP_BETA_CONDITIONAL_OK) {
            if(scalar_sign == 0)
              other_x_sign = 0;
            else if(other_x_sign != 0) {
              other_x_sign *= scalar_sign;
              other_x_log += scalar_log;
            }
          }
          if(status != NP_BETA_CONDITIONAL_OK) {
            free(workspace);
            if(bad_evaluation != NULL) *bad_evaluation = evaluation_index;
            if(bad_dimension != NULL) *bad_dimension = dimension;
            return status;
          }
        }

        if(other_x_sign == 0 || derivative_level_sign == 0) {
          level_log[observation_index] = -INFINITY;
          level_sign[observation_index] = 0.0;
        } else {
          level_log[observation_index] =
            other_x_log + derivative_level_log;
          level_sign[observation_index] =
            (double)(other_x_sign * derivative_level_sign);
          maximum_x_log = fmax(maximum_x_log, level_log[observation_index]);
        }
        if(other_x_sign == 0 || derivative.regular_sign == 0) {
          regular_log[observation_index] = -INFINITY;
          regular_sign[observation_index] = 0.0;
        } else {
          regular_log[observation_index] =
            other_x_log + derivative.regular_log_absolute;
          regular_sign[observation_index] =
            (double)(other_x_sign * derivative.regular_sign);
          maximum_x_log = fmax(maximum_x_log, regular_log[observation_index]);
        }
        if(other_x_sign == 0 || derivative.jump_sign == 0) {
          jump_log[observation_index] = -INFINITY;
          jump_sign[observation_index] = 0.0;
        } else {
          jump_log[observation_index] =
            other_x_log + derivative.jump_log_absolute;
          jump_sign[observation_index] =
            (double)(other_x_sign * derivative.jump_sign);
          maximum_x_log = fmax(maximum_x_log, jump_log[observation_index]);
        }

        {
          double product_y_log = 0.0;
          int product_y_sign = 1;
          for(dimension = 0; dimension < num_y; ++dimension) {
            const double evaluation = train_is_eval ?
              train_y[dimension * num_train + evaluation_index] :
              eval_y[dimension * num_eval + evaluation_index];
            const double observation =
              train_y[dimension * num_train + observation_index];
            const double bandwidth = np_beta_conditional_bandwidth(
              bandwidth_mode, bandwidth_eval_y, bandwidth_train_y,
              dimension, evaluation_index, observation_index,
              num_eval, num_train);
            double scalar_log = -INFINITY;
            int scalar_sign = 0;
            const np_beta_conditional_status status =
              np_beta_conditional_scalar_log(
                family_y, kernel_code_y, order_y, do_distribution,
                evaluation, observation, bandwidth,
                lower_y[dimension], upper_y[dimension],
                &scalar_log, &scalar_sign, kernel_status);
            if(status != NP_BETA_CONDITIONAL_OK) {
              free(workspace);
              if(bad_evaluation != NULL) *bad_evaluation = evaluation_index;
              if(bad_dimension != NULL) *bad_dimension = num_x + dimension;
              return status;
            }
            if(scalar_sign == 0) {
              product_y_sign = 0;
              product_y_log = -INFINITY;
              break;
            }
            product_y_sign *= scalar_sign;
            product_y_log += scalar_log;
          }
          y_log[observation_index] = product_y_log;
          y_sign[observation_index] = (double)product_y_sign;
        }
      }

      if(maximum_x_log == -INFINITY) {
        free(workspace);
        if(bad_evaluation != NULL) *bad_evaluation = evaluation_index;
        return NP_BETA_CONDITIONAL_ERR_ZERO_WEIGHT;
      }
      for(observation_index = 0; observation_index < num_train;
          ++observation_index) {
        double q = 0.0;
        const double w = (level_sign[observation_index] == 0.0) ? 0.0 :
          level_sign[observation_index] *
          exp(level_log[observation_index] - maximum_x_log);
        const double d = (regular_sign[observation_index] == 0.0) ? 0.0 :
          regular_sign[observation_index] *
          exp(regular_log[observation_index] - maximum_x_log);
        const double j = (jump_sign[observation_index] == 0.0) ? 0.0 :
          jump_sign[observation_index] *
          exp(jump_log[observation_index] - maximum_x_log);
        if(y_sign[observation_index] != 0.0) {
          if(y_log[observation_index] > log(DBL_MAX)) {
            free(workspace);
            if(bad_evaluation != NULL) *bad_evaluation = evaluation_index;
            return NP_BETA_CONDITIONAL_ERR_NUMERIC;
          }
          q = y_sign[observation_index] * exp(y_log[observation_index]);
        }
        total_weight += w;
        weighted_y += w * q;
        regular_total += d;
        regular_y += d * q;
        jump_total += j;
        jump_y += j * q;
      }

      {
        const double evaluation = train_is_eval ?
          train_x[derivative_dimension * num_train + evaluation_index] :
          eval_x[derivative_dimension * num_eval + evaluation_index];
        const int at_lower = family_x == NP_CKERNEL_FAMILY_BETA &&
          evaluation == lower_x[derivative_dimension];
        const int at_upper = family_x == NP_CKERNEL_FAMILY_BETA &&
          evaluation == upper_x[derivative_dimension];
        const double side_orientation = at_upper ? -1.0 : 1.0;
        double side_weight = total_weight;
        double side_y = weighted_y;
        double base_value;
        double side_value;
        double gradient_value;
        int q_is_constant = 1;
        int have_active = 0;
        double first_q = 0.0;

        if(!R_FINITE(total_weight) || total_weight == 0.0 ||
           !R_FINITE(weighted_y)) {
          free(workspace);
          if(bad_evaluation != NULL) *bad_evaluation = evaluation_index;
          return NP_BETA_CONDITIONAL_ERR_NUMERIC;
        }
        base_value = weighted_y / total_weight;
        if(at_lower || at_upper) {
          side_weight += side_orientation * jump_total;
          side_y += side_orientation * jump_y;
        }
        if(!R_FINITE(side_weight) || side_weight == 0.0 ||
           !R_FINITE(side_y)) {
          gradient[derivative_dimension * num_eval + evaluation_index] = NA_REAL;
          gradient_stderr[derivative_dimension * num_eval + evaluation_index] = NA_REAL;
          if(undefined_count != NULL) ++*undefined_count;
          continue;
        }
        side_value = side_y / side_weight;

        for(observation_index = 0; observation_index < num_train;
            ++observation_index) {
          double q = 0.0;
          const double w = (level_sign[observation_index] == 0.0) ? 0.0 :
            level_sign[observation_index] *
            exp(level_log[observation_index] - maximum_x_log);
          const double j = (jump_sign[observation_index] == 0.0) ? 0.0 :
            jump_sign[observation_index] *
            exp(jump_log[observation_index] - maximum_x_log);
          if(y_sign[observation_index] != 0.0)
            q = y_sign[observation_index] * exp(y_log[observation_index]);
          if(w != 0.0 || j != 0.0) {
            if(!have_active) {
              first_q = q;
              have_active = 1;
            } else if(q != first_q) {
              q_is_constant = 0;
            }
          }
        }

        if(q_is_constant) {
          gradient[derivative_dimension * num_eval + evaluation_index] = 0.0;
          gradient_stderr[derivative_dimension * num_eval + evaluation_index] = 0.0;
          continue;
        }

        if((at_lower || at_upper) && side_value != base_value) {
          const double jump_in_ratio = side_orientation *
            (side_value - base_value);
          const double tolerance = 128.0 * DBL_EPSILON *
            fmax(1.0, fmax(fabs(side_value), fabs(base_value)));
          if(fabs(jump_in_ratio) <= tolerance) {
            gradient[derivative_dimension * num_eval + evaluation_index] = NA_REAL;
            gradient_stderr[derivative_dimension * num_eval + evaluation_index] = NA_REAL;
            if(undefined_count != NULL) ++*undefined_count;
          } else {
            gradient[derivative_dimension * num_eval + evaluation_index] =
              (jump_in_ratio > 0.0) ? INFINITY : -INFINITY;
            gradient_stderr[derivative_dimension * num_eval + evaluation_index] = NA_REAL;
            if(infinite_count != NULL) ++*infinite_count;
          }
          continue;
        }

        gradient_value = (regular_y - side_value * regular_total) /
          side_weight;
        gradient[derivative_dimension * num_eval + evaluation_index] =
          gradient_value;
        if(num_train <= 1) {
          gradient_stderr[derivative_dimension * num_eval + evaluation_index] = 0.0;
        } else {
          double influence_square_sum = 0.0;
          for(observation_index = 0; observation_index < num_train;
              ++observation_index) {
            double q = 0.0;
            const double w = (level_sign[observation_index] == 0.0) ? 0.0 :
              level_sign[observation_index] *
              exp(level_log[observation_index] - maximum_x_log);
            const double d = (regular_sign[observation_index] == 0.0) ? 0.0 :
              regular_sign[observation_index] *
              exp(regular_log[observation_index] - maximum_x_log);
            const double j = (jump_sign[observation_index] == 0.0) ? 0.0 :
              jump_sign[observation_index] *
              exp(jump_log[observation_index] - maximum_x_log);
            const double side_w = w +
              ((at_lower || at_upper) ? side_orientation * j : 0.0);
            const double a = side_w / side_weight;
            const double aprime =
              (d * side_weight - side_w * regular_total) /
              (side_weight * side_weight);
            double influence;
            if(y_sign[observation_index] != 0.0)
              q = y_sign[observation_index] * exp(y_log[observation_index]);
            influence = aprime * (q - side_value) - gradient_value * a;
            influence_square_sum += influence * influence;
          }
          gradient_stderr[derivative_dimension * num_eval + evaluation_index] =
            sqrt(influence_square_sum / (double)(num_train - 1));
        }
        if(!R_FINITE(gradient_value) ||
           !R_FINITE(gradient_stderr[derivative_dimension * num_eval +
                                     evaluation_index])) {
          gradient[derivative_dimension * num_eval + evaluation_index] = NA_REAL;
          gradient_stderr[derivative_dimension * num_eval + evaluation_index] = NA_REAL;
          if(undefined_count != NULL) ++*undefined_count;
        }
      }
    }
  }

  free(workspace);
  return NP_BETA_CONDITIONAL_OK;
}
