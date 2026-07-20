#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include <R.h>
#include <R_ext/Arith.h>

#include "headers.h"
#include "beta_conditional.h"
#include "beta_kernel.h"
#include "beta_objectives.h"

#define NP_BETA_CVLS_GRID_1D 81
#define NP_BETA_CVLS_GRID_2D 31

typedef struct {
  double positive;
  double negative;
} np_beta_log_sum;

static double np_beta_objective_log_add(double accumulator, double term)
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

static void np_beta_log_sum_init(np_beta_log_sum *sum)
{
  sum->positive = -INFINITY;
  sum->negative = -INFINITY;
}

static void np_beta_log_sum_add(np_beta_log_sum *sum,
                                double log_absolute,
                                int sign)
{
  if(sign > 0)
    sum->positive = np_beta_objective_log_add(sum->positive, log_absolute);
  else if(sign < 0)
    sum->negative = np_beta_objective_log_add(sum->negative, log_absolute);
}

static int np_beta_log_sum_finish(const np_beta_log_sum *sum,
                                  double *log_absolute,
                                  int *sign)
{
  return np_beta_signed_log_absolute(sum->positive,
                                     sum->negative,
                                     log_absolute,
                                     sign) == NP_BETA_OK ? 0 : 1;
}

static int np_beta_log_signed_value(double log_absolute,
                                    int sign,
                                    double *value)
{
  double magnitude;

  if(value == NULL || ISNAN(log_absolute) || log_absolute == INFINITY)
    return 1;
  if(sign == 0 || log_absolute == -INFINITY) {
    *value = 0.0;
    return 0;
  }
  if(log_absolute > log(DBL_MAX))
    return 1;
  magnitude = exp(log_absolute);
  if(!R_FINITE(magnitude))
    return 1;
  *value = (sign > 0) ? magnitude : -magnitude;
  return 0;
}

static int np_beta_objective_bandwidth_rows(int bandwidth_mode,
                                            int num_train,
                                            int num_eval)
{
  if(bandwidth_mode == BW_FIXED)
    return 1;
  if(bandwidth_mode == BW_GEN_NN)
    return num_eval;
  if(bandwidth_mode == BW_ADAP_NN)
    return num_train;
  return 0;
}

static int np_beta_objective_prepare_bandwidth(
  int bandwidth_mode,
  int num_train,
  int num_eval,
  int num_continuous,
  double **train_continuous,
  double **eval_continuous,
  const double *candidate,
  double **matrix_bandwidth)
{
  if(matrix_bandwidth == NULL || train_continuous == NULL ||
     eval_continuous == NULL || candidate == NULL || num_train <= 0 ||
     num_eval <= 0 || num_continuous <= 0)
    return 1;

  return kernel_bandwidth_mean(
    0,
    bandwidth_mode,
    num_train,
    num_eval,
    0, 0, 0,
    num_continuous, 0, 0,
    0,
    (double *)candidate,
    NULL, NULL,
    train_continuous,
    eval_continuous,
    NULL,
    matrix_bandwidth,
    NULL);
}

static int np_beta_objective_prepare_bandwidth_offset(
  int bandwidth_mode,
  int num_train,
  int num_eval,
  int num_continuous,
  int standard_deviation_offset,
  double **train_continuous,
  double **eval_continuous,
  const double *candidate,
  double **matrix_bandwidth)
{
  extern double *vector_continuous_stddev_extern;
  double *saved_standard_deviation = vector_continuous_stddev_extern;
  int status;

  if(standard_deviation_offset < 0)
    return 1;
  if(saved_standard_deviation != NULL)
    vector_continuous_stddev_extern =
      saved_standard_deviation + standard_deviation_offset;
  status = np_beta_objective_prepare_bandwidth(
    bandwidth_mode, num_train, num_eval, num_continuous,
    train_continuous, eval_continuous, candidate, matrix_bandwidth);
  vector_continuous_stddev_extern = saved_standard_deviation;
  return status;
}

static double np_beta_objective_bandwidth(
  int bandwidth_mode,
  double **matrix_bandwidth,
  int dimension,
  int evaluation_index,
  int observation_index)
{
  if(bandwidth_mode == BW_FIXED)
    return matrix_bandwidth[dimension][0];
  if(bandwidth_mode == BW_GEN_NN)
    return matrix_bandwidth[dimension][evaluation_index];
  return matrix_bandwidth[dimension][observation_index];
}

static int np_beta_objective_product_log(
  int do_cdf,
  int bandwidth_mode,
  int order,
  int num_train,
  int num_eval,
  int num_continuous,
  double **train_continuous,
  double **eval_continuous,
  double **matrix_bandwidth,
  const double *lower,
  const double *upper,
  int evaluation_index,
  int observation_index,
  double *log_absolute,
  int *sign)
{
  double log_product = 0.0;
  int product_sign = 1;
  int dimension;

  (void)num_train;
  (void)num_eval;
  if(log_absolute == NULL || sign == NULL)
    return 1;

  for(dimension = 0; dimension < num_continuous; ++dimension) {
    const double evaluation = eval_continuous[dimension][evaluation_index];
    const double observation = train_continuous[dimension][observation_index];
    const double bandwidth = np_beta_objective_bandwidth(
      bandwidth_mode, matrix_bandwidth, dimension,
      evaluation_index, observation_index);
    np_beta_status scalar_status = NP_BETA_OK;
    double scalar_log = -INFINITY;
    int scalar_sign = 0;

    if(do_cdf) {
      const double scalar_value = np_beta_cdf_order(
        evaluation, observation, bandwidth,
        lower[dimension], upper[dimension], order, &scalar_status);
      if(scalar_status != NP_BETA_OK || !R_FINITE(scalar_value))
        return 1;
      if(scalar_value != 0.0) {
        scalar_log = log(fabs(scalar_value));
        scalar_sign = (scalar_value > 0.0) ? 1 : -1;
      }
    } else {
      scalar_log = np_beta_log_abs_pdf_order(
        evaluation, observation, bandwidth,
        lower[dimension], upper[dimension], order,
        &scalar_sign, &scalar_status);
      if(scalar_status != NP_BETA_OK || ISNAN(scalar_log) ||
         scalar_log == INFINITY)
        return 1;
    }

    if(scalar_sign == 0 || scalar_log == -INFINITY) {
      *log_absolute = -INFINITY;
      *sign = 0;
      return 0;
    }
    log_product += scalar_log;
    product_sign *= scalar_sign;
    if(!R_FINITE(log_product))
      return 1;
  }

  *log_absolute = log_product;
  *sign = product_sign;
  return 0;
}

static int np_beta_objective_kernel_sum_log(
  int do_cdf,
  int bandwidth_mode,
  int order,
  int num_train,
  int num_eval,
  int num_continuous,
  double **train_continuous,
  double **eval_continuous,
  double **matrix_bandwidth,
  const double *lower,
  const double *upper,
  int evaluation_index,
  int omitted_observation,
  const double *multiplier,
  double *log_absolute,
  int *sign)
{
  np_beta_log_sum sum;
  int observation_index;

  np_beta_log_sum_init(&sum);
  for(observation_index = 0; observation_index < num_train;
      ++observation_index) {
    double log_product = -INFINITY;
    int product_sign = 0;

    if(observation_index == omitted_observation)
      continue;
    if(np_beta_objective_product_log(
         do_cdf, bandwidth_mode, order,
         num_train, num_eval, num_continuous,
         train_continuous, eval_continuous, matrix_bandwidth,
         lower, upper, evaluation_index, observation_index,
         &log_product, &product_sign) != 0)
      return 1;
    if(multiplier != NULL) {
      const double scale = multiplier[observation_index];
      if(!R_FINITE(scale))
        return 1;
      if(scale == 0.0)
        product_sign = 0;
      else {
        log_product += log(fabs(scale));
        product_sign *= (scale > 0.0) ? 1 : -1;
      }
    }
    np_beta_log_sum_add(&sum, log_product, product_sign);
  }

  return np_beta_log_sum_finish(&sum, log_absolute, sign);
}

static double np_beta_objective_cvml_contribution(double log_sum,
                                                   int sign,
                                                   int denominator)
{
  const double log_dbl_min = log(DBL_MIN);
  const double log_fit = log_sum - log((double)denominator);

  if(sign > 0 && log_fit > log_dbl_min)
    return -log_fit;

  np_guarded_cvml_hit();
  if(sign < 0 && log_fit > log_dbl_min)
    return log_fit - 2.0 * log_dbl_min;
  return -log_dbl_min;
}

double np_beta_objective_density_ml(
  int bandwidth_mode,
  int order,
  int num_obs,
  int num_continuous,
  double **train_continuous,
  const double *candidate,
  const double *lower,
  const double *upper)
{
  double **matrix_bandwidth = NULL;
  double objective = 0.0;
  int bandwidth_rows;
  int evaluation_index;

  if(num_obs < 2 || num_continuous <= 0 || train_continuous == NULL ||
     candidate == NULL || lower == NULL || upper == NULL ||
     !np_beta_order_supported(order))
    return DBL_MAX;

  bandwidth_rows = np_beta_objective_bandwidth_rows(
    bandwidth_mode, num_obs, num_obs);
  if(bandwidth_rows <= 0)
    return DBL_MAX;
  matrix_bandwidth = alloc_matd(bandwidth_rows, num_continuous);
  if(matrix_bandwidth == NULL)
    return DBL_MAX;
  if(np_beta_objective_prepare_bandwidth(
       bandwidth_mode, num_obs, num_obs, num_continuous,
       train_continuous, train_continuous, candidate,
       matrix_bandwidth) != 0) {
    free_mat(matrix_bandwidth, num_continuous);
    return DBL_MAX;
  }

  for(evaluation_index = 0; evaluation_index < num_obs; ++evaluation_index) {
    double log_sum = -INFINITY;
    int sum_sign = 0;

    if((evaluation_index & 31) == 0)
      np_progress_bandwidth_loop_step();
    if(np_beta_objective_kernel_sum_log(
         0, bandwidth_mode, order,
         num_obs, num_obs, num_continuous,
         train_continuous, train_continuous, matrix_bandwidth,
         lower, upper, evaluation_index, evaluation_index,
         NULL, &log_sum, &sum_sign) != 0) {
      objective = DBL_MAX;
      break;
    }
    objective += np_beta_objective_cvml_contribution(
      log_sum, sum_sign, num_obs - 1);
    if(!R_FINITE(objective)) {
      objective = DBL_MAX;
      break;
    }
  }

  free_mat(matrix_bandwidth, num_continuous);
  return objective;
}

static int np_beta_objective_fill_grid(int num_continuous,
                                       const double *lower,
                                       const double *upper,
                                       int points,
                                       double **grid,
                                       double *weights)
{
  int total_points = 1;
  int evaluation_index;
  int dimension;

  for(dimension = 0; dimension < num_continuous; ++dimension) {
    if(!R_FINITE(lower[dimension]) || !R_FINITE(upper[dimension]) ||
       !(upper[dimension] > lower[dimension]) ||
       total_points > INT_MAX / points)
      return 1;
    total_points *= points;
  }

  for(evaluation_index = 0; evaluation_index < total_points;
      ++evaluation_index) {
    int index = evaluation_index;
    double weight = 1.0;

    for(dimension = 0; dimension < num_continuous; ++dimension) {
      const int grid_index = index % points;
      const double step = (upper[dimension] - lower[dimension]) /
        (double)(points - 1);
      grid[dimension][evaluation_index] = lower[dimension] +
        (double)grid_index * step;
      weight *= step * ((grid_index == 0 || grid_index == points - 1) ?
                        0.5 : 1.0);
      index /= points;
    }
    weights[evaluation_index] = weight;
  }
  return 0;
}

double np_beta_objective_density_ls(
  int bandwidth_mode,
  int order,
  int num_obs,
  int num_continuous,
  double **train_continuous,
  const double *candidate,
  const double *lower,
  const double *upper)
{
  const int points = (num_continuous == 1) ?
    NP_BETA_CVLS_GRID_1D : NP_BETA_CVLS_GRID_2D;
  const int total_points = (num_continuous == 1) ? points : points * points;
  double **grid = NULL;
  double **bandwidth_grid = NULL;
  double **bandwidth_train = NULL;
  double *weights = NULL;
  double integrated_square = 0.0;
  double cross_term = 0.0;
  double objective = DBL_MAX;
  int grid_bandwidth_rows;
  int train_bandwidth_rows;
  int evaluation_index;

  if(num_obs < 2 || num_continuous < 1 || num_continuous > 2 ||
     train_continuous == NULL || candidate == NULL ||
     lower == NULL || upper == NULL || !np_beta_order_supported(order))
    return DBL_MAX;

  grid = alloc_matd(total_points, num_continuous);
  weights = alloc_vecd(total_points);
  grid_bandwidth_rows = np_beta_objective_bandwidth_rows(
    bandwidth_mode, num_obs, total_points);
  train_bandwidth_rows = np_beta_objective_bandwidth_rows(
    bandwidth_mode, num_obs, num_obs);
  bandwidth_grid = alloc_matd(grid_bandwidth_rows, num_continuous);
  bandwidth_train = alloc_matd(train_bandwidth_rows, num_continuous);
  if(grid == NULL || weights == NULL || bandwidth_grid == NULL ||
     bandwidth_train == NULL)
    goto cleanup_density_ls;
  if(np_beta_objective_fill_grid(num_continuous, lower, upper,
                                 points, grid, weights) != 0)
    goto cleanup_density_ls;
  if(np_beta_objective_prepare_bandwidth(
       bandwidth_mode, num_obs, total_points, num_continuous,
       train_continuous, grid, candidate, bandwidth_grid) != 0 ||
     np_beta_objective_prepare_bandwidth(
       bandwidth_mode, num_obs, num_obs, num_continuous,
       train_continuous, train_continuous, candidate,
       bandwidth_train) != 0)
    goto cleanup_density_ls;

  for(evaluation_index = 0; evaluation_index < total_points;
      ++evaluation_index) {
    double log_sum = -INFINITY;
    double density = 0.0;
    int sum_sign = 0;

    if((evaluation_index & 31) == 0)
      np_progress_bandwidth_loop_step();
    if(np_beta_objective_kernel_sum_log(
         0, bandwidth_mode, order,
         num_obs, total_points, num_continuous,
         train_continuous, grid, bandwidth_grid,
         lower, upper, evaluation_index, -1,
         NULL, &log_sum, &sum_sign) != 0 ||
       np_beta_log_signed_value(log_sum - log((double)num_obs),
                                sum_sign, &density) != 0)
      goto cleanup_density_ls;
    integrated_square += weights[evaluation_index] * density * density;
    if(!R_FINITE(integrated_square))
      goto cleanup_density_ls;
  }

  for(evaluation_index = 0; evaluation_index < num_obs; ++evaluation_index) {
    double log_sum = -INFINITY;
    double density = 0.0;
    int sum_sign = 0;

    if(np_beta_objective_kernel_sum_log(
         0, bandwidth_mode, order,
         num_obs, num_obs, num_continuous,
         train_continuous, train_continuous, bandwidth_train,
         lower, upper, evaluation_index, evaluation_index,
         NULL, &log_sum, &sum_sign) != 0 ||
       np_beta_log_signed_value(log_sum - log((double)(num_obs - 1)),
                                sum_sign, &density) != 0)
      goto cleanup_density_ls;
    cross_term += density;
    if(!R_FINITE(cross_term))
      goto cleanup_density_ls;
  }

  objective = integrated_square - 2.0 * cross_term / (double)num_obs;
  if(!R_FINITE(objective))
    objective = DBL_MAX;

cleanup_density_ls:
  if(grid != NULL)
    free_mat(grid, num_continuous);
  if(bandwidth_grid != NULL)
    free_mat(bandwidth_grid, num_continuous);
  if(bandwidth_train != NULL)
    free_mat(bandwidth_train, num_continuous);
  if(weights != NULL)
    free(weights);
  return objective;
}

double np_beta_objective_distribution_ls(
  int bandwidth_mode,
  int order,
  int num_obs_train,
  int num_obs_eval,
  int num_continuous,
  int cdf_on_train,
  double **train_continuous,
  double **eval_continuous,
  const double *candidate,
  const double *lower,
  const double *upper)
{
  double **matrix_bandwidth = NULL;
  double *contribution = NULL;
  double objective = 0.0;
  int bandwidth_rows;
  int evaluation_index;

  if(num_obs_train < 2 || num_obs_eval <= 0 || num_continuous <= 0 ||
     train_continuous == NULL || eval_continuous == NULL ||
     candidate == NULL || lower == NULL || upper == NULL ||
     !np_beta_order_supported(order))
    return DBL_MAX;

  bandwidth_rows = np_beta_objective_bandwidth_rows(
    bandwidth_mode, num_obs_train, num_obs_eval);
  matrix_bandwidth = alloc_matd(bandwidth_rows, num_continuous);
  contribution = alloc_vecd(num_obs_train);
  if(matrix_bandwidth == NULL || contribution == NULL)
    goto distribution_failure;
  if(np_beta_objective_prepare_bandwidth(
       bandwidth_mode, num_obs_train, num_obs_eval, num_continuous,
       train_continuous, eval_continuous, candidate,
       matrix_bandwidth) != 0)
    goto distribution_failure;

  for(evaluation_index = 0; evaluation_index < num_obs_eval;
      ++evaluation_index) {
    double total = 0.0;
    double correction = 0.0;
    int observation_index;

    if((evaluation_index & 31) == 0)
      np_progress_bandwidth_loop_step();
    for(observation_index = 0; observation_index < num_obs_train;
        ++observation_index) {
      double log_product = -INFINITY;
      double value = 0.0;
      double adjusted;
      double temporary;
      int product_sign = 0;

      if(np_beta_objective_product_log(
           1, bandwidth_mode, order,
           num_obs_train, num_obs_eval, num_continuous,
           train_continuous, eval_continuous, matrix_bandwidth,
           lower, upper, evaluation_index, observation_index,
           &log_product, &product_sign) != 0 ||
         np_beta_log_signed_value(log_product, product_sign, &value) != 0)
        goto distribution_failure;
      contribution[observation_index] = value;
      adjusted = value - correction;
      temporary = total + adjusted;
      correction = (temporary - total) - adjusted;
      total = temporary;
    }

    for(observation_index = 0; observation_index < num_obs_train;
        ++observation_index) {
      double indicator = 1.0;
      const double fitted = (total - contribution[observation_index]) /
        (double)(num_obs_train - 1);
      double difference;
      int dimension;

      if(cdf_on_train && observation_index == evaluation_index)
        continue;

      for(dimension = 0; dimension < num_continuous; ++dimension) {
        if(train_continuous[dimension][observation_index] >
           eval_continuous[dimension][evaluation_index]) {
          indicator = 0.0;
          break;
        }
      }
      difference = indicator - fitted;
      objective += difference * difference;
      if(!R_FINITE(objective))
        goto distribution_failure;
    }
  }

  objective /= (double)num_obs_train * (double)num_obs_eval;
  free_mat(matrix_bandwidth, num_continuous);
  free(contribution);
  return R_FINITE(objective) ? objective : DBL_MAX;

distribution_failure:
  if(matrix_bandwidth != NULL)
    free_mat(matrix_bandwidth, num_continuous);
  if(contribution != NULL)
    free(contribution);
  return DBL_MAX;
}

static int np_beta_objective_ratio(double numerator_log,
                                   int numerator_sign,
                                   double denominator_log,
                                   int denominator_sign,
                                   double *ratio)
{
  if(ratio == NULL || denominator_sign == 0 ||
     denominator_log == -INFINITY)
    return 1;
  if(numerator_sign == 0 || numerator_log == -INFINITY) {
    *ratio = 0.0;
    return 0;
  }
  return np_beta_log_signed_value(
    numerator_log - denominator_log,
    numerator_sign * denominator_sign,
    ratio);
}

double np_beta_objective_regression_lc(
  int bandwidth_mode,
  int order,
  int use_aic,
  int num_obs,
  int num_continuous,
  double **train_continuous,
  const double *response,
  const double *candidate,
  const double *lower,
  const double *upper)
{
  double **matrix_bandwidth = NULL;
  double loss = 0.0;
  double trace_hat = 0.0;
  double objective = DBL_MAX;
  int bandwidth_rows;
  int evaluation_index;

  if(num_obs < 2 || num_continuous <= 0 || train_continuous == NULL ||
     response == NULL || candidate == NULL || lower == NULL || upper == NULL ||
     !np_beta_order_supported(order))
    return DBL_MAX;

  bandwidth_rows = np_beta_objective_bandwidth_rows(
    bandwidth_mode, num_obs, num_obs);
  matrix_bandwidth = alloc_matd(bandwidth_rows, num_continuous);
  if(matrix_bandwidth == NULL)
    return DBL_MAX;
  if(np_beta_objective_prepare_bandwidth(
       bandwidth_mode, num_obs, num_obs, num_continuous,
       train_continuous, train_continuous, candidate,
       matrix_bandwidth) != 0)
    goto cleanup_regression;

  for(evaluation_index = 0; evaluation_index < num_obs; ++evaluation_index) {
    const int omitted = use_aic ? -1 : evaluation_index;
    double denominator_log = -INFINITY;
    double numerator_log = -INFINITY;
    double fitted = 0.0;
    double residual;
    int denominator_sign = 0;
    int numerator_sign = 0;

    if((evaluation_index & 31) == 0)
      np_progress_bandwidth_loop_step();
    if(np_beta_objective_kernel_sum_log(
         0, bandwidth_mode, order,
         num_obs, num_obs, num_continuous,
         train_continuous, train_continuous, matrix_bandwidth,
         lower, upper, evaluation_index, omitted,
         NULL, &denominator_log, &denominator_sign) != 0 ||
       np_beta_objective_kernel_sum_log(
         0, bandwidth_mode, order,
         num_obs, num_obs, num_continuous,
         train_continuous, train_continuous, matrix_bandwidth,
         lower, upper, evaluation_index, omitted,
         response, &numerator_log, &numerator_sign) != 0 ||
       np_beta_objective_ratio(numerator_log, numerator_sign,
                               denominator_log, denominator_sign,
                               &fitted) != 0)
      goto cleanup_regression;

    residual = response[evaluation_index] - fitted;
    loss += residual * residual;
    if(!R_FINITE(loss))
      goto cleanup_regression;

    if(use_aic) {
      double self_log = -INFINITY;
      double diagonal = 0.0;
      int self_sign = 0;

      if(np_beta_objective_product_log(
           0, bandwidth_mode, order,
           num_obs, num_obs, num_continuous,
           train_continuous, train_continuous, matrix_bandwidth,
           lower, upper, evaluation_index, evaluation_index,
           &self_log, &self_sign) != 0 ||
         np_beta_objective_ratio(self_log, self_sign,
                                 denominator_log, denominator_sign,
                                 &diagonal) != 0)
        goto cleanup_regression;
      trace_hat += diagonal;
      if(!R_FINITE(trace_hat))
        goto cleanup_regression;
    }
  }

  loss /= (double)num_obs;
  if(!use_aic) {
    objective = loss;
  } else {
    const double denominator = 1.0 - (trace_hat + 2.0) / (double)num_obs;
    const double numerator = 1.0 + trace_hat / (double)num_obs;
    const double penalty = numerator / denominator;

    if(loss > 0.0 && R_FINITE(loss) && R_FINITE(penalty))
      objective = log(loss) + penalty;
  }

cleanup_regression:
  free_mat(matrix_bandwidth, num_continuous);
  return R_FINITE(objective) ? objective : DBL_MAX;
}

static np_continuous_kernel_family
np_beta_objective_conditional_family(int kernel_code)
{
  return (kernel_code == NP_CKERNEL_COORDINATE_CODE) ?
    NP_CKERNEL_FAMILY_BETA : NP_CKERNEL_FAMILY_LEGACY;
}

static int np_beta_objective_conditional_product_log(
  int do_cdf,
  int bandwidth_mode,
  int kernel_code,
  int order,
  int num_dimensions,
  double **train,
  double **evaluation,
  double **matrix_bandwidth,
  const double *lower,
  const double *upper,
  int evaluation_index,
  int observation_index,
  double *log_absolute,
  int *sign)
{
  const np_continuous_kernel_family family =
    np_beta_objective_conditional_family(kernel_code);
  double log_product = 0.0;
  int product_sign = 1;
  int dimension;

  if(num_dimensions <= 0 || train == NULL || evaluation == NULL ||
     matrix_bandwidth == NULL || lower == NULL || upper == NULL ||
     log_absolute == NULL || sign == NULL)
    return 1;

  for(dimension = 0; dimension < num_dimensions; ++dimension) {
    const double bandwidth = np_beta_objective_bandwidth(
      bandwidth_mode, matrix_bandwidth, dimension,
      evaluation_index, observation_index);
    double scalar_log = -INFINITY;
    int scalar_sign = 0;
    np_beta_status beta_status = NP_BETA_OK;

    if(np_beta_conditional_scalar_log(
         family, kernel_code, order, do_cdf,
         evaluation[dimension][evaluation_index],
         train[dimension][observation_index], bandwidth,
         lower[dimension], upper[dimension],
         &scalar_log, &scalar_sign, &beta_status) !=
       NP_BETA_CONDITIONAL_OK)
      return 1;
    if(scalar_sign == 0 || scalar_log == -INFINITY) {
      *log_absolute = -INFINITY;
      *sign = 0;
      return 0;
    }
    log_product += scalar_log;
    product_sign *= scalar_sign;
    if(!R_FINITE(log_product))
      return 1;
  }

  *log_absolute = log_product;
  *sign = product_sign;
  return 0;
}

static int np_beta_objective_conditional_sums(
  int do_cdf_y,
  int bandwidth_mode,
  int kernel_code_x,
  int order_x,
  int kernel_code_y,
  int order_y,
  int num_obs,
  int num_x,
  int num_y,
  double **train_x,
  double **train_y,
  double **eval_x,
  double **eval_y,
  double **bandwidth_x,
  double **bandwidth_y,
  const double *lower_x,
  const double *upper_x,
  const double *lower_y,
  const double *upper_y,
  int evaluation_x,
  int evaluation_y,
  int omitted_observation,
  double *denominator_log,
  int *denominator_sign,
  double *numerator_log,
  int *numerator_sign)
{
  np_beta_log_sum denominator;
  np_beta_log_sum numerator;
  int observation_index;

  np_beta_log_sum_init(&denominator);
  np_beta_log_sum_init(&numerator);
  for(observation_index = 0; observation_index < num_obs;
      ++observation_index) {
    double x_log = -INFINITY;
    double y_log = -INFINITY;
    int x_sign = 0;
    int y_sign = 0;

    if(observation_index == omitted_observation)
      continue;
    if(np_beta_objective_conditional_product_log(
         0, bandwidth_mode, kernel_code_x, order_x, num_x,
         train_x, eval_x, bandwidth_x, lower_x, upper_x,
         evaluation_x, observation_index, &x_log, &x_sign) != 0)
      return 1;
    np_beta_log_sum_add(&denominator, x_log, x_sign);
    if(x_sign == 0)
      continue;
    if(np_beta_objective_conditional_product_log(
         do_cdf_y, bandwidth_mode, kernel_code_y, order_y, num_y,
         train_y, eval_y, bandwidth_y, lower_y, upper_y,
         evaluation_y, observation_index, &y_log, &y_sign) != 0)
      return 1;
    np_beta_log_sum_add(&numerator, x_log + y_log, x_sign * y_sign);
  }

  if(np_beta_log_sum_finish(&denominator,
                            denominator_log, denominator_sign) != 0 ||
     np_beta_log_sum_finish(&numerator,
                            numerator_log, numerator_sign) != 0 ||
     *denominator_sign == 0)
    return 1;
  return 0;
}

static int np_beta_objective_conditional_layout_ok(int bandwidth_mode,
                                                    int order_x,
                                                    int order_y,
                                                    int num_obs,
                                                    int num_x,
                                                    int num_y,
                                                    double **train_x,
                                                    double **train_y,
                                                    const double *candidate,
                                                    const double *lower_x,
                                                    const double *upper_x,
                                                    const double *lower_y,
                                                    const double *upper_y)
{
  return num_obs >= 2 && num_x > 0 && num_y > 0 &&
    train_x != NULL && train_y != NULL && candidate != NULL &&
    lower_x != NULL && upper_x != NULL &&
    lower_y != NULL && upper_y != NULL &&
    np_beta_order_supported(order_x) && np_beta_order_supported(order_y) &&
    np_beta_objective_bandwidth_rows(bandwidth_mode, num_obs, num_obs) > 0;
}

double np_beta_objective_conditional_density_ml(
  int bandwidth_mode,
  int kernel_code_x,
  int order_x,
  int kernel_code_y,
  int order_y,
  int num_obs,
  int num_x,
  int num_y,
  double **train_x,
  double **train_y,
  const double *candidate,
  const double *lower_x,
  const double *upper_x,
  const double *lower_y,
  const double *upper_y)
{
  double **bandwidth_x = NULL;
  double **bandwidth_y = NULL;
  double objective = DBL_MAX;
  double sum = 0.0;
  int bandwidth_rows;
  int evaluation_index;

  if(!np_beta_objective_conditional_layout_ok(
       bandwidth_mode, order_x, order_y, num_obs, num_x, num_y,
       train_x, train_y, candidate, lower_x, upper_x, lower_y, upper_y))
    return DBL_MAX;
  bandwidth_rows = np_beta_objective_bandwidth_rows(
    bandwidth_mode, num_obs, num_obs);
  bandwidth_x = alloc_matd(bandwidth_rows, num_x);
  bandwidth_y = alloc_matd(bandwidth_rows, num_y);
  if(bandwidth_x == NULL || bandwidth_y == NULL)
    goto cleanup_conditional_density_ml;
  if(np_beta_objective_prepare_bandwidth_offset(
       bandwidth_mode, num_obs, num_obs, num_x,
       0, train_x, train_x, candidate, bandwidth_x) != 0 ||
     np_beta_objective_prepare_bandwidth_offset(
       bandwidth_mode, num_obs, num_obs, num_y, num_x,
       train_y, train_y, candidate + num_x, bandwidth_y) != 0)
    goto cleanup_conditional_density_ml;

  for(evaluation_index = 0; evaluation_index < num_obs; ++evaluation_index) {
    double denominator_log = -INFINITY;
    double numerator_log = -INFINITY;
    int denominator_sign = 0;
    int numerator_sign = 0;
    int fit_sign;

    if((evaluation_index & 31) == 0)
      np_progress_bandwidth_loop_step();
    if(np_beta_objective_conditional_sums(
         0, bandwidth_mode, kernel_code_x, order_x,
         kernel_code_y, order_y, num_obs, num_x, num_y,
         train_x, train_y, train_x, train_y,
         bandwidth_x, bandwidth_y,
         lower_x, upper_x, lower_y, upper_y,
         evaluation_index, evaluation_index, evaluation_index,
         &denominator_log, &denominator_sign,
         &numerator_log, &numerator_sign) != 0)
      goto cleanup_conditional_density_ml;
    fit_sign = denominator_sign * numerator_sign;
    sum += np_beta_objective_cvml_contribution(
      numerator_log - denominator_log, fit_sign, 1);
    if(!R_FINITE(sum))
      goto cleanup_conditional_density_ml;
  }
  objective = sum;

cleanup_conditional_density_ml:
  if(bandwidth_x != NULL)
    free_mat(bandwidth_x, num_x);
  if(bandwidth_y != NULL)
    free_mat(bandwidth_y, num_y);
  return R_FINITE(objective) ? objective : DBL_MAX;
}

double np_beta_objective_conditional_distribution_ls(
  int bandwidth_mode,
  int kernel_code_x,
  int order_x,
  int kernel_code_y,
  int order_y,
  int num_obs_train,
  int num_obs_eval,
  int num_x,
  int num_y,
  int cdf_on_train,
  double **train_x,
  double **train_y,
  double **eval_y,
  const double *candidate,
  const double *lower_x,
  const double *upper_x,
  const double *lower_y,
  const double *upper_y)
{
  double **bandwidth_x = NULL;
  double **bandwidth_y = NULL;
  double objective = DBL_MAX;
  double sum = 0.0;
  int bandwidth_rows_x;
  int bandwidth_rows_y;
  int observation_index;
  int evaluation_index;

  if(!np_beta_objective_conditional_layout_ok(
       bandwidth_mode, order_x, order_y, num_obs_train, num_x, num_y,
       train_x, train_y, candidate, lower_x, upper_x, lower_y, upper_y) ||
     num_obs_eval <= 0 || eval_y == NULL)
    return DBL_MAX;
  bandwidth_rows_x = np_beta_objective_bandwidth_rows(
    bandwidth_mode, num_obs_train, num_obs_train);
  bandwidth_rows_y = np_beta_objective_bandwidth_rows(
    bandwidth_mode, num_obs_train, num_obs_eval);
  bandwidth_x = alloc_matd(bandwidth_rows_x, num_x);
  bandwidth_y = alloc_matd(bandwidth_rows_y, num_y);
  if(bandwidth_x == NULL || bandwidth_y == NULL)
    goto cleanup_conditional_distribution_ls;
  if(np_beta_objective_prepare_bandwidth_offset(
       bandwidth_mode, num_obs_train, num_obs_train, num_x,
       0, train_x, train_x, candidate, bandwidth_x) != 0 ||
     np_beta_objective_prepare_bandwidth_offset(
       bandwidth_mode, num_obs_train, num_obs_eval, num_y, num_x,
       train_y, eval_y, candidate + num_x, bandwidth_y) != 0)
    goto cleanup_conditional_distribution_ls;

  for(observation_index = 0; observation_index < num_obs_train;
      ++observation_index) {
    if((observation_index & 15) == 0)
      np_progress_bandwidth_loop_step();
    for(evaluation_index = 0; evaluation_index < num_obs_eval;
        ++evaluation_index) {
      double denominator_log = -INFINITY;
      double numerator_log = -INFINITY;
      double fitted = 0.0;
      double indicator = 1.0;
      double difference;
      int denominator_sign = 0;
      int numerator_sign = 0;
      int dimension;

      if(cdf_on_train && evaluation_index == observation_index)
        continue;
      if(np_beta_objective_conditional_sums(
           1, bandwidth_mode, kernel_code_x, order_x,
           kernel_code_y, order_y, num_obs_train, num_x, num_y,
           train_x, train_y, train_x, eval_y,
           bandwidth_x, bandwidth_y,
           lower_x, upper_x, lower_y, upper_y,
           observation_index, evaluation_index, observation_index,
           &denominator_log, &denominator_sign,
           &numerator_log, &numerator_sign) != 0 ||
         np_beta_objective_ratio(
           numerator_log, numerator_sign,
           denominator_log, denominator_sign, &fitted) != 0)
        goto cleanup_conditional_distribution_ls;
      for(dimension = 0; dimension < num_y; ++dimension)
        indicator *= train_y[dimension][observation_index] <=
          eval_y[dimension][evaluation_index] ? 1.0 : 0.0;
      difference = indicator - fitted;
      sum += difference * difference;
      if(!R_FINITE(sum))
        goto cleanup_conditional_distribution_ls;
    }
  }
  objective = sum / ((double)num_obs_train * (double)num_obs_eval);

cleanup_conditional_distribution_ls:
  if(bandwidth_x != NULL)
    free_mat(bandwidth_x, num_x);
  if(bandwidth_y != NULL)
    free_mat(bandwidth_y, num_y);
  return R_FINITE(objective) ? objective : DBL_MAX;
}

double np_beta_objective_conditional_density_ls(
  int bandwidth_mode,
  int kernel_code_x,
  int order_x,
  int kernel_code_y,
  int order_y,
  int quadrature_points,
  int num_obs,
  int num_x,
  int num_y,
  double **train_x,
  double **train_y,
  const double *candidate,
  const double *lower_x,
  const double *upper_x,
  const double *lower_y,
  const double *upper_y)
{
  double **grid_y = NULL;
  double **bandwidth_x = NULL;
  double **bandwidth_y_train = NULL;
  double **bandwidth_y_grid = NULL;
  double *weights = NULL;
  double objective = DBL_MAX;
  double sum = 0.0;
  int total_points = 1;
  int bandwidth_rows_train;
  int bandwidth_rows_grid;
  int actual_points;
  int observation_index;
  int dimension;

  if(!np_beta_objective_conditional_layout_ok(
       bandwidth_mode, order_x, order_y, num_obs, num_x, num_y,
       train_x, train_y, candidate, lower_x, upper_x, lower_y, upper_y) ||
     num_y > 2 || quadrature_points < 2)
    return DBL_MAX;
  for(dimension = 0; dimension < num_y; ++dimension) {
    if(total_points > INT_MAX / quadrature_points)
      return DBL_MAX;
    total_points *= quadrature_points;
  }

  bandwidth_rows_train = np_beta_objective_bandwidth_rows(
    bandwidth_mode, num_obs, num_obs);
  grid_y = alloc_matd(total_points, num_y);
  weights = alloc_vecd(total_points);
  bandwidth_x = alloc_matd(bandwidth_rows_train, num_x);
  bandwidth_y_train = alloc_matd(bandwidth_rows_train, num_y);
  if(grid_y == NULL || weights == NULL || bandwidth_x == NULL ||
     bandwidth_y_train == NULL)
    goto cleanup_conditional_density_ls;
  actual_points = total_points;
  if(((num_y == 1) ?
       np_bounded_cvls_build_conditional_grid_1d_extern(
         train_y[0], num_obs, lower_y[0], upper_y[0],
         quadrature_points, grid_y[0], weights, &actual_points) :
       np_beta_objective_fill_grid(num_y, lower_y, upper_y,
                                   quadrature_points, grid_y, weights)) != 0)
    goto cleanup_conditional_density_ls;

  total_points = actual_points;
  bandwidth_rows_grid = np_beta_objective_bandwidth_rows(
    bandwidth_mode, num_obs, total_points);
  bandwidth_y_grid = alloc_matd(bandwidth_rows_grid, num_y);
  if(bandwidth_y_grid == NULL ||
     np_beta_objective_prepare_bandwidth_offset(
       bandwidth_mode, num_obs, num_obs, num_x,
       0, train_x, train_x, candidate, bandwidth_x) != 0 ||
     np_beta_objective_prepare_bandwidth_offset(
       bandwidth_mode, num_obs, num_obs, num_y, num_x,
       train_y, train_y, candidate + num_x, bandwidth_y_train) != 0 ||
     np_beta_objective_prepare_bandwidth_offset(
       bandwidth_mode, num_obs, total_points, num_y, num_x,
       train_y, grid_y, candidate + num_x, bandwidth_y_grid) != 0)
    goto cleanup_conditional_density_ls;
  for(observation_index = 0; observation_index < num_obs;
      ++observation_index) {
    double denominator_log = -INFINITY;
    double cross_numerator_log = -INFINITY;
    double cross_fit = 0.0;
    double integrated_square = 0.0;
    int denominator_sign = 0;
    int cross_numerator_sign = 0;
    int grid_index;

    if((observation_index & 7) == 0)
      np_progress_bandwidth_loop_step();
    if(np_beta_objective_conditional_sums(
         0, bandwidth_mode, kernel_code_x, order_x,
         kernel_code_y, order_y, num_obs, num_x, num_y,
         train_x, train_y, train_x, train_y,
         bandwidth_x, bandwidth_y_train,
         lower_x, upper_x, lower_y, upper_y,
         observation_index, observation_index, observation_index,
         &denominator_log, &denominator_sign,
         &cross_numerator_log, &cross_numerator_sign) != 0 ||
       np_beta_objective_ratio(
         cross_numerator_log, cross_numerator_sign,
         denominator_log, denominator_sign, &cross_fit) != 0)
      goto cleanup_conditional_density_ls;

    for(grid_index = 0; grid_index < total_points; ++grid_index) {
      double ignored_denominator_log = -INFINITY;
      double grid_numerator_log = -INFINITY;
      double fitted = 0.0;
      int ignored_denominator_sign = 0;
      int grid_numerator_sign = 0;

      if(np_beta_objective_conditional_sums(
           0, bandwidth_mode, kernel_code_x, order_x,
           kernel_code_y, order_y, num_obs, num_x, num_y,
           train_x, train_y, train_x, grid_y,
           bandwidth_x, bandwidth_y_grid,
           lower_x, upper_x, lower_y, upper_y,
           observation_index, grid_index, observation_index,
           &ignored_denominator_log, &ignored_denominator_sign,
           &grid_numerator_log, &grid_numerator_sign) != 0 ||
         np_beta_objective_ratio(
           grid_numerator_log, grid_numerator_sign,
           denominator_log, denominator_sign, &fitted) != 0)
        goto cleanup_conditional_density_ls;
      integrated_square += weights[grid_index] * fitted * fitted;
      if(!R_FINITE(integrated_square))
        goto cleanup_conditional_density_ls;
    }
    sum += integrated_square - 2.0 * cross_fit;
    if(!R_FINITE(sum))
      goto cleanup_conditional_density_ls;
  }
  objective = sum / (double)num_obs;

cleanup_conditional_density_ls:
  if(grid_y != NULL)
    free_mat(grid_y, num_y);
  if(bandwidth_x != NULL)
    free_mat(bandwidth_x, num_x);
  if(bandwidth_y_train != NULL)
    free_mat(bandwidth_y_train, num_y);
  if(bandwidth_y_grid != NULL)
    free_mat(bandwidth_y_grid, num_y);
  if(weights != NULL)
    free(weights);
  return R_FINITE(objective) ? objective : DBL_MAX;
}
