/*
 * Associated beta-kernel primitives for bounded continuous variables.
 *
 * The public metric bandwidth h is converted to the dimensionless
 * concentration tau = ((upper - lower) / h)^2.  The order-2 PDF is
 *
 *   Beta(observation_unit;
 *        1 + target_unit * tau,
 *        1 + (1 - target_unit) * tau) / (upper - lower).
 *
 * The complete original-scale normalization is returned here.  Callers must
 * not apply either a legacy external 1/h factor or generic bounded-kernel
 * truncation normalization.
 */

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include <R_ext/Arith.h>
#include <Rmath.h>

#include "beta_kernel.h"

static void np_beta_set_status(np_beta_status *status,
                               np_beta_status value)
{
  if(status != NULL)
    *status = value;
}

static double np_beta_unit_coordinate(double value,
                                      double lower,
                                      double upper,
                                      double support_length)
{
  const double midpoint = lower + 0.5 * support_length;

  /* Use the nearer endpoint to avoid losing most of a small upper-tail gap. */
  if(value <= midpoint)
    return (value - lower) / support_length;

  return 1.0 - (upper - value) / support_length;
}

/* Stable log-gamma differences used by the analytic beta-product overlap.
 * Directly subtracting three lbeta values loses several digits when the
 * concentration is large because the desired log overlap is small relative
 * to the individual log beta functions. */
static double np_beta_log_gamma_backward_half(double x)
{
  double power = -0.5;
  double factorial = 1.0;
  double result = 0.0;
  int order;

  if(!R_FINITE(x) || x < 1.0)
    return NAN;
  if(x < 32.0)
    return lgammafn(x - 0.5) - lgammafn(x);

  /* Taylor-expand log Gamma(x-1/2) about x. */
  for(order = 1; order <= 8; ++order) {
    const double derivative = psigamma(x, (double)(order - 1));
    factorial *= (double)order;
    if(!R_FINITE(derivative))
      return NAN;
    result += power * derivative / factorial;
    power *= -0.5;
  }
  return result;
}

static double np_beta_log_gamma_centered(double x, double displacement)
{
  const double absolute_displacement = fabs(displacement);

  if(!R_FINITE(x) || !R_FINITE(displacement) ||
     x < 1.0 || absolute_displacement >= x)
    return NAN;
  if(displacement == 0.0)
    return 0.0;

  if(absolute_displacement <= 0.05 * x) {
    double squared_power = displacement * displacement;
    double factorial = 2.0;
    double result = 0.0;
    int term;

    /* f(x+d)+f(x-d)-2f(x), f=log Gamma. Odd powers cancel. */
    for(term = 1; term <= 6; ++term) {
      const int derivative_order = 2 * term - 1;
      const double derivative = psigamma(x, (double)derivative_order);
      if(!R_FINITE(derivative))
        return NAN;
      result += 2.0 * squared_power * derivative / factorial;
      squared_power *= displacement * displacement;
      factorial *= (double)(2 * term + 1) * (double)(2 * term + 2);
    }
    return result;
  }

  return lgammafn(x + displacement) + lgammafn(x - displacement) -
    2.0 * lgammafn(x);
}

static int np_beta_order_coefficients(int order,
                                      const int **coefficients)
{
  static const int order_two[] = {1};
  static const int order_four[] = {2, -1};
  static const int order_six[] = {3, -3, 1};
  static const int order_eight[] = {4, -6, 4, -1};

  if(coefficients == NULL)
    return 0;

  switch(order) {
  case 2:
    *coefficients = order_two;
    return 1;
  case 4:
    *coefficients = order_four;
    return 2;
  case 6:
    *coefficients = order_six;
    return 3;
  case 8:
    *coefficients = order_eight;
    return 4;
  default:
    *coefficients = NULL;
    return 0;
  }
}

int np_beta_order_supported(int order)
{
  const int *coefficients = NULL;
  return np_beta_order_coefficients(order, &coefficients) > 0;
}

static double np_beta_log_add(double accumulator, double term)
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

np_beta_status np_beta_signed_log_absolute(double positive_log,
                                           double negative_log,
                                           double *log_absolute,
                                           int *sign)
{
  double scaled_difference;

  if(log_absolute == NULL || sign == NULL ||
     ISNAN(positive_log) || ISNAN(negative_log))
    return NP_BETA_ERR_NUMERIC;
  if(positive_log == -INFINITY && negative_log == -INFINITY) {
    *log_absolute = -INFINITY;
    *sign = 0;
    return NP_BETA_OK;
  }
  if(negative_log == -INFINITY) {
    *log_absolute = positive_log;
    *sign = 1;
  } else if(positive_log == -INFINITY) {
    *log_absolute = negative_log;
    *sign = -1;
  } else if(positive_log == negative_log) {
    *log_absolute = -INFINITY;
    *sign = 0;
    return NP_BETA_OK;
  } else if(positive_log > negative_log) {
    scaled_difference = -expm1(negative_log - positive_log);
    if(!R_FINITE(scaled_difference) || scaled_difference <= 0.0)
      return NP_BETA_ERR_NUMERIC;
    *log_absolute = positive_log + log(scaled_difference);
    *sign = 1;
  } else {
    scaled_difference = -expm1(positive_log - negative_log);
    if(!R_FINITE(scaled_difference) || scaled_difference <= 0.0)
      return NP_BETA_ERR_NUMERIC;
    *log_absolute = negative_log + log(scaled_difference);
    *sign = -1;
  }

  if(*log_absolute == -INFINITY) {
    *sign = 0;
    return NP_BETA_OK;
  }
  if(!R_FINITE(*log_absolute))
    return NP_BETA_ERR_RANGE;

  return NP_BETA_OK;
}

static double np_beta_signed_log_difference(double positive_log,
                                            double negative_log,
                                            np_beta_status *status)
{
  double log_absolute = -INFINITY;
  double value;
  int sign = 0;
  const np_beta_status difference_status = np_beta_signed_log_absolute(
    positive_log, negative_log, &log_absolute, &sign);

  if(difference_status != NP_BETA_OK) {
    np_beta_set_status(status, difference_status);
    return (difference_status == NP_BETA_ERR_RANGE && sign != 0) ?
      ((sign > 0) ? INFINITY : -INFINITY) : NAN;
  }
  if(sign == 0) {
    np_beta_set_status(status, NP_BETA_OK);
    return 0.0;
  }
  if(log_absolute > log(DBL_MAX)) {
    np_beta_set_status(status, NP_BETA_ERR_RANGE);
    return (sign > 0) ? INFINITY : -INFINITY;
  }

  value = exp(log_absolute);
  if(!R_FINITE(value)) {
    np_beta_set_status(status, NP_BETA_ERR_RANGE);
    return (sign > 0) ? INFINITY : -INFINITY;
  }

  np_beta_set_status(status, NP_BETA_OK);
  return (sign > 0) ? value : -value;
}

const char *np_beta_status_message(np_beta_status status)
{
  switch(status) {
  case NP_BETA_OK:
    return "success";
  case NP_BETA_ERR_NONFINITE:
    return "beta-kernel arguments must be finite";
  case NP_BETA_ERR_BOUNDS:
    return "beta-kernel bounds must satisfy lower < upper with finite width";
  case NP_BETA_ERR_BANDWIDTH:
    return "beta-kernel bandwidth must be finite and strictly positive";
  case NP_BETA_ERR_SCALE:
    return "beta-kernel concentration scale must be a positive integer";
  case NP_BETA_ERR_OBSERVATION:
    return "beta-kernel observations must lie within the declared bounds";
  case NP_BETA_ERR_CONCENTRATION:
    return "beta-kernel concentration is outside the representable range";
  case NP_BETA_ERR_NUMERIC:
    return "beta-kernel evaluation produced an undefined numeric result";
  case NP_BETA_ERR_RANGE:
    return "beta-kernel value exceeds the representable numeric range";
  default:
    return "unknown beta-kernel status";
  }
}

np_beta_status np_beta_shape_init(double evaluation,
                                  double observation,
                                  double bandwidth,
                                  double lower,
                                  double upper,
                                  int scale,
                                  np_beta_shape *shape)
{
  double concentration_limit;
  double ratio;

  if(shape == NULL)
    return NP_BETA_ERR_NUMERIC;

  if(!R_FINITE(evaluation) || !R_FINITE(observation) ||
     !R_FINITE(lower) || !R_FINITE(upper))
    return NP_BETA_ERR_NONFINITE;

  if(!R_FINITE(bandwidth) || bandwidth <= 0.0)
    return NP_BETA_ERR_BANDWIDTH;

  if(scale <= 0)
    return NP_BETA_ERR_SCALE;

  shape->lower = lower;
  shape->upper = upper;
  shape->support_length = upper - lower;
  if(!R_FINITE(shape->support_length) || shape->support_length <= 0.0)
    return NP_BETA_ERR_BOUNDS;

  if(observation < lower || observation > upper)
    return NP_BETA_ERR_OBSERVATION;
  shape->observation_unit = np_beta_unit_coordinate(observation,
                                                    lower,
                                                    upper,
                                                    shape->support_length);
  shape->observation_complement_unit =
    (upper - observation) / shape->support_length;
  if(!R_FINITE(shape->observation_unit) ||
     shape->observation_unit < 0.0 || shape->observation_unit > 1.0)
    return NP_BETA_ERR_OBSERVATION;
  if(!R_FINITE(shape->observation_complement_unit) ||
     shape->observation_complement_unit < 0.0 ||
     shape->observation_complement_unit > 1.0)
    return NP_BETA_ERR_OBSERVATION;

  if(observation == lower) {
    shape->observation_endpoint = -1;
    shape->log_observation_unit = -INFINITY;
    shape->log_observation_complement_unit = 0.0;
  } else if(observation == upper) {
    shape->observation_endpoint = 1;
    shape->log_observation_unit = 0.0;
    shape->log_observation_complement_unit = -INFINITY;
  } else {
    shape->observation_endpoint = 0;
    shape->log_observation_unit =
      log(observation - lower) - log(shape->support_length);
    shape->log_observation_complement_unit =
      log(upper - observation) - log(shape->support_length);
  }

  if(evaluation < lower) {
    shape->eval_location = NP_BETA_EVAL_BELOW;
    shape->target_unit = 0.0;
    shape->target_complement_unit = 1.0;
  } else if(evaluation > upper) {
    shape->eval_location = NP_BETA_EVAL_ABOVE;
    shape->target_unit = 1.0;
    shape->target_complement_unit = 0.0;
  } else {
    shape->eval_location = NP_BETA_EVAL_INSIDE;
    shape->target_unit = (evaluation - lower) / shape->support_length;
    shape->target_complement_unit =
      (upper - evaluation) / shape->support_length;
    if(!R_FINITE(shape->target_unit) ||
       shape->target_unit < 0.0 || shape->target_unit > 1.0)
      return NP_BETA_ERR_NUMERIC;
    if(!R_FINITE(shape->target_complement_unit) ||
       shape->target_complement_unit < 0.0 ||
       shape->target_complement_unit > 1.0)
      return NP_BETA_ERR_NUMERIC;
  }

  ratio = shape->support_length / bandwidth;
  concentration_limit = sqrt(DBL_MAX) * sqrt((double) scale);
  if(!R_FINITE(ratio) || !R_FINITE(concentration_limit) ||
     ratio > concentration_limit)
    return NP_BETA_ERR_CONCENTRATION;

  shape->concentration = (ratio / (double) scale) * ratio;
  if(!R_FINITE(shape->concentration))
    return NP_BETA_ERR_CONCENTRATION;

  return NP_BETA_OK;
}

double np_beta_log_pdf_scale(const np_beta_shape *shape,
                             np_beta_status *status)
{
  double alpha;
  double beta;
  double alpha_minus_one;
  double beta_minus_one;
  double log_value;

  if(shape == NULL) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }

  if(shape->eval_location != NP_BETA_EVAL_INSIDE) {
    np_beta_set_status(status, NP_BETA_OK);
    return -INFINITY;
  }

  alpha_minus_one = shape->target_unit * shape->concentration;
  beta_minus_one = shape->target_complement_unit * shape->concentration;
  alpha = 1.0 + alpha_minus_one;
  beta = 1.0 + beta_minus_one;
  if(!R_FINITE(alpha) || !R_FINITE(beta) || alpha < 1.0 || beta < 1.0) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }

  if(shape->observation_endpoint < 0) {
    if(alpha_minus_one == 0.0) {
      log_value = log(beta) - log(shape->support_length);
    } else {
      np_beta_set_status(status, NP_BETA_OK);
      return -INFINITY;
    }
  } else if(shape->observation_endpoint > 0) {
    if(beta_minus_one == 0.0) {
      log_value = log(alpha) - log(shape->support_length);
    } else {
      np_beta_set_status(status, NP_BETA_OK);
      return -INFINITY;
    }
  } else {
    log_value = -log(shape->support_length) +
      alpha_minus_one * shape->log_observation_unit +
      beta_minus_one * shape->log_observation_complement_unit -
      lbeta(alpha, beta);
  }

  if(ISNAN(log_value)) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }
  if(log_value == INFINITY) {
    np_beta_set_status(status, NP_BETA_ERR_RANGE);
    return INFINITY;
  }

  np_beta_set_status(status, NP_BETA_OK);
  return log_value;
}

double np_beta_pdf_scale(const np_beta_shape *shape,
                         np_beta_status *status)
{
  np_beta_status log_status = NP_BETA_OK;
  const double log_value = np_beta_log_pdf_scale(shape, &log_status);

  if(log_status != NP_BETA_OK) {
    np_beta_set_status(status, log_status);
    return (log_status == NP_BETA_ERR_RANGE) ? INFINITY : NAN;
  }
  if(log_value == -INFINITY) {
    np_beta_set_status(status, NP_BETA_OK);
    return 0.0;
  }
  if(log_value > log(DBL_MAX)) {
    np_beta_set_status(status, NP_BETA_ERR_RANGE);
    return INFINITY;
  }

  np_beta_set_status(status, NP_BETA_OK);
  return exp(log_value);
}

double np_beta_pdf_order2(double evaluation,
                          double observation,
                          double bandwidth,
                          double lower,
                          double upper,
                          np_beta_status *status)
{
  np_beta_shape shape;
  np_beta_status shape_status = np_beta_shape_init(evaluation,
                                                   observation,
                                                   bandwidth,
                                                   lower,
                                                   upper,
                                                   1,
                                                   &shape);
  if(shape_status != NP_BETA_OK) {
    np_beta_set_status(status, shape_status);
    return NAN;
  }

  return np_beta_pdf_scale(&shape, status);
}

double np_beta_pdf_order(double evaluation,
                         double observation,
                         double bandwidth,
                         double lower,
                         double upper,
                         int order,
                         np_beta_status *status)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  double positive_log = -INFINITY;
  double negative_log = -INFINITY;
  int component;

  if(component_count == 0) {
    np_beta_set_status(status, NP_BETA_ERR_SCALE);
    return NAN;
  }
  if(order == 2)
    return np_beta_pdf_order2(evaluation, observation, bandwidth,
                              lower, upper, status);

  for(component = 0; component < component_count; ++component) {
    np_beta_shape shape;
    np_beta_status component_status = np_beta_shape_init(
      evaluation, observation, bandwidth, lower, upper,
      component + 1, &shape);
    double log_value;
    double log_term;

    if(component_status != NP_BETA_OK) {
      np_beta_set_status(status, component_status);
      return NAN;
    }
    log_value = np_beta_log_pdf_scale(&shape, &component_status);
    if(component_status != NP_BETA_OK) {
      np_beta_set_status(status, component_status);
      return (component_status == NP_BETA_ERR_RANGE) ? INFINITY : NAN;
    }
    log_term = (log_value == -INFINITY) ? -INFINITY :
      log_value + log((double)abs(coefficients[component]));
    if(coefficients[component] > 0)
      positive_log = np_beta_log_add(positive_log, log_term);
    else
      negative_log = np_beta_log_add(negative_log, log_term);
  }

  return np_beta_signed_log_difference(positive_log, negative_log, status);
}

double np_beta_log_abs_pdf_order(double evaluation,
                                 double observation,
                                 double bandwidth,
                                 double lower,
                                 double upper,
                                 int order,
                                 int *sign,
                                 np_beta_status *status)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  double positive_log = -INFINITY;
  double negative_log = -INFINITY;
  double log_absolute = -INFINITY;
  np_beta_status difference_status;
  int component;

  if(sign == NULL || component_count == 0) {
    np_beta_set_status(status, NP_BETA_ERR_SCALE);
    return NAN;
  }
  *sign = 0;

  for(component = 0; component < component_count; ++component) {
    np_beta_shape shape;
    np_beta_status component_status = np_beta_shape_init(
      evaluation, observation, bandwidth, lower, upper,
      component + 1, &shape);
    double log_value;
    double log_term;

    if(component_status != NP_BETA_OK) {
      np_beta_set_status(status, component_status);
      return NAN;
    }
    log_value = np_beta_log_pdf_scale(&shape, &component_status);
    if(component_status != NP_BETA_OK) {
      np_beta_set_status(status, component_status);
      return (component_status == NP_BETA_ERR_RANGE) ? INFINITY : NAN;
    }
    log_term = (log_value == -INFINITY) ? -INFINITY :
      log_value + log((double)abs(coefficients[component]));
    if(coefficients[component] > 0)
      positive_log = np_beta_log_add(positive_log, log_term);
    else
      negative_log = np_beta_log_add(negative_log, log_term);
  }

  difference_status = np_beta_signed_log_absolute(
    positive_log, negative_log, &log_absolute, sign);
  np_beta_set_status(status, difference_status);
  return (difference_status == NP_BETA_OK) ? log_absolute : NAN;
}

double np_beta_log_pdf_order2(double evaluation,
                              double observation,
                              double bandwidth,
                              double lower,
                              double upper,
                              np_beta_status *status)
{
  np_beta_shape shape;
  np_beta_status shape_status = np_beta_shape_init(evaluation,
                                                   observation,
                                                   bandwidth,
                                                   lower,
                                                   upper,
                                                   1,
                                                   &shape);
  if(shape_status != NP_BETA_OK) {
    np_beta_set_status(status, shape_status);
    return NAN;
  }

  return np_beta_log_pdf_scale(&shape, status);
}

static double np_beta_cdf_scale(double evaluation,
                                double observation,
                                double bandwidth,
                                double lower,
                                double upper,
                                int scale,
                                np_beta_status *status)
{
  np_beta_shape shape;
  np_beta_status shape_status;
  double alpha;
  double beta;
  double tail_coordinate;
  double midpoint;
  double value;

  if(!R_FINITE(evaluation)) {
    np_beta_set_status(status, NP_BETA_ERR_NONFINITE);
    return NAN;
  }

  /* The CDF is observation-centred: the observation supplies the beta
   * shapes, while the distribution target is the incomplete-beta argument.
   * Initializing with the observation in both scalar positions validates the
   * support/bandwidth contract and prepares exactly those shapes. */
  shape_status = np_beta_shape_init(observation,
                                    observation,
                                    bandwidth,
                                    lower,
                                    upper,
                                    scale,
                                    &shape);
  if(shape_status != NP_BETA_OK) {
    np_beta_set_status(status, shape_status);
    return NAN;
  }

  if(evaluation <= lower) {
    np_beta_set_status(status, NP_BETA_OK);
    return 0.0;
  }
  if(evaluation >= upper) {
    np_beta_set_status(status, NP_BETA_OK);
    return 1.0;
  }

  alpha = 1.0 + shape.target_unit * shape.concentration;
  beta = 1.0 + shape.target_complement_unit * shape.concentration;
  midpoint = lower + 0.5 * shape.support_length;
  tail_coordinate = (evaluation <= midpoint) ?
    (evaluation - lower) / shape.support_length :
    (upper - evaluation) / shape.support_length;
  if(!R_FINITE(tail_coordinate) || tail_coordinate <= 0.0 ||
     tail_coordinate >= 1.0 || !R_FINITE(alpha) || !R_FINITE(beta) ||
     alpha < 1.0 || beta < 1.0) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }

  /* Evaluate the nearer tail directly.  Forming 1-(upper-y)/L near the upper
   * endpoint can lose several bits before a concentrated beta CDF amplifies
   * that coordinate error.  The symmetry identity
   * I_s(a,b) = P(Beta(b,a) > 1-s) lets Rmath evaluate the upper branch without
   * constructing the cancellation-prone unit coordinate. */
  value = (evaluation <= midpoint) ?
    pbeta(tail_coordinate, alpha, beta, 1, 0) :
    pbeta(tail_coordinate, beta, alpha, 0, 0);
  if(!R_FINITE(value) || value < 0.0 || value > 1.0) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }

  np_beta_set_status(status, NP_BETA_OK);
  return value;
}

double np_beta_cdf_order2(double evaluation,
                          double observation,
                          double bandwidth,
                          double lower,
                          double upper,
                          np_beta_status *status)
{
  return np_beta_cdf_scale(evaluation, observation, bandwidth,
                           lower, upper, 1, status);
}

double np_beta_cdf_order(double evaluation,
                         double observation,
                         double bandwidth,
                         double lower,
                         double upper,
                         int order,
                         np_beta_status *status)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  double positive_log = -INFINITY;
  double negative_log = -INFINITY;
  int component;

  if(component_count == 0) {
    np_beta_set_status(status, NP_BETA_ERR_SCALE);
    return NAN;
  }
  if(order == 2)
    return np_beta_cdf_order2(evaluation, observation, bandwidth,
                              lower, upper, status);
  if(evaluation <= lower || evaluation >= upper)
    return np_beta_cdf_scale(evaluation, observation, bandwidth,
                             lower, upper, 1, status);

  for(component = 0; component < component_count; ++component) {
    np_beta_status component_status = NP_BETA_OK;
    const double component_value = np_beta_cdf_scale(
      evaluation, observation, bandwidth, lower, upper,
      component + 1, &component_status);
    double log_term;

    if(component_status != NP_BETA_OK) {
      np_beta_set_status(status, component_status);
      return NAN;
    }
    log_term = (component_value == 0.0) ? -INFINITY :
      log(component_value) + log((double)abs(coefficients[component]));
    if(coefficients[component] > 0)
      positive_log = np_beta_log_add(positive_log, log_term);
    else
      negative_log = np_beta_log_add(negative_log, log_term);
  }

  return np_beta_signed_log_difference(positive_log, negative_log, status);
}

static double np_beta_log_overlap_scale(double center_one,
                                        double bandwidth_one,
                                        int scale_one,
                                        double center_two,
                                        double bandwidth_two,
                                        int scale_two,
                                        double lower,
                                        double upper,
                                        np_beta_status *status)
{
  np_beta_shape shape_one;
  np_beta_shape shape_two;
  np_beta_status shape_status;
  double alpha_one;
  double beta_one;
  double alpha_two;
  double beta_two;
  double overlap_alpha;
  double overlap_beta;
  double alpha_mean;
  double beta_mean;
  double total_mean;
  double alpha_displacement;
  double beta_displacement;
  double total_displacement;
  double alpha_centered;
  double beta_centered;
  double total_centered;
  double alpha_half;
  double beta_half;
  double total_half;
  double log_value;

  /* Each center defines one associated beta density in the common
   * observation coordinate.  Passing the center in both scalar positions
   * validates it against the support while preparing its shape parameters. */
  shape_status = np_beta_shape_init(center_one,
                                    center_one,
                                    bandwidth_one,
                                    lower,
                                    upper,
                                    scale_one,
                                    &shape_one);
  if(shape_status != NP_BETA_OK) {
    np_beta_set_status(status, shape_status);
    return NAN;
  }
  shape_status = np_beta_shape_init(center_two,
                                    center_two,
                                    bandwidth_two,
                                    lower,
                                    upper,
                                    scale_two,
                                    &shape_two);
  if(shape_status != NP_BETA_OK) {
    np_beta_set_status(status, shape_status);
    return NAN;
  }

  alpha_one = 1.0 + shape_one.target_unit * shape_one.concentration;
  beta_one = 1.0 +
    shape_one.target_complement_unit * shape_one.concentration;
  alpha_two = 1.0 + shape_two.target_unit * shape_two.concentration;
  beta_two = 1.0 +
    shape_two.target_complement_unit * shape_two.concentration;
  overlap_alpha = alpha_one + alpha_two - 1.0;
  overlap_beta = beta_one + beta_two - 1.0;
  if(!R_FINITE(alpha_one) || !R_FINITE(beta_one) ||
     !R_FINITE(alpha_two) || !R_FINITE(beta_two) ||
     !R_FINITE(overlap_alpha) || !R_FINITE(overlap_beta) ||
     alpha_one < 1.0 || beta_one < 1.0 ||
     alpha_two < 1.0 || beta_two < 1.0 ||
     overlap_alpha < 1.0 || overlap_beta < 1.0) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }

  /* Re-express the beta-function ratio with the Legendre duplication formula.
   * The resulting half-step gamma ratios and centered gamma differences can
   * be evaluated without subtracting O(tau log tau) quantities.  The identity
   * is exact and also covers unequal bandwidths. */
  alpha_mean = 0.5 * (alpha_one + alpha_two);
  beta_mean = 0.5 * (beta_one + beta_two);
  total_mean = alpha_mean + beta_mean;
  alpha_displacement = 0.5 * (alpha_one - alpha_two);
  beta_displacement = 0.5 * (beta_one - beta_two);
  total_displacement = 0.5 *
    ((alpha_one + beta_one) - (alpha_two + beta_two));
  alpha_centered = np_beta_log_gamma_centered(alpha_mean,
                                               alpha_displacement);
  beta_centered = np_beta_log_gamma_centered(beta_mean,
                                              beta_displacement);
  total_centered = np_beta_log_gamma_centered(total_mean,
                                               total_displacement);
  alpha_half = np_beta_log_gamma_backward_half(alpha_mean);
  beta_half = np_beta_log_gamma_backward_half(beta_mean);
  total_half = np_beta_log_gamma_backward_half(total_mean);
  if(!R_FINITE(alpha_centered) || !R_FINITE(beta_centered) ||
     !R_FINITE(total_centered) || !R_FINITE(alpha_half) ||
     !R_FINITE(beta_half) || !R_FINITE(total_half) ||
     total_mean <= 1.0) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }
  log_value = -log(shape_one.support_length) - log(2.0) -
    0.57236494292470008707 +
    alpha_half - alpha_centered + beta_half - beta_centered +
    total_centered + log(total_mean - 1.0) - total_half;
  if(ISNAN(log_value)) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }
  if(log_value == INFINITY) {
    np_beta_set_status(status, NP_BETA_ERR_RANGE);
    return INFINITY;
  }

  np_beta_set_status(status, NP_BETA_OK);
  return log_value;
}

double np_beta_log_overlap_order2(double center_one,
                                  double bandwidth_one,
                                  double center_two,
                                  double bandwidth_two,
                                  double lower,
                                  double upper,
                                  np_beta_status *status)
{
  return np_beta_log_overlap_scale(center_one, bandwidth_one, 1,
                                   center_two, bandwidth_two, 1,
                                   lower, upper, status);
}

double np_beta_overlap_order2(double center_one,
                              double bandwidth_one,
                              double center_two,
                              double bandwidth_two,
                              double lower,
                              double upper,
                              np_beta_status *status)
{
  np_beta_status log_status = NP_BETA_OK;
  const double log_value = np_beta_log_overlap_order2(center_one,
                                                       bandwidth_one,
                                                       center_two,
                                                       bandwidth_two,
                                                       lower,
                                                       upper,
                                                       &log_status);

  if(log_status != NP_BETA_OK) {
    np_beta_set_status(status, log_status);
    return (log_status == NP_BETA_ERR_RANGE) ? INFINITY : NAN;
  }
  if(log_value == -INFINITY) {
    np_beta_set_status(status, NP_BETA_OK);
    return 0.0;
  }
  if(log_value > log(DBL_MAX)) {
    np_beta_set_status(status, NP_BETA_ERR_RANGE);
    return INFINITY;
  }

  np_beta_set_status(status, NP_BETA_OK);
  return exp(log_value);
}

double np_beta_overlap_order(double center_one,
                             double bandwidth_one,
                             double center_two,
                             double bandwidth_two,
                             double lower,
                             double upper,
                             int order,
                             np_beta_status *status)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  double positive_log = -INFINITY;
  double negative_log = -INFINITY;
  int component_one;
  int component_two;

  if(component_count == 0) {
    np_beta_set_status(status, NP_BETA_ERR_SCALE);
    return NAN;
  }
  if(order == 2)
    return np_beta_overlap_order2(center_one, bandwidth_one,
                                  center_two, bandwidth_two,
                                  lower, upper, status);

  for(component_one = 0; component_one < component_count; ++component_one) {
    for(component_two = 0; component_two < component_count; ++component_two) {
      const int coefficient = coefficients[component_one] *
        coefficients[component_two];
      np_beta_status component_status = NP_BETA_OK;
      const double log_value = np_beta_log_overlap_scale(
        center_one, bandwidth_one, component_one + 1,
        center_two, bandwidth_two, component_two + 1,
        lower, upper, &component_status);
      const double log_term = (log_value == -INFINITY) ? -INFINITY :
        log_value + log((double)abs(coefficient));

      if(component_status != NP_BETA_OK) {
        np_beta_set_status(status, component_status);
        return (component_status == NP_BETA_ERR_RANGE) ? INFINITY : NAN;
      }
      if(coefficient > 0)
        positive_log = np_beta_log_add(positive_log, log_term);
      else
        negative_log = np_beta_log_add(negative_log, log_term);
    }
  }

  return np_beta_signed_log_difference(positive_log, negative_log, status);
}
