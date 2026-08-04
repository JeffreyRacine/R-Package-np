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

#if defined(__clang__) || defined(__GNUC__)
# define NP_BETA_KERNEL_ALWAYS_INLINE \
  static inline __attribute__((always_inline))
#else
# define NP_BETA_KERNEL_ALWAYS_INLINE static inline
#endif

NP_BETA_KERNEL_ALWAYS_INLINE np_beta_status
np_beta_shape_target_coordinate_init(
  double evaluation,
  double lower,
  double upper,
  double support_length,
  np_beta_shape *shape)
{
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
    shape->target_unit =
      (evaluation - lower) / support_length;
    shape->target_complement_unit =
      (upper - evaluation) / support_length;
    if(!R_FINITE(shape->target_unit) ||
       shape->target_unit < 0.0 || shape->target_unit > 1.0)
      return NP_BETA_ERR_NUMERIC;
    if(!R_FINITE(shape->target_complement_unit) ||
       shape->target_complement_unit < 0.0 ||
       shape->target_complement_unit > 1.0)
      return NP_BETA_ERR_NUMERIC;
  }
  return NP_BETA_OK;
}

NP_BETA_KERNEL_ALWAYS_INLINE np_beta_status
np_beta_shape_concentration_init(
  double bandwidth,
  double support_length,
  int scale,
  np_beta_shape *shape)
{
  double concentration_limit;
  double ratio;

  ratio = support_length / bandwidth;
  concentration_limit = sqrt(DBL_MAX) * sqrt((double)scale);
  if(!R_FINITE(ratio) || !R_FINITE(concentration_limit) ||
     ratio > concentration_limit)
    return NP_BETA_ERR_CONCENTRATION;

  shape->concentration = (ratio / (double)scale) * ratio;
  if(!R_FINITE(shape->concentration))
    return NP_BETA_ERR_CONCENTRATION;
  return NP_BETA_OK;
}

/* Complete the evaluation/bandwidth-owned part of a shape whose support has
 * already been validated.  Scalar and prepared PDF consumers share this
 * exact arithmetic; keeping it inline avoids adding a call boundary to the
 * public scalar primitive. */
NP_BETA_KERNEL_ALWAYS_INLINE np_beta_status np_beta_shape_target_init(
  double evaluation,
  double bandwidth,
  double lower,
  double upper,
  double support_length,
  int scale,
  np_beta_shape *shape)
{
  const np_beta_status target_status =
    np_beta_shape_target_coordinate_init(
      evaluation, lower, upper, support_length, shape);

  if(target_status != NP_BETA_OK)
    return target_status;
  return np_beta_shape_concentration_init(
    bandwidth, support_length, scale, shape);
}

#undef NP_BETA_KERNEL_ALWAYS_INLINE

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
  shape->observation_unit = np_beta_observation_unit_coordinate(
    observation, lower, upper, shape->support_length);
  shape->observation_complement_unit =
    np_beta_observation_complement_unit_coordinate(
      observation, upper, shape->support_length);
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
    shape->log_observation_unit = np_beta_observation_log_unit(
      observation, lower, shape->support_length);
    shape->log_observation_complement_unit =
      np_beta_observation_log_complement_unit(
        observation, upper, shape->support_length);
  }

  return np_beta_shape_target_init(
    evaluation, bandwidth, lower, upper, shape->support_length,
    scale, shape);
}

void np_beta_pdf_observation_from_shape(
  const np_beta_shape *shape,
  np_beta_pdf_observation *observation)
{
  if(shape == NULL || observation == NULL)
    return;
  observation->log_unit = shape->log_observation_unit;
  observation->log_complement_unit =
    shape->log_observation_complement_unit;
  observation->endpoint = shape->observation_endpoint;
}

np_beta_status np_beta_pdf_component_from_shape(
  const np_beta_shape *shape,
  np_beta_pdf_component *component)
{
  if(shape == NULL || component == NULL)
    return NP_BETA_ERR_NUMERIC;

  component->log_support_length = log(shape->support_length);
  component->eval_location = shape->eval_location;
  component->normalizer_ready = 0;
  component->log_beta = 0.0;
  component->log_abs_coefficient = 0.0;
  component->coefficient_sign = 1;
  component->alpha_minus_one = 0.0;
  component->beta_minus_one = 0.0;
  component->alpha = 1.0;
  component->beta = 1.0;
  if(shape->eval_location != NP_BETA_EVAL_INSIDE)
    return NP_BETA_OK;

  component->alpha_minus_one =
    shape->target_unit * shape->concentration;
  component->beta_minus_one =
    shape->target_complement_unit * shape->concentration;
  component->alpha = 1.0 + component->alpha_minus_one;
  component->beta = 1.0 + component->beta_minus_one;
  if(!R_FINITE(component->alpha) || !R_FINITE(component->beta) ||
     component->alpha < 1.0 || component->beta < 1.0)
    return NP_BETA_ERR_NUMERIC;

  return NP_BETA_OK;
}

np_beta_status np_beta_pdf_component_prepare_coefficient(
  np_beta_pdf_component *component,
  int order,
  int component_index)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  int coefficient;

  if(component == NULL || component_index < 0 ||
     component_index >= component_count)
    return NP_BETA_ERR_SCALE;
  coefficient = coefficients[component_index];
  if(coefficient == 0)
    return NP_BETA_ERR_SCALE;
  component->log_abs_coefficient =
    log((double)abs(coefficient));
  component->coefficient_sign = coefficient > 0 ? 1 : -1;
  return R_FINITE(component->log_abs_coefficient) ?
    NP_BETA_OK : NP_BETA_ERR_NUMERIC;
}

np_beta_status np_beta_pdf_component_prepare_normalizer(
  np_beta_pdf_component *component)
{
  if(component == NULL)
    return NP_BETA_ERR_NUMERIC;
  /* An outside-support target is an exact zero PDF row.  No observation can
   * consume a normalizer, so retain the unprepared state and succeed. */
  if(component->eval_location != NP_BETA_EVAL_INSIDE)
    return NP_BETA_OK;
  component->log_beta = lbeta(component->alpha, component->beta);
  if(ISNAN(component->log_beta))
    return NP_BETA_ERR_NUMERIC;
  if(component->log_beta == INFINITY)
    return NP_BETA_ERR_RANGE;
  component->normalizer_ready = 1;
  return NP_BETA_OK;
}

double np_beta_log_pdf_component_prepared(
  const np_beta_pdf_component *component,
  const np_beta_pdf_observation *observation,
  np_beta_status *status)
{
  double log_value;

  if(component == NULL || observation == NULL) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }
  if(component->eval_location != NP_BETA_EVAL_INSIDE) {
    np_beta_set_status(status, NP_BETA_OK);
    return -INFINITY;
  }

  if(observation->endpoint < 0) {
    if(component->alpha_minus_one == 0.0) {
      log_value = log(component->beta) - component->log_support_length;
    } else {
      np_beta_set_status(status, NP_BETA_OK);
      return -INFINITY;
    }
  } else if(observation->endpoint > 0) {
    if(component->beta_minus_one == 0.0) {
      log_value = log(component->alpha) - component->log_support_length;
    } else {
      np_beta_set_status(status, NP_BETA_OK);
      return -INFINITY;
    }
  } else {
    if(!component->normalizer_ready) {
      np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
      return NAN;
    }
    log_value = -component->log_support_length +
      component->alpha_minus_one * observation->log_unit +
      component->beta_minus_one * observation->log_complement_unit -
      component->log_beta;
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

double np_beta_log_abs_pdf_prepared(
  const np_beta_pdf_component *components,
  const np_beta_pdf_observation *observation,
  int component_count,
  int *sign,
  np_beta_status *status)
{
  double positive_log = -INFINITY;
  double negative_log = -INFINITY;
  double log_absolute = -INFINITY;
  np_beta_status difference_status;
  int component;

  if(components == NULL || observation == NULL || sign == NULL ||
     component_count <= 0 ||
     component_count > NP_BETA_ORDER_MAX_COMPONENTS) {
    np_beta_set_status(status, NP_BETA_ERR_SCALE);
    return NAN;
  }
  *sign = 0;
  for(component = 0; component < component_count; ++component) {
    np_beta_status component_status = NP_BETA_OK;
    const double log_value = np_beta_log_pdf_component_prepared(
      &components[component], observation, &component_status);
    double log_term;

    if(component_status != NP_BETA_OK) {
      np_beta_set_status(status, component_status);
      return (component_status == NP_BETA_ERR_RANGE) ? INFINITY : NAN;
    }
    if(components[component].coefficient_sign != -1 &&
       components[component].coefficient_sign != 1) {
      np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
      return NAN;
    }
    log_term = (log_value == -INFINITY) ? -INFINITY :
      log_value + components[component].log_abs_coefficient;
    if(components[component].coefficient_sign > 0)
      positive_log = np_beta_log_add(positive_log, log_term);
    else
      negative_log = np_beta_log_add(negative_log, log_term);
  }

  difference_status = np_beta_signed_log_absolute(
    positive_log, negative_log, &log_absolute, sign);
  np_beta_set_status(status, difference_status);
  return (difference_status == NP_BETA_OK) ? log_absolute : NAN;
}

/* Adaptive-NN PDF rows own their bandwidth by observation.  The observation
 * coordinate and coefficient algebra are nevertheless invocation-invariant.
 * Prepare only the bandwidth-dependent target components here, preserving
 * the scalar shape-init arithmetic and its error precedence. */
double np_beta_log_abs_pdf_order_prepared_observation(
  double evaluation,
  double bandwidth,
  double lower,
  double upper,
  int order,
  const np_beta_pdf_observation *observation,
  np_beta_status observation_status,
  const double *log_abs_coefficient,
  const signed char *coefficient_sign,
  int *sign,
  np_beta_status *status)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  np_beta_pdf_component components[NP_BETA_ORDER_MAX_COMPONENTS];
  np_beta_shape shape_template;
  np_beta_status component_status;
  int component;

  (void)coefficients;
  if(sign == NULL || observation == NULL ||
     log_abs_coefficient == NULL || coefficient_sign == NULL ||
     component_count == 0) {
    np_beta_set_status(status, NP_BETA_ERR_SCALE);
    return NAN;
  }
  *sign = 0;

  if(!R_FINITE(evaluation) || !R_FINITE(lower) || !R_FINITE(upper) ||
     observation_status == NP_BETA_ERR_NONFINITE) {
    np_beta_set_status(status, NP_BETA_ERR_NONFINITE);
    return NAN;
  }
  if(!R_FINITE(bandwidth) || bandwidth <= 0.0) {
    np_beta_set_status(status, NP_BETA_ERR_BANDWIDTH);
    return NAN;
  }
  shape_template.lower = lower;
  shape_template.upper = upper;
  shape_template.support_length = upper - lower;
  if(!R_FINITE(shape_template.support_length) ||
     shape_template.support_length <= 0.0) {
    np_beta_set_status(status, NP_BETA_ERR_BOUNDS);
    return NAN;
  }
  if(observation_status != NP_BETA_OK) {
    np_beta_set_status(status, observation_status);
    return NAN;
  }

  shape_template.log_observation_unit = observation->log_unit;
  shape_template.log_observation_complement_unit =
    observation->log_complement_unit;
  shape_template.observation_endpoint = observation->endpoint;
  component_status = np_beta_shape_target_coordinate_init(
    evaluation, lower, upper, shape_template.support_length,
    &shape_template);
  if(component_status != NP_BETA_OK) {
    np_beta_set_status(status, component_status);
    return NAN;
  }

  for(component = 0; component < component_count; ++component) {
    np_beta_shape shape = shape_template;

    component_status = np_beta_shape_concentration_init(
      bandwidth, shape.support_length, component + 1, &shape);
    if(component_status != NP_BETA_OK) {
      np_beta_set_status(status, component_status);
      return NAN;
    }

    component_status = np_beta_pdf_component_from_shape(
      &shape, &components[component]);
    if(component_status == NP_BETA_OK && observation->endpoint == 0)
      component_status = np_beta_pdf_component_prepare_normalizer(
        &components[component]);
    if(component_status != NP_BETA_OK) {
      np_beta_set_status(status, component_status);
      return (component_status == NP_BETA_ERR_RANGE) ? INFINITY : NAN;
    }
    components[component].log_abs_coefficient =
      log_abs_coefficient[component];
    components[component].coefficient_sign = coefficient_sign[component];
  }

  return np_beta_log_abs_pdf_prepared(
    components, observation, component_count, sign, status);
}

double np_beta_log_pdf_scale(const np_beta_shape *shape,
                             np_beta_status *status)
{
  np_beta_pdf_observation observation;
  np_beta_pdf_component component;
  np_beta_status component_status;

  if(shape == NULL) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }

  np_beta_pdf_observation_from_shape(shape, &observation);
  component_status = np_beta_pdf_component_from_shape(shape, &component);
  if(component_status == NP_BETA_OK && observation.endpoint == 0)
    component_status = np_beta_pdf_component_prepare_normalizer(&component);
  if(component_status != NP_BETA_OK) {
    np_beta_set_status(status, component_status);
    return (component_status == NP_BETA_ERR_RANGE) ? INFINITY : NAN;
  }
  return np_beta_log_pdf_component_prepared(
    &component, &observation, status);
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

static void np_beta_derivative_zero(np_beta_derivative *derivative)
{
  derivative->regular_log_absolute = -INFINITY;
  derivative->jump_log_absolute = -INFINITY;
  derivative->regular_sign = 0;
  derivative->jump_sign = 0;
}

static np_beta_status np_beta_derivative_accumulate(double log_term,
                                                    int sign,
                                                    double *positive_log,
                                                    double *negative_log)
{
  if(ISNAN(log_term) || (sign != -1 && sign != 1))
    return NP_BETA_ERR_NUMERIC;
  if(log_term == -INFINITY)
    return NP_BETA_OK;
  if(log_term == INFINITY)
    return NP_BETA_ERR_RANGE;
  if(sign > 0)
    *positive_log = np_beta_log_add(*positive_log, log_term);
  else
    *negative_log = np_beta_log_add(*negative_log, log_term);
  return NP_BETA_OK;
}

static np_beta_status np_beta_derivative_component_accumulate(
  const np_beta_shape *shape,
  double evaluation,
  double observation,
  double log_pdf,
  double log_abs_coefficient,
  int coefficient_sign,
  double *regular_positive_log,
  double *regular_negative_log,
  double *jump_positive_log,
  double *jump_negative_log)
{
  double score;
  double log_term;
  int term_sign;

  if(shape == NULL || coefficient_sign == 0 ||
     !R_FINITE(log_abs_coefficient))
    return NP_BETA_ERR_NUMERIC;

  /* Exact matching endpoints are jump discontinuities.  Their regular
   * one-sided derivative is zero; the signed jump is retained so ratio
   * estimators can cancel it structurally. */
  if((evaluation == shape->lower && observation == shape->lower) ||
     (evaluation == shape->upper && observation == shape->upper)) {
    log_term = log_pdf + log_abs_coefficient;
    term_sign = coefficient_sign;
    if(evaluation == shape->lower)
      term_sign = -term_sign;
    return np_beta_derivative_accumulate(
      log_term, term_sign, jump_positive_log, jump_negative_log);
  }

  if(log_pdf == -INFINITY)
    return NP_BETA_OK;

  score = (shape->concentration / shape->support_length) *
    (shape->log_observation_unit -
     shape->log_observation_complement_unit -
     digamma(1.0 + shape->target_unit * shape->concentration) +
     digamma(1.0 + shape->target_complement_unit * shape->concentration));
  if(ISNAN(score))
    return NP_BETA_ERR_NUMERIC;
  if(score == 0.0)
    return NP_BETA_OK;
  if(!R_FINITE(score))
    return NP_BETA_ERR_RANGE;

  log_term = log_pdf + log(fabs(score)) + log_abs_coefficient;
  term_sign = coefficient_sign * ((score > 0.0) ? 1 : -1);
  return np_beta_derivative_accumulate(
    log_term, term_sign, regular_positive_log, regular_negative_log);
}

np_beta_status np_beta_pdf_derivative_order(double evaluation,
                                            double observation,
                                            double bandwidth,
                                            double lower,
                                            double upper,
                                            int order,
                                            np_beta_derivative *derivative)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  double regular_positive_log = -INFINITY;
  double regular_negative_log = -INFINITY;
  double jump_positive_log = -INFINITY;
  double jump_negative_log = -INFINITY;
  int component;

  if(derivative == NULL)
    return NP_BETA_ERR_NUMERIC;
  np_beta_derivative_zero(derivative);
  if(component_count == 0)
    return NP_BETA_ERR_SCALE;

  for(component = 0; component < component_count; ++component) {
    np_beta_shape shape;
    np_beta_status status = np_beta_shape_init(
      evaluation, observation, bandwidth, lower, upper,
      component + 1, &shape);
    double log_pdf;
    double log_abs_coefficient;
    int coefficient_sign;

    if(status != NP_BETA_OK)
      return status;
    if(shape.eval_location != NP_BETA_EVAL_INSIDE ||
       shape.concentration == 0.0)
      continue;
    log_pdf = np_beta_log_pdf_scale(&shape, &status);
    if(status != NP_BETA_OK)
      return status;
    log_abs_coefficient = log((double)abs(coefficients[component]));
    coefficient_sign = coefficients[component] > 0 ? 1 : -1;
    status = np_beta_derivative_component_accumulate(
      &shape, evaluation, observation, log_pdf,
      log_abs_coefficient, coefficient_sign,
      &regular_positive_log, &regular_negative_log,
      &jump_positive_log, &jump_negative_log);
    if(status != NP_BETA_OK)
      return status;
  }

  {
    np_beta_status status = np_beta_signed_log_absolute(
      regular_positive_log, regular_negative_log,
      &derivative->regular_log_absolute, &derivative->regular_sign);
    if(status != NP_BETA_OK)
      return status;
    status = np_beta_signed_log_absolute(
      jump_positive_log, jump_negative_log,
      &derivative->jump_log_absolute, &derivative->jump_sign);
    return status;
  }
}

np_beta_status np_beta_log_abs_pdf_derivative_order(
  double evaluation,
  double observation,
  double bandwidth,
  double lower,
  double upper,
  int order,
  double *level_log_absolute,
  int *level_sign,
  np_beta_derivative *derivative)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  np_beta_shape shape[NP_BETA_ORDER_MAX_COMPONENTS];
  double log_pdf[NP_BETA_ORDER_MAX_COMPONENTS];
  double log_abs_coefficient[NP_BETA_ORDER_MAX_COMPONENTS];
  double level_positive_log = -INFINITY;
  double level_negative_log = -INFINITY;
  double regular_positive_log = -INFINITY;
  double regular_negative_log = -INFINITY;
  double jump_positive_log = -INFINITY;
  double jump_negative_log = -INFINITY;
  np_beta_status status;
  int component;

  if(level_log_absolute == NULL || level_sign == NULL ||
     derivative == NULL || component_count == 0)
    return NP_BETA_ERR_SCALE;
  *level_log_absolute = -INFINITY;
  *level_sign = 0;
  np_beta_derivative_zero(derivative);

  /* Complete the level phase before derivative work.  Besides preserving
   * arithmetic order, this retains the incumbent first-failure contract of
   * np_beta_log_abs_pdf_order() followed by np_beta_pdf_derivative_order(). */
  for(component = 0; component < component_count; ++component) {
    double log_term;

    status = np_beta_shape_init(
      evaluation, observation, bandwidth, lower, upper,
      component + 1, &shape[component]);
    if(status != NP_BETA_OK)
      return status;
    log_pdf[component] = np_beta_log_pdf_scale(
      &shape[component], &status);
    if(status != NP_BETA_OK)
      return status;
    log_abs_coefficient[component] =
      log((double)abs(coefficients[component]));
    log_term = log_pdf[component] == -INFINITY ? -INFINITY :
      log_pdf[component] + log_abs_coefficient[component];
    if(coefficients[component] > 0)
      level_positive_log = np_beta_log_add(level_positive_log, log_term);
    else
      level_negative_log = np_beta_log_add(level_negative_log, log_term);
  }
  status = np_beta_signed_log_absolute(
    level_positive_log, level_negative_log,
    level_log_absolute, level_sign);
  if(status != NP_BETA_OK)
    return status;

  for(component = 0; component < component_count; ++component) {
    if(shape[component].eval_location != NP_BETA_EVAL_INSIDE ||
       shape[component].concentration == 0.0)
      continue;
    status = np_beta_derivative_component_accumulate(
      &shape[component], evaluation, observation, log_pdf[component],
      log_abs_coefficient[component],
      coefficients[component] > 0 ? 1 : -1,
      &regular_positive_log, &regular_negative_log,
      &jump_positive_log, &jump_negative_log);
    if(status != NP_BETA_OK)
      return status;
  }

  status = np_beta_signed_log_absolute(
    regular_positive_log, regular_negative_log,
    &derivative->regular_log_absolute, &derivative->regular_sign);
  if(status != NP_BETA_OK)
    return status;
  return np_beta_signed_log_absolute(
    jump_positive_log, jump_negative_log,
    &derivative->jump_log_absolute, &derivative->jump_sign);
}

np_beta_status np_beta_log_abs_pdf_derivative_prepared(
  const np_beta_pdf_component *components,
  const np_beta_pdf_derivative_component *derivative_components,
  const np_beta_pdf_observation *observation,
  int component_count,
  double evaluation,
  double observed,
  double lower,
  double upper,
  double *level_log_absolute,
  int *level_sign,
  np_beta_derivative *derivative)
{
  double regular_positive_log = -INFINITY;
  double regular_negative_log = -INFINITY;
  double jump_positive_log = -INFINITY;
  double jump_negative_log = -INFINITY;
  np_beta_status status = NP_BETA_OK;
  int component;

  if(components == NULL || derivative_components == NULL ||
     observation == NULL ||
     level_log_absolute == NULL || level_sign == NULL ||
     derivative == NULL || component_count <= 0 ||
     component_count > NP_BETA_ORDER_MAX_COMPONENTS)
    return NP_BETA_ERR_SCALE;
  *level_log_absolute = -INFINITY;
  *level_sign = 0;
  np_beta_derivative_zero(derivative);

  /* The canonical prepared level owner remains responsible for the complete
   * level phase.  Only after it succeeds do we consume the same component
   * state for the derivative, preserving the fused scalar failure order. */
  *level_log_absolute = np_beta_log_abs_pdf_prepared(
    components, observation, component_count, level_sign, &status);
  if(status != NP_BETA_OK)
    return status;

  for(component = 0; component < component_count; ++component) {
    np_beta_shape shape;
    double log_pdf;

    if(components[component].eval_location != NP_BETA_EVAL_INSIDE ||
       derivative_components[component].concentration == 0.0)
      continue;
    log_pdf = np_beta_log_pdf_component_prepared(
      &components[component], observation, &status);
    if(status != NP_BETA_OK)
      return status;
    shape.lower = lower;
    shape.upper = upper;
    shape.support_length =
      derivative_components[component].support_length;
    shape.concentration =
      derivative_components[component].concentration;
    shape.target_unit = derivative_components[component].target_unit;
    shape.target_complement_unit =
      derivative_components[component].target_complement_unit;
    shape.log_observation_unit = observation->log_unit;
    shape.log_observation_complement_unit =
      observation->log_complement_unit;
    status = np_beta_derivative_component_accumulate(
      &shape, evaluation, observed, log_pdf,
      components[component].log_abs_coefficient,
      components[component].coefficient_sign,
      &regular_positive_log, &regular_negative_log,
      &jump_positive_log, &jump_negative_log);
    if(status != NP_BETA_OK)
      return status;
  }

  status = np_beta_signed_log_absolute(
    regular_positive_log, regular_negative_log,
    &derivative->regular_log_absolute, &derivative->regular_sign);
  if(status != NP_BETA_OK)
    return status;
  return np_beta_signed_log_absolute(
    jump_positive_log, jump_negative_log,
    &derivative->jump_log_absolute, &derivative->jump_sign);
}

np_beta_status
np_beta_log_abs_pdf_derivative_order_prepared_observation(
  double evaluation,
  double observed,
  double bandwidth,
  double lower,
  double upper,
  int order,
  const np_beta_pdf_observation *observation,
  np_beta_status observation_status,
  const double *log_abs_coefficient,
  const signed char *coefficient_sign,
  double *level_log_absolute,
  int *level_sign,
  np_beta_derivative *derivative)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  np_beta_shape shape[NP_BETA_ORDER_MAX_COMPONENTS];
  np_beta_pdf_component components[NP_BETA_ORDER_MAX_COMPONENTS];
  double log_pdf[NP_BETA_ORDER_MAX_COMPONENTS];
  double level_positive_log = -INFINITY;
  double level_negative_log = -INFINITY;
  double regular_positive_log = -INFINITY;
  double regular_negative_log = -INFINITY;
  double jump_positive_log = -INFINITY;
  double jump_negative_log = -INFINITY;
  np_beta_shape shape_template;
  np_beta_status status;
  int component;

  (void)coefficients;
  if(observation == NULL || log_abs_coefficient == NULL ||
     coefficient_sign == NULL || level_log_absolute == NULL ||
     level_sign == NULL || derivative == NULL || component_count == 0)
    return NP_BETA_ERR_SCALE;
  *level_log_absolute = -INFINITY;
  *level_sign = 0;
  np_beta_derivative_zero(derivative);

  if(!R_FINITE(evaluation) || !R_FINITE(lower) || !R_FINITE(upper) ||
     observation_status == NP_BETA_ERR_NONFINITE)
    return NP_BETA_ERR_NONFINITE;
  if(!R_FINITE(bandwidth) || bandwidth <= 0.0)
    return NP_BETA_ERR_BANDWIDTH;
  shape_template.lower = lower;
  shape_template.upper = upper;
  shape_template.support_length = upper - lower;
  if(!R_FINITE(shape_template.support_length) ||
     shape_template.support_length <= 0.0)
    return NP_BETA_ERR_BOUNDS;
  if(observation_status != NP_BETA_OK)
    return observation_status;

  shape_template.log_observation_unit = observation->log_unit;
  shape_template.log_observation_complement_unit =
    observation->log_complement_unit;
  shape_template.observation_endpoint = observation->endpoint;
  status = np_beta_shape_target_coordinate_init(
    evaluation, lower, upper, shape_template.support_length,
    &shape_template);
  if(status != NP_BETA_OK)
    return status;

  /* Complete the level phase before derivative work, matching the scalar
   * fused primitive while consuming invocation-owned observation state. */
  for(component = 0; component < component_count; ++component) {
    double log_term;

    if(!R_FINITE(log_abs_coefficient[component]) ||
       (coefficient_sign[component] != -1 &&
        coefficient_sign[component] != 1))
      return NP_BETA_ERR_NUMERIC;
    shape[component] = shape_template;
    status = np_beta_shape_concentration_init(
      bandwidth, shape[component].support_length,
      component + 1, &shape[component]);
    if(status != NP_BETA_OK)
      return status;
    status = np_beta_pdf_component_from_shape(
      &shape[component], &components[component]);
    if(status == NP_BETA_OK && observation->endpoint == 0)
      status = np_beta_pdf_component_prepare_normalizer(
        &components[component]);
    if(status != NP_BETA_OK)
      return status;
    components[component].log_abs_coefficient =
      log_abs_coefficient[component];
    components[component].coefficient_sign =
      coefficient_sign[component];
    log_pdf[component] = np_beta_log_pdf_component_prepared(
      &components[component], observation, &status);
    if(status != NP_BETA_OK)
      return status;
    log_term = log_pdf[component] == -INFINITY ? -INFINITY :
      log_pdf[component] + log_abs_coefficient[component];
    if(coefficient_sign[component] > 0)
      level_positive_log = np_beta_log_add(level_positive_log, log_term);
    else
      level_negative_log = np_beta_log_add(level_negative_log, log_term);
  }
  status = np_beta_signed_log_absolute(
    level_positive_log, level_negative_log,
    level_log_absolute, level_sign);
  if(status != NP_BETA_OK)
    return status;

  for(component = 0; component < component_count; ++component) {
    if(shape[component].eval_location != NP_BETA_EVAL_INSIDE ||
       shape[component].concentration == 0.0)
      continue;
    status = np_beta_derivative_component_accumulate(
      &shape[component], evaluation, observed, log_pdf[component],
      log_abs_coefficient[component], coefficient_sign[component],
      &regular_positive_log, &regular_negative_log,
      &jump_positive_log, &jump_negative_log);
    if(status != NP_BETA_OK)
      return status;
  }
  status = np_beta_signed_log_absolute(
    regular_positive_log, regular_negative_log,
    &derivative->regular_log_absolute, &derivative->regular_sign);
  if(status != NP_BETA_OK)
    return status;
  return np_beta_signed_log_absolute(
    jump_positive_log, jump_negative_log,
    &derivative->jump_log_absolute, &derivative->jump_sign);
}

static double np_beta_signed_log_value(double log_absolute,
                                       int sign,
                                       np_beta_status *status)
{
  if(sign == 0 || log_absolute == -INFINITY) {
    np_beta_set_status(status, NP_BETA_OK);
    return 0.0;
  }
  if((sign != -1 && sign != 1) || ISNAN(log_absolute)) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }
  if(log_absolute > log(DBL_MAX)) {
    np_beta_set_status(status, NP_BETA_ERR_RANGE);
    return (sign > 0) ? INFINITY : -INFINITY;
  }
  np_beta_set_status(status, NP_BETA_OK);
  return (sign > 0) ? exp(log_absolute) : -exp(log_absolute);
}

double np_beta_derivative_regular_value(
  const np_beta_derivative *derivative,
  np_beta_status *status)
{
  if(derivative == NULL) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }
  return np_beta_signed_log_value(derivative->regular_log_absolute,
                                  derivative->regular_sign, status);
}

double np_beta_derivative_public_value(
  const np_beta_derivative *derivative,
  np_beta_status *status)
{
  if(derivative == NULL) {
    np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
    return NAN;
  }
  if(derivative->jump_sign != 0) {
    np_beta_set_status(status, NP_BETA_OK);
    return (derivative->jump_sign > 0) ? INFINITY : -INFINITY;
  }
  return np_beta_derivative_regular_value(derivative, status);
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

/* Scalar and prepared CDF consumers must share the same incomplete-beta
 * boundary, tail, and numeric contract.  The caller owns shape preparation;
 * this small primitive owns the one canonical CDF calculation. */
#if defined(__clang__) || defined(__GNUC__)
# define NP_BETA_CDF_ALWAYS_INLINE \
  static inline __attribute__((always_inline))
#else
# define NP_BETA_CDF_ALWAYS_INLINE static inline
#endif

NP_BETA_CDF_ALWAYS_INLINE double np_beta_cdf_from_shape(
  double evaluation,
  const np_beta_shape *shape,
  np_beta_status *status)
{
  double alpha;
  double beta;
  double tail_coordinate;
  double midpoint;
  double value;

  if(evaluation <= shape->lower) {
    np_beta_set_status(status, NP_BETA_OK);
    return 0.0;
  }
  if(evaluation >= shape->upper) {
    np_beta_set_status(status, NP_BETA_OK);
    return 1.0;
  }

  alpha = 1.0 + shape->target_unit * shape->concentration;
  beta = 1.0 + shape->target_complement_unit * shape->concentration;
  midpoint = shape->lower + 0.5 * shape->support_length;
  tail_coordinate = (evaluation <= midpoint) ?
    (evaluation - shape->lower) / shape->support_length :
    (shape->upper - evaluation) / shape->support_length;
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

#undef NP_BETA_CDF_ALWAYS_INLINE

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

  return np_beta_cdf_from_shape(evaluation, &shape, status);
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

static np_beta_status np_beta_cdf_order_log_parts(
  double evaluation,
  double observation,
  double bandwidth,
  double lower,
  double upper,
  int order,
  double *positive_log,
  double *negative_log)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  int component;

  if(component_count <= 1 || positive_log == NULL || negative_log == NULL)
    return NP_BETA_ERR_SCALE;
  *positive_log = -INFINITY;
  *negative_log = -INFINITY;

  for(component = 0; component < component_count; ++component) {
    np_beta_status component_status = NP_BETA_OK;
    const double component_value = np_beta_cdf_scale(
      evaluation, observation, bandwidth, lower, upper,
      component + 1, &component_status);
    double log_term;

    if(component_status != NP_BETA_OK)
      return component_status;
    log_term = (component_value == 0.0) ? -INFINITY :
      log(component_value) + log((double)abs(coefficients[component]));
    if(coefficients[component] > 0)
      *positive_log = np_beta_log_add(*positive_log, log_term);
    else
      *negative_log = np_beta_log_add(*negative_log, log_term);
  }

  return NP_BETA_OK;
}

double np_beta_cdf_order(double evaluation,
                         double observation,
                         double bandwidth,
                         double lower,
                         double upper,
                         int order,
                         np_beta_status *status)
{
  double positive_log = -INFINITY;
  double negative_log = -INFINITY;
  np_beta_status parts_status;

  if(!np_beta_order_supported(order)) {
    np_beta_set_status(status, NP_BETA_ERR_SCALE);
    return NAN;
  }
  if(order == 2)
    return np_beta_cdf_order2(evaluation, observation, bandwidth,
                              lower, upper, status);
  if(evaluation <= lower || evaluation >= upper)
    return np_beta_cdf_scale(evaluation, observation, bandwidth,
                             lower, upper, 1, status);

  parts_status = np_beta_cdf_order_log_parts(
    evaluation, observation, bandwidth, lower, upper, order,
    &positive_log, &negative_log);
  if(parts_status != NP_BETA_OK) {
    np_beta_set_status(status, parts_status);
    return NAN;
  }

  return np_beta_signed_log_difference(positive_log, negative_log, status);
}

double np_beta_log_abs_cdf_order(double evaluation,
                                 double observation,
                                 double bandwidth,
                                 double lower,
                                 double upper,
                                 int order,
                                 int *sign,
                                 np_beta_status *status)
{
  double positive_log = -INFINITY;
  double negative_log = -INFINITY;
  double log_absolute = -INFINITY;
  np_beta_status parts_status;

  if(sign == NULL || !np_beta_order_supported(order)) {
    np_beta_set_status(status, NP_BETA_ERR_SCALE);
    return NAN;
  }
  *sign = 0;

  if(order == 2 || evaluation <= lower || evaluation >= upper) {
    np_beta_status value_status = NP_BETA_OK;
    const double value = np_beta_cdf_order(
      evaluation, observation, bandwidth, lower, upper,
      order, &value_status);

    if(value_status != NP_BETA_OK || !R_FINITE(value)) {
      np_beta_set_status(status, value_status);
      return NAN;
    }
    if(value != 0.0) {
      log_absolute = log(fabs(value));
      *sign = (value > 0.0) ? 1 : -1;
    }
    np_beta_set_status(status, NP_BETA_OK);
    return log_absolute;
  }

  parts_status = np_beta_cdf_order_log_parts(
    evaluation, observation, bandwidth, lower, upper, order,
    &positive_log, &negative_log);
  if(parts_status == NP_BETA_OK)
    parts_status = np_beta_signed_log_absolute(
      positive_log, negative_log, &log_absolute, sign);
  np_beta_set_status(status, parts_status);
  return (parts_status == NP_BETA_OK) ? log_absolute : NAN;
}

np_beta_status np_beta_cdf_observation_init(
  double observation,
  double lower,
  double upper,
  double support_length,
  np_beta_cdf_observation *prepared)
{
  np_beta_pdf_observation pdf_observation;
  np_beta_shape shape = {0};
  double observation_unit;
  double observation_complement_unit;
  np_beta_status status;

  if(prepared == NULL)
    return NP_BETA_ERR_NUMERIC;
  status = np_beta_pdf_observation_init(
    observation, lower, upper, support_length,
    &observation_unit, &observation_complement_unit,
    &pdf_observation.log_unit, &pdf_observation.log_complement_unit,
    &pdf_observation.endpoint);
  if(status != NP_BETA_OK)
    return status;
  status = np_beta_shape_target_coordinate_init(
    observation, lower, upper, support_length, &shape);
  if(status != NP_BETA_OK)
    return status;
  prepared->target_unit = shape.target_unit;
  prepared->target_complement_unit = shape.target_complement_unit;
  return NP_BETA_OK;
}

double np_beta_log_abs_cdf_order_prepared_observation(
  double evaluation,
  double bandwidth,
  double lower,
  double upper,
  int order,
  const np_beta_cdf_observation *observation,
  np_beta_status observation_status,
  const double *log_abs_coefficient,
  const signed char *coefficient_sign,
  int *sign,
  np_beta_status *status)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  double positive_log = -INFINITY;
  double negative_log = -INFINITY;
  double support_length;
  int component_limit;
  int component;

  (void)coefficients;
  if(sign == NULL || observation == NULL ||
     log_abs_coefficient == NULL || coefficient_sign == NULL ||
     component_count == 0) {
    np_beta_set_status(status, NP_BETA_ERR_SCALE);
    return NAN;
  }
  *sign = 0;
  if(!R_FINITE(evaluation) || !R_FINITE(lower) || !R_FINITE(upper) ||
     observation_status == NP_BETA_ERR_NONFINITE) {
    np_beta_set_status(status, NP_BETA_ERR_NONFINITE);
    return NAN;
  }
  if(!R_FINITE(bandwidth) || bandwidth <= 0.0) {
    np_beta_set_status(status, NP_BETA_ERR_BANDWIDTH);
    return NAN;
  }
  support_length = upper - lower;
  if(!R_FINITE(support_length) || support_length <= 0.0) {
    np_beta_set_status(status, NP_BETA_ERR_BOUNDS);
    return NAN;
  }
  if(observation_status != NP_BETA_OK) {
    np_beta_set_status(status, observation_status);
    return NAN;
  }

  component_limit = (evaluation <= lower || evaluation >= upper) ?
    1 : component_count;

  for(component = 0; component < component_limit; ++component) {
    np_beta_shape shape = {0};
    double component_value;
    double log_term;
    np_beta_status component_status;

    if(!R_FINITE(log_abs_coefficient[component]) ||
       (coefficient_sign[component] != -1 &&
        coefficient_sign[component] != 1)) {
      np_beta_set_status(status, NP_BETA_ERR_NUMERIC);
      return NAN;
    }
    shape.lower = lower;
    shape.upper = upper;
    shape.support_length = support_length;
    shape.target_unit = observation->target_unit;
    shape.target_complement_unit = observation->target_complement_unit;
    component_status = np_beta_shape_concentration_init(
      bandwidth, support_length, component + 1, &shape);
    if(component_status != NP_BETA_OK) {
      np_beta_set_status(status, component_status);
      return NAN;
    }
    component_value = np_beta_cdf_from_shape(
      evaluation, &shape, &component_status);
    if(component_status != NP_BETA_OK) {
      np_beta_set_status(status, component_status);
      return NAN;
    }
    if(evaluation <= lower) {
      np_beta_set_status(status, NP_BETA_OK);
      return -INFINITY;
    }
    if(evaluation >= upper) {
      *sign = 1;
      np_beta_set_status(status, NP_BETA_OK);
      return 0.0;
    }
    if(component_count == 1) {
      if(component_value != 0.0) {
        *sign = 1;
        positive_log = log(fabs(component_value));
      }
      np_beta_set_status(status, NP_BETA_OK);
      return positive_log;
    }
    log_term = component_value == 0.0 ? -INFINITY :
      log(component_value) + log_abs_coefficient[component];
    if(coefficient_sign[component] > 0)
      positive_log = np_beta_log_add(positive_log, log_term);
    else
      negative_log = np_beta_log_add(negative_log, log_term);
  }

  {
    double log_absolute = -INFINITY;
    const np_beta_status parts_status = np_beta_signed_log_absolute(
      positive_log, negative_log, &log_absolute, sign);

    np_beta_set_status(status, parts_status);
    return parts_status == NP_BETA_OK ? log_absolute : NAN;
  }
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

static np_beta_status np_beta_overlap_order_log_parts(
  double center_one,
  double bandwidth_one,
  double center_two,
  double bandwidth_two,
  double lower,
  double upper,
  int order,
  double *positive_log,
  double *negative_log)
{
  const int *coefficients = NULL;
  const int component_count =
    np_beta_order_coefficients(order, &coefficients);
  int component_one;
  int component_two;

  if(component_count <= 1 || positive_log == NULL || negative_log == NULL)
    return NP_BETA_ERR_SCALE;
  *positive_log = -INFINITY;
  *negative_log = -INFINITY;

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

      if(component_status != NP_BETA_OK)
        return component_status;
      if(coefficient > 0)
        *positive_log = np_beta_log_add(*positive_log, log_term);
      else
        *negative_log = np_beta_log_add(*negative_log, log_term);
    }
  }

  return NP_BETA_OK;
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
  double positive_log = -INFINITY;
  double negative_log = -INFINITY;
  np_beta_status parts_status;

  if(!np_beta_order_supported(order)) {
    np_beta_set_status(status, NP_BETA_ERR_SCALE);
    return NAN;
  }
  if(order == 2)
    return np_beta_overlap_order2(center_one, bandwidth_one,
                                  center_two, bandwidth_two,
                                  lower, upper, status);

  parts_status = np_beta_overlap_order_log_parts(
    center_one, bandwidth_one, center_two, bandwidth_two,
    lower, upper, order, &positive_log, &negative_log);
  if(parts_status != NP_BETA_OK) {
    np_beta_set_status(status, parts_status);
    return (parts_status == NP_BETA_ERR_RANGE) ? INFINITY : NAN;
  }

  return np_beta_signed_log_difference(positive_log, negative_log, status);
}

double np_beta_log_abs_overlap_order(double center_one,
                                     double bandwidth_one,
                                     double center_two,
                                     double bandwidth_two,
                                     double lower,
                                     double upper,
                                     int order,
                                     int *sign,
                                     np_beta_status *status)
{
  double positive_log = -INFINITY;
  double negative_log = -INFINITY;
  double log_absolute = -INFINITY;
  np_beta_status parts_status;

  if(sign == NULL || !np_beta_order_supported(order)) {
    np_beta_set_status(status, NP_BETA_ERR_SCALE);
    return NAN;
  }
  *sign = 0;

  if(order == 2) {
    log_absolute = np_beta_log_overlap_order2(
      center_one, bandwidth_one, center_two, bandwidth_two,
      lower, upper, &parts_status);
    if(parts_status == NP_BETA_OK && log_absolute != -INFINITY)
      *sign = 1;
    np_beta_set_status(status, parts_status);
    return (parts_status == NP_BETA_OK) ? log_absolute : NAN;
  }

  parts_status = np_beta_overlap_order_log_parts(
    center_one, bandwidth_one, center_two, bandwidth_two,
    lower, upper, order, &positive_log, &negative_log);
  if(parts_status == NP_BETA_OK)
    parts_status = np_beta_signed_log_absolute(
      positive_log, negative_log, &log_absolute, sign);
  np_beta_set_status(status, parts_status);
  return (parts_status == NP_BETA_OK) ? log_absolute : NAN;
}
