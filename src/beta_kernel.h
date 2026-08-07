#ifndef NP_BETA_KERNEL_H
#define NP_BETA_KERNEL_H

#include <math.h>
#include <stddef.h>

#include <R_ext/Arith.h>

#define NP_BETA_ORDER_MAX_COMPONENTS 4

/*
 * Coordinate-aware associated beta-kernel primitives.
 *
 * These functions intentionally do not use the legacy unary continuous-kernel
 * registry.  A beta kernel needs the evaluation coordinate, observation
 * coordinate, finite support, and public metric bandwidth separately.
 */

typedef enum {
  NP_BETA_OK = 0,
  NP_BETA_ERR_NONFINITE = 1,
  NP_BETA_ERR_BOUNDS = 2,
  NP_BETA_ERR_BANDWIDTH = 3,
  NP_BETA_ERR_SCALE = 4,
  NP_BETA_ERR_OBSERVATION = 5,
  NP_BETA_ERR_CONCENTRATION = 6,
  NP_BETA_ERR_NUMERIC = 7,
  NP_BETA_ERR_RANGE = 8
} np_beta_status;

typedef enum {
  NP_BETA_EVAL_BELOW = -1,
  NP_BETA_EVAL_INSIDE = 0,
  NP_BETA_EVAL_ABOVE = 1
} np_beta_eval_location;

typedef struct {
  double lower;
  double upper;
  double support_length;
  double target_unit;
  double target_complement_unit;
  double observation_unit;
  double observation_complement_unit;
  double log_observation_unit;
  double log_observation_complement_unit;
  double concentration;
  np_beta_eval_location eval_location;
  int observation_endpoint;
} np_beta_shape;

/* PDF-only decomposition of np_beta_shape.  Observation state is immutable
 * within one native invocation; component state is immutable within one
 * evaluation row for fixed and generalized-NN bandwidths. */
typedef struct {
  double log_unit;
  double log_complement_unit;
  int endpoint;
} np_beta_pdf_observation;

/* CDF-only observation state.  The observation centres the beta shape, so
 * these two coordinates are invariant for every evaluation row and bandwidth
 * topology. */
typedef struct {
  double target_unit;
  double target_complement_unit;
} np_beta_cdf_observation;

#if defined(__clang__) || defined(__GNUC__)
# define NP_BETA_ALWAYS_INLINE static inline __attribute__((always_inline))
#else
# define NP_BETA_ALWAYS_INLINE static inline
#endif

NP_BETA_ALWAYS_INLINE double np_beta_observation_unit_coordinate(
  double value,
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

NP_BETA_ALWAYS_INLINE double np_beta_observation_complement_unit_coordinate(
  double observation,
  double upper,
  double support_length)
{
  return (upper - observation) / support_length;
}

NP_BETA_ALWAYS_INLINE double np_beta_observation_log_unit(
  double observation,
  double lower,
  double support_length)
{
  return log(observation - lower) - log(support_length);
}

NP_BETA_ALWAYS_INLINE double np_beta_observation_log_complement_unit(
  double observation,
  double upper,
  double support_length)
{
  return log(upper - observation) - log(support_length);
}

NP_BETA_ALWAYS_INLINE np_beta_status np_beta_pdf_observation_init(
  double observation,
  double lower,
  double upper,
  double support_length,
  double *observation_unit,
  double *observation_complement_unit,
  double *log_observation_unit,
  double *log_observation_complement_unit,
  int *observation_endpoint)
{
  if(observation_unit == NULL || observation_complement_unit == NULL ||
     log_observation_unit == NULL ||
     log_observation_complement_unit == NULL ||
     observation_endpoint == NULL)
    return NP_BETA_ERR_NUMERIC;
  if(!R_FINITE(observation) || !R_FINITE(lower) || !R_FINITE(upper))
    return NP_BETA_ERR_NONFINITE;
  if(!R_FINITE(support_length) || support_length <= 0.0)
    return NP_BETA_ERR_BOUNDS;
  if(observation < lower || observation > upper)
    return NP_BETA_ERR_OBSERVATION;

  *observation_unit = np_beta_observation_unit_coordinate(
    observation, lower, upper, support_length);
  *observation_complement_unit =
    np_beta_observation_complement_unit_coordinate(
      observation, upper, support_length);
  if(!R_FINITE(*observation_unit) ||
     *observation_unit < 0.0 || *observation_unit > 1.0 ||
     !R_FINITE(*observation_complement_unit) ||
     *observation_complement_unit < 0.0 ||
     *observation_complement_unit > 1.0)
    return NP_BETA_ERR_OBSERVATION;
  if(observation == lower) {
    *observation_endpoint = -1;
    *log_observation_unit = -INFINITY;
    *log_observation_complement_unit = 0.0;
  } else if(observation == upper) {
    *observation_endpoint = 1;
    *log_observation_unit = 0.0;
    *log_observation_complement_unit = -INFINITY;
  } else {
    *observation_endpoint = 0;
    *log_observation_unit = np_beta_observation_log_unit(
      observation, lower, support_length);
    *log_observation_complement_unit =
      np_beta_observation_log_complement_unit(
        observation, upper, support_length);
  }
  return NP_BETA_OK;
}

#undef NP_BETA_ALWAYS_INLINE

typedef struct {
  double log_support_length;
  double alpha_minus_one;
  double beta_minus_one;
  double alpha;
  double beta;
  double log_beta;
  double log_abs_coefficient;
  np_beta_eval_location eval_location;
  int coefficient_sign;
  int normalizer_ready;
} np_beta_pdf_component;

/* Derivative-only state prepared once per fixed/GNN evaluation row.  Keeping
 * it separate prevents the scalar/adaptive PDF component from carrying cold
 * fields through its O(n_train * n_eval) stack path. */
typedef struct {
  double support_length;
  double concentration;
  double digamma_alpha;
  double digamma_beta;
} np_beta_pdf_derivative_component;

/*
 * Derivatives of an observation-centred beta kernel can have an ordinary
 * one-sided derivative and, only when an observation exactly matches a
 * support endpoint, a jump at that endpoint.  Keeping both pieces in signed
 * log form lets ratio estimators cancel endpoint jumps before deciding
 * whether their public derivative is finite or infinite.
 */
typedef struct {
  double regular_log_absolute;
  double jump_log_absolute;
  int regular_sign;
  int jump_sign;
} np_beta_derivative;

np_beta_status np_beta_shape_init(double evaluation,
                                  double observation,
                                  double bandwidth,
                                  double lower,
                                  double upper,
                                  int scale,
                                  np_beta_shape *shape);

void np_beta_pdf_observation_from_shape(
  const np_beta_shape *shape,
  np_beta_pdf_observation *observation);

np_beta_status np_beta_pdf_component_from_shape(
  const np_beta_shape *shape,
  np_beta_pdf_component *component);

np_beta_status np_beta_pdf_component_prepare_normalizer(
  np_beta_pdf_component *component);

np_beta_status np_beta_pdf_component_prepare_coefficient(
  np_beta_pdf_component *component,
  int order,
  int component_index);

double np_beta_log_pdf_component_prepared(
  const np_beta_pdf_component *component,
  const np_beta_pdf_observation *observation,
  np_beta_status *status);

double np_beta_log_abs_pdf_prepared(
  const np_beta_pdf_component *components,
  const np_beta_pdf_observation *observation,
  int component_count,
  int *sign,
  np_beta_status *status);

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
  np_beta_derivative *derivative);

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
  np_beta_derivative *derivative);

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
  np_beta_status *status);

double np_beta_log_pdf_scale(const np_beta_shape *shape,
                             np_beta_status *status);

double np_beta_log_abs_pdf_order(double evaluation,
                                 double observation,
                                 double bandwidth,
                                 double lower,
                                 double upper,
                                 int order,
                                 int *sign,
                                 np_beta_status *status);

np_beta_status np_beta_pdf_derivative_order(double evaluation,
                                            double observation,
                                            double bandwidth,
                                            double lower,
                                            double upper,
                                            int order,
                                            np_beta_derivative *derivative);

/* Form the signed-log level and its regular/jump target derivative from one
 * component preparation.  The level phase is completed first so native
 * failure precedence remains identical to calling the two scalar primitives
 * in sequence. */
np_beta_status np_beta_log_abs_pdf_derivative_order(
  double evaluation,
  double observation,
  double bandwidth,
  double lower,
  double upper,
  int order,
  double *level_log_absolute,
  int *level_sign,
  np_beta_derivative *derivative);

double np_beta_cdf_order2(double evaluation,
                          double observation,
                          double bandwidth,
                          double lower,
                          double upper,
                          np_beta_status *status);

double np_beta_cdf_order(double evaluation,
                         double observation,
                         double bandwidth,
                         double lower,
                         double upper,
                         int order,
                         np_beta_status *status);

double np_beta_log_abs_cdf_order(double evaluation,
                                 double observation,
                                 double bandwidth,
                                 double lower,
                                 double upper,
                                 int order,
                                 int *sign,
                                 np_beta_status *status);

np_beta_status np_beta_cdf_observation_init(
  double observation,
  double lower,
  double upper,
  double support_length,
  np_beta_cdf_observation *prepared);

np_beta_status np_beta_concentration_prepare(
  double bandwidth,
  double support_length,
  int scale,
  double *concentration);

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
  np_beta_status *status);

double np_beta_log_abs_cdf_order_prepared_concentration(
  double evaluation,
  double bandwidth,
  double lower,
  double upper,
  int order,
  const np_beta_cdf_observation *observation,
  np_beta_status observation_status,
  const double *log_abs_coefficient,
  const signed char *coefficient_sign,
  const double *prepared_concentration,
  const np_beta_status *prepared_concentration_status,
  int *sign,
  np_beta_status *status);

double np_beta_log_overlap_order2(double center_one,
                                  double bandwidth_one,
                                  double center_two,
                                  double bandwidth_two,
                                  double lower,
                                  double upper,
                                  np_beta_status *status);

double np_beta_log_abs_overlap_order(double center_one,
                                     double bandwidth_one,
                                     double center_two,
                                     double bandwidth_two,
                                     double lower,
                                     double upper,
                                     int order,
                                     int *sign,
                                     np_beta_status *status);

int np_beta_order_supported(int order);

np_beta_status np_beta_signed_log_absolute(double positive_log,
                                           double negative_log,
                                           double *log_absolute,
                                           int *sign);

const char *np_beta_status_message(np_beta_status status);

#endif
