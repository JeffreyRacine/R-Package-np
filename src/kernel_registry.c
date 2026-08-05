#include <float.h>
#include <math.h>
#include <stddef.h>

#include <R_ext/Arith.h>

#include "beta_kernel.h"
#include "kernel_registry.h"

/* Legacy scalar implementations are owned by kernel.c. */
extern double kernel(int kernel_code, double value);
extern double cdf_kernel(int kernel_code, double value);

/* Compile-time uniqueness and ABI sentinels without requiring C11. */
typedef char np_ckernel_family_codes_must_differ[
  (NP_CKERNEL_FAMILY_LEGACY != NP_CKERNEL_FAMILY_BETA) ? 1 : -1];
typedef char np_ckernel_coordinate_code_must_be_nonlegacy[
  (NP_CKERNEL_COORDINATE_CODE < NP_CKERNEL_LEGACY_CODE_MIN) ? 1 : -1];

static int np_continuous_kernel_finite_bound(double value)
{
  return R_FINITE(value) && fabs(value) < 0.5 * DBL_MAX;
}

static double np_continuous_kernel_legacy_pdf(int kernel_code,
                                              double evaluation,
                                              double observation,
                                              double bandwidth,
                                              double lower,
                                              double upper)
{
  const int finite_lower = np_continuous_kernel_finite_bound(lower);
  const int finite_upper = np_continuous_kernel_finite_bound(upper);
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

static double np_continuous_kernel_legacy_cdf(int kernel_code,
                                              double evaluation,
                                              double observation,
                                              double bandwidth,
                                              double lower,
                                              double upper)
{
  const int finite_lower = np_continuous_kernel_finite_bound(lower);
  const int finite_upper = np_continuous_kernel_finite_bound(upper);
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

NPContinuousKernelScalarStatus
np_continuous_kernel_scalar_log(np_continuous_kernel_family family,
                                int kernel_code,
                                int order,
                                int do_cdf,
                                double evaluation,
                                double observation,
                                double bandwidth,
                                double lower,
                                double upper,
                                double *log_absolute,
                                int *sign)
{
  double value;

  if(log_absolute == NULL || sign == NULL)
    return NP_CONTINUOUS_KERNEL_SCALAR_ERR_LAYOUT;
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
    if(scalar_status != NP_BETA_OK)
      return NP_CONTINUOUS_KERNEL_SCALAR_ERR_KERNEL;
  } else if(family == NP_CKERNEL_FAMILY_LEGACY) {
    value = do_cdf ?
      np_continuous_kernel_legacy_cdf(kernel_code, evaluation, observation,
                                      bandwidth, lower, upper) :
      np_continuous_kernel_legacy_pdf(kernel_code, evaluation, observation,
                                      bandwidth, lower, upper);
    if(!R_FINITE(value))
      return NP_CONTINUOUS_KERNEL_SCALAR_ERR_KERNEL;
    if(value != 0.0) {
      *log_absolute = log(fabs(value));
      *sign = (value > 0.0) ? 1 : -1;
    }
  } else {
    return NP_CONTINUOUS_KERNEL_SCALAR_ERR_LAYOUT;
  }

  if(ISNAN(*log_absolute) || *log_absolute == INFINITY)
    return NP_CONTINUOUS_KERNEL_SCALAR_ERR_NUMERIC;
  return NP_CONTINUOUS_KERNEL_SCALAR_OK;
}

np_continuous_kernel_descriptor_status
np_continuous_kernel_descriptor_init(int family,
                                     int code,
                                     int order,
                                     np_continuous_kernel_descriptor *descriptor)
{
  if(descriptor == NULL)
    return NP_CKERNEL_DESCRIPTOR_ERR_FAMILY;

  if(order != 2 && order != 4 && order != 6 && order != 8)
    return NP_CKERNEL_DESCRIPTOR_ERR_ORDER;

  if(family == NP_CKERNEL_FAMILY_LEGACY) {
    if(code < NP_CKERNEL_LEGACY_CODE_MIN ||
       code > NP_CKERNEL_LEGACY_CODE_MAX)
      return NP_CKERNEL_DESCRIPTOR_ERR_CODE;
  } else if(family == NP_CKERNEL_FAMILY_BETA) {
    if(code != NP_CKERNEL_COORDINATE_CODE)
      return NP_CKERNEL_DESCRIPTOR_ERR_CODE;
  } else {
    return NP_CKERNEL_DESCRIPTOR_ERR_FAMILY;
  }

  descriptor->family = (np_continuous_kernel_family) family;
  descriptor->legacy_code = code;
  descriptor->order = order;
  return NP_CKERNEL_DESCRIPTOR_OK;
}

const char *np_continuous_kernel_descriptor_status_message(
  np_continuous_kernel_descriptor_status status)
{
  switch(status) {
  case NP_CKERNEL_DESCRIPTOR_OK:
    return "success";
  case NP_CKERNEL_DESCRIPTOR_ERR_FAMILY:
    return "unknown continuous-kernel family";
  case NP_CKERNEL_DESCRIPTOR_ERR_CODE:
    return "continuous-kernel family/code mismatch";
  case NP_CKERNEL_DESCRIPTOR_ERR_ORDER:
    return "continuous-kernel order must be one of 2, 4, 6, or 8";
  default:
    return "unknown continuous-kernel descriptor status";
  }
}

np_continuous_kernel_route_status
np_continuous_kernel_route_validate(const NPContinuousKernelRoute *route,
                                    int coordinate_count)
{
  int coordinate_end = 0;
  int segment_index;

  if(route == NULL)
    return NP_CKERNEL_ROUTE_ERR_NULL;
  if(route->segment_count < 0 ||
     route->segment_count > NP_CKERNEL_ROUTE_MAX_SEGMENTS)
    return NP_CKERNEL_ROUTE_ERR_SEGMENT_COUNT;
  if(coordinate_count < 0)
    return NP_CKERNEL_ROUTE_ERR_COORDINATES;
  if(route->segment_count == 0)
    return (coordinate_count == 0) ? NP_CKERNEL_ROUTE_OK :
      NP_CKERNEL_ROUTE_ERR_COORDINATES;

  for(segment_index = 0; segment_index < route->segment_count;
      segment_index++) {
    const NPContinuousKernelSegment * const segment =
      &route->segment[segment_index];
    np_continuous_kernel_descriptor descriptor;

    if(segment->coordinate_offset != coordinate_end ||
       segment->coordinate_count <= 0 ||
       segment->coordinate_count > coordinate_count - coordinate_end)
      return NP_CKERNEL_ROUTE_ERR_COORDINATES;
    if(np_continuous_kernel_descriptor_init(
         segment->descriptor.family,
         segment->descriptor.legacy_code,
         segment->descriptor.order,
         &descriptor) != NP_CKERNEL_DESCRIPTOR_OK)
      return NP_CKERNEL_ROUTE_ERR_DESCRIPTOR;
    if(descriptor.family == NP_CKERNEL_FAMILY_BETA) {
      int coordinate;

      if(segment->lower == NULL || segment->upper == NULL)
        return NP_CKERNEL_ROUTE_ERR_BOUNDS;
      for(coordinate = 0; coordinate < segment->coordinate_count;
          coordinate++) {
        const double lower = segment->lower[coordinate];
        const double upper = segment->upper[coordinate];

        if(!isfinite(lower) || !isfinite(upper) || !(upper > lower))
          return NP_CKERNEL_ROUTE_ERR_BOUNDS;
      }
    }
    coordinate_end += segment->coordinate_count;
  }

  return (coordinate_end == coordinate_count) ? NP_CKERNEL_ROUTE_OK :
    NP_CKERNEL_ROUTE_ERR_COORDINATES;
}

const NPContinuousKernelSegment *
np_continuous_kernel_route_segment(const NPContinuousKernelRoute *route,
                                   int coordinate)
{
  int segment_index;

  if(route == NULL || coordinate < 0)
    return NULL;
  for(segment_index = 0; segment_index < route->segment_count;
      segment_index++) {
    const NPContinuousKernelSegment * const segment =
      &route->segment[segment_index];
    const int coordinate_end = segment->coordinate_offset +
      segment->coordinate_count;

    if(coordinate >= segment->coordinate_offset &&
       coordinate < coordinate_end)
      return segment;
  }
  return NULL;
}

int np_continuous_kernel_route_has_beta(
  const NPContinuousKernelRoute *route)
{
  int segment_index;

  if(route == NULL)
    return 0;
  for(segment_index = 0; segment_index < route->segment_count;
      segment_index++)
    if(route->segment[segment_index].descriptor.family ==
       NP_CKERNEL_FAMILY_BETA)
      return 1;
  return 0;
}
