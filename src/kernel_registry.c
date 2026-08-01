#include <math.h>
#include <stddef.h>

#include "kernel_registry.h"

/* Compile-time uniqueness and ABI sentinels without requiring C11. */
typedef char np_ckernel_family_codes_must_differ[
  (NP_CKERNEL_FAMILY_LEGACY != NP_CKERNEL_FAMILY_BETA) ? 1 : -1];
typedef char np_ckernel_coordinate_code_must_be_nonlegacy[
  (NP_CKERNEL_COORDINATE_CODE < NP_CKERNEL_LEGACY_CODE_MIN) ? 1 : -1];

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

const char *np_continuous_kernel_route_status_message(
  np_continuous_kernel_route_status status)
{
  switch(status) {
  case NP_CKERNEL_ROUTE_OK:
    return "success";
  case NP_CKERNEL_ROUTE_ERR_NULL:
    return "continuous-kernel route is missing";
  case NP_CKERNEL_ROUTE_ERR_SEGMENT_COUNT:
    return "continuous-kernel route has an invalid segment count";
  case NP_CKERNEL_ROUTE_ERR_COORDINATES:
    return "continuous-kernel route does not cover its coordinates exactly";
  case NP_CKERNEL_ROUTE_ERR_DESCRIPTOR:
    return "continuous-kernel route has an invalid descriptor";
  case NP_CKERNEL_ROUTE_ERR_BOUNDS:
    return "beta continuous-kernel route has invalid bounds";
  default:
    return "unknown continuous-kernel route status";
  }
}
