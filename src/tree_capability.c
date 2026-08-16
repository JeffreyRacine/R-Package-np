#include <float.h>
#include <math.h>
#include <stddef.h>

#include "tree_capability.h"

static int np_tree_support_endpoint_unbounded(const double value)
{
  return !isfinite(value) || fabs(value) == DBL_MAX;
}

NPPruneSupportDescriptor np_prune_support_descriptor(
  const int effective_kernel,
  const double lower,
  const double upper)
{
  NPPruneSupportDescriptor descriptor = {
    NP_PRUNE_SUPPORT_UNKNOWN,
    lower,
    upper
  };
  int lower_unbounded;
  int upper_unbounded;

  if(effective_kernel < 0 || effective_kernel >= 40 ||
     effective_kernel % 10 == 9 || isnan(lower) || isnan(upper) ||
     lower > upper)
    return descriptor;

  if(effective_kernel == 28 || (lower == 0.0 && upper == 0.0)) {
    descriptor.kind = NP_PRUNE_SUPPORT_UNSUPPORTED;
    return descriptor;
  }

  lower_unbounded = np_tree_support_endpoint_unbounded(lower);
  upper_unbounded = np_tree_support_endpoint_unbounded(upper);

  if(lower_unbounded && upper_unbounded) {
    descriptor.kind = NP_PRUNE_SUPPORT_UNBOUNDED;
    return descriptor;
  }

  if(lower_unbounded != upper_unbounded) {
    descriptor.kind = NP_PRUNE_SUPPORT_ONE_SIDED_TAIL;
    return descriptor;
  }

  if(lower < 0.0 && upper > 0.0)
    descriptor.kind = NP_PRUNE_SUPPORT_COMPACT_ZERO;
  else
    descriptor.kind = NP_PRUNE_SUPPORT_UNSUPPORTED;

  return descriptor;
}

NPTreeCapability np_tree_fixed_coordinate_capability(
  const NPPruneSupportDescriptor descriptor,
  const double bandwidth,
  const NPTreeCoordinateEnvelope bounds)
{
  double inverse_bandwidth;
  double lower_difference;
  double upper_difference;
  double standardized_lower;
  double standardized_upper;
  int lower_unbounded;
  int upper_unbounded;

  if(descriptor.kind == NP_PRUNE_SUPPORT_UNKNOWN ||
     descriptor.kind == NP_PRUNE_SUPPORT_UNSUPPORTED ||
     !isfinite(bandwidth) || bandwidth <= 0.0 ||
     !isfinite(bounds.train_min) || !isfinite(bounds.train_max) ||
     !isfinite(bounds.eval_min) || !isfinite(bounds.eval_max) ||
     bounds.train_min > bounds.train_max ||
     bounds.eval_min > bounds.eval_max)
    return NP_TREE_CAPABILITY_UNKNOWN;

  if(descriptor.kind == NP_PRUNE_SUPPORT_UNBOUNDED)
    return NP_TREE_CANNOT_PRUNE;

  inverse_bandwidth = 1.0/bandwidth;
  if(!isfinite(inverse_bandwidth) || inverse_bandwidth <= 0.0)
    return NP_TREE_CAPABILITY_UNKNOWN;

  lower_difference = bounds.eval_min - bounds.train_max;
  upper_difference = bounds.eval_max - bounds.train_min;
  standardized_lower = lower_difference*inverse_bandwidth;
  standardized_upper = upper_difference*inverse_bandwidth;
  lower_unbounded = np_tree_support_endpoint_unbounded(descriptor.lower);
  upper_unbounded = np_tree_support_endpoint_unbounded(descriptor.upper);

  if(descriptor.kind == NP_PRUNE_SUPPORT_COMPACT_ZERO) {
    if(!isfinite(standardized_lower) ||
       !isfinite(standardized_upper) ||
       standardized_lower <= descriptor.lower ||
       standardized_upper >= descriptor.upper)
      return NP_TREE_CAN_PRUNE_ZERO_SUPPORT;
    return NP_TREE_CANNOT_PRUNE;
  }

  if(descriptor.kind != NP_PRUNE_SUPPORT_ONE_SIDED_TAIL ||
     lower_unbounded == upper_unbounded)
    return NP_TREE_CAPABILITY_UNKNOWN;

  if((!lower_unbounded &&
      (!isfinite(standardized_lower) ||
       standardized_lower <= descriptor.lower)) ||
     (!upper_unbounded &&
      (!isfinite(standardized_upper) ||
       standardized_upper >= descriptor.upper)))
    return NP_TREE_CAN_USE_ONE_SIDED_TAIL;

  return NP_TREE_CANNOT_PRUNE;
}

NPTreeCapability np_tree_capability_combine(
  const NPTreeCapability current,
  const NPTreeCapability coordinate)
{
  if(current == NP_TREE_CAPABILITY_UNKNOWN ||
     coordinate == NP_TREE_CAPABILITY_UNKNOWN)
    return NP_TREE_CAPABILITY_UNKNOWN;
  if(current == NP_TREE_CAN_PRUNE_ZERO_SUPPORT ||
     coordinate == NP_TREE_CAN_PRUNE_ZERO_SUPPORT)
    return NP_TREE_CAN_PRUNE_ZERO_SUPPORT;
  if(current == NP_TREE_CAN_USE_ONE_SIDED_TAIL ||
     coordinate == NP_TREE_CAN_USE_ONE_SIDED_TAIL)
    return NP_TREE_CAN_USE_ONE_SIDED_TAIL;
  return NP_TREE_CANNOT_PRUNE;
}

NPTreeCapability np_tree_fixed_capability(
  const NPPruneSupportDescriptor * const support,
  const double * const bandwidth,
  const NPTreeCoordinateEnvelope * const envelope,
  const int dimensions)
{
  int coordinate;
  NPTreeCapability capability = NP_TREE_CANNOT_PRUNE;

  if(support == NULL || bandwidth == NULL || envelope == NULL ||
     dimensions <= 0)
    return NP_TREE_CAPABILITY_UNKNOWN;

  for(coordinate = 0; coordinate < dimensions; ++coordinate) {
    capability = np_tree_capability_combine(
      capability,
      np_tree_fixed_coordinate_capability(support[coordinate],
                                          bandwidth[coordinate],
                                          envelope[coordinate]));
    if(capability == NP_TREE_CAPABILITY_UNKNOWN ||
       capability == NP_TREE_CAN_PRUNE_ZERO_SUPPORT)
      return capability;
  }

  return capability;
}
