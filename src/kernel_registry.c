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
