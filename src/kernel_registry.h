#ifndef NP_KERNEL_REGISTRY_H
#define NP_KERNEL_REGISTRY_H

/* Legacy continuous-kernel numeric codes remain 0:9. */
#define NP_CKERNEL_LEGACY_CODE_MIN 0
#define NP_CKERNEL_LEGACY_CODE_MAX 9
#define NP_CKERNEL_COORDINATE_CODE (-1)

typedef enum {
  NP_CKERNEL_FAMILY_LEGACY = 0,
  NP_CKERNEL_FAMILY_BETA = 1
} np_continuous_kernel_family;

typedef enum {
  NP_CKERNEL_DESCRIPTOR_OK = 0,
  NP_CKERNEL_DESCRIPTOR_ERR_FAMILY = 1,
  NP_CKERNEL_DESCRIPTOR_ERR_CODE = 2,
  NP_CKERNEL_DESCRIPTOR_ERR_ORDER = 3
} np_continuous_kernel_descriptor_status;

typedef enum {
  NP_CONTINUOUS_KERNEL_SCALAR_OK = 0,
  NP_CONTINUOUS_KERNEL_SCALAR_ERR_LAYOUT = 1,
  NP_CONTINUOUS_KERNEL_SCALAR_ERR_KERNEL = 2,
  NP_CONTINUOUS_KERNEL_SCALAR_ERR_NUMERIC = 3
} NPContinuousKernelScalarStatus;

typedef struct {
  np_continuous_kernel_family family;
  int legacy_code;
  int order;
} np_continuous_kernel_descriptor;

#define NP_CKERNEL_ROUTE_MAX_SEGMENTS 2

typedef struct {
  np_continuous_kernel_descriptor descriptor;
  int coordinate_offset;
  int coordinate_count;
  const double *lower;
  const double *upper;
} NPContinuousKernelSegment;

typedef struct {
  int segment_count;
  NPContinuousKernelSegment segment[NP_CKERNEL_ROUTE_MAX_SEGMENTS];
} NPContinuousKernelRoute;

/* Optional estimator-owned observation-progress hook for canonical rows. */
typedef void (*NPContinuousKernelProgressFunction)(int completed, int total);

typedef enum {
  NP_CKERNEL_ROUTE_OK = 0,
  NP_CKERNEL_ROUTE_ERR_NULL = 1,
  NP_CKERNEL_ROUTE_ERR_SEGMENT_COUNT = 2,
  NP_CKERNEL_ROUTE_ERR_COORDINATES = 3,
  NP_CKERNEL_ROUTE_ERR_DESCRIPTOR = 4,
  NP_CKERNEL_ROUTE_ERR_BOUNDS = 5
} np_continuous_kernel_route_status;

np_continuous_kernel_descriptor_status
np_continuous_kernel_descriptor_init(int family,
                                     int code,
                                     int order,
                                     np_continuous_kernel_descriptor *descriptor);

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
                                int *sign);

const char *np_continuous_kernel_descriptor_status_message(
  np_continuous_kernel_descriptor_status status);

np_continuous_kernel_route_status
np_continuous_kernel_route_validate(const NPContinuousKernelRoute *route,
                                    int coordinate_count);

const NPContinuousKernelSegment *
np_continuous_kernel_route_segment(const NPContinuousKernelRoute *route,
                                   int coordinate);

int np_continuous_kernel_route_has_beta(
  const NPContinuousKernelRoute *route);

const char *np_continuous_kernel_route_status_message(
  np_continuous_kernel_route_status status);

#endif
