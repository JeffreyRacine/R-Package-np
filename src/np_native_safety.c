/*
 * Out-of-line native-safety operations used by MPI translation units.
 *
 * Keep partition arithmetic here rather than inlining it into jksum.c: that
 * file is exceptionally large and performance-sensitive, so changing its
 * inlining and code layout can create timing drift unrelated to the checked
 * arithmetic itself.
 */

#include "np_native_safety.h"

int
np_int_padded_count_nonnegative(int count,
                                int partitions,
                                int minimum_one_per_partition,
                                int *stride,
                                int *padded_count)
{
  int local_stride;
  size_t padded_size;

  if((stride == NULL) || (padded_count == NULL) ||
     ((minimum_one_per_partition != 0) &&
      (minimum_one_per_partition != 1)) ||
     !np_int_ceil_div_nonnegative(count, partitions, &local_stride))
    return 0;
  if(minimum_one_per_partition && (local_stride == 0))
    local_stride = 1;
  if(!np_size_mul_checked((size_t)local_stride,
                          (size_t)partitions,
                          &padded_size) ||
     !np_size_to_int_checked(padded_size, padded_count))
    return 0;
  *stride = local_stride;
  return 1;
}

int
np_int64_padded_count_nonnegative(int64_t count,
                                  int partitions,
                                  int minimum_one_per_partition,
                                  int64_t *stride,
                                  int64_t *padded_count)
{
  int64_t local_stride;

  if((stride == NULL) || (padded_count == NULL) ||
     (count < 0) || (partitions <= 0) ||
     ((minimum_one_per_partition != 0) &&
      (minimum_one_per_partition != 1)))
    return 0;
  local_stride = count/(int64_t)partitions;
  if(count%(int64_t)partitions != 0)
    local_stride++;
  if(minimum_one_per_partition && (local_stride == 0))
    local_stride = 1;
  if(local_stride > INT64_MAX/(int64_t)partitions)
    return 0;
  *stride = local_stride;
  *padded_count = local_stride*(int64_t)partitions;
  return 1;
}
