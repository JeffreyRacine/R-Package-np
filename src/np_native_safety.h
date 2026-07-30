#ifndef NP_NATIVE_SAFETY_H
#define NP_NATIVE_SAFETY_H

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

typedef enum {
  NP_NATIVE_ALLOC_OK = 0,
  NP_NATIVE_ALLOC_OVERFLOW = 1,
  NP_NATIVE_ALLOC_UNAVAILABLE = 2,
  NP_NATIVE_ALLOC_INVALID = 3
} NPNativeAllocStatus;

static inline int
np_size_add_checked(size_t left, size_t right, size_t *result)
{
  if((result == NULL) || (left > SIZE_MAX - right))
    return 0;
  *result = left + right;
  return 1;
}

static inline int
np_size_mul_checked(size_t left, size_t right, size_t *result)
{
  if((result == NULL) || ((left != 0U) && (right > SIZE_MAX/left)))
    return 0;
  *result = left*right;
  return 1;
}

static inline int
np_size_mul3_checked(size_t first,
                     size_t second,
                     size_t third,
                     size_t *result)
{
  size_t partial;

  return np_size_mul_checked(first, second, &partial) &&
    np_size_mul_checked(partial, third, result);
}

static inline int
np_size_array_bytes_checked(size_t count,
                            size_t element_size,
                            size_t *result)
{
  return np_size_mul_checked(count, element_size, result);
}

static inline int
np_size_accumulate_array_bytes(size_t *total,
                               size_t count,
                               size_t element_size)
{
  size_t bytes;
  size_t updated;

  if((total == NULL) ||
     !np_size_array_bytes_checked(count, element_size, &bytes) ||
     !np_size_add_checked(*total, bytes, &updated))
    return 0;
  *total = updated;
  return 1;
}

static inline int
np_size_to_int_checked(size_t value, int *result)
{
  if((result == NULL) || (value > (size_t)INT_MAX))
    return 0;
  *result = (int)value;
  return 1;
}

static inline int
np_int_ceil_div_nonnegative(int numerator, int denominator, int *result)
{
  int quotient;

  if((result == NULL) || (numerator < 0) || (denominator <= 0))
    return 0;
  quotient = numerator/denominator;
  if(numerator%denominator != 0) {
    if(quotient == INT_MAX)
      return 0;
    quotient++;
  }
  *result = quotient;
  return 1;
}

static inline int
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

static inline NPNativeAllocStatus
np_native_malloc_array(void **result, size_t count, size_t element_size)
{
  size_t bytes;
  void *allocation;

  if(result == NULL)
    return NP_NATIVE_ALLOC_INVALID;
  *result = NULL;
  if(!np_size_array_bytes_checked(count, element_size, &bytes))
    return NP_NATIVE_ALLOC_OVERFLOW;
  if(bytes == 0U)
    return NP_NATIVE_ALLOC_OK;
  allocation = malloc(bytes);
  if(allocation == NULL)
    return NP_NATIVE_ALLOC_UNAVAILABLE;
  *result = allocation;
  return NP_NATIVE_ALLOC_OK;
}

static inline NPNativeAllocStatus
np_native_calloc_array(void **result, size_t count, size_t element_size)
{
  size_t bytes;
  void *allocation;

  if(result == NULL)
    return NP_NATIVE_ALLOC_INVALID;
  *result = NULL;
  if(!np_size_array_bytes_checked(count, element_size, &bytes))
    return NP_NATIVE_ALLOC_OVERFLOW;
  if(bytes == 0U)
    return NP_NATIVE_ALLOC_OK;
  allocation = calloc(count, element_size);
  if(allocation == NULL)
    return NP_NATIVE_ALLOC_UNAVAILABLE;
  *result = allocation;
  return NP_NATIVE_ALLOC_OK;
}

static inline NPNativeAllocStatus
np_native_realloc_array(void **target, size_t count, size_t element_size)
{
  size_t bytes;
  void *replacement;

  if(target == NULL)
    return NP_NATIVE_ALLOC_INVALID;
  if(!np_size_array_bytes_checked(count, element_size, &bytes))
    return NP_NATIVE_ALLOC_OVERFLOW;
  if(bytes == 0U) {
    free(*target);
    *target = NULL;
    return NP_NATIVE_ALLOC_OK;
  }
  replacement = realloc(*target, bytes);
  if(replacement == NULL)
    return NP_NATIVE_ALLOC_UNAVAILABLE;
  *target = replacement;
  return NP_NATIVE_ALLOC_OK;
}

#endif
