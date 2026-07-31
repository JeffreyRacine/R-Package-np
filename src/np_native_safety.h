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
np_native_cache_table_bytes(size_t capacity,
                            size_t key_width,
                            size_t key_element_size,
                            size_t fixed_bytes,
                            size_t *result)
{
  size_t key_count;
  size_t total = fixed_bytes;

  if((result == NULL) || (capacity == 0U) || (key_width == 0U) ||
     (key_element_size == 0U) ||
     !np_size_mul_checked(capacity, key_width, &key_count) ||
     !np_size_accumulate_array_bytes(&total, key_count, key_element_size) ||
     !np_size_accumulate_array_bytes(&total, capacity, sizeof(double)) ||
     !np_size_accumulate_array_bytes(&total,
                                     capacity,
                                     sizeof(unsigned char)))
    return 0;
  *result = total;
  return 1;
}

static inline int
np_native_cache_rehash_peak_bytes(size_t old_capacity,
                                  size_t new_capacity,
                                  size_t key_width,
                                  size_t key_element_size,
                                  size_t fixed_bytes,
                                  size_t *result)
{
  size_t old_bytes;
  size_t new_bytes;
  size_t peak_bytes;

  if((result == NULL) || (old_capacity == 0U) ||
     !np_native_cache_table_bytes(old_capacity,
                                  key_width,
                                  key_element_size,
                                  0U,
                                  &old_bytes) ||
     !np_native_cache_table_bytes(new_capacity,
                                  key_width,
                                  key_element_size,
                                  0U,
                                  &new_bytes) ||
     !np_size_add_checked(old_bytes, new_bytes, &peak_bytes) ||
     !np_size_add_checked(peak_bytes, fixed_bytes, result))
    return 0;
  return 1;
}

static inline int
np_native_cache_growth_size(size_t capacity, size_t *result)
{
  size_t capacity_scaled;
  size_t prospective_size;

  if((result == NULL) || (capacity == 0U) ||
     !np_size_mul_checked(capacity, 7U, &capacity_scaled))
    return 0;
  /*
   * Preserve `(size + 1) * 10 >= capacity * 7`, but compute the
   * overflow-checked current-size threshold once outside the insertion path.
   */
  prospective_size = capacity_scaled/10U;
  if(capacity_scaled%10U != 0U) {
    if(!np_size_add_checked(prospective_size, 1U, &prospective_size))
      return 0;
  }
  if((prospective_size == 0U) || (prospective_size > capacity))
    return 0;
  *result = prospective_size - 1U;
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

int np_int_padded_count_nonnegative(int count,
                                    int partitions,
                                    int minimum_one_per_partition,
                                    int *stride,
                                    int *padded_count);

int np_int64_padded_count_nonnegative(int64_t count,
                                      int partitions,
                                      int minimum_one_per_partition,
                                      int64_t *stride,
                                      int64_t *padded_count);

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

static inline NPNativeAllocStatus
np_native_malloc_column_matrix(double ***result, int nrows, int ncols)
{
  double **columns = NULL;
  double *data = NULL;
  NPNativeAllocStatus status;
  size_t cells;
  int column;

  if(result == NULL)
    return NP_NATIVE_ALLOC_INVALID;
  *result = NULL;
  if((nrows <= 0) || (ncols <= 0))
    return NP_NATIVE_ALLOC_INVALID;
  if(!np_size_mul_checked((size_t)nrows, (size_t)ncols, &cells))
    return NP_NATIVE_ALLOC_OVERFLOW;

  status = np_native_malloc_array((void **)&columns,
                                  (size_t)ncols,
                                  sizeof(*columns));
  if(status != NP_NATIVE_ALLOC_OK)
    return status;
  status = np_native_malloc_array((void **)&data, cells, sizeof(*data));
  if(status != NP_NATIVE_ALLOC_OK) {
    free(columns);
    return status;
  }

  columns[0] = data;
  for(column = 1; column < ncols; column++)
    columns[column] = columns[column - 1] + nrows;
  *result = columns;
  return NP_NATIVE_ALLOC_OK;
}

static inline int
np_native_tile_workspace_bytes(size_t nrows,
                               size_t width,
                               size_t linear_matrix_count,
                               size_t square_matrix_count,
                               size_t fixed_bytes,
                               size_t *result)
{
  size_t linear_cells;
  size_t square_cells;
  size_t data_cells;
  size_t data_bytes;
  size_t pointer_count;
  size_t pointer_bytes;
  size_t total;

  if((result == NULL) || (nrows == 0U) || (width == 0U))
    return 0;
  if(!np_size_mul3_checked(nrows,
                           width,
                           linear_matrix_count,
                           &linear_cells) ||
     !np_size_mul3_checked(width,
                           width,
                           square_matrix_count,
                           &square_cells) ||
     !np_size_add_checked(linear_cells, square_cells, &data_cells) ||
     !np_size_array_bytes_checked(data_cells, sizeof(double), &data_bytes) ||
     !np_size_mul_checked(width,
                          linear_matrix_count,
                          &pointer_count) ||
     !np_size_array_bytes_checked(pointer_count,
                                  sizeof(double *),
                                  &pointer_bytes) ||
     !np_size_add_checked(data_bytes, pointer_bytes, &total) ||
     !np_size_add_checked(total, fixed_bytes, result))
    return 0;
  return 1;
}

static inline int
np_native_bounded_tile_width(size_t nrows,
                             int preferred_width,
                             size_t linear_matrix_count,
                             size_t square_matrix_count,
                             size_t fixed_bytes,
                             size_t byte_ceiling)
{
  int width;

  if((nrows == 0U) || (preferred_width <= 0) ||
     ((linear_matrix_count == 0U) && (square_matrix_count == 0U)) ||
     (byte_ceiling == 0U))
    return 0;

  for(width = preferred_width; width > 0; width /= 2) {
    size_t bytes;

    if(np_native_tile_workspace_bytes(nrows,
                                      (size_t)width,
                                      linear_matrix_count,
                                      square_matrix_count,
                                      fixed_bytes,
                                      &bytes) &&
       (bytes <= byte_ceiling))
      return width;
  }
  return 0;
}

#endif
