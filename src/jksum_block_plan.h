#ifndef NP_JKSUM_BLOCK_PLAN_H
#define NP_JKSUM_BLOCK_PLAN_H

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "np_native_safety.h"

typedef enum {
  NP_JKSUM_BLOCK_PLAN_OK = 0,
  NP_JKSUM_BLOCK_PLAN_INVALID = 1,
  NP_JKSUM_BLOCK_PLAN_OVERFLOW = 2,
  NP_JKSUM_BLOCK_PLAN_CAPACITY = 3
} NPJksumBlockPlanStatus;

typedef struct {
  int64_t x_width;
  int64_t x_alloc;
  int64_t y_width;
  int64_t y_alloc;
  size_t total_cells;
} NPJksumBlockPlan;

static inline NPJksumBlockPlanStatus
np_jksum_memfac_cells(double memfac, size_t *cells)
{
  const size_t ceiling = SIZE_MAX/10U;
  long double scaled;

  if((cells == NULL) || !isfinite(memfac) || (memfac <= 0.0))
    return NP_JKSUM_BLOCK_PLAN_INVALID;
  scaled = ceill((long double)memfac*300000.0L);
  if(!isfinite(scaled) || (scaled >= (long double)ceiling)) {
    *cells = ceiling;
    return NP_JKSUM_BLOCK_PLAN_OK;
  }
  if(scaled < 1.0L)
    return NP_JKSUM_BLOCK_PLAN_INVALID;
  *cells = (size_t)scaled;
  return NP_JKSUM_BLOCK_PLAN_OK;
}

static inline int
np_jksum_block_matrix_fits(size_t train,
                            size_t width,
                            size_t matrix_byte_limit)
{
  size_t cells;
  size_t bytes;

  return np_size_mul_checked(train, width, &cells) &&
    np_size_array_bytes_checked(cells, sizeof(double), &bytes) &&
    (bytes <= matrix_byte_limit);
}

static inline int
np_jksum_distribution_one_full_plan_safe(size_t train,
                                          size_t eval,
                                          size_t budget_cells,
                                          size_t matrix_byte_limit,
                                          size_t *total_cells)
{
  size_t row_cells;

  return np_size_add_checked(train, 1U, &row_cells) &&
    np_size_mul_checked(row_cells, eval, total_cells) &&
    (*total_cells <= budget_cells) &&
    np_jksum_block_matrix_fits(train, eval, matrix_byte_limit);
}

static inline int
np_jksum_distribution_two_full_plan_safe(size_t train,
                                          size_t eval,
                                          size_t budget_cells,
                                          size_t matrix_byte_limit,
                                          size_t *total_cells)
{
  size_t x_cells;
  size_t y_cells;
  size_t x_total;
  size_t y_total;

  return np_size_mul_checked(train, 2U, &y_cells) &&
    np_size_add_checked(y_cells, 1U, &x_cells) &&
    np_size_mul_checked(x_cells, train, &x_total) &&
    np_size_mul_checked(y_cells, eval, &y_total) &&
    np_size_add_checked(x_total, y_total, total_cells) &&
    (*total_cells <= budget_cells) &&
    np_jksum_block_matrix_fits(train, train, matrix_byte_limit) &&
    np_jksum_block_matrix_fits(train, eval, matrix_byte_limit);
}

NPJksumBlockPlanStatus
np_jksum_distribution_one_block_plan(int64_t train_alloc,
                                      int64_t eval_alloc,
                                      int partitions,
                                      size_t budget_cells,
                                      size_t matrix_byte_limit,
                                      NPJksumBlockPlan *plan);

NPJksumBlockPlanStatus
np_jksum_distribution_two_block_plan(int64_t train_alloc,
                                      int64_t eval_alloc,
                                      int partitions,
                                      size_t budget_cells,
                                      size_t matrix_byte_limit,
                                      NPJksumBlockPlan *plan);

const char *
np_jksum_block_plan_status_message(NPJksumBlockPlanStatus status);

#endif
