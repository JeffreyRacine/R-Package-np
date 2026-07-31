#include "jksum_block_plan.h"

#include <math.h>
#include <stdint.h>

#include "np_native_safety.h"

static void
np_jksum_block_plan_clear(NPJksumBlockPlan *plan)
{
  if(plan == NULL)
    return;
  plan->x_width = 0;
  plan->x_alloc = 0;
  plan->y_width = 0;
  plan->y_alloc = 0;
  plan->total_cells = 0U;
}

static int
np_jksum_positive_i64_to_size(int64_t value, size_t *result)
{
  if((result == NULL) || (value <= 0) ||
     ((uint64_t)value > (uint64_t)SIZE_MAX))
    return 0;
  *result = (size_t)value;
  return 1;
}

static int
np_jksum_padded_width(int64_t width,
                      int partitions,
                      int64_t *allocated)
{
  int64_t stride;

  return np_int64_padded_count_nonnegative(width,
                                           partitions,
                                           1,
                                           &stride,
                                           allocated);
}

static int
np_jksum_one_plan_cells(size_t train,
                        size_t allocated_width,
                        size_t *cells)
{
  size_t row_cells;

  return np_size_add_checked(train, 1U, &row_cells) &&
    np_size_mul_checked(row_cells, allocated_width, cells);
}

static int
np_jksum_two_plan_denominators(size_t train,
                               size_t *x_cells,
                               size_t *y_cells)
{
  return np_size_mul_checked(train, 2U, y_cells) &&
    np_size_add_checked(*y_cells, 1U, x_cells);
}

static int
np_jksum_two_plan_cells(size_t train,
                        size_t x_alloc,
                        size_t y_alloc,
                        size_t *cells)
{
  size_t x_cells;
  size_t y_cells;
  size_t x_total;
  size_t y_total;

  return np_jksum_two_plan_denominators(train, &x_cells, &y_cells) &&
    np_size_mul_checked(x_cells, x_alloc, &x_total) &&
    np_size_mul_checked(y_cells, y_alloc, &y_total) &&
    np_size_add_checked(x_total, y_total, cells);
}

static int
np_jksum_aligned_capacity(size_t capacity,
                          int partitions,
                          size_t *aligned)
{
  const size_t partition_count = (size_t)partitions;

  if((aligned == NULL) || (partitions <= 0))
    return 0;
  *aligned = (capacity/partition_count)*partition_count;
  return 1;
}

NPJksumBlockPlanStatus
np_jksum_distribution_one_block_plan(int64_t train_alloc,
                                      int64_t eval_alloc,
                                      int partitions,
                                      size_t budget_cells,
                                      size_t matrix_byte_limit,
                                      NPJksumBlockPlan *plan)
{
  size_t train;
  size_t eval;
  size_t full_cells;
  size_t legacy_width;
  size_t legacy_alloc;
  size_t actual_cells;
  size_t matrix_cell_limit;
  size_t memory_cap;
  size_t matrix_cap;
  size_t safe_cap;
  size_t aligned_cap;
  int64_t legacy_alloc_i64;

  np_jksum_block_plan_clear(plan);
  if((plan == NULL) || (partitions <= 0) ||
     (matrix_byte_limit == 0U) ||
     !np_jksum_positive_i64_to_size(train_alloc, &train) ||
     !np_jksum_positive_i64_to_size(eval_alloc, &eval) ||
     (train_alloc%(int64_t)partitions != 0) ||
     (eval_alloc%(int64_t)partitions != 0))
    return NP_JKSUM_BLOCK_PLAN_INVALID;

  if(np_jksum_distribution_one_full_plan_safe(train,
                                               eval,
                                               budget_cells,
                                               matrix_byte_limit,
                                               &full_cells)) {
    plan->x_width = eval_alloc;
    plan->x_alloc = eval_alloc;
    plan->total_cells = full_cells;
    return NP_JKSUM_BLOCK_PLAN_OK;
  }

  if(!np_size_add_checked(train, 1U, &memory_cap))
    return NP_JKSUM_BLOCK_PLAN_OVERFLOW;
  legacy_width = budget_cells/memory_cap;
  if(legacy_width > eval)
    legacy_width = eval;
  if((legacy_width > 0U) &&
     (legacy_width <= (size_t)INT64_MAX) &&
     np_jksum_padded_width((int64_t)legacy_width,
                           partitions,
                           &legacy_alloc_i64) &&
     (legacy_alloc_i64 == (int64_t)legacy_width) &&
     np_jksum_positive_i64_to_size(legacy_alloc_i64, &legacy_alloc) &&
     np_jksum_one_plan_cells(train, legacy_alloc, &actual_cells) &&
     (actual_cells <= budget_cells) &&
     np_jksum_block_matrix_fits(train,
                                legacy_alloc,
                                matrix_byte_limit)) {
    plan->x_width = (int64_t)legacy_width;
    plan->x_alloc = legacy_alloc_i64;
    plan->total_cells = actual_cells;
    return NP_JKSUM_BLOCK_PLAN_OK;
  }

  memory_cap = budget_cells/memory_cap;
  matrix_cell_limit = matrix_byte_limit/sizeof(double);
  matrix_cap = matrix_cell_limit/train;
  safe_cap = memory_cap < matrix_cap ? memory_cap : matrix_cap;
  if(!np_jksum_aligned_capacity(safe_cap, partitions, &aligned_cap) ||
     (aligned_cap < (size_t)partitions))
    return NP_JKSUM_BLOCK_PLAN_CAPACITY;
  if(aligned_cap > eval)
    aligned_cap = eval;
  if((aligned_cap == 0U) || (aligned_cap > (size_t)INT64_MAX) ||
     !np_jksum_one_plan_cells(train, aligned_cap, &actual_cells) ||
     (actual_cells > budget_cells) ||
     !np_jksum_block_matrix_fits(train,
                                 aligned_cap,
                                 matrix_byte_limit))
    return NP_JKSUM_BLOCK_PLAN_OVERFLOW;

  plan->x_width = (int64_t)aligned_cap;
  plan->x_alloc = (int64_t)aligned_cap;
  plan->total_cells = actual_cells;
  return NP_JKSUM_BLOCK_PLAN_OK;
}

NPJksumBlockPlanStatus
np_jksum_distribution_two_block_plan(int64_t train_alloc,
                                      int64_t eval_alloc,
                                      int partitions,
                                      size_t budget_cells,
                                      size_t matrix_byte_limit,
                                      NPJksumBlockPlan *plan)
{
  size_t train;
  size_t eval;
  size_t x_cells;
  size_t y_cells;
  size_t full_cells;
  size_t y_width;
  size_t x_width;
  size_t y_alloc;
  size_t x_alloc;
  size_t y_used;
  size_t remaining;
  size_t actual_cells;
  size_t matrix_cell_limit;
  size_t matrix_width_cap;
  size_t x_reserve;
  size_t y_cap;
  size_t x_cap;
  size_t aligned_cap;
  int64_t y_alloc_i64;
  int64_t x_alloc_i64;

  np_jksum_block_plan_clear(plan);
  if((plan == NULL) || (partitions <= 0) ||
     (matrix_byte_limit == 0U) ||
     !np_jksum_positive_i64_to_size(train_alloc, &train) ||
     !np_jksum_positive_i64_to_size(eval_alloc, &eval) ||
     (train_alloc%(int64_t)partitions != 0) ||
     (eval_alloc%(int64_t)partitions != 0))
    return NP_JKSUM_BLOCK_PLAN_INVALID;
  if(!np_jksum_two_plan_denominators(train, &x_cells, &y_cells))
    return NP_JKSUM_BLOCK_PLAN_OVERFLOW;

  if(np_jksum_distribution_two_full_plan_safe(train,
                                               eval,
                                               budget_cells,
                                               matrix_byte_limit,
                                               &full_cells)) {
    plan->x_width = train_alloc;
    plan->x_alloc = train_alloc;
    plan->y_width = eval_alloc;
    plan->y_alloc = eval_alloc;
    plan->total_cells = full_cells;
    return NP_JKSUM_BLOCK_PLAN_OK;
  }

  if(budget_cells > x_cells) {
    y_width = (budget_cells - x_cells)/y_cells;
    if(y_width > eval)
      y_width = eval;
    if((y_width > 0U) &&
       np_size_mul_checked(y_cells, y_width, &y_used) &&
       (y_used <= budget_cells)) {
      remaining = budget_cells - y_used;
      x_width = remaining/x_cells;
      if((x_width > 0U) &&
         (x_width <= (size_t)INT64_MAX) &&
         (y_width <= (size_t)INT64_MAX) &&
         np_jksum_padded_width((int64_t)x_width,
                               partitions,
                               &x_alloc_i64) &&
         np_jksum_padded_width((int64_t)y_width,
                               partitions,
                               &y_alloc_i64) &&
         np_jksum_positive_i64_to_size(x_alloc_i64, &x_alloc) &&
         np_jksum_positive_i64_to_size(y_alloc_i64, &y_alloc) &&
         np_jksum_two_plan_cells(train,
                                  x_alloc,
                                  y_alloc,
                                  &actual_cells) &&
         (actual_cells <= budget_cells) &&
         np_jksum_block_matrix_fits(train,
                                    x_alloc,
                                    matrix_byte_limit) &&
         np_jksum_block_matrix_fits(train,
                                    y_alloc,
                                    matrix_byte_limit)) {
        plan->x_width = (int64_t)x_width;
        plan->x_alloc = x_alloc_i64;
        plan->y_width = (int64_t)y_width;
        plan->y_alloc = y_alloc_i64;
        plan->total_cells = actual_cells;
        return NP_JKSUM_BLOCK_PLAN_OK;
      }
    }
  }

  matrix_cell_limit = matrix_byte_limit/sizeof(double);
  matrix_width_cap = matrix_cell_limit/train;
  if(!np_size_mul_checked(x_cells, (size_t)partitions, &x_reserve) ||
     (x_reserve > budget_cells))
    return NP_JKSUM_BLOCK_PLAN_CAPACITY;
  y_cap = (budget_cells - x_reserve)/y_cells;
  if(y_cap > matrix_width_cap)
    y_cap = matrix_width_cap;
  if(!np_jksum_aligned_capacity(y_cap, partitions, &aligned_cap) ||
     (aligned_cap < (size_t)partitions))
    return NP_JKSUM_BLOCK_PLAN_CAPACITY;
  y_alloc = aligned_cap > eval ? eval : aligned_cap;

  if(!np_size_mul_checked(y_cells, y_alloc, &y_used) ||
     (y_used > budget_cells))
    return NP_JKSUM_BLOCK_PLAN_OVERFLOW;
  remaining = budget_cells - y_used;
  x_cap = remaining/x_cells;
  if(x_cap > matrix_width_cap)
    x_cap = matrix_width_cap;
  if(!np_jksum_aligned_capacity(x_cap, partitions, &aligned_cap) ||
     (aligned_cap < (size_t)partitions))
    return NP_JKSUM_BLOCK_PLAN_CAPACITY;
  x_alloc = aligned_cap > train ? train : aligned_cap;

  if((x_alloc == 0U) || (y_alloc == 0U) ||
     (x_alloc > (size_t)INT64_MAX) ||
     (y_alloc > (size_t)INT64_MAX) ||
     !np_jksum_two_plan_cells(train, x_alloc, y_alloc, &actual_cells) ||
     (actual_cells > budget_cells) ||
     !np_jksum_block_matrix_fits(train, x_alloc, matrix_byte_limit) ||
     !np_jksum_block_matrix_fits(train, y_alloc, matrix_byte_limit))
    return NP_JKSUM_BLOCK_PLAN_OVERFLOW;

  plan->x_width = (int64_t)x_alloc;
  plan->x_alloc = (int64_t)x_alloc;
  plan->y_width = (int64_t)y_alloc;
  plan->y_alloc = (int64_t)y_alloc;
  plan->total_cells = actual_cells;
  return NP_JKSUM_BLOCK_PLAN_OK;
}

const char *
np_jksum_block_plan_status_message(NPJksumBlockPlanStatus status)
{
  switch(status) {
    case NP_JKSUM_BLOCK_PLAN_OK:
      return "success";
    case NP_JKSUM_BLOCK_PLAN_INVALID:
      return "invalid block-planning input";
    case NP_JKSUM_BLOCK_PLAN_OVERFLOW:
      return "block-planning size overflow";
    case NP_JKSUM_BLOCK_PLAN_CAPACITY:
      return "memory capacity cannot hold one rank-symmetric block row";
  }
  return "unknown block-planning failure";
}
