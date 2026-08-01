#ifndef NP_CONTINUOUS_KERNEL_ROW_H
#define NP_CONTINUOUS_KERNEL_ROW_H

#include <stddef.h>

#include "beta_kernel.h"
#include "kernel_registry.h"

typedef enum {
  NP_CONTINUOUS_ROW_OK = 0,
  NP_CONTINUOUS_ROW_ERR_LAYOUT = 1,
  NP_CONTINUOUS_ROW_ERR_ROUTE = 2,
  NP_CONTINUOUS_ROW_ERR_MEMORY = 3,
  NP_CONTINUOUS_ROW_ERR_KERNEL = 4,
  NP_CONTINUOUS_ROW_ERR_NUMERIC = 5
} NPContinuousKernelRowStatus;

typedef struct {
  double *primary_log_absolute;
  double *secondary_log_absolute;
  signed char *primary_sign;
  signed char *secondary_sign;
  size_t capacity;
  size_t secondary_capacity;
} NPContinuousKernelRowWorkspace;

typedef struct {
  const NPContinuousKernelRoute *route;
  int bandwidth_mode;
  int num_train;
  int num_eval;
  int num_continuous;
  int train_is_eval;
  double * const *train;
  double * const *evaluation;
  double * const *bandwidth_eval;
  double * const *bandwidth_train;
  const int *operator;
} NPContinuousKernelRowPlan;

typedef struct {
  double *row;
  double segment_log_scale[NP_CKERNEL_ROUTE_MAX_SEGMENTS];
  double total_log_scale;
  int bad_coordinate;
  int bad_observation;
  np_beta_status beta_status;
} NPContinuousKernelRowResult;

typedef struct {
  double *regular_row;
  double *jump_row;
  double regular_segment_log_scale[NP_CKERNEL_ROUTE_MAX_SEGMENTS];
  double jump_segment_log_scale[NP_CKERNEL_ROUTE_MAX_SEGMENTS];
  double regular_total_log_scale;
  double jump_total_log_scale;
  int bad_coordinate;
  int bad_observation;
  np_beta_status beta_status;
} NPContinuousKernelDerivativeRowResult;

void np_continuous_kernel_row_workspace_init(
  NPContinuousKernelRowWorkspace *workspace);

void np_continuous_kernel_row_workspace_release(
  NPContinuousKernelRowWorkspace *workspace);

NPContinuousKernelRowStatus np_continuous_kernel_row_workspace_reserve(
  NPContinuousKernelRowWorkspace *workspace,
  size_t observation_count,
  int need_secondary);

NPContinuousKernelRowStatus np_continuous_kernel_beta_factor_row(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *result);

NPContinuousKernelRowStatus np_continuous_kernel_beta_derivative_factor_row(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  int derivative_coordinate,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelDerivativeRowResult *result);

NPContinuousKernelRowStatus np_continuous_kernel_scaled_restore(
  double scaled_value,
  double log_scale,
  int power,
  double *value);

NPContinuousKernelRowStatus np_continuous_kernel_signed_log_restore(
  double log_absolute,
  int sign,
  double *value);

const char *np_continuous_kernel_row_status_message(
  NPContinuousKernelRowStatus status);

#endif
