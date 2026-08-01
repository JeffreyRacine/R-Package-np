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

/*
 * Optional O(n) factor supplied in signed-log form for one evaluation row.
 * The provider owns neither output buffer and must fill every non-omitted
 * entry through observation_count.  Leave-one-out exclusion remains the row
 * engine's responsibility; omitted_observation is supplied so a provider may
 * leave that entry untouched.
 */
typedef NPContinuousKernelRowStatus (*NPContinuousKernelLogFactorFunction)(
  const void *context,
  int evaluation_index,
  int omitted_observation,
  int observation_count,
  double *log_absolute,
  signed char *sign);

typedef struct {
  NPContinuousKernelLogFactorFunction function;
  const void *context;
} NPContinuousKernelLogFactorProvider;

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

/*
 * Scratch owned by an absolute derivative-row consumer.  The four signed-log
 * channels are sized only by the response/weight tensor width, never by the
 * observation count.  They are reusable across evaluation rows.
 */
typedef struct {
  double *regular_positive_log;
  double *regular_negative_log;
  double *jump_positive_log;
  double *jump_negative_log;
  size_t capacity;
} NPContinuousKernelDerivativeAccumulator;

typedef struct NPContinuousKernelDerivativeDiagnostics {
  int bad_coordinate;
  int bad_observation;
  int undefined_count;
  np_beta_status beta_status;
} NPContinuousKernelDerivativeDiagnostics;

void np_continuous_kernel_row_workspace_init(
  NPContinuousKernelRowWorkspace *workspace);

void np_continuous_kernel_row_workspace_release(
  NPContinuousKernelRowWorkspace *workspace);

NPContinuousKernelRowStatus np_continuous_kernel_row_workspace_reserve(
  NPContinuousKernelRowWorkspace *workspace,
  size_t observation_count,
  int need_secondary);

NPContinuousKernelRowStatus np_continuous_kernel_row_plan_validate(
  const NPContinuousKernelRowPlan *plan);

NPContinuousKernelRowStatus np_continuous_kernel_beta_factor_row(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *result);

/*
 * Compose the canonical beta factor with one family-neutral signed-log
 * factor.  On success, workspace's primary channel is the complete product,
 * result->row is scaled by its complete-row maximum, and
 * result->total_log_scale records that maximum.  Per-segment scales retain
 * their beta-only diagnostic meaning.  The incumbent beta-only entry point
 * remains separate so existing hot paths incur no provider branch.
 */
NPContinuousKernelRowStatus
np_continuous_kernel_beta_factor_row_with_log_factor(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  const NPContinuousKernelLogFactorProvider *provider,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *result);

NPContinuousKernelRowStatus np_continuous_kernel_beta_derivative_factor_row(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  int derivative_coordinate,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelDerivativeRowResult *result);

void np_continuous_kernel_derivative_accumulator_init(
  NPContinuousKernelDerivativeAccumulator *accumulator);

void np_continuous_kernel_derivative_accumulator_release(
  NPContinuousKernelDerivativeAccumulator *accumulator);

NPContinuousKernelRowStatus
np_continuous_kernel_beta_derivative_absolute_rows_validated(
  const NPContinuousKernelRowPlan *plan,
  int leave_one_out,
  int leave_one_out_offset,
  int derivative_coordinate,
  double * const *response,
  int response_columns,
  double * const *case_weights,
  int weight_columns,
  NPContinuousKernelDerivativeAccumulator *accumulator,
  double *weighted_sum,
  double *kernel_weights,
  NPContinuousKernelDerivativeDiagnostics *diagnostics);

NPContinuousKernelRowStatus
np_continuous_kernel_beta_derivative_absolute_rows_with_log_factor_validated(
  const NPContinuousKernelRowPlan *plan,
  int leave_one_out,
  int leave_one_out_offset,
  int derivative_coordinate,
  const NPContinuousKernelLogFactorProvider *provider,
  double * const *response,
  int response_columns,
  double * const *case_weights,
  int weight_columns,
  NPContinuousKernelDerivativeAccumulator *accumulator,
  double *weighted_sum,
  double *kernel_weights,
  NPContinuousKernelDerivativeDiagnostics *diagnostics);

NPContinuousKernelRowStatus
np_continuous_kernel_beta_derivative_powered_rows_validated(
  const NPContinuousKernelRowPlan *plan,
  int leave_one_out,
  int leave_one_out_offset,
  int derivative_coordinate,
  int kernel_power,
  const NPContinuousKernelLogFactorProvider *provider,
  double * const *response,
  int response_columns,
  double * const *case_weights,
  int weight_columns,
  NPContinuousKernelDerivativeAccumulator *accumulator,
  double *weighted_sum,
  double *kernel_weights,
  NPContinuousKernelDerivativeDiagnostics *diagnostics);

NPContinuousKernelRowStatus
np_continuous_kernel_beta_dual_power_rows_validated(
  const NPContinuousKernelRowPlan *plan,
  int leave_one_out,
  int leave_one_out_offset,
  const NPContinuousKernelLogFactorProvider *provider,
  double * const *response,
  int response_columns,
  double * const *case_weights,
  int weight_columns,
  double * const *power2_response,
  int power2_response_columns,
  double * const *power2_case_weights,
  int power2_weight_columns,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *row_result,
  double *weighted_sum,
  double *weighted_sum_power2,
  NPContinuousKernelDerivativeDiagnostics *diagnostics,
  NPContinuousKernelProgressFunction progress);

/*
 * Accumulate the complete canonical row and its online centered second
 * moment in training order.  This is the stable moment contract required by
 * distribution-function consumers; it deliberately does not reconstruct the
 * centered moment by subtracting two raw moments.
 */
NPContinuousKernelRowStatus
np_continuous_kernel_beta_centered_moment_rows_validated(
  const NPContinuousKernelRowPlan *plan,
  int leave_one_out,
  int leave_one_out_offset,
  const NPContinuousKernelLogFactorProvider *provider,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *row_result,
  double *sum,
  double *centered_m2,
  NPContinuousKernelDerivativeDiagnostics *diagnostics,
  NPContinuousKernelProgressFunction progress);

NPContinuousKernelRowStatus np_continuous_kernel_scaled_restore(
  double scaled_value,
  double log_scale,
  int power,
  double *value);

NPContinuousKernelRowStatus np_continuous_kernel_signed_log_restore(
  double log_absolute,
  int sign,
  double *value);

NPContinuousKernelRowStatus np_continuous_kernel_signed_log_power_restore(
  double log_absolute,
  int sign,
  int power,
  double *value);

NPContinuousKernelRowStatus np_continuous_kernel_signed_log_power(
  double log_absolute,
  int sign,
  int power,
  double *powered_log_absolute,
  int *powered_sign);

const char *np_continuous_kernel_row_status_message(
  NPContinuousKernelRowStatus status);

#endif
