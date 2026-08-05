#ifndef NP_CONTINUOUS_KERNEL_ROW_H
#define NP_CONTINUOUS_KERNEL_ROW_H

#include <stddef.h>

#include "beta_kernel.h"
#include "kernel_registry.h"
#include "regression_contract.h"

typedef enum {
  NP_CONTINUOUS_ROW_OK = 0,
  NP_CONTINUOUS_ROW_ERR_LAYOUT = 1,
  NP_CONTINUOUS_ROW_ERR_ROUTE = 2,
  NP_CONTINUOUS_ROW_ERR_MEMORY = 3,
  NP_CONTINUOUS_ROW_ERR_KERNEL = 4,
  NP_CONTINUOUS_ROW_ERR_NUMERIC = 5,
  NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT = 6
} NPContinuousKernelRowStatus;

#define NP_BETA_PREPARED_MAX_COMPONENTS NP_BETA_ORDER_MAX_COMPONENTS

typedef struct {
  /* R transient allocation scope: nested owners must release in LIFO order.
   * coordinate_slot maps global continuous coordinates to dense beta-only
   * storage and is -1 for every non-beta coordinate.  pdf_active owns the
   * invocation-invariant observation plane for every normal/PDF bandwidth
   * topology.  pdf_row_component_active additionally identifies fixed/GNN
   * rows whose complete component state is evaluation-row invariant.
   * cdf_active owns the invocation-invariant observation plane for every CDF
   * bandwidth topology.  cdf_concentration_active additionally owns one
   * coordinate/component concentration table for fixed bandwidths only; it
   * never owns an observation or evaluation plane. */
  int pdf_active;
  int pdf_row_component_active;
  int cdf_active;
  int cdf_concentration_active;
  int allocation_active;
  int num_train;
  int num_continuous;
  int num_beta_coordinates;
  int num_cdf_coordinates;
  const void *allocation_marker;
  int *coordinate_slot;
  int *cdf_coordinate_slot;
  np_beta_pdf_observation *pdf_observation;
  np_beta_status *pdf_observation_status;
  np_beta_pdf_component *pdf_row_component;
  np_beta_pdf_derivative_component *pdf_row_derivative_component;
  double *pdf_log_abs_coefficient;
  signed char *pdf_coefficient_sign;
  int *pdf_first_interior;
  int *pdf_second_interior;
  np_beta_cdf_observation *cdf_observation;
  np_beta_status *cdf_observation_status;
  double *cdf_log_abs_coefficient;
  signed char *cdf_coefficient_sign;
  double *cdf_concentration;
  np_beta_status *cdf_concentration_status;
} NPContinuousKernelBetaPreparedContext;

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
  NPContinuousKernelBetaPreparedContext *beta_prepared;
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

/*
 * Reusable O(n) signed-log channels for a level row and the regular/jump
 * parts of one target derivative.  The row owner fills all three channels
 * against one common maximum so ratio consumers never restore absolute beta
 * scale or multiply separately rounded scaled rows.
 */
typedef struct {
  double *level_log_absolute;
  double *regular_log_absolute;
  double *jump_log_absolute;
  signed char *level_sign;
  signed char *regular_sign;
  signed char *jump_sign;
  size_t capacity;
} NPContinuousKernelLevelDerivativeWorkspace;

typedef struct NPContinuousKernelDerivativeDiagnostics {
  int bad_coordinate;
  int bad_observation;
  int undefined_count;
  np_beta_status beta_status;
} NPContinuousKernelDerivativeDiagnostics;

void np_continuous_kernel_row_workspace_init(
  NPContinuousKernelRowWorkspace *workspace);

void np_continuous_kernel_beta_prepared_context_init(
  NPContinuousKernelBetaPreparedContext *context);

void np_continuous_kernel_beta_prepared_context_release(
  NPContinuousKernelBetaPreparedContext *context);

NPContinuousKernelRowStatus np_continuous_kernel_beta_prepared_context_prepare(
  NPContinuousKernelBetaPreparedContext *context,
  const NPContinuousKernelRowPlan *plan);

NPContinuousKernelRowStatus
np_continuous_kernel_beta_prepared_derivative_context_prepare(
  NPContinuousKernelBetaPreparedContext *context);

NPContinuousKernelRowStatus
np_continuous_kernel_beta_prepared_pdf_row_prepare(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int evaluation_index,
  int omitted_observation,
  NPContinuousKernelRowResult *result);

NPContinuousKernelRowStatus
np_continuous_kernel_beta_prepared_derivative_row_prepare(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int evaluation_index,
  int omitted_observation,
  NPContinuousKernelRowResult *result);

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

void np_continuous_kernel_level_derivative_workspace_init(
  NPContinuousKernelLevelDerivativeWorkspace *workspace);

void np_continuous_kernel_level_derivative_workspace_release(
  NPContinuousKernelLevelDerivativeWorkspace *workspace);

NPContinuousKernelRowStatus
np_continuous_kernel_beta_level_derivative_log_row_validated(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  int derivative_coordinate,
  const NPContinuousKernelLogFactorProvider *provider,
  NPContinuousKernelLevelDerivativeWorkspace *workspace,
  double *common_log_scale,
  NPContinuousKernelDerivativeDiagnostics *diagnostics);

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
  NPContinuousKernelRowWorkspace *factor_workspace,
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
  NPContinuousKernelRowWorkspace *workspace,
  double *row_storage,
  double *factor_log_absolute,
  signed char *factor_sign,
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
  int retain_common_scale,
  double *scaled_kernel_weights,
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

/*
 * Consume one normalized width-one regression row without restoring its
 * common log scale.  Positive rows use weighted Welford accumulation; signed
 * rows use a two-pass centered moment.  Ordinary regression retains the
 * exact squared-weight effective-sample-size factor.
 */
NPContinuousKernelRowStatus
np_continuous_kernel_beta_regression_moment_rows_validated(
  const NPContinuousKernelRowPlan *plan,
  int leave_one_out,
  int leave_one_out_offset,
  const NPContinuousKernelLogFactorProvider *provider,
  const double *response,
  int positive_weights,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *row_result,
  double *mean,
  double *mean_stderr,
  NPContinuousKernelDerivativeDiagnostics *diagnostics,
  NPContinuousKernelProgressFunction progress);

/* Conditional scalar sibling with the established observation-influence
 * standard-error contract.  Row construction and moment arithmetic remain
 * canonical; the separate entry keeps ordinary regression's hot body exact. */
NPContinuousKernelRowStatus
np_continuous_kernel_beta_conditional_moment_rows_validated(
  const NPContinuousKernelRowPlan *plan,
  int leave_one_out,
  int leave_one_out_offset,
  const NPContinuousKernelLogFactorProvider *provider,
  const double *response,
  int positive_weights,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *row_result,
  double *mean,
  double *mean_stderr,
  NPContinuousKernelDerivativeDiagnostics *diagnostics,
  NPContinuousKernelProgressFunction progress);

NPContinuousKernelRowStatus np_continuous_kernel_scaled_restore(
  double scaled_value,
  double log_scale,
  int power,
  double *value);

/* Derivative outputs may have mathematically defined signed-infinite endpoint
 * limits or an explicit R missing value for an undefined derivative
 * cancellation/standard error. Keep this contract separate from finite-only
 * row/objective restoration so no non-derivative consumer can admit them. */
NPContinuousKernelRowStatus np_continuous_kernel_scaled_derivative_restore(
  double scaled_value,
  double log_scale,
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
