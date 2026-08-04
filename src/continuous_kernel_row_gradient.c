#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include <R_ext/Arith.h>

#include "headers.h"
#include "continuous_kernel_row.h"

#if defined(__clang__) || defined(__GNUC__)
# define NP_CONTINUOUS_ROW_GRADIENT_ALWAYS_INLINE \
  static inline __attribute__((always_inline))
#else
# define NP_CONTINUOUS_ROW_GRADIENT_ALWAYS_INLINE static inline
#endif
void np_continuous_kernel_level_derivative_workspace_init(
  NPContinuousKernelLevelDerivativeWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  workspace->level_log_absolute = NULL;
  workspace->regular_log_absolute = NULL;
  workspace->jump_log_absolute = NULL;
  workspace->level_sign = NULL;
  workspace->regular_sign = NULL;
  workspace->jump_sign = NULL;
  workspace->capacity = 0;
}

void np_continuous_kernel_level_derivative_workspace_release(
  NPContinuousKernelLevelDerivativeWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  free(workspace->level_log_absolute);
  free(workspace->level_sign);
  np_continuous_kernel_level_derivative_workspace_init(workspace);
}

static NPContinuousKernelRowStatus
np_continuous_kernel_level_derivative_workspace_reserve(
  NPContinuousKernelLevelDerivativeWorkspace *workspace,
  size_t observation_count)
{
  double *log_storage;
  signed char *sign_storage;

  if(workspace == NULL || observation_count == 0 ||
     observation_count > SIZE_MAX / (3U * sizeof(double)) ||
     observation_count > SIZE_MAX / (3U * sizeof(signed char)))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(workspace->capacity >= observation_count &&
     workspace->level_log_absolute != NULL &&
     workspace->regular_log_absolute != NULL &&
     workspace->jump_log_absolute != NULL &&
     workspace->level_sign != NULL && workspace->regular_sign != NULL &&
     workspace->jump_sign != NULL)
    return NP_CONTINUOUS_ROW_OK;

  log_storage = (double *)malloc(
    3U * observation_count * sizeof(double));
  sign_storage = (signed char *)malloc(
    3U * observation_count * sizeof(signed char));
  if(log_storage == NULL || sign_storage == NULL) {
    free(log_storage);
    free(sign_storage);
    return NP_CONTINUOUS_ROW_ERR_MEMORY;
  }

  free(workspace->level_log_absolute);
  free(workspace->level_sign);
  workspace->level_log_absolute = log_storage;
  workspace->regular_log_absolute = log_storage + observation_count;
  workspace->jump_log_absolute = log_storage + 2U * observation_count;
  workspace->level_sign = sign_storage;
  workspace->regular_sign = sign_storage + observation_count;
  workspace->jump_sign = sign_storage + 2U * observation_count;
  workspace->capacity = observation_count;
  return NP_CONTINUOUS_ROW_OK;
}
/*
 * Form one level/derivative observation after the caller has selected the
 * factor and omission topology.  Keeping this arithmetic in one inlined
 * owner lets pure fits avoid categorical-provider and leave-one-out tests in
 * their O(n_train * n_eval) loop without duplicating beta-kernel algebra.
 */
NP_CONTINUOUS_ROW_GRADIENT_ALWAYS_INLINE NPContinuousKernelRowStatus
np_continuous_kernel_beta_level_derivative_observation_bound(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int evaluation_index,
  int observation,
  int derivative_coordinate,
  double other_log,
  int other_sign,
  NPContinuousKernelLevelDerivativeWorkspace *workspace,
  double *common_log_scale,
  NPContinuousKernelDerivativeDiagnostics *diagnostics)
{
  double level_log = -INFINITY;
  int level_sign = 0;
  np_beta_derivative derivative;
  int coordinate;

  for(coordinate = 0; coordinate < plan->num_continuous; ++coordinate) {
    const int local_coordinate =
      coordinate - segment->coordinate_offset;
    const double evaluation = plan->train_is_eval ?
      plan->train[coordinate][evaluation_index] :
      plan->evaluation[coordinate][evaluation_index];
    const double bandwidth = plan->bandwidth_mode == BW_FIXED ?
      plan->bandwidth_eval[coordinate][0] :
      (plan->bandwidth_mode == BW_GEN_NN ?
       plan->bandwidth_eval[coordinate][evaluation_index] :
       plan->bandwidth_train[coordinate][observation]);
    const double observed = plan->train[coordinate][observation];
    np_beta_status beta_status = NP_BETA_OK;

    if(coordinate == derivative_coordinate) {
      beta_status = np_beta_log_abs_pdf_derivative_order(
        evaluation, observed, bandwidth,
        segment->lower[local_coordinate],
        segment->upper[local_coordinate], segment->descriptor.order,
        &level_log, &level_sign, &derivative);
    } else {
      int scalar_sign = 0;
      const double scalar_log = np_beta_log_abs_pdf_order(
        evaluation, observed, bandwidth,
        segment->lower[local_coordinate],
        segment->upper[local_coordinate], segment->descriptor.order,
        &scalar_sign, &beta_status);

      if(beta_status == NP_BETA_OK) {
        if(scalar_sign == 0)
          other_sign = 0;
        else if(other_sign != 0) {
          other_sign *= scalar_sign;
          other_log += scalar_log;
        }
      }
    }
    if(beta_status != NP_BETA_OK) {
      if(diagnostics != NULL) {
        diagnostics->bad_coordinate = coordinate;
        diagnostics->bad_observation = observation;
        diagnostics->beta_status = beta_status;
      }
      return NP_CONTINUOUS_ROW_ERR_KERNEL;
    }
  }

  if(other_sign == 0 || level_sign == 0) {
    workspace->level_log_absolute[observation] = -INFINITY;
    workspace->level_sign[observation] = 0;
  } else {
    workspace->level_log_absolute[observation] = other_log + level_log;
    workspace->level_sign[observation] =
      (signed char)(other_sign * level_sign);
    *common_log_scale = fmax(
      *common_log_scale, workspace->level_log_absolute[observation]);
  }
  if(other_sign == 0 || derivative.regular_sign == 0) {
    workspace->regular_log_absolute[observation] = -INFINITY;
    workspace->regular_sign[observation] = 0;
  } else {
    workspace->regular_log_absolute[observation] =
      other_log + derivative.regular_log_absolute;
    workspace->regular_sign[observation] =
      (signed char)(other_sign * derivative.regular_sign);
    *common_log_scale = fmax(
      *common_log_scale, workspace->regular_log_absolute[observation]);
  }
  if(other_sign == 0 || derivative.jump_sign == 0) {
    workspace->jump_log_absolute[observation] = -INFINITY;
    workspace->jump_sign[observation] = 0;
  } else {
    workspace->jump_log_absolute[observation] =
      other_log + derivative.jump_log_absolute;
    workspace->jump_sign[observation] =
      (signed char)(other_sign * derivative.jump_sign);
    *common_log_scale = fmax(
      *common_log_scale, workspace->jump_log_absolute[observation]);
  }
  return NP_CONTINUOUS_ROW_OK;
}

NPContinuousKernelRowStatus
np_continuous_kernel_beta_level_derivative_log_row_validated(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  int derivative_coordinate,
  const NPContinuousKernelLogFactorProvider *provider,
  NPContinuousKernelLevelDerivativeWorkspace *workspace,
  double *common_log_scale,
  NPContinuousKernelDerivativeDiagnostics *diagnostics)
{
  const NPContinuousKernelSegment *segment;
  NPContinuousKernelRowStatus row_status;
  int observation;
  int coordinate;

  if(diagnostics != NULL) {
    diagnostics->bad_coordinate = -1;
    diagnostics->bad_observation = -1;
    diagnostics->undefined_count = 0;
    diagnostics->beta_status = NP_BETA_OK;
  }
  if(plan == NULL || plan->route == NULL || workspace == NULL ||
     common_log_scale == NULL || plan->num_train <= 0 ||
     plan->num_eval <= 0 || plan->num_continuous <= 0 ||
     plan->route->segment_count != 1 ||
     plan->route->segment[0].descriptor.family != NP_CKERNEL_FAMILY_BETA ||
     evaluation_index < 0 || evaluation_index >= plan->num_eval ||
     omitted_observation < -1 || omitted_observation >= plan->num_train ||
     derivative_coordinate < 0 ||
     derivative_coordinate >= plan->num_continuous ||
     (provider != NULL && provider->function == NULL))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  row_status = np_continuous_kernel_row_plan_validate(plan);
  if(row_status != NP_CONTINUOUS_ROW_OK)
    return row_status;
  for(coordinate = 0; coordinate < plan->num_continuous; ++coordinate)
    if(plan->operator[coordinate] != OP_NORMAL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  row_status = np_continuous_kernel_level_derivative_workspace_reserve(
    workspace, (size_t)plan->num_train);
  if(row_status != NP_CONTINUOUS_ROW_OK)
    return row_status;
  if(provider != NULL) {
    row_status = provider->function(
      provider->context, evaluation_index, omitted_observation,
      plan->num_train, workspace->level_log_absolute,
      workspace->level_sign);
    if(row_status != NP_CONTINUOUS_ROW_OK)
      return row_status;
  }

  segment = &plan->route->segment[0];
  *common_log_scale = -INFINITY;
  if(provider == NULL && omitted_observation == -1) {
    for(observation = 0; observation < plan->num_train; ++observation) {
      row_status =
        np_continuous_kernel_beta_level_derivative_observation_bound(
          plan, segment, evaluation_index, observation,
          derivative_coordinate, 0.0, 1, workspace, common_log_scale,
          diagnostics);
      if(row_status != NP_CONTINUOUS_ROW_OK)
        return row_status;
    }
    return NP_CONTINUOUS_ROW_OK;
  }

  for(observation = 0; observation < plan->num_train; ++observation) {
    double other_log = 0.0;
    int other_sign = 1;

    if(observation == omitted_observation) {
      workspace->level_log_absolute[observation] = -INFINITY;
      workspace->regular_log_absolute[observation] = -INFINITY;
      workspace->jump_log_absolute[observation] = -INFINITY;
      workspace->level_sign[observation] = 0;
      workspace->regular_sign[observation] = 0;
      workspace->jump_sign[observation] = 0;
      continue;
    }
    if(provider != NULL) {
      other_log = workspace->level_log_absolute[observation];
      other_sign = workspace->level_sign[observation];
      if(ISNAN(other_log) || other_log == INFINITY ||
         (other_sign != -1 && other_sign != 0 && other_sign != 1) ||
         ((other_sign == 0) != (other_log == -INFINITY))) {
        if(diagnostics != NULL)
          diagnostics->bad_observation = observation;
        return NP_CONTINUOUS_ROW_ERR_NUMERIC;
      }
    }
    row_status =
      np_continuous_kernel_beta_level_derivative_observation_bound(
        plan, segment, evaluation_index, observation,
        derivative_coordinate, other_log, other_sign, workspace,
        common_log_scale, diagnostics);
    if(row_status != NP_CONTINUOUS_ROW_OK)
      return row_status;
  }
  return NP_CONTINUOUS_ROW_OK;
}
