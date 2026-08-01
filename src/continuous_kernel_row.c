#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include <R_ext/Arith.h>

#include "headers.h"
#include "continuous_kernel_row.h"

static void np_continuous_kernel_row_reset_diagnostics(
  int *bad_coordinate,
  int *bad_observation,
  np_beta_status *beta_status)
{
  if(bad_coordinate != NULL)
    *bad_coordinate = -1;
  if(bad_observation != NULL)
    *bad_observation = -1;
  if(beta_status != NULL)
    *beta_status = NP_BETA_OK;
}

void np_continuous_kernel_row_workspace_init(
  NPContinuousKernelRowWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  workspace->primary_log_absolute = NULL;
  workspace->secondary_log_absolute = NULL;
  workspace->primary_sign = NULL;
  workspace->secondary_sign = NULL;
  workspace->capacity = 0;
  workspace->secondary_capacity = 0;
}

void np_continuous_kernel_row_workspace_release(
  NPContinuousKernelRowWorkspace *workspace)
{
  if(workspace == NULL)
    return;
  free(workspace->primary_log_absolute);
  free(workspace->secondary_log_absolute);
  free(workspace->primary_sign);
  free(workspace->secondary_sign);
  np_continuous_kernel_row_workspace_init(workspace);
}

NPContinuousKernelRowStatus np_continuous_kernel_row_workspace_reserve(
  NPContinuousKernelRowWorkspace *workspace,
  size_t observation_count,
  int need_secondary)
{
  const int need_primary_allocation =
    workspace != NULL && workspace->capacity < observation_count;
  const int need_secondary_allocation =
    workspace != NULL && need_secondary &&
    workspace->secondary_capacity < observation_count;
  double *primary_log_absolute = NULL;
  signed char *primary_sign = NULL;
  double *secondary_log_absolute = NULL;
  signed char *secondary_sign = NULL;

  if(workspace == NULL || observation_count == 0 ||
     observation_count > SIZE_MAX / sizeof(double))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(!need_primary_allocation &&
     workspace->primary_log_absolute != NULL &&
     workspace->primary_sign != NULL &&
     !need_secondary_allocation &&
     (!need_secondary ||
      (workspace->secondary_log_absolute != NULL &&
       workspace->secondary_sign != NULL)))
    return NP_CONTINUOUS_ROW_OK;

  if(need_primary_allocation || workspace->primary_log_absolute == NULL ||
     workspace->primary_sign == NULL) {
    primary_log_absolute = (double *)malloc(
      observation_count * sizeof(double));
    primary_sign = (signed char *)malloc(
      observation_count * sizeof(signed char));
  }
  if(need_secondary_allocation ||
     (need_secondary &&
      (workspace->secondary_log_absolute == NULL ||
       workspace->secondary_sign == NULL))) {
    secondary_log_absolute = (double *)malloc(
      observation_count * sizeof(double));
    secondary_sign = (signed char *)malloc(
      observation_count * sizeof(signed char));
  }
  if(((need_primary_allocation || workspace->primary_log_absolute == NULL ||
       workspace->primary_sign == NULL) &&
      (primary_log_absolute == NULL || primary_sign == NULL)) ||
     ((need_secondary_allocation ||
       (need_secondary &&
        (workspace->secondary_log_absolute == NULL ||
         workspace->secondary_sign == NULL))) &&
      (secondary_log_absolute == NULL || secondary_sign == NULL))) {
    free(primary_log_absolute);
    free(primary_sign);
    free(secondary_log_absolute);
    free(secondary_sign);
    return NP_CONTINUOUS_ROW_ERR_MEMORY;
  }

  if(primary_log_absolute != NULL) {
    free(workspace->primary_log_absolute);
    free(workspace->primary_sign);
    workspace->primary_log_absolute = primary_log_absolute;
    workspace->primary_sign = primary_sign;
    workspace->capacity = observation_count;
  }
  if(secondary_log_absolute != NULL) {
    free(workspace->secondary_log_absolute);
    free(workspace->secondary_sign);
    workspace->secondary_log_absolute = secondary_log_absolute;
    workspace->secondary_sign = secondary_sign;
    workspace->secondary_capacity = observation_count;
  }
  return NP_CONTINUOUS_ROW_OK;
}

static NPContinuousKernelRowStatus np_continuous_kernel_row_plan_status(
  const NPContinuousKernelRowPlan *plan)
{
  int coordinate;

  if(plan == NULL || plan->route == NULL || plan->num_train <= 0 ||
     plan->num_eval <= 0 || plan->num_continuous <= 0 ||
     plan->train == NULL || plan->operator == NULL ||
     (!plan->train_is_eval && plan->evaluation == NULL) ||
     (plan->train_is_eval != 0 && plan->train_is_eval != 1) ||
     (plan->bandwidth_mode != BW_FIXED &&
      plan->bandwidth_mode != BW_GEN_NN &&
      plan->bandwidth_mode != BW_ADAP_NN))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(np_continuous_kernel_route_validate(
       plan->route, plan->num_continuous) != NP_CKERNEL_ROUTE_OK)
    return NP_CONTINUOUS_ROW_ERR_ROUTE;
  for(coordinate = 0; coordinate < plan->num_continuous; ++coordinate)
    if(plan->train[coordinate] == NULL ||
       (!plan->train_is_eval && plan->evaluation[coordinate] == NULL))
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  return NP_CONTINUOUS_ROW_OK;
}

static NPContinuousKernelRowStatus np_continuous_kernel_row_bandwidths(
  const NPContinuousKernelRowPlan *plan,
  int coordinate,
  int evaluation_index,
  int observation_index,
  int convolution,
  double *evaluation_bandwidth,
  double *observation_bandwidth)
{
  if(evaluation_bandwidth == NULL || observation_bandwidth == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;

  if(plan->bandwidth_mode == BW_FIXED) {
    if(plan->bandwidth_eval == NULL ||
       plan->bandwidth_eval[coordinate] == NULL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
    *evaluation_bandwidth = plan->bandwidth_eval[coordinate][0];
    *observation_bandwidth = *evaluation_bandwidth;
  } else if(plan->bandwidth_mode == BW_GEN_NN) {
    if(plan->bandwidth_eval == NULL ||
       plan->bandwidth_eval[coordinate] == NULL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
    *evaluation_bandwidth =
      plan->bandwidth_eval[coordinate][evaluation_index];
    *observation_bandwidth = *evaluation_bandwidth;
    if(convolution && plan->bandwidth_train != NULL &&
       plan->bandwidth_train[coordinate] != NULL)
      *observation_bandwidth =
        plan->bandwidth_train[coordinate][observation_index];
  } else {
    if(plan->bandwidth_train == NULL ||
       plan->bandwidth_train[coordinate] == NULL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
    *observation_bandwidth =
      plan->bandwidth_train[coordinate][observation_index];
    *evaluation_bandwidth = *observation_bandwidth;
    if(convolution) {
      if(plan->bandwidth_eval == NULL ||
         plan->bandwidth_eval[coordinate] == NULL)
        return NP_CONTINUOUS_ROW_ERR_LAYOUT;
      *evaluation_bandwidth =
        plan->bandwidth_eval[coordinate][evaluation_index];
    }
  }
  return NP_CONTINUOUS_ROW_OK;
}

static NPContinuousKernelRowStatus np_continuous_kernel_beta_log_value(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int local_coordinate,
  int coordinate,
  int evaluation_index,
  int observation_index,
  double *log_absolute,
  int *sign,
  np_beta_status *beta_status)
{
  const int operator_code = plan->operator[coordinate];
  const double evaluation = plan->train_is_eval ?
    plan->train[coordinate][evaluation_index] :
    plan->evaluation[coordinate][evaluation_index];
  const double observation =
    plan->train[coordinate][observation_index];
  double evaluation_bandwidth;
  double observation_bandwidth;
  NPContinuousKernelRowStatus row_status;

  row_status = np_continuous_kernel_row_bandwidths(
    plan, coordinate, evaluation_index, observation_index,
    operator_code == OP_CONVOLUTION,
    &evaluation_bandwidth, &observation_bandwidth);
  if(row_status != NP_CONTINUOUS_ROW_OK)
    return row_status;

  if(operator_code == OP_NORMAL) {
    *log_absolute = np_beta_log_abs_pdf_order(
      evaluation, observation, evaluation_bandwidth,
      segment->lower[local_coordinate], segment->upper[local_coordinate],
      segment->descriptor.order, sign, beta_status);
  } else if(operator_code == OP_INTEGRAL) {
    *log_absolute = np_beta_log_abs_cdf_order(
      evaluation, observation, evaluation_bandwidth,
      segment->lower[local_coordinate], segment->upper[local_coordinate],
      segment->descriptor.order, sign, beta_status);
  } else if(operator_code == OP_CONVOLUTION) {
    *log_absolute = np_beta_log_abs_overlap_order(
      evaluation, evaluation_bandwidth, observation, observation_bandwidth,
      segment->lower[local_coordinate], segment->upper[local_coordinate],
      segment->descriptor.order, sign, beta_status);
  } else {
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  }

  if(*beta_status != NP_BETA_OK)
    return NP_CONTINUOUS_ROW_ERR_KERNEL;
  if(ISNAN(*log_absolute) || *log_absolute == INFINITY ||
     (*sign != -1 && *sign != 0 && *sign != 1) ||
     ((*sign == 0) != (*log_absolute == -INFINITY)))
    return NP_CONTINUOUS_ROW_ERR_NUMERIC;
  return NP_CONTINUOUS_ROW_OK;
}

static double np_continuous_kernel_row_maximum(
  const double *log_absolute,
  const signed char *sign,
  int observation_count)
{
  double maximum = -INFINITY;
  int observation;

  for(observation = 0; observation < observation_count; ++observation)
    if(sign[observation] != 0 && log_absolute[observation] > maximum)
      maximum = log_absolute[observation];
  return maximum;
}

static NPContinuousKernelRowStatus np_continuous_kernel_row_multiply_scaled(
  double *row,
  const double *log_absolute,
  const signed char *sign,
  int observation_count,
  double log_scale)
{
  int observation;

  if(row == NULL || ISNAN(log_scale) || log_scale == INFINITY)
    return NP_CONTINUOUS_ROW_ERR_NUMERIC;
  if(log_scale == -INFINITY) {
    for(observation = 0; observation < observation_count; ++observation)
      row[observation] = 0.0;
    return NP_CONTINUOUS_ROW_OK;
  }

  for(observation = 0; observation < observation_count; ++observation) {
    const double factor = (sign[observation] == 0) ? 0.0 :
      (double)sign[observation] *
      exp(log_absolute[observation] - log_scale);

    if(!R_FINITE(factor))
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
    row[observation] *= factor;
    if(!R_FINITE(row[observation]))
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
  }
  return NP_CONTINUOUS_ROW_OK;
}

static NPContinuousKernelRowStatus np_continuous_kernel_row_scale_sum(
  const double *segment_log_scale,
  int segment_count,
  double *total_log_scale)
{
  double total = 0.0;
  int segment;

  if(segment_log_scale == NULL || total_log_scale == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(segment = 0; segment < segment_count; ++segment) {
    if(ISNAN(segment_log_scale[segment]) ||
       segment_log_scale[segment] == INFINITY)
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
    if(segment_log_scale[segment] == -INFINITY) {
      *total_log_scale = -INFINITY;
      return NP_CONTINUOUS_ROW_OK;
    }
    total += segment_log_scale[segment];
    if(!R_FINITE(total))
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
  }
  *total_log_scale = total;
  return NP_CONTINUOUS_ROW_OK;
}

NPContinuousKernelRowStatus np_continuous_kernel_beta_factor_row(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *result)
{
  NPContinuousKernelRowStatus status;
  int segment_index;
  int observation;

  if(result == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  np_continuous_kernel_row_reset_diagnostics(
    &result->bad_coordinate, &result->bad_observation, &result->beta_status);
  result->total_log_scale = 0.0;
  for(segment_index = 0; segment_index < NP_CKERNEL_ROUTE_MAX_SEGMENTS;
      ++segment_index)
    result->segment_log_scale[segment_index] = 0.0;

  status = np_continuous_kernel_row_plan_status(plan);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  if(result->row == NULL ||
     evaluation_index < 0 || evaluation_index >= plan->num_eval ||
     omitted_observation < -1 || omitted_observation >= plan->num_train)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  status = np_continuous_kernel_row_workspace_reserve(
    workspace, (size_t)plan->num_train, 0);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;

  for(observation = 0; observation < plan->num_train; ++observation)
    result->row[observation] =
      (observation == omitted_observation) ? 0.0 : 1.0;

  for(segment_index = 0; segment_index < plan->route->segment_count;
      ++segment_index) {
    const NPContinuousKernelSegment * const segment =
      &plan->route->segment[segment_index];
    int coordinate;

    if(segment->descriptor.family != NP_CKERNEL_FAMILY_BETA)
      continue;
    for(observation = 0; observation < plan->num_train; ++observation) {
      double log_product = 0.0;
      int product_sign = 1;

      if(observation == omitted_observation) {
        workspace->primary_log_absolute[observation] = -INFINITY;
        workspace->primary_sign[observation] = 0;
        continue;
      }
      for(coordinate = segment->coordinate_offset;
          coordinate < segment->coordinate_offset + segment->coordinate_count;
          ++coordinate) {
        const int local_coordinate = coordinate - segment->coordinate_offset;
        double scalar_log = -INFINITY;
        int scalar_sign = 0;

        status = np_continuous_kernel_beta_log_value(
          plan, segment, local_coordinate, coordinate,
          evaluation_index, observation,
          &scalar_log, &scalar_sign, &result->beta_status);
        if(status != NP_CONTINUOUS_ROW_OK) {
          result->bad_coordinate = coordinate;
          result->bad_observation = observation;
          return status;
        }
        if(scalar_sign == 0) {
          product_sign = 0;
          log_product = -INFINITY;
          break;
        }
        log_product += scalar_log;
        product_sign *= scalar_sign;
        if(ISNAN(log_product) || log_product == INFINITY) {
          result->bad_coordinate = coordinate;
          result->bad_observation = observation;
          return NP_CONTINUOUS_ROW_ERR_NUMERIC;
        }
      }
      workspace->primary_log_absolute[observation] = log_product;
      workspace->primary_sign[observation] = (signed char)product_sign;
    }

    result->segment_log_scale[segment_index] =
      np_continuous_kernel_row_maximum(
        workspace->primary_log_absolute, workspace->primary_sign,
        plan->num_train);
    status = np_continuous_kernel_row_multiply_scaled(
      result->row, workspace->primary_log_absolute,
      workspace->primary_sign, plan->num_train,
      result->segment_log_scale[segment_index]);
    if(status != NP_CONTINUOUS_ROW_OK)
      return status;
  }

  return np_continuous_kernel_row_scale_sum(
    result->segment_log_scale, plan->route->segment_count,
    &result->total_log_scale);
}

NPContinuousKernelRowStatus np_continuous_kernel_beta_derivative_factor_row(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  int derivative_coordinate,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelDerivativeRowResult *result)
{
  NPContinuousKernelRowStatus status;
  int segment_index;
  int observation;

  if(result == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  np_continuous_kernel_row_reset_diagnostics(
    &result->bad_coordinate, &result->bad_observation, &result->beta_status);
  result->regular_total_log_scale = 0.0;
  result->jump_total_log_scale = 0.0;
  for(segment_index = 0; segment_index < NP_CKERNEL_ROUTE_MAX_SEGMENTS;
      ++segment_index) {
    result->regular_segment_log_scale[segment_index] = 0.0;
    result->jump_segment_log_scale[segment_index] = 0.0;
  }

  status = np_continuous_kernel_row_plan_status(plan);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  if(result->regular_row == NULL || result->jump_row == NULL ||
     evaluation_index < 0 || evaluation_index >= plan->num_eval ||
     derivative_coordinate < 0 ||
     derivative_coordinate >= plan->num_continuous ||
     omitted_observation < -1 || omitted_observation >= plan->num_train)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  status = np_continuous_kernel_row_workspace_reserve(
    workspace, (size_t)plan->num_train, 1);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;

  for(observation = 0; observation < plan->num_train; ++observation) {
    const double initial =
      (observation == omitted_observation) ? 0.0 : 1.0;
    result->regular_row[observation] = initial;
    result->jump_row[observation] = initial;
  }

  for(segment_index = 0; segment_index < plan->route->segment_count;
      ++segment_index) {
    const NPContinuousKernelSegment * const segment =
      &plan->route->segment[segment_index];
    const int segment_has_derivative =
      segment->descriptor.family == NP_CKERNEL_FAMILY_BETA &&
      derivative_coordinate >= segment->coordinate_offset &&
      derivative_coordinate <
        segment->coordinate_offset + segment->coordinate_count;
    int coordinate;

    if(segment->descriptor.family != NP_CKERNEL_FAMILY_BETA)
      continue;
    for(observation = 0; observation < plan->num_train; ++observation) {
      double regular_log = 0.0;
      double jump_log = 0.0;
      int regular_sign = 1;
      int jump_sign = 1;

      if(observation == omitted_observation) {
        workspace->primary_log_absolute[observation] = -INFINITY;
        workspace->secondary_log_absolute[observation] = -INFINITY;
        workspace->primary_sign[observation] = 0;
        workspace->secondary_sign[observation] = 0;
        continue;
      }

      for(coordinate = segment->coordinate_offset;
          coordinate < segment->coordinate_offset + segment->coordinate_count;
          ++coordinate) {
        const int local_coordinate = coordinate - segment->coordinate_offset;

        if(coordinate == derivative_coordinate) {
          const double evaluation = plan->train_is_eval ?
            plan->train[coordinate][evaluation_index] :
            plan->evaluation[coordinate][evaluation_index];
          const double observed = plan->train[coordinate][observation];
          double evaluation_bandwidth;
          double observation_bandwidth;
          np_beta_derivative derivative;

          if(plan->operator[coordinate] != OP_DERIVATIVE)
            return NP_CONTINUOUS_ROW_ERR_LAYOUT;
          status = np_continuous_kernel_row_bandwidths(
            plan, coordinate, evaluation_index, observation, 0,
            &evaluation_bandwidth, &observation_bandwidth);
          if(status != NP_CONTINUOUS_ROW_OK)
            return status;
          (void)observation_bandwidth;
          result->beta_status = np_beta_pdf_derivative_order(
            evaluation, observed, evaluation_bandwidth,
            segment->lower[local_coordinate],
            segment->upper[local_coordinate],
            segment->descriptor.order, &derivative);
          if(result->beta_status != NP_BETA_OK) {
            result->bad_coordinate = coordinate;
            result->bad_observation = observation;
            return NP_CONTINUOUS_ROW_ERR_KERNEL;
          }
          regular_log += derivative.regular_log_absolute;
          regular_sign *= derivative.regular_sign;
          jump_log += derivative.jump_log_absolute;
          jump_sign *= derivative.jump_sign;
          if(derivative.regular_sign == 0) {
            regular_log = -INFINITY;
            regular_sign = 0;
          }
          if(derivative.jump_sign == 0) {
            jump_log = -INFINITY;
            jump_sign = 0;
          }
        } else {
          double scalar_log = -INFINITY;
          int scalar_sign = 0;

          status = np_continuous_kernel_beta_log_value(
            plan, segment, local_coordinate, coordinate,
            evaluation_index, observation,
            &scalar_log, &scalar_sign, &result->beta_status);
          if(status != NP_CONTINUOUS_ROW_OK) {
            result->bad_coordinate = coordinate;
            result->bad_observation = observation;
            return status;
          }
          if(scalar_sign == 0) {
            regular_log = -INFINITY;
            jump_log = -INFINITY;
            regular_sign = 0;
            jump_sign = 0;
            break;
          }
          if(regular_sign != 0) {
            regular_log += scalar_log;
            regular_sign *= scalar_sign;
          }
          if(jump_sign != 0) {
            jump_log += scalar_log;
            jump_sign *= scalar_sign;
          }
        }
      }

      workspace->primary_log_absolute[observation] = regular_log;
      workspace->secondary_log_absolute[observation] =
        segment_has_derivative ? jump_log : regular_log;
      workspace->primary_sign[observation] = (signed char)regular_sign;
      workspace->secondary_sign[observation] = (signed char)
        (segment_has_derivative ? jump_sign : regular_sign);
    }

    result->regular_segment_log_scale[segment_index] =
      np_continuous_kernel_row_maximum(
        workspace->primary_log_absolute, workspace->primary_sign,
        plan->num_train);
    result->jump_segment_log_scale[segment_index] =
      np_continuous_kernel_row_maximum(
        workspace->secondary_log_absolute, workspace->secondary_sign,
        plan->num_train);
    status = np_continuous_kernel_row_multiply_scaled(
      result->regular_row, workspace->primary_log_absolute,
      workspace->primary_sign, plan->num_train,
      result->regular_segment_log_scale[segment_index]);
    if(status == NP_CONTINUOUS_ROW_OK)
      status = np_continuous_kernel_row_multiply_scaled(
        result->jump_row, workspace->secondary_log_absolute,
        workspace->secondary_sign, plan->num_train,
        result->jump_segment_log_scale[segment_index]);
    if(status != NP_CONTINUOUS_ROW_OK)
      return status;
  }

  status = np_continuous_kernel_row_scale_sum(
    result->regular_segment_log_scale, plan->route->segment_count,
    &result->regular_total_log_scale);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  return np_continuous_kernel_row_scale_sum(
    result->jump_segment_log_scale, plan->route->segment_count,
    &result->jump_total_log_scale);
}

NPContinuousKernelRowStatus np_continuous_kernel_scaled_restore(
  double scaled_value,
  double log_scale,
  int power,
  double *value)
{
  double log_absolute;

  if(value == NULL || power <= 0 || !R_FINITE(scaled_value) ||
     ISNAN(log_scale) || log_scale == INFINITY)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(scaled_value == 0.0 || log_scale == -INFINITY) {
    *value = 0.0;
    return NP_CONTINUOUS_ROW_OK;
  }
  if(log_scale == 0.0) {
    *value = scaled_value;
    return NP_CONTINUOUS_ROW_OK;
  }
  log_absolute = log(fabs(scaled_value)) + (double)power * log_scale;
  if(ISNAN(log_absolute) || log_absolute > log(DBL_MAX))
    return NP_CONTINUOUS_ROW_ERR_NUMERIC;
  *value = copysign(exp(log_absolute), scaled_value);
  return R_FINITE(*value) ? NP_CONTINUOUS_ROW_OK :
    NP_CONTINUOUS_ROW_ERR_NUMERIC;
}

NPContinuousKernelRowStatus np_continuous_kernel_signed_log_restore(
  double log_absolute,
  int sign,
  double *value)
{
  if(value == NULL || ISNAN(log_absolute) || log_absolute == INFINITY ||
     (sign != -1 && sign != 0 && sign != 1) ||
     ((sign == 0) != (log_absolute == -INFINITY)))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(sign == 0) {
    *value = 0.0;
    return NP_CONTINUOUS_ROW_OK;
  }
  if(log_absolute > log(DBL_MAX))
    return NP_CONTINUOUS_ROW_ERR_NUMERIC;
  *value = (double)sign * exp(log_absolute);
  return R_FINITE(*value) ? NP_CONTINUOUS_ROW_OK :
    NP_CONTINUOUS_ROW_ERR_NUMERIC;
}

const char *np_continuous_kernel_row_status_message(
  NPContinuousKernelRowStatus status)
{
  switch(status) {
  case NP_CONTINUOUS_ROW_OK:
    return "success";
  case NP_CONTINUOUS_ROW_ERR_LAYOUT:
    return "continuous-kernel row has an invalid layout";
  case NP_CONTINUOUS_ROW_ERR_ROUTE:
    return "continuous-kernel row has invalid route metadata";
  case NP_CONTINUOUS_ROW_ERR_MEMORY:
    return "continuous-kernel row could not allocate its workspace";
  case NP_CONTINUOUS_ROW_ERR_KERNEL:
    return "continuous-kernel scalar evaluation failed";
  case NP_CONTINUOUS_ROW_ERR_NUMERIC:
    return "continuous-kernel row produced an invalid numeric result";
  default:
    return "unknown continuous-kernel row status";
  }
}
