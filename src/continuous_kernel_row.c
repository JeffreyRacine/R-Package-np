#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include <R_ext/Arith.h>
#include <R_ext/Memory.h>
#include <R_ext/Utils.h>

#include "headers.h"
#include "continuous_kernel_row.h"

#if defined(__clang__) || defined(__GNUC__)
# define NP_CONTINUOUS_ROW_ALWAYS_INLINE \
  static inline __attribute__((always_inline))
# define NP_CONTINUOUS_ROW_NOINLINE __attribute__((noinline))
#else
# define NP_CONTINUOUS_ROW_ALWAYS_INLINE static inline
# define NP_CONTINUOUS_ROW_NOINLINE
#endif

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

void np_continuous_kernel_beta_prepared_context_init(
  NPContinuousKernelBetaPreparedContext *context)
{
  if(context == NULL)
    return;
  context->pdf_active = 0;
  context->pdf_row_component_active = 0;
  context->cdf_active = 0;
  context->cdf_concentration_active = 0;
  context->allocation_active = 0;
  context->num_train = 0;
  context->num_continuous = 0;
  context->num_beta_coordinates = 0;
  context->num_cdf_coordinates = 0;
  context->allocation_marker = NULL;
  context->coordinate_slot = NULL;
  context->cdf_coordinate_slot = NULL;
  context->pdf_observation = NULL;
  context->pdf_observation_status = NULL;
  context->pdf_row_component = NULL;
  context->pdf_row_derivative_component = NULL;
  context->pdf_log_abs_coefficient = NULL;
  context->pdf_coefficient_sign = NULL;
  context->pdf_first_interior = NULL;
  context->pdf_second_interior = NULL;
  context->cdf_observation = NULL;
  context->cdf_observation_status = NULL;
  context->cdf_log_abs_coefficient = NULL;
  context->cdf_coefficient_sign = NULL;
  context->cdf_concentration = NULL;
  context->cdf_concentration_status = NULL;
}

void np_continuous_kernel_beta_prepared_context_release(
  NPContinuousKernelBetaPreparedContext *context)
{
  if(context == NULL)
    return;
  if(context->allocation_active)
    vmaxset(context->allocation_marker);
  np_continuous_kernel_beta_prepared_context_init(context);
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

NPContinuousKernelRowStatus np_continuous_kernel_row_plan_validate(
  const NPContinuousKernelRowPlan *plan)
{
  int coordinate;

  if(plan == NULL || plan->route == NULL || plan->num_train <= 0 ||
     plan->num_eval <= 0 || plan->num_continuous <= 0 ||
     plan->train == NULL || plan->operator == NULL ||
     (!plan->train_is_eval && plan->evaluation == NULL) ||
     (plan->train_is_eval != 0 && plan->train_is_eval != 1) ||
     (plan->train_is_eval && plan->num_train != plan->num_eval) ||
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

static NPContinuousKernelRowStatus
np_continuous_kernel_beta_prepared_cdf_context_prepare(
  NPContinuousKernelBetaPreparedContext *context,
  const NPContinuousKernelRowPlan *plan)
{
  size_t observation_count;
  size_t component_count;
  int cdf_coordinate_count = 0;
  int concentration_eligible = plan->bandwidth_mode == BW_FIXED;
  int coordinate;

  for(coordinate = 0; coordinate < plan->num_continuous; ++coordinate) {
    const NPContinuousKernelSegment * const segment =
      np_continuous_kernel_route_segment(plan->route, coordinate);

    if(segment == NULL)
      return NP_CONTINUOUS_ROW_ERR_ROUTE;
    if(segment->descriptor.family == NP_CKERNEL_FAMILY_BETA &&
       plan->operator[coordinate] == OP_INTEGRAL) {
      ++cdf_coordinate_count;
      if(plan->bandwidth_eval == NULL ||
         plan->bandwidth_eval[coordinate] == NULL)
        concentration_eligible = 0;
    }
  }
  if(cdf_coordinate_count == 0)
    return NP_CONTINUOUS_ROW_OK;
  if((size_t)cdf_coordinate_count >
       SIZE_MAX / (size_t)plan->num_train ||
     (size_t)cdf_coordinate_count >
       SIZE_MAX / NP_BETA_PREPARED_MAX_COMPONENTS)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  observation_count = (size_t)cdf_coordinate_count *
    (size_t)plan->num_train;
  component_count = (size_t)cdf_coordinate_count *
    NP_BETA_PREPARED_MAX_COMPONENTS;

  context->allocation_marker = vmaxget();
  context->allocation_active = 1;
  context->cdf_coordinate_slot = (int *)R_alloc(
    (size_t)plan->num_continuous, (int)sizeof(int));
  context->cdf_observation = (np_beta_cdf_observation *)R_alloc(
    observation_count, (int)sizeof(np_beta_cdf_observation));
  context->cdf_observation_status = (np_beta_status *)R_alloc(
    observation_count, (int)sizeof(np_beta_status));
  context->cdf_log_abs_coefficient = (double *)R_alloc(
    component_count, (int)sizeof(double));
  context->cdf_coefficient_sign = (signed char *)R_alloc(
    component_count, (int)sizeof(signed char));
  if(concentration_eligible) {
    context->cdf_concentration = (double *)R_alloc(
      component_count, (int)sizeof(double));
    context->cdf_concentration_status = (np_beta_status *)R_alloc(
      component_count, (int)sizeof(np_beta_status));
  }

  cdf_coordinate_count = 0;
  for(coordinate = 0; coordinate < plan->num_continuous; ++coordinate) {
    const NPContinuousKernelSegment * const segment =
      np_continuous_kernel_route_segment(plan->route, coordinate);
    const int local_coordinate = segment == NULL ? -1 :
      coordinate - segment->coordinate_offset;
    int observation;

    context->cdf_coordinate_slot[coordinate] = -1;
    if(segment == NULL ||
       segment->descriptor.family != NP_CKERNEL_FAMILY_BETA ||
       plan->operator[coordinate] != OP_INTEGRAL)
      continue;
    context->cdf_coordinate_slot[coordinate] = cdf_coordinate_count;
    for(int component = 0;
        component < NP_BETA_PREPARED_MAX_COMPONENTS;
        ++component) {
      const size_t component_offset =
        (size_t)cdf_coordinate_count * NP_BETA_PREPARED_MAX_COMPONENTS +
        (size_t)component;
      np_beta_pdf_component coefficient_component = {0};

      context->cdf_log_abs_coefficient[component_offset] = 0.0;
      context->cdf_coefficient_sign[component_offset] = 0;
      if(concentration_eligible) {
        context->cdf_concentration[component_offset] = 0.0;
        context->cdf_concentration_status[component_offset] =
          NP_BETA_ERR_SCALE;
      }
      if(component < segment->descriptor.order / 2) {
        const np_beta_status coefficient_status =
          np_beta_pdf_component_prepare_coefficient(
            &coefficient_component, segment->descriptor.order, component);

        if(coefficient_status != NP_BETA_OK)
          return NP_CONTINUOUS_ROW_ERR_ROUTE;
        context->cdf_log_abs_coefficient[component_offset] =
          coefficient_component.log_abs_coefficient;
        context->cdf_coefficient_sign[component_offset] =
          (signed char)coefficient_component.coefficient_sign;
        if(concentration_eligible)
          context->cdf_concentration_status[component_offset] =
            np_beta_concentration_prepare(
              plan->bandwidth_eval[coordinate][0],
              segment->upper[local_coordinate] -
                segment->lower[local_coordinate],
              component + 1,
              &context->cdf_concentration[component_offset]);
      }
    }
    for(observation = 0; observation < plan->num_train; ++observation) {
      const size_t observation_offset =
        (size_t)cdf_coordinate_count * (size_t)plan->num_train +
        (size_t)observation;

      context->cdf_observation_status[observation_offset] =
        np_beta_cdf_observation_init(
          plan->train[coordinate][observation],
          segment->lower[local_coordinate],
          segment->upper[local_coordinate],
          segment->upper[local_coordinate] -
            segment->lower[local_coordinate],
          &context->cdf_observation[observation_offset]);
    }
    ++cdf_coordinate_count;
  }
  context->num_train = plan->num_train;
  context->num_continuous = plan->num_continuous;
  context->num_cdf_coordinates = cdf_coordinate_count;
  context->cdf_active = 1;
  context->cdf_concentration_active = concentration_eligible;
  return NP_CONTINUOUS_ROW_OK;
}

NPContinuousKernelRowStatus np_continuous_kernel_beta_prepared_context_prepare(
  NPContinuousKernelBetaPreparedContext *context,
  const NPContinuousKernelRowPlan *plan)
{
  size_t observation_count;
  size_t component_count;
  int beta_coordinate_count = 0;
  int coordinate;

  if(context == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(context->allocation_active)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  np_continuous_kernel_beta_prepared_context_init(context);
  if(np_continuous_kernel_row_plan_validate(plan) !=
     NP_CONTINUOUS_ROW_OK)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(coordinate = 0; coordinate < plan->num_continuous; ++coordinate) {
    const NPContinuousKernelSegment * const segment =
      np_continuous_kernel_route_segment(plan->route, coordinate);

    if(segment == NULL)
      return NP_CONTINUOUS_ROW_ERR_ROUTE;
    if(segment->descriptor.family == NP_CKERNEL_FAMILY_BETA &&
       plan->operator[coordinate] == OP_INTEGRAL)
      return np_continuous_kernel_beta_prepared_cdf_context_prepare(
        context, plan);
  }
  for(coordinate = 0; coordinate < plan->num_continuous; ++coordinate) {
    const NPContinuousKernelSegment * const segment =
      np_continuous_kernel_route_segment(plan->route, coordinate);

    if(segment == NULL)
      return NP_CONTINUOUS_ROW_ERR_ROUTE;
    if(segment->descriptor.family != NP_CKERNEL_FAMILY_BETA)
      continue;
    ++beta_coordinate_count;
    if(plan->operator[coordinate] != OP_NORMAL)
      return NP_CONTINUOUS_ROW_OK;
    if((plan->bandwidth_mode == BW_FIXED ||
        plan->bandwidth_mode == BW_GEN_NN) &&
       (plan->bandwidth_eval == NULL ||
        plan->bandwidth_eval[coordinate] == NULL))
      return NP_CONTINUOUS_ROW_OK;
    if(plan->bandwidth_mode == BW_ADAP_NN &&
       (plan->bandwidth_train == NULL ||
        plan->bandwidth_train[coordinate] == NULL))
      return NP_CONTINUOUS_ROW_OK;
  }
  if(beta_coordinate_count == 0)
    return NP_CONTINUOUS_ROW_OK;
  if((size_t)beta_coordinate_count >
       SIZE_MAX / (size_t)plan->num_train ||
     (size_t)beta_coordinate_count >
       SIZE_MAX / NP_BETA_PREPARED_MAX_COMPONENTS)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  observation_count = (size_t)beta_coordinate_count *
    (size_t)plan->num_train;
  component_count = (size_t)beta_coordinate_count *
    NP_BETA_PREPARED_MAX_COMPONENTS;
  context->allocation_marker = vmaxget();
  context->allocation_active = 1;
  context->coordinate_slot = (int *)R_alloc(
    (size_t)plan->num_continuous, (int)sizeof(int));
  context->pdf_observation = (np_beta_pdf_observation *)R_alloc(
    observation_count, (int)sizeof(np_beta_pdf_observation));
  context->pdf_observation_status = (np_beta_status *)R_alloc(
    observation_count, (int)sizeof(np_beta_status));
  if(plan->bandwidth_mode == BW_ADAP_NN) {
    context->pdf_log_abs_coefficient = (double *)R_alloc(
      component_count, (int)sizeof(double));
    context->pdf_coefficient_sign = (signed char *)R_alloc(
      component_count, (int)sizeof(signed char));
  } else {
    context->pdf_row_component = (np_beta_pdf_component *)R_alloc(
      component_count, (int)sizeof(np_beta_pdf_component));
  }
  context->pdf_first_interior = (int *)R_alloc(
    (size_t)beta_coordinate_count, (int)sizeof(int));
  context->pdf_second_interior = (int *)R_alloc(
    (size_t)beta_coordinate_count, (int)sizeof(int));

  beta_coordinate_count = 0;
  for(coordinate = 0; coordinate < plan->num_continuous; ++coordinate) {
    const NPContinuousKernelSegment * const segment =
      np_continuous_kernel_route_segment(plan->route, coordinate);
    const int local_coordinate = segment == NULL ? -1 :
      coordinate - segment->coordinate_offset;
    int observation;

    context->coordinate_slot[coordinate] = -1;
    if(segment == NULL ||
       segment->descriptor.family != NP_CKERNEL_FAMILY_BETA)
      continue;
    context->coordinate_slot[coordinate] = beta_coordinate_count;
    context->pdf_first_interior[beta_coordinate_count] = -1;
    context->pdf_second_interior[beta_coordinate_count] = -1;
    if(plan->bandwidth_mode == BW_ADAP_NN) {
      for(int component = 0;
          component < NP_BETA_PREPARED_MAX_COMPONENTS;
          ++component) {
        const size_t component_offset =
          (size_t)beta_coordinate_count *
            NP_BETA_PREPARED_MAX_COMPONENTS +
          (size_t)component;
        np_beta_pdf_component coefficient_component = {0};

        context->pdf_log_abs_coefficient[component_offset] = 0.0;
        context->pdf_coefficient_sign[component_offset] = 0;
        if(component < segment->descriptor.order / 2) {
          const np_beta_status coefficient_status =
            np_beta_pdf_component_prepare_coefficient(
              &coefficient_component, segment->descriptor.order, component);

          if(coefficient_status != NP_BETA_OK)
            return NP_CONTINUOUS_ROW_ERR_ROUTE;
          context->pdf_log_abs_coefficient[component_offset] =
            coefficient_component.log_abs_coefficient;
          context->pdf_coefficient_sign[component_offset] =
            (signed char)coefficient_component.coefficient_sign;
        }
      }
    }
    for(observation = 0; observation < plan->num_train; ++observation) {
      const double observed = plan->train[coordinate][observation];
      const double support_length =
        segment->upper[local_coordinate] -
        segment->lower[local_coordinate];
      const size_t observation_offset =
        (size_t)beta_coordinate_count * (size_t)plan->num_train +
        (size_t)observation;
      np_beta_pdf_observation * const pdf_observation =
        &context->pdf_observation[observation_offset];
      double observation_unit;
      double observation_complement_unit;
      const np_beta_status beta_status =
        np_beta_pdf_observation_init(
          observed,
          segment->lower[local_coordinate],
          segment->upper[local_coordinate], support_length,
          &observation_unit, &observation_complement_unit,
          &pdf_observation->log_unit,
          &pdf_observation->log_complement_unit,
          &pdf_observation->endpoint);

      context->pdf_observation_status[observation_offset] = beta_status;
      if(beta_status != NP_BETA_OK)
        continue;
      if(pdf_observation->endpoint == 0) {
        if(context->pdf_first_interior[beta_coordinate_count] < 0)
          context->pdf_first_interior[beta_coordinate_count] = observation;
        else if(context->pdf_second_interior[beta_coordinate_count] < 0)
          context->pdf_second_interior[beta_coordinate_count] = observation;
      }
    }
    ++beta_coordinate_count;
  }
  context->num_train = plan->num_train;
  context->num_continuous = plan->num_continuous;
  context->num_beta_coordinates = beta_coordinate_count;
  context->pdf_active = 1;
  context->pdf_row_component_active =
    plan->bandwidth_mode == BW_FIXED ||
    plan->bandwidth_mode == BW_GEN_NN;
  return NP_CONTINUOUS_ROW_OK;
}

NPContinuousKernelRowStatus
np_continuous_kernel_beta_prepared_derivative_context_prepare(
  NPContinuousKernelBetaPreparedContext *context)
{
  size_t component_count;

  if(context == NULL || !context->allocation_active ||
     !context->pdf_active || !context->pdf_row_component_active ||
     context->num_beta_coordinates <= 0 ||
     context->pdf_row_component == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(context->pdf_row_derivative_component != NULL)
    return NP_CONTINUOUS_ROW_OK;
  if((size_t)context->num_beta_coordinates >
       SIZE_MAX / NP_BETA_PREPARED_MAX_COMPONENTS)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  component_count = (size_t)context->num_beta_coordinates *
    NP_BETA_PREPARED_MAX_COMPONENTS;
  context->pdf_row_derivative_component =
    (np_beta_pdf_derivative_component *)R_alloc(
      component_count,
      (int)sizeof(np_beta_pdf_derivative_component));
  return NP_CONTINUOUS_ROW_OK;
}

NPContinuousKernelRowStatus
np_continuous_kernel_beta_prepared_pdf_row_prepare(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int evaluation_index,
  int omitted_observation,
  NPContinuousKernelRowResult *result)
{
  NPContinuousKernelBetaPreparedContext * const context =
    plan == NULL ? NULL : plan->beta_prepared;
  const int first_observation = omitted_observation == 0 ? 1 : 0;
  int coordinate;

  if(context == NULL || !context->pdf_active || segment == NULL ||
     result == NULL ||
     context->num_train != plan->num_train ||
     context->num_continuous != plan->num_continuous)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(first_observation >= plan->num_train)
    return NP_CONTINUOUS_ROW_OK;
  for(coordinate = segment->coordinate_offset;
      coordinate < segment->coordinate_offset + segment->coordinate_count;
      ++coordinate) {
    const int local_coordinate = coordinate - segment->coordinate_offset;
    const double evaluation = plan->train_is_eval ?
      plan->train[coordinate][evaluation_index] :
      plan->evaluation[coordinate][evaluation_index];
    const double bandwidth = plan->bandwidth_mode == BW_FIXED ?
      plan->bandwidth_eval[coordinate][0] :
      plan->bandwidth_eval[coordinate][evaluation_index];
    const int beta_slot = context->coordinate_slot[coordinate];
    const int component_count = segment->descriptor.order / 2;
    int normalizer_observation =
      beta_slot < 0 ? -1 : context->pdf_first_interior[beta_slot];
    int component;

    if(beta_slot < 0 || beta_slot >= context->num_beta_coordinates)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
    if(normalizer_observation == omitted_observation)
      normalizer_observation = context->pdf_second_interior[beta_slot];

    if(component_count <= 0 ||
       component_count > NP_BETA_PREPARED_MAX_COMPONENTS)
      return NP_CONTINUOUS_ROW_ERR_ROUTE;
    for(component = 0; component < component_count; ++component) {
      np_beta_pdf_component * const row_component =
        &context->pdf_row_component[
          (size_t)beta_slot * NP_BETA_PREPARED_MAX_COMPONENTS +
          (size_t)component];
      np_beta_shape shape;
      np_beta_status beta_status = np_beta_shape_init(
        evaluation, segment->lower[local_coordinate], bandwidth,
        segment->lower[local_coordinate],
        segment->upper[local_coordinate], component + 1, &shape);

      if(beta_status != NP_BETA_OK) {
        result->beta_status = beta_status;
        result->bad_coordinate = coordinate;
        result->bad_observation = first_observation;
        return NP_CONTINUOUS_ROW_ERR_KERNEL;
      }
      beta_status = np_beta_pdf_component_from_shape(
        &shape, row_component);
      if(beta_status != NP_BETA_OK) {
        result->beta_status = beta_status;
        result->bad_coordinate = coordinate;
        result->bad_observation = first_observation;
        return NP_CONTINUOUS_ROW_ERR_KERNEL;
      }
      beta_status = np_beta_pdf_component_prepare_coefficient(
        row_component, segment->descriptor.order, component);
      if(beta_status != NP_BETA_OK) {
        result->beta_status = beta_status;
        result->bad_coordinate = coordinate;
        result->bad_observation = first_observation;
        return NP_CONTINUOUS_ROW_ERR_KERNEL;
      }
      if(normalizer_observation >= 0) {
        beta_status = np_beta_pdf_component_prepare_normalizer(
          row_component);
        if(beta_status != NP_BETA_OK) {
          result->beta_status = beta_status;
          result->bad_coordinate = coordinate;
          result->bad_observation = normalizer_observation;
          return NP_CONTINUOUS_ROW_ERR_KERNEL;
        }
      }
    }
  }
  return NP_CONTINUOUS_ROW_OK;
}

NPContinuousKernelRowStatus
np_continuous_kernel_beta_prepared_derivative_row_prepare(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int evaluation_index,
  int omitted_observation,
  NPContinuousKernelRowResult *result)
{
  NPContinuousKernelBetaPreparedContext * const context =
    plan == NULL ? NULL : plan->beta_prepared;
  const int first_observation = omitted_observation == 0 ? 1 : 0;
  int coordinate;

  if(context == NULL || !context->pdf_active ||
     !context->pdf_row_component_active ||
     context->pdf_row_derivative_component == NULL ||
     segment == NULL || result == NULL ||
     context->num_train != plan->num_train ||
     context->num_continuous != plan->num_continuous)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(first_observation >= plan->num_train)
    return NP_CONTINUOUS_ROW_OK;

  for(coordinate = segment->coordinate_offset;
      coordinate < segment->coordinate_offset + segment->coordinate_count;
      ++coordinate) {
    const int local_coordinate = coordinate - segment->coordinate_offset;
    const double evaluation = plan->train_is_eval ?
      plan->train[coordinate][evaluation_index] :
      plan->evaluation[coordinate][evaluation_index];
    const double bandwidth = plan->bandwidth_mode == BW_FIXED ?
      plan->bandwidth_eval[coordinate][0] :
      plan->bandwidth_eval[coordinate][evaluation_index];
    const int beta_slot = context->coordinate_slot[coordinate];
    const int component_count = segment->descriptor.order / 2;
    int component;

    if(beta_slot < 0 || beta_slot >= context->num_beta_coordinates)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
    if(component_count <= 0 ||
       component_count > NP_BETA_PREPARED_MAX_COMPONENTS)
      return NP_CONTINUOUS_ROW_ERR_ROUTE;
    for(component = 0; component < component_count; ++component) {
      const size_t component_offset =
        (size_t)beta_slot * NP_BETA_PREPARED_MAX_COMPONENTS +
        (size_t)component;
      np_beta_pdf_derivative_component * const derivative_component =
        &context->pdf_row_derivative_component[component_offset];
      np_beta_shape shape;
      const np_beta_status beta_status = np_beta_shape_init(
        evaluation, segment->lower[local_coordinate], bandwidth,
        segment->lower[local_coordinate],
        segment->upper[local_coordinate], component + 1, &shape);

      if(beta_status != NP_BETA_OK) {
        result->beta_status = beta_status;
        result->bad_coordinate = coordinate;
        result->bad_observation = first_observation;
        return NP_CONTINUOUS_ROW_ERR_KERNEL;
      }
      derivative_component->support_length = shape.support_length;
      derivative_component->concentration = shape.concentration;
      derivative_component->target_unit = shape.target_unit;
      derivative_component->target_complement_unit =
        shape.target_complement_unit;
    }
  }
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

NP_CONTINUOUS_ROW_ALWAYS_INLINE NPContinuousKernelRowStatus
np_continuous_kernel_beta_log_value(
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

NP_CONTINUOUS_ROW_ALWAYS_INLINE NPContinuousKernelRowStatus
np_continuous_kernel_beta_cdf_prepared_segment_log_fill_core(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int evaluation_index,
  int omitted_observation,
  double *log_absolute,
  signed char *sign,
  double *maximum,
  NPContinuousKernelRowResult *result,
  int use_prepared_concentration)
{
  NPContinuousKernelBetaPreparedContext * const prepared =
    plan == NULL ? NULL : plan->beta_prepared;
  int observation;

  if(prepared == NULL || !prepared->cdf_active ||
     prepared->cdf_coordinate_slot == NULL ||
     prepared->cdf_observation == NULL ||
     prepared->cdf_observation_status == NULL ||
     prepared->cdf_log_abs_coefficient == NULL ||
     prepared->cdf_coefficient_sign == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(use_prepared_concentration &&
     (!prepared->cdf_concentration_active ||
      prepared->cdf_concentration == NULL ||
      prepared->cdf_concentration_status == NULL))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  *maximum = -INFINITY;
  for(observation = 0; observation < plan->num_train; ++observation) {
    double log_product = 0.0;
    int product_sign = 1;

    if(observation == omitted_observation) {
      log_absolute[observation] = -INFINITY;
      sign[observation] = 0;
      continue;
    }
    for(int coordinate = segment->coordinate_offset;
        coordinate < segment->coordinate_offset + segment->coordinate_count;
        ++coordinate) {
      const int local_coordinate = coordinate - segment->coordinate_offset;
      double scalar_log = -INFINITY;
      int scalar_sign = 0;
      NPContinuousKernelRowStatus row_status;

      if(plan->operator[coordinate] == OP_INTEGRAL) {
        const int cdf_slot = prepared->cdf_coordinate_slot[coordinate];
        size_t observation_offset;
        size_t component_offset;
        double evaluation_bandwidth;
        double observation_bandwidth;

        if(cdf_slot < 0 || cdf_slot >= prepared->num_cdf_coordinates) {
          result->bad_coordinate = coordinate;
          result->bad_observation = observation;
          return NP_CONTINUOUS_ROW_ERR_LAYOUT;
        }
        observation_offset =
          (size_t)cdf_slot * (size_t)plan->num_train +
          (size_t)observation;
        component_offset =
          (size_t)cdf_slot * NP_BETA_PREPARED_MAX_COMPONENTS;
        row_status = np_continuous_kernel_row_bandwidths(
          plan, coordinate, evaluation_index, observation, 0,
          &evaluation_bandwidth, &observation_bandwidth);
        if(row_status != NP_CONTINUOUS_ROW_OK) {
          result->bad_coordinate = coordinate;
          result->bad_observation = observation;
          return row_status;
        }
        (void)observation_bandwidth;
        if(use_prepared_concentration) {
          scalar_log =
            np_beta_log_abs_cdf_order_prepared_concentration(
              plan->train_is_eval ?
                plan->train[coordinate][evaluation_index] :
                plan->evaluation[coordinate][evaluation_index],
              evaluation_bandwidth,
              segment->lower[local_coordinate],
              segment->upper[local_coordinate],
              segment->descriptor.order,
              &prepared->cdf_observation[observation_offset],
              prepared->cdf_observation_status[observation_offset],
              &prepared->cdf_log_abs_coefficient[component_offset],
              &prepared->cdf_coefficient_sign[component_offset],
              &prepared->cdf_concentration[component_offset],
              &prepared->cdf_concentration_status[component_offset],
              &scalar_sign, &result->beta_status);
        } else {
          scalar_log = np_beta_log_abs_cdf_order_prepared_observation(
            plan->train_is_eval ?
              plan->train[coordinate][evaluation_index] :
              plan->evaluation[coordinate][evaluation_index],
            evaluation_bandwidth,
            segment->lower[local_coordinate],
            segment->upper[local_coordinate], segment->descriptor.order,
            &prepared->cdf_observation[observation_offset],
            prepared->cdf_observation_status[observation_offset],
            &prepared->cdf_log_abs_coefficient[component_offset],
            &prepared->cdf_coefficient_sign[component_offset],
            &scalar_sign, &result->beta_status);
        }
        row_status = result->beta_status == NP_BETA_OK ?
          NP_CONTINUOUS_ROW_OK : NP_CONTINUOUS_ROW_ERR_KERNEL;
        if(row_status == NP_CONTINUOUS_ROW_OK &&
           (ISNAN(scalar_log) || scalar_log == INFINITY ||
            (scalar_sign != -1 && scalar_sign != 0 && scalar_sign != 1) ||
            ((scalar_sign == 0) != (scalar_log == -INFINITY))))
          row_status = NP_CONTINUOUS_ROW_ERR_NUMERIC;
      } else {
        row_status = np_continuous_kernel_beta_log_value(
          plan, segment, local_coordinate, coordinate,
          evaluation_index, observation,
          &scalar_log, &scalar_sign, &result->beta_status);
      }
      if(row_status != NP_CONTINUOUS_ROW_OK) {
        result->bad_coordinate = coordinate;
        result->bad_observation = observation;
        return row_status;
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
    log_absolute[observation] = log_product;
    sign[observation] = (signed char)product_sign;
    if(product_sign != 0 && log_product > *maximum)
      *maximum = log_product;
  }
  return NP_CONTINUOUS_ROW_OK;
}

static NPContinuousKernelRowStatus NP_CONTINUOUS_ROW_NOINLINE
np_continuous_kernel_beta_cdf_prepared_segment_log_fill(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int evaluation_index,
  int omitted_observation,
  double *log_absolute,
  signed char *sign,
  double *maximum,
  NPContinuousKernelRowResult *result)
{
  return np_continuous_kernel_beta_cdf_prepared_segment_log_fill_core(
    plan, segment, evaluation_index, omitted_observation,
    log_absolute, sign, maximum, result, 0);
}

static NPContinuousKernelRowStatus NP_CONTINUOUS_ROW_NOINLINE
np_continuous_kernel_beta_cdf_prepared_concentration_segment_log_fill(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int evaluation_index,
  int omitted_observation,
  double *log_absolute,
  signed char *sign,
  double *maximum,
  NPContinuousKernelRowResult *result)
{
  return np_continuous_kernel_beta_cdf_prepared_segment_log_fill_core(
    plan, segment, evaluation_index, omitted_observation,
    log_absolute, sign, maximum, result, 1);
}

typedef struct {
  double other_log;
  double regular_log;
  double jump_log;
  int other_sign;
  int regular_sign;
  int jump_sign;
  int has_derivative;
} NPContinuousKernelDerivativeComponents;

NP_CONTINUOUS_ROW_ALWAYS_INLINE NPContinuousKernelRowStatus
np_continuous_kernel_beta_derivative_value(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int local_coordinate,
  int coordinate,
  int evaluation_index,
  int observation,
  np_beta_derivative *derivative,
  np_beta_status *beta_status)
{
  const double evaluation = plan->train_is_eval ?
    plan->train[coordinate][evaluation_index] :
    plan->evaluation[coordinate][evaluation_index];
  const double observed = plan->train[coordinate][observation];
  double evaluation_bandwidth;
  double observation_bandwidth;
  NPContinuousKernelRowStatus status;

  if(plan->operator[coordinate] != OP_DERIVATIVE)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  status = np_continuous_kernel_row_bandwidths(
    plan, coordinate, evaluation_index, observation, 0,
    &evaluation_bandwidth, &observation_bandwidth);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  (void)observation_bandwidth;
  *beta_status = np_beta_pdf_derivative_order(
    evaluation, observed, evaluation_bandwidth,
    segment->lower[local_coordinate], segment->upper[local_coordinate],
    segment->descriptor.order, derivative);
  return *beta_status == NP_BETA_OK ? NP_CONTINUOUS_ROW_OK :
    NP_CONTINUOUS_ROW_ERR_KERNEL;
}

NP_CONTINUOUS_ROW_ALWAYS_INLINE NPContinuousKernelRowStatus
np_continuous_kernel_beta_derivative_segment_components(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int evaluation_index,
  int observation,
  int derivative_coordinate,
  NPContinuousKernelDerivativeComponents *components,
  int *bad_coordinate,
  np_beta_status *beta_status)
{
  const int segment_has_derivative =
    segment->descriptor.family == NP_CKERNEL_FAMILY_BETA &&
    derivative_coordinate >= segment->coordinate_offset &&
    derivative_coordinate <
      segment->coordinate_offset + segment->coordinate_count;
  NPContinuousKernelRowStatus status;
  int coordinate;

  components->other_log = 0.0;
  components->regular_log = -INFINITY;
  components->jump_log = -INFINITY;
  components->other_sign = 1;
  components->regular_sign = 0;
  components->jump_sign = 0;
  components->has_derivative = segment_has_derivative;
  if(segment->descriptor.family != NP_CKERNEL_FAMILY_BETA)
    return NP_CONTINUOUS_ROW_OK;

  for(coordinate = segment->coordinate_offset;
      coordinate < segment->coordinate_offset + segment->coordinate_count;
      ++coordinate) {
    const int local_coordinate = coordinate - segment->coordinate_offset;

    if(coordinate == derivative_coordinate) {
      np_beta_derivative derivative;

      status = np_continuous_kernel_beta_derivative_value(
        plan, segment, local_coordinate, coordinate, evaluation_index,
        observation, &derivative, beta_status);
      if(status != NP_CONTINUOUS_ROW_OK) {
        *bad_coordinate = coordinate;
        return status;
      }
      components->regular_log = derivative.regular_log_absolute;
      components->regular_sign = derivative.regular_sign;
      components->jump_log = derivative.jump_log_absolute;
      components->jump_sign = derivative.jump_sign;
    } else {
      double scalar_log = -INFINITY;
      int scalar_sign = 0;

      status = np_continuous_kernel_beta_log_value(
        plan, segment, local_coordinate, coordinate,
        evaluation_index, observation,
        &scalar_log, &scalar_sign, beta_status);
      if(status != NP_CONTINUOUS_ROW_OK) {
        *bad_coordinate = coordinate;
        return status;
      }
      if(scalar_sign == 0) {
        components->other_log = -INFINITY;
        components->other_sign = 0;
        break;
      }
      components->other_log += scalar_log;
      components->other_sign *= scalar_sign;
    }
  }

  if(!segment_has_derivative) {
    components->regular_log = 0.0;
    components->jump_log = 0.0;
    components->regular_sign = 1;
    components->jump_sign = 1;
  }
  if(ISNAN(components->other_log) || components->other_log == INFINITY ||
     ISNAN(components->regular_log) ||
     components->regular_log == INFINITY ||
     ISNAN(components->jump_log) || components->jump_log == INFINITY)
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

static NPContinuousKernelRowStatus
np_continuous_kernel_beta_segment_log_fill(
  const NPContinuousKernelRowPlan *plan,
  const NPContinuousKernelSegment *segment,
  int evaluation_index,
  int omitted_observation,
  double *log_absolute,
  signed char *sign,
  double *maximum,
  NPContinuousKernelRowResult *result)
{
  int normal_operator = 1;
  int observation;

  if(plan == NULL || segment == NULL || log_absolute == NULL ||
     sign == NULL || maximum == NULL || result == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(plan->beta_prepared != NULL &&
     plan->beta_prepared->cdf_concentration_active)
    return
      np_continuous_kernel_beta_cdf_prepared_concentration_segment_log_fill(
        plan, segment, evaluation_index, omitted_observation,
        log_absolute, sign, maximum, result);
  if(plan->beta_prepared != NULL && plan->beta_prepared->cdf_active)
    return np_continuous_kernel_beta_cdf_prepared_segment_log_fill(
      plan, segment, evaluation_index, omitted_observation,
      log_absolute, sign, maximum, result);
  *maximum = -INFINITY;
  for(int coordinate = segment->coordinate_offset;
      coordinate < segment->coordinate_offset + segment->coordinate_count;
      ++coordinate) {
    normal_operator &= plan->operator[coordinate] == OP_NORMAL;
    if((plan->bandwidth_mode == BW_FIXED ||
        plan->bandwidth_mode == BW_GEN_NN) &&
       (plan->bandwidth_eval == NULL ||
        plan->bandwidth_eval[coordinate] == NULL))
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
    if(plan->bandwidth_mode == BW_ADAP_NN &&
       (plan->bandwidth_train == NULL ||
        plan->bandwidth_train[coordinate] == NULL))
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  }

  if(normal_operator) {
    NPContinuousKernelBetaPreparedContext * const prepared =
      plan->beta_prepared;

    if(prepared != NULL && prepared->pdf_active &&
       prepared->pdf_row_component_active) {
      const int prepared_component_count =
        segment->descriptor.order / 2;
          const NPContinuousKernelRowStatus prepare_status =
        np_continuous_kernel_beta_prepared_pdf_row_prepare(
          plan, segment, evaluation_index, omitted_observation, result);

      if(prepare_status != NP_CONTINUOUS_ROW_OK)
        return prepare_status;
      for(observation = 0; observation < plan->num_train; ++observation) {
        double log_product = 0.0;
        int product_sign = 1;

        if(observation == omitted_observation) {
          log_absolute[observation] = -INFINITY;
          sign[observation] = 0;
          continue;
        }
        for(int coordinate = segment->coordinate_offset;
            coordinate < segment->coordinate_offset +
              segment->coordinate_count;
            ++coordinate) {
          const int beta_slot = prepared->coordinate_slot[coordinate];
          const np_beta_pdf_component *components;
          const np_beta_pdf_observation *prepared_observation;
          const size_t observation_offset =
            (size_t)beta_slot * (size_t)plan->num_train +
            (size_t)observation;
          int scalar_sign = 0;

          if(beta_slot < 0 ||
             beta_slot >= prepared->num_beta_coordinates) {
            result->bad_coordinate = coordinate;
            result->bad_observation = observation;
            return NP_CONTINUOUS_ROW_ERR_LAYOUT;
          }
          components = &prepared->pdf_row_component[
            (size_t)beta_slot * NP_BETA_PREPARED_MAX_COMPONENTS];
          if(prepared->pdf_observation_status[observation_offset] !=
               NP_BETA_OK) {
            result->beta_status =
              prepared->pdf_observation_status[observation_offset];
            result->bad_coordinate = coordinate;
            result->bad_observation = observation;
            return NP_CONTINUOUS_ROW_ERR_KERNEL;
          }
          prepared_observation =
            &prepared->pdf_observation[observation_offset];

          log_absolute[observation] = np_beta_log_abs_pdf_prepared(
            components, prepared_observation, prepared_component_count,
            &scalar_sign, &result->beta_status);
          if(result->beta_status != NP_BETA_OK) {
            result->bad_coordinate = coordinate;
            result->bad_observation = observation;
            return NP_CONTINUOUS_ROW_ERR_KERNEL;
          }
          if(ISNAN(log_absolute[observation]) ||
             log_absolute[observation] == INFINITY ||
             (scalar_sign != -1 && scalar_sign != 0 && scalar_sign != 1) ||
             ((scalar_sign == 0) !=
              (log_absolute[observation] == -INFINITY))) {
            result->bad_coordinate = coordinate;
            result->bad_observation = observation;
            return NP_CONTINUOUS_ROW_ERR_NUMERIC;
          }
          if(scalar_sign == 0) {
            product_sign = 0;
            log_product = -INFINITY;
            break;
          }
          log_product += log_absolute[observation];
          product_sign *= scalar_sign;
          if(ISNAN(log_product) || log_product == INFINITY) {
            result->bad_coordinate = coordinate;
            result->bad_observation = observation;
            return NP_CONTINUOUS_ROW_ERR_NUMERIC;
          }
        }
        log_absolute[observation] = log_product;
        sign[observation] = (signed char)product_sign;
        if(product_sign != 0 && log_product > *maximum)
          *maximum = log_product;
      }
      return NP_CONTINUOUS_ROW_OK;
    }

    if(prepared != NULL && prepared->pdf_active &&
       plan->bandwidth_mode == BW_ADAP_NN) {
      for(observation = 0; observation < plan->num_train; ++observation) {
        double log_product = 0.0;
        int product_sign = 1;

        if(observation == omitted_observation) {
          log_absolute[observation] = -INFINITY;
          sign[observation] = 0;
          continue;
        }
        for(int coordinate = segment->coordinate_offset;
            coordinate < segment->coordinate_offset +
              segment->coordinate_count;
            ++coordinate) {
          const int local_coordinate =
            coordinate - segment->coordinate_offset;
          const int beta_slot = prepared->coordinate_slot[coordinate];
          const double evaluation = plan->train_is_eval ?
            plan->train[coordinate][evaluation_index] :
            plan->evaluation[coordinate][evaluation_index];
          size_t observation_offset;
          size_t component_offset;
          int scalar_sign = 0;

          if(beta_slot < 0 ||
             beta_slot >= prepared->num_beta_coordinates) {
            result->bad_coordinate = coordinate;
            result->bad_observation = observation;
            return NP_CONTINUOUS_ROW_ERR_LAYOUT;
          }
          observation_offset =
            (size_t)beta_slot * (size_t)plan->num_train +
            (size_t)observation;
          component_offset =
            (size_t)beta_slot * NP_BETA_PREPARED_MAX_COMPONENTS;
          log_absolute[observation] =
            np_beta_log_abs_pdf_order_prepared_observation(
              evaluation,
              plan->bandwidth_train[coordinate][observation],
              segment->lower[local_coordinate],
              segment->upper[local_coordinate],
              segment->descriptor.order,
              &prepared->pdf_observation[observation_offset],
              prepared->pdf_observation_status[observation_offset],
              &prepared->pdf_log_abs_coefficient[component_offset],
              &prepared->pdf_coefficient_sign[component_offset],
              &scalar_sign, &result->beta_status);
          if(result->beta_status != NP_BETA_OK) {
            result->bad_coordinate = coordinate;
            result->bad_observation = observation;
            return NP_CONTINUOUS_ROW_ERR_KERNEL;
          }
          if(ISNAN(log_absolute[observation]) ||
             log_absolute[observation] == INFINITY ||
             (scalar_sign != -1 && scalar_sign != 0 &&
              scalar_sign != 1) ||
             ((scalar_sign == 0) !=
              (log_absolute[observation] == -INFINITY))) {
            result->bad_coordinate = coordinate;
            result->bad_observation = observation;
            return NP_CONTINUOUS_ROW_ERR_NUMERIC;
          }
          if(scalar_sign == 0) {
            product_sign = 0;
            log_product = -INFINITY;
            break;
          }
          log_product += log_absolute[observation];
          product_sign *= scalar_sign;
          if(ISNAN(log_product) || log_product == INFINITY) {
            result->bad_coordinate = coordinate;
            result->bad_observation = observation;
            return NP_CONTINUOUS_ROW_ERR_NUMERIC;
          }
        }
        log_absolute[observation] = log_product;
        sign[observation] = (signed char)product_sign;
        if(product_sign != 0 && log_product > *maximum)
          *maximum = log_product;
      }
      return NP_CONTINUOUS_ROW_OK;
    }

    for(observation = 0; observation < plan->num_train; ++observation) {
      double log_product = 0.0;
      int product_sign = 1;

      if(observation == omitted_observation) {
        log_absolute[observation] = -INFINITY;
        sign[observation] = 0;
        continue;
      }
      for(int coordinate = segment->coordinate_offset;
          coordinate < segment->coordinate_offset + segment->coordinate_count;
          ++coordinate) {
        const int local_coordinate = coordinate - segment->coordinate_offset;
        const double evaluation = plan->train_is_eval ?
          plan->train[coordinate][evaluation_index] :
          plan->evaluation[coordinate][evaluation_index];
        const double bandwidth = plan->bandwidth_mode == BW_FIXED ?
          plan->bandwidth_eval[coordinate][0] :
          (plan->bandwidth_mode == BW_GEN_NN ?
           plan->bandwidth_eval[coordinate][evaluation_index] :
           plan->bandwidth_train[coordinate][observation]);
        int scalar_sign = 0;

        log_absolute[observation] = np_beta_log_abs_pdf_order(
          evaluation, plan->train[coordinate][observation], bandwidth,
          segment->lower[local_coordinate], segment->upper[local_coordinate],
          segment->descriptor.order, &scalar_sign, &result->beta_status);
        if(result->beta_status != NP_BETA_OK) {
          result->bad_coordinate = coordinate;
          result->bad_observation = observation;
          return NP_CONTINUOUS_ROW_ERR_KERNEL;
        }
        if(ISNAN(log_absolute[observation]) ||
           log_absolute[observation] == INFINITY ||
           (scalar_sign != -1 && scalar_sign != 0 && scalar_sign != 1) ||
           ((scalar_sign == 0) !=
            (log_absolute[observation] == -INFINITY))) {
          result->bad_coordinate = coordinate;
          result->bad_observation = observation;
          return NP_CONTINUOUS_ROW_ERR_NUMERIC;
        }
        if(scalar_sign == 0) {
          product_sign = 0;
          log_product = -INFINITY;
          break;
        }
        log_product += log_absolute[observation];
        product_sign *= scalar_sign;
        if(ISNAN(log_product) || log_product == INFINITY) {
          result->bad_coordinate = coordinate;
          result->bad_observation = observation;
          return NP_CONTINUOUS_ROW_ERR_NUMERIC;
        }
      }
      log_absolute[observation] = log_product;
      sign[observation] = (signed char)product_sign;
      if(product_sign != 0 && log_product > *maximum)
        *maximum = log_product;
    }
    return NP_CONTINUOUS_ROW_OK;
  }

  for(observation = 0; observation < plan->num_train; ++observation) {
    double log_product = 0.0;
    int product_sign = 1;

    if(observation == omitted_observation) {
      log_absolute[observation] = -INFINITY;
      sign[observation] = 0;
      continue;
    }
    for(int coordinate = segment->coordinate_offset;
        coordinate < segment->coordinate_offset + segment->coordinate_count;
        ++coordinate) {
      const int local_coordinate = coordinate - segment->coordinate_offset;
      double scalar_log = -INFINITY;
      int scalar_sign = 0;
      NPContinuousKernelRowStatus status;

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
    log_absolute[observation] = log_product;
    sign[observation] = (signed char)product_sign;
    if(product_sign != 0 && log_product > *maximum)
      *maximum = log_product;
  }
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

  status = np_continuous_kernel_row_plan_validate(plan);
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

    if(segment->descriptor.family != NP_CKERNEL_FAMILY_BETA)
      continue;
    status = np_continuous_kernel_beta_segment_log_fill(
      plan, segment, evaluation_index, omitted_observation,
      workspace->primary_log_absolute, workspace->primary_sign,
      &result->segment_log_scale[segment_index], result);
    if(status != NP_CONTINUOUS_ROW_OK)
      return status;

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

static NPContinuousKernelRowStatus
np_continuous_kernel_log_channels_multiply(
  const NPContinuousKernelRowPlan *plan,
  int omitted_observation,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *result)
{
  int observation;

  if(plan == NULL || workspace == NULL || result == NULL ||
     workspace->primary_log_absolute == NULL ||
     workspace->primary_sign == NULL ||
     workspace->secondary_log_absolute == NULL ||
     workspace->secondary_sign == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;

  for(observation = 0; observation < plan->num_train; ++observation) {
    const double factor_log =
      workspace->secondary_log_absolute[observation];
    const int factor_sign = workspace->secondary_sign[observation];
    const double beta_log = workspace->primary_log_absolute[observation];
    const int beta_sign = workspace->primary_sign[observation];

    if(observation == omitted_observation) {
      workspace->primary_log_absolute[observation] = -INFINITY;
      workspace->primary_sign[observation] = 0;
      continue;
    }
    if(ISNAN(factor_log) || factor_log == INFINITY ||
       (factor_sign != -1 && factor_sign != 0 && factor_sign != 1) ||
       ((factor_sign == 0) != (factor_log == -INFINITY)) ||
       ISNAN(beta_log) || beta_log == INFINITY ||
       (beta_sign != -1 && beta_sign != 0 && beta_sign != 1) ||
       ((beta_sign == 0) != (beta_log == -INFINITY))) {
      result->bad_observation = observation;
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
    }
    if(beta_sign == 0 || factor_sign == 0) {
      workspace->primary_log_absolute[observation] = -INFINITY;
      workspace->primary_sign[observation] = 0;
      continue;
    }
    workspace->primary_log_absolute[observation] = beta_log + factor_log;
    workspace->primary_sign[observation] =
      (signed char)(beta_sign * factor_sign);
    if(ISNAN(workspace->primary_log_absolute[observation]) ||
       workspace->primary_log_absolute[observation] == INFINITY) {
      result->bad_observation = observation;
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
    }
  }
  return NP_CONTINUOUS_ROW_OK;
}

static NPContinuousKernelRowStatus
np_continuous_kernel_log_factor_compose(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  const NPContinuousKernelLogFactorProvider *provider,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *result)
{
  NPContinuousKernelRowStatus status;

  if(plan == NULL || provider == NULL || provider->function == NULL ||
     workspace == NULL || result == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  status = np_continuous_kernel_row_workspace_reserve(
    workspace, (size_t)plan->num_train, 1);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  status = provider->function(
    provider->context, evaluation_index, omitted_observation,
    plan->num_train, workspace->secondary_log_absolute,
    workspace->secondary_sign);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  return np_continuous_kernel_log_channels_multiply(
    plan, omitted_observation, workspace, result);
}

static NPContinuousKernelRowStatus NP_CONTINUOUS_ROW_NOINLINE
np_continuous_kernel_beta_log_factor_row_multi_prevalidated(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  const NPContinuousKernelLogFactorProvider *provider,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *result)
{
  NPContinuousKernelRowStatus status;
  int beta_segment_count = 0;
  int segment_index;

  if(evaluation_index < 0 || evaluation_index >= plan->num_eval ||
     omitted_observation < -1 || omitted_observation >= plan->num_train)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(segment_index = 0; segment_index < plan->route->segment_count;
      ++segment_index)
    beta_segment_count +=
      plan->route->segment[segment_index].descriptor.family ==
      NP_CKERNEL_FAMILY_BETA;
  if(beta_segment_count == 0)
    return NP_CONTINUOUS_ROW_ERR_ROUTE;
  status = np_continuous_kernel_row_workspace_reserve(
    workspace, (size_t)plan->num_train,
    provider != NULL || beta_segment_count > 1);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;

  beta_segment_count = 0;
  for(segment_index = 0; segment_index < plan->route->segment_count;
      ++segment_index) {
    const NPContinuousKernelSegment * const segment =
      &plan->route->segment[segment_index];
    double * const segment_log = beta_segment_count == 0 ?
      workspace->primary_log_absolute :
      workspace->secondary_log_absolute;
    signed char * const segment_sign = beta_segment_count == 0 ?
      workspace->primary_sign : workspace->secondary_sign;

    if(segment->descriptor.family != NP_CKERNEL_FAMILY_BETA)
      continue;
    status = np_continuous_kernel_beta_segment_log_fill(
      plan, segment, evaluation_index, omitted_observation,
      segment_log, segment_sign,
      &result->segment_log_scale[segment_index], result);
    if(status != NP_CONTINUOUS_ROW_OK)
      return status;
    if(beta_segment_count > 0) {
      status = np_continuous_kernel_log_channels_multiply(
        plan, omitted_observation, workspace, result);
      if(status != NP_CONTINUOUS_ROW_OK)
        return status;
    }
    ++beta_segment_count;
  }
  if(provider == NULL) {
    result->total_log_scale = np_continuous_kernel_row_maximum(
      workspace->primary_log_absolute, workspace->primary_sign,
      plan->num_train);
    return NP_CONTINUOUS_ROW_OK;
  }
  status = np_continuous_kernel_log_factor_compose(
    plan, evaluation_index, omitted_observation,
    provider, workspace, result);
  if(status == NP_CONTINUOUS_ROW_OK)
    result->total_log_scale = np_continuous_kernel_row_maximum(
      workspace->primary_log_absolute, workspace->primary_sign,
      plan->num_train);
  return status;
}

static NPContinuousKernelRowStatus
np_continuous_kernel_beta_log_factor_row(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  const NPContinuousKernelLogFactorProvider *provider,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *result)
{
  NPContinuousKernelRowStatus status;
  int segment_index;

  if(plan == NULL || plan->route == NULL || workspace == NULL ||
     result == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  np_continuous_kernel_row_reset_diagnostics(
    &result->bad_coordinate, &result->bad_observation, &result->beta_status);
  result->total_log_scale = 0.0;
  for(segment_index = 0; segment_index < NP_CKERNEL_ROUTE_MAX_SEGMENTS;
      ++segment_index)
    result->segment_log_scale[segment_index] = 0.0;

  status = np_continuous_kernel_row_plan_validate(plan);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  if(plan->route->segment_count != 1 ||
     plan->route->segment[0].descriptor.family != NP_CKERNEL_FAMILY_BETA)
    return np_continuous_kernel_beta_log_factor_row_multi_prevalidated(
      plan, evaluation_index, omitted_observation,
      provider, workspace, result);
  if(evaluation_index < 0 || evaluation_index >= plan->num_eval ||
     omitted_observation < -1 || omitted_observation >= plan->num_train)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  status = np_continuous_kernel_row_workspace_reserve(
    workspace, (size_t)plan->num_train, provider != NULL);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  status = np_continuous_kernel_beta_segment_log_fill(
    plan, &plan->route->segment[0], evaluation_index,
    omitted_observation, workspace->primary_log_absolute,
    workspace->primary_sign, &result->segment_log_scale[0], result);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  if(provider == NULL) {
    result->total_log_scale = result->segment_log_scale[0];
    return NP_CONTINUOUS_ROW_OK;
  }
  status = np_continuous_kernel_log_factor_compose(
    plan, evaluation_index, omitted_observation,
    provider, workspace, result);
  if(status == NP_CONTINUOUS_ROW_OK)
    result->total_log_scale = np_continuous_kernel_row_maximum(
      workspace->primary_log_absolute, workspace->primary_sign,
      plan->num_train);
  return status;
}

NPContinuousKernelRowStatus
np_continuous_kernel_beta_factor_row_with_log_factor(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  const NPContinuousKernelLogFactorProvider *provider,
  NPContinuousKernelRowWorkspace *workspace,
  NPContinuousKernelRowResult *result)
{
  NPContinuousKernelRowStatus status;
  double complete_log_scale;
  int observation;

  if(result == NULL || result->row == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  status = np_continuous_kernel_beta_log_factor_row(
    plan, evaluation_index, omitted_observation,
    provider, workspace, result);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;

  complete_log_scale = result->total_log_scale;
  for(observation = 0; observation < plan->num_train; ++observation)
    result->row[observation] = 1.0;
  status = np_continuous_kernel_row_multiply_scaled(
    result->row, workspace->primary_log_absolute,
    workspace->primary_sign, plan->num_train, complete_log_scale);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  result->total_log_scale = complete_log_scale;
  return NP_CONTINUOUS_ROW_OK;
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

  status = np_continuous_kernel_row_plan_validate(plan);
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

    if(segment->descriptor.family != NP_CKERNEL_FAMILY_BETA)
      continue;
    for(observation = 0; observation < plan->num_train; ++observation) {
      NPContinuousKernelDerivativeComponents components;
      double regular_log;
      double jump_log;
      int regular_sign;
      int jump_sign;

      if(observation == omitted_observation) {
        workspace->primary_log_absolute[observation] = -INFINITY;
        workspace->secondary_log_absolute[observation] = -INFINITY;
        workspace->primary_sign[observation] = 0;
        workspace->secondary_sign[observation] = 0;
        continue;
      }
      status = np_continuous_kernel_beta_derivative_segment_components(
        plan, segment, evaluation_index, observation, derivative_coordinate,
        &components,
        &result->bad_coordinate, &result->beta_status);
      if(status != NP_CONTINUOUS_ROW_OK) {
        result->bad_observation = observation;
        return status;
      }
      regular_sign = components.other_sign * components.regular_sign;
      jump_sign = components.other_sign * components.jump_sign;
      regular_log = regular_sign == 0 ? -INFINITY :
        components.other_log + components.regular_log;
      jump_log = jump_sign == 0 ? -INFINITY :
        components.other_log + components.jump_log;

      workspace->primary_log_absolute[observation] = regular_log;
      workspace->secondary_log_absolute[observation] = jump_log;
      workspace->primary_sign[observation] = (signed char)regular_sign;
      workspace->secondary_sign[observation] = (signed char)jump_sign;
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

void np_continuous_kernel_derivative_accumulator_init(
  NPContinuousKernelDerivativeAccumulator *accumulator)
{
  if(accumulator == NULL)
    return;
  accumulator->regular_positive_log = NULL;
  accumulator->regular_negative_log = NULL;
  accumulator->jump_positive_log = NULL;
  accumulator->jump_negative_log = NULL;
  accumulator->capacity = 0;
}

void np_continuous_kernel_derivative_accumulator_release(
  NPContinuousKernelDerivativeAccumulator *accumulator)
{
  if(accumulator == NULL)
    return;
  free(accumulator->regular_positive_log);
  free(accumulator->regular_negative_log);
  free(accumulator->jump_positive_log);
  free(accumulator->jump_negative_log);
  np_continuous_kernel_derivative_accumulator_init(accumulator);
}

static NPContinuousKernelRowStatus
np_continuous_kernel_derivative_accumulator_reserve(
  NPContinuousKernelDerivativeAccumulator *accumulator,
  size_t width)
{
  double *regular_positive_log;
  double *regular_negative_log;
  double *jump_positive_log;
  double *jump_negative_log;

  if(accumulator == NULL || width == 0 || width > SIZE_MAX / sizeof(double))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(accumulator->capacity >= width &&
     accumulator->regular_positive_log != NULL &&
     accumulator->regular_negative_log != NULL &&
     accumulator->jump_positive_log != NULL &&
     accumulator->jump_negative_log != NULL)
    return NP_CONTINUOUS_ROW_OK;

  regular_positive_log = (double *)malloc(width * sizeof(double));
  regular_negative_log = (double *)malloc(width * sizeof(double));
  jump_positive_log = (double *)malloc(width * sizeof(double));
  jump_negative_log = (double *)malloc(width * sizeof(double));
  if(regular_positive_log == NULL || regular_negative_log == NULL ||
     jump_positive_log == NULL || jump_negative_log == NULL) {
    free(regular_positive_log);
    free(regular_negative_log);
    free(jump_positive_log);
    free(jump_negative_log);
    return NP_CONTINUOUS_ROW_ERR_MEMORY;
  }

  free(accumulator->regular_positive_log);
  free(accumulator->regular_negative_log);
  free(accumulator->jump_positive_log);
  free(accumulator->jump_negative_log);
  accumulator->regular_positive_log = regular_positive_log;
  accumulator->regular_negative_log = regular_negative_log;
  accumulator->jump_positive_log = jump_positive_log;
  accumulator->jump_negative_log = jump_negative_log;
  accumulator->capacity = width;
  return NP_CONTINUOUS_ROW_OK;
}

static double np_continuous_kernel_log_add(double accumulator, double term)
{
  double maximum;
  double minimum;

  if(accumulator == -INFINITY)
    return term;
  if(term == -INFINITY)
    return accumulator;
  maximum = fmax(accumulator, term);
  minimum = fmin(accumulator, term);
  return maximum + log1p(exp(minimum - maximum));
}

static void np_continuous_kernel_signed_log_add(
  double log_absolute,
  int sign,
  double *positive_log,
  double *negative_log)
{
  if(sign > 0)
    *positive_log = np_continuous_kernel_log_add(
      *positive_log, log_absolute);
  else if(sign < 0)
    *negative_log = np_continuous_kernel_log_add(
      *negative_log, log_absolute);
}

static double np_continuous_kernel_signed_log_value(
  double positive_log,
  double negative_log,
  int *undefined)
{
  double log_absolute = -INFINITY;
  int sign = 0;
  const np_beta_status status = np_beta_signed_log_absolute(
    positive_log, negative_log, &log_absolute, &sign);

  if(status != NP_BETA_OK || ISNAN(log_absolute)) {
    if(undefined != NULL)
      *undefined = 1;
    return NA_REAL;
  }
  if(sign == 0)
    return 0.0;
  if(log_absolute > log(DBL_MAX))
    return (sign > 0) ? INFINITY : -INFINITY;
  return (sign > 0) ? exp(log_absolute) : -exp(log_absolute);
}

NP_CONTINUOUS_ROW_ALWAYS_INLINE NPContinuousKernelRowStatus
np_continuous_kernel_beta_derivative_absolute_row_bound(
  const NPContinuousKernelRowPlan *plan,
  int evaluation_index,
  int omitted_observation,
  int derivative_coordinate,
  double * const *response,
  int response_columns,
  double * const *case_weights,
  int weight_columns,
  const double *factor_log_absolute,
  const signed char *factor_sign,
  NPContinuousKernelDerivativeAccumulator *accumulator,
  double *weighted_sum,
  double *kernel_weights,
  NPContinuousKernelDerivativeDiagnostics *diagnostics)
{
  const int response_extent = response_columns > 0 ? response_columns : 1;
  const int weight_extent = weight_columns > 0 ? weight_columns : 1;
  const NPContinuousKernelSegment *single_beta_segment = NULL;
  size_t sum_extent;
  NPContinuousKernelRowStatus status;
  int observation;
  int response_column;
  int weight_column;

  if((factor_log_absolute == NULL) != (factor_sign == NULL))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  sum_extent = (size_t)response_extent * (size_t)weight_extent;
  if(plan->route->segment_count == 1 &&
     plan->route->segment[0].descriptor.family == NP_CKERNEL_FAMILY_BETA)
    single_beta_segment = &plan->route->segment[0];
  for(size_t output = 0; output < sum_extent; ++output) {
    accumulator->regular_positive_log[output] = -INFINITY;
    accumulator->regular_negative_log[output] = -INFINITY;
    accumulator->jump_positive_log[output] = -INFINITY;
    accumulator->jump_negative_log[output] = -INFINITY;
  }

  for(observation = 0; observation < plan->num_train; ++observation) {
    double other_log = 0.0;
    double derivative_regular_log = -INFINITY;
    double derivative_jump_log = -INFINITY;
    int other_sign = 1;
    int derivative_regular_sign = 0;
    int derivative_jump_sign = 0;
    int derivative_found = 0;
    int segment_index;

    if(observation == omitted_observation) {
      if(kernel_weights != NULL)
        kernel_weights[observation] = 0.0;
      continue;
    }
    if(factor_log_absolute != NULL) {
      const double factor_log = factor_log_absolute[observation];
      const int sign = factor_sign[observation];

      if(ISNAN(factor_log) || factor_log == INFINITY ||
         (sign != -1 && sign != 0 && sign != 1) ||
         ((sign == 0) != (factor_log == -INFINITY)))
        return NP_CONTINUOUS_ROW_ERR_NUMERIC;
      other_log = factor_log;
      other_sign = sign;
    }
    if(single_beta_segment != NULL) {
      int coordinate;

      for(coordinate = single_beta_segment->coordinate_offset;
          coordinate < single_beta_segment->coordinate_offset +
            single_beta_segment->coordinate_count;
          ++coordinate) {
        const int local_coordinate =
          coordinate - single_beta_segment->coordinate_offset;
        np_beta_status beta_status = NP_BETA_OK;

        if(coordinate == derivative_coordinate) {
          np_beta_derivative derivative;

          status = np_continuous_kernel_beta_derivative_value(
            plan, single_beta_segment, local_coordinate, coordinate,
            evaluation_index, observation, &derivative, &beta_status);
          if(status == NP_CONTINUOUS_ROW_OK) {
            derivative_regular_log = derivative.regular_log_absolute;
            derivative_regular_sign = derivative.regular_sign;
            derivative_jump_log = derivative.jump_log_absolute;
            derivative_jump_sign = derivative.jump_sign;
            derivative_found = 1;
          }
        } else {
          double scalar_log = -INFINITY;
          int scalar_sign = 0;

          status = np_continuous_kernel_beta_log_value(
            plan, single_beta_segment, local_coordinate, coordinate,
            evaluation_index, observation,
            &scalar_log, &scalar_sign, &beta_status);
          if(status == NP_CONTINUOUS_ROW_OK) {
            if(scalar_sign == 0) {
              other_log = -INFINITY;
              other_sign = 0;
            } else if(other_sign != 0) {
              other_log += scalar_log;
              other_sign *= scalar_sign;
            }
          }
        }
        if(status != NP_CONTINUOUS_ROW_OK) {
          if(diagnostics != NULL) {
            diagnostics->bad_coordinate = coordinate;
            diagnostics->bad_observation = observation;
            diagnostics->beta_status = beta_status;
          }
          return status;
        }
      }
    } else {
      for(segment_index = 0; segment_index < plan->route->segment_count;
          ++segment_index) {
        NPContinuousKernelDerivativeComponents components;
        int bad_coordinate = -1;
        np_beta_status beta_status = NP_BETA_OK;

        status = np_continuous_kernel_beta_derivative_segment_components(
          plan, &plan->route->segment[segment_index], evaluation_index,
          observation, derivative_coordinate, &components,
          &bad_coordinate, &beta_status);
        if(status != NP_CONTINUOUS_ROW_OK) {
          if(diagnostics != NULL) {
            diagnostics->bad_coordinate = bad_coordinate;
            diagnostics->bad_observation = observation;
            diagnostics->beta_status = beta_status;
          }
          return status;
        }
        if(other_sign != 0 && components.other_sign != 0) {
          other_log += components.other_log;
          other_sign *= components.other_sign;
        } else {
          other_log = -INFINITY;
          other_sign = 0;
        }
        if(components.has_derivative) {
          if(derivative_found)
            return NP_CONTINUOUS_ROW_ERR_LAYOUT;
          derivative_found = 1;
          derivative_regular_log = components.regular_log;
          derivative_regular_sign = components.regular_sign;
          derivative_jump_log = components.jump_log;
          derivative_jump_sign = components.jump_sign;
        }
      }
    }
    if(!derivative_found)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;

    if(kernel_weights != NULL) {
      if(other_sign != 0 && derivative_jump_sign != 0)
        kernel_weights[observation] =
          other_sign * derivative_jump_sign > 0 ? INFINITY : -INFINITY;
      else if(other_sign == 0 || derivative_regular_sign == 0)
        kernel_weights[observation] = 0.0;
      else if(other_log + derivative_regular_log > log(DBL_MAX))
        kernel_weights[observation] =
          other_sign * derivative_regular_sign > 0 ? INFINITY : -INFINITY;
      else
        kernel_weights[observation] =
          other_sign * derivative_regular_sign > 0 ?
          exp(other_log + derivative_regular_log) :
          -exp(other_log + derivative_regular_log);
    }

    if(other_sign == 0)
      continue;

    if(response_columns == 0 && weight_columns == 0) {
      if(derivative_regular_sign != 0)
        np_continuous_kernel_signed_log_add(
          other_log + derivative_regular_log,
          other_sign * derivative_regular_sign,
          &accumulator->regular_positive_log[0],
          &accumulator->regular_negative_log[0]);
      if(derivative_jump_sign != 0)
        np_continuous_kernel_signed_log_add(
          other_log + derivative_jump_log,
          other_sign * derivative_jump_sign,
          &accumulator->jump_positive_log[0],
          &accumulator->jump_negative_log[0]);
      continue;
    }

    for(response_column = 0; response_column < response_extent;
        ++response_column) {
      const double response_value = response_columns > 0 ?
        response[response_column][observation] : 1.0;

      for(weight_column = 0; weight_column < weight_extent;
          ++weight_column) {
        const double weight_value = weight_columns > 0 ?
          case_weights[weight_column][observation] : 1.0;
        const size_t output =
          (size_t)response_column * (size_t)weight_extent +
          (size_t)weight_column;
        double multiplier_log;
        int multiplier_sign;

        if(!R_FINITE(response_value) || !R_FINITE(weight_value))
          return NP_CONTINUOUS_ROW_ERR_NUMERIC;
        if(response_value == 0.0 || weight_value == 0.0)
          continue;
        multiplier_sign = other_sign *
          ((response_value > 0.0) ? 1 : -1) *
          ((weight_value > 0.0) ? 1 : -1);
        multiplier_log = other_log +
          log(fabs(response_value)) + log(fabs(weight_value));
        if(derivative_regular_sign != 0)
          np_continuous_kernel_signed_log_add(
            multiplier_log + derivative_regular_log,
            multiplier_sign * derivative_regular_sign,
            &accumulator->regular_positive_log[output],
            &accumulator->regular_negative_log[output]);
        if(derivative_jump_sign != 0)
          np_continuous_kernel_signed_log_add(
            multiplier_log + derivative_jump_log,
            multiplier_sign * derivative_jump_sign,
            &accumulator->jump_positive_log[output],
            &accumulator->jump_negative_log[output]);
      }
    }
  }

  for(size_t output = 0; output < sum_extent; ++output) {
    int undefined = 0;
    const double jump = np_continuous_kernel_signed_log_value(
      accumulator->jump_positive_log[output],
      accumulator->jump_negative_log[output], &undefined);

    if(undefined) {
      weighted_sum[output] = NA_REAL;
      if(diagnostics != NULL)
        ++diagnostics->undefined_count;
    } else if(jump != 0.0) {
      weighted_sum[output] = jump > 0.0 ? INFINITY : -INFINITY;
    } else {
      weighted_sum[output] = np_continuous_kernel_signed_log_value(
        accumulator->regular_positive_log[output],
        accumulator->regular_negative_log[output], &undefined);
      if(undefined && diagnostics != NULL)
        ++diagnostics->undefined_count;
    }
  }
  return NP_CONTINUOUS_ROW_OK;
}

static NPContinuousKernelRowStatus
np_continuous_kernel_beta_derivative_absolute_contract_validate(
  const NPContinuousKernelRowPlan *plan,
  int leave_one_out,
  int leave_one_out_offset,
  int derivative_coordinate,
  double * const *response,
  int response_columns,
  double * const *case_weights,
  int weight_columns,
  NPContinuousKernelDerivativeAccumulator *accumulator,
  const double *weighted_sum,
  double *kernel_weights,
  size_t *sum_extent)
{
  const int response_extent = response_columns > 0 ? response_columns : 1;
  const int weight_extent = weight_columns > 0 ? weight_columns : 1;
  NPContinuousKernelRowStatus status;
  int derivative_count = 0;
  int index;

  if(plan == NULL || plan->route == NULL || plan->num_train <= 0 ||
     plan->num_eval <= 0 || plan->num_continuous <= 0 ||
     plan->operator == NULL ||
     (leave_one_out != 0 && leave_one_out != 1) ||
     leave_one_out_offset < 0 || derivative_coordinate < 0 ||
     derivative_coordinate >= plan->num_continuous ||
     response_columns < 0 || weight_columns < 0 ||
     (response_columns > 0 && response == NULL) ||
     (weight_columns > 0 && case_weights == NULL) ||
     accumulator == NULL || weighted_sum == NULL || sum_extent == NULL ||
     (size_t)response_extent > SIZE_MAX / (size_t)weight_extent)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(index = 0; index < response_columns; ++index)
    if(response[index] == NULL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(index = 0; index < weight_columns; ++index)
    if(case_weights[index] == NULL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  *sum_extent = (size_t)response_extent * (size_t)weight_extent;
  if((size_t)plan->num_eval > SIZE_MAX / *sum_extent ||
     (kernel_weights != NULL &&
      (size_t)plan->num_eval > SIZE_MAX / (size_t)plan->num_train) ||
     (leave_one_out &&
      (plan->num_eval > plan->num_train ||
       leave_one_out_offset > plan->num_train - plan->num_eval)))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(index = 0; index < plan->route->segment_count; ++index)
    if(plan->route->segment[index].descriptor.family !=
       NP_CKERNEL_FAMILY_BETA)
      return NP_CONTINUOUS_ROW_ERR_ROUTE;
  for(index = 0; index < plan->num_continuous; ++index)
    derivative_count += plan->operator[index] == OP_DERIVATIVE;
  if(derivative_count != 1 ||
     plan->operator[derivative_coordinate] != OP_DERIVATIVE)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  status = np_continuous_kernel_derivative_accumulator_reserve(
    accumulator, *sum_extent);
  return status;
}

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
  NPContinuousKernelDerivativeDiagnostics *diagnostics)
{
  size_t sum_extent;
  NPContinuousKernelRowStatus status;
  int evaluation;

  if(diagnostics != NULL) {
    diagnostics->bad_coordinate = -1;
    diagnostics->bad_observation = -1;
    diagnostics->undefined_count = 0;
    diagnostics->beta_status = NP_BETA_OK;
  }
  status = np_continuous_kernel_beta_derivative_absolute_contract_validate(
    plan, leave_one_out, leave_one_out_offset, derivative_coordinate,
    response, response_columns, case_weights, weight_columns, accumulator,
    weighted_sum, kernel_weights, &sum_extent);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;

  for(evaluation = 0; evaluation < plan->num_eval; ++evaluation) {
    const int omitted_observation = leave_one_out ?
      evaluation + leave_one_out_offset : -1;
    NPContinuousKernelDerivativeDiagnostics row_diagnostics = {
      -1, -1, 0, NP_BETA_OK
    };

    if(omitted_observation >= plan->num_train)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
    status = np_continuous_kernel_beta_derivative_absolute_row_bound(
      plan, evaluation, omitted_observation, derivative_coordinate,
      response, response_columns, case_weights, weight_columns,
      NULL, NULL, accumulator,
      weighted_sum + (size_t)evaluation * sum_extent,
      kernel_weights == NULL ? NULL :
        kernel_weights + (size_t)evaluation * (size_t)plan->num_train,
      &row_diagnostics);
    if(diagnostics != NULL) {
      diagnostics->undefined_count += row_diagnostics.undefined_count;
      if(status != NP_CONTINUOUS_ROW_OK) {
        diagnostics->bad_coordinate = row_diagnostics.bad_coordinate;
        diagnostics->bad_observation = row_diagnostics.bad_observation;
        diagnostics->beta_status = row_diagnostics.beta_status;
      }
    }
    if(status != NP_CONTINUOUS_ROW_OK)
      return status;
    if((evaluation & 255) == 0)
      R_CheckUserInterrupt();
  }
  return NP_CONTINUOUS_ROW_OK;
}

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
  NPContinuousKernelDerivativeDiagnostics *diagnostics)
{
  size_t sum_extent;
  NPContinuousKernelRowStatus status;
  int evaluation;

  if(diagnostics != NULL) {
    diagnostics->bad_coordinate = -1;
    diagnostics->bad_observation = -1;
    diagnostics->undefined_count = 0;
    diagnostics->beta_status = NP_BETA_OK;
  }
  if(provider == NULL || provider->function == NULL ||
     factor_workspace == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  status = np_continuous_kernel_beta_derivative_absolute_contract_validate(
    plan, leave_one_out, leave_one_out_offset, derivative_coordinate,
    response, response_columns, case_weights, weight_columns, accumulator,
    weighted_sum, kernel_weights, &sum_extent);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;

  status = np_continuous_kernel_row_workspace_reserve(
    factor_workspace, (size_t)plan->num_train, 0);
  if(status != NP_CONTINUOUS_ROW_OK)
    goto cleanup;
  for(evaluation = 0; evaluation < plan->num_eval; ++evaluation) {
    const int omitted_observation = leave_one_out ?
      evaluation + leave_one_out_offset : -1;
    NPContinuousKernelDerivativeDiagnostics row_diagnostics = {
      -1, -1, 0, NP_BETA_OK
    };

    if(omitted_observation >= plan->num_train) {
      status = NP_CONTINUOUS_ROW_ERR_LAYOUT;
      goto cleanup;
    }
    status = provider->function(
      provider->context, evaluation, omitted_observation,
      plan->num_train, factor_workspace->primary_log_absolute,
      factor_workspace->primary_sign);
    if(status != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
    status = np_continuous_kernel_beta_derivative_absolute_row_bound(
      plan, evaluation, omitted_observation, derivative_coordinate,
      response, response_columns, case_weights, weight_columns,
      factor_workspace->primary_log_absolute,
      factor_workspace->primary_sign, accumulator,
      weighted_sum + (size_t)evaluation * sum_extent,
      kernel_weights == NULL ? NULL :
        kernel_weights + (size_t)evaluation * (size_t)plan->num_train,
      &row_diagnostics);
    if(diagnostics != NULL) {
      diagnostics->undefined_count += row_diagnostics.undefined_count;
      if(status != NP_CONTINUOUS_ROW_OK) {
        diagnostics->bad_coordinate = row_diagnostics.bad_coordinate;
        diagnostics->bad_observation = row_diagnostics.bad_observation;
        diagnostics->beta_status = row_diagnostics.beta_status;
      }
    }
    if(status != NP_CONTINUOUS_ROW_OK)
      goto cleanup;
    if((evaluation & 255) == 0)
      R_CheckUserInterrupt();
  }
  status = NP_CONTINUOUS_ROW_OK;

cleanup:
  return status;
}

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
  NPContinuousKernelDerivativeDiagnostics *diagnostics)
{
  const int response_extent = response_columns > 0 ? response_columns : 1;
  const int weight_extent = weight_columns > 0 ? weight_columns : 1;
  NPContinuousKernelDerivativeRowResult row_result;
  size_t sum_extent;
  NPContinuousKernelRowStatus status = NP_CONTINUOUS_ROW_ERR_LAYOUT;
  int derivative_count = 0;
  int evaluation;
  int index;

  if(diagnostics != NULL) {
    diagnostics->bad_coordinate = -1;
    diagnostics->bad_observation = -1;
    diagnostics->undefined_count = 0;
    diagnostics->beta_status = NP_BETA_OK;
  }
  if(plan == NULL || plan->route == NULL || plan->num_train <= 0 ||
     plan->num_eval <= 0 || plan->num_continuous <= 0 ||
     plan->operator == NULL || plan->route->segment_count != 1 ||
     plan->route->segment[0].descriptor.family != NP_CKERNEL_FAMILY_BETA ||
     (leave_one_out != 0 && leave_one_out != 1) ||
     leave_one_out_offset < 0 || derivative_coordinate < 0 ||
     derivative_coordinate >= plan->num_continuous || kernel_power == 1 ||
     response_columns < 0 || weight_columns < 0 ||
     (response_columns > 0 && response == NULL) ||
     (weight_columns > 0 && case_weights == NULL) ||
     accumulator == NULL || workspace == NULL || row_storage == NULL ||
     weighted_sum == NULL ||
     (provider != NULL &&
      (factor_log_absolute == NULL || factor_sign == NULL)) ||
     (size_t)response_extent > SIZE_MAX / (size_t)weight_extent)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(index = 0; index < response_columns; ++index)
    if(response[index] == NULL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(index = 0; index < weight_columns; ++index)
    if(case_weights[index] == NULL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  sum_extent = (size_t)response_extent * (size_t)weight_extent;
  if((size_t)plan->num_eval > SIZE_MAX / sum_extent ||
     (kernel_weights != NULL &&
      (size_t)plan->num_eval > SIZE_MAX / (size_t)plan->num_train) ||
     (leave_one_out &&
      (plan->num_eval > plan->num_train ||
       leave_one_out_offset > plan->num_train - plan->num_eval)) ||
     (size_t)plan->num_train > SIZE_MAX / 2 ||
     (size_t)plan->num_train * 2 > SIZE_MAX / sizeof(double))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(index = 0; index < plan->num_continuous; ++index)
    derivative_count += plan->operator[index] == OP_DERIVATIVE;
  if(derivative_count != 1 ||
     plan->operator[derivative_coordinate] != OP_DERIVATIVE)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  status = np_continuous_kernel_derivative_accumulator_reserve(
    accumulator, sum_extent);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;

  row_result.regular_row = row_storage;
  row_result.jump_row = row_storage + plan->num_train;
  if(provider != NULL) {
    if(provider->function == NULL) {
      status = NP_CONTINUOUS_ROW_ERR_LAYOUT;
      goto cleanup;
    }
  }

  for(evaluation = 0; evaluation < plan->num_eval; ++evaluation) {
    const int omitted_observation = leave_one_out ?
      evaluation + leave_one_out_offset : -1;
    int observation;

    for(size_t output = 0; output < sum_extent; ++output) {
      accumulator->regular_positive_log[output] = -INFINITY;
      accumulator->regular_negative_log[output] = -INFINITY;
      accumulator->jump_positive_log[output] = -INFINITY;
      accumulator->jump_negative_log[output] = -INFINITY;
    }
    status = np_continuous_kernel_beta_derivative_factor_row(
      plan, evaluation, omitted_observation, derivative_coordinate,
      workspace, &row_result);
    if(status != NP_CONTINUOUS_ROW_OK) {
      if(diagnostics != NULL) {
        diagnostics->bad_coordinate = row_result.bad_coordinate;
        diagnostics->bad_observation = row_result.bad_observation;
        diagnostics->beta_status = row_result.beta_status;
      }
      goto cleanup;
    }
    if(provider != NULL) {
      status = provider->function(
        provider->context, evaluation, omitted_observation,
        plan->num_train, factor_log_absolute, factor_sign);
      if(status != NP_CONTINUOUS_ROW_OK)
        goto cleanup;
      for(observation = 0; observation < plan->num_train; ++observation) {
        const double factor_log = factor_log_absolute[observation];
        const int sign = factor_sign[observation];

        if(observation == omitted_observation)
          continue;
        if(ISNAN(factor_log) || factor_log == INFINITY ||
           (sign != -1 && sign != 0 && sign != 1) ||
           ((sign == 0) != (factor_log == -INFINITY))) {
          status = NP_CONTINUOUS_ROW_ERR_NUMERIC;
          goto cleanup;
        }
        if(sign == 0) {
          workspace->primary_log_absolute[observation] = -INFINITY;
          workspace->secondary_log_absolute[observation] = -INFINITY;
          workspace->primary_sign[observation] = 0;
          workspace->secondary_sign[observation] = 0;
        } else {
          if(workspace->primary_sign[observation] != 0) {
            workspace->primary_log_absolute[observation] += factor_log;
            workspace->primary_sign[observation] = (signed char)(
              workspace->primary_sign[observation] * sign);
          }
          if(workspace->secondary_sign[observation] != 0) {
            workspace->secondary_log_absolute[observation] += factor_log;
            workspace->secondary_sign[observation] = (signed char)(
              workspace->secondary_sign[observation] * sign);
          }
          if(ISNAN(workspace->primary_log_absolute[observation]) ||
             workspace->primary_log_absolute[observation] == INFINITY ||
             ISNAN(workspace->secondary_log_absolute[observation]) ||
             workspace->secondary_log_absolute[observation] == INFINITY) {
            status = NP_CONTINUOUS_ROW_ERR_NUMERIC;
            goto cleanup;
          }
        }
      }
    }

    for(observation = 0; observation < plan->num_train; ++observation) {
      const double regular_log = workspace->primary_log_absolute[observation];
      const double jump_log = workspace->secondary_log_absolute[observation];
      const int regular_sign = workspace->primary_sign[observation];
      const int jump_sign = workspace->secondary_sign[observation];
      double powered_log = -INFINITY;
      int powered_sign = 0;
      int powered_jump = 0;

      if(kernel_weights != NULL) {
        if(jump_sign != 0)
          kernel_weights[(size_t)evaluation * (size_t)plan->num_train +
                         (size_t)observation] =
            jump_sign > 0 ? INFINITY : -INFINITY;
        else if(regular_sign == 0)
          kernel_weights[(size_t)evaluation * (size_t)plan->num_train +
                         (size_t)observation] = 0.0;
        else if(regular_log > log(DBL_MAX))
          kernel_weights[(size_t)evaluation * (size_t)plan->num_train +
                         (size_t)observation] =
            regular_sign > 0 ? INFINITY : -INFINITY;
        else
          kernel_weights[(size_t)evaluation * (size_t)plan->num_train +
                         (size_t)observation] =
            regular_sign > 0 ? exp(regular_log) : -exp(regular_log);
      }
      if(observation == omitted_observation)
        continue;

      if(jump_sign != 0) {
        if(kernel_power > 0) {
          status = np_continuous_kernel_signed_log_power(
            jump_log, jump_sign, kernel_power,
            &powered_log, &powered_sign);
          powered_jump = 1;
        } else if(kernel_power == 0) {
          powered_log = 0.0;
          powered_sign = 1;
        } else {
          continue;
        }
      } else {
        status = np_continuous_kernel_signed_log_power(
          regular_log, regular_sign, kernel_power,
          &powered_log, &powered_sign);
      }
      if(status != NP_CONTINUOUS_ROW_OK)
        goto cleanup;
      if(powered_sign == 0)
        continue;

      for(int response_column = 0; response_column < response_extent;
          ++response_column) {
        const double response_value = response_columns > 0 ?
          response[response_column][observation] : 1.0;

        for(int weight_column = 0; weight_column < weight_extent;
            ++weight_column) {
          const double weight_value = weight_columns > 0 ?
            case_weights[weight_column][observation] : 1.0;
          const size_t output =
            (size_t)response_column * (size_t)weight_extent +
            (size_t)weight_column;
          double term_log;
          int term_sign;

          if(!R_FINITE(response_value) || !R_FINITE(weight_value)) {
            status = NP_CONTINUOUS_ROW_ERR_NUMERIC;
            goto cleanup;
          }
          if(response_value == 0.0 || weight_value == 0.0)
            continue;
          term_log = powered_log + log(fabs(response_value)) +
            log(fabs(weight_value));
          term_sign = powered_sign *
            (response_value > 0.0 ? 1 : -1) *
            (weight_value > 0.0 ? 1 : -1);
          if(powered_jump)
            np_continuous_kernel_signed_log_add(
              term_log, term_sign,
              &accumulator->jump_positive_log[output],
              &accumulator->jump_negative_log[output]);
          else
            np_continuous_kernel_signed_log_add(
              term_log, term_sign,
              &accumulator->regular_positive_log[output],
              &accumulator->regular_negative_log[output]);
        }
      }
    }

    for(size_t output = 0; output < sum_extent; ++output) {
      int undefined = 0;
      const double jump = np_continuous_kernel_signed_log_value(
        accumulator->jump_positive_log[output],
        accumulator->jump_negative_log[output], &undefined);
      double * const destination = weighted_sum +
        (size_t)evaluation * sum_extent + output;

      if(undefined) {
        *destination = NA_REAL;
        if(diagnostics != NULL)
          ++diagnostics->undefined_count;
      } else if(jump != 0.0) {
        *destination = jump > 0.0 ? INFINITY : -INFINITY;
      } else {
        *destination = np_continuous_kernel_signed_log_value(
          accumulator->regular_positive_log[output],
          accumulator->regular_negative_log[output], &undefined);
        if(undefined && diagnostics != NULL)
          ++diagnostics->undefined_count;
      }
    }
    if((evaluation & 255) == 0)
      R_CheckUserInterrupt();
  }
  status = NP_CONTINUOUS_ROW_OK;

cleanup:
  return status;
}

/*
 * Consume power one and, when requested, power two from one complete ordinary
 * beta route.  Retaining the common row scale is also useful to power-one-only
 * LP hat/apply consumers, which must not reconstruct absolute weights.  The
 * two output tensors may bind different response/weight columns, matching the
 * shared dual-power context contract.  Every beta segment is composed in
 * signed-log form; an optional provider supplies the product of non-beta
 * factors.  Derivatives retain their separate regular/jump owner.
 */
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
  NPContinuousKernelProgressFunction progress)
{
  const int compute_power2 = weighted_sum_power2 != NULL;
  const int response_extent = response_columns > 0 ? response_columns : 1;
  const int weight_extent = weight_columns > 0 ? weight_columns : 1;
  const int power2_response_extent = power2_response_columns > 0 ?
    power2_response_columns : 1;
  const int power2_weight_extent = power2_weight_columns > 0 ?
    power2_weight_columns : 1;
  size_t sum_extent;
  size_t power2_sum_extent;
  NPContinuousKernelRowStatus status;
  int evaluation;
  int index;

  if(diagnostics != NULL) {
    diagnostics->bad_coordinate = -1;
    diagnostics->bad_observation = -1;
    diagnostics->undefined_count = 0;
    diagnostics->beta_status = NP_BETA_OK;
  }
  if(plan == NULL || plan->route == NULL || plan->num_train <= 0 ||
     plan->num_eval <= 0 || plan->num_continuous <= 0 ||
     (leave_one_out != 0 && leave_one_out != 1) ||
     leave_one_out_offset < 0 || response_columns < 0 ||
     weight_columns < 0 || power2_response_columns < 0 ||
     power2_weight_columns < 0 ||
     (response_columns > 0 && response == NULL) ||
     (weight_columns > 0 && case_weights == NULL) ||
     (power2_response_columns > 0 && power2_response == NULL) ||
     (power2_weight_columns > 0 && power2_case_weights == NULL) ||
     workspace == NULL || row_result == NULL || weighted_sum == NULL ||
     (!compute_power2 &&
      (power2_response != NULL || power2_response_columns != 0 ||
       power2_case_weights != NULL || power2_weight_columns != 0)) ||
     (retain_common_scale != 0 && retain_common_scale != 1) ||
     (!retain_common_scale && scaled_kernel_weights != NULL) ||
     (retain_common_scale && row_result->row == NULL) ||
     (size_t)response_extent > SIZE_MAX / (size_t)weight_extent ||
     (compute_power2 && (size_t)power2_response_extent >
       SIZE_MAX / (size_t)power2_weight_extent))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(index = 0; index < response_columns; ++index)
    if(response[index] == NULL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(index = 0; index < weight_columns; ++index)
    if(case_weights[index] == NULL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(compute_power2) {
    for(index = 0; index < power2_response_columns; ++index)
      if(power2_response[index] == NULL)
        return NP_CONTINUOUS_ROW_ERR_LAYOUT;
    for(index = 0; index < power2_weight_columns; ++index)
      if(power2_case_weights[index] == NULL)
        return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  }
  for(index = 0; index < plan->num_continuous; ++index)
    if(plan->operator[index] != OP_NORMAL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;

  sum_extent = (size_t)response_extent * (size_t)weight_extent;
  power2_sum_extent = compute_power2 ?
    (size_t)power2_response_extent * (size_t)power2_weight_extent : 0;
  if((size_t)plan->num_eval > SIZE_MAX / sum_extent ||
     (compute_power2 &&
      (size_t)plan->num_eval > SIZE_MAX / power2_sum_extent) ||
     (leave_one_out &&
      (plan->num_eval > plan->num_train ||
       leave_one_out_offset > plan->num_train - plan->num_eval)))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;

  for(evaluation = 0; evaluation < plan->num_eval; ++evaluation) {
    const int omitted_observation = leave_one_out ?
      evaluation + leave_one_out_offset : -1;
    const size_t output_base = (size_t)evaluation * sum_extent;
    const size_t power2_output_base =
      (size_t)evaluation * power2_sum_extent;
    int observation;

    status = retain_common_scale ?
      np_continuous_kernel_beta_factor_row_with_log_factor(
        plan, evaluation, omitted_observation, provider,
        workspace, row_result) :
      np_continuous_kernel_beta_log_factor_row(
        plan, evaluation, omitted_observation, provider,
        workspace, row_result);
    if(status != NP_CONTINUOUS_ROW_OK) {
      if(diagnostics != NULL) {
        diagnostics->bad_coordinate = row_result->bad_coordinate;
        diagnostics->bad_observation = row_result->bad_observation;
        diagnostics->beta_status = row_result->beta_status;
      }
      return status;
    }
    for(size_t output = 0; output < sum_extent; ++output)
      weighted_sum[output_base + output] = 0.0;
    if(compute_power2)
      for(size_t output = 0; output < power2_sum_extent; ++output)
        weighted_sum_power2[power2_output_base + output] = 0.0;

    for(observation = 0; observation < plan->num_train; ++observation) {
      double value;
      double value_power2;

      if(retain_common_scale) {
        value = row_result->row[observation];
        status = R_FINITE(value) ? NP_CONTINUOUS_ROW_OK :
          NP_CONTINUOUS_ROW_ERR_NUMERIC;
      } else {
        status = np_continuous_kernel_signed_log_restore(
          workspace->primary_log_absolute[observation],
          workspace->primary_sign[observation], &value);
      }
      if(compute_power2 && status == NP_CONTINUOUS_ROW_OK) {
        /* A dual consumer already owns K in representable arithmetic.  Square
         * it directly instead of paying for a second exp(2 * log(K)). */
        value_power2 = value * value;
        if(!R_FINITE(value_power2))
          status = NP_CONTINUOUS_ROW_ERR_NUMERIC;
      }
      if(status != NP_CONTINUOUS_ROW_OK)
        return status;
      if(scaled_kernel_weights != NULL)
        scaled_kernel_weights[
          (size_t)evaluation * (size_t)plan->num_train +
          (size_t)observation
        ] = value;
      if(observation == omitted_observation)
        continue;

      for(int response_column = 0; response_column < response_extent;
          ++response_column) {
        const double response_value = response_columns > 0 ?
          response[response_column][observation] : 1.0;

        for(int weight_column = 0; weight_column < weight_extent;
            ++weight_column) {
          const double weight_value = weight_columns > 0 ?
            case_weights[weight_column][observation] : 1.0;
          const size_t output = output_base +
            (size_t)response_column * (size_t)weight_extent +
            (size_t)weight_column;

          if(!R_FINITE(response_value) || !R_FINITE(weight_value))
            return NP_CONTINUOUS_ROW_ERR_NUMERIC;
          weighted_sum[output] += value * response_value * weight_value;
          if(!R_FINITE(weighted_sum[output]))
            return NP_CONTINUOUS_ROW_ERR_NUMERIC;
        }
      }
      if(compute_power2) {
        for(int response_column = 0;
            response_column < power2_response_extent; ++response_column) {
          const double response_value = power2_response_columns > 0 ?
            power2_response[response_column][observation] : 1.0;

          for(int weight_column = 0;
              weight_column < power2_weight_extent; ++weight_column) {
            const double weight_value = power2_weight_columns > 0 ?
              power2_case_weights[weight_column][observation] : 1.0;
            const size_t output = power2_output_base +
              (size_t)response_column * (size_t)power2_weight_extent +
              (size_t)weight_column;

            if(!R_FINITE(response_value) || !R_FINITE(weight_value))
              return NP_CONTINUOUS_ROW_ERR_NUMERIC;
            weighted_sum_power2[output] +=
              value_power2 * response_value * weight_value;
            if(!R_FINITE(weighted_sum_power2[output]))
              return NP_CONTINUOUS_ROW_ERR_NUMERIC;
          }
        }
      }
    }
    if(progress != NULL)
      progress(evaluation + 1, plan->num_eval);
    if((evaluation & 31) == 0)
      R_CheckUserInterrupt();
  }
  return NP_CONTINUOUS_ROW_OK;
}

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
  NPContinuousKernelProgressFunction progress)
{
  NPContinuousKernelRowStatus status;
  int evaluation;
  int coordinate;

  if(diagnostics != NULL) {
    diagnostics->bad_coordinate = -1;
    diagnostics->bad_observation = -1;
    diagnostics->undefined_count = 0;
    diagnostics->beta_status = NP_BETA_OK;
  }
  if(plan == NULL || plan->route == NULL || plan->num_train <= 0 ||
     plan->num_eval <= 0 || plan->num_continuous <= 0 ||
     (leave_one_out != 0 && leave_one_out != 1) ||
     leave_one_out_offset < 0 || plan->operator == NULL ||
     workspace == NULL || row_result == NULL || sum == NULL ||
     centered_m2 == NULL ||
     (leave_one_out &&
      (plan->num_eval > plan->num_train ||
       leave_one_out_offset > plan->num_train - plan->num_eval)))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(coordinate = 0; coordinate < plan->num_continuous; ++coordinate)
    if(plan->operator[coordinate] != OP_NORMAL &&
       plan->operator[coordinate] != OP_INTEGRAL &&
       plan->operator[coordinate] != OP_CONVOLUTION)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;

  for(evaluation = 0; evaluation < plan->num_eval; ++evaluation) {
    const int omitted_observation = leave_one_out ?
      evaluation + leave_one_out_offset : -1;
    double running_mean = 0.0;
    double running_m2 = 0.0;
    int running_count = 0;
    int observation;

    sum[evaluation] = 0.0;
    centered_m2[evaluation] = 0.0;
    status = np_continuous_kernel_beta_log_factor_row(
      plan, evaluation, omitted_observation, provider,
      workspace, row_result);
    if(status != NP_CONTINUOUS_ROW_OK) {
      if(diagnostics != NULL) {
        diagnostics->bad_coordinate = row_result->bad_coordinate;
        diagnostics->bad_observation = row_result->bad_observation;
        diagnostics->beta_status = row_result->beta_status;
      }
      return status;
    }

    for(observation = 0; observation < plan->num_train; ++observation) {
      double value;
      double delta;

      status = np_continuous_kernel_signed_log_restore(
        workspace->primary_log_absolute[observation],
        workspace->primary_sign[observation], &value);
      if(status != NP_CONTINUOUS_ROW_OK) {
        if(diagnostics != NULL)
          diagnostics->bad_observation = observation;
        return status;
      }
      if(observation == omitted_observation)
        continue;

      sum[evaluation] += value;
      delta = value - running_mean;
      ++running_count;
      running_mean += delta / (double)running_count;
      running_m2 += delta * (value - running_mean);
      if(!R_FINITE(sum[evaluation]) || !R_FINITE(running_mean) ||
         !R_FINITE(running_m2)) {
        if(diagnostics != NULL)
          diagnostics->bad_observation = observation;
        return NP_CONTINUOUS_ROW_ERR_NUMERIC;
      }
    }
    if(running_m2 < 0.0)
      running_m2 = 0.0;
    centered_m2[evaluation] = running_m2;
    if(progress != NULL)
      progress(evaluation + 1, plan->num_eval);
    if((evaluation & 31) == 0)
      R_CheckUserInterrupt();
  }
  return NP_CONTINUOUS_ROW_OK;
}

NPContinuousKernelRowStatus np_continuous_kernel_scaled_restore(
  double scaled_value,
  double log_scale,
  int power,
  double *value)
{
  double log_absolute;
  const int sign = scaled_value < 0.0 ? -1 : 1;

  if(value == NULL || !R_FINITE(scaled_value) ||
     ISNAN(log_scale) || log_scale == INFINITY)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(scaled_value == 0.0 || log_scale == -INFINITY) {
    *value = 0.0;
    return NP_CONTINUOUS_ROW_OK;
  }
  if(power == 1 && log_scale == 0.0) {
    *value = scaled_value;
    return NP_CONTINUOUS_ROW_OK;
  }
  log_absolute = log(fabs(scaled_value)) + log_scale;
  if(ISNAN(log_absolute) || log_absolute == INFINITY)
    return NP_CONTINUOUS_ROW_ERR_NUMERIC;
  return np_continuous_kernel_signed_log_power_restore(
    log_absolute, sign, power, value);
}

NPContinuousKernelRowStatus np_continuous_kernel_scaled_derivative_restore(
  double scaled_value,
  double log_scale,
  double *value)
{
  if(value == NULL || ISNAN(log_scale) || log_scale == INFINITY)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(ISNA(scaled_value)) {
    *value = NA_REAL;
    return NP_CONTINUOUS_ROW_OK;
  }
  if(ISNAN(scaled_value))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(R_FINITE(scaled_value))
    return np_continuous_kernel_scaled_restore(
      scaled_value, log_scale, 1, value);
  if(log_scale == -INFINITY)
    return NP_CONTINUOUS_ROW_ERR_NUMERIC;
  *value = scaled_value;
  return NP_CONTINUOUS_ROW_OK;
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

NPContinuousKernelRowStatus np_continuous_kernel_signed_log_power_restore(
  double log_absolute,
  int sign,
  int power,
  double *value)
{
  double powered_log_absolute;
  int powered_sign;
  NPContinuousKernelRowStatus status;

  if(value == NULL)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  status = np_continuous_kernel_signed_log_power(
    log_absolute, sign, power, &powered_log_absolute, &powered_sign);
  if(status != NP_CONTINUOUS_ROW_OK)
    return status;
  return np_continuous_kernel_signed_log_restore(
    powered_log_absolute, powered_sign, value);
}

NPContinuousKernelRowStatus np_continuous_kernel_signed_log_power(
  double log_absolute,
  int sign,
  int power,
  double *powered_log_absolute,
  int *powered_sign)
{
  if(powered_log_absolute == NULL || powered_sign == NULL ||
     ISNAN(log_absolute) || log_absolute == INFINITY ||
     (sign != -1 && sign != 0 && sign != 1) ||
     ((sign == 0) != (log_absolute == -INFINITY)))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(sign == 0) {
    *powered_log_absolute = -INFINITY;
    *powered_sign = 0;
    return NP_CONTINUOUS_ROW_OK;
  }
  if(power == 0) {
    *powered_log_absolute = 0.0;
    *powered_sign = 1;
    return NP_CONTINUOUS_ROW_OK;
  }

  *powered_log_absolute = (double)power * log_absolute;
  if(ISNAN(*powered_log_absolute) || *powered_log_absolute == INFINITY)
    return NP_CONTINUOUS_ROW_ERR_NUMERIC;
  if(*powered_log_absolute == -INFINITY) {
    *powered_sign = 0;
    return NP_CONTINUOUS_ROW_OK;
  }
  *powered_sign = ((power % 2) == 0) ? 1 : sign;
  return NP_CONTINUOUS_ROW_OK;
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
  case NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT:
    return "continuous-kernel row has zero total weight";
  default:
    return "unknown continuous-kernel row status";
  }
}

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
  NPContinuousKernelProgressFunction progress)
{
  NPContinuousKernelRowStatus status;
  int evaluation;
  int coordinate;
  int observation;

  if(diagnostics != NULL) {
    diagnostics->bad_coordinate = -1;
    diagnostics->bad_observation = -1;
    diagnostics->undefined_count = 0;
    diagnostics->beta_status = NP_BETA_OK;
  }
  if(plan == NULL || plan->route == NULL || plan->num_train <= 0 ||
     plan->num_eval <= 0 || plan->num_continuous <= 0 ||
     plan->route->segment_count != 1 ||
     plan->route->segment[0].descriptor.family != NP_CKERNEL_FAMILY_BETA ||
     (leave_one_out != 0 && leave_one_out != 1) ||
     leave_one_out_offset < 0 ||
     (positive_weights != 0 && positive_weights != 1) ||
     plan->operator == NULL || response == NULL || workspace == NULL ||
     row_result == NULL || mean == NULL || mean_stderr == NULL ||
     (leave_one_out &&
      (plan->num_eval > plan->num_train ||
       leave_one_out_offset > plan->num_train - plan->num_eval)))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(coordinate = 0; coordinate < plan->num_continuous; ++coordinate)
    if(plan->operator[coordinate] != OP_NORMAL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(observation = 0; observation < plan->num_train; ++observation)
    if(!R_FINITE(response[observation])) {
      if(diagnostics != NULL)
        diagnostics->bad_observation = observation;
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
    }

  for(evaluation = 0; evaluation < plan->num_eval; ++evaluation) {
    const int omitted_observation = leave_one_out ?
      evaluation + leave_one_out_offset : -1;
    double total_weight = 0.0;
    double weighted_mean = 0.0;
    double weighted_m2 = 0.0;
    double squared_weight_sum = 0.0;

    status = np_continuous_kernel_beta_log_factor_row(
      plan, evaluation, omitted_observation, provider,
      workspace, row_result);
    if(status != NP_CONTINUOUS_ROW_OK) {
      if(diagnostics != NULL) {
        diagnostics->bad_coordinate = row_result->bad_coordinate;
        diagnostics->bad_observation = row_result->bad_observation;
        diagnostics->beta_status = row_result->beta_status;
      }
      return status;
    }
    if(row_result->total_log_scale == -INFINITY)
      return NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT;

    if(positive_weights) {
      for(observation = 0; observation < plan->num_train; ++observation) {
        const double weight = workspace->primary_sign[observation] == 0 ?
          0.0 : (double)workspace->primary_sign[observation] * exp(
            workspace->primary_log_absolute[observation] -
            row_result->total_log_scale);

        if(observation == omitted_observation)
          continue;
        if(weight < 0.0 || !R_FINITE(weight)) {
          if(diagnostics != NULL)
            diagnostics->bad_observation = observation;
          return NP_CONTINUOUS_ROW_ERR_NUMERIC;
        }
        if(weight > 0.0) {
          const double new_total_weight = total_weight + weight;
          const double delta = response[observation] - weighted_mean;
          const double new_mean = weighted_mean +
            (weight / new_total_weight) * delta;

          weighted_m2 += weight * delta *
            (response[observation] - new_mean);
          squared_weight_sum += weight * weight;
          total_weight = new_total_weight;
          weighted_mean = new_mean;
        }
      }
    } else {
      double weighted_response_sum = 0.0;

      for(observation = 0; observation < plan->num_train; ++observation) {
        const double weight = workspace->primary_sign[observation] == 0 ?
          0.0 : (double)workspace->primary_sign[observation] * exp(
            workspace->primary_log_absolute[observation] -
            row_result->total_log_scale);

        if(observation == omitted_observation)
          continue;
        total_weight += weight;
        weighted_response_sum += weight * response[observation];
      }
      if(!R_FINITE(total_weight) || total_weight == 0.0 ||
         !R_FINITE(weighted_response_sum))
        return total_weight == 0.0 ? NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT :
          NP_CONTINUOUS_ROW_ERR_NUMERIC;
      weighted_mean = weighted_response_sum / total_weight;

      for(observation = 0; observation < plan->num_train; ++observation) {
        const double weight = workspace->primary_sign[observation] == 0 ?
          0.0 : (double)workspace->primary_sign[observation] * exp(
            workspace->primary_log_absolute[observation] -
            row_result->total_log_scale);
        const double residual = response[observation] - weighted_mean;

        if(observation == omitted_observation)
          continue;
        weighted_m2 += weight * residual * residual;
        squared_weight_sum += weight * weight;
      }
    }

    if(!R_FINITE(total_weight) ||
       (positive_weights ? total_weight <= 0.0 : total_weight == 0.0))
      return total_weight == 0.0 ? NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT :
        NP_CONTINUOUS_ROW_ERR_NUMERIC;
    if(!R_FINITE(weighted_mean) || !R_FINITE(weighted_m2) ||
       !R_FINITE(squared_weight_sum))
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
    if((positive_weights && weighted_m2 < 0.0) ||
       (!positive_weights && weighted_m2 / total_weight < 0.0))
      weighted_m2 = 0.0;

    mean[evaluation] = weighted_mean;
    mean_stderr[evaluation] = sqrt(
      (weighted_m2 / total_weight) *
      (squared_weight_sum / (total_weight * total_weight)));
    if(!R_FINITE(mean_stderr[evaluation]))
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
    if(progress != NULL)
      progress(evaluation + 1, plan->num_eval);
    if((evaluation & 31) == 0)
      R_CheckUserInterrupt();
  }
  return NP_CONTINUOUS_ROW_OK;
}

/*
 * Conditional scalar estimators use an observation-influence variance that
 * is intentionally distinct from ordinary regression's local-residual
 * covariance.  Keep it out of line so adding this policy cannot enlarge or
 * perturb the ordinary regression moment hot path.
 */
static NPContinuousKernelRowStatus NP_CONTINUOUS_ROW_NOINLINE
np_continuous_kernel_beta_conditional_influence_stderr(
  const NPContinuousKernelRowWorkspace *workspace,
  const double *response,
  const int num_train,
  const int omitted_observation,
  const double total_log_scale,
  const double weighted_mean,
  const double total_weight,
  double *mean_stderr)
{
  const int variance_count =
    num_train - (omitted_observation >= 0 ? 1 : 0);
  double influence_scale = 0.0;
  double influence_sum_squares = 1.0;
  int observation;

  if(workspace == NULL || response == NULL || mean_stderr == NULL ||
     num_train <= 0 || omitted_observation < -1 ||
     omitted_observation >= num_train || ISNAN(total_log_scale) ||
     total_log_scale == INFINITY || !R_FINITE(weighted_mean) ||
     !R_FINITE(total_weight) || total_weight == 0.0)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(variance_count <= 1) {
    *mean_stderr = 0.0;
    return NP_CONTINUOUS_ROW_OK;
  }

  for(observation = 0; observation < num_train; ++observation) {
    const double weight = workspace->primary_sign[observation] == 0 ?
      0.0 : (double)workspace->primary_sign[observation] * exp(
        workspace->primary_log_absolute[observation] - total_log_scale);
    const double influence =
      weight * (response[observation] - weighted_mean);
    const double absolute_influence = fabs(influence);

    if(observation == omitted_observation)
      continue;
    if(!R_FINITE(influence))
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
    if(absolute_influence != 0.0) {
      if(influence_scale < absolute_influence) {
        const double ratio = influence_scale / absolute_influence;

        influence_sum_squares =
          1.0 + influence_sum_squares * ratio * ratio;
        influence_scale = absolute_influence;
      } else {
        const double ratio = absolute_influence / influence_scale;

        influence_sum_squares += ratio * ratio;
      }
    }
  }
  *mean_stderr = influence_scale * sqrt(influence_sum_squares) /
    (fabs(total_weight) * sqrt((double)(variance_count - 1)));
  return R_FINITE(*mean_stderr) ? NP_CONTINUOUS_ROW_OK :
    NP_CONTINUOUS_ROW_ERR_NUMERIC;
}

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
  NPContinuousKernelProgressFunction progress)
{
  NPContinuousKernelRowStatus status;
  int evaluation;
  int coordinate;
  int observation;

  if(diagnostics != NULL) {
    diagnostics->bad_coordinate = -1;
    diagnostics->bad_observation = -1;
    diagnostics->undefined_count = 0;
    diagnostics->beta_status = NP_BETA_OK;
  }
  if(plan == NULL || plan->route == NULL || plan->num_train <= 0 ||
     plan->num_eval <= 0 || plan->num_continuous <= 0 ||
     plan->route->segment_count != 1 ||
     plan->route->segment[0].descriptor.family != NP_CKERNEL_FAMILY_BETA ||
     (leave_one_out != 0 && leave_one_out != 1) ||
     leave_one_out_offset < 0 ||
     (positive_weights != 0 && positive_weights != 1) ||
     plan->operator == NULL || response == NULL || workspace == NULL ||
     row_result == NULL || mean == NULL || mean_stderr == NULL ||
     (leave_one_out &&
      (plan->num_eval > plan->num_train ||
       leave_one_out_offset > plan->num_train - plan->num_eval)))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(coordinate = 0; coordinate < plan->num_continuous; ++coordinate)
    if(plan->operator[coordinate] != OP_NORMAL)
      return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  for(observation = 0; observation < plan->num_train; ++observation)
    if(!R_FINITE(response[observation])) {
      if(diagnostics != NULL)
        diagnostics->bad_observation = observation;
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
    }

  for(evaluation = 0; evaluation < plan->num_eval; ++evaluation) {
    const int omitted_observation = leave_one_out ?
      evaluation + leave_one_out_offset : -1;
    double total_weight = 0.0;
    double weighted_mean = 0.0;
    double weighted_m2 = 0.0;
    double squared_weight_sum = 0.0;

    status = np_continuous_kernel_beta_log_factor_row(
      plan, evaluation, omitted_observation, provider,
      workspace, row_result);
    if(status != NP_CONTINUOUS_ROW_OK) {
      if(diagnostics != NULL) {
        diagnostics->bad_coordinate = row_result->bad_coordinate;
        diagnostics->bad_observation = row_result->bad_observation;
        diagnostics->beta_status = row_result->beta_status;
      }
      return status;
    }
    if(row_result->total_log_scale == -INFINITY)
      return NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT;

    if(positive_weights) {
      for(observation = 0; observation < plan->num_train; ++observation) {
        const double weight = workspace->primary_sign[observation] == 0 ?
          0.0 : (double)workspace->primary_sign[observation] * exp(
            workspace->primary_log_absolute[observation] -
            row_result->total_log_scale);

        if(observation == omitted_observation)
          continue;
        if(weight < 0.0 || !R_FINITE(weight)) {
          if(diagnostics != NULL)
            diagnostics->bad_observation = observation;
          return NP_CONTINUOUS_ROW_ERR_NUMERIC;
        }
        if(weight > 0.0) {
          const double new_total_weight = total_weight + weight;
          const double delta = response[observation] - weighted_mean;
          const double new_mean = weighted_mean +
            (weight / new_total_weight) * delta;

          weighted_m2 += weight * delta *
            (response[observation] - new_mean);
          squared_weight_sum += weight * weight;
          total_weight = new_total_weight;
          weighted_mean = new_mean;
        }
      }
    } else {
      double weighted_response_sum = 0.0;

      for(observation = 0; observation < plan->num_train; ++observation) {
        const double weight = workspace->primary_sign[observation] == 0 ?
          0.0 : (double)workspace->primary_sign[observation] * exp(
            workspace->primary_log_absolute[observation] -
            row_result->total_log_scale);

        if(observation == omitted_observation)
          continue;
        total_weight += weight;
        weighted_response_sum += weight * response[observation];
      }
      if(!R_FINITE(total_weight) || total_weight == 0.0 ||
         !R_FINITE(weighted_response_sum))
        return total_weight == 0.0 ? NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT :
          NP_CONTINUOUS_ROW_ERR_NUMERIC;
      weighted_mean = weighted_response_sum / total_weight;

      for(observation = 0; observation < plan->num_train; ++observation) {
        const double weight = workspace->primary_sign[observation] == 0 ?
          0.0 : (double)workspace->primary_sign[observation] * exp(
            workspace->primary_log_absolute[observation] -
            row_result->total_log_scale);
        const double residual = response[observation] - weighted_mean;

        if(observation == omitted_observation)
          continue;
        weighted_m2 += weight * residual * residual;
        squared_weight_sum += weight * weight;
      }
    }

    if(!R_FINITE(total_weight) ||
       (positive_weights ? total_weight <= 0.0 : total_weight == 0.0))
      return total_weight == 0.0 ? NP_CONTINUOUS_ROW_ERR_ZERO_WEIGHT :
        NP_CONTINUOUS_ROW_ERR_NUMERIC;
    if(!R_FINITE(weighted_mean) || !R_FINITE(weighted_m2) ||
       !R_FINITE(squared_weight_sum))
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
    if((positive_weights && weighted_m2 < 0.0) ||
       (!positive_weights && weighted_m2 / total_weight < 0.0))
      weighted_m2 = 0.0;

    mean[evaluation] = weighted_mean;
    status = np_continuous_kernel_beta_conditional_influence_stderr(
      workspace, response, plan->num_train, omitted_observation,
      row_result->total_log_scale, weighted_mean, total_weight,
      &mean_stderr[evaluation]);
    if(status != NP_CONTINUOUS_ROW_OK)
      return status;
    if(!R_FINITE(mean_stderr[evaluation]))
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
    if(progress != NULL)
      progress(evaluation + 1, plan->num_eval);
    if((evaluation & 31) == 0)
      R_CheckUserInterrupt();
  }
  return NP_CONTINUOUS_ROW_OK;
}
