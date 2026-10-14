static NPContinuousKernelRowStatus np_beta_categorical_log_factor(
  const void *state,
  int evaluation_index,
  int omitted_observation,
  int observation_count,
  double *log_absolute,
  signed char *sign)
{
  const NPBetaScaledRowCategoricalContext * const context =
    (const NPBetaScaledRowCategoricalContext *)state;
  NPCategoricalProfileTileStatus tile_status;
  int observation;

  if(context == NULL || observation_count != context->dense_spec.ntrain)
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;
  if(!context->use_compressed) {
    tile_status = np_categorical_profile_log_row_fill_prevalidated(
      &context->dense_spec, evaluation_index, omitted_observation,
      context->scratch, context->scratch_capacity, log_absolute, sign);
    return np_beta_categorical_tile_status(tile_status);
  }

  tile_status = np_categorical_profile_log_row_fill_prevalidated(
    &context->compressed_spec, evaluation_index, -1,
    context->scratch, context->scratch_capacity,
    context->compressed_log_absolute, context->compressed_sign);
  if(tile_status != NP_PROFILE_TILE_OK)
    return np_beta_categorical_tile_status(tile_status);
  for(observation = 0; observation < observation_count; ++observation) {
    const int profile = context->profile_id[observation];

    log_absolute[observation] =
      context->compressed_log_absolute[profile];
    sign[observation] = context->compressed_sign[profile];
  }
  return NP_CONTINUOUS_ROW_OK;
}
