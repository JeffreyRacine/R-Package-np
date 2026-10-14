NPCategoricalProfileTileStatus
np_categorical_profile_log_row_fill_prevalidated(
  const NPCategoricalProfileKernelSpec *spec,
  int eval_index,
  int omitted_observation,
  double *scratch,
  size_t scratch_elements,
  double *log_absolute,
  signed char *sign)
{
  NPCategoricalProfileTileStatus status;
  int coordinate;
  int observation;

  if(spec == NULL || scratch == NULL || log_absolute == NULL || sign == NULL)
    return NP_PROFILE_TILE_ERR_ARGUMENT;
  if(eval_index < 0 || eval_index >= spec->neval ||
     omitted_observation < -1 || omitted_observation >= spec->ntrain)
    return NP_PROFILE_TILE_ERR_RANGE;
  if(scratch_elements < (size_t)spec->ntrain)
    return NP_PROFILE_TILE_ERR_CAPACITY;
  status = np_categorical_profile_evaluation_data_validate(
    spec, eval_index, 1);
  if(status != NP_PROFILE_TILE_OK)
    return status;

  for(observation = 0; observation < spec->ntrain; ++observation) {
    log_absolute[observation] =
      observation == omitted_observation ? -INFINITY : 0.0;
    sign[observation] = observation == omitted_observation ? 0 : 1;
  }
  for(coordinate = 0;
      coordinate < spec->nunordered + spec->nordered;
      ++coordinate) {
    np_categorical_profile_coordinate_fill_prevalidated(
      spec, eval_index, coordinate, 0, scratch);
    for(observation = 0; observation < spec->ntrain; ++observation) {
      const double factor = scratch[observation];

      if(observation == omitted_observation)
        continue;
      if(!isfinite(factor))
        return NP_PROFILE_TILE_ERR_NONFINITE;
      if(sign[observation] == 0)
        continue;
      if(factor == 0.0) {
        log_absolute[observation] = -INFINITY;
        sign[observation] = 0;
      } else {
        log_absolute[observation] += log(fabs(factor));
        sign[observation] =
          (signed char)(sign[observation] * (factor < 0.0 ? -1 : 1));
        if(!isfinite(log_absolute[observation]))
          return NP_PROFILE_TILE_ERR_NONFINITE;
      }
    }
  }
  return NP_PROFILE_TILE_OK;
}
