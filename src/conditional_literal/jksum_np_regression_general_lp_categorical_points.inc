static NP_ALWAYS_INLINE int np_regression_general_lp_categorical_points(
  const int prepared_ann,
  const NPRegressionGeneralLPFitCall *call,
  NPRegressionGeneralLPFitOwner *owner,
  const int moment_stride,
  const int response_y_offset,
  const int response_basis_offset,
  const double epsilon,
  const double fitted_point,
  const int row,
  const int compute_hc0)
{
  double base_point = fitted_point;
  int coordinate;

  if(call == NULL || owner == NULL || owner->categorical_point == NULL)
    return 0;
  const int compute_pair = call->ordinary_response;

  if(call->categorical_point_invariant) {
    for(coordinate = call->categorical_range.begin - call->num_reg_continuous;
        coordinate < call->categorical_range.end - call->num_reg_continuous;
        ++coordinate)
      owner->categorical_point[coordinate] = 0.0;
    return 1;
  }

  if(call->categorical_base_requires_refit &&
     !np_regression_general_lp_point_at_frame(
       prepared_ann, call, owner, moment_stride, response_y_offset,
       response_basis_offset, epsilon, row, -1, &base_point,
       compute_pair ? owner->categorical_base_kernel_row : NULL,
       compute_pair ? owner->power2_projection : NULL))
    return 0;

  for(coordinate = call->unordered_range.begin - call->num_reg_continuous;
      coordinate < call->unordered_range.end - call->num_reg_continuous;
      ++coordinate) {
    const double current = owner->eval_unordered[coordinate][0];
    const double alternate = call->matrix_categorical_vals[coordinate][0];
    double alternate_point;

    if(current == alternate) {
      owner->categorical_point[coordinate] = 0.0;
      if(compute_hc0) {
        owner->categorical_stderr[coordinate] = 0.0;
        if(call->hc0_context->gradient_structural_zero != NULL)
          call->hc0_context->gradient_structural_zero[
            call->num_reg_continuous + coordinate][row] = 1;
      }
      continue;
    }
    owner->eval_unordered[coordinate][0] = alternate;
    if(!np_regression_general_lp_point_at_frame(
         prepared_ann, call, owner, moment_stride, response_y_offset,
         response_basis_offset, epsilon, row,
         call->num_reg_continuous + coordinate, &alternate_point,
         compute_pair ? owner->categorical_alternate_kernel_row : NULL,
         compute_pair ? owner->coefficient : NULL)) {
      owner->eval_unordered[coordinate][0] = current;
      return 0;
    }
    owner->eval_unordered[coordinate][0] = current;
    owner->categorical_point[coordinate] = base_point - alternate_point;
    if(ISNA(alternate_point)) {
      if(compute_hc0) owner->categorical_stderr[coordinate] = NA_REAL;
      continue;
    }
    if(compute_pair &&
       !np_regression_lp_categorical_contrast(
         owner->basis,
         owner->nterms,
         call->num_obs_train,
         owner->power2_projection,
         owner->coefficient,
         owner->categorical_base_kernel_row,
         owner->categorical_alternate_kernel_row,
         call->vector_Y, &owner->categorical_point[coordinate],
         compute_hc0 ? &owner->categorical_stderr[coordinate] : NULL,
         compute_hc0 ? call->hc0_context : NULL,
         compute_hc0 && call->hc0_context->gradient_unavailable != NULL ?
           &call->hc0_context->gradient_unavailable[
             call->num_reg_continuous + coordinate][row] : NULL,
         compute_hc0 && call->hc0_context->gradient_structural_zero != NULL ?
           &call->hc0_context->gradient_structural_zero[
             call->num_reg_continuous + coordinate][row] : NULL))
      return 0;
  }

  for(coordinate = call->ordered_range.begin - call->num_reg_continuous - call->num_reg_unordered;
      coordinate < call->ordered_range.end - call->num_reg_continuous - call->num_reg_unordered;
      ++coordinate) {
    const int category = call->num_reg_unordered + coordinate;
    const int count = call->num_categories[category];
    const double current = owner->eval_ordered[coordinate][0];
    double alternate_point;
    double alternate;
    double sign;
    int level = -1;
    int index;

    if(count <= 1) {
      owner->categorical_point[category] = 0.0;
      if(compute_hc0) {
        owner->categorical_stderr[category] = 0.0;
        if(call->hc0_context->gradient_structural_zero != NULL)
          call->hc0_context->gradient_structural_zero[
            call->num_reg_continuous + category][row] = 1;
      }
      continue;
    }
    for(index = 0; index < count; ++index)
      if(call->matrix_categorical_vals[category][index] == current) {
        level = index;
        break;
      }
    if(level < 0)
      return 0;
    if(level == 0) {
      alternate = call->matrix_categorical_vals[category][1];
      sign = -1.0;
    } else {
      alternate = call->matrix_categorical_vals[category][level - 1];
      sign = 1.0;
    }
    owner->eval_ordered[coordinate][0] = alternate;
    if(!np_regression_general_lp_point_at_frame(
         prepared_ann, call, owner, moment_stride, response_y_offset,
         response_basis_offset, epsilon, row,
         call->num_reg_continuous + category, &alternate_point,
         compute_pair ? owner->categorical_alternate_kernel_row : NULL,
         compute_pair ? owner->coefficient : NULL)) {
      owner->eval_ordered[coordinate][0] = current;
      return 0;
    }
    owner->eval_ordered[coordinate][0] = current;
    owner->categorical_point[category] =
      sign*(base_point - alternate_point);
    if(ISNA(alternate_point)) {
      if(compute_hc0) owner->categorical_stderr[category] = NA_REAL;
      continue;
    }
    if(compute_pair &&
       !np_regression_lp_categorical_contrast(
         owner->basis,
         owner->nterms,
         call->num_obs_train,
         owner->power2_projection,
         owner->coefficient,
         owner->categorical_base_kernel_row,
         owner->categorical_alternate_kernel_row,
         call->vector_Y, &owner->categorical_point[category],
         compute_hc0 ? &owner->categorical_stderr[category] : NULL,
         compute_hc0 ? call->hc0_context : NULL,
         compute_hc0 && call->hc0_context->gradient_unavailable != NULL ?
           &call->hc0_context->gradient_unavailable[
             call->num_reg_continuous + category][row] : NULL,
         compute_hc0 && call->hc0_context->gradient_structural_zero != NULL ?
           &call->hc0_context->gradient_structural_zero[
             call->num_reg_continuous + category][row] : NULL))
      return 0;
    if(compute_pair) owner->categorical_point[category] *= sign;
  }
  return 1;
}
