static int np_regression_categorical_profile_fit_body(
  const NPRegressionCategoricalProfileFitCall *call,
  NPCategoricalProfileOwner *owner)
{
  int * const kernel_c = call->kernel_c;
  int * const kernel_u = call->kernel_u;
  int * const kernel_o = call->kernel_o;
  const int lp_engine = call->lp_engine;
  const int BANDWIDTH_reg = call->bandwidth_mode;
  const int num_obs_train = call->num_obs_train;
  const int num_obs_eval = call->num_obs_eval;
  const int num_reg_unordered = call->num_reg_unordered;
  const int num_reg_ordered = call->num_reg_ordered;
  const int num_reg_continuous = call->num_reg_continuous;
  double ** const matrix_X_unordered_train =
    call->matrix_X_unordered_train;
  double ** const matrix_X_ordered_train = call->matrix_X_ordered_train;
  double ** const matrix_X_unordered_eval = call->matrix_X_unordered_eval;
  double ** const matrix_X_ordered_eval = call->matrix_X_ordered_eval;
  double * const vector_Y = call->vector_Y;
  double * const vector_scale_factor = call->vector_scale_factor;
  double * const lambda = call->lambda;
  int * const num_categories = call->num_categories;
  int * const operator = call->operator;
  double ** const matrix_categorical_vals = call->matrix_categorical_vals;
  double * const mean = call->mean;
  double * const mean_stderr = call->mean_stderr;
  const int do_merr = mean_stderr != NULL;
  const int response_column_count = do_merr ? 3 : 2;
  const int emit_fit_progress = call->emit_fit_progress;

  int i, j, g, status;
  int *train_prof_id = NULL, *train_prof_rep = NULL;
  int *eval_prof_id = NULL, *eval_prof_rep = NULL;
  int nprof_train = 0, nprof_eval = 0;
  double **profile_unordered_train = NULL, **profile_ordered_train = NULL;
  double **profile_unordered_eval = NULL, **profile_ordered_eval = NULL;
  double *counts = NULL, *sums = NULL, *sums2 = NULL, *weighted_sum = NULL;
  double *profile_y[3];
  double *profile_mean = NULL, *profile_stderr = NULL;
  NPRegressionResidualKernelContext residual_context = {0};
  long double *profile_anchor = NULL, *profile_centered_sum = NULL;
  int *profile_eval_to_train = NULL;
  int ok = 0;

  if((int_TREE_PROFILE_X != NP_TREE_TRUE) ||
     (BANDWIDTH_reg != BW_FIXED) ||
     (num_reg_continuous != 0) ||
     ((num_reg_unordered + num_reg_ordered) <= 0) ||
     (vector_Y == NULL) ||
     (mean == NULL) ||
     (lambda == NULL) ||
     (operator == NULL))
    return 0;

  if((lp_engine != NP_LP_ENGINE_SCALAR) && (lp_engine != NP_LP_ENGINE_GENERAL))
    return 0;

  for(i = 0; i < (num_reg_unordered + num_reg_ordered); i++)
    if(operator[i] != OP_NORMAL)
      return 0;

  for(i = 0; i < num_reg_unordered; i++)
    if((lambda[i] > 0.5) ||
       np_disc_near_upper(kernel_u[i], lambda[i], num_categories[i]))
      return 0;

  for(i = 0; i < num_reg_ordered; i++)
    if((lambda[num_reg_unordered + i] > 0.5) ||
       np_disc_ordered_near_upper(kernel_o[i], lambda[num_reg_unordered + i]))
      return 0;

  if(!np_build_discrete_profile_index(num_obs_train,
                                      num_reg_unordered,
                                      num_reg_ordered,
                                      matrix_X_unordered_train,
                                      matrix_X_ordered_train,
                                      &train_prof_id,
                                      &train_prof_rep,
                                      &nprof_train))
    return 0;
  np_categorical_profile_owner_take_index(owner, train_prof_id);
  np_categorical_profile_owner_take_index(owner, train_prof_rep);

  if(!np_build_discrete_profile_index(num_obs_eval,
                                      num_reg_unordered,
                                      num_reg_ordered,
                                      matrix_X_unordered_eval,
                                      matrix_X_ordered_eval,
                                      &eval_prof_id,
                                      &eval_prof_rep,
                                      &nprof_eval))
    goto cleanup;
  np_categorical_profile_owner_take_index(owner, eval_prof_id);
  np_categorical_profile_owner_take_index(owner, eval_prof_rep);

  if((nprof_train <= 0) || (nprof_eval <= 0) ||
     (4*nprof_train > 3*num_obs_train) ||
     (4*nprof_eval > 3*num_obs_eval))
    goto cleanup;

  if(call->residual_preparation != NULL) {
    if(!np_regression_residual_preparation_valid(call->residual_preparation))
      error("regression profile residual preparation has an invalid identity map");
    profile_anchor = R_allocLD((size_t)nprof_train);
    profile_centered_sum = R_allocLD((size_t)nprof_train);
    profile_eval_to_train = (int *)R_alloc((size_t)nprof_eval, sizeof(int));
    memset(profile_centered_sum, 0, (size_t)nprof_train*sizeof(long double));
    residual_context.output = call->residual_preparation;
    residual_context.state = np_regression_residual_states(nprof_eval);
    residual_context.num_rows = nprof_eval;
    residual_context.profile_eval_to_train = profile_eval_to_train;
    residual_context.profile_anchor = profile_anchor;
    residual_context.profile_centered_sum = profile_centered_sum;
    residual_context.profile_self_weight = R_allocLD((size_t)nprof_eval);
    memset(residual_context.profile_self_weight, 0,
           (size_t)nprof_eval*sizeof(long double));
  }

  profile_unordered_train = alloc_tmatd(nprof_train, num_reg_unordered);
  np_categorical_profile_owner_take_matrix(owner, profile_unordered_train);
  profile_ordered_train = alloc_tmatd(nprof_train, num_reg_ordered);
  np_categorical_profile_owner_take_matrix(owner, profile_ordered_train);
  profile_unordered_eval = alloc_tmatd(nprof_eval, num_reg_unordered);
  np_categorical_profile_owner_take_matrix(owner, profile_unordered_eval);
  profile_ordered_eval = alloc_tmatd(nprof_eval, num_reg_ordered);
  np_categorical_profile_owner_take_matrix(owner, profile_ordered_eval);
  counts = alloc_vecd(nprof_train);
  np_categorical_profile_owner_take_vector(owner, counts);
  sums = alloc_vecd(nprof_train);
  np_categorical_profile_owner_take_vector(owner, sums);
  if(do_merr) {
    sums2 = alloc_vecd(nprof_train);
    np_categorical_profile_owner_take_vector(owner, sums2);
  }
  weighted_sum = alloc_vecd(response_column_count*nprof_eval);
  np_categorical_profile_owner_take_vector(owner, weighted_sum);
  profile_mean = alloc_vecd(nprof_eval);
  np_categorical_profile_owner_take_vector(owner, profile_mean);
  if(do_merr) {
    profile_stderr = alloc_vecd(nprof_eval);
    np_categorical_profile_owner_take_vector(owner, profile_stderr);
  }

  if(((num_reg_unordered > 0) &&
      (profile_unordered_train == NULL || profile_unordered_eval == NULL)) ||
     ((num_reg_ordered > 0) &&
      (profile_ordered_train == NULL || profile_ordered_eval == NULL)) ||
     counts == NULL || sums == NULL ||
     (do_merr && sums2 == NULL) || weighted_sum == NULL ||
     profile_mean == NULL || (do_merr && profile_stderr == NULL))
    goto cleanup;

  for(g = 0; g < nprof_train; g++){
    const int rep = train_prof_rep[g];
    counts[g] = 0.0;
    sums[g] = 0.0;
    if(profile_anchor != NULL)
      profile_anchor[g] = (long double)vector_Y[rep];
    if(do_merr)
      sums2[g] = 0.0;
    for(j = 0; j < num_reg_unordered; j++)
      profile_unordered_train[j][g] = matrix_X_unordered_train[j][rep];
    for(j = 0; j < num_reg_ordered; j++)
      profile_ordered_train[j][g] = matrix_X_ordered_train[j][rep];
  }

  for(g = 0; g < nprof_eval; g++){
    const int rep = eval_prof_rep[g];
    if(profile_eval_to_train != NULL)
      profile_eval_to_train[g] = train_prof_id[
        call->residual_preparation->evaluation_to_donor[rep]];
    for(j = 0; j < num_reg_unordered; j++)
      profile_unordered_eval[j][g] = matrix_X_unordered_eval[j][rep];
    for(j = 0; j < num_reg_ordered; j++)
      profile_ordered_eval[j][g] = matrix_X_ordered_eval[j][rep];
  }

  for(i = 0; i < num_obs_train; i++){
    g = train_prof_id[i];
    counts[g] += 1.0;
    sums[g] += vector_Y[i];
    if(profile_centered_sum != NULL)
      profile_centered_sum[g] += (long double)vector_Y[i] - profile_anchor[g];
    if(do_merr)
      sums2[g] += vector_Y[i]*vector_Y[i];
    if(emit_fit_progress)
      np_progress_fit_loop_step(i + 1, num_obs_train);
  }

  profile_y[0] = sums;
  profile_y[1] = counts;
  profile_y[2] = sums2;

  if(call->residual_preparation != NULL) {
    residual_context.profile_count = counts;
    status = np_regression_profile_residual_kernel_sum(
      call, nprof_train, nprof_eval,
      profile_unordered_train, profile_ordered_train,
      profile_unordered_eval, profile_ordered_eval,
      profile_y, response_column_count, weighted_sum, &residual_context);
  } else {
    status = kernel_weighted_sum_np(kernel_c,
                                  kernel_u,
                                  kernel_o,
                                  BANDWIDTH_reg,
                                  nprof_train,
                                  nprof_eval,
                                  num_reg_unordered,
                                  num_reg_ordered,
                                  0,
                                  0,
                                  0,
                                  1,
                                  1,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  operator,
                                  OP_NOOP,
                                  0,
                                  0,
                                  NULL,
                                  1,
                                  response_column_count,
                                  0,
                                  NP_TREE_FALSE,
                                  0,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  profile_unordered_train,
                                  profile_ordered_train,
                                  NULL,
                                  profile_unordered_eval,
                                  profile_ordered_eval,
                                  NULL,
                                  profile_y,
                                  NULL,
                                  NULL,
                                  vector_scale_factor,
                                  1,
                                  NULL,
                                  NULL,
                                  lambda,
                                  num_categories,
                                  matrix_categorical_vals,
                                  NULL,
                                  weighted_sum,
                                  NULL,
                                  NULL,
                                  NULL);
  }

  if(status != 0)
    goto cleanup;

  for(g = 0; g < nprof_eval; g++){
    const double A = weighted_sum[response_column_count*g];
    const double B = weighted_sum[response_column_count*g + 1];

    if(!(fabs(B) > DBL_MIN) || !R_FINITE(B))
      goto cleanup;

    profile_mean[g] = A/B;
    if(do_merr) {
      const double C = weighted_sum[response_column_count*g + 2];
      const double v = C/B - profile_mean[g]*profile_mean[g];
      profile_stderr[g] = (v <= 0.0 || !R_FINITE(v)) ? 0.0 : sqrt(v/B);
    }
  }

  for(i = 0; i < num_obs_eval; i++){
    g = eval_prof_id[i];
    mean[i] = profile_mean[g];
    if(do_merr)
      mean_stderr[i] = profile_stderr[g];
    if(call->residual_preparation != NULL) {
      NPRegressionResidualPreparation * const preparation =
        call->residual_preparation;
      const int donor = preparation->evaluation_to_donor[i];
      const int self_profile = profile_eval_to_train[g];
      NPResidualOffDiagonal state = residual_context.state[g];
      long double residual = 0.0L;
      if(state.scale > 0.0L) {
        const long double normalized_total = state.sum +
          residual_context.profile_self_weight[g] / state.scale;
        state.numerator += normalized_total *
          ((long double)preparation->response[donor] -
           profile_anchor[self_profile]);
      }
      NPResidualInformation information = np_residual_offdiag_finish(
        &state, (long double)weighted_sum[response_column_count*g + 1],
        &residual);
      if(information != NP_RESIDUAL_INVALID &&
         !R_FINITE((double)residual))
        information = NP_RESIDUAL_INVALID;
      preparation->normalized_residual[donor] = (double)residual;
      preparation->information[donor] = (int)information;
    }
  }

  ok = 1;

cleanup:
  np_categorical_profile_owner_clear(owner);

  return ok;
}
