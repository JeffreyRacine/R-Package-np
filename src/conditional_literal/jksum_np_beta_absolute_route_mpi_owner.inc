static int np_beta_absolute_route_mpi_owner(
  const NPBetaAbsoluteRouteCall * const call)
{
  const int stride =
    call->num_obs_eval / iNum_Processors +
    ((call->num_obs_eval % iNum_Processors) != 0);
  const size_t evaluation_start_size = MIN(
    (size_t)my_rank * (size_t)stride, (size_t)call->num_obs_eval);
  const int evaluation_start = (int)evaluation_start_size;
  const int evaluation_count = MIN(
    stride, call->num_obs_eval - evaluation_start);
  const int original_train_is_eval =
    np_beta_absolute_route_train_is_eval(call);
  NPBetaAbsoluteRouteCall local_call = *call;
  NPNNKernelFold local_fold;
  NPContinuousKernelDerivativeDiagnostics local_diagnostics = {
    .bad_coordinate = -1,
    .bad_observation = -1,
    .undefined_count = 0,
    .beta_status = NP_BETA_OK
  };
  double **continuous_eval_view = NULL;
  double **unordered_eval_view = NULL;
  double **ordered_eval_view = NULL;
  double **bandwidth_eval_view = NULL;
  int *recvcounts = NULL;
  int *displs = NULL;
  int primary_width = 0;
  int power2_width = 0;
  int status = KWSNP_ERR_BADINVOC;
  int route_ownership_consumed = 0;
  int failing_rank;
  int failure_payload[5];
  int undefined_count;

  if(np_mpi_rank_failure_injected(
       "NP_RMPI_INJECT_BETA_NPKSUM_FAIL_RANK") ||
     call->num_obs_train <= 0 || call->num_obs_eval <= 0 ||
     call->num_reg_continuous <= 0 ||
     call->leave_one_out_offset < 0 ||
     evaluation_start > INT_MAX - call->leave_one_out_offset ||
     call->weighted_sum == NULL || stride <= 0 ||
     !np_beta_absolute_route_output_width(
       call->ncol_Y, call->ncol_W, &primary_width) ||
     (call->weighted_sum_power2 != NULL &&
      !np_beta_absolute_route_output_width(
        call->ncol_Y_power2, call->ncol_W_power2, &power2_width)) ||
     (size_t)iNum_Processors > SIZE_MAX / sizeof(*recvcounts))
    goto rendezvous;

  recvcounts = (int *)malloc(
    (size_t)iNum_Processors * sizeof(*recvcounts));
  displs = (int *)malloc((size_t)iNum_Processors * sizeof(*displs));
  if(recvcounts == NULL || displs == NULL)
    goto rendezvous;

  status = 0;
  if(evaluation_count > 0) {
    continuous_eval_view = np_beta_absolute_route_matrix_view(
      call->matrix_X_continuous_eval,
      call->num_reg_continuous, evaluation_start);
    unordered_eval_view = np_beta_absolute_route_matrix_view(
      call->matrix_X_unordered_eval,
      call->num_reg_unordered, evaluation_start);
    ordered_eval_view = np_beta_absolute_route_matrix_view(
      call->matrix_X_ordered_eval,
      call->num_reg_ordered, evaluation_start);
    if(continuous_eval_view == NULL ||
       (call->num_reg_unordered > 0 && unordered_eval_view == NULL) ||
       (call->num_reg_ordered > 0 && ordered_eval_view == NULL))
      status = KWSNP_ERR_BADINVOC;

    if(status == 0 && call->bandwidth_mode == BW_GEN_NN) {
      bandwidth_eval_view = np_beta_absolute_route_matrix_view(
        call->matrix_bandwidth_eval,
        call->num_reg_continuous, evaluation_start);
      if(bandwidth_eval_view == NULL)
        status = KWSNP_ERR_BADINVOC;
    }

    if(status == 0) {
      local_call.num_obs_eval = evaluation_count;
      if(call->fold != NULL) {
        local_fold = *call->fold;
        local_fold.num_eval = evaluation_count;
        local_fold.eval_to_train += evaluation_start;
        local_call.fold = &local_fold;
      }
      local_call.leave_one_out_train_is_eval = original_train_is_eval;
      local_call.leave_one_out_offset += evaluation_start;
      local_call.matrix_X_continuous_eval = continuous_eval_view;
      local_call.matrix_X_unordered_eval = unordered_eval_view;
      local_call.matrix_X_ordered_eval = ordered_eval_view;
      if(call->bandwidth_mode == BW_GEN_NN)
        local_call.matrix_bandwidth_eval = bandwidth_eval_view;
      local_call.weighted_sum = call->weighted_sum +
        (size_t)evaluation_start * (size_t)primary_width;
      if(call->weighted_sum_power2 != NULL)
        local_call.weighted_sum_power2 = call->weighted_sum_power2 +
          (size_t)evaluation_start * (size_t)power2_width;
      if(call->centered_m2 != NULL)
        local_call.centered_m2 = call->centered_m2 + evaluation_start;
      if(call->kw != NULL)
        local_call.kw = call->kw +
          (size_t)evaluation_start * (size_t)call->num_obs_train;
      local_call.route_diagnostics = &local_diagnostics;
      status = np_beta_absolute_route(&local_call);
      route_ownership_consumed = 1;
    }
  }

rendezvous:
  if(!route_ownership_consumed) {
    free(call->row);
    free_tmat(call->owned_bandwidth_tmatrix);
  }
  free(continuous_eval_view);
  free(unordered_eval_view);
  free(ordered_eval_view);
  free(bandwidth_eval_view);

  failing_rank = status == 0 ? iNum_Processors : my_rank;
  MPI_Allreduce(MPI_IN_PLACE, &failing_rank, 1, MPI_INT, MPI_MIN, comm[1]);
  if(failing_rank < iNum_Processors) {
    if(my_rank == failing_rank) {
      failure_payload[0] = status;
      failure_payload[1] = local_diagnostics.bad_coordinate;
      failure_payload[2] = local_diagnostics.bad_observation;
      failure_payload[3] = local_diagnostics.undefined_count;
      failure_payload[4] = (int)local_diagnostics.beta_status;
    }
    MPI_Bcast(failure_payload, 5, MPI_INT, failing_rank, comm[1]);
    if(call->route_diagnostics != NULL) {
      call->route_diagnostics->bad_coordinate = failure_payload[1];
      call->route_diagnostics->bad_observation = failure_payload[2];
      call->route_diagnostics->undefined_count = failure_payload[3];
      call->route_diagnostics->beta_status =
        (np_beta_status)failure_payload[4];
    }
    free(recvcounts);
    free(displs);
    return failure_payload[0];
  }

  if(!np_beta_absolute_route_gather_rows(
       call->num_obs_eval, stride, primary_width,
       call->weighted_sum, recvcounts, displs,
       "beta kernel sum MPI_Allgatherv") ||
     (call->weighted_sum_power2 != NULL &&
      !np_beta_absolute_route_gather_rows(
        call->num_obs_eval, stride, power2_width,
        call->weighted_sum_power2, recvcounts, displs,
        "beta dual-power kernel sum MPI_Allgatherv")) ||
     (call->centered_m2 != NULL &&
      !np_beta_absolute_route_gather_rows(
        call->num_obs_eval, stride, 1, call->centered_m2,
        recvcounts, displs,
        "beta centered moment MPI_Allgatherv")) ||
     (call->kw != NULL &&
      !np_beta_absolute_route_gather_rows(
        call->num_obs_eval, stride, call->num_obs_train, call->kw,
        recvcounts, displs,
        "beta kernel weights MPI_Allgatherv"))) {
    free(recvcounts);
    free(displs);
    return KWSNP_ERR_BADINVOC;
  }

  undefined_count = local_diagnostics.undefined_count;
  MPI_Allreduce(MPI_IN_PLACE, &undefined_count, 1, MPI_INT, MPI_SUM, comm[1]);
  if(call->route_diagnostics != NULL) {
    call->route_diagnostics->bad_coordinate = -1;
    call->route_diagnostics->bad_observation = -1;
    call->route_diagnostics->undefined_count = undefined_count;
    call->route_diagnostics->beta_status = NP_BETA_OK;
  }
  free(recvcounts);
  free(displs);
  return 0;
}
