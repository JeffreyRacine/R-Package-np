static NP_NOINLINE int np_beta_scalar_regression_fit_canonical(
  const int BANDWIDTH_reg,
  const int KERNEL_unordered_reg,
  const int KERNEL_ordered_reg,
  const int num_obs_train,
  const int num_obs_eval,
  const int num_reg_unordered,
  const int num_reg_ordered,
  const int num_reg_continuous,
  double **matrix_X_unordered_train,
  double **matrix_X_ordered_train,
  double **matrix_X_continuous_train,
  double **matrix_X_unordered_eval,
  double **matrix_X_ordered_eval,
  double **matrix_X_continuous_eval,
  const double *vector_Y,
  double *vector_scale_factor,
  int *num_categories,
  double **matrix_categorical_vals,
  double *mean,
  double **gradient,
  double *mean_stderr,
  double **gradient_stderr,
  const NPContinuousKernelRoute *kernel_route,
  NPContinuousKernelDerivativeDiagnostics *kernel_route_diagnostics,
  const int categorical_compress,
  const NPRegressionStandardErrorMode standard_error_mode,
  const int ordinary_response,
  const NPContinuousPreparedBandwidthView *prepared_bandwidth,
  const NPNNGeometryContext *nn_geometry_context,
  const NPRegressionHC0Context *hc0_context,
  NPRegressionFailure *failure,
  NPRegressionLPEmptyRows *empty_rows)
{
  const int original_train_is_eval =
    matrix_X_continuous_train == matrix_X_continuous_eval;
  np_beta_bandwidth_prepare_status bandwidth_status =
    NP_BETA_BANDWIDTH_PREPARE_OK;
  NPContinuousKernelRowStatus row_status = NP_CONTINUOUS_ROW_OK;
  int status;

#ifdef MPI2
  if(iNum_Processors > 1 && !np_mpi_local_regression_active()) {
    const int num_predictors =
      num_reg_continuous + num_reg_unordered + num_reg_ordered;
    const int stride =
      num_obs_eval/iNum_Processors +
      ((num_obs_eval%iNum_Processors) != 0);
    size_t evaluation_start_size;
    int evaluation_start;
    int evaluation_count;
    int *recvcounts;
    int *displs;
    int equal_counts = 0;
    double **unordered_eval_view = NULL;
    double **ordered_eval_view = NULL;
    double **continuous_eval_view = NULL;
    double **gradient_view = NULL;
    double **gradient_stderr_view = NULL;
    double **owned_prepared_matrix = NULL;
    NPContinuousPreparedBandwidthView owned_prepared_view;
    NPContinuousPreparedBandwidthView prepared_slice;
    const NPContinuousPreparedBandwidthView *prepared_source =
      prepared_bandwidth;
    const NPContinuousPreparedBandwidthView *prepared_view = NULL;
    NPRegressionHC0Context hc0_slice;
    NPRegressionResidualPreparation preparation_slice;
    const NPRegressionHC0Context *hc0_view = hc0_context;
    NPRegressionResidualPreparation * const preparation =
      hc0_context != NULL &&
      hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_PREPARING ?
        hc0_context->preparation : NULL;
    double *preparation_values = NULL;
    int *preparation_information = NULL;
    NPContinuousKernelDerivativeDiagnostics local_diagnostics = {
      .bad_coordinate = -1,
      .bad_observation = -1,
      .undefined_count = 0,
      .beta_status = NP_BETA_OK
    };
    int failing_rank;
    int payload[7];
    int coordinate;

    if(num_obs_eval <= 0 || stride <= 0 || num_predictors < 0) {
      np_beta_scalar_regression_fit_error(
        NP_BETA_SCALAR_REGRESSION_FIT_ERR_MPI_LAYOUT,
        bandwidth_status, row_status, kernel_route_diagnostics, NULL);
    }
    evaluation_start_size = MIN(
      (size_t)my_rank*(size_t)stride, (size_t)num_obs_eval);
    evaluation_start = (int)evaluation_start_size;
    evaluation_count = MIN(stride, num_obs_eval - evaluation_start);
    recvcounts = (int *)R_alloc((size_t)iNum_Processors,
                                (int)sizeof(*recvcounts));
    displs = (int *)R_alloc((size_t)iNum_Processors,
                            (int)sizeof(*displs));
    if(!np_jksum_mpi_contiguous_row_layout(
         num_obs_eval, stride, 1, recvcounts, displs, &equal_counts)) {
      np_beta_scalar_regression_fit_error(
        NP_BETA_SCALAR_REGRESSION_FIT_ERR_MPI_LAYOUT,
        bandwidth_status, row_status, kernel_route_diagnostics, NULL);
    }
    (void)equal_counts;

    if(preparation != NULL) {
      /* The source is a validated whole-training map.  Rank-local beta rows
       * consume a slice of that map, but write by global donor identity. */
      preparation_slice = *preparation;
      preparation_slice.evaluation_to_donor += evaluation_start;
      preparation_slice.num_obs_eval = evaluation_count;
      preparation_values = (double *)R_alloc(
        (size_t)num_obs_eval, (int)sizeof(*preparation_values));
      preparation_information = (int *)R_alloc(
        (size_t)num_obs_eval, (int)sizeof(*preparation_information));
    }
    if(hc0_context != NULL) {
      hc0_slice = *hc0_context;
      if(preparation != NULL)
        hc0_slice.preparation = &preparation_slice;
      if(hc0_context->mean_unavailable != NULL)
        hc0_slice.mean_unavailable += evaluation_start;
      hc0_slice.gradient_unavailable = np_beta_scalar_regression_mask_view(
        hc0_context->gradient_unavailable, num_predictors, evaluation_start);
      hc0_slice.gradient_structural_zero = np_beta_scalar_regression_mask_view(
        hc0_context->gradient_structural_zero, num_predictors, evaluation_start);
      hc0_view = &hc0_slice;
    }

    /*
     * Nearest-neighbour realization contains legacy rank-symmetric MPI
     * collectives.  Complete it once on the unsliced invocation, then give
     * every rank an immutable slice.  Entering it from unequal rank-local
     * blocks would mismatch those collectives and would also repeat setup.
     */
    status = NP_BETA_SCALAR_REGRESSION_FIT_OK;
    if(BANDWIDTH_reg != BW_FIXED && prepared_source == NULL) {
      const int bandwidth_rows =
        BANDWIDTH_reg == BW_GEN_NN ? num_obs_eval : num_obs_train;

      owned_prepared_matrix = alloc_tmatd(
        bandwidth_rows, num_reg_continuous);
      bandwidth_status = np_beta_bandwidth_prepare_matrix_ctx(
        BANDWIDTH_reg == BW_GEN_NN ?
          NP_BETA_BANDWIDTH_GENERALIZED_NN :
          NP_BETA_BANDWIDTH_ADAPTIVE_NN,
        matrix_X_continuous_train, matrix_X_continuous_eval,
        vector_scale_factor,
        num_obs_train, num_obs_eval, num_reg_continuous,
        original_train_is_eval,
        BANDWIDTH_reg == BW_GEN_NN, BANDWIDTH_reg == BW_ADAP_NN,
        0, owned_prepared_matrix, owned_prepared_matrix, nn_geometry_context);
      if(bandwidth_status != NP_BETA_BANDWIDTH_PREPARE_OK) {
        status = NP_BETA_SCALAR_REGRESSION_FIT_ERR_BANDWIDTH;
      } else {
        owned_prepared_view.bandwidth_mode = BANDWIDTH_reg;
        owned_prepared_view.num_train = num_obs_train;
        owned_prepared_view.num_eval_total = num_obs_eval;
        owned_prepared_view.num_continuous = num_reg_continuous;
        owned_prepared_view.evaluation_offset = 0;
        owned_prepared_view.evaluation_count = num_obs_eval;
        owned_prepared_view.bandwidth_eval =
          BANDWIDTH_reg == BW_GEN_NN ? owned_prepared_matrix : NULL;
        owned_prepared_view.bandwidth_train =
          BANDWIDTH_reg == BW_ADAP_NN ? owned_prepared_matrix : NULL;
        prepared_source = &owned_prepared_view;
      }
    }

    if(status == NP_BETA_SCALAR_REGRESSION_FIT_OK && evaluation_count > 0) {
      unordered_eval_view = np_beta_scalar_regression_fit_matrix_view(
        matrix_X_unordered_eval, num_reg_unordered, evaluation_start);
      ordered_eval_view = np_beta_scalar_regression_fit_matrix_view(
        matrix_X_ordered_eval, num_reg_ordered, evaluation_start);
      continuous_eval_view = np_beta_scalar_regression_fit_matrix_view(
        matrix_X_continuous_eval, num_reg_continuous, evaluation_start);
      if((num_reg_unordered > 0 && unordered_eval_view == NULL) ||
         (num_reg_ordered > 0 && ordered_eval_view == NULL) ||
         (num_reg_continuous > 0 && continuous_eval_view == NULL))
        status = NP_BETA_SCALAR_REGRESSION_FIT_ERR_MPI_LAYOUT;

      if(status == NP_BETA_SCALAR_REGRESSION_FIT_OK && gradient != NULL) {
        gradient_view = np_beta_scalar_regression_fit_matrix_view(
          gradient, num_predictors, evaluation_start);
        if(gradient_view == NULL)
          status = NP_BETA_SCALAR_REGRESSION_FIT_ERR_MPI_LAYOUT;
      }
      if(status == NP_BETA_SCALAR_REGRESSION_FIT_OK &&
         gradient_stderr != NULL) {
        gradient_stderr_view = np_beta_scalar_regression_fit_matrix_view(
          gradient_stderr, num_predictors, evaluation_start);
        if(gradient_stderr_view == NULL)
          status = NP_BETA_SCALAR_REGRESSION_FIT_ERR_MPI_LAYOUT;
      }
      if(status == NP_BETA_SCALAR_REGRESSION_FIT_OK &&
         prepared_source != NULL) {
        const int prepared_status =
          np_beta_scalar_regression_prepared_view_slice(
            prepared_source, BANDWIDTH_reg,
            num_obs_train, num_obs_eval, num_reg_continuous,
            evaluation_start, evaluation_count, &prepared_slice);
        if(prepared_status < 0)
          status = NP_BETA_SCALAR_REGRESSION_FIT_ERR_PREPARED_BANDWIDTH;
        else if(prepared_status > 0)
          prepared_view = &prepared_slice;
      }
      if(status == NP_BETA_SCALAR_REGRESSION_FIT_OK)
        status = np_beta_scalar_regression_fit_block_canonical(
          BANDWIDTH_reg, KERNEL_unordered_reg, KERNEL_ordered_reg,
          num_obs_train, evaluation_count,
          num_reg_unordered, num_reg_ordered, num_reg_continuous,
          matrix_X_unordered_train, matrix_X_ordered_train,
          matrix_X_continuous_train,
          unordered_eval_view, ordered_eval_view, continuous_eval_view,
          vector_Y, vector_scale_factor,
          num_categories, matrix_categorical_vals,
          mean + evaluation_start, gradient_view,
          mean_stderr == NULL ? NULL : mean_stderr + evaluation_start,
          gradient_stderr_view,
          kernel_route, &local_diagnostics, categorical_compress,
          standard_error_mode, ordinary_response, prepared_view, original_train_is_eval,
          NULL, hc0_view,
          &bandwidth_status, &row_status, NULL);
    }
    if(owned_prepared_matrix != NULL)
      free_tmat(owned_prepared_matrix);

    failing_rank =
      status == NP_BETA_SCALAR_REGRESSION_FIT_OK ?
      iNum_Processors : my_rank;
    MPI_Allreduce(MPI_IN_PLACE, &failing_rank, 1, MPI_INT, MPI_MIN, comm[1]);
    if(failing_rank < iNum_Processors) {
      if(my_rank == failing_rank) {
        payload[0] = status;
        payload[1] = (int)bandwidth_status;
        payload[2] = (int)row_status;
        payload[3] = local_diagnostics.bad_coordinate;
        payload[4] = local_diagnostics.bad_observation;
        payload[5] = local_diagnostics.undefined_count;
        payload[6] = (int)local_diagnostics.beta_status;
      }
      MPI_Bcast(payload, 7, MPI_INT, failing_rank, comm[1]);
      if(kernel_route_diagnostics != NULL) {
        kernel_route_diagnostics->bad_coordinate = payload[3];
        kernel_route_diagnostics->bad_observation = payload[4];
        kernel_route_diagnostics->undefined_count = payload[5];
        kernel_route_diagnostics->beta_status = (np_beta_status)payload[6];
      }
      if(payload[0] == NP_BETA_SCALAR_REGRESSION_FIT_ERR_BANDWIDTH &&
         payload[1] == NP_BETA_BANDWIDTH_PREPARE_ERR_ZERO_RADIUS) {
        const NPNNZeroRadiusInfo info = np_nn_zero_radius_info(
          BANDWIDTH_reg, num_obs_train, num_obs_eval, num_reg_continuous,
          matrix_X_continuous_train, matrix_X_continuous_eval,
          vector_scale_factor, nn_geometry_context, 0);
        np_nn_zero_radius_error(&info);
      }
      np_beta_scalar_regression_fit_error(
        payload[0], (np_beta_bandwidth_prepare_status)payload[1],
        (NPContinuousKernelRowStatus)payload[2], kernel_route_diagnostics, NULL);
    }

    if(preparation != NULL) {
      /* Gather in evaluation order, then restore the explicit donor order.
       * This remains valid for non-identity maps and empty rank chunks. */
      int evaluation;
      for(evaluation = evaluation_start;
          evaluation < evaluation_start + evaluation_count; ++evaluation) {
        const int donor = preparation->evaluation_to_donor[evaluation];
        preparation_values[evaluation] = preparation->normalized_residual[donor];
        preparation_information[evaluation] = preparation->information[donor];
      }
      np_mpi_allgatherv_in_place_double(
        evaluation_count, preparation_values, recvcounts, displs,
        "beta scalar regression residual preparation MPI_Allgatherv");
      MPI_Allgatherv(MPI_IN_PLACE, evaluation_count, MPI_INT,
        preparation_information, recvcounts, displs, MPI_INT, comm[1]);
      for(evaluation = 0; evaluation < num_obs_eval; ++evaluation) {
        const int donor = preparation->evaluation_to_donor[evaluation];
        preparation->normalized_residual[donor] = preparation_values[evaluation];
        preparation->information[donor] = preparation_information[evaluation];
      }
    }
    np_mpi_allgatherv_in_place_double(
      evaluation_count, mean, recvcounts, displs,
      "beta scalar regression mean MPI_Allgatherv");
    if(mean_stderr != NULL)
      np_mpi_allgatherv_in_place_double(
        evaluation_count, mean_stderr, recvcounts, displs,
        "beta scalar regression standard error MPI_Allgatherv");
    if(hc0_context != NULL && hc0_context->mean_unavailable != NULL)
      MPI_Allgatherv(MPI_IN_PLACE, evaluation_count, MPI_INT,
        hc0_context->mean_unavailable, recvcounts, displs, MPI_INT, comm[1]);
    for(coordinate = 0; coordinate < num_predictors; ++coordinate) {
      if(gradient != NULL)
        np_mpi_allgatherv_in_place_double(
          evaluation_count, gradient[coordinate], recvcounts, displs,
          "beta scalar regression gradient MPI_Allgatherv");
      if(gradient_stderr != NULL)
        np_mpi_allgatherv_in_place_double(
          evaluation_count, gradient_stderr[coordinate], recvcounts, displs,
          "beta scalar regression gradient standard error MPI_Allgatherv");
      if(hc0_context != NULL &&
         hc0_context->gradient_unavailable != NULL &&
         hc0_context->gradient_unavailable[coordinate] != NULL)
        MPI_Allgatherv(MPI_IN_PLACE, evaluation_count, MPI_INT,
          hc0_context->gradient_unavailable[coordinate], recvcounts, displs,
          MPI_INT, comm[1]);
      if(hc0_context != NULL &&
         hc0_context->gradient_structural_zero != NULL &&
         hc0_context->gradient_structural_zero[coordinate] != NULL)
        MPI_Allgatherv(MPI_IN_PLACE, evaluation_count, MPI_INT,
          hc0_context->gradient_structural_zero[coordinate], recvcounts, displs,
          MPI_INT, comm[1]);
    }

    {
      const int diagnostics_rank =
        MIN(iNum_Processors - 1, (num_obs_eval - 1)/stride);
      if(my_rank == diagnostics_rank) {
        payload[0] = local_diagnostics.bad_coordinate;
        payload[1] = local_diagnostics.bad_observation;
        payload[2] = local_diagnostics.undefined_count;
        payload[3] = (int)local_diagnostics.beta_status;
      }
      MPI_Bcast(payload, 4, MPI_INT, diagnostics_rank, comm[1]);
      if(kernel_route_diagnostics != NULL) {
        kernel_route_diagnostics->bad_coordinate = payload[0];
        kernel_route_diagnostics->bad_observation = payload[1];
        kernel_route_diagnostics->undefined_count = payload[2];
        kernel_route_diagnostics->beta_status = (np_beta_status)payload[3];
      }
    }
    np_progress_fit_loop_step(num_obs_eval, num_obs_eval);
    return 0;
  }
#endif

  status = np_beta_scalar_regression_fit_block_canonical(
    BANDWIDTH_reg, KERNEL_unordered_reg, KERNEL_ordered_reg,
    num_obs_train, num_obs_eval,
    num_reg_unordered, num_reg_ordered, num_reg_continuous,
    matrix_X_unordered_train, matrix_X_ordered_train,
    matrix_X_continuous_train,
    matrix_X_unordered_eval, matrix_X_ordered_eval,
    matrix_X_continuous_eval,
    vector_Y, vector_scale_factor,
    num_categories, matrix_categorical_vals,
    mean, gradient, mean_stderr, gradient_stderr,
    kernel_route, kernel_route_diagnostics, categorical_compress,
    standard_error_mode, ordinary_response, prepared_bandwidth, original_train_is_eval,
    nn_geometry_context, hc0_context,
    &bandwidth_status, &row_status, empty_rows);
  if(status == NP_BETA_SCALAR_REGRESSION_FIT_ERR_BANDWIDTH &&
     bandwidth_status == NP_BETA_BANDWIDTH_PREPARE_ERR_ZERO_RADIUS) {
    const NPNNZeroRadiusInfo info = np_nn_zero_radius_info(
      BANDWIDTH_reg, num_obs_train, num_obs_eval, num_reg_continuous,
      matrix_X_continuous_train, matrix_X_continuous_eval,
      vector_scale_factor, nn_geometry_context, 0);
    if(failure != NULL) {
      np_regression_failure_zero_radius(failure, &info);
      return NP_REGRESSION_FIT_ERR_ZERO_NN_RADIUS;
    }
    np_nn_zero_radius_error(&info);
  }
  if(status != NP_BETA_SCALAR_REGRESSION_FIT_OK)
    return np_beta_scalar_regression_fit_error(
      status, bandwidth_status, row_status, kernel_route_diagnostics, failure);
  return 0;
}
