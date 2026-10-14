static SEXP np_regression_scalar_fit_execute(void *data)
{
  NPRegressionScalarFitExecution * const execution =
    (NPRegressionScalarFitExecution *)data;
  const NPRegressionScalarFitCall * const call = execution->call;
  NPRegressionScalarFitOwner * const owner = &execution->owner;
  const int conditional_influence =
    call->standard_error_mode == NP_REGRESSION_STDERR_CONDITIONAL_INFLUENCE;
  const int conditional_se = conditional_influence && call->do_merr;
  const int hc0_residual_preparing = call->hc0_context != NULL &&
    call->hc0_context->status == NP_REGRESSION_HC0_RESIDUAL_PREPARING;
  const int ordinary_hc0 = call->hc0_context != NULL &&
    !hc0_residual_preparing;
  const int point_already_computed = ordinary_hc0 &&
    call->hc0_context->point_already_computed;
  const int defer_hc0_mean_stderr = ordinary_hc0 && call->do_gerr &&
    call->bandwidth_mode == BW_ADAP_NN;
  const int response_column_count = ordinary_hc0 ?
    (point_already_computed ? 1 : 2) : (call->do_merr ? 3 : 2);
  const int p_nvar = call->do_grad ?
    (call->num_reg_continuous + call->num_reg_unordered +
     call->num_reg_ordered) : 0;
  const NPRegressionGradientRequest *gradient_request =
    call->hc0_context != NULL ? call->hc0_context->gradient_request : NULL;
  const NPRegressionGradientRange continuous_range =
    np_regression_gradient_range(gradient_request, 0, call->num_reg_continuous);
  const NPRegressionGradientRange unordered_range = np_regression_gradient_range(
    gradient_request, call->num_reg_continuous,
    call->num_reg_continuous + call->num_reg_unordered);
  const NPRegressionGradientRange ordered_range = np_regression_gradient_range(
    gradient_request, call->num_reg_continuous + call->num_reg_unordered, p_nvar);
  const NPRegressionGradientRange categorical_range = np_regression_gradient_range(
    gradient_request,call->num_reg_continuous,p_nvar);
  const int stable_scalar_contrast = call->ordinary_response && call->do_grad &&
    categorical_range.begin < categorical_range.end;
  double *response_columns[NP_REGRESSION_SCALAR_RESPONSE_COLUMNS_MAX] = {
    NULL, NULL, NULL
  };
  double *hc0_residual_columns[1] = {NULL};
  NP_DualPowerCtx hc0_dual_power_ctx = {0};
  NP_DualPowerCtx residual_dual_power_ctx = {0};
  NP_DualPowerCtx contrast_dual_power_ctx = {0};
  NPRegressionScalarContrastContext scalar_contrast_context = {0};
  NPRegressionScalarANNContrast scalar_ann_contrast = {0};
  NPRegressionResidualKernelContext residual_kernel_context = {0};
  NPRegressionHC0DerivativeMomentCtx hc0_derivative_ctx = {0};
  int permutation_operator = OP_NOOP;
  size_t allocation_count;
  size_t allocation_bytes;
  int i;
  int predictor;
  const NPEmptyXRowContext empty_context = {
    .n = call->num_obs_train, .bw = call->bandwidth_mode,
    .nu = call->num_reg_unordered, .no = call->num_reg_ordered, .nc = call->num_reg_continuous,
    .kc = call->kernel_c, .ku = call->kernel_u, .ko = call->kernel_o, .op = call->operator,
    .train = {call->matrix_X_unordered_train, call->matrix_X_ordered_train, call->matrix_X_continuous_train},
    .eval = {call->matrix_X_unordered_eval, call->matrix_X_ordered_eval, call->matrix_X_continuous_eval},
    .bandwidth = call->matrix_bandwidth, .lambda = call->lambda, .vsf = call->vector_scale_factor,
    .categories = call->num_categories, .catvals = call->matrix_categorical_vals
  };

  allocation_count =
    (size_t)response_column_count *
    (size_t)call->num_obs_eval_alloc;
  if(!np_size_array_bytes_checked(
       allocation_count, sizeof(double), &allocation_bytes)) {
    execution->status = NP_REGRESSION_SCALAR_FIT_ERR_WORKSPACE_DIMENSION;
    return R_NilValue;
  }
  owner->mean_columns = (double *)malloc(allocation_bytes);

  if(conditional_se) {
    owner->conditional_weights = (double *)malloc(
      (size_t)call->num_obs_train * sizeof(double));
    if(p_nvar > 0 &&
       (!np_size_mul_checked((size_t)call->num_obs_train,
                             (size_t)p_nvar, &allocation_count) ||
        !np_size_array_bytes_checked(
          allocation_count, sizeof(double), &allocation_bytes))) {
      execution->status = NP_REGRESSION_SCALAR_FIT_ERR_DIMENSION;
      return R_NilValue;
    }
    if(p_nvar > 0)
      owner->conditional_permutation_weights =
        (double *)malloc(allocation_bytes);
    if(owner->conditional_weights == NULL ||
       (p_nvar > 0 && owner->conditional_permutation_weights == NULL)) {
      execution->status = NP_REGRESSION_SCALAR_FIT_ERR_CONDITIONAL_ALLOC;
      return R_NilValue;
    }
  }

  if(call->do_grad) {
    if(!np_size_mul_checked(
         (size_t)response_column_count,
         (size_t)call->num_obs_eval_alloc, &allocation_count) ||
       !np_size_mul_checked(
         allocation_count, (size_t)p_nvar, &allocation_count) ||
       !np_size_array_bytes_checked(
         allocation_count, sizeof(double), &allocation_bytes)) {
      execution->status = NP_REGRESSION_SCALAR_FIT_ERR_WORKSPACE_DIMENSION;
      return R_NilValue;
    }
    owner->permutation_columns = (double *)malloc(allocation_bytes);
    if(owner->permutation_columns == NULL) {
      execution->status = NP_REGRESSION_SCALAR_FIT_ERR_ALLOC;
      return R_NilValue;
    }
    permutation_operator = OP_DERIVATIVE;
  }

  if(owner->mean_columns == NULL) {
    execution->status = NP_REGRESSION_SCALAR_FIT_ERR_ALLOC;
    return R_NilValue;
  }

  owner->unit_response = (double *)malloc(
    (size_t)call->num_obs_train * sizeof(double));
  if(call->do_merr && !ordinary_hc0)
    owner->squared_response = (double *)malloc(
      (size_t)call->num_obs_train * sizeof(double));
  if(owner->unit_response == NULL ||
     (call->do_merr && !ordinary_hc0 && owner->squared_response == NULL)) {
    execution->status = NP_REGRESSION_SCALAR_FIT_ERR_ALLOC;
    return R_NilValue;
  }

  if(point_already_computed) {
    response_columns[0] = owner->unit_response;
  } else {
    response_columns[0] = call->vector_Y;
    response_columns[1] = owner->unit_response;
    if(call->do_merr && !ordinary_hc0)
      response_columns[2] = owner->squared_response;
  }

  for(i = 0; i < call->num_obs_train; i++) {
    owner->unit_response[i] = 1.0;
    if(call->do_merr && !ordinary_hc0)
      owner->squared_response[i] = call->vector_Y[i] * call->vector_Y[i];
  }

  if(ordinary_hc0) {
    hc0_residual_columns[0] =
      (double *)call->hc0_context->scaled_residual;
    hc0_dual_power_ctx.weighted_sum = call->mean_stderr;
    hc0_dual_power_ctx.kernel_pow = 2;
    hc0_dual_power_ctx.matrix_Y = hc0_residual_columns;
    hc0_dual_power_ctx.matrix_W = hc0_residual_columns;
    hc0_dual_power_ctx.ncol_Y = 1;
    hc0_dual_power_ctx.ncol_W = 1;
    hc0_dual_power_ctx.progress = NULL;
    hc0_dual_power_ctx.retain_common_scale = 0;
    hc0_dual_power_ctx.uncertainty_context = call->hc0_context;
    if(call->do_gerr && p_nvar > 0) {
      hc0_derivative_ctx.scaled_residual =
        call->hc0_context->scaled_residual;
      hc0_derivative_ctx.gradient = call->gradient;
      hc0_derivative_ctx.gradient_stderr = call->gradient_stderr;
      hc0_derivative_ctx.num_continuous = call->num_reg_continuous;
      hc0_derivative_ctx.num_predictors = p_nvar;
      hc0_derivative_ctx.num_eval = call->num_obs_eval;
      hc0_derivative_ctx.num_eval_alloc = call->num_obs_eval_alloc;
      hc0_derivative_ctx.response_column_count = response_column_count;
      hc0_derivative_ctx.denominator_column =
        point_already_computed ? 0 : 1;
      hc0_derivative_ctx.continuous_range = continuous_range;
      hc0_derivative_ctx.categorical_range = np_regression_gradient_range(
        gradient_request, call->num_reg_continuous, p_nvar);
      if(stable_scalar_contrast)
        hc0_derivative_ctx.categorical_range.end = hc0_derivative_ctx.categorical_range.begin;
      hc0_derivative_ctx.uncertainty_context = call->hc0_context;
      hc0_dual_power_ctx.regression_derivative = &hc0_derivative_ctx;
    }
  }

  if((!ordinary_hc0) && (!conditional_influence) && (!call->do_grad) &&
     np_regression_categorical_profile_fit(
       call->kernel_c, call->kernel_u, call->kernel_o,
       call->lp_engine, call->bandwidth_mode,
       call->num_obs_train, call->num_obs_eval,
       call->num_reg_unordered, call->num_reg_ordered,
       call->num_reg_continuous,
       call->matrix_X_unordered_train, call->matrix_X_ordered_train,
       call->matrix_X_unordered_eval, call->matrix_X_ordered_eval,
       call->vector_Y, call->vector_scale_factor, call->lambda,
       call->num_categories, call->operator,
       call->matrix_categorical_vals, call->mean, call->mean_stderr,
       hc0_residual_preparing,
       hc0_residual_preparing ? call->hc0_context->preparation : NULL)) {
    execution->status = NP_REGRESSION_SCALAR_FIT_PROFILE;
    return R_NilValue;
  }

  if(hc0_residual_preparing) {
    NPRegressionResidualPreparation * const preparation =
      call->hc0_context->preparation;
    if(!np_regression_residual_preparation_valid(preparation) ||
       preparation->num_obs_train != call->num_obs_train ||
       preparation->num_obs_eval != call->num_obs_eval ||
       call->do_merr || conditional_influence) {
      execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
      return R_NilValue;
    }
    residual_kernel_context.output = preparation;
    residual_kernel_context.state =
      np_regression_residual_states(call->num_obs_eval);
    residual_kernel_context.num_rows = call->num_obs_eval;
    residual_dual_power_ctx.residual_preparation = &residual_kernel_context;
  }

  if(stable_scalar_contrast) {
    scalar_contrast_context.response = call->vector_Y;
    scalar_contrast_context.gradient = call->gradient;
    scalar_contrast_context.gradient_stderr = call->gradient_stderr;
    scalar_contrast_context.uncertainty = ordinary_hc0 ? call->hc0_context : NULL;
    scalar_contrast_context.categorical_range = np_regression_gradient_range(
      gradient_request,call->num_reg_continuous,p_nvar);
    scalar_contrast_context.num_eval = call->num_obs_eval;
    scalar_contrast_context.num_eval_alloc = call->num_obs_eval_alloc;
    scalar_contrast_context.want_variance = ordinary_hc0 && call->do_gerr;
    if(call->bandwidth_mode==BW_ADAP_NN &&
       !np_regression_scalar_ann_contrast_allocate(&scalar_contrast_context,&scalar_ann_contrast)) {
      execution->status=NP_REGRESSION_SCALAR_FIT_ERR_WORKSPACE_DIMENSION;
      return R_NilValue;
    }
    hc0_dual_power_ctx.scalar_contrast = &scalar_contrast_context;
    residual_dual_power_ctx.scalar_contrast = &scalar_contrast_context;
    contrast_dual_power_ctx.scalar_contrast = &scalar_contrast_context;
  }

  {
    NPPermutationWeightOutput conditional_pkw_output = np_pkw_output_make(
      owner->conditional_permutation_weights, p_nvar);
    const int replay_information = ordinary_hc0 && call->do_gerr && p_nvar > 0 &&
      call->bandwidth_mode == BW_ADAP_NN && call->hc0_context->unknown_count > 0;
    const int replay_contrast=stable_scalar_contrast && call->bandwidth_mode==BW_ADAP_NN;
    for(int information_pass = 0;
        information_pass <= (replay_information || replay_contrast); ++information_pass) {
    NPRegressionScalarANNInformation information = {0};
    NP_DualPowerCtx information_dual = {0};
    if(replay_contrast) {
      scalar_ann_contrast.stage=information_pass;
      information_dual.scalar_contrast=&scalar_contrast_context;
    }
    if(information_pass && replay_information) {
      size_t cells;
      if(!np_size_mul_checked((size_t)p_nvar,(size_t)call->num_obs_eval,&cells) ||
         cells > (size_t)INT_MAX/5U) {
        execution->status = NP_REGRESSION_SCALAR_FIT_ERR_WORKSPACE_DIMENSION;
        return R_NilValue;
      }
      double * const denominators = (double *)R_alloc(
        (size_t)call->num_obs_eval,sizeof(double));
      double * const alternate_denominators = (double *)R_alloc(cells,sizeof(double));
      information.derivative = &hc0_derivative_ctx;
      information.level_denominator = denominators;
      information.alternate_denominator = alternate_denominators;
      information.certificate = R_allocLD(5U*cells);
      memset(information.certificate,0,5U*cells*sizeof(long double));
      const int denominator_column = point_already_computed ? 0 : 1;
      for(i = 0; i < call->num_obs_eval; ++i) {
        denominators[i] = owner->mean_columns[response_column_count*i+denominator_column];
        for(predictor = 0; predictor < p_nvar; ++predictor) {
          const size_t cell = (size_t)predictor*(size_t)call->num_obs_eval+(size_t)i;
          alternate_denominators[cell] =
            owner->permutation_columns[(size_t)response_column_count*cell+
                                       (size_t)denominator_column];
        }
      }
      information_dual.ann_information = &information;
    }
    const int weighted_sum_status = kernel_weighted_sum_np_ctx_ex(
      call->kernel_c, call->kernel_u, call->kernel_o,
      call->bandwidth_mode,
      call->num_obs_train, call->num_obs_eval,
      call->num_reg_unordered, call->num_reg_ordered,
      call->num_reg_continuous,
      0, 0, 1, 1, conditional_se, 0, 0, 0, 0,
      call->operator, permutation_operator,
      0, call->do_grad, NULL, 0,
      response_column_count, 0,
      call->tree_enabled, 0, kdt_extern_X,
      NULL, NULL, NULL,
      call->matrix_X_unordered_train,
      call->matrix_X_ordered_train,
      call->matrix_X_continuous_train,
      call->matrix_X_unordered_eval,
      call->matrix_X_ordered_eval,
      call->matrix_X_continuous_eval,
      response_columns, NULL, NULL,
      call->vector_scale_factor, 1,
      call->matrix_bandwidth, call->matrix_bandwidth,
      call->lambda, call->num_categories,
      call->matrix_categorical_vals, call->matrix_ordered_indices,
      owner->mean_columns, owner->permutation_columns,
      owner->conditional_weights, call->gate_context,
      information_pass ? &information_dual : ordinary_hc0 ? &hc0_dual_power_ctx :
        (hc0_residual_preparing ? &residual_dual_power_ctx :
         stable_scalar_contrast ? &contrast_dual_power_ctx : NULL),
      NULL, NULL, NULL,
      owner->conditional_permutation_weights != NULL ?
        &conditional_pkw_output : NULL,
        NULL);

    if(weighted_sum_status != 0) {
      execution->status = NP_REGRESSION_SCALAR_FIT_ERR_TRAVERSAL;
      return R_NilValue;
    }
    }
  }

  if(hc0_residual_preparing) {
    NPRegressionResidualPreparation * const preparation =
      call->hc0_context->preparation;
    for(i = 0; i < call->num_obs_eval; ++i) {
      const int donor = preparation->evaluation_to_donor[i];
      long double residual = 0.0L;
      NPResidualInformation information = np_residual_offdiag_finish(
        residual_kernel_context.state + i,
        (long double)owner->mean_columns[response_column_count*i + 1],
        &residual);
      if(information != NP_RESIDUAL_INVALID && !R_FINITE((double)residual))
        information = NP_RESIDUAL_INVALID;
      preparation->normalized_residual[donor] = (double)residual;
      preparation->information[donor] = (int)information;
    }
  }

  for(i = 0; i < call->num_obs_eval; i++) {
    const int response_offset = response_column_count * i;
    const int denominator_column = point_already_computed ? 0 : 1;
    const double denominator =
      owner->mean_columns[response_offset + denominator_column];

    if(!R_FINITE(denominator) || denominator == 0.0) {
      if(call->empty_rows != NULL && denominator == 0.0)
        np_empty_x_row_mark(&empty_context, i, -1, call->empty_rows);
      call->mean[i] = NA_REAL;
      if(call->do_merr)
        call->mean_stderr[i] = NA_REAL;
      continue;
    }
    if(!point_already_computed) {
      call->mean[i] = owner->mean_columns[response_offset] / denominator;
      if(!R_FINITE(call->mean[i])) {
        call->mean[i] = NA_REAL;
        if(call->do_merr)
          call->mean_stderr[i] = NA_REAL;
        continue;
      }
    }
    if(ordinary_hc0 && !defer_hc0_mean_stderr) {
      const double moment = call->mean_stderr[i];
      double se;

      if(!R_FINITE(moment) || moment < 0.0) {
        execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
        return R_NilValue;
      }
      se = call->hc0_context->residual_scale *
        sqrt(moment) / fabs(denominator);
      call->mean_stderr[i] = R_FINITE(se) ? se : NA_REAL;
    } else if(!ordinary_hc0 && call->do_merr) {
      call->mean_stderr[i] =
        owner->mean_columns[response_offset + 2] / denominator -
        call->mean[i] * call->mean[i];
      call->mean_stderr[i] = (call->mean_stderr[i] <= 0.0) ? 0.0 :
        sqrt(call->mean_stderr[i] * call->kernel_squared_integral /
             (denominator * call->bandwidth_product));
      if(!R_FINITE(call->mean_stderr[i]))
        call->mean_stderr[i] = NA_REAL;
    }
  }

  if(call->do_grad) {
    for(predictor = continuous_range.begin;
        predictor < continuous_range.end; predictor++) {
      for(i = 0; i < call->num_obs_eval; i++) {
        const int response_offset =
          response_column_count * i;
        const int permutation_offset =
          predictor * call->num_obs_eval *
          response_column_count + response_offset;
        const int derivative_denominator_column =
          point_already_computed ? 0 : 1;
        const double denominator =
          owner->mean_columns[
            response_offset + derivative_denominator_column];
        const double derivative_denominator =
          owner->permutation_columns[
            permutation_offset + derivative_denominator_column];
        const double derivative_numerator = point_already_computed ? 0.0 :
          owner->permutation_columns[permutation_offset];
        double hc0_cross_moment = 0.0;

        if(!R_FINITE(call->mean[i]) || !R_FINITE(denominator) ||
           denominator == 0.0) {
          call->gradient[predictor][i] = NA_REAL;
          if(call->do_gerr)
            call->gradient_stderr[predictor][i] = NA_REAL;
          continue;
        }
        if(ordinary_hc0 && call->bandwidth_mode == BW_ADAP_NN)
          hc0_cross_moment = call->gradient[predictor][i];
        if(!point_already_computed)
          call->gradient[predictor][i] =
            (derivative_numerator -
             call->mean[i] * derivative_denominator) / denominator;
        if(!R_FINITE(call->gradient[predictor][i])) {
          call->gradient[predictor][i] = NA_REAL;
          if(call->do_gerr)
            call->gradient_stderr[predictor][i] = NA_REAL;
          continue;
        }
        if(call->do_gerr && ordinary_hc0) {
          double quadratic = call->gradient_stderr[predictor][i];

          if(call->hc0_context->gradient_structural_zero != NULL &&
             call->hc0_context->gradient_structural_zero[predictor][i]) {
            quadratic = 0.0;
          } else if(call->bandwidth_mode == BW_ADAP_NN &&
             !np_regression_hc0_derivative_quadratic(
               call->mean_stderr[i], quadratic, hc0_cross_moment,
               derivative_denominator / denominator, &quadratic)) {
            execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
            return R_NilValue;
          }
          if(!R_FINITE(quadratic) || quadratic < 0.0) {
            execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
            return R_NilValue;
          }
          call->gradient_stderr[predictor][i] =
            call->hc0_context->residual_scale * sqrt(quadratic) /
            fabs(denominator);
          if(!R_FINITE(call->gradient_stderr[predictor][i])) {
            execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
            return R_NilValue;
          }
        } else if(call->do_gerr) {
          const double bandwidth =
            call->bandwidth_mode == BW_ADAP_NN ? 1.0 :
            (call->bandwidth_mode == BW_GEN_NN ?
             call->matrix_bandwidth[predictor][i] :
             call->matrix_bandwidth[predictor][0]);
          call->gradient_stderr[predictor][i] =
            call->gradient_stderr_factor * call->mean_stderr[i] / bandwidth;
        }
      }
    }

    for(predictor = unordered_range.begin;
        predictor < unordered_range.end;
        predictor++) {
      for(i = 0; i < call->num_obs_eval; i++) {
        const int response_offset =
          response_column_count * i;
        const int permutation_offset =
          predictor * call->num_obs_eval *
          response_column_count + response_offset;
        const int denominator_column = point_already_computed ? 0 : 1;
        const double level_denominator =
          owner->mean_columns[response_offset + denominator_column];
        const double alternate_denominator =
          owner->permutation_columns[
            permutation_offset + denominator_column];
        double hc0_cross_moment = 0.0;

        if(!R_FINITE(call->mean[i]) || !R_FINITE(level_denominator) ||
           level_denominator == 0.0 ||
           !R_FINITE(alternate_denominator) ||
           alternate_denominator == 0.0) {
          if(call->empty_rows != NULL && R_FINITE(call->mean[i]) &&
             alternate_denominator == 0.0)
            np_empty_x_row_mark(&empty_context, i, predictor, call->empty_rows);
          call->gradient[predictor][i] = NA_REAL;
          if(call->do_gerr)
            call->gradient_stderr[predictor][i] = NA_REAL;
          continue;
        }
        if(ordinary_hc0 && call->bandwidth_mode == BW_ADAP_NN)
          hc0_cross_moment = call->gradient[predictor][i];
        if(!point_already_computed && !stable_scalar_contrast)
          call->gradient[predictor][i] = call->mean[i] -
            owner->permutation_columns[permutation_offset] /
              alternate_denominator;
        if(!R_FINITE(call->gradient[predictor][i])) {
          call->gradient[predictor][i] = NA_REAL;
          if(call->do_gerr)
            call->gradient_stderr[predictor][i] = NA_REAL;
          continue;
        }
        if(call->do_gerr && ordinary_hc0 && stable_scalar_contrast) {
          call->gradient_stderr[predictor][i] *= call->hc0_context->residual_scale;
          if(!R_FINITE(call->gradient_stderr[predictor][i])) {
            execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
            return R_NilValue;
          }
        } else if(call->do_gerr && ordinary_hc0) {
          double quadratic = call->gradient_stderr[predictor][i];

          if(call->hc0_context->gradient_structural_zero != NULL &&
             call->hc0_context->gradient_structural_zero[predictor][i]) {
            quadratic = 0.0;
          } else if(call->bandwidth_mode == BW_ADAP_NN &&
             !np_regression_hc0_derivative_quadratic(
               quadratic, call->mean_stderr[i], hc0_cross_moment,
               level_denominator / alternate_denominator, &quadratic)) {
            execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
            return R_NilValue;
          }
          if(!R_FINITE(quadratic) || quadratic < 0.0) {
            execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
            return R_NilValue;
          }
          call->gradient_stderr[predictor][i] =
            call->hc0_context->residual_scale * sqrt(quadratic) /
            fabs(level_denominator);
          if(!R_FINITE(call->gradient_stderr[predictor][i])) {
            execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
            return R_NilValue;
          }
        } else if(call->do_gerr && call->num_reg_continuous > 0) {
          const double alternate_mean =
            owner->permutation_columns[permutation_offset] /
            alternate_denominator;
          const double variance =
            owner->permutation_columns[permutation_offset + 2] /
              alternate_denominator -
            alternate_mean * alternate_mean;
          const double nonnegative_variance = variance <= 0.0 ? 0.0 : variance;
          call->gradient_stderr[predictor][i] = sqrt(
            call->mean_stderr[i] * call->mean_stderr[i] +
            nonnegative_variance * call->kernel_squared_integral /
            (alternate_denominator * call->bandwidth_product));
        } else if(call->do_gerr) {
          call->gradient_stderr[predictor][i] = 0.0;
        }
      }
    }

    for(predictor = ordered_range.begin;
        predictor < ordered_range.end; predictor++) {
      const int ordered_coordinate = predictor -
        call->num_reg_continuous - call->num_reg_unordered;

      if(call->num_categories[
           call->num_reg_unordered + ordered_coordinate] == 1) {
        for(i = 0; i < call->num_obs_eval; ++i) {
          if(!R_FINITE(call->mean[i])) {
            call->gradient[predictor][i] = NA_REAL;
            if(call->do_gerr)
              call->gradient_stderr[predictor][i] = NA_REAL;
          } else if(!point_already_computed) {
            call->gradient[predictor][i] = 0.0;
          }
          if(call->do_gerr && R_FINITE(call->mean[i]))
            call->gradient_stderr[predictor][i] = 0.0;
        }
        continue;
      }
      for(i = 0; i < call->num_obs_eval; i++) {
        const int response_offset =
          response_column_count * i;
        const int permutation_offset =
          predictor * call->num_obs_eval *
          response_column_count + response_offset;
        const int denominator_column = point_already_computed ? 0 : 1;
        const double level_denominator =
          owner->mean_columns[response_offset + denominator_column];
        const double alternate_denominator =
          owner->permutation_columns[
            permutation_offset + denominator_column];
        double hc0_cross_moment = 0.0;

        if(!R_FINITE(call->mean[i]) || !R_FINITE(level_denominator) ||
           level_denominator == 0.0 ||
           !R_FINITE(alternate_denominator) ||
           alternate_denominator == 0.0) {
          if(call->empty_rows != NULL && R_FINITE(call->mean[i]) &&
             alternate_denominator == 0.0)
            np_empty_x_row_mark(&empty_context, i, predictor, call->empty_rows);
          call->gradient[predictor][i] = NA_REAL;
          if(call->do_gerr)
            call->gradient_stderr[predictor][i] = NA_REAL;
          continue;
        }
        if(ordinary_hc0 && call->bandwidth_mode == BW_ADAP_NN)
          hc0_cross_moment = call->gradient[predictor][i];
        if(stable_scalar_contrast)
          call->gradient[predictor][i] *=
            call->matrix_ordered_indices[ordered_coordinate][i] != 0 ? 1.0 : -1.0;
        else if(!point_already_computed)
          call->gradient[predictor][i] =
            (call->mean[i] -
             owner->permutation_columns[permutation_offset] /
               alternate_denominator) *
            (call->matrix_ordered_indices[ordered_coordinate][i] != 0 ?
             1.0 : -1.0);
        if(!R_FINITE(call->gradient[predictor][i])) {
          call->gradient[predictor][i] = NA_REAL;
          if(call->do_gerr)
            call->gradient_stderr[predictor][i] = NA_REAL;
          continue;
        }
        if(call->do_gerr && ordinary_hc0 && stable_scalar_contrast) {
          call->gradient_stderr[predictor][i] *= call->hc0_context->residual_scale;
          if(!R_FINITE(call->gradient_stderr[predictor][i])) {
            execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
            return R_NilValue;
          }
        } else if(call->do_gerr && ordinary_hc0) {
          double quadratic = call->gradient_stderr[predictor][i];

          if(call->hc0_context->gradient_structural_zero != NULL &&
             call->hc0_context->gradient_structural_zero[predictor][i]) {
            quadratic = 0.0;
          } else if(call->bandwidth_mode == BW_ADAP_NN &&
             !np_regression_hc0_derivative_quadratic(
               quadratic, call->mean_stderr[i], hc0_cross_moment,
               level_denominator / alternate_denominator, &quadratic)) {
            execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
            return R_NilValue;
          }
          if(!R_FINITE(quadratic) || quadratic < 0.0) {
            execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
            return R_NilValue;
          }
          call->gradient_stderr[predictor][i] =
            call->hc0_context->residual_scale * sqrt(quadratic) /
            fabs(level_denominator);
          if(!R_FINITE(call->gradient_stderr[predictor][i])) {
            execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
            return R_NilValue;
          }
        } else if(call->do_gerr && call->num_reg_continuous > 0) {
          const double alternate_mean =
            owner->permutation_columns[permutation_offset] /
            alternate_denominator;
          const double variance =
            owner->permutation_columns[permutation_offset + 2] /
              alternate_denominator -
            alternate_mean * alternate_mean;
          const double nonnegative_variance = variance <= 0.0 ? 0.0 : variance;
          call->gradient_stderr[predictor][i] = sqrt(
            call->mean_stderr[i] * call->mean_stderr[i] +
            nonnegative_variance * call->kernel_squared_integral /
            (alternate_denominator * call->bandwidth_product));
        } else if(call->do_gerr) {
          call->gradient_stderr[predictor][i] = 0.0;
        }
      }
    }
  }

  if(defer_hc0_mean_stderr) {
    for(i = 0; i < call->num_obs_eval; ++i) {
      const int response_offset = response_column_count * i;
      const int denominator_column = point_already_computed ? 0 : 1;
      const double denominator =
        owner->mean_columns[response_offset + denominator_column];
      const double moment = call->mean_stderr[i];

      if(!R_FINITE(call->mean[i]) || !R_FINITE(denominator) ||
         denominator == 0.0) {
        call->mean_stderr[i] = NA_REAL;
        continue;
      }
      if(!R_FINITE(moment) || moment < 0.0) {
        execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
        return R_NilValue;
      }
      call->mean_stderr[i] = call->hc0_context->residual_scale *
        sqrt(moment) / fabs(denominator);
      if(!R_FINITE(call->mean_stderr[i])) {
        execution->status = NP_REGRESSION_SCALAR_FIT_ERR_HC0;
        return R_NilValue;
      }
    }
  }

  if(conditional_se) {
    if(!R_FINITE(call->mean[0])) {
      call->mean_stderr[0] = NA_REAL;
      if(call->do_gerr)
        for(predictor = 0; predictor < p_nvar; ++predictor)
          call->gradient_stderr[predictor][0] = NA_REAL;
    } else if(!np_regression_conditional_influence_finish(
                call->num_obs_train, p_nvar,
                call->do_grad ? call->num_reg_continuous : 0,
                call->do_grad ? call->num_reg_unordered : 0,
                call->vector_Y, owner->conditional_weights,
                owner->conditional_permutation_weights,
                owner->mean_columns, owner->permutation_columns,
                call->mean, call->gradient,
                call->mean_stderr, call->gradient_stderr,
                call->empty_rows != NULL ? call->empty_rows->component_flags : NULL)) {
      execution->status = NP_REGRESSION_SCALAR_FIT_ERR_INFLUENCE;
      return R_NilValue;
    }
  }

  execution->status = NP_REGRESSION_SCALAR_FIT_OK;
  return R_NilValue;
}
