#line 13631 "np.c"
SEXP C_np_kernelsum(SEXP tuno,
                    SEXP tord,
                    SEXP tcon,
                    SEXP ty,
                    SEXP weights,
                    SEXP euno,
                    SEXP eord,
                    SEXP econ,
                    SEXP bw,
                    SEXP mcv,
                    SEXP padnum,
                    SEXP op,
                    SEXP myopti,
                    SEXP kpow,
                    SEXP ksum_len,
                    SEXP pksum_len,
                    SEXP kw_len,
                    SEXP ckerlb,
                    SEXP ckerub)
{
  SEXP tuno_r=R_NilValue, tord_r=R_NilValue, tcon_r=R_NilValue, ty_r=R_NilValue, weights_r=R_NilValue;
  SEXP euno_r=R_NilValue, eord_r=R_NilValue, econ_r=R_NilValue, bw_r=R_NilValue, mcv_r=R_NilValue;
  SEXP padnum_r=R_NilValue, op_i=R_NilValue, myopti_i=R_NilValue, kpow_r=R_NilValue, ckerlb_r=R_NilValue, ckerub_r=R_NilValue;
  SEXP out=R_NilValue, out_names=R_NilValue, out_ksum=R_NilValue, out_pksum=R_NilValue, out_kw=R_NilValue, out_pkw=R_NilValue;
  int n_ksum = asInteger(ksum_len);
  int n_pksum = asInteger(pksum_len);
  int n_kw = asInteger(kw_len);
  int ncon = 0, nuno = 0, nord = 0, p_operator = 0, do_score = 0, do_ocg = 0, return_kernel_weights = 0, p_nvar = 0;
  R_xlen_t n_pkw = 0;
  const int *eval_to_train = NULL;
  double * ckerlb_p = NULL;
  double * ckerub_p = NULL;
  np_continuous_kernel_descriptor descriptor;

  if(n_ksum < 0) n_ksum = 0;
  if(n_pksum < 0) n_pksum = 0;
  if(n_kw < 0) n_kw = 0;

  PROTECT(tuno_r = coerceVector(tuno, REALSXP));
  PROTECT(tord_r = coerceVector(tord, REALSXP));
  PROTECT(tcon_r = coerceVector(tcon, REALSXP));
  PROTECT(ty_r = coerceVector(ty, REALSXP));
  PROTECT(weights_r = coerceVector(weights, REALSXP));
  PROTECT(euno_r = coerceVector(euno, REALSXP));
  PROTECT(eord_r = coerceVector(eord, REALSXP));
  PROTECT(econ_r = coerceVector(econ, REALSXP));
  PROTECT(bw_r = coerceVector(bw, REALSXP));
  PROTECT(mcv_r = coerceVector(mcv, REALSXP));
  PROTECT(padnum_r = coerceVector(padnum, REALSXP));
  PROTECT(op_i = coerceVector(op, INTSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(kpow_r = coerceVector(kpow, REALSXP));
  PROTECT(ckerlb_r = coerceVector(ckerlb, REALSXP));
  PROTECT(ckerub_r = coerceVector(ckerub, REALSXP));

  if(XLENGTH(myopti_i) < KWS_OPTIONS_COUNT)
    error("C_np_kernelsum: invalid internal option vector");
  if(XLENGTH(kpow_r) != 1 || !R_FINITE(REAL(kpow_r)[0]) ||
     REAL(kpow_r)[0] != floor(REAL(kpow_r)[0]) ||
     fabs(REAL(kpow_r)[0]) > INT_MAX)
    error("C_np_kernelsum: kernel.pow must be one finite integer");

  descriptor = np_kernelsum_descriptor_or_error(myopti_i, "C_np_kernelsum");

  ncon = (int)INTEGER(myopti_i)[KWS_NCONI];
  nuno = (int)INTEGER(myopti_i)[KWS_NUNOI];
  nord = (int)INTEGER(myopti_i)[KWS_NORDI];
  p_operator = (int)INTEGER(myopti_i)[KWS_POPI];
  do_score = (int)INTEGER(myopti_i)[KWS_PSCOREI];
  do_ocg = (int)INTEGER(myopti_i)[KWS_POCGI];
  return_kernel_weights = (int)INTEGER(myopti_i)[KWS_RKWI];
  p_nvar = ((p_operator != OP_NOOP) ? ncon : 0) + ((do_score || do_ocg) ? (nuno + nord) : 0);
  n_pkw = (return_kernel_weights && (p_nvar > 0)) ? ((R_xlen_t)n_kw * (R_xlen_t)p_nvar) : 0;
  resolve_bounds_or_default(ckerlb_r, ckerub_r, ncon, &ckerlb_p, &ckerub_p);

  const int fold_requested =
    getAttrib(myopti, install("np.fold.train.index")) != R_NilValue;
  if(INTEGER(myopti_i)[KWS_BWI] == BW_GEN_NN || fold_requested) {
    SEXP map = getAttrib(myopti, install(fold_requested ?
                        "np.fold.train.index" : "np.eval.train.index"));
    if(map != R_NilValue) {
      const int nt = INTEGER(myopti_i)[KWS_TNOBSI];
      const int ne = INTEGER(myopti_i)[KWS_ENOBSI];
      if(TYPEOF(map) != INTSXP || nt <= 0 || ne <= 0 ||
         XLENGTH(map) != ne || ncon <= 0 ||
         (!fold_requested && INTEGER(myopti_i)[KWS_TISEI]) ||
         INTEGER(myopti_i)[KWS_LOOI] ||
         (fold_requested && INTEGER(myopti_i)[KWS_BWI] != BW_GEN_NN &&
          INTEGER(myopti_i)[KWS_BWI] != BW_ADAP_NN) ||
         !np_real_buffer_has_matrix(tcon_r, nt, ncon) ||
         (!INTEGER(myopti_i)[KWS_TISEI] &&
          !np_real_buffer_has_matrix(econ_r, ne, ncon)))
        error("C_np_kernelsum: invalid internal NN training map");
      eval_to_train = INTEGER(map);
      for(int j = 0; j < ne; ++j)
        if(eval_to_train[j] < 0 || eval_to_train[j] >= nt)
          error("C_np_kernelsum: invalid internal NN training index");
    }
  }

  PROTECT(out_ksum = allocVector(REALSXP, n_ksum));
  PROTECT(out_pksum = allocVector(REALSXP, n_pksum));
  PROTECT(out_kw = allocVector(REALSXP, n_kw));
  PROTECT(out_pkw = allocVector(REALSXP, n_pkw));

  if(fold_requested) {
    if(p_operator != OP_NOOP || do_score || do_ocg ||
       (REAL(kpow_r)[0] != 1.0 && REAL(kpow_r)[0] != 2.0) ||
       XLENGTH(op_i) != ncon + nuno + nord)
      error("invalid operator for raw delete-one kernel weights");
    for(int d = 0; d < ncon + nuno + nord; ++d)
      if(INTEGER(op_i)[d] != OP_NORMAL)
        error("raw delete-one kernel weights require normal operators");
  }

  if(descriptor.family == NP_CKERNEL_FAMILY_BETA) {
    const int num_train = INTEGER(myopti_i)[KWS_TNOBSI];
    const int num_eval = INTEGER(myopti_i)[KWS_ENOBSI];
    const int num_response_columns = INTEGER(myopti_i)[KWS_YNCOLI];
    const int num_weight_columns = INTEGER(myopti_i)[KWS_WNCOLI];
    const int train_is_eval = INTEGER(myopti_i)[KWS_TISEI];
    const int leave_one_out = INTEGER(myopti_i)[KWS_LOOI];
    const int categorical_compress =
      XLENGTH(myopti_i) > KWS_CCOMPRESSI ?
      INTEGER(myopti_i)[KWS_CCOMPRESSI] : 0;
    const int response_extent = (num_response_columns > 0) ?
      num_response_columns : 1;
    const int weight_extent = (num_weight_columns > 0) ?
      num_weight_columns : 1;
    const R_xlen_t expected_sum = (R_xlen_t)num_eval *
      (R_xlen_t)response_extent * (R_xlen_t)weight_extent;
    const R_xlen_t expected_weights = (R_xlen_t)num_train *
      (R_xlen_t)num_eval;
    np_beta_bandwidth_mode beta_bandwidth_mode = NP_BETA_BANDWIDTH_FIXED;
    const double *beta_bandwidth_eval = REAL(bw_r);
    const double *beta_bandwidth_train = REAL(bw_r);
    int beta_bandwidth_code = INTEGER(myopti_i)[KWS_BWI];
    const int beta_kernel_power = (int)REAL(kpow_r)[0];
    int has_overlap = 0;
    int derivative_dimension = -1;
    int undefined_count = 0;
    int ncat;
    int i;

    if(ncon < 0 || nuno < 0 || nord < 0 ||
       nuno > INT_MAX - nord)
      error("C_np_kernelsum: invalid beta kernel-sum dimensions");
    ncat = nuno + nord;
    if(ncat > 0 && XLENGTH(myopti_i) <= KWS_CCOMPRESSI)
      error("C_np_kernelsum: categorical-compression state is missing");
    if(categorical_compress != 0 && categorical_compress != 1)
      error("C_np_kernelsum: invalid categorical-compression option");
    if(beta_bandwidth_code != BW_FIXED &&
       beta_bandwidth_code != BW_GEN_NN &&
       beta_bandwidth_code != BW_ADAP_NN)
      error("C_np_kernelsum: invalid beta bandwidth mode");
    if((p_operator != OP_NOOP && p_operator != OP_DERIVATIVE) ||
       do_score || do_ocg)
      error("C_np_kernelsum: beta kernels support only derivative permutation operators");
    if(ncon > INT_MAX - ncat || XLENGTH(op_i) != ncon + ncat)
      error("C_np_kernelsum: beta operator vector has the wrong length");
    for(i = 0; i < ncon; ++i) {
      const int operator_code = INTEGER(op_i)[i];

      if(operator_code == OP_CONVOLUTION) {
        has_overlap = 1;
      } else if(operator_code == OP_DERIVATIVE) {
        if(derivative_dimension >= 0)
          error("C_np_kernelsum: beta kernels support one direct derivative dimension at a time");
        derivative_dimension = i;
      } else if(operator_code != OP_NORMAL && operator_code != OP_INTEGRAL) {
        error("C_np_kernelsum: unsupported beta operator");
      }
    }
    for(i = ncon; i < ncon + ncat; ++i)
      if(INTEGER(op_i)[i] != OP_NORMAL &&
         INTEGER(op_i)[i] != OP_CONVOLUTION &&
         INTEGER(op_i)[i] != OP_INTEGRAL)
        error("C_np_kernelsum: unsupported categorical operator");
    if(derivative_dimension >= 0 && p_operator == OP_DERIVATIVE)
      error("C_np_kernelsum: direct and permutation beta derivatives cannot be combined");
    if(num_train <= 0 || num_eval <= 0 || ncon <= 0 ||
       (train_is_eval != 0 && train_is_eval != 1) ||
       (train_is_eval && num_train != num_eval) ||
       num_response_columns < 0 || num_weight_columns < 0)
      error("C_np_kernelsum: invalid beta kernel-sum dimensions");
    if(!np_real_buffer_has_matrix(tcon_r, num_train, ncon) ||
       (!train_is_eval && !np_real_buffer_has_matrix(
          econ_r, num_eval, ncon)) ||
       !np_real_buffer_has_matrix(tuno_r, num_train, nuno) ||
       !np_real_buffer_has_matrix(tord_r, num_train, nord) ||
       (!train_is_eval &&
        (!np_real_buffer_has_matrix(euno_r, num_eval, nuno) ||
         !np_real_buffer_has_matrix(eord_r, num_eval, nord))) ||
       XLENGTH(bw_r) < ncon + ncat || XLENGTH(ckerlb_r) < ncon ||
       XLENGTH(ckerub_r) < ncon ||
       !np_real_buffer_has_matrix(ty_r, num_train,
                                  num_response_columns) ||
       !np_real_buffer_has_matrix(weights_r, num_train,
                                  num_weight_columns))
      error("C_np_kernelsum: beta kernel-sum input buffer is too short");
    if(ncat > 0) {
      const int max_levels = INTEGER(myopti_i)[KWS_MLEVI];

      if(max_levels <= 0 || XLENGTH(padnum_r) < 1 ||
         !R_FINITE(REAL(padnum_r)[0]) ||
         !np_real_buffer_has_matrix(mcv_r, max_levels, ncat))
        error("C_np_kernelsum: invalid categorical support metadata");
    }
    if(beta_bandwidth_code != BW_FIXED && !fold_requested) {
      const int need_eval = (beta_bandwidth_code == BW_GEN_NN) || has_overlap;
      const int need_train = (beta_bandwidth_code == BW_ADAP_NN) || has_overlap;
      double *bandwidth_eval_storage = need_eval ?
        (double *)R_alloc((size_t)ncon * (size_t)num_eval,
                          sizeof(double)) : NULL;
      double *bandwidth_train_storage = need_train ?
        (double *)R_alloc((size_t)ncon * (size_t)num_train,
                          sizeof(double)) : NULL;
      np_beta_bandwidth_prepare_status bandwidth_status;

      beta_bandwidth_mode = (beta_bandwidth_code == BW_GEN_NN) ?
        NP_BETA_BANDWIDTH_GENERALIZED_NN : NP_BETA_BANDWIDTH_ADAPTIVE_NN;
      if(eval_to_train != NULL) {
        if(need_train || beta_bandwidth_code != BW_GEN_NN)
          error("C_np_kernelsum: mapped beta queries require ordinary GNN rows");
        /* Private beta GNN hat queries retain observation identity through
         * evaluation chunks. Public kernel sums keep their external-query
         * radius contract; fixed and ANN owners never enter this branch. */
        NPNNGeometryContext geometry = {
          .mode = NP_NN_QUERY_TRAINING_MAP,
          .eval_to_train = eval_to_train,
          .adaptive_successor = NULL
        };
        double **train_columns = (double **)R_alloc((size_t)ncon, sizeof(double *));
        double **eval_columns = (double **)R_alloc((size_t)ncon, sizeof(double *));
        double **bandwidth_columns = (double **)R_alloc((size_t)ncon, sizeof(double *));
        for(int d = 0; d < ncon; ++d) {
          train_columns[d] = REAL(tcon_r) + (size_t)d*num_train;
          eval_columns[d] = (train_is_eval ? REAL(tcon_r) : REAL(econ_r)) +
            (size_t)d*num_eval;
          bandwidth_columns[d] = bandwidth_eval_storage + (size_t)d*num_eval;
        }
        bandwidth_status = np_beta_bandwidth_prepare_matrix_ctx(
          beta_bandwidth_mode, train_columns, eval_columns, REAL(bw_r),
          num_train, num_eval, ncon, 0, 1, 0, 0,
          bandwidth_columns, NULL, &geometry);
      } else {
        bandwidth_status = np_beta_bandwidth_prepare(
          beta_bandwidth_mode,
          REAL(tcon_r), train_is_eval ? NULL : REAL(econ_r), REAL(bw_r),
          num_train, num_eval, ncon, train_is_eval,
          need_eval, need_train, 0,
          bandwidth_eval_storage, bandwidth_train_storage);
      }
      if(bandwidth_status == NP_BETA_BANDWIDTH_PREPARE_ERR_ZERO_RADIUS) {
        const NPNNGeometryContext failure_geometry = {
          .mode = NP_NN_QUERY_TRAINING_MAP,
          .eval_to_train = eval_to_train,
          .adaptive_successor = NULL
        };
        const NPNNZeroRadiusInfo info = np_nn_zero_radius_info_flat(
          beta_bandwidth_code, num_train, num_eval, ncon,
          REAL(tcon_r), train_is_eval ? NULL : REAL(econ_r), REAL(bw_r),
          eval_to_train != NULL ? &failure_geometry : NULL, 0);
        np_nn_zero_radius_error(&info);
      }
      if(bandwidth_status != NP_BETA_BANDWIDTH_PREPARE_OK)
        error("C_np_kernelsum: %s",
              np_beta_bandwidth_prepare_status_message(bandwidth_status));

      beta_bandwidth_eval = bandwidth_eval_storage;
      beta_bandwidth_train = bandwidth_train_storage;
    }
    if(n_ksum != expected_sum)
      error("C_np_kernelsum: beta kernel-sum output buffer has the wrong length");
    if(return_kernel_weights && n_kw != expected_weights)
      error("C_np_kernelsum: beta kernel-weight buffer has the wrong length");
    if(p_operator == OP_DERIVATIVE &&
       n_pksum != (R_xlen_t)ncon * expected_sum)
      error("C_np_kernelsum: beta derivative-sum buffer has the wrong length");
    if(p_operator == OP_NOOP && n_pksum != 0)
      error("C_np_kernelsum: unexpected beta permutation-sum buffer");
    if(p_operator == OP_DERIVATIVE && return_kernel_weights &&
       n_pkw != (R_xlen_t)ncon * expected_weights)
      error("C_np_kernelsum: beta derivative-weight buffer has the wrong length");

    {
      NPContinuousKernelRoute route;
      NPContinuousKernelDerivativeDiagnostics route_diagnostics;
      double **train_columns = (double **)R_alloc(
        (size_t)ncon, sizeof(double *));
      double **evaluation_columns = (double **)R_alloc(
        (size_t)ncon, sizeof(double *));
      double **bandwidth_eval_columns = NULL;
      double **bandwidth_train_columns = NULL;
      double **response_columns = NULL;
      double **weight_columns = NULL;
      NPBetaCategoricalIngress categorical_ingress;
      int route_status;
      NPNNKernelFold fold = {0};

      route.segment_count = 1;
      route.segment[0].descriptor = descriptor;
      route.segment[0].coordinate_offset = 0;
      route.segment[0].coordinate_count = ncon;
      route.segment[0].lower = ckerlb_p;
      route.segment[0].upper = ckerub_p;
      for(i = 0; i < ncon; ++i) {
        train_columns[i] = REAL(tcon_r) + (size_t)i * (size_t)num_train;
        evaluation_columns[i] = train_is_eval ? train_columns[i] :
          REAL(econ_r) + (size_t)i * (size_t)num_eval;
      }
      np_beta_categorical_ingress_prepare_or_error(
        &categorical_ingress, tuno_r, tord_r, euno_r, eord_r,
        mcv_r, padnum_r, myopti_i, num_train, num_eval,
        nuno, nord, train_is_eval, "C_np_kernelsum");
      if(beta_bandwidth_eval != NULL &&
         beta_bandwidth_code != BW_FIXED) {
        bandwidth_eval_columns = (double **)R_alloc(
          (size_t)ncon, sizeof(double *));
        for(i = 0; i < ncon; ++i)
          bandwidth_eval_columns[i] = (double *)beta_bandwidth_eval +
            (size_t)i * (size_t)num_eval;
      }
      if(beta_bandwidth_train != NULL &&
         beta_bandwidth_code != BW_FIXED) {
        bandwidth_train_columns = (double **)R_alloc(
          (size_t)ncon, sizeof(double *));
        for(i = 0; i < ncon; ++i)
          bandwidth_train_columns[i] = (double *)beta_bandwidth_train +
            (size_t)i * (size_t)num_train;
      }
      if(num_response_columns > 0) {
        response_columns = (double **)R_alloc(
          (size_t)num_response_columns, sizeof(double *));
        for(i = 0; i < num_response_columns; ++i)
          response_columns[i] = REAL(ty_r) +
            (size_t)i * (size_t)num_train;
      }
      if(num_weight_columns > 0) {
        weight_columns = (double **)R_alloc(
          (size_t)num_weight_columns, sizeof(double *));
        for(i = 0; i < num_weight_columns; ++i)
          weight_columns[i] = REAL(weights_r) +
            (size_t)i * (size_t)num_train;
      }

      if(fold_requested) {
        if(np_kernelsum_fold_prepare(&fold, beta_bandwidth_code,
          num_train, num_eval, ncon, 0, 0, train_columns, evaluation_columns,
          eval_to_train, REAL(bw_r), NULL, 1, 0) != 0)
          error("raw delete-one kernel geometry failed");
        if(beta_bandwidth_code == BW_GEN_NN)
          bandwidth_eval_columns = fold.primary;
        else
          bandwidth_train_columns = fold.selected;
      }
      route_status = kernel_weighted_sum_np_fold_route(
        NULL, categorical_ingress.kernel_unordered,
        categorical_ingress.kernel_ordered,
        beta_bandwidth_code, num_train, num_eval,
        nuno, nord, ncon, leave_one_out, 0,
        beta_kernel_power, 0, 0, 0, 0, 0, 0,
        INTEGER(op_i), OP_NOOP, 0, 0, NULL, 0,
        num_response_columns, num_weight_columns,
        NP_TREE_FALSE, 0, NULL, NULL, NULL, NULL,
        categorical_ingress.train_unordered,
        categorical_ingress.train_ordered, train_columns,
        categorical_ingress.evaluation_unordered,
        categorical_ingress.evaluation_ordered,
        evaluation_columns,
        response_columns, weight_columns, NULL,
        (beta_bandwidth_code == BW_FIXED) ? REAL(bw_r) : NULL,
        beta_bandwidth_code != BW_FIXED,
        bandwidth_train_columns, bandwidth_eval_columns,
        ncat > 0 ? REAL(bw_r) + ncon : NULL,
        categorical_ingress.num_categories,
        categorical_ingress.category_values, NULL,
        REAL(out_ksum), REAL(out_pksum),
        return_kernel_weights ? REAL(out_kw) : NULL, NULL,
        categorical_compress, &route,
        &route_diagnostics, fold_requested ? &fold : NULL);
      undefined_count = route_diagnostics.undefined_count;
      if(route_status != 0 &&
         route_diagnostics.beta_status != NP_BETA_OK)
        error("C_np_kernelsum: beta scalar operator failed in continuous dimension %d: %s",
              route_diagnostics.bad_coordinate + 1,
              np_beta_status_message(route_diagnostics.beta_status));
      if(route_status != 0)
        error("C_np_kernelsum: canonical beta kernel row failed with code %d",
              route_status);
      if(p_operator == OP_DERIVATIVE) {
        int *direct_operators = (int *)R_alloc(
          (size_t)(ncon + ncat), sizeof(int));

        memcpy(direct_operators, INTEGER(op_i),
               (size_t)(ncon + ncat) * sizeof(int));
        for(i = 0; i < ncon; ++i) {
          if(i > 0)
            direct_operators[i - 1] = INTEGER(op_i)[i - 1];
          direct_operators[i] = OP_DERIVATIVE;
          route_status = kernel_weighted_sum_np_route(
            NULL, categorical_ingress.kernel_unordered,
            categorical_ingress.kernel_ordered,
            beta_bandwidth_code, num_train, num_eval,
            nuno, nord, ncon, leave_one_out, 0,
            beta_kernel_power, 0, 0, 0, 0, 0, 0,
            direct_operators, OP_NOOP, 0, 0, NULL, 0,
            num_response_columns, num_weight_columns,
            NP_TREE_FALSE, 0, NULL, NULL, NULL, NULL,
            categorical_ingress.train_unordered,
            categorical_ingress.train_ordered, train_columns,
            categorical_ingress.evaluation_unordered,
            categorical_ingress.evaluation_ordered,
            evaluation_columns,
            response_columns, weight_columns, NULL,
            (beta_bandwidth_code == BW_FIXED) ? REAL(bw_r) : NULL,
            beta_bandwidth_code != BW_FIXED,
            bandwidth_train_columns, bandwidth_eval_columns,
            ncat > 0 ? REAL(bw_r) + ncon : NULL,
            categorical_ingress.num_categories,
            categorical_ingress.category_values, NULL,
            REAL(out_pksum) + (R_xlen_t)i * expected_sum, NULL,
            return_kernel_weights ?
              REAL(out_pkw) + (R_xlen_t)i * expected_weights : NULL,
            NULL, categorical_compress, &route, &route_diagnostics);
          undefined_count += route_diagnostics.undefined_count;
          if(route_status != 0 &&
             route_diagnostics.beta_status != NP_BETA_OK)
            error("C_np_kernelsum: beta derivative failed in continuous dimension %d: %s",
                  route_diagnostics.bad_coordinate + 1,
                  np_beta_status_message(route_diagnostics.beta_status));
          if(route_status != 0)
            error("C_np_kernelsum: canonical beta derivative row failed with code %d",
                  route_status);
        }
      }
    }

    if(derivative_dimension >= 0 || p_operator == OP_DERIVATIVE) {
      int infinite_count = 0;
      R_xlen_t aggregate_index;
      const R_xlen_t aggregate_length =
        (R_xlen_t)n_ksum + (R_xlen_t)n_pksum;
      for(aggregate_index = 0; aggregate_index < aggregate_length;
          ++aggregate_index) {
        const double value = (aggregate_index < n_ksum) ?
          REAL(out_ksum)[aggregate_index] :
          REAL(out_pksum)[aggregate_index - n_ksum];
        if(!ISNA(value) && !R_FINITE(value))
          ++infinite_count;
      }
      if(infinite_count > 0 || undefined_count > 0)
        warning("beta derivative produced %d infinite endpoint value(s) and %d undefined cancellation(s)",
                infinite_count, undefined_count);
    }
  } else {
    np_kernelsum_common(REAL(tuno_r), REAL(tord_r), REAL(tcon_r),
                 REAL(ty_r), REAL(weights_r),
                 REAL(euno_r), REAL(eord_r), REAL(econ_r),
                 REAL(bw_r),
                 REAL(mcv_r), REAL(padnum_r),
                 INTEGER(op_i), INTEGER(myopti_i), REAL(kpow_r),
                 REAL(out_ksum), NULL, REAL(out_pksum), REAL(out_kw), REAL(out_pkw),
                 ckerlb_p, ckerub_p, eval_to_train, fold_requested);
  }

  PROTECT(out = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, out_ksum);
  SET_VECTOR_ELT(out, 1, out_pksum);
  SET_VECTOR_ELT(out, 2, out_kw);
  SET_VECTOR_ELT(out, 3, out_pkw);

  PROTECT(out_names = allocVector(STRSXP, 4));
  SET_STRING_ELT(out_names, 0, mkChar("ksum"));
  SET_STRING_ELT(out_names, 1, mkChar("p.ksum"));
  SET_STRING_ELT(out_names, 2, mkChar("kernel.weights"));
  SET_STRING_ELT(out_names, 3, mkChar("p.kernel.weights"));
  setAttrib(out, R_NamesSymbol, out_names);

  UNPROTECT(22);
  return out;
}
