#line 9578 "np.c"
SEXP C_np_regression(SEXP tuno,
                     SEXP tord,
                     SEXP tcon,
                     SEXP ty,
                     SEXP euno,
                     SEXP eord,
                     SEXP econ,
                     SEXP ey,
                     SEXP rbw,
                     SEXP mcv,
                     SEXP padnum,
                     SEXP nconfac,
                     SEXP ncatfac,
                     SEXP mysd,
                     SEXP myopti,
                     SEXP glp_degree,
                     SEXP glp_gradient_order,
                     SEXP glp_bernstein,
                     SEXP glp_basis,
                     SEXP enrow,
                     SEXP ncol,
                     SEXP output_request,
                     SEXP ckerlb,
                     SEXP ckerub)
{
  SEXP tuno_r = R_NilValue, tord_r = R_NilValue, tcon_r = R_NilValue, ty_r = R_NilValue;
  SEXP euno_r = R_NilValue, eord_r = R_NilValue, econ_r = R_NilValue, ey_r = R_NilValue;
  SEXP rbw_r = R_NilValue, mcv_r = R_NilValue, padnum_r = R_NilValue;
  SEXP nconfac_r = R_NilValue, ncatfac_r = R_NilValue, mysd_r = R_NilValue, myopti_i = R_NilValue;
  SEXP degree_i = R_NilValue, gradient_order_i = R_NilValue, ckerlb_r = R_NilValue, ckerub_r = R_NilValue;
  SEXP out = R_NilValue, out_names = R_NilValue;
  SEXP out_mean = R_NilValue, out_merr = R_NilValue, out_g = R_NilValue, out_gerr = R_NilValue, out_xtra = R_NilValue;
  int bern = asInteger(glp_bernstein);
  int basis = asInteger(glp_basis);
  int en = asInteger(enrow);
  int nc = asInteger(ncol);
  int request = asInteger(output_request);
  SEXP gradient_coordinate = PROTECT(getAttrib(output_request, install(".np.gradient.coordinate")));
  NPRegressionGradientRequest gradient_request;
  const NPRegressionGradientRequest *gradient_request_ptr = NULL;
  int variance_unavailable = 0;
  int *gradient_zero_out = NULL;
  SEXP empty_flags = R_NilValue;
  NPRegressionLPEmptyRows empty_rows = {NULL, 0, NULL, NULL, NULL};
  NPRegressionLPEmptyRows *empty_rows_ptr = NULL;
  NPRegressionFitProjection fit_projection = {0};
  NPRegressionFitProjection *fit_projection_ptr = NULL;
  SEXP projection_basis = R_NilValue, projection_rows = R_NilValue;
  int extra_protect = 0;
  int do_merr;
  int do_grad;
  int do_gerr;
  int ncon = 0;
  int nunordered = 0;
  int nordered = 0;
  int num_categorical = 0;
  int num_predictors = 0;
  int num_train = 0;
  int num_eval = 0;
  int train_is_eval = 0;
  int ey_is_ty = 0;
  int categorical_compress = 0;
  double * ckerlb_p = NULL;
  double * ckerub_p = NULL;
  R_xlen_t gsize;
  np_continuous_kernel_descriptor descriptor;
  NPContinuousKernelRoute beta_route;
  NPContinuousKernelDerivativeDiagnostics beta_diagnostics;
  const NPContinuousKernelRoute *active_route = NULL;
  NPContinuousKernelDerivativeDiagnostics *active_diagnostics = NULL;

  if (en < 0) en = 0;
  if (nc < 0) nc = 0;
  if(!np_regression_output_request_valid(request))
    error("C_np_regression: invalid output request");
  do_merr = np_regression_output_requests_errors(request);
  do_grad = np_regression_output_requests_gradients(request);
  do_gerr = np_regression_output_requests_gradient_errors(request);
  gsize = do_grad ? (R_xlen_t)en * (R_xlen_t)nc : 1;

  PROTECT(tuno_r = coerceVector(tuno, REALSXP));
  PROTECT(tord_r = coerceVector(tord, REALSXP));
  PROTECT(tcon_r = coerceVector(tcon, REALSXP));
  PROTECT(ty_r = coerceVector(ty, REALSXP));
  PROTECT(euno_r = coerceVector(euno, REALSXP));
  PROTECT(eord_r = coerceVector(eord, REALSXP));
  PROTECT(econ_r = coerceVector(econ, REALSXP));
  PROTECT(ey_r = coerceVector(ey, REALSXP));
  PROTECT(rbw_r = coerceVector(rbw, REALSXP));
  PROTECT(mcv_r = coerceVector(mcv, REALSXP));
  PROTECT(padnum_r = coerceVector(padnum, REALSXP));
  PROTECT(nconfac_r = coerceVector(nconfac, REALSXP));
  PROTECT(ncatfac_r = coerceVector(ncatfac, REALSXP));
  PROTECT(mysd_r = coerceVector(mysd, REALSXP));
  PROTECT(myopti_i = coerceVector(myopti, INTSXP));
  PROTECT(degree_i = coerceVector(glp_degree, INTSXP));
  PROTECT(gradient_order_i = coerceVector(glp_gradient_order, INTSXP));
  PROTECT(ckerlb_r = coerceVector(ckerlb, REALSXP));
  PROTECT(ckerub_r = coerceVector(ckerub, REALSXP));

  if(XLENGTH(myopti_i) <= REG_CKORDERI)
    error("C_np_regression: continuous-kernel descriptor is missing");
  descriptor = np_continuous_kernel_descriptor_or_error(
    INTEGER(myopti_i)[REG_CKFAMILYI],
    INTEGER(myopti_i)[REG_CKRNEVI],
    INTEGER(myopti_i)[REG_CKORDERI],
    "C_np_regression");

  ncon = (int)INTEGER(myopti_i)[REG_NCONI];
  nunordered = (int)INTEGER(myopti_i)[REG_NUNOI];
  nordered = (int)INTEGER(myopti_i)[REG_NORDI];
  if(gradient_coordinate != R_NilValue) {
    if(TYPEOF(gradient_coordinate) != INTSXP || XLENGTH(gradient_coordinate) != 1 ||
       request != NP_REGRESSION_OUTPUT_FULL || ncon < 0 ||
       nunordered < 0 || nordered < 0 ||
       (R_xlen_t)ncon + nunordered + nordered != nc ||
       INTEGER(gradient_coordinate)[0] <= ncon ||
       INTEGER(gradient_coordinate)[0] > nc)
      error("C_np_regression: invalid categorical gradient coordinate request");
    gradient_request.coordinate = INTEGER(gradient_coordinate)[0] - 1;
    gradient_request_ptr = &gradient_request;
  }
  num_train = (int)INTEGER(myopti_i)[REG_TNOBSI];
  num_eval = (int)INTEGER(myopti_i)[REG_ENOBSI];
  train_is_eval = (int)INTEGER(myopti_i)[REG_TISEI];
  ey_is_ty = (int)INTEGER(myopti_i)[REG_EY];
  if(XLENGTH(ty_r) < num_train)
    error("C_np_regression: training-response buffer is too short");
  /* An absent evaluation response aliases training y only at training
   * locations; external prediction legitimately supplies no response. */
  if(!ey_is_ty && XLENGTH(ey_r) < num_eval)
    error("C_np_regression: evaluation-response buffer is too short");
  categorical_compress = INTEGER(myopti_i)[REG_CATCOMPI];
  if(categorical_compress != 0 && categorical_compress != 1)
    error("C_np_regression: categorical compression must be TRUE or FALSE");
  resolve_bounds_or_default(ckerlb_r, ckerub_r, ncon, &ckerlb_p, &ckerub_p);

  PROTECT(out_mean = allocVector(REALSXP, en));
  PROTECT(out_merr = do_merr ? allocVector(REALSXP, en) : R_NilValue);
  PROTECT(out_g = do_grad ? allocVector(REALSXP, gsize) : R_NilValue);
  PROTECT(out_gerr = do_gerr ? allocVector(REALSXP, gsize) : R_NilValue);
  PROTECT(out_xtra = allocVector(REALSXP, 6));
  if(do_gerr)
    for(R_xlen_t ii = 0; ii < gsize; ii++)
      REAL(out_gerr)[ii] = 0.0;

  if(descriptor.family == NP_CKERNEL_FAMILY_BETA) {
    const int beta_bandwidth_code = INTEGER(myopti_i)[REG_BWI];
    const int beta_lp_engine = INTEGER(myopti_i)[REG_LL];
    const int max_levels = INTEGER(myopti_i)[REG_MLEVI];

    if(nunordered < 0 || nordered < 0 || ncon < 0 ||
       nunordered > INT_MAX - nordered)
      error("C_np_regression: invalid beta categorical dimensions");
    num_categorical = nunordered + nordered;
    if(ncon > INT_MAX - num_categorical)
      error("C_np_regression: invalid beta total dimension");
    num_predictors = ncon + num_categorical;
    if(beta_bandwidth_code != BW_FIXED &&
       beta_bandwidth_code != BW_GEN_NN &&
       beta_bandwidth_code != BW_ADAP_NN)
      error("C_np_regression: invalid beta bandwidth mode");
    if(beta_lp_engine != NP_LP_ENGINE_SCALAR &&
       beta_lp_engine != NP_LP_ENGINE_GENERAL)
      error("C_np_regression: invalid beta regression engine");
    if(beta_lp_engine == NP_LP_ENGINE_GENERAL &&
       (XLENGTH(degree_i) < ncon || XLENGTH(gradient_order_i) < ncon))
      error("C_np_regression: beta LP degree metadata is too short");
    if(do_grad != (INTEGER(myopti_i)[REG_GRAD] != 0))
      error("C_np_regression: inconsistent beta gradient flags");
    if(num_train <= 0 || num_eval <= 0 || en != num_eval ||
       ncon <= 0 || nc != num_predictors ||
       (train_is_eval != 0 && train_is_eval != 1) ||
       (ey_is_ty != 0 && ey_is_ty != 1))
      error("C_np_regression: invalid beta regression dimensions");
    if(num_categorical > 0 &&
       (max_levels <= 0 ||
        !np_real_buffer_has_matrix(
          mcv_r, max_levels, num_categorical)))
      error("C_np_regression: beta categorical level buffer is too short");
    if(!np_real_buffer_has_matrix(tcon_r, num_train, ncon) ||
       (!train_is_eval &&
        !np_real_buffer_has_matrix(econ_r, num_eval, ncon)) ||
       !np_real_buffer_has_matrix(tuno_r, num_train, nunordered) ||
       !np_real_buffer_has_matrix(tord_r, num_train, nordered) ||
       (!train_is_eval &&
        (!np_real_buffer_has_matrix(euno_r, num_eval, nunordered) ||
         !np_real_buffer_has_matrix(eord_r, num_eval, nordered))) ||
       XLENGTH(rbw_r) < num_predictors ||
       XLENGTH(padnum_r) < 1 || XLENGTH(nconfac_r) < 1 ||
       XLENGTH(ncatfac_r) < 1 || XLENGTH(mysd_r) < ncon ||
       XLENGTH(ckerlb_r) < ncon || XLENGTH(ckerub_r) < ncon)
      error("C_np_regression: beta regression input buffer is too short");
    if(train_is_eval && num_eval != num_train)
      error("C_np_regression: beta train/evaluation dimensions are inconsistent");

    beta_route.segment_count = 1;
    beta_route.segment[0].descriptor = descriptor;
    beta_route.segment[0].coordinate_offset = 0;
    beta_route.segment[0].coordinate_count = ncon;
    beta_route.segment[0].lower = ckerlb_p;
    beta_route.segment[0].upper = ckerub_p;
    beta_diagnostics.bad_coordinate = -1;
    beta_diagnostics.bad_observation = -1;
    beta_diagnostics.undefined_count = 0;
    beta_diagnostics.beta_status = NP_BETA_OK;
    active_route = &beta_route;
    active_diagnostics = &beta_diagnostics;
  }

  if(TYPEOF(output_request) == INTSXP && XLENGTH(output_request) == 3) {
    const int width = INTEGER(output_request)[2];
    if(width < 1 || do_merr || do_grad || ncon < 1 || num_train < 1 ||
       num_eval != en || INTEGER(myopti_i)[REG_LL] != NP_LP_ENGINE_GENERAL ||
       (double)width * ((double)num_train + en) > 8 * 1024 * 1024)
      error("invalid or oversized plot fit projection request");
#ifdef MPI2
    if(!np_mpi_local_regression_active() && iNum_Processors > 1)
      error("plot fit projection requires an explicit local regression scope");
#endif
    PROTECT(projection_basis = allocMatrix(REALSXP, num_train, width));
    PROTECT(projection_rows = allocMatrix(REALSXP, en, width));
    extra_protect += 2;
    fit_projection = (NPRegressionFitProjection){
      .ntrain = num_train, .neval = en, .width = width,
      .basis = REAL(projection_basis), .projection = REAL(projection_rows)};
    fit_projection_ptr = &fit_projection;
  }

  if(TYPEOF(output_request) == INTSXP && XLENGTH(output_request) >= 2 &&
     INTEGER(output_request)[1] == 1 && !train_is_eval) {
    PROTECT(empty_flags = allocVector(INTSXP, en));
    ++extra_protect;
    memset(INTEGER(empty_flags), 0, (size_t)en*sizeof(int));
    empty_rows.flags = INTEGER(empty_flags);
    empty_rows_ptr = &empty_rows;
  }

  np_regression(REAL(tuno_r), REAL(tord_r), REAL(tcon_r), REAL(ty_r),
                REAL(euno_r), REAL(eord_r), REAL(econ_r), REAL(ey_r),
                REAL(rbw_r), REAL(mcv_r), REAL(padnum_r),
                REAL(nconfac_r), REAL(ncatfac_r), REAL(mysd_r),
                INTEGER(myopti_i),
                INTEGER(degree_i), INTEGER(gradient_order_i), &bern, &basis,
                REAL(out_mean), do_merr ? REAL(out_merr) : NULL,
                do_grad ? REAL(out_g) : NULL,
                do_gerr ? REAL(out_gerr) : NULL,
                REAL(out_xtra), ckerlb_p, ckerub_p,
                active_route, active_diagnostics, categorical_compress, empty_rows_ptr,
                gradient_request_ptr, &variance_unavailable, &gradient_zero_out, fit_projection_ptr);

  if(gradient_request_ptr != NULL) {
    for(int column = 0; column < nc; ++column)
      if(column != gradient_request.coordinate)
        for(int row = 0; row < en; ++row) {
          REAL(out_g)[(R_xlen_t)column*en + row] = NA_REAL;
          REAL(out_gerr)[(R_xlen_t)column*en + row] = NA_REAL;
        }
  }

  PROTECT(out = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, out_mean);
  SET_VECTOR_ELT(out, 1, out_merr);
  SET_VECTOR_ELT(out, 2, out_g);
  SET_VECTOR_ELT(out, 3, out_gerr);
  SET_VECTOR_ELT(out, 4, out_xtra);

  PROTECT(out_names = allocVector(STRSXP, 5));
  SET_STRING_ELT(out_names, 0, mkChar("mean"));
  SET_STRING_ELT(out_names, 1, mkChar("merr"));
  SET_STRING_ELT(out_names, 2, mkChar("g"));
  SET_STRING_ELT(out_names, 3, mkChar("gerr"));
  SET_STRING_ELT(out_names, 4, mkChar("xtra"));
  setAttrib(out, R_NamesSymbol, out_names);

  if(empty_rows.count > 0)
    setAttrib(out, install(".np.empty.rows"), empty_flags);
  if(variance_unavailable)
    setAttrib(out, install(".np.variance.unavailable"), ScalarLogical(TRUE));
  if(gradient_zero_out != NULL) {
    SEXP certificate = PROTECT(allocVector(LGLSXP, gsize));
    ++extra_protect;
    memcpy(LOGICAL(certificate), gradient_zero_out, (size_t)gsize*sizeof(int));
    setAttrib(out, install(".np.gradient.structural.zero"), certificate);
  }
  if(fit_projection_ptr != NULL) {
    if(fit_projection.emitted != en)
      error("plot fit projection did not visit every fitted row");
    SEXP projection_result = PROTECT(allocVector(VECSXP, 3));
    ++extra_protect;
    SET_VECTOR_ELT(projection_result, 0, projection_basis);
    SET_VECTOR_ELT(projection_result, 1, projection_rows);
    SET_VECTOR_ELT(projection_result, 2, ScalarLogical(fit_projection.constant));
    setAttrib(out, install(".np.fit.projection"), projection_result);
  }
  UNPROTECT(27 + extra_protect);
  return out;
}
