static int np_regression_general_lp_fit(
  const NPRegressionGeneralLPFitCall *call)
{
  NPRegressionGeneralLPFitExecution execution;

  execution.call = call;
  execution.status = NP_REGRESSION_GENERAL_LP_FIT_ERR_ALLOC;
  np_regression_general_lp_fit_owner_init(&execution.owner);
  R_UnwindProtect(
    (call->prepared_ann_row != NULL ||
     (call->bandwidth_mode == BW_ADAP_NN && call->ordinary_response &&
      call->tree_enabled == NP_TREE_TRUE && call->kernel_route == NULL &&
      !int_cker_bound_extern && call->num_reg_continuous > 0 &&
      (call->kernel_c[0] == 4 || call->kernel_c[0] == 8))) ?
      np_regression_general_lp_fit_ann_execute :
      np_regression_general_lp_fit_execute, &execution,
    np_regression_general_lp_fit_owner_cleanup, &execution, NULL);
  return execution.status;
}
