static int np_regression_scalar_fit(
  const NPRegressionScalarFitCall *call)
{
  NPRegressionScalarFitExecution execution;

  execution.call = call;
  execution.status = NP_REGRESSION_SCALAR_FIT_ERR_ALLOC;
  np_regression_scalar_fit_owner_init(
    &execution.owner, call->enclosing_owner);
  R_UnwindProtect(
    np_regression_scalar_fit_execute, &execution,
    np_regression_scalar_fit_owner_cleanup, &execution.owner, NULL);
  return execution.status;
}
