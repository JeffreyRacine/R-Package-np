static SEXP np_regression_categorical_profile_fit_execute(void *data)
{
  NPRegressionCategoricalProfileFitExecution * const execution =
    (NPRegressionCategoricalProfileFitExecution *)data;

  execution->status = np_regression_categorical_profile_fit_body(
    execution->call, &execution->owner);
  return R_NilValue;
}
