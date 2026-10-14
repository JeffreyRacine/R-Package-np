static int np_regression_categorical_profile_fit(
int *kernel_c,
int *kernel_u,
int *kernel_o,
const int lp_engine,
const int BANDWIDTH_reg,
const int num_obs_train,
const int num_obs_eval,
const int num_reg_unordered,
const int num_reg_ordered,
const int num_reg_continuous,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double *vector_Y,
double *vector_scale_factor,
double *lambda,
int *num_categories,
int *operator,
double **matrix_categorical_vals,
double *mean,
double *mean_stderr,
const int emit_fit_progress,
NPRegressionResidualPreparation *residual_preparation)
{
  const NPRegressionCategoricalProfileFitCall call = {
    .kernel_c = kernel_c,
    .kernel_u = kernel_u,
    .kernel_o = kernel_o,
    .lp_engine = lp_engine,
    .bandwidth_mode = BANDWIDTH_reg,
    .num_obs_train = num_obs_train,
    .num_obs_eval = num_obs_eval,
    .num_reg_unordered = num_reg_unordered,
    .num_reg_ordered = num_reg_ordered,
    .num_reg_continuous = num_reg_continuous,
    .matrix_X_unordered_train = matrix_X_unordered_train,
    .matrix_X_ordered_train = matrix_X_ordered_train,
    .matrix_X_unordered_eval = matrix_X_unordered_eval,
    .matrix_X_ordered_eval = matrix_X_ordered_eval,
    .vector_Y = vector_Y,
    .vector_scale_factor = vector_scale_factor,
    .lambda = lambda,
    .num_categories = num_categories,
    .operator = operator,
    .matrix_categorical_vals = matrix_categorical_vals,
    .mean = mean,
    .mean_stderr = mean_stderr,
    .emit_fit_progress = emit_fit_progress,
    .residual_preparation = residual_preparation
  };
  NPRegressionCategoricalProfileFitExecution execution;

  execution.call = &call;
  execution.status = 0;
  np_categorical_profile_owner_init(&execution.owner);
  R_UnwindProtect(
    np_regression_categorical_profile_fit_execute, &execution,
    np_categorical_profile_owner_cleanup, &execution.owner, NULL);
  return execution.status;
}
