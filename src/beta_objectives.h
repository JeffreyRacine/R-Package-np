#ifndef NP_BETA_OBJECTIVES_H
#define NP_BETA_OBJECTIVES_H

/* C-only cross-validation objectives for the coordinate-aware beta family.
 * Candidate bandwidth vectors are zero-based and contain one entry per
 * continuous variable.  Fixed-bandwidth scale-factor conversion and the
 * package's generalized/adaptive nearest-neighbor rules are inherited from
 * kernel_bandwidth_mean().  Each routine returns DBL_MAX when the candidate
 * cannot be evaluated. */

double np_beta_objective_density_ml(
  int bandwidth_mode,
  int order,
  int num_obs,
  int num_continuous,
  double **train_continuous,
  const double *candidate,
  const double *lower,
  const double *upper);

double np_beta_objective_density_ls(
  int bandwidth_mode,
  int order,
  int num_obs,
  int num_continuous,
  double **train_continuous,
  const double *candidate,
  const double *lower,
  const double *upper);

double np_beta_objective_distribution_ls(
  int bandwidth_mode,
  int order,
  int num_obs_train,
  int num_obs_eval,
  int num_continuous,
  int cdf_on_train,
  double **train_continuous,
  double **eval_continuous,
  const double *candidate,
  const double *lower,
  const double *upper);

double np_beta_objective_regression_lc(
  int bandwidth_mode,
  int order,
  int use_aic,
  int num_obs,
  int num_continuous,
  double **train_continuous,
  const double *response,
  const double *candidate,
  const double *lower,
  const double *upper);

#endif
