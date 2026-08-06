#ifndef NP_REGRESSION_CONTRACT_H
#define NP_REGRESSION_CONTRACT_H

/*
 * Standard-error algebra belongs to the regression response consumer, not to
 * a continuous-kernel family.  Ordinary regression retains its local-
 * residual covariance contract.  A scalar conditional density/distribution
 * fit instead treats the dependent-side kernel row as the response and uses
 * the established observation-influence variance.
 */
typedef enum {
  NP_REGRESSION_STDERR_LOCAL_RESIDUAL = 0,
  NP_REGRESSION_STDERR_CONDITIONAL_INFLUENCE = 1
} NPRegressionStandardErrorMode;

/*
 * Public regression fits select one explicit output contract at the native
 * boundary.  The values are intentionally bit flags so all owners can derive
 * the same request without route-specific boolean combinations or silent
 * promotion to a more expensive mode.
 */
typedef enum {
  NP_REGRESSION_OUTPUT_MEAN_ONLY = 0,
  NP_REGRESSION_OUTPUT_MEAN_AND_SE = 1,
  NP_REGRESSION_OUTPUT_MEAN_AND_GRADIENT = 2,
  NP_REGRESSION_OUTPUT_FULL = 3
} NPRegressionOutputRequest;

static inline int np_regression_output_request_valid(const int request)
{
  return request >= NP_REGRESSION_OUTPUT_MEAN_ONLY &&
    request <= NP_REGRESSION_OUTPUT_FULL;
}

static inline int np_regression_output_requests_errors(const int request)
{
  return (request & NP_REGRESSION_OUTPUT_MEAN_AND_SE) != 0;
}

static inline int np_regression_output_requests_gradients(const int request)
{
  return (request & NP_REGRESSION_OUTPUT_MEAN_AND_GRADIENT) != 0;
}

/*
 * Read-only consumer view of a call-scoped realized continuous-bandwidth
 * layout.  The owning conditional loop advances only evaluation_offset;
 * bandwidth storage remains immutable after realization.
 * Conditional response-row consumers can reuse one full X realization while
 * invoking the common regression owner on successive response rows.  The
 * owner copies only its current view into its existing bounded workspace; it
 * never assumes ownership of these pointers.
 */
typedef struct {
  int bandwidth_mode;
  int num_train;
  int num_eval_total;
  int num_continuous;
  int evaluation_offset;
  int evaluation_count;
  double **bandwidth_eval;
  double **bandwidth_train;
} NPContinuousPreparedBandwidthView;

#endif
