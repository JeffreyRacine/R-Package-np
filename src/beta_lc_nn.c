/* Univariate second-order beta density for frozen-radius LC derivative rows.
 * The caller owns the physical NN geometry, row partition, sums and jumps. */
#include <float.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include "headers.h"
#include "continuous_kernel_row.h"

typedef struct {
  double near_unit, log_unit, log_complement, concentration;
  int upper;
} NPBetaLCNNDonor;

typedef struct NPBetaLCNNRows {
  NPBetaLCNNDonor *donors;
  int evaluation, shared_radius;
  double lower, upper, span, log_span;
  double query, u, v, lower_gap, upper_gap, log_u, log_v;
  double a1, b1, mode_log;
  np_beta_pdf_derivative_component score;
} NPBetaLCNNRows;

static np_beta_status prepare_donors(const NPContinuousKernelRowPlan *p,
                                    NPBetaLCNNRows **out)
{
  const NPContinuousKernelSegment *segment = &p->route->segment[0];
  /* The derivative-only kernel-sum plan may have no PDF preparation.
   * Give its rank-local scratch storage the same unwind lifetime. */
  if(!p->beta_prepared->allocation_active) {
    p->beta_prepared->allocation_marker = vmaxget();
    p->beta_prepared->allocation_active = 1;
  }
  NPBetaLCNNRows *r = (NPBetaLCNNRows *)R_alloc(1, sizeof(*r));
  r->donors = (NPBetaLCNNDonor *)R_alloc(p->num_train, sizeof(*r->donors));
  r->lower = segment->lower[0];
  r->upper = segment->upper[0];
  r->span = r->upper - r->lower;
  r->log_span = log(r->span);
  r->shared_radius = p->bandwidth_mode == BW_GEN_NN;
  r->evaluation = -1;
  r->score.concentration = 0.0;
  for(int i = 0; i < p->num_train; ++i) {
    const double x = p->train[0][i];
    np_beta_shape s;
    const double h = r->shared_radius ? r->span : p->bandwidth_train[0][i];
    const np_beta_status status = np_beta_shape_init(r->lower, x, h,
      r->lower, r->upper, 1, &s);
    if(status != NP_BETA_OK) return status;
    NPBetaLCNNDonor *d = &r->donors[i];
    d->upper = s.observation_unit > 0.5;
    d->near_unit = d->upper ? s.observation_complement_unit : s.observation_unit;
    d->log_unit = s.log_observation_unit;
    d->log_complement = s.log_observation_complement_unit;
    d->concentration = s.concentration;
  }
  *out = r;
  return NP_BETA_OK;
}

static np_beta_status prepare_query(const NPContinuousKernelRowPlan *p,
                                   NPBetaLCNNRows *r, int j)
{
  r->query = p->train_is_eval ? p->train[0][j] : p->evaluation[0][j];
  r->lower_gap = r->query - r->lower;
  r->upper_gap = r->upper - r->query;
  r->u = r->lower_gap / r->span;
  r->v = r->upper_gap / r->span;
  if(r->shared_radius && r->query >= r->lower && r->query <= r->upper) {
    np_beta_shape s;
    const np_beta_status status = np_beta_shape_init(r->query, r->query,
      p->bandwidth_eval[0][j], r->lower, r->upper, 1, &s);
    if(status != NP_BETA_OK) return status;
    r->log_u = s.log_observation_unit;
    r->log_v = s.log_observation_complement_unit;
    r->a1 = s.target_unit * s.concentration;
    r->b1 = s.target_complement_unit * s.concentration;
    r->score.concentration = s.concentration;
    r->score.support_length = r->span;
    r->score.digamma_alpha = digamma(1.0 + r->a1);
    r->score.digamma_beta = digamma(1.0 + r->b1);
    r->mode_log = (s.observation_unit <= 0.5 ?
      dbeta(s.observation_unit, 1.0+r->a1, 1.0+r->b1, 1) :
      dbeta(s.observation_complement_unit, 1.0+r->b1, 1.0+r->a1, 1)) - r->log_span;
  }
  r->evaluation = j;
  return NP_BETA_OK;
}

/* The GNN normalizer is common to the row. Measure the log density from
 * its mode using physical gaps. Use the complementary log when a direct
 * distance ratio loses the smaller endpoint. */
static double shared_log_density(const NPBetaLCNNRows *r, const NPBetaLCNNDonor *d,
                                 double x)
{
  if(r->a1 == 0.0)
    return r->mode_log + r->b1 * (d->upper ? d->log_complement : log1p(-d->near_unit));
  if(r->b1 == 0.0)
    return r->mode_log + r->a1 * (d->upper ? log1p(-d->near_unit) : d->log_unit);
  const double dx = x - r->query;
  const double du = dx / r->lower_gap, dv = -dx / r->upper_gap;
  if(fmax(fabs(du),fabs(dv)) < sqrt(DBL_EPSILON))
    return r->mode_log + r->a1 * Rf_log1pmx(du) + r->b1 * Rf_log1pmx(dv);
  const double lu = fabs(du) < 0.5 ? log1p(du) : d->log_unit - r->log_u;
  const double lv = fabs(dv) < 0.5 ? log1p(dv) : d->log_complement - r->log_v;
  return r->mode_log + fma(r->a1, lu, r->b1 * lv);
}

np_beta_status np_beta_lc_nn_terms(const NPContinuousKernelRowPlan *plan,
  int evaluation_index, int observation, double *level, np_beta_derivative *derivative)
{
  NPBetaLCNNRows *r = plan->beta_prepared->lc_nn_rows;
  if(r == NULL) {
    np_beta_status status = prepare_donors(plan, &r);
    if(status != NP_BETA_OK) return status;
    plan->beta_prepared->lc_nn_rows = r;
  }
  if(evaluation_index != r->evaluation) {
    np_beta_status status = prepare_query(plan, r, evaluation_index);
    if(status != NP_BETA_OK) return status;
  }
  const double x = plan->train[0][observation];
  const NPBetaLCNNDonor *d = &r->donors[observation];
  const double tau = r->shared_radius ? r->score.concentration : d->concentration;
  if(x == r->lower || x == r->upper || r->query < r->lower ||
     r->query > r->upper || tau == 0.0) {
    const double h = r->shared_radius ? plan->bandwidth_eval[0][evaluation_index] :
                                       plan->bandwidth_train[0][observation];
    np_beta_status status;
    int sign;
    if(derivative != NULL)
      return np_beta_log_abs_pdf_derivative_order(r->query, x, h, r->lower, r->upper,
        2, level, &sign, derivative);
    *level = np_beta_log_abs_pdf_order(r->query, x, h, r->lower, r->upper, 2, &sign, &status);
    return status;
  }
  if(r->shared_radius) {
    *level = shared_log_density(r, d, x);
  } else {
    const double a = 1.0 + r->u * tau, b = 1.0 + r->v * tau;
    *level = (d->upper ? dbeta(d->near_unit,b,a,1) : dbeta(d->near_unit,a,b,1)) -
      r->log_span;
  }
  if(ISNAN(*level) || *level == INFINITY) return NP_BETA_ERR_NUMERIC;
  if(derivative != NULL) {
    np_beta_shape s = {.lower=r->lower, .upper=r->upper, .support_length=r->span,
      .target_unit=r->u, .target_complement_unit=r->v, .concentration=tau,
      .log_observation_unit=d->log_unit, .log_observation_complement_unit=d->log_complement};
    return np_beta_single_component_derivative(&s, r->shared_radius ? &r->score : NULL,
      r->query, x, *level, derivative);
  }
  return NP_BETA_OK;
}
