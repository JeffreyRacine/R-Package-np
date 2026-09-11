/* Requested-only empirical LP endpoint covariance. No point/search owner calls
 * this entry. Kernel rows and the realized basis are supplied by the existing
 * conditional adapters; solves use the canonical accepted ridge policy. */
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include "jksum_lp_solve.h"

typedef struct {
  SEXP basis, evaluation, weights, response, slope;
  int n, p, m;
  NPLPSolveWorkspace solve;
} NPLPPairCall;

static int np_lp_pair_dims(SEXP x, int nr, int nc)
{
  SEXP d = getAttrib(x, R_DimSymbol);
  return TYPEOF(x) == REALSXP && TYPEOF(d) == INTSXP &&
    XLENGTH(d) == 2 && INTEGER(d)[0] == nr && INTEGER(d)[1] == nc;
}

static void np_lp_pair_cleanup(void *data, Rboolean jump)
{
  (void)jump;
  np_lp_solve_workspace_clear(&((NPLPPairCall *)data)->solve);
}

/* 0=available, 1=empty support, 2=nondifferentiable ridge scale,
 * 3=nonfinite influence, 4=nonpositive quantile density. */
static int np_lp_pair_endpoint(NPLPPairCall *call, int query,
  double *weighted, double *theta, double *direction, double *score,
  double *point, double *slope)
{
  const int n = call->n, p = call->p, rows = 2*call->m, one = 1;
  const char trans_t = 'T', trans_n = 'N';
  const double alpha = 1.0, zero = 0.0;
  const double *z = REAL(call->basis);
  const double *w = REAL(call->weights) + (size_t)query*n;
  const double *r = REAL(call->response) + (size_t)query*n;
  NPLPSolveWorkspace *s = &call->solve;
  NPLPSolvePolicyDiagnostics diagnostics = {0, 0.0};
  double anchor, b0, max_abs = 0.0, ridge_sign = 1.0;
  int active = 0, scale_a = 0, scale_b = 0;
  *point = *slope = NA_REAL;
  for(int i = 0; i < n; ++i) active += w[i] != 0.0;
  if(!active) return 1;
  for(int j = 0; j < p; ++j)
    for(int i = 0; i < n; ++i)
      weighted[i + (size_t)n*j] = z[i + (size_t)n*j]*w[i];
  F77_CALL(dgemm)(&trans_t, &trans_n, &p, &p, &n, &alpha, z, &n,
    weighted, &n, &zero, s->gram_source, &p FCONE FCONE);
  F77_CALL(dgemv)(&trans_t, &n, &p, &alpha, weighted, &n, r, &one,
    &zero, s->rhs_source, &one FCONE);
  anchor = s->gram_source[0];
  b0 = s->rhs_source[0];
  for(int j = 0; j < p; ++j) {
    const double a = fabs(s->gram_source[j + (size_t)p*j]);
    if(a > max_abs) { max_abs = a; scale_a = scale_b = j; }
  }
  /* Match the canonical scale's off-diagonal contingency for signed rows. */
  const int diagonal_scale = max_abs > 0.0;
  if(!diagonal_scale)
    for(int j = 0; j < p; ++j)
      for(int k = 0; k < p; ++k) {
        const double a = fabs(s->gram_source[k + (size_t)p*j]);
        if(a > max_abs) { max_abs = a; scale_a = k; scale_b = j; }
      }
  ridge_sign = s->gram_source[scale_a + (size_t)p*scale_b] < 0.0 ? -1.0 : 1.0;
  /* Freeze the pristine maximum before the solve mutates the diagonal.
   * A tie is smooth only when its empirical-mass derivatives agree. */
  int scale_nonsmooth = 0;
  for(int j = 0; j < p; ++j)
    for(int k = diagonal_scale ? j : 0; k < (diagonal_scale ? j+1 : p); ++k)
      if((k != scale_a || j != scale_b) &&
         fabs(s->gram_source[k + (size_t)p*j]) == max_abs) {
        const double sign = s->gram_source[k + (size_t)p*j] < 0.0 ? -1.0 : 1.0;
        for(int i = 0; i < n && !scale_nonsmooth; ++i)
          if(w[i] != 0.0 &&
             sign*z[i + (size_t)n*k]*z[i + (size_t)n*j] !=
             ridge_sign*z[i + (size_t)n*scale_a]*z[i + (size_t)n*scale_b])
            scale_nonsmooth = 1;
      }

  if(np_lp_solve_workspace_solve_response_ranked(s, p, 1, 1.0/n,
       NP_LP_RANK_UPPER_BOUND_UNKNOWN, &diagnostics) != NP_LP_SOLVE_POLICY_OK)
    error("conditional LP contrast covariance: accepted endpoint solve failed");
  memcpy(theta, s->rhs_work, (size_t)p*sizeof(double));
  *point = 0.0;
  for(int j = 0; j < p; ++j) {
    direction[j] = REAL(call->evaluation)[query + (size_t)rows*j];
    *point += direction[j]*theta[j];
    s->rhs_source[j] = direction[j];
  }
  /* The derivative needs S^{-T} d. Use the retained factor, with no second
   * admission/ridge selection or assumption of bitwise Gram symmetry. */
  if(!np_lp_solve_workspace_solve_transpose_factored(s, p, 1))
    error("conditional LP contrast covariance: retained direction solve failed");
  memcpy(direction, s->rhs_work, (size_t)p*sizeof(double));
  const double ridge = diagnostics.ridge_total;
  const double c0 = ridge == 0.0 ? 1.0 : 1.0 + ridge/anchor;
  if(!R_FINITE(*point) || !R_FINITE(c0)) return 3;

  if(ridge > 0.0 && scale_nonsmooth) return 2;
  double vt = 0.0;
  for(int j = 0; j < p; ++j) vt += direction[j]*theta[j];
  double density = 0.0;
  for(int i = 0; i < n; ++i) {
    double vz = 0.0, zt = 0.0;
    for(int j = 0; j < p; ++j) {
      vz += direction[j]*z[i + (size_t)n*j];
      zt += theta[j]*z[i + (size_t)n*j];
    }
    double value = w[i]*vz*(r[i]-zt);
    if(ridge > 0.0) {
      const double da = w[i]*z[i]*z[i];
      const double dridge = ((double)diagnostics.ridge_steps/n)*ridge_sign*
        w[i]*z[i + (size_t)n*scale_a]*z[i + (size_t)n*scale_b];
      value += w[i]*direction[0]*(ridge/anchor)*z[i]*r[i] +
        direction[0]*b0*(dridge/anchor - (ridge/anchor)*(da/anchor)) -
        dridge*vt;
    }
    if(!R_FINITE(value)) return 3;
    score[i] = value;
    if(call->slope != R_NilValue)
      density += w[i]*(vz + (c0-1.0)*direction[0]*z[i])*
        REAL(call->slope)[i + (size_t)n*query];
  }
  if(call->slope != R_NilValue) {
    *slope = density;
    if(!R_FINITE(density) || density <= 0.0) return 4;
    for(int i = 0; i < n; ++i) {
      score[i] = -score[i]/density;
      if(!R_FINITE(score[i])) return 3;
    }
  }
  return 0;
}

static SEXP np_lp_pair_run(void *data)
{
  NPLPPairCall *call = (NPLPPairCall *)data;
  const int n = call->n, p = call->p, m = call->m;
  double *weighted = (double *)R_alloc((size_t)n*p, sizeof(double));
  double *theta = (double *)R_alloc(p, sizeof(double));
  double *direction = (double *)R_alloc(p, sizeof(double));
  double *upper = (double *)R_alloc(n, sizeof(double));
  double *lower = (double *)R_alloc(n, sizeof(double));
  SEXP out = PROTECT(allocMatrix(REALSXP, m, 6));
  for(int j = 0; j < m; ++j) {
    double pu, pl, fu, fl;
    R_CheckUserInterrupt();
    int status = np_lp_pair_endpoint(call, j, weighted, theta, direction,
                                     upper, &pu, &fu);
    int status_l = np_lp_pair_endpoint(call, j+m, weighted, theta, direction,
                                       lower, &pl, &fl);
    if(!status) status = status_l;
    double se = NA_REAL;
    if(!status) {
      double mean = 0.0, scale = 0.0, sumsq = 1.0;
      for(int i = 0; i < n; ++i) {
        upper[i] -= lower[i];
        mean += (upper[i]-mean)/(i+1.0);
      }
      for(int i = 0; i < n; ++i) {
        const double value = fabs(upper[i]-mean);
        if(value > scale) {
          const double ratio = scale/value;
          sumsq = 1.0 + sumsq*ratio*ratio;
          scale = value;
        } else if(value != 0.0) {
          const double ratio = value/scale;
          sumsq += ratio*ratio;
        }
      }
      se = scale*sqrt(sumsq*((double)n/(n-1.0)));
      if(!R_FINITE(se)) { se = NA_REAL; status = 3; }
    }
    REAL(out)[j] = se;
    REAL(out)[j+m] = pu;
    REAL(out)[j+(size_t)2*m] = pl;
    REAL(out)[j+(size_t)3*m] = fu;
    REAL(out)[j+(size_t)4*m] = fl;
    REAL(out)[j+(size_t)5*m] = status;
  }
  UNPROTECT(1);
  return out;
}

SEXP C_np_conditional_lp_pair_se(SEXP basis, SEXP evaluation, SEXP weights,
                                SEXP response, SEXP slope)
{
  SEXP d = getAttrib(basis, R_DimSymbol);
  SEXP e = getAttrib(evaluation, R_DimSymbol);
  if(TYPEOF(basis) != REALSXP || TYPEOF(evaluation) != REALSXP ||
     TYPEOF(d) != INTSXP || XLENGTH(d) != 2 ||
     TYPEOF(e) != INTSXP || XLENGTH(e) != 2)
    error("conditional LP contrast covariance requires matrix inputs");
  const int n = INTEGER(d)[0], p = INTEGER(d)[1], ne = INTEGER(e)[0];
  if(n < 2 || p < 1 || ne < 2 || ne % 2 != 0 ||
     INTEGER(e)[1] != p || !np_lp_pair_dims(weights, n, ne) ||
     !np_lp_pair_dims(response, n, ne) ||
     (slope != R_NilValue && !np_lp_pair_dims(slope, n, ne)) ||
     (size_t)n > SIZE_MAX/sizeof(double)/(size_t)p)
    error("conditional LP contrast covariance has inconsistent dimensions");
  SEXP inputs[5] = {basis, evaluation, weights, response, slope};
  for(int j = 0; j < (slope == R_NilValue ? 4 : 5); ++j)
    for(R_xlen_t i = 0; i < XLENGTH(inputs[j]); ++i)
      if(!R_FINITE(REAL(inputs[j])[i]))
        error("conditional LP contrast covariance received non-finite input");
  NPLPPairCall call = {basis, evaluation, weights, response, slope, n, p, ne/2, {0}};
  np_lp_solve_workspace_init(&call.solve);
  if(!np_lp_solve_workspace_reserve(&call.solve, p, 1))
    error("conditional LP contrast covariance workspace allocation failed");
  return R_UnwindProtect(np_lp_pair_run, &call, np_lp_pair_cleanup, &call, NULL);
}
