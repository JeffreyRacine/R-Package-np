#ifndef NP_CONDITIONAL_LOCAL_QR_H
#define NP_CONDITIONAL_LOCAL_QR_H

/* Conditional-only local-design rank and QR owner. The caller supplies original X,
 * complete-product weights and the existing downward-closed term table.
 * Public basis preparation and all regression dependencies remain unchanged.
 *
 * fast owns QR storage and may reserve only when dimensions change. policy
 * lends gram_work, rank_values and rank_work; no pristine Gram/RHS is changed.
 * No row is published unless classification and reconstruction both succeed.
 * Deficiency/ambiguity are policy results, never finite objective values. */
typedef enum {
  NP_CQR_LOCAL_FULL = 0,
  NP_CQR_LOCAL_DEFICIENT,
  NP_CQR_LOCAL_AMBIGUOUS,
  NP_CQR_LOCAL_EMPTY,
  NP_CQR_LOCAL_NONFINITE,
  NP_CQR_LOCAL_FAILED
} NPConditionalLocalQRStatus;

static double np_cqr_local_power(double x, int exponent)
{
  double value = 1.0;
  while(exponent > 0){
    if(exponent & 1) value *= x;
    exponent >>= 1;
    if(exponent) x *= x;
  }
  return value;
}

static NPConditionalLocalQRStatus np_cqr_local_row(
  NPConditionalQRFast *fast, NPLPSolveWorkspace *policy,
  int n, int dimensions, int terms_count, const int *terms,
  double * const *x, const double *weights, int position,
  double *row, double *ratio_out)
{
  const int one = 1;
  const int k = terms_count;
  int info = 0;
  int positive = 0;
  double maxw = 0.0;

  if(ratio_out) *ratio_out = NA_REAL;
  if(!fast || !policy || !terms || !x || !weights || !row ||
     n < 2 || dimensions < 1 || k < 1 || position < 0 || position >= n ||
     k > INT_MAX/5 || policy->p_capacity < k ||
     !policy->gram_work || !policy->rank_values || !policy->rank_work ||
     policy->rank_work_capacity < 5U*(size_t)k)
    return NP_CQR_LOCAL_FAILED;
  for(int t = 0; t < k; ++t){
    for(int d = 0; d < dimensions; ++d){
      const int power = terms[(size_t)t*dimensions+d];
      if(power < 0 || (t == 0 && power != 0)) return NP_CQR_LOCAL_FAILED;
    }
  }
  for(int i = 0; i < n; ++i){
    if(!R_FINITE(weights[i])) return NP_CQR_LOCAL_NONFINITE;
    if(weights[i] < 0.0) return NP_CQR_LOCAL_FAILED;
    if(i != position && weights[i] > 0.0){
      ++positive;
      maxw = fmax(maxw, weights[i]);
    }
  }
  if(!positive) return NP_CQR_LOCAL_EMPTY;
  if(positive < k) return NP_CQR_LOCAL_DEFICIENT;
  if(np_cqr_fast_init(fast, n, k)) return NP_CQR_LOCAL_FAILED;

  for(int i = 0; i < n; ++i)
    fast->sqrtw[i] = i == position ? 0.0 : sqrt(weights[i]/maxw);
  for(int t = 0; t < k; ++t){
    double *column = fast->a + (size_t)fast->capacity*t;
    memcpy(column, fast->sqrtw, (size_t)n*sizeof(double));
  }
  /* A scale multiplies each monomial column by a single positive constant,
   * which cancels under column normalization. No rounded raw basis is shifted. */
  for(int d = 0; d < dimensions; ++d){
    int active = 0;
    for(int t = 0; t < k; ++t)
      if(terms[(size_t)t*dimensions+d]) active = 1;
    if(!active) continue;
    if(!x[d] || !R_FINITE(x[d][position])) return NP_CQR_LOCAL_NONFINITE;
    const double center = x[d][position];
    double radius = 0.0;
    for(int i = 0; i < n; ++i){
      if(i == position || weights[i] == 0.0) continue;
      const double delta = x[d][i] - center;
      if(!R_FINITE(delta)) return NP_CQR_LOCAL_NONFINITE;
      radius = fmax(radius, fabs(delta));
    }
    if(!(radius > 0.0)) return NP_CQR_LOCAL_DEFICIENT;
    for(int t = 0; t < k; ++t){
      const int power = terms[(size_t)t*dimensions+d];
      if(!power) continue;
      double *column = fast->a + (size_t)fast->capacity*t;
      for(int i = 0; i < n; ++i){
        if(fast->sqrtw[i] == 0.0) continue;
        column[i] *= np_cqr_local_power((x[d][i]-center)/radius, power);
      }
    }
  }
  for(int t = 0; t < k; ++t){
    double *column = fast->a + (size_t)fast->capacity*t;
    fast->scale[t] = F77_CALL(dnrm2)(&n, column, &one);
    if(!R_FINITE(fast->scale[t])) return NP_CQR_LOCAL_NONFINITE;
    if(!(fast->scale[t] > 0.0)) return NP_CQR_LOCAL_DEFICIENT;
    for(int i = 0; i < n; ++i) column[i] /= fast->scale[t];
    fast->pivot[t] = 0;
  }
  F77_CALL(dgeqp3)(&n, &k, fast->a, &fast->capacity, fast->pivot,
                   fast->tau, fast->work, &fast->lwork, &info);
  if(info) return NP_CQR_LOCAL_FAILED;
  for(int j = 0; j < k; ++j)
    for(int i = 0; i < k; ++i)
      policy->gram_work[i+(size_t)k*j] = i <= j ?
        fast->a[i+(size_t)fast->capacity*j] : 0.0;
  const char no = 'N';
  const int lwork = 5*k;
  double dummy = 0.0;
  F77_CALL(dgesvd)(&no, &no, &k, &k, policy->gram_work, &k,
                   policy->rank_values, &dummy, &one, &dummy, &one,
                   policy->rank_work, &lwork, &info FCONE FCONE);
  if(info) return NP_CQR_LOCAL_FAILED;
  const double ratio = policy->rank_values[k-1]/policy->rank_values[0];
  const double neps = ((double)n+k)*DBL_EPSILON;
  const double threshold = neps/(1.0-neps);
  if(!R_FINITE(ratio)) return NP_CQR_LOCAL_NONFINITE;
  if(ratio_out) *ratio_out = ratio;
  if(ratio <= threshold/10.0) return NP_CQR_LOCAL_DEFICIENT;
  if(ratio <= 10.0*threshold) return NP_CQR_LOCAL_AMBIGUOUS;

  memset(fast->v, 0, (size_t)n*sizeof(double));
  for(int i = 0; i < k; ++i){
    const int pivot = fast->pivot[i]-1;
    if(pivot < 0 || pivot >= k) return NP_CQR_LOCAL_FAILED;
    double value = pivot == 0 ? 1.0/fast->scale[0] : 0.0;
    for(int j = 0; j < i; ++j)
      value -= fast->a[j+(size_t)fast->capacity*i]*fast->v[j];
    fast->v[i] = value/fast->a[i+(size_t)fast->capacity*i];
    if(!R_FINITE(fast->v[i])) return NP_CQR_LOCAL_NONFINITE;
  }
  F77_CALL(dormqr)("L", "N", &n, &one, &k, fast->a, &fast->capacity,
                   fast->tau, fast->v, &fast->capacity, fast->work,
                   &fast->lwork, &info FCONE FCONE);
  if(info) return NP_CQR_LOCAL_FAILED;
  for(int i = 0; i < n; ++i)
    if(!R_FINITE(fast->sqrtw[i]*fast->v[i])) return NP_CQR_LOCAL_NONFINITE;
  for(int i = 0; i < n; ++i) row[i] = fast->sqrtw[i]*fast->v[i];
  return NP_CQR_LOCAL_FULL;
}
#endif
