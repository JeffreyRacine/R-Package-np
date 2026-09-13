#ifndef NP_REGRESSION_RESIDUAL_H
#define NP_REGRESSION_RESIDUAL_H

/* Shared uncertainty arithmetic, not a replacement fitting owner.
 * A caller supplies the actual influence and certifies constant reproduction
 * before using the off-diagonal identity. No bandwidth or magnitude threshold
 * is an information certificate. Initialize each state with {0}. INVALID is
 * an execution failure; only UNIDENTIFIED denotes unavailable uncertainty.
 */
#include <math.h>
#include <stddef.h>

typedef enum {
  NP_RESIDUAL_IDENTIFIED = 0,
  NP_RESIDUAL_UNIDENTIFIED = 1,
  NP_RESIDUAL_INVALID = 2
} NPResidualInformation;

typedef struct {
  long double scale;
  long double sum;
  long double sumsq;
  long double numerator;
  int invalid;
} NPResidualOffDiagonal;

typedef struct {
  long double scale;
  long double inverse_norm;
  NPResidualInformation information;
} NPResidualNormalization;

/* Geometry is prepared once, then reused across a bounded response tile.
 * Divide a raw weight by its scale before multiplying by inverse_norm, so
 * tiny raw weights never require forming their reciprocal or square. */
static inline NPResidualInformation np_residual_offdiag_normalization(
  const NPResidualOffDiagonal *state, long double denominator_sign,
  NPResidualNormalization *normalization)
{
  if(normalization == NULL) return NP_RESIDUAL_INVALID;
  normalization->scale = 0.0L;
  normalization->inverse_norm = 0.0L;
  normalization->information = NP_RESIDUAL_INVALID;
  if(state == NULL || state->invalid || !isfinite(state->scale) ||
     !isfinite(state->sum) || !isfinite(state->sumsq) ||
     state->scale < 0.0L || state->sumsq < 0.0L ||
     !isfinite(denominator_sign) || denominator_sign == 0.0L)
    return normalization->information;
  if(state->scale == 0.0L) {
    normalization->information = NP_RESIDUAL_UNIDENTIFIED;
    return normalization->information;
  }
  const long double norm = hypotl(state->sum, sqrtl(state->sumsq));
  if(!isfinite(norm) || norm <= 0.0L) return normalization->information;
  normalization->scale = state->scale;
  normalization->inverse_norm = copysignl(1.0L, denominator_sign)/norm;
  normalization->information = NP_RESIDUAL_IDENTIFIED;
  return normalization->information;
}

static inline long double np_residual_normalized_weight(
  const NPResidualNormalization *normalization, long double weight)
{
  if(normalization == NULL ||
     normalization->information != NP_RESIDUAL_IDENTIFIED ||
     !isfinite(weight)) return NAN;
  return (weight/normalization->scale)*normalization->inverse_norm;
}

/* Exact equality of two products of finite doubles. Scaling their mantissas
 * first prevents overflow/underflow; FMA retains the multiplication roundoff
 * rather than mistaking two rounded products for a zero influence. This is
 * an uncertainty-only dependency certificate, not a numerical tolerance. */
static inline int np_residual_products_equal(
  double a, double b, double c, double d)
{
  if(!isfinite(a) || !isfinite(b) || !isfinite(c) || !isfinite(d)) return 0;
  const int left_zero = a == 0.0 || b == 0.0;
  const int right_zero = c == 0.0 || d == 0.0;
  if(left_zero || right_zero) return left_zero && right_zero;
  int ea, eb, ec, ed;
  double ma=frexp(a,&ea), mb=frexp(b,&eb);
  double mc=frexp(c,&ec), md=frexp(d,&ed);
  const int difference=(ea+eb)-(ec+ed);
  if(difference > 1 || difference < -1) return 0;
  if(difference == 1) mc *= 0.5;
  if(difference == -1) ma *= 0.5;
  const double left=ma*mb, right=mc*md;
  return left == right && fma(ma,mb,-left) == fma(mc,md,-right);
}

/* Weights may be signed. LC may supply unnormalized kernel weights; the
 * final denominator sign restores the sign of e/sqrt(q). General smoothers
 * supply their canonical constant-preserving influence directly. */
static inline void np_residual_offdiag_add(
  NPResidualOffDiagonal *state, long double weight, long double difference)
{
  const long double absolute = fabsl(weight);
  if(!isfinite(weight) || !isfinite(difference)) {
    state->invalid = 1;
    return;
  }
  if(absolute == 0.0L) return;
  if(absolute > state->scale) {
    const long double ratio = state->scale / absolute;
    state->sum *= ratio;
    state->sumsq *= ratio * ratio;
    state->numerator *= ratio;
    state->scale = absolute;
  }
  const long double scaled = weight / state->scale;
  state->sum += scaled;
  state->sumsq += scaled * scaled;
  state->numerator += scaled * difference;
}

static inline NPResidualInformation np_residual_offdiag_finish(
  const NPResidualOffDiagonal *state, long double denominator_sign,
  long double *residual)
{
  if(state->invalid || !isfinite(denominator_sign) ||
     denominator_sign == 0.0L || residual == NULL)
    return NP_RESIDUAL_INVALID;
  if(state->scale == 0.0L) {
    *residual = 0.0L; /* finite storage only; information governs use */
    return NP_RESIDUAL_UNIDENTIFIED;
  }
  const long double norm = hypotl(state->sum, sqrtl(state->sumsq));
  *residual = copysignl(1.0L, denominator_sign) * state->numerator / norm;
  return isfinite(*residual) ? NP_RESIDUAL_IDENTIFIED : NP_RESIDUAL_INVALID;
}

/* A compressed categorical profile represents count equal weights. Its
 * centered differences are already summed without expanding donors. */
static inline void np_residual_offdiag_add_group(
  NPResidualOffDiagonal *state, long double weight, long double count,
  long double difference_sum)
{
  const long double absolute = fabsl(weight);
  if(!isfinite(weight) || !isfinite(count) || count < 0.0L ||
     !isfinite(difference_sum) || (count == 0.0L && difference_sum != 0.0L)) {
    state->invalid = 1;
    return;
  }
  if(absolute == 0.0L || count == 0.0L) return;
  if(absolute > state->scale) {
    const long double ratio = state->scale / absolute;
    state->sum *= ratio;
    state->sumsq *= ratio * ratio;
    state->numerator *= ratio;
    state->scale = absolute;
  }
  const long double scaled = weight / state->scale;
  state->sum += count * scaled;
  state->sumsq += count * scaled * scaled;
  state->numerator += scaled * difference_sum;
}

typedef struct {
  long double scale;
  long double sumsq;
  size_t nonzero_influences;
  int unavailable;
  int invalid;
} NPResidualVariance;

static inline void np_residual_variance_add(
  NPResidualVariance *state, long double influence,
  int certified_zero, NPResidualInformation information, long double residual)
{
  if(!isfinite(influence) || information == NP_RESIDUAL_INVALID ||
     (certified_zero && influence != 0.0L)) {
    state->invalid = 1;
    return;
  }
  /* A producer-owned structural certificate, not rounded subtraction,
   * proves independence from an unidentified donor. */
  if(certified_zero) return;
  if(influence != 0.0L) ++state->nonzero_influences;
  if(information == NP_RESIDUAL_UNIDENTIFIED) {
    state->unavailable = 1;
    return;
  }
  if(information != NP_RESIDUAL_IDENTIFIED || !isfinite(residual)) {
    state->invalid = 1;
    return;
  }
  if(influence == 0.0L) return;
  const long double value = fabsl(influence * residual);
  if(!isfinite(value)) {
    state->invalid = 1;
  } else if(value > state->scale) {
    const long double ratio = state->scale / value;
    state->sumsq = 1.0L + state->sumsq * ratio * ratio;
    state->scale = value;
  } else if(value != 0.0L) {
    const long double ratio = value / state->scale;
    state->sumsq += ratio * ratio;
  }
}

static inline NPResidualInformation np_residual_variance_finish(
  const NPResidualVariance *state, long double *standard_error)
{
  if(state->invalid || standard_error == NULL) return NP_RESIDUAL_INVALID;
  if(state->unavailable) return NP_RESIDUAL_UNIDENTIFIED;
  *standard_error = state->scale * sqrtl(state->sumsq);
  return isfinite(*standard_error) ? NP_RESIDUAL_IDENTIFIED : NP_RESIDUAL_INVALID;
}
#endif
