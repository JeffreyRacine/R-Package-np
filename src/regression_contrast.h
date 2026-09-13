#ifndef NP_REGRESSION_CONTRAST_H
#define NP_REGRESSION_CONTRAST_H

#include <math.h>
#include "regression_residual.h"

/* Two-double arithmetic is intentional: long double is not wider than
 * double on every supported platform. No magnitude threshold is used. */
typedef struct { double hi, lo; } NPContrastNumber;

static inline NPContrastNumber np_contrast_number(double value)
{
  const NPContrastNumber result = {value,0.0};
  return result;
}

static inline int np_contrast_finite(NPContrastNumber value)
{
  return isfinite(value.hi) && isfinite(value.lo);
}

static inline NPContrastNumber np_contrast_two_sum(double a, double b)
{
  const double sum=a+b, bv=sum-a;
  const NPContrastNumber result={sum,(a-(sum-bv))+(b-bv)};
  return result;
}

static inline NPContrastNumber np_contrast_add(
  NPContrastNumber a, NPContrastNumber b)
{
  NPContrastNumber high=np_contrast_two_sum(a.hi,b.hi);
  const NPContrastNumber low=np_contrast_two_sum(a.lo,b.lo);
  high=np_contrast_two_sum(high.hi,high.lo+low.hi);
  return np_contrast_two_sum(high.hi,high.lo+low.lo);
}

static inline NPContrastNumber np_contrast_negative(NPContrastNumber value)
{
  const NPContrastNumber result={-value.hi,-value.lo};
  return result;
}

static inline NPContrastNumber np_contrast_subtract(
  NPContrastNumber a, NPContrastNumber b)
{
  return np_contrast_add(a,np_contrast_negative(b));
}

static inline NPContrastNumber np_contrast_multiply(
  NPContrastNumber a, NPContrastNumber b)
{
  const double product=a.hi*b.hi;
  const double cross=a.hi*b.lo+a.lo*b.hi;
  const double remainder=(fma(a.hi,b.hi,-product)+cross)+a.lo*b.lo;
  return np_contrast_two_sum(product,remainder);
}

static inline NPContrastNumber np_contrast_two_product(double a, double b)
{
  const double product=a*b;
  const NPContrastNumber result={product,fma(a,b,-product)};
  return result;
}

static inline NPContrastNumber np_contrast_scale(NPContrastNumber a, double factor)
{
  const double product=a.hi*factor;
  return np_contrast_two_sum(product,fma(a.hi,factor,-product)+a.lo*factor);
}

static inline NPContrastNumber np_contrast_divide(
  NPContrastNumber numerator, NPContrastNumber denominator)
{
  if(!np_contrast_finite(numerator) || !np_contrast_finite(denominator) ||
     denominator.hi == 0.0) return np_contrast_number(NAN);
  const NPContrastNumber first=np_contrast_number(numerator.hi/denominator.hi);
  NPContrastNumber remainder=np_contrast_subtract(numerator,
    np_contrast_multiply(denominator,first));
  NPContrastNumber result=np_contrast_add(first,
    np_contrast_number(remainder.hi/denominator.hi));
  remainder=np_contrast_subtract(numerator,np_contrast_multiply(denominator,result));
  return np_contrast_add(result,np_contrast_number(remainder.hi/denominator.hi));
}

static inline NPContrastNumber np_contrast_ratio_difference(
  double left, NPContrastNumber left_sum,
  double right, NPContrastNumber right_sum)
{
  return np_contrast_subtract(
    np_contrast_divide(np_contrast_number(left),left_sum),
    np_contrast_divide(np_contrast_number(right),right_sum));
}

typedef struct {
  NPContrastNumber left_sum,right_sum,denominator;
  double reciprocal;
  int direct;
} NPContrastRatio;

static inline NPContrastRatio np_contrast_ratio_prepare(
  NPContrastNumber left_sum,NPContrastNumber right_sum)
{
  NPContrastRatio ratio={0};
  ratio.left_sum=left_sum;ratio.right_sum=right_sum;
  ratio.denominator=np_contrast_multiply(left_sum,right_sum);
  ratio.reciprocal=1.0/ratio.denominator.hi;
  ratio.direct=np_contrast_finite(ratio.denominator) &&
    isfinite(ratio.reciprocal) && ratio.reciprocal!=0.0;
  return ratio;
}

/* Normalizers are fixed over a resident row. Cache their complete product
 * once rather than performing repeated corrected divisions for every donor.
 * The uncommon exponent-overflow/underflow cases keep the precise ratio
 * division path; this is arithmetic selection, not a kernel/owner fallback. */
static inline NPContrastNumber np_contrast_ratio_apply(
  const NPContrastRatio *ratio,double left,double right)
{
  const double lp=left*ratio->right_sum.hi,rp=right*ratio->left_sum.hi;
  if(!ratio->direct || !isfinite(lp) || !isfinite(rp) ||
     (left!=0.0 && lp==0.0) || (right!=0.0 && rp==0.0))
    return np_contrast_ratio_difference(left,ratio->left_sum,right,ratio->right_sum);
  NPContrastNumber numerator=np_contrast_two_sum(lp,-rp);
  numerator.lo+=(fma(left,ratio->right_sum.hi,-lp)-fma(right,ratio->left_sum.hi,-rp))+
    (left*ratio->right_sum.lo-right*ratio->left_sum.lo);
  const double high=numerator.hi*ratio->reciprocal;
  const double remainder=(fma(-high,ratio->denominator.hi,numerator.hi)+numerator.lo)-
    high*ratio->denominator.lo;
  return np_contrast_two_sum(high,remainder*ratio->reciprocal);
}

static inline double np_contrast_value(NPContrastNumber value)
{
  return value.hi+value.lo;
}

/* Optional exact normalizer provenance for unidentified-donor dependency
 * certificates. An expansion is an exact sum of its input doubles, unlike
 * an ordinary or two-double rounded normalizer. This cold helper may decline
 * a certificate if its finite capacity or representability is exhausted;
 * it must never turn that into a false zero influence. */
enum { NP_CONTRAST_EXACT_CAPACITY=64 };
typedef struct {
  double component[NP_CONTRAST_EXACT_CAPACITY];
  int size;
  int invalid;
} NPContrastExactSum;

static inline void np_contrast_exact_add(NPContrastExactSum *sum,double value)
{
  if(sum->invalid) return;
  if(!isfinite(value)) { sum->invalid=1; return; }
  double result[NP_CONTRAST_EXACT_CAPACITY];
  int count=0;
  for(int i=0;i<sum->size;++i) {
    const NPContrastNumber part=np_contrast_two_sum(value,sum->component[i]);
    if(!np_contrast_finite(part)) { sum->invalid=1; return; }
    if(part.lo != 0.0) {
      if(count==NP_CONTRAST_EXACT_CAPACITY) {sum->invalid=1;return;}
      result[count++]=part.lo;
    }
    value=part.hi;
  }
  if(value != 0.0) {
    if(count==NP_CONTRAST_EXACT_CAPACITY) {sum->invalid=1;return;}
    result[count++]=value;
  }
  for(int i=0;i<count;++i) sum->component[i]=result[i];
  sum->size=count;
}

static inline int np_contrast_exact_products_equal(
  double left,const NPContrastExactSum *left_sum,
  double right,const NPContrastExactSum *right_sum)
{
  if(!isfinite(left) || !isfinite(right) || left_sum->invalid || right_sum->invalid)
    return 0;
  int minimum=0,maximum=0,any=0;
  for(int side=0;side<2;++side) {
    const double factor=side==0 ? left : right;
    const NPContrastExactSum *sum=side==0 ? left_sum : right_sum;
    if(factor==0.0) continue;
    int factor_exponent;
    const double factor_mantissa=frexp(factor,&factor_exponent);
    for(int i=0;i<sum->size;++i) {
      int exponent;
      const double mantissa=frexp(sum->component[i],&exponent);
      const double high=factor_mantissa*mantissa;
      const double low=fma(factor_mantissa,mantissa,-high);
      for(int part=0;part<2;++part) {
        const double value=part==0 ? high : low;
        if(value==0.0) continue;
        int part_exponent;
        (void)frexp(value,&part_exponent);
        const int total=factor_exponent+exponent+part_exponent;
        if(!any || total<minimum) minimum=total;
        if(!any || total>maximum) maximum=total;
        any=1;
      }
    }
  }
  if(!any) return 1;
  const int origin=minimum+(maximum-minimum)/2;
  NPContrastExactSum difference={0};
  for(int side=0;side<2;++side) {
    const double factor=side==0 ? left : right;
    const NPContrastExactSum *sum=side==0 ? left_sum : right_sum;
    if(factor==0.0) continue;
    int factor_exponent;
    const double factor_mantissa=frexp(factor,&factor_exponent);
    for(int i=0;i<sum->size;++i) {
      int exponent;
      const double mantissa=frexp(sum->component[i],&exponent);
      const double high=factor_mantissa*mantissa;
      const double low=fma(factor_mantissa,mantissa,-high);
      const int shift=factor_exponent+exponent-origin;
      for(int part=0;part<2;++part) {
        const double value=part==0 ? high : low;
        if(value==0.0) continue;
        const double scaled=ldexp(value,shift);
        if(!isfinite(scaled) || scaled==0.0 || ldexp(scaled,-shift)!=value) return 0;
        np_contrast_exact_add(&difference,side==0 ? scaled : -scaled);
      }
    }
  }
  return !difference.invalid && difference.size==0;
}

typedef struct {
  NPContrastNumber off_sum;
  NPContrastNumber centered_response;
  NPResidualVariance variance;
  double anchor_response;
  int want_variance;
  int all_off_certified_zero;
  int invalid;
} NPRegressionContrastAccumulator;

typedef struct {
  double contrast;
  double standard_error;
  int structural_zero;
  NPResidualInformation information;
} NPRegressionContrastResult;

/* Compensated streaming sum: retain each high-addition error and each
 * product remainder. Renormalizing four times per donor adds no useful
 * precision to this first-order influence calculation. */
static inline void np_contrast_accumulate(NPContrastNumber *sum,NPContrastNumber term)
{
  const NPContrastNumber high=np_contrast_two_sum(sum->hi,term.hi);
  sum->hi=high.hi;
  sum->lo+=high.lo+term.lo;
}

static inline NPRegressionContrastAccumulator np_regression_contrast_begin(
  double anchor_response, int want_variance)
{
  NPRegressionContrastAccumulator result={0};
  result.anchor_response=anchor_response;
  result.want_variance=want_variance;
  result.all_off_certified_zero=1;
  result.invalid=!isfinite(anchor_response);
  return result;
}

/* The owner supplies an independently justified zero certificate. A rounded
 * zero coefficient alone does not erase unavailable donor variance. */
static inline void np_regression_contrast_add(
  NPRegressionContrastAccumulator *state, NPContrastNumber coefficient,
  double response, int certified_zero,
  NPResidualInformation information, double residual)
{
  if(!np_contrast_finite(coefficient) || !isfinite(response) ||
     (certified_zero && (coefficient.hi != 0.0 || coefficient.lo != 0.0))) {
    state->invalid=1;
    return;
  }
  state->all_off_certified_zero &= certified_zero;
  np_contrast_accumulate(&state->off_sum,coefficient);
  const NPContrastNumber centered=np_contrast_two_sum(response,-state->anchor_response);
  np_contrast_accumulate(&state->centered_response,np_contrast_multiply(coefficient,centered));
  if(state->want_variance)
    np_residual_variance_add(&state->variance,
      (long double)coefficient.hi+(long double)coefficient.lo,
      certified_zero,information,(long double)residual);
}

static inline NPRegressionContrastResult np_regression_contrast_finish(
  NPRegressionContrastAccumulator *state, NPResidualInformation anchor_information,
  double anchor_residual, int anchor_certified_zero)
{
  NPRegressionContrastResult result={NAN,NAN,0,NP_RESIDUAL_INVALID};
  if(state->invalid || !np_contrast_finite(state->off_sum) ||
     !np_contrast_finite(state->centered_response)) return result;
  result.contrast=np_contrast_value(state->centered_response);
  if(!isfinite(result.contrast)) return result;
  result.structural_zero=state->all_off_certified_zero;
  result.information=NP_RESIDUAL_IDENTIFIED;
  if(state->want_variance) {
    const NPContrastNumber anchor=anchor_certified_zero ? np_contrast_number(0.0) :
      np_contrast_negative(state->off_sum);
    np_residual_variance_add(&state->variance,
      (long double)anchor.hi+(long double)anchor.lo,
      state->all_off_certified_zero || anchor_certified_zero,
      anchor_information,(long double)anchor_residual);
    long double standard_error=0.0L;
    result.information=np_residual_variance_finish(&state->variance,&standard_error);
    if(result.information == NP_RESIDUAL_IDENTIFIED) {
      result.standard_error=(double)standard_error;
      if(!isfinite(result.standard_error)) result.information=NP_RESIDUAL_INVALID;
    }
  }
  return result;
}

#endif
