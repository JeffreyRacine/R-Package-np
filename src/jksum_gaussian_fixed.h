#ifndef NP_JKSUM_GAUSSIAN_FIXED_H
#define NP_JKSUM_GAUSSIAN_FIXED_H

#include <R_ext/Visibility.h>
#include <math.h>

/* Gaussian4/6 convolution has one arithmetic owner.  Writing the polynomial
 * in squared standardized separation makes translation invariance explicit.
 * The derivative representation K4=phi-phi''/2 and
 * K6=phi-phi''/2+phi''''/8 gives these Hermite coefficients after convolution.
 * Native callers supply unnormalized kernels; bandwidth division stays with
 * the caller.  Fixed rows prepare once; adaptive scalar calls use the same
 * coefficients for each actual pair of bandwidths. */
typedef struct {
  double variance;
  double scale;
  double coefficient[5];
} np_gaussian_convolution_polynomial;

static inline np_gaussian_convolution_polynomial
np_gaussian_convolution_prepare(const int kernel, const double hx,
                                const double hy)
{
  np_gaussian_convolution_polynomial p;
  const double hx2 = hx*hx;
  const double hy2 = hy*hy;
  const double variance = hx2 + hy2;
  const double ab = (hx2/variance)*(hy2/variance);
  p.variance = variance;
  p.scale = 0.39894228040143267794*hx*hy/sqrt(variance);
  if(kernel == 1){
    p.coefficient[0] = 1.5 + 0.75*ab;
    p.coefficient[1] = -0.5 - 1.5*ab;
    p.coefficient[2] = 0.25*ab;
    p.coefficient[3] = p.coefficient[4] = 0.0;
  } else {
    const double c6 = ab/16.0;
    const double c8 = ab*ab/64.0;
    p.coefficient[0] = 1.875 + 15.0*c6 + 105.0*c8;
    p.coefficient[1] = -1.25 - 45.0*c6 - 420.0*c8;
    p.coefficient[2] = 0.125 + 15.0*c6 + 210.0*c8;
    p.coefficient[3] = -c6 - 28.0*c8;
    p.coefficient[4] = c8;
  }
  return p;
}

static inline double np_gaussian_convolution_evaluate(
  const np_gaussian_convolution_polynomial *p, const int kernel,
  const double delta)
{
  const double t = delta*delta/p->variance;
  const double *c = p->coefficient;
  const double polynomial = kernel == 1 ? (c[2]*t + c[1])*t + c[0] :
    (((c[4]*t + c[3])*t + c[2])*t + c[1])*t + c[0];
  return p->scale*exp(-0.5*t)*polynomial;
}

attribute_hidden int np_fixed_gaussian_convolution_row_try(
  int kernel,
  const double *xt,
  int num_xt,
  int do_xw,
  double x,
  double hy,
  double h,
  double *result,
  int power);

attribute_hidden int np_fixed_gaussian_convolution_product_try(
  const int *kernels,
  double * const *xt,
  double * const *xeval,
  double * const *bandwidth,
  double * const *alt_bandwidth,
  const int *bandwidth_power,
  int ndim,
  int n,
  int eval_index,
  int bandwidth_index,
  double *result);

#endif
