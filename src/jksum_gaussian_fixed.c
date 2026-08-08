#include <math.h>

#include "headers.h"
#include "jksum_gaussian_fixed.h"

attribute_hidden int np_fixed_gaussian_convolution_row_try(
  const int kernel,
  const double * const xt,
  const int num_xt,
  const int do_xw,
  const double x,
  const double hy,
  const double h,
  double * const result,
  const int power)
{
  const int multiply_existing = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = multiply_existing ? result : &unit_weight;
  int i;
  int j;

  if((xt == NULL) || (result == NULL) || (num_xt < 0) ||
     ((kernel != 1) && (kernel != 2)))
    return 0;

  if(kernel == 1){
    const double hx2 = h*h;
    const double hy2 = hy*hy;
    const double hxy2 = hx2+hy2;
    const double x2 = x*x;
    const double hx3 = hx2*h;
    const double hy3 = hy2*hy;
    const double hx5 = hx3*hx2;
    const double hy5 = hy3*hy2;
    const double hx7 = hx5*hx2;
    const double hy7 = hy5*hy2;
    const double hx9 = hx7*hx2;
    const double hy9 = hy7*hy2;
    const double denominator =
      sqrt(hy2 + hx2)*4*hxy2*hxy2*hxy2*hxy2;
    const double hy_power = ipow(hy, power);
    const double scale = ONE_OVER_SQRT_TWO_PI/(denominator*hy_power);

    for(i = 0, j = 0; i < num_xt; i++, j += multiply_existing){
      const double y = xt[i];
      const double y2 = y*y;
      const double delta = x-y;
      double kval;

      if(xw[j] == 0.0)
        continue;

      kval =
        (hx3*hy3*(y2*y2 - 4*x*y*y2 + x2*x2)
         + (6*hx3*hy3*x2 - 2*h*hy7 - 6*hx3*hy5 -
            12*hx5*hy3 - 2*hx7*hy)*y2
         + ((4*h*hy7 + 24*hx3*hy5 + 24*hx5*hy3 +
             4*hx7*hy)*x - 4*hx3*hy3*x2*x)*y
         + (-2*h*hy7 - 12*hx3*hy5 - 12*hx5*hy3 -
            2*hx7*hy)*x2
         + 6*h*hy9 + 27*hx3*hy7 + 42*hx5*hy5 +
           27*hx7*hy3 + 6*hx9*hy)*
        exp(-0.5*delta*delta/hxy2)*scale;
      result[i] = xw[j]*kval;
    }
    return 1;
  }

  {
    const double hx2 = h*h;
    const double hx4 = hx2*hx2;
    const double hx6 = hx4*hx2;
    const double hx8 = hx4*hx4;
    const double hx10 = hx8*hx2;
    const double hx12 = hx10*hx2;
    const double hx14 = hx12*hx2;
    const double hx16 = hx8*hx8;
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x2*x2;
    const double x5 = x*x4;
    const double x6 = x3*x3;
    const double x7 = x*x6;
    const double x8 = x4*x4;
    const double hy2 = hy*hy;
    const double hy4 = hy2*hy2;
    const double hy6 = hy4*hy2;
    const double hy8 = hy4*hy4;
    const double hy10 = hy8*hy2;
    const double hy12 = hy10*hy2;
    const double hy14 = hy12*hy2;
    const double hy16 = hy8*hy8;
    const double hxy2 = hx2+hy2;
    const double hxy4 = hxy2*hxy2;
    const double hxy8 = hxy4*hxy4;
    const double denominator = sqrt(hxy2)*64*hxy8;
    const double hy_power = ipow(hy, power);
    const double scale = ONE_OVER_SQRT_TWO_PI/(denominator*hy_power);

    for(i = 0, j = 0; i < num_xt; i++, j += multiply_existing){
      const double y = xt[i];
      const double y2 = y*y;
      const double y3 = y*y2;
      const double y4 = y2*y2;
      const double y5 = y*y4;
      const double y6 = y3*y3;
      const double y7 = y*y6;
      const double y8 = y4*y4;
      const double delta = x-y;
      double kval;

      if(xw[j] == 0.0)
        continue;

      kval =
        h*hy*
        (hx4*hy4*y8-8*hx4*hy4*x*y7+28*hx4*hy4*x2*y6
         -4*hx2*hy8*y6-40*hx4*hy6*y6-40*hx6*hy4*y6
         -4*hx8*hy2*y6-56*hx4*hy4*x3*y5+24*hx2*hy8*x*y5
         +240*hx4*hy6*x*y5+240*hx6*hy4*x*y5+24*hx8*hy2*x*y5
         +70*hx4*hy4*x4*y4-60*hx2*hy8*x2*y4
         -600*hx4*hy6*x2*y4-600*hx6*hy4*x2*y4
         -60*hx8*hy2*x2*y4+8*hy12*y4+108*hx2*hy10*y4
         +570*hx4*hy8*y4+940*hx6*hy6*y4+570*hx8*hy4*y4
         +108*hx10*hy2*y4+8*hx12*y4-56*hx4*hy4*x5*y3
         +80*hx2*hy8*x3*y3+800*hx4*hy6*x3*y3
         +800*hx6*hy4*x3*y3+80*hx8*hy2*x3*y3
         -32*hy12*x*y3-432*hx2*hy10*x*y3
         -2280*hx4*hy8*x*y3-3760*hx6*hy6*x*y3
         -2280*hx8*hy4*x*y3-432*hx10*hy2*x*y3
         -32*hx12*x*y3+28*hx4*hy4*x6*y2
         -60*hx2*hy8*x4*y2-600*hx4*hy6*x4*y2
         -600*hx6*hy4*x4*y2-60*hx8*hy2*x4*y2
         +48*hy12*x2*y2+648*hx2*hy10*x2*y2
         +3420*hx4*hy8*x2*y2+5640*hx6*hy6*x2*y2
         +3420*hx8*hy4*x2*y2+648*hx10*hy2*x2*y2
         +48*hx12*x2*y2-80*hy14*y2-740*hx2*hy12*y2
         -3000*hx4*hy10*y2-5860*hx6*hy8*y2
         -5860*hx8*hy6*y2-3000*hx10*hy4*y2
         -740*hx12*hy2*y2-80*hx14*y2-8*hx4*hy4*x7*y
         +24*hx2*hy8*x5*y+240*hx4*hy6*x5*y
         +240*hx6*hy4*x5*y+24*hx8*hy2*x5*y
         -32*hy12*x3*y-432*hx2*hy10*x3*y
         -2280*hx4*hy8*x3*y-3760*hx6*hy6*x3*y
         -2280*hx8*hy4*x3*y-432*hx10*hy2*x3*y
         -32*hx12*x3*y+160*hy14*x*y+1480*hx2*hy12*x*y
         +6000*hx4*hy10*x*y+11720*hx6*hy8*x*y
         +11720*hx8*hy6*x*y+6000*hx10*hy4*x*y
         +1480*hx12*hy2*x*y+160*hx14*x*y+hx4*hy4*x8
         -4*hx2*hy8*x6-40*hx4*hy6*x6-40*hx6*hy4*x6
         -4*hx8*hy2*x6+8*hy12*x4+108*hx2*hy10*x4
         +570*hx4*hy8*x4+940*hx6*hy6*x4+570*hx8*hy4*x4
         +108*hx10*hy2*x4+8*hx12*x4-80*hy14*x2
         -740*hx2*hy12*x2-3000*hx4*hy10*x2
         -5860*hx6*hy8*x2-5860*hx8*hy6*x2
         -3000*hx10*hy4*x2-740*hx12*hy2*x2
         -80*hx14*x2+120*hy16+1020*hx2*hy14
         +3825*hx4*hy12+8040*hx6*hy10+10230*hx8*hy8
         +8040*hx10*hy6+3825*hx12*hy4+1020*hx14*hy2
         +120*hx16)*
        exp(-0.5*delta*delta/hxy2)*scale;
      result[i] = xw[j]*kval;
    }
  }

  return 1;
}

attribute_hidden int np_fixed_gaussian_convolution_product_try(
  const int * const kernels,
  double * const * const xt,
  double * const * const xeval,
  double * const * const bandwidth,
  double * const * const alt_bandwidth,
  const int * const bandwidth_power,
  const int ndim,
  const int n,
  const int eval_index,
  const int bandwidth_index,
  double * const result)
{
  int dimension;

  if((kernels == NULL) || (xt == NULL) || (xeval == NULL) ||
     (bandwidth == NULL) || (alt_bandwidth == NULL) ||
     (bandwidth_power == NULL) || (result == NULL) ||
     (ndim <= 0) || (n < 0))
    return 0;

  for(dimension = 0; dimension < ndim; dimension++){
    if(((kernels[dimension] != 1) && (kernels[dimension] != 2)) ||
       (xt[dimension] == NULL) || (xeval[dimension] == NULL) ||
       (bandwidth[dimension] == NULL) ||
       (alt_bandwidth[dimension] == NULL))
      return 0;
  }

  for(dimension = 0; dimension < ndim; dimension++){
    if(!np_fixed_gaussian_convolution_row_try(
         kernels[dimension], xt[dimension], n, dimension != 0,
         xeval[dimension][eval_index], alt_bandwidth[dimension][0],
         bandwidth[dimension][bandwidth_index], result,
         bandwidth_power[dimension]))
      return 0;
  }

  return 1;
}
