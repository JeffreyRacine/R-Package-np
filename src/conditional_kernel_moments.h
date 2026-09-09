#ifndef NP_CONDITIONAL_KERNEL_MOMENTS_H
#define NP_CONDITIONAL_KERNEL_MOMENTS_H

/* Scalar conditional leading variance uses independent X/Y kernel products.
 * Keep the legacy initializer (including its finite-shift constants) intact
 * for every other consumer. This helper allocates nothing and runs once/fit.
 * headers.h supplies initialize_kernel_regression_asymptotic_constants(). */
static double np_conditional_kernel_square_product(
  const int kernel_x, const int kernel_y,
  const int dimension_x, const int dimension_y, const int density)
{
  double integral, product_x = 1.0, product_y = 1.0, half, difference;

  if(density && kernel_x == kernel_y) {
    if(dimension_x + dimension_y > 0)
      initialize_kernel_regression_asymptotic_constants(
        kernel_x, dimension_x + dimension_y,
        &integral, &product_x, &half, &difference);
    return product_x;
  }
  if(dimension_x > 0)
    initialize_kernel_regression_asymptotic_constants(
      kernel_x, dimension_x, &integral, &product_x, &half, &difference);
  if(density && dimension_y > 0)
    initialize_kernel_regression_asymptotic_constants(
      kernel_y, dimension_y, &integral, &product_y, &half, &difference);
  return product_x * product_y;
}

/* Integral of the square of the analytic derivative implemented by
 * np_deriv_gauss{2,4,6,8}/np_deriv_epan{2,4,6,8}. Values use the actual
 * polynomial coefficients and Epanechnikov support (-sqrt(5),sqrt(5)).
 * Independent polynomial/quadrature checks are in the SE theory packet.
 * Uniform has no ordinary smooth-kernel derivative theorem and is excluded.
 */
static double np_conditional_kernel_derivative_square_integral(const int kernel)
{
  static const double moment[8] = {
    0.14104739588693907, 0.48485042336135320,
    0.97865405056318600, 1.5972429562844692,
    0.13416407864998742, 0.83852549282614142,
    2.5679843227599548, 5.7779647262104845
  };
  return kernel >= 0 && kernel < 8 ? moment[kernel] : NAN;
}

#endif
