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

#endif
