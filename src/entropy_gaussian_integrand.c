#include <R.h>
#include <Rinternals.h>
#include <math.h>

#define NP_ENTROPY_INTERRUPT_INTERVAL 65536

/*
 * Evaluate the Gaussian bivariate-entropy integrand used by npdeptest() and
 * npsdeptest().  The joint product is evaluated as one bivariate Gaussian:
 * phi(z_x) * phi(z_y) = (2*pi)^-1 * exp(-(z_x^2 + z_y^2)/2).
 *
 * The caller owns all inputs.  Apart from the returned 1 x n_eval matrix, this
 * routine uses constant storage; in particular, it never materializes an
 * n_train x n_eval kernel or distance matrix.
 */
SEXP C_np_entropy_gaussian_integrand(SEXP xy,
                                     SEXP x,
                                     SEXP y,
                                     SEXP bandwidths)
{
  const double normal_constant = 0.398942280401432677939946059934;
  const double bivariate_normal_constant = normal_constant * normal_constant;
  SEXP xy_dim;
  SEXP answer;
  const double *xy_ptr;
  const double *x_ptr;
  const double *y_ptr;
  const double *bw;
  double *answer_ptr;
  double bw_joint_product;
  R_xlen_t n;
  int n_eval;
  int interrupt_countdown = NP_ENTROPY_INTERRUPT_INTERVAL;

  if (TYPEOF(xy) != REALSXP || !Rf_isMatrix(xy))
    Rf_error("entropy integration points must be a numeric matrix");

  xy_dim = Rf_getAttrib(xy, R_DimSymbol);
  if (INTEGER(xy_dim)[0] != 2)
    Rf_error("entropy integration points must have two rows");

  if (TYPEOF(x) != REALSXP || TYPEOF(y) != REALSXP)
    Rf_error("entropy integration data must be numeric vectors");

  n = XLENGTH(x);
  if (n < 1 || XLENGTH(y) != n)
    Rf_error("entropy integration data must have equal positive lengths");

  if (TYPEOF(bandwidths) != REALSXP || XLENGTH(bandwidths) != 4)
    Rf_error("entropy integration requires four numeric bandwidths");

  bw = REAL(bandwidths);
  for (int k = 0; k < 4; ++k) {
    if (!R_FINITE(bw[k]) || bw[k] <= 0.0)
      Rf_error("entropy integration bandwidths must be finite and positive");
  }

  n_eval = INTEGER(xy_dim)[1];
  answer = PROTECT(Rf_allocMatrix(REALSXP, 1, n_eval));
  answer_ptr = REAL(answer);
  xy_ptr = REAL(xy);
  x_ptr = REAL(x);
  y_ptr = REAL(y);
  bw_joint_product = bw[2] * bw[3];

  for (int j = 0; j < n_eval; ++j) {
    const double eval_x = xy_ptr[2 * j];
    const double eval_y = xy_ptr[2 * j + 1];
    long double sum_x = 0.0L;
    long double sum_y = 0.0L;
    long double sum_joint = 0.0L;

    for (R_xlen_t i = 0; i < n; ++i) {
      const double z_x = (eval_x - x_ptr[i]) / bw[0];
      const double z_y = (eval_y - y_ptr[i]) / bw[1];
      const double z_joint_x = (eval_x - x_ptr[i]) / bw[2];
      const double z_joint_y = (eval_y - y_ptr[i]) / bw[3];
      const double kernel_x = normal_constant * exp(-0.5 * z_x * z_x);
      const double kernel_y = normal_constant * exp(-0.5 * z_y * z_y);
      const double kernel_joint = bivariate_normal_constant * exp(
        -0.5 * (z_joint_x * z_joint_x + z_joint_y * z_joint_y)
      );

      sum_x += kernel_x;
      sum_y += kernel_y;
      sum_joint += kernel_joint / bw_joint_product;

      if (--interrupt_countdown == 0) {
        R_CheckUserInterrupt();
        interrupt_countdown = NP_ENTROPY_INTERRUPT_INTERVAL;
      }
    }

    {
      const double density_x =
        ((double)(sum_x / (long double)n)) / bw[0];
      const double density_y =
        ((double)(sum_y / (long double)n)) / bw[1];
      const double density_joint =
        (double)(sum_joint / (long double)n);
      const double delta = sqrt(density_joint) -
        sqrt(density_x) * sqrt(density_y);
      answer_ptr[j] = delta * delta;
    }
  }

  UNPROTECT(1);
  return answer;
}
