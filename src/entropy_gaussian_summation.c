#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <math.h>

#define NP_ENTROPY_SUMMATION_INTERRUPT_INTERVAL 65536

/*
 * Evaluate fixed, second-order Gaussian symmetry statistics from bootstrap
 * multiplicities.  The support is shared by all columns of counts.  Work is
 * streamed over active support points, so storage is O(n * chunk), never
 * O(n^2).  Gaussian normalizing constants and the common bandwidth cancel in
 * the density ratio.
 */
SEXP C_np_entropy_symmetric_summation_counts(SEXP support,
                                              SEXP counts,
                                              SEXP bandwidth)
{
  SEXP dim;
  SEXP answer;
  const double *x;
  const double *weights;
  const double h = Rf_asReal(bandwidth);
  double *out;
  int n_support;
  int n_replications;
  int interrupt_countdown = NP_ENTROPY_SUMMATION_INTERRUPT_INTERVAL;

  if (TYPEOF(support) != REALSXP || TYPEOF(counts) != REALSXP ||
      !Rf_isMatrix(counts))
    Rf_error("entropy symmetry counts require numeric inputs");

  dim = Rf_getAttrib(counts, R_DimSymbol);
  n_support = INTEGER(dim)[0];
  n_replications = INTEGER(dim)[1];
  if (n_support < 1 || n_replications < 1 ||
      XLENGTH(support) != n_support)
    Rf_error("entropy symmetry counts have incompatible dimensions");
  if (!R_FINITE(h) || h <= 0.0)
    Rf_error("entropy symmetry bandwidth must be finite and positive");

  x = REAL(support);
  for (int i = 0; i < n_support; ++i)
    if (!R_FINITE(x[i]))
      Rf_error("entropy symmetry support must be finite");

  answer = PROTECT(Rf_allocVector(REALSXP, n_replications));
  weights = REAL(counts);
  out = REAL(answer);

  for (int replication = 0; replication < n_replications; ++replication) {
    const double *replication_weights =
      weights + ((R_xlen_t)n_support * replication);
    long double total_weight = 0.0L;
    long double weighted_sum = 0.0L;
    long double statistic_sum = 0.0L;

    for (int i = 0; i < n_support; ++i) {
      const double weight = replication_weights[i];
      if (!R_FINITE(weight) || weight < 0.0)
        Rf_error("entropy symmetry counts must be finite and nonnegative");
      if (weight > 0.0) {
        total_weight += weight;
        weighted_sum += weight * x[i];
      }
    }
    if (total_weight <= 0.0L)
      Rf_error("entropy symmetry counts must have positive column sums");

    {
      const double location = (double)(weighted_sum / total_weight);

      for (int evaluation = 0; evaluation < n_support; ++evaluation) {
        long double density_sum = 0.0L;
        long double rotated_sum = 0.0L;
        const double evaluation_weight = replication_weights[evaluation];

        if (evaluation_weight <= 0.0)
          continue;

        for (int training = 0; training < n_support; ++training) {
          const double training_weight = replication_weights[training];
          if (training_weight > 0.0) {
            const double z =
              (x[evaluation] - x[training]) / h;
            const double z_rotated =
              (2.0 * location - x[evaluation] - x[training]) / h;
            density_sum += training_weight * exp(-0.5 * z * z);
            rotated_sum +=
              training_weight * exp(-0.5 * z_rotated * z_rotated);
          }

          if (--interrupt_countdown == 0) {
            R_CheckUserInterrupt();
            interrupt_countdown = NP_ENTROPY_SUMMATION_INTERRUPT_INTERVAL;
          }
        }

        {
          const double ratio = (double)(rotated_sum / density_sum);
          const double term = 1.0 - sqrt(ratio);
          statistic_sum += evaluation_weight * term * term;
        }
      }
    }

    out[replication] =
      0.5 * (double)(statistic_sum / total_weight);
  }

  UNPROTECT(1);
  return answer;
}

/*
 * Fuse the three Gaussian density sums in the bivariate Srho summation
 * statistic.  This avoids three general npksum() setups and makes one streamed
 * pass over each training/evaluation pair with constant auxiliary storage.
 */
SEXP C_np_entropy_bivariate_summation(SEXP x,
                                      SEXP y,
                                      SEXP bandwidths)
{
  const double *x_ptr;
  const double *y_ptr;
  const double *bw;
  const R_xlen_t n = XLENGTH(x);
  long double statistic_sum = 0.0L;
  int interrupt_countdown = NP_ENTROPY_SUMMATION_INTERRUPT_INTERVAL;

  if (TYPEOF(x) != REALSXP || TYPEOF(y) != REALSXP ||
      TYPEOF(bandwidths) != REALSXP || n < 1 || XLENGTH(y) != n ||
      XLENGTH(bandwidths) != 4)
    Rf_error("bivariate entropy summation requires compatible numeric inputs");

  x_ptr = REAL(x);
  y_ptr = REAL(y);
  bw = REAL(bandwidths);
  for (int k = 0; k < 4; ++k)
    if (!R_FINITE(bw[k]) || bw[k] <= 0.0)
      Rf_error("bivariate entropy bandwidths must be finite and positive");
  for (R_xlen_t i = 0; i < n; ++i)
    if (!R_FINITE(x_ptr[i]) || !R_FINITE(y_ptr[i]))
      Rf_error("bivariate entropy data must be finite");

  for (R_xlen_t evaluation = 0; evaluation < n; ++evaluation) {
    long double sum_x = 0.0L;
    long double sum_y = 0.0L;
    long double sum_joint = 0.0L;

    for (R_xlen_t training = 0; training < n; ++training) {
      const double z_x =
        (x_ptr[evaluation] - x_ptr[training]) / bw[0];
      const double z_y =
        (y_ptr[evaluation] - y_ptr[training]) / bw[1];
      const double z_joint_x =
        (x_ptr[evaluation] - x_ptr[training]) / bw[2];
      const double z_joint_y =
        (y_ptr[evaluation] - y_ptr[training]) / bw[3];

      sum_x += exp(-0.5 * z_x * z_x);
      sum_y += exp(-0.5 * z_y * z_y);
      sum_joint += exp(
        -0.5 * (z_joint_x * z_joint_x + z_joint_y * z_joint_y)
      );

      if (--interrupt_countdown == 0) {
        R_CheckUserInterrupt();
        interrupt_countdown = NP_ENTROPY_SUMMATION_INTERRUPT_INTERVAL;
      }
    }

    {
      const long double ratio =
        (sum_x * sum_y * bw[2] * bw[3]) /
        ((long double)n * sum_joint * bw[0] * bw[1]);
      const double term = 1.0 - sqrt((double)ratio);
      statistic_sum += term * term;
    }
  }

  return Rf_ScalarReal(0.5 * (double)(statistic_sum / n));
}

/*
 * Evaluate a bounded chunk of dependence-test bootstrap assignments.  Each
 * integer row maps the fixed y positions to resampled x observations.
 * Looping over position pairs first reuses the invariant y marginal work
 * across the chunk while retaining the scalar training order for every
 * replication.  Workspace is O(chunk), in addition to the caller-owned
 * O(n * chunk) index plan.
 */
SEXP C_np_entropy_bivariate_summation_xindex(SEXP x,
                                             SEXP y,
                                             SEXP index,
                                             SEXP bandwidths)
{
  SEXP dim;
  SEXP answer;
  const double *x_ptr;
  const double *y_ptr;
  const double *bw;
  const int *index_ptr;
  double *out;
  long double *sum_x;
  long double *sum_joint;
  long double *statistic_sum;
  const R_xlen_t n = XLENGTH(x);
  int n_replications;
  int interrupt_countdown = NP_ENTROPY_SUMMATION_INTERRUPT_INTERVAL;

  if (TYPEOF(x) != REALSXP || TYPEOF(y) != REALSXP ||
      TYPEOF(index) != INTSXP || !Rf_isMatrix(index) ||
      TYPEOF(bandwidths) != REALSXP || n < 1 || XLENGTH(y) != n ||
      XLENGTH(bandwidths) != 4 || n > INT_MAX)
    Rf_error("bivariate entropy index batch requires compatible inputs");

  dim = Rf_getAttrib(index, R_DimSymbol);
  if (INTEGER(dim)[0] < 1 || (R_xlen_t)INTEGER(dim)[1] != n)
    Rf_error("bivariate entropy index batch has incompatible dimensions");
  n_replications = INTEGER(dim)[0];

  x_ptr = REAL(x);
  y_ptr = REAL(y);
  index_ptr = INTEGER(index);
  bw = REAL(bandwidths);
  for (int k = 0; k < 4; ++k)
    if (!R_FINITE(bw[k]) || bw[k] <= 0.0)
      Rf_error("bivariate entropy bandwidths must be finite and positive");
  for (R_xlen_t i = 0; i < n; ++i)
    if (!R_FINITE(x_ptr[i]) || !R_FINITE(y_ptr[i]))
      Rf_error("bivariate entropy data must be finite");
  for (R_xlen_t i = 0; i < XLENGTH(index); ++i)
    if (index_ptr[i] < 1 || index_ptr[i] > n)
      Rf_error("bivariate entropy bootstrap indices are out of range");

  answer = PROTECT(Rf_allocVector(REALSXP, n_replications));
  out = REAL(answer);
  sum_x = (long double *)R_alloc((size_t)n_replications,
                                 sizeof(long double));
  sum_joint = (long double *)R_alloc((size_t)n_replications,
                                     sizeof(long double));
  statistic_sum = (long double *)R_alloc((size_t)n_replications,
                                         sizeof(long double));
  for (int replication = 0; replication < n_replications; ++replication)
    statistic_sum[replication] = 0.0L;

  for (R_xlen_t evaluation = 0; evaluation < n; ++evaluation) {
    long double sum_y = 0.0L;

    for (int replication = 0; replication < n_replications; ++replication) {
      sum_x[replication] = 0.0L;
      sum_joint[replication] = 0.0L;
    }

    for (R_xlen_t training = 0; training < n; ++training) {
      const double z_y =
        (y_ptr[evaluation] - y_ptr[training]) / bw[1];
      const double z_joint_y =
        (y_ptr[evaluation] - y_ptr[training]) / bw[3];
      const double z_joint_y_squared = z_joint_y * z_joint_y;

      sum_y += exp(-0.5 * z_y * z_y);

      for (int replication = 0; replication < n_replications; ++replication) {
        const double x_evaluation =
          x_ptr[index_ptr[replication +
                          (R_xlen_t)n_replications * evaluation] - 1];
        const double x_training =
          x_ptr[index_ptr[replication +
                          (R_xlen_t)n_replications * training] - 1];
        const double z_x = (x_evaluation - x_training) / bw[0];
        const double z_joint_x =
          (x_evaluation - x_training) / bw[2];

        sum_x[replication] += exp(-0.5 * z_x * z_x);
        sum_joint[replication] += exp(
          -0.5 * (z_joint_x * z_joint_x + z_joint_y_squared)
        );

        if (--interrupt_countdown == 0) {
          R_CheckUserInterrupt();
          interrupt_countdown = NP_ENTROPY_SUMMATION_INTERRUPT_INTERVAL;
        }
      }
    }

    for (int replication = 0; replication < n_replications; ++replication) {
      const long double ratio =
        (sum_x[replication] * sum_y * bw[2] * bw[3]) /
        ((long double)n * sum_joint[replication] * bw[0] * bw[1]);
      const double term = 1.0 - sqrt((double)ratio);
      statistic_sum[replication] += term * term;
    }
  }

  for (int replication = 0; replication < n_replications; ++replication)
    out[replication] =
      0.5 * (double)(statistic_sum[replication] / n);

  UNPROTECT(1);
  return answer;
}

/*
 * Evaluate a chunk of two-sample univariate Srho statistics from bootstrap
 * multiplicities on a common support.  The x multiplicities also identify the
 * evaluation sample.  Storage is owned by the caller and bounded by its chunk.
 */
SEXP C_np_entropy_univariate_summation_counts(SEXP support,
                                               SEXP counts_x,
                                               SEXP counts_y,
                                               SEXP bandwidths)
{
  SEXP dim_x;
  SEXP dim_y;
  SEXP answer;
  const double *support_ptr;
  const double *weights_x;
  const double *weights_y;
  const double *bw;
  double *out;
  int n_support;
  int n_replications;
  int interrupt_countdown = NP_ENTROPY_SUMMATION_INTERRUPT_INTERVAL;

  if (TYPEOF(support) != REALSXP || TYPEOF(counts_x) != REALSXP ||
      TYPEOF(counts_y) != REALSXP || TYPEOF(bandwidths) != REALSXP ||
      !Rf_isMatrix(counts_x) || !Rf_isMatrix(counts_y) ||
      XLENGTH(bandwidths) != 2)
    Rf_error("univariate entropy counts require numeric inputs");

  dim_x = PROTECT(Rf_getAttrib(counts_x, R_DimSymbol));
  dim_y = PROTECT(Rf_getAttrib(counts_y, R_DimSymbol));
  n_support = INTEGER(dim_x)[0];
  n_replications = INTEGER(dim_x)[1];
  if (n_support < 1 || n_replications < 1 ||
      INTEGER(dim_y)[0] != n_support ||
      INTEGER(dim_y)[1] != n_replications ||
      XLENGTH(support) != n_support)
    Rf_error("univariate entropy counts have incompatible dimensions");

  support_ptr = REAL(support);
  weights_x = REAL(counts_x);
  weights_y = REAL(counts_y);
  bw = REAL(bandwidths);
  for (int k = 0; k < 2; ++k)
    if (!R_FINITE(bw[k]) || bw[k] <= 0.0)
      Rf_error("univariate entropy bandwidths must be finite and positive");
  for (int i = 0; i < n_support; ++i)
    if (!R_FINITE(support_ptr[i]))
      Rf_error("univariate entropy support must be finite");

  answer = PROTECT(Rf_allocVector(REALSXP, n_replications));
  out = REAL(answer);

  for (int replication = 0; replication < n_replications; ++replication) {
    const double *replication_x =
      weights_x + ((R_xlen_t)n_support * replication);
    const double *replication_y =
      weights_y + ((R_xlen_t)n_support * replication);
    long double size_x = 0.0L;
    long double size_y = 0.0L;
    long double statistic_sum = 0.0L;

    for (int i = 0; i < n_support; ++i) {
      if (!R_FINITE(replication_x[i]) || replication_x[i] < 0.0 ||
          !R_FINITE(replication_y[i]) || replication_y[i] < 0.0)
        Rf_error("univariate entropy counts must be finite and nonnegative");
      size_x += replication_x[i];
      size_y += replication_y[i];
    }
    if (size_x <= 0.0L || size_y <= 0.0L)
      Rf_error("univariate entropy counts must have positive column sums");

    for (int evaluation = 0; evaluation < n_support; ++evaluation) {
      long double density_x = 0.0L;
      long double density_y = 0.0L;
      const double evaluation_weight = replication_x[evaluation];

      if (evaluation_weight <= 0.0)
        continue;

      for (int training = 0; training < n_support; ++training) {
        if (replication_x[training] > 0.0) {
          const double z =
            (support_ptr[evaluation] - support_ptr[training]) / bw[0];
          density_x += replication_x[training] * exp(-0.5 * z * z);
        }
        if (replication_y[training] > 0.0) {
          const double z =
            (support_ptr[evaluation] - support_ptr[training]) / bw[1];
          density_y += replication_y[training] * exp(-0.5 * z * z);
        }

        if (--interrupt_countdown == 0) {
          R_CheckUserInterrupt();
          interrupt_countdown = NP_ENTROPY_SUMMATION_INTERRUPT_INTERVAL;
        }
      }

      {
        const long double ratio =
          (density_y * size_x * bw[0]) /
          (density_x * size_y * bw[1]);
        const double term = 1.0 - sqrt((double)ratio);
        statistic_sum += evaluation_weight * term * term;
      }
    }

    out[replication] = 0.5 * (double)(statistic_sum / size_x);
  }

  UNPROTECT(3);
  return answer;
}
