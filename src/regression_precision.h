#ifndef NP_REGRESSION_PRECISION_H
#define NP_REGRESSION_PRECISION_H
#include "jksum_lp_solve.h"
/* Invocation-owned accepted-policy QR. No rank choice, retry or ridge change. */
typedef struct {
  int n, p, m, work_capacity, work_m, work_p;
  size_t matrix_capacity, vector_capacity;
  double *a, *vectors, *work;
  /* Optional exact-zero donor map for the Epanechnikov HC0 variance pass. */
  int compact_zero, active_count, active_capacity;
  int *active_donors;
  int *pivot;
  int pivot_capacity;
  double *tau, *scale, *sw, *u, *row;
  double max_weight, root_ridge, anchor;
  /* The seven-column cubic-by-linear GLP may have one algebraic alias
   * on three rounded abscissae. All other precision owners leave this off. */
  int check_cubic_support, cubic_alias;
  int local_cubic;
  double cubic_center, cubic_scale, *local_columns;
  size_t local_capacity;
  int cubic_columns[4];
  double cubic_levels[3][4];
  const double *weights;
} NPRegressionPrecision;
void np_regression_precision_init(NPRegressionPrecision *q);
void np_regression_precision_clear(NPRegressionPrecision *q);
int np_regression_precision_needed(NPLPSolveWorkspace *s, int p);
int np_regression_precision_factor(NPRegressionPrecision *q,
  double *const *basis, const double *weights, int n, int p,
  double ridge, double anchor);
int np_regression_precision_row(NPRegressionPrecision *q,
  const double *direction);
double np_regression_precision_apply(const NPRegressionPrecision *q,
  const double *response);
#endif
