#ifndef NP_REGRESSION_PRECISION_H
#define NP_REGRESSION_PRECISION_H
#include "jksum_lp_solve.h"
/* Invocation-owned accepted-policy QR. No rank choice, retry or ridge change. */
typedef struct {
  int n, p, m, work_capacity, work_m, work_p;
  size_t matrix_capacity, vector_capacity;
  double *a, *vectors, *work;
  int *pivot;
  int pivot_capacity;
  double *tau, *scale, *sw, *u, *row;
  double max_weight, root_ridge, anchor;
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
