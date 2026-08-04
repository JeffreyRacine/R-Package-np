#ifndef NP_BETA_SCALED_ROW_H
#define NP_BETA_SCALED_ROW_H

#include <stddef.h>

#include "categorical_profile_tile.h"
#include "continuous_kernel_row.h"

/*
 * Persistent O(n) owner for one canonical beta kernel row.  The row is kept
 * on a single positive common scale and may include the package's ordinary
 * categorical product.  Estimator consumers own the supplied row buffer and
 * restore absolute units only after completing their scale-equivariant
 * algebra.
 */
typedef struct {
  NPCategoricalProfileKernelSpec dense_spec;
  NPCategoricalProfileKernelSpec compressed_spec;
  int use_compressed;
  int *profile_id;
  double **compressed_unordered;
  double **compressed_ordered;
  double *compressed_unordered_storage;
  double *compressed_ordered_storage;
  double *compressed_log_absolute;
  signed char *compressed_sign;
  double *scratch;
  size_t scratch_capacity;
} NPBetaScaledRowCategoricalContext;

typedef struct {
  int ready;
  int has_categories;
  NPContinuousKernelBetaPreparedContext beta_prepared;
  NPContinuousKernelRowWorkspace row_workspace;
  NPContinuousKernelRowPlan plan;
  NPContinuousKernelRowResult row_result;
  NPBetaScaledRowCategoricalContext categorical_context;
  NPContinuousKernelLogFactorProvider categorical_provider;
  NPContinuousKernelDerivativeDiagnostics *diagnostics;
} NPBetaScaledRowContext;

void np_beta_scaled_row_context_init(NPBetaScaledRowContext *context);
void np_beta_scaled_row_context_clear(NPBetaScaledRowContext *context);

NPContinuousKernelRowStatus np_beta_scaled_row_context_prepare(
  NPBetaScaledRowContext *context,
  const NPContinuousKernelRoute *route,
  NPContinuousKernelDerivativeDiagnostics *diagnostics,
  int bandwidth_mode,
  int ntrain,
  int neval,
  int ncontinuous,
  int nunordered,
  int nordered,
  double **train_continuous,
  double **eval_continuous,
  double **train_unordered,
  double **eval_unordered,
  double **train_ordered,
  double **eval_ordered,
  double **bandwidth_train,
  double **bandwidth_eval,
  const int *operator_code,
  const int *kernel_unordered,
  const int *kernel_ordered,
  const double *lambda,
  const int *num_categories,
  double **category_values,
  int categorical_compress,
  double *row);

NPContinuousKernelRowStatus np_beta_scaled_row_context_fill(
  NPBetaScaledRowContext *context,
  int evaluation,
  double *sum,
  double *common_log_scale);

NPContinuousKernelRowStatus np_beta_scaled_row_context_fill_omitting(
  NPBetaScaledRowContext *context,
  int evaluation,
  int omitted_observation,
  double *sum,
  double *common_log_scale);

#endif
