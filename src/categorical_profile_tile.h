#ifndef NP_CATEGORICAL_PROFILE_TILE_H
#define NP_CATEGORICAL_PROFILE_TILE_H

#include <stddef.h>

#define NP_CATEGORICAL_PROFILE_TILE_MAX_BYTES \
  ((size_t)64U*(size_t)1024U*(size_t)1024U)

typedef enum {
  NP_PROFILE_TILE_OK = 0,
  NP_PROFILE_TILE_ERR_ARGUMENT = 1,
  NP_PROFILE_TILE_ERR_DIMENSION = 2,
  NP_PROFILE_TILE_ERR_KERNEL = 3,
  NP_PROFILE_TILE_ERR_OPERATOR = 4,
  NP_PROFILE_TILE_ERR_NONFINITE = 5,
  NP_PROFILE_TILE_ERR_SUPPORT = 6,
  NP_PROFILE_TILE_ERR_RANGE = 7,
  NP_PROFILE_TILE_ERR_CAPACITY = 8
} NPCategoricalProfileTileStatus;

typedef struct {
  int ntrain;
  int neval;
  int nunordered;
  int nordered;
  double * const *train_unordered;
  double * const *train_ordered;
  double * const *eval_unordered;
  double * const *eval_ordered;
  const int *kernel_unordered;
  const int *kernel_ordered;
  const int *operator_code;
  const double *lambda;
  const int *num_categories;
  double * const *category_values;
} NPCategoricalProfileKernelSpec;

NPCategoricalProfileTileStatus
np_categorical_profile_tile_bytes(size_t ntrain,
                                  size_t nrows,
                                  size_t *elements,
                                  size_t *bytes);

/*
 * Validate immutable metadata and training profiles once before a tiled
 * traversal. A caller using the prevalidated fill must not mutate the
 * specification or its training storage until the traversal is complete.
 */
NPCategoricalProfileTileStatus
np_categorical_profile_spec_validate(
  const NPCategoricalProfileKernelSpec *spec);

/*
 * Fill one tile after np_categorical_profile_spec_validate() succeeds.
 * Evaluation values are still checked for the requested slice, and capacity
 * and range checks remain mandatory on every call.
 */
NPCategoricalProfileTileStatus
np_categorical_profile_tile_fill_prevalidated(
  const NPCategoricalProfileKernelSpec *spec,
  int eval_start,
  int eval_count,
  double *output,
  size_t output_elements);

/* Fully checked convenience entry for isolated calls and proof hooks. */
NPCategoricalProfileTileStatus
np_categorical_profile_tile_fill(
  const NPCategoricalProfileKernelSpec *spec,
  int eval_start,
  int eval_count,
  double *output,
  size_t output_elements);

const char *np_categorical_profile_tile_status_message(
  NPCategoricalProfileTileStatus status);

#endif
