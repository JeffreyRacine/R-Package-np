#ifndef NP_REGRESSION_INFERENCE_REUSE_H
#define NP_REGRESSION_INFERENCE_REUSE_H

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Successful uncertainty results only, scoped to one immutable basis/residual
 * invocation. Mutable retained rows and accepted adjoints are value-keyed.
 * Each direction keeps its last result; no cross-response or MPI state. */
typedef struct {
  double kernel_value;
  double divisor;
  double standard_error;
  int unavailable;
  int valid;
} NPInferenceReuseEntry;

typedef struct {
  double *projection;
  NPInferenceReuseEntry *entry;
  double **basis;
  const void *context;
  int nterms;
  int ndirections;
  int ntrain;
  int row_valid;
  double row_value;
  double row_divisor;
} NPInferenceReuse;

static inline void np_inference_reuse_clear(NPInferenceReuse *reuse)
{
  free(reuse->projection);
  free(reuse->entry);
  memset(reuse, 0, sizeof(*reuse));
}

static inline int np_inference_reuse_reserve(
  NPInferenceReuse *reuse, int nterms, int ndirections, int ntrain,
  double **basis, const void *context)
{
  if(reuse == NULL || reuse->projection != NULL || reuse->entry != NULL ||
     nterms <= 0 || ndirections <= 0 || ntrain <= 0 ||
     basis == NULL || context == NULL ||
     (size_t)nterms > SIZE_MAX / sizeof(double) / (size_t)ndirections ||
     (size_t)ndirections > SIZE_MAX / sizeof(NPInferenceReuseEntry))
    return 0;
  reuse->projection = (double *)malloc(
    (size_t)nterms * (size_t)ndirections * sizeof(double));
  reuse->entry = (NPInferenceReuseEntry *)calloc(
    (size_t)ndirections, sizeof(NPInferenceReuseEntry));
  if(reuse->projection == NULL || reuse->entry == NULL) {
    np_inference_reuse_clear(reuse);
    return 0;
  }
  reuse->basis = basis;
  reuse->context = context;
  reuse->nterms = nterms;
  reuse->ndirections = ndirections;
  reuse->ntrain = ntrain;
  return 1;
}

/* Called once per actual evaluation row, outside the direction loop. Failure
 * to prove a constant row is a cache miss, not an estimator failure/fallback.
 * The unchanged donor owner still validates and computes that row. */
static inline void np_inference_reuse_begin_row(
  NPInferenceReuse *reuse, double **basis, const void *context,
  const double *row, int ntrain, double divisor, void (*activity)(void))
{
  reuse->row_valid = 0;
  if(reuse->entry == NULL || basis != reuse->basis ||
     context != reuse->context || ntrain != reuse->ntrain || row == NULL ||
     !isfinite(divisor) || divisor <= 0.0 || !isfinite(row[0]))
    return;
  for(int donor = 1; donor < ntrain; ++donor) {
    if(memcmp(&row[donor], &row[0], sizeof(double)) != 0)
      return;
    if(activity != NULL && (donor & 1023) == 0) activity();
  }
  reuse->row_value = row[0];
  reuse->row_divisor = divisor;
  reuse->row_valid = 1;
}

static inline int np_inference_reuse_get(
  const NPInferenceReuse *reuse, int direction, const double *projection,
  double *standard_error, int *unavailable)
{
  if(!reuse->row_valid || direction < 0 ||
     direction >= reuse->ndirections || projection == NULL)
    return 0;
  const NPInferenceReuseEntry *entry = &reuse->entry[direction];
  if(!entry->valid ||
     memcmp(&entry->kernel_value, &reuse->row_value, sizeof(double)) != 0 ||
     memcmp(&entry->divisor, &reuse->row_divisor, sizeof(double)) != 0 ||
     memcmp(reuse->projection + (size_t)direction * (size_t)reuse->nterms,
            projection, (size_t)reuse->nterms * sizeof(double)) != 0)
    return 0;
  *standard_error = entry->standard_error;
  if(unavailable != NULL) *unavailable = entry->unavailable;
  return 1;
}

/* Only invoke after the existing donor calculation returned success. */
static inline void np_inference_reuse_put(
  NPInferenceReuse *reuse, int direction, const double *projection,
  double standard_error, int unavailable)
{
  if(!reuse->row_valid || direction < 0 ||
     direction >= reuse->ndirections || projection == NULL ||
     (unavailable != 0 && unavailable != 1) ||
     (!unavailable && !isfinite(standard_error)))
    return;
  for(int term = 0; term < reuse->nterms; ++term)
    if(!isfinite(projection[term])) return;
  NPInferenceReuseEntry *entry = &reuse->entry[direction];
  memcpy(reuse->projection + (size_t)direction * (size_t)reuse->nterms,
         projection, (size_t)reuse->nterms * sizeof(double));
  entry->kernel_value = reuse->row_value;
  entry->divisor = reuse->row_divisor;
  entry->standard_error = standard_error;
  entry->unavailable = unavailable;
  entry->valid = 1;
}

#endif
