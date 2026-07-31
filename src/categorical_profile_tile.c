#include <math.h>
#include <limits.h>
#include <stddef.h>

#include "categorical_profile_tile.h"
#include "headers.h"
#include "np_native_safety.h"

/*
 * These rank-local vector primitives own the package's categorical kernel
 * formulas. The bounded profile engine composes them but does not duplicate
 * estimator state or formula implementations.
 */
extern void np_ukernelv(const int KERNEL,
                        const double * const xt,
                        const int num_xt,
                        const int do_xw,
                        const double x,
                        const double lambda,
                        const int ncat,
                        double * const result,
                        const XL * const xl,
                        const int skip_upper_gate);
extern void np_okernelv(const int KERNEL,
                        const double * const xt,
                        const int num_xt,
                        const int do_xw,
                        const double x,
                        const double lambda,
                        const double *cats,
                        const int ncat,
                        double * const result,
                        const XL * const xl,
                        const int swap_xxt);
extern void np_convol_okernelv(const int KERNEL,
                               const double * const xt,
                               const int num_xt,
                               const int do_xw,
                               const double x,
                               const double lambda,
                               int ncat,
                               double *cat,
                               double * const result,
                               const int swap_xxt);

static NPCategoricalProfileTileStatus
np_categorical_profile_metadata_validate(
  const NPCategoricalProfileKernelSpec *spec)
{
  int i;
  int nvar;

  if(spec == NULL)
    return NP_PROFILE_TILE_ERR_ARGUMENT;
  if((spec->ntrain <= 0) || (spec->neval <= 0) ||
     (spec->nunordered < 0) || (spec->nordered < 0) ||
     (spec->nunordered > INT_MAX - spec->nordered))
    return NP_PROFILE_TILE_ERR_DIMENSION;
  nvar = spec->nunordered + spec->nordered;
  if(nvar <= 0)
    return NP_PROFILE_TILE_ERR_DIMENSION;
  if((spec->operator_code == NULL) || (spec->lambda == NULL) ||
     (spec->num_categories == NULL) ||
     (spec->category_values == NULL))
    return NP_PROFILE_TILE_ERR_ARGUMENT;
  if((spec->nunordered > 0) &&
     ((spec->train_unordered == NULL) ||
      (spec->eval_unordered == NULL) ||
      (spec->kernel_unordered == NULL)))
    return NP_PROFILE_TILE_ERR_ARGUMENT;
  if((spec->nordered > 0) &&
     ((spec->train_ordered == NULL) ||
      (spec->eval_ordered == NULL) ||
      (spec->kernel_ordered == NULL)))
    return NP_PROFILE_TILE_ERR_ARGUMENT;

  for(i = 0; i < spec->nunordered; i++) {
    if((spec->train_unordered[i] == NULL) ||
       (spec->eval_unordered[i] == NULL))
      return NP_PROFILE_TILE_ERR_ARGUMENT;
    if((spec->kernel_unordered[i] < 0) ||
       (spec->kernel_unordered[i] > 1))
      return NP_PROFILE_TILE_ERR_KERNEL;
    if((spec->operator_code[i] < OP_NORMAL) ||
       (spec->operator_code[i] > OP_INTEGRAL))
      return NP_PROFILE_TILE_ERR_OPERATOR;
    if(!isfinite(spec->lambda[i]))
      return NP_PROFILE_TILE_ERR_NONFINITE;
    if((spec->lambda[i] < 0.0) || (spec->lambda[i] > 1.0) ||
       (spec->num_categories[i] <= 0))
      return NP_PROFILE_TILE_ERR_SUPPORT;
  }

  for(i = 0; i < spec->nordered; i++) {
    const int offset = spec->nunordered + i;
    const int ncat = spec->num_categories[offset];
    double *cats = spec->category_values[offset];
    int j;

    if((spec->train_ordered[i] == NULL) ||
       (spec->eval_ordered[i] == NULL))
      return NP_PROFILE_TILE_ERR_ARGUMENT;
    if((spec->kernel_ordered[i] < 0) ||
       (spec->kernel_ordered[i] > 3))
      return NP_PROFILE_TILE_ERR_KERNEL;
    if((spec->operator_code[offset] < OP_NORMAL) ||
       (spec->operator_code[offset] > OP_INTEGRAL))
      return NP_PROFILE_TILE_ERR_OPERATOR;
    if(!isfinite(spec->lambda[offset]))
      return NP_PROFILE_TILE_ERR_NONFINITE;
    if((spec->lambda[offset] < 0.0) ||
       (spec->lambda[offset] > 1.0) ||
       (ncat <= 0) || (cats == NULL))
      return NP_PROFILE_TILE_ERR_SUPPORT;
    for(j = 0; j < ncat; j++) {
      if(!isfinite(cats[j]))
        return NP_PROFILE_TILE_ERR_NONFINITE;
      if((j > 0) && !(cats[j] > cats[j - 1]))
        return NP_PROFILE_TILE_ERR_SUPPORT;
    }
  }

  return NP_PROFILE_TILE_OK;
}

static NPCategoricalProfileTileStatus
np_categorical_profile_training_data_validate(
  const NPCategoricalProfileKernelSpec *spec)
{
  int i;

  for(i = 0; i < spec->nunordered; i++) {
    int j;

    for(j = 0; j < spec->ntrain; j++)
      if(!isfinite(spec->train_unordered[i][j]))
        return NP_PROFILE_TILE_ERR_NONFINITE;
  }

  for(i = 0; i < spec->nordered; i++) {
    int j;

    for(j = 0; j < spec->ntrain; j++)
      if(!isfinite(spec->train_ordered[i][j]))
        return NP_PROFILE_TILE_ERR_NONFINITE;
  }

  return NP_PROFILE_TILE_OK;
}

static NPCategoricalProfileTileStatus
np_categorical_profile_evaluation_data_validate(
  const NPCategoricalProfileKernelSpec *spec,
  int eval_start,
  int eval_count)
{
  int i;

  for(i = 0; i < spec->nunordered; i++) {
    int j;

    for(j = eval_start; j < eval_start + eval_count; j++)
      if(!isfinite(spec->eval_unordered[i][j]))
        return NP_PROFILE_TILE_ERR_NONFINITE;
  }

  for(i = 0; i < spec->nordered; i++) {
    int j;

    for(j = eval_start; j < eval_start + eval_count; j++)
      if(!isfinite(spec->eval_ordered[i][j]))
        return NP_PROFILE_TILE_ERR_NONFINITE;
  }

  return NP_PROFILE_TILE_OK;
}

NPCategoricalProfileTileStatus
np_categorical_profile_spec_validate(
  const NPCategoricalProfileKernelSpec *spec)
{
  NPCategoricalProfileTileStatus status =
    np_categorical_profile_metadata_validate(spec);

  if(status != NP_PROFILE_TILE_OK)
    return status;
  return np_categorical_profile_training_data_validate(spec);
}

NPCategoricalProfileTileStatus
np_categorical_profile_tile_bytes(size_t ntrain,
                                  size_t nrows,
                                  size_t *elements,
                                  size_t *bytes)
{
  size_t required_elements;
  size_t required_bytes;

  if((elements == NULL) || (bytes == NULL) ||
     (ntrain == 0U) || (nrows == 0U))
    return NP_PROFILE_TILE_ERR_ARGUMENT;
  if(!np_size_mul_checked(ntrain, nrows, &required_elements) ||
     !np_size_array_bytes_checked(required_elements,
                                  sizeof(double),
                                  &required_bytes))
    return NP_PROFILE_TILE_ERR_CAPACITY;
  *elements = required_elements;
  *bytes = required_bytes;
  return NP_PROFILE_TILE_OK;
}

NPCategoricalProfileTileStatus
np_categorical_profile_tile_fill_prevalidated(
  const NPCategoricalProfileKernelSpec *spec,
  int eval_start,
  int eval_count,
  double *output,
  size_t output_elements)
{
  NPCategoricalProfileTileStatus status;
  size_t required_elements;
  size_t required_bytes;
  int row;

  if(spec == NULL)
    return NP_PROFILE_TILE_ERR_ARGUMENT;
  if((eval_start < 0) || (eval_count <= 0) ||
     (eval_start > spec->neval) ||
     (eval_count > spec->neval - eval_start))
    return NP_PROFILE_TILE_ERR_RANGE;
  status = np_categorical_profile_evaluation_data_validate(spec,
                                                           eval_start,
                                                           eval_count);
  if(status != NP_PROFILE_TILE_OK)
    return status;
  status = np_categorical_profile_tile_bytes((size_t)spec->ntrain,
                                             (size_t)eval_count,
                                             &required_elements,
                                             &required_bytes);
  if(status != NP_PROFILE_TILE_OK)
    return status;
  if((output == NULL) || (output_elements < required_elements) ||
     (required_bytes > NP_CATEGORICAL_PROFILE_TILE_MAX_BYTES))
    return NP_PROFILE_TILE_ERR_CAPACITY;

  for(row = 0; row < eval_count; row++) {
    const int eval_index = eval_start + row;
    double *weights = output + (size_t)row*(size_t)spec->ntrain;
    int have_weights = 0;
    int i;

    for(i = 0; i < spec->nunordered; i++) {
      const int kernel =
        spec->kernel_unordered[i] +
        OP_UFUN_OFFSETS[spec->operator_code[i]];
      np_ukernelv(kernel,
                  spec->train_unordered[i],
                  spec->ntrain,
                  have_weights,
                  spec->eval_unordered[i][eval_index],
                  spec->lambda[i],
                  spec->num_categories[i],
                  weights,
                  NULL,
                  0);
      have_weights = 1;
    }

    for(i = 0; i < spec->nordered; i++) {
      const int offset = spec->nunordered + i;
      const int base_kernel = spec->kernel_ordered[i];
      const int operator_code = spec->operator_code[offset];

      /*
       * Preserve the incumbent ordered Li--Racine convolution owner. Other
       * ordered kernels use the established operator-code offsets.
       */
      if((operator_code == OP_CONVOLUTION) &&
         (spec->kernel_ordered[0] == 1)) {
        np_convol_okernelv(base_kernel,
                           spec->train_ordered[i],
                           spec->ntrain,
                           have_weights,
                           spec->eval_ordered[i][eval_index],
                           spec->lambda[offset],
                           spec->num_categories[offset],
                           spec->category_values[offset],
                           weights,
                           0);
      } else {
        const int kernel =
          base_kernel + OP_OFUN_OFFSETS[operator_code];
        np_okernelv(kernel,
                    spec->train_ordered[i],
                    spec->ntrain,
                    have_weights,
                    spec->eval_ordered[i][eval_index],
                    spec->lambda[offset],
                    spec->category_values[offset],
                    spec->num_categories[offset],
                    weights,
                    NULL,
                    0);
      }
      have_weights = 1;
    }
  }

  return NP_PROFILE_TILE_OK;
}

NPCategoricalProfileTileStatus
np_categorical_profile_tile_fill(
  const NPCategoricalProfileKernelSpec *spec,
  int eval_start,
  int eval_count,
  double *output,
  size_t output_elements)
{
  NPCategoricalProfileTileStatus status =
    np_categorical_profile_spec_validate(spec);

  if(status != NP_PROFILE_TILE_OK)
    return status;
  return np_categorical_profile_tile_fill_prevalidated(spec,
                                                       eval_start,
                                                       eval_count,
                                                       output,
                                                       output_elements);
}

const char *np_categorical_profile_tile_status_message(
  NPCategoricalProfileTileStatus status)
{
  switch(status) {
  case NP_PROFILE_TILE_OK:
    return "success";
  case NP_PROFILE_TILE_ERR_ARGUMENT:
    return "invalid categorical-profile tile argument";
  case NP_PROFILE_TILE_ERR_DIMENSION:
    return "invalid categorical-profile dimensions";
  case NP_PROFILE_TILE_ERR_KERNEL:
    return "unsupported categorical-profile kernel code";
  case NP_PROFILE_TILE_ERR_OPERATOR:
    return "unsupported categorical-profile operator code";
  case NP_PROFILE_TILE_ERR_NONFINITE:
    return "categorical-profile inputs must be finite";
  case NP_PROFILE_TILE_ERR_SUPPORT:
    return "invalid categorical-profile support or bandwidth";
  case NP_PROFILE_TILE_ERR_RANGE:
    return "invalid categorical-profile evaluation range";
  case NP_PROFILE_TILE_ERR_CAPACITY:
    return "categorical-profile tile exceeds checked capacity";
  default:
    return "unknown categorical-profile tile status";
  }
}
