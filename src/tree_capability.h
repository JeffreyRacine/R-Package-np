#ifndef NP_TREE_CAPABILITY_H
#define NP_TREE_CAPABILITY_H

typedef enum {
  NP_PRUNE_SUPPORT_UNKNOWN = 0,
  NP_PRUNE_SUPPORT_UNBOUNDED,
  NP_PRUNE_SUPPORT_COMPACT_ZERO,
  NP_PRUNE_SUPPORT_ONE_SIDED_TAIL,
  NP_PRUNE_SUPPORT_UNSUPPORTED
} NPPruneSupportKind;

typedef struct {
  NPPruneSupportKind kind;
  double lower;
  double upper;
} NPPruneSupportDescriptor;

typedef struct {
  double train_min;
  double train_max;
  double eval_min;
  double eval_max;
} NPTreeCoordinateEnvelope;

typedef enum {
  NP_TREE_CAPABILITY_UNKNOWN = 0,
  NP_TREE_CANNOT_PRUNE,
  NP_TREE_CAN_PRUNE_ZERO_SUPPORT,
  NP_TREE_CAN_USE_ONE_SIDED_TAIL
} NPTreeCapability;

NPPruneSupportDescriptor np_prune_support_descriptor(
  int effective_kernel,
  double lower,
  double upper);

NPTreeCapability np_tree_fixed_coordinate_capability(
  NPPruneSupportDescriptor support,
  double bandwidth,
  NPTreeCoordinateEnvelope envelope);

NPTreeCapability np_tree_capability_combine(
  NPTreeCapability current,
  NPTreeCapability coordinate);

NPTreeCapability np_tree_fixed_capability(
  const NPPruneSupportDescriptor *support,
  const double *bandwidth,
  const NPTreeCoordinateEnvelope *envelope,
  int dimensions);

#endif
